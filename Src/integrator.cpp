#include "pch.h"
#include "integrator.h"
#include "integration_params.h"
#include "JKFit.h"
#include "libCintMain.h"
#include "nos_math.h"

#if defined(__APPLE__)
// On macOS we are using Accelerate for BLAS/LAPACK
#include <Accelerate/Accelerate.h>
#else
// Linux/Windows with oneMKL
#include <mkl.h>
#endif

vec einsum_ijk_ij_p(const dMatrix3 &v1, const dMatrix2 &v2)
{
    const int I = (int)v1.extent(0);
    const int J = (int)v1.extent(1);
    const int P = (int)v1.extent(2);
    // Initialize the result vector
    vec rho(P, 0.0);

    // Perform the summation
    for (int p = 0; p < P; ++p)
    {
        for (int i = 0; i < I; ++i)
        {
            for (int j = 0; j < J; ++j)
            {
                rho[p] += v1(i, j, p) * v2(i, j);
            }
        }
    }
    return rho;
}

// Reorder p-Orbitals to the SALTED convention of:
//  L = 1 components following the - 1, 0, +1 convention
// Meaning for every p-Orbital we swap the first and last component
vec reorder_p(vec coefs_in, WFN aux_basis)
{
    vec coefs_out = coefs_in;
    int coef_idx = 0;
    for (int atm_idx = 0; atm_idx < aux_basis.get_ncen(); atm_idx++)
    {
        for (int shell = 0; shell < aux_basis.get_atom_shell_count(atm_idx); shell++)
        {
            int type = aux_basis.get_shell_type(atm_idx, shell); // Guessing only NoSpherA2 basis sets are used! l starts at 0!!!!!
            if (type != 1)
            {
                coef_idx += 2 * type + 1;
                continue;
            }
            coefs_out[coef_idx] = coefs_in[coef_idx + 1];
            coefs_out[coef_idx + 1] = coefs_in[coef_idx + 2];
            coefs_out[coef_idx + 2] = coefs_in[coef_idx];
            coef_idx += 3;
        }
    }

    return coefs_out;
}



// Enhanced electron restraint with adaptive weighting for s-orbitals only
// Only s-orbitals are restrained as other orbitals don't contribute to 
// spherically averaged electron density used in electron counting
void add_electron_restraint(vec& eri2c, vec& rho, const WFN& wavy_aux, 
                          double base_restraint_coef, 
                          bool adaptive_weighting,
                          const vec& expected_charges)
{
    double radial;
    basis_set_entry bf;
    
    size_t original_size = rho.size();
    size_t n_atoms = wavy_aux.get_ncen();
    
    // Adaptive restraint coefficient based on system size and basis set quality
    double restraint_coef = base_restraint_coef;
    if (adaptive_weighting) {
        // Scale restraint based on auxiliary basis size - larger basis needs less restraint
        double basis_scaling = std::max(0.1, 1.0 - std::log10(original_size) * 0.1);
        restraint_coef *= basis_scaling;
        
        // Scale based on number of atoms - more atoms need stronger restraint per atom
        double atom_scaling = std::min(2.0, 1.0 + std::sqrt(n_atoms) * 0.1);
        restraint_coef *= atom_scaling;
    }
    
    std::cout << "Setting adaptive restraint coefficient to: " << std::fixed 
              << std::showpoint << std::setprecision(6) << restraint_coef << std::endl;

    // Resize matrices for restraint equations
    size_t new_size = original_size * original_size + n_atoms * original_size;
    eri2c.resize(new_size);
    
    // Create matrix view for restraint block
    dMatrixRef2 eri2c_restraint(eri2c.data() + original_size * original_size, n_atoms, original_size);
    
    rho.resize(original_size + n_atoms);

    // Calculate electron integration matrix for all basis functions
    vec electron_integrals(original_size, 0.0);
    int coef_idx = 0;
    
    for (int atm_idx = 0; atm_idx < n_atoms; atm_idx++) {
        atom current_atom = wavy_aux.get_atoms()[atm_idx];
        
        // Determine expected electron count for this atom
        double expected_electrons = current_atom.get_charge(); // Default to nuclear charge
        if (!expected_charges.empty() && atm_idx < expected_charges.size()) {
            expected_electrons = expected_charges[atm_idx];
        }
        
        // Set target electron count with adaptive weighting
        double atom_weight = restraint_coef;
        if (adaptive_weighting) {
            // Heavier atoms (higher Z) get stronger restraint
            atom_weight *= std::min(2.0, 1.0 + current_atom.get_charge() * 0.02);
        }
        
        rho[original_size + atm_idx] = expected_electrons * atom_weight;
        
        // Reset coefficient index for this atom
        int type = -1, prim = 0;
        
        for (unsigned int shell = 0; shell < current_atom.get_shellcount().size(); shell++) {
            type = current_atom.get_basis_set_entry(prim).get_type();
            
            // Only apply restraints to s-orbitals (type == 0)
            // Other orbitals don't contribute to spherically averaged electron density
            if (type == 0) { // s-orbital only
                radial = 0.0;
                
                // Sum over primitives in this s-shell
                for (unsigned int e = 0; e < current_atom.get_shellcount()[shell]; e++) {
                    bf = current_atom.get_basis_set_entry(prim + e);
                    primitive p(0, bf.get_type(), bf.get_exponent(), bf.get_coefficient());
                    
                    // Electron-nucleus attraction integral for s-orbital: <χ|1|χ>
                    radial += constants::PI / (2.0 * std::pow(p.get_exp(), 1.5)) 
                            * p.normalization_constant() * p.get_coef();
                }
                
                // Store the electron integral for this s-orbital
                eri2c_restraint(atm_idx, coef_idx) = radial * atom_weight;
                coef_idx++;
            } else {
                // Skip non-s orbitals but still increment indices appropriately
                coef_idx += (2 * type + 1);
            }
            
            prim += current_atom.get_shellcount()[shell];
        }
    }
    
    std::cout << "Added electron restraints for " << n_atoms << " atoms." << std::endl;
}


vec density_fit_restrain(const WFN &wavy, const WFN& wavy_aux, const double max_mem, const char metric,
                 double restraint_strength, bool adaptive_restraint, 
                 const std::string& charge_scheme, bool analyze_quality)
{
    vec eri2c;
    vec eri3c;
    vec rho;

    // Initialize basis functions (qmBasis and auxBasis)
    Int_Params normal_basis(wavy);
    Int_Params aux_basis(wavy_aux);
    dMatrix2 dm = wavy.get_dm();

    std::cout << "\n=== Density Fitting with Enhanced Restraints ===" << std::endl;
    std::cout << "Normal basis functions: " << normal_basis.get_nao() << std::endl;
    std::cout << "Auxiliary basis functions: " << aux_basis.get_nao() << std::endl;
    std::cout << "Metric: " << (metric == 'C' ? "Coulomb" : "Overlap") << std::endl;

    // Compute integrals
    if (metric == 'C')
    {
        computeEri2c(aux_basis, eri2c);
        computeRho_Coulomb(normal_basis, aux_basis, dm, rho, max_mem);
    }
    else if (metric == 'O')
    {
        compute2c_Overlap(aux_basis, eri2c);
        computeRho_Overlap(normal_basis, aux_basis, dm, rho, max_mem);
    }

    // Calculate expected charges based on chosen scheme
    vec expected_charges = calculate_expected_charges(wavy, wavy_aux, charge_scheme);
    std::cout << "Using charge scheme: " << charge_scheme << std::endl;
    
    // Apply enhanced electron restraints
    add_electron_restraint(eri2c, rho, wavy_aux, restraint_strength, 
                          adaptive_restraint, expected_charges);
    
    // Solve the regularized system
    std::cout << "Solving regularized linear system..." << std::endl;
    solve_linear_system(eri2c, rho.size(), aux_basis.get_nao(), rho);

    // Reorder p-orbitals to SALTED convention
    rho = reorder_p(rho, wavy_aux);

    // Analyze the quality of the fit if requested
    if (analyze_quality) {
        analyze_density_fit_quality(rho, wavy_aux, expected_charges);
    }

    std::cout << "===============================================\n" << std::endl;
    return rho;
}

// Hybrid approach combining thikonov regularization and electron restraints
vec density_fit_hybrid(const WFN &wavy, const WFN& wavy_aux, const double max_mem, 
                      const char metric, double restraint_strength, 
                      double tikhonov_lambda, const std::string& charge_scheme, bool analyze_quality)
{
    vec eri2c;
    vec eri3c;
    vec rho;

    // Initialize basis functions
    Int_Params normal_basis(wavy);
    Int_Params aux_basis(wavy_aux);
    dMatrix2 dm = wavy.get_dm();

    std::cout << "\n=== Hybrid Density Fitting (Restraints + Tikhonov) ===" << std::endl;

    // Compute integrals
    if (metric == 'C')
    {
        computeEri2c(aux_basis, eri2c);
        computeRho_Coulomb(normal_basis, aux_basis, dm, rho, max_mem);
    }
    else if (metric == 'O')
    {
        compute2c_Overlap(aux_basis, eri2c);
        computeRho_Overlap(normal_basis, aux_basis, dm, rho, max_mem);
    }

    // First apply Tikhonov regularization to the original matrix
    size_t n = aux_basis.get_nao();
    for (size_t i = 0; i < n; i++) {
        eri2c[i * n + i] += tikhonov_lambda;
    }

    // Then add electron restraints
    vec expected_charges = calculate_expected_charges(wavy, wavy_aux, charge_scheme);
    add_electron_restraint(eri2c, rho, wavy_aux, restraint_strength, true, expected_charges);
    
    std::cout << "Solving hybrid regularized system..." << std::endl;
    solve_linear_system(eri2c, rho.size(), n, rho);

    rho = reorder_p(rho, wavy_aux);
    
    // Analyze the quality of the fit if requested
    if (analyze_quality) {
        analyze_density_fit_quality(rho, wavy_aux, expected_charges);
    }

    std::cout << "====================================================\n" << std::endl;
    return rho;
}

// Calculate expected atomic charges based on different partitioning schemes
vec calculate_expected_charges(const WFN& wavy, const WFN& wavy_aux, const std::string& scheme)
{
    vec expected_charges(wavy_aux.get_ncen());
    
    if (scheme == "nuclear") {
        // Simple nuclear charge
        for (int i = 0; i < wavy_aux.get_ncen(); i++) {
            expected_charges[i] = wavy_aux.get_atoms()[i].get_charge();
        }
    }
    else if (scheme == "mulliken_estimate") {
        // Rough Mulliken-like estimate based on electronegativity
        // This is a simplified approach - in practice you'd want proper Mulliken analysis
        double total_electrons = 0.0;
        vec electronegativities(wavy_aux.get_ncen());
        
        for (int i = 0; i < wavy_aux.get_ncen(); i++) {
            double Z = wavy_aux.get_atoms()[i].get_charge();
            total_electrons += Z;
            
            // Approximate electronegativity (Pauling scale)
            if (Z == 1) electronegativities[i] = 2.20;      // H
            else if (Z == 6) electronegativities[i] = 2.55;  // C
            else if (Z == 7) electronegativities[i] = 3.04;  // N
            else if (Z == 8) electronegativities[i] = 3.44;  // O
            else if (Z == 9) electronegativities[i] = 3.98;  // F
            else electronegativities[i] = 2.0 + Z * 0.1;    // Rough estimate
        }
        
        // Redistribute electrons based on electronegativity differences
        double avg_electronegativity = 0.0;
        for (double en : electronegativities) avg_electronegativity += en;
        avg_electronegativity /= electronegativities.size();
        
        for (int i = 0; i < wavy_aux.get_ncen(); i++) {
            double Z = wavy_aux.get_atoms()[i].get_charge();
            double en_diff = electronegativities[i] - avg_electronegativity;
            // More electronegative atoms get slightly more electrons
            expected_charges[i] = Z + en_diff * 0.1;
        }
    }
    else if (scheme == "formal_charges") {
        // Use formal charges if available (would need to be passed in)
        // For now, default to nuclear charges
        for (int i = 0; i < wavy_aux.get_ncen(); i++) {
            expected_charges[i] = wavy_aux.get_atoms()[i].get_charge();
        }
    }
    
    return expected_charges;
}

// Analyze the quality of density fitting and detect problematic charges
void analyze_density_fit_quality(const vec& coefficients, const WFN& wavy_aux, 
                                const vec& expected_charges)
{
    std::cout << "\n=== Density Fitting Quality Analysis ===" << std::endl;
    
    vec atomic_populations(wavy_aux.get_ncen(), 0.0);
    int coef_idx = 0;
    
    // Calculate atomic populations from coefficients using only s-orbitals
    for (int atm_idx = 0; atm_idx < wavy_aux.get_ncen(); atm_idx++) {
        atom current_atom = wavy_aux.get_atoms()[atm_idx];
        
        int type = -1, prim = 0;
        for (unsigned int shell = 0; shell < current_atom.get_shellcount().size(); shell++) {
            type = current_atom.get_basis_set_entry(prim).get_type();
            
            // Only calculate population from s-orbitals (type == 0)
            if (type == 0) { // s-orbital
                double radial_integral = 0.0;
                basis_set_entry bf;
                
                // Calculate the radial integration for s-orbital
                for (unsigned int e = 0; e < current_atom.get_shellcount()[shell]; e++) {
                    bf = current_atom.get_basis_set_entry(prim + e);
                    primitive p(0, bf.get_type(), bf.get_exponent(), bf.get_coefficient());
                    
                    // Radial integral for s-orbital: <χ|1|χ>
                    radial_integral += constants::PI / (2.0 * std::pow(p.get_exp(), 1.5)) 
                                     * p.normalization_constant() * p.get_coef();
                }
                
                // Multiply radial part with the coefficient to get electron population contribution
                atomic_populations[atm_idx] += coefficients[coef_idx] * radial_integral;
                coef_idx++;
            } else {
                // Skip non-s orbitals but still increment coefficient index
                coef_idx += (2 * type + 1);
            }
            
            prim += current_atom.get_shellcount()[shell];
        }
        
        double n_electrons = current_atom.get_charge();
        double computed_charge = n_electrons - atomic_populations[atm_idx];
        double expected_charge = 0.0;
        
        if (!expected_charges.empty() && atm_idx < expected_charges.size()) {
            expected_charge = n_electrons - expected_charges[atm_idx];
        }
        
        std::cout << "Atom " << atm_idx + 1 << std::fixed << std::setprecision(3) << " (Z=" << n_electrons << "): "
                  << "Population = " << atomic_populations[atm_idx]
                  << ", Charge = " << computed_charge
                  << ", Expected = " << expected_charge
                  << ", Deviation = " << std::abs(computed_charge - expected_charge) << std::endl;
                  if (std::abs(computed_charge - expected_charge) > 1) {
                      std::cout << "Warning: Significant deviation for atom " << atm_idx + 1 << "!" << std::endl;
                  }
    }
}

// Example usage function demonstrating the enhanced density fitting approaches
void demonstrate_enhanced_density_fitting(const WFN& wavy, const WFN& wavy_aux)
{
    std::cout << "\n=== Enhanced Density Fitting Demonstration ===" << std::endl;
    
    double max_mem = 1000.0; // MB
    char metric = 'C'; // Coulomb metric
    
    // Method 0: Unrestrained (original approach)
    std::cout << "\n--- Method 0: Unrestrained (Baseline) ---" << std::endl;
    vec coeff_unrestrained = density_fit_unrestrained(wavy, wavy_aux, max_mem, metric, true);
    
    // Method 1: Enhanced restraints with adaptive weighting
    std::cout << "\n--- Method 1: Enhanced Adaptive Restraints ---" << std::endl;
    vec coeff_enhanced = density_fit_restrain(wavy, wavy_aux, max_mem, metric,
                                   0.0001,     // restraint strength
                                   true,      // adaptive weighting
                                   "mulliken_estimate", // charge scheme
                                   true);     // analyze quality
    
    // Method 2: Hybrid approach
    std::cout << "\n--- Method 2: Hybrid Regularization ---" << std::endl;
    vec coeff_hybrid = density_fit_hybrid(wavy, wavy_aux, max_mem, metric,
                                         0.0001,     // restraint strength
                                         1e-6,      // tikhonov lambda
                                         "mulliken_estimate"); // charge scheme
    
    // Compare results
    std::cout << "\n=== Method Comparison ===" << std::endl;
    std::cout << "Unrestrained        - Max coeff: " << std::fixed << std::setprecision(6) 
              << *std::max_element(coeff_unrestrained.begin(), coeff_unrestrained.end()) 
              << ", Min coeff: " << *std::min_element(coeff_unrestrained.begin(), coeff_unrestrained.end()) << std::endl;
    std::cout << "Enhanced restraints - Max coeff: " << std::fixed << std::setprecision(6) 
              << *std::max_element(coeff_enhanced.begin(), coeff_enhanced.end()) 
              << ", Min coeff: " << *std::min_element(coeff_enhanced.begin(), coeff_enhanced.end()) << std::endl;
    std::cout << "Hybrid approach     - Max coeff: " << std::fixed << std::setprecision(6) 
              << *std::max_element(coeff_hybrid.begin(), coeff_hybrid.end()) 
              << ", Min coeff: " << *std::min_element(coeff_hybrid.begin(), coeff_hybrid.end()) << std::endl;
    
    // Calculate RMS differences from unrestrained baseline
    double rms_unrestrained_vs_enhanced = 0.0;
    double rms_unrestrained_vs_hybrid = 0.0;
    double rms_enhanced_vs_hybrid = 0.0;
    
    size_t n_coeff = std::min({coeff_unrestrained.size(), coeff_enhanced.size(), coeff_hybrid.size()});
    
    for (size_t i = 0; i < n_coeff; i++) {
        double diff1 = coeff_unrestrained[i] - coeff_enhanced[i];
        double diff2 = coeff_unrestrained[i] - coeff_hybrid[i];
        double diff3 = coeff_enhanced[i] - coeff_hybrid[i];
        
        rms_unrestrained_vs_enhanced += diff1 * diff1;
        rms_unrestrained_vs_hybrid += diff2 * diff2;
        rms_enhanced_vs_hybrid += diff3 * diff3;
    }
    
    rms_unrestrained_vs_enhanced = std::sqrt(rms_unrestrained_vs_enhanced / n_coeff);
    rms_unrestrained_vs_hybrid = std::sqrt(rms_unrestrained_vs_hybrid / n_coeff);
    rms_enhanced_vs_hybrid = std::sqrt(rms_enhanced_vs_hybrid / n_coeff);
    
    std::cout << "\n=== RMS Coefficient Differences ===" << std::endl;
    std::cout << "Unrestrained vs Enhanced: " << std::fixed << std::setprecision(6) << rms_unrestrained_vs_enhanced << std::endl;
    std::cout << "Unrestrained vs Hybrid:   " << std::fixed << std::setprecision(6) << rms_unrestrained_vs_hybrid << std::endl;
    std::cout << "Enhanced vs Hybrid:       " << std::fixed << std::setprecision(6) << rms_enhanced_vs_hybrid << std::endl;
    
    // Calculate coefficient magnitude statistics
    double avg_unrestrained = 0.0, avg_enhanced = 0.0, avg_hybrid = 0.0;
    for (size_t i = 0; i < n_coeff; i++) {
        avg_unrestrained += std::abs(coeff_unrestrained[i]);
        avg_enhanced += std::abs(coeff_enhanced[i]);
        avg_hybrid += std::abs(coeff_hybrid[i]);
    }
    avg_unrestrained /= n_coeff;
    avg_enhanced /= n_coeff;
    avg_hybrid /= n_coeff;
    
    std::cout << "\n=== Average Coefficient Magnitudes ===" << std::endl;
    std::cout << "Unrestrained: " << std::fixed << std::setprecision(6) << avg_unrestrained << std::endl;
    std::cout << "Enhanced:     " << std::fixed << std::setprecision(6) << avg_enhanced << std::endl;
    std::cout << "Hybrid:       " << std::fixed << std::setprecision(6) << avg_hybrid << std::endl;
    
    std::cout << "=================================================\n" << std::endl;
}

// Basic unrestrained density fitting (original approach)
vec density_fit_unrestrained(const WFN &wavy, const WFN& wavy_aux, const double max_mem, 
                             const char metric, bool analyze_quality)
{
    vec eri2c;
    vec eri3c;
    vec rho;

    // Initialize basis functions
    Int_Params normal_basis(wavy);
    Int_Params aux_basis(wavy_aux);
    dMatrix2 dm = wavy.get_dm();

    std::cout << "\n=== Unrestrained Density Fitting (Original) ===" << std::endl;
    std::cout << "Normal basis functions: " << normal_basis.get_nao() << std::endl;
    std::cout << "Auxiliary basis functions: " << aux_basis.get_nao() << std::endl;
    std::cout << "Metric: " << (metric == 'C' ? "Coulomb" : "Overlap") << std::endl;

    // Compute integrals
    if (metric == 'C')
    {
        computeEri2c(aux_basis, eri2c);
        computeRho_Coulomb(normal_basis, aux_basis, dm, rho, max_mem);
    }
    else if (metric == 'O')
    {
        compute2c_Overlap(aux_basis, eri2c);
        computeRho_Overlap(normal_basis, aux_basis, dm, rho, max_mem);
    }

    std::cout << "Solving unrestrained linear system..." << std::endl;
    solve_linear_system(eri2c, aux_basis.get_nao(), rho);

    // Reorder p-orbitals to SALTED convention
    rho = reorder_p(rho, wavy_aux);

    // Analyze quality if requested
    if (analyze_quality) {
        vec empty_charges; // No expected charges for unrestrained case
        analyze_density_fit_quality(rho, wavy_aux, empty_charges);
    }

    std::cout << "==============================================\n" << std::endl;
    return rho;
}

