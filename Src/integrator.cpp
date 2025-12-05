#include "pch.h"
#include "integrator.h"
#include "integration_params.h"
#include "JKFit.h"
#include "libCintMain.h"
#include "nos_math.h"
#include "GridManager.h"

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

// Basic unrestrained density fitting (original approach)
vec density_fit_unrestrained(const WFN& wavy, const WFN& wavy_aux, const double max_mem,
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
    switch (metric) {
    case 'C':
        compute2C<Coulomb2C>(aux_basis, eri2c);
        computeRho<Coulomb3C>(normal_basis, aux_basis, dm, rho, max_mem);
        break;
    case 'O':
        compute2C<Overlap2C>(aux_basis, eri2c);
        computeRho<Overlap3C>(normal_basis, aux_basis, dm, rho, max_mem);
        break;
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


    switch (metric) {
        case 'C':
            compute2C<Coulomb2C>(aux_basis, eri2c);
            computeRho<Coulomb3C>(normal_basis, aux_basis, dm, rho, max_mem);
            break;
        case 'O':
            compute2C<Overlap2C>(aux_basis, eri2c);
            computeRho<Overlap3C>(normal_basis, aux_basis, dm, rho, max_mem);
            break;
    }

    // Calculate expected charges based on chosen scheme
    std::cout << "Using charge scheme: " << charge_scheme << std::endl;
    vec expected_populations = calculate_expected_populations(wavy, wavy_aux, charge_scheme);
    
    // Apply enhanced electron restraints
    add_electron_restraint(eri2c, rho, wavy_aux, restraint_strength, 
                          adaptive_restraint, expected_populations);
    
    // Solve the regularized 
    std::cout << "Solving regularized linear system..." << std::endl;
    solve_linear_system(eri2c, rho.size(), aux_basis.get_nao(), rho);

    // Reorder p-orbitals to SALTED convention
    rho = reorder_p(rho, wavy_aux);

    // Analyze the quality of the fit if requested
    if (analyze_quality) {
        analyze_density_fit_quality(rho, wavy_aux, expected_populations);
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
    switch (metric) {
    case 'C':
        compute2C<Coulomb2C>(aux_basis, eri2c);
        computeRho<Coulomb3C>(normal_basis, aux_basis, dm, rho, max_mem);
        break;
    case 'O':
        compute2C<Overlap2C>(aux_basis, eri2c);
        computeRho<Overlap3C>(normal_basis, aux_basis, dm, rho, max_mem);
        break;
    }

    // First apply Tikhonov regularization to the original matrix
    size_t n = aux_basis.get_nao();
    for (size_t i = 0; i < n; i++) {
        eri2c[i * n + i] += tikhonov_lambda;
    }

    // Then add electron restraints
    vec expected_populations = calculate_expected_populations(wavy, wavy_aux, charge_scheme);
    add_electron_restraint(eri2c, rho, wavy_aux, restraint_strength, true, expected_populations);
    
    std::cout << "Solving hybrid regularized system..." << std::endl;
    solve_linear_system(eri2c, rho.size(), n, rho);

    rho = reorder_p(rho, wavy_aux);
    
    // Analyze the quality of the fit if requested
    if (analyze_quality) {
        analyze_density_fit_quality(rho, wavy_aux, expected_populations);
    }

    std::cout << "====================================================\n" << std::endl;
    return rho;
}

// Calculate expected atomic populations based on different partitioning schemes
vec calculate_expected_populations(const WFN& wavy, const WFN& wavy_aux, const std::string& scheme)
{
    vec expected_populations(wavy_aux.get_ncen());
    
    if (scheme == "nuclear") {
        // Simple nuclear populations
        for (int i = 0; i < wavy_aux.get_ncen(); i++) {
            expected_populations[i] = wavy_aux.get_atoms()[i].get_charge();
        }
    }
    // https://pubs.acs.org/doi/10.1021/ed065p227
    else if (scheme == "sanderson_estimate") {
        double compound_electronegativity = 1.0;
        for (const auto& atom : wavy.get_atoms()) {
            compound_electronegativity *= constants::allen_electronegativities[atom.get_charge()-1];
        }
        compound_electronegativity = std::pow(compound_electronegativity, 1.0 / wavy.get_ncen());

        for (int iat = 0; iat < wavy_aux.get_ncen(); iat++) {
            double atom_electronegativity = constants::allen_electronegativities[wavy_aux.get_atoms()[iat].get_charge() - 1];
            expected_populations[iat] = wavy_aux.get_atoms()[iat].get_charge() + (compound_electronegativity - atom_electronegativity) / (1.57 * std::sqrt(atom_electronegativity));
        }
    }
    else if (scheme == "mulliken") {
        dMatrix2 dm = wavy.get_dm();
        vec eri2c;
        Int_Params normal_basis(wavy);
        compute2C<Overlap2C>(normal_basis, eri2c);
        dMatrixRef2 eri2c_ref(eri2c.data(), normal_basis.get_nao(), normal_basis.get_nao());
        const size_t nao = dm.extent(1);

        // optional: symmetric Mulliken operator M = 1/2 (P S + S P)
        // (avoid forming full matrices if memory is tight; just accumulate the diagonal)
        vec mulliken_pop(wavy.get_ncen(), 0.0);

        size_t mu_begin = 0;
        for (unsigned iat = 0; iat < wavy.get_ncen(); ++iat) {
            const atom A = wavy.get_atoms()[iat];

            // determine how many AOs belong to this atom (use contracted shells, not primitives!)
            size_t nAO_A = 0;
            int prim = 0;
            for (size_t sh = 0; sh < A.get_shellcount().size(); ++sh) {
                int l = A.get_basis_set_entry(prim).get_type() - 1;
                nAO_A += size_t(2 * l + 1);
                prim += A.get_shellcount()[sh];
            }
            size_t mu_end = mu_begin + nAO_A;

            double GA = 0.0;
            for (size_t m = mu_begin; m < mu_end; ++m) {
                double diag_PS = 0.0;
                for (size_t n = 0; n < nao; ++n)
                    diag_PS += dm(m, n) * eri2c_ref(n, m);
                GA += diag_PS;
            }
            expected_populations[iat] = GA;//A.get_charge() - GA;        // Mulliken charge
            mu_begin = mu_end;
        }
    }
    else if (scheme == "TFVC" || scheme == "Hirsh") {
        PartitionType type = (scheme == "TFVC") ? PartitionType::TFVC : PartitionType::Hirshfeld;
        GridConfiguration config;
        config.accuracy = 0;
        config.partition_type = type;
        config.pbc = 0;
        config.debug = false;

        cell unit_cell = cell();
        ivec asym_atom_list(wavy.get_ncen());
        for (int atom_nr = 0; atom_nr < wavy.get_ncen(); atom_nr++) {
            asym_atom_list[atom_nr] = atom_nr;
        }
        svec labels(wavy.get_ncen());
        for (int i = 0; i < wavy.get_ncen(); i++) {
            labels[i] = wavy.get_atoms()[i].get_label();
        }

        GridManager grid_manager(config);

        WFN temp = wavy;
        temp.delete_unoccupied_MOs();
        // Setup grids for the molecule
        bvec needs_grid = GridManager::determineAtomsNeedingGrids(temp, asym_atom_list);
        grid_manager.setup3DGridsForMolecule(temp, needs_grid, asym_atom_list, unit_cell);
        

        // Calculate partitioned charges
        auto results = grid_manager.calculatePartitionedCharges(temp, unit_cell);
        //results.printChargeTable(labels, temp, std::cout);
        for (int i = 0; i < temp.get_ncen(); i++) {
            expected_populations[i] = results.atom_charges[static_cast<int>(type)][i];
        }
    }
    else {
        std::cerr << "Warning: Unknown charge scheme '" << scheme 
            << "'. Defaulting to nuclear charges." << std::endl;
        expected_populations = calculate_expected_populations(wavy, wavy_aux, "nuclear");
    }

    return expected_populations;
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


#include "SALTED_utilities.h"
//#include "test_functions.h"
// Example usage function demonstrating the enhanced density fitting approaches
void demonstrate_enhanced_density_fitting(const WFN& wavy, const WFN& wavy_aux)
{
    //#include "test_functions.h"
    std::cout << "\n=== Enhanced Density Fitting Demonstration ===" << std::endl;

    double max_mem = 4000.0; // MB
    char metric = 'C'; // Coulomb metric

    // Method 0: Unrestrained (original approach)
    std::cout << "\n--- Method 0: Unrestrained (Baseline) ---" << std::endl;
    vec coeff_unrestrained = density_fit_unrestrained(wavy, wavy_aux, max_mem, metric, true);

    // Method 1: Enhanced restraints with adaptive weighting
    std::cout << "\n--- Method 1: Enhanced Adaptive Restraints ---" << std::endl;
    vec coeff_enhanced = density_fit_restrain(wavy, wavy_aux, max_mem, metric,
        2e-4,     // restraint strength
        true,      // adaptive weighting
        "Hirsh", // charge scheme
        true);     // analyze quality

    // Method 2: Hybrid approach
    std::cout << "\n--- Method 2: Hybrid Regularization ---" << std::endl;
    vec coeff_hybrid = density_fit_hybrid(wavy, wavy_aux, max_mem, metric,
        2e-4,     // restraint strength
        1e-6,      // tikhonov lambda
        "Hirsh", true); // charge scheme
    const double radius = 3.0;
    const double increment = 0.1;

    //std::vector<vec> all_coeffs = { coeff_unrestrained , coeff_enhanced, coeff_hybrid};

    GridConfiguration config;
    config.accuracy = 2;
    config.partition_type = PartitionType::Hirshfeld;
    config.pbc = 0;
    config.debug = false;
    config.all_charges = true;
    const int weight_index = GridData::HIRSH_WEIGHT;

    cell unit_cell = cell();
    ivec asym_atom_list(wavy.get_ncen());
    for (int atom_nr = 0; atom_nr < wavy.get_ncen(); atom_nr++) {
        asym_atom_list[atom_nr] = atom_nr;
    }

    GridManager grid_manager(config);
    WFN temp = wavy;

    temp.delete_unoccupied_MOs();
    // Setup grids for the molecule
    bvec needs_grid = GridManager::determineAtomsNeedingGrids(temp, asym_atom_list);
    grid_manager.setup3DGridsForMolecule(temp, needs_grid, asym_atom_list, unit_cell);
    GridData grid_data = grid_manager.getGridData();

    enum DiffDensityIndex { DIFF_UNRESTRAINED = 0, DIFF_ENHANCED = 1, DIFF_HYBRID = 2 };
    enum SumIndex { SUM_NO_DIFF = 0, RRS = 1, ABS_SUM = 2, SUM = 3 };
    vec3 diff_densities(wavy.get_ncen(), vec2(3, vec(4)));
    std::vector<atom> atoms = wavy_aux.get_atoms();
    vec partitioned_densities(wavy.get_ncen(), 0.0);

    for (int i = 0; i < wavy.get_ncen(); i++) {


        auto calc_density_unrestrained = [&](double x, double y, double z) {
            return calc_density_ML(x, y, z, coeff_unrestrained, atoms, i);
            };
        auto calc_density_enhanced = [&](double x, double y, double z) {
            return calc_density_ML(x, y, z, coeff_enhanced, atoms, i);
            };
        auto calc_density_hybrid = [&](double x, double y, double z) {
            return calc_density_ML(x, y, z, coeff_hybrid, atoms, i);
            };


        double diff_pos = 0, diff_neg = 0;
        const int natom_points = grid_data.num_points_per_atom[i];
        vec riDensity = grid_manager.evaluateFunctionOnGrid(grid_data.atomic_grids[i], calc_density_unrestrained);
        for (int p = 0; p < natom_points; p++) {
            diff_densities[i][DiffDensityIndex::DIFF_UNRESTRAINED][SumIndex::SUM_NO_DIFF] += riDensity[p];

            double val = riDensity[p] - (grid_data.atomic_grids[i][weight_index][p] * grid_data.atomic_grids[i][GridData::WFN_DENSITY][p]);
            diff_pos += std::abs(riDensity[p] + (grid_data.atomic_grids[i][weight_index][p] * grid_data.atomic_grids[i][GridData::WFN_DENSITY][p]));
            diff_neg += std::abs(val);

            diff_densities[i][DiffDensityIndex::DIFF_UNRESTRAINED][SumIndex::SUM] += val;
            diff_densities[i][DiffDensityIndex::DIFF_UNRESTRAINED][SumIndex::ABS_SUM] += std::abs(val) / 2;
            partitioned_densities[i] += grid_data.atomic_grids[i][weight_index][p] * grid_data.atomic_grids[i][GridData::WFN_DENSITY][p];

        }
        diff_densities[i][DiffDensityIndex::DIFF_UNRESTRAINED][SumIndex::RRS] = diff_neg / diff_pos;

        riDensity = grid_manager.evaluateFunctionOnGrid(grid_data.atomic_grids[i], calc_density_enhanced);
        diff_pos = 0, diff_neg = 0;
        for (int p = 0; p < natom_points; p++) {
            diff_densities[i][DiffDensityIndex::DIFF_ENHANCED][SumIndex::SUM_NO_DIFF] += riDensity[p];

            double val = riDensity[p] - (grid_data.atomic_grids[i][weight_index][p] * grid_data.atomic_grids[i][GridData::WFN_DENSITY][p]);
            diff_pos += std::abs(riDensity[p] + (grid_data.atomic_grids[i][weight_index][p] * grid_data.atomic_grids[i][GridData::WFN_DENSITY][p]));
            diff_neg += std::abs(val);

            diff_densities[i][DiffDensityIndex::DIFF_ENHANCED][SumIndex::SUM] += val;
            diff_densities[i][DiffDensityIndex::DIFF_ENHANCED][SumIndex::ABS_SUM] += std::abs(val) / 2;

        }
        diff_densities[i][DiffDensityIndex::DIFF_ENHANCED][SumIndex::RRS] = diff_neg / diff_pos;

        riDensity = grid_manager.evaluateFunctionOnGrid(grid_data.atomic_grids[i], calc_density_hybrid);
        diff_pos = 0, diff_neg = 0;
        for (int p = 0; p < natom_points; p++) {
            diff_densities[i][DiffDensityIndex::DIFF_HYBRID][SumIndex::SUM_NO_DIFF] += riDensity[p];

            double val = riDensity[p] - (grid_data.atomic_grids[i][weight_index][p] * grid_data.atomic_grids[i][GridData::WFN_DENSITY][p]);
            diff_pos += std::abs(riDensity[p] + (grid_data.atomic_grids[i][weight_index][p] * grid_data.atomic_grids[i][GridData::WFN_DENSITY][p]));
            diff_neg += std::abs(val);

            diff_densities[i][DiffDensityIndex::DIFF_HYBRID][SumIndex::SUM] += val;
            diff_densities[i][DiffDensityIndex::DIFF_HYBRID][SumIndex::ABS_SUM] += std::abs(val) / 2;

        }
        diff_densities[i][DiffDensityIndex::DIFF_HYBRID][SumIndex::RRS] = diff_neg / diff_pos;
    }
    std::cout << "\n=======================Unrestrained==========================" << std::endl;
    for (int i = 0; i < wavy.get_ncen(); i++) {
        std::string label = wavy.get_atoms()[i].get_label();
        //std::string out = std::format("Atom {:<3}({:<3}) | Unrestrained: RRS = {:>10.6f}, DiffSum = {:>10.6f}, |DiffSum/2| = {:>10.6f}, Sum = {:>10.6f}, Partitioned Sum = {:>10.6f}",
        //    label,
        //    i,
        //    diff_densities[i][DiffDensityIndex::DIFF_UNRESTRAINED][SumIndex::RRS],
        //    diff_densities[i][DiffDensityIndex::DIFF_UNRESTRAINED][SumIndex::SUM],
        //    diff_densities[i][DiffDensityIndex::DIFF_UNRESTRAINED][SumIndex::ABS_SUM],
        //    diff_densities[i][DiffDensityIndex::DIFF_UNRESTRAINED][SumIndex::SUM_NO_DIFF],
        //    partitioned_densities[i]
        //);
        std::string out = "Atom " + label + "(" + std::to_string(i) + ") | Unrestrained: RRS = " 
            + std::to_string(diff_densities[i][DiffDensityIndex::DIFF_UNRESTRAINED][SumIndex::RRS]) 
            + ", DiffSum = " + std::to_string(diff_densities[i][DiffDensityIndex::DIFF_UNRESTRAINED][SumIndex::SUM]) 
            + ", |DiffSum/2| = " + std::to_string(diff_densities[i][DiffDensityIndex::DIFF_UNRESTRAINED][SumIndex::ABS_SUM]) 
            + ", Sum = " + std::to_string(diff_densities[i][DiffDensityIndex::DIFF_UNRESTRAINED][SumIndex::SUM_NO_DIFF]);
        std::cout << out << std::endl;
    }
    std::cout << "\n=======================Enhanced==========================" << std::endl;
    for (int i = 0; i < wavy.get_ncen(); i++) {
        std::string label = wavy.get_atoms()[i].get_label();
        //std::string out = std::format("Atom {:<3}({:<3}) | Enhanced: RRS = {:>10.6f}, DiffSum = {:>10.6f}, |DiffSum/2| = {:>10.6f}, Sum = {:>10.6f}",
        //    label,
        //    i,
        //    diff_densities[i][DiffDensityIndex::DIFF_ENHANCED][SumIndex::RRS],
        //    diff_densities[i][DiffDensityIndex::DIFF_ENHANCED][SumIndex::SUM],
        //    diff_densities[i][DiffDensityIndex::DIFF_ENHANCED][SumIndex::ABS_SUM],
        //    diff_densities[i][DiffDensityIndex::DIFF_ENHANCED][SumIndex::SUM_NO_DIFF]
        //);
        std::string out = "Atom " + label + "(" + std::to_string(i) + ") | Enhanced: RRS = "
            + std::to_string(diff_densities[i][DiffDensityIndex::DIFF_ENHANCED][SumIndex::RRS])
            + ", DiffSum = " + std::to_string(diff_densities[i][DiffDensityIndex::DIFF_ENHANCED][SumIndex::SUM])
            + ", |DiffSum/2| = " + std::to_string(diff_densities[i][DiffDensityIndex::DIFF_ENHANCED][SumIndex::ABS_SUM])
            + ", Sum = " + std::to_string(diff_densities[i][DiffDensityIndex::DIFF_ENHANCED][SumIndex::SUM_NO_DIFF]);
        std::cout << out << std::endl;
    }
    std::cout << "\n=======================Hybrid==========================" << std::endl;
    for (int i = 0; i < wavy.get_ncen(); i++) {
        std::string label = wavy.get_atoms()[i].get_label();
        //std::string out = std::format("Atom {:<3}({:<3}) | Hybrid: RRS = {:>10.6f}, DiffSum = {:>10.6f}, |DiffSum/2| = {:>10.6f}, Sum = {:>10.6f}",
        //    label,
        //    i,
        //    diff_densities[i][DiffDensityIndex::DIFF_HYBRID][SumIndex::RRS],
        //    diff_densities[i][DiffDensityIndex::DIFF_HYBRID][SumIndex::SUM],
        //    diff_densities[i][DiffDensityIndex::DIFF_HYBRID][SumIndex::ABS_SUM],
        //    diff_densities[i][DiffDensityIndex::DIFF_HYBRID][SumIndex::SUM_NO_DIFF]
        //);
        std::string out = "Atom " + label + "(" + std::to_string(i) + ") | Hybrid: RRS = "
            + std::to_string(diff_densities[i][DiffDensityIndex::DIFF_HYBRID][SumIndex::RRS])
            + ", DiffSum = " + std::to_string(diff_densities[i][DiffDensityIndex::DIFF_HYBRID][SumIndex::SUM])
            + ", |DiffSum/2| = " + std::to_string(diff_densities[i][DiffDensityIndex::DIFF_HYBRID][SumIndex::ABS_SUM])
            + ", Sum = " + std::to_string(diff_densities[i][DiffDensityIndex::DIFF_HYBRID][SumIndex::SUM_NO_DIFF]);
        std::cout << out << std::endl;
    }

    //WFN dummy = wavy;
    //WFN dummy_aux = wavy_aux;
    //double MinMax[6]{ 0, 0, 0, 0, 0, 0 };
    //int steps[3]{ 0, 0, 0 };
    //readxyzMinMax_fromWFN(dummy, MinMax, steps, radius, increment, true);
    //dummy.delete_unoccupied_MOs();
    ////Cubes: 0=WFN, 1=RI_Unrestrained, 2=RI_Enhanced, 3=RI_Hybrid
    //std::vector<cube> cubes;
    //for (int a = 0; a < 4; a++) {
    //    cubes.push_back(cube(steps[0], steps[1], steps[2], dummy.get_ncen(), true));
    //    for (int i = 0; i < 3; i++)
    //    {
    //        cubes[a].set_origin(i, MinMax[i]);
    //        cubes[a].set_vector(i, i, (MinMax[i + 3] - MinMax[i]) / steps[i]);
    //    }
    //}
    //cubes[0].give_parent_wfn(dummy);
    //cubes[1].give_parent_wfn(dummy_aux);
    //cubes[2].give_parent_wfn(dummy_aux);
    //cubes[3].give_parent_wfn(dummy_aux);
    //std::cout << "Starting work..." << std::endl;
    //cubes[0].set_comment1("Calculated density using NoSpherA2 for WFN");
    //cubes[1].set_comment1("Calculated density using SALTED RI (Unrestrained)");
    //cubes[2].set_comment1("Calculated density using SALTED RI (Enhanced Restraints)");
    //cubes[3].set_comment1("Calculated density using SALTED RI (Hybrid Approach)");
    //cubes[0].set_path((dummy.get_path().parent_path() / dummy.get_path().stem()).string() + "_rho_WFN.cube");
    //cubes[1].set_path((dummy.get_path().parent_path() / dummy.get_path().stem()).string() + "_rho_RI_Unrestrained.cube");
    //cubes[2].set_path((dummy.get_path().parent_path() / dummy.get_path().stem()).string() + "_rho_RI_Enhanced.cube");
    //cubes[3].set_path((dummy.get_path().parent_path() / dummy.get_path().stem()).string() + "_rho_RI_Hybrid.cube");
    //std::vector<vec> all_coeffs = { coeff_unrestrained , coeff_enhanced, coeff_hybrid };
    //for (int a = 0; a < 4; a++) {
    //    //See if the cubefile already exists, if yes skip the calculation and just load it
    //    if (std::filesystem::exists(cubes[a].get_path())) {
    //        cubes[a].read_file(true, true);
    //    }
    //    else {
    //        std::cout << "Calculating cube for density..." << std::endl;
    //        if (a == 0)
    //            wavy.calc_rho_cube(cubes[a]);
    //        else
    //            calc_cube_ML(all_coeffs[a-1], dummy_aux, cubes[a]);
    //        cubes[a].write_file(false, false);
    //    }
    //}
    ////Calculate difference cubes
    //for (int a = 1; a < 4; a++) {
    //    cube diff = cubes[0] - cubes[a];
    //    diff.calc_dv();
    //    double rrs = cubes[0].rrs(cubes[a]);
    //    vec sums = diff.double_sum();
    //    std::cout << "Difference WFN - Method " << a - 1 << ": RRS = " << std::scientific << std::setprecision(6) << rrs
    //        << ", Sum = " << std::fixed << std::setprecision(6) << sums[0]
    //        << ", |Sum| = " << std::fixed << std::setprecision(6) << sums[1] << std::endl;
    //    diff.set_comment1("Difference cube: WFN - Method " + std::to_string(a - 1));
    //    diff.give_parent_wfn(dummy);
    //    diff.set_path((dummy.get_path().parent_path() / dummy.get_path().stem()).string() + "_rho_diff_Method" + std::to_string(a - 1) + ".cube");
    //    diff.write_file(true, false);
    //}


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

    size_t n_coeff = std::min({ coeff_unrestrained.size(), coeff_enhanced.size(), coeff_hybrid.size() });

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

