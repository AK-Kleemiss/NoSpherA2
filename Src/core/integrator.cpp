#include "pch.h"
#include "integrator.h"
#include "libCintMain.h"
#undef I // I had to include complex.h, so I have to undefine I for it to work in Linux.
#include "nos_math.h"
#include "GridManager.h"
#include "basis_set.h"
#include "SALTED_utilities.h"
#include <occ/disp/dftd4.h>
#include <occ/interaction/polarization.h>


vec einsum_ijk_ij_p(const dMatrix3& v1, const dMatrix2& v2)
{
    const int I = (int)v1.extent(0);
    const int J = (int)v1.extent(1);
    const int P = (int)v1.extent(2);
    // Initialize the result vector
    vec rho(P, 0.0);

    // Perform the summation
    for (int p = 0; p < P; p++)
    {
        for (int i = 0; i < I; i++)
        {
            for (int j = 0; j < J; j++)
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

vec DensityFitting::density_fit(const WFN& wavy, const WFN& wavy_aux, const CONFIG& config) {

    //if (wavy.get_origin() == e_origin::xtb || wavy.get_origin() == e_origin::ptb) {
    //    dMatrix2 cart2sph_matrix = get_cart2sph_matrix(wavy, true);
    //    dMatrix2 temp = dot(wavy.get_dm(), cart2sph_matrix, false, false);
    //    temp = dot(cart2sph_matrix, temp, true, false);
    //    wavy.set_dm(temp);
    //}

    vec eri2c;
    vec rho;

    // Initialize basis functions
    Int_Params normal_basis(wavy);
    Int_Params aux_basis(wavy_aux);
    dMatrix2 dm = wavy.get_dm();

    //Write out dm
    //for (size_t i = 0; i < dm.extent(0); i++) {
    //    std::cout << "[";
    //    for (size_t j = 0; j < dm.extent(1); j++) {
    //        std::cout << std::fixed << std::setprecision(6) << dm(i, j) << ", ";
    //    }
    //    std::cout <<"]," << std::endl;
    //}

    std::cout << "\n=== Density Fitting ===" << std::endl;
    std::cout << "Normal basis functions: " << normal_basis.get_nao() << std::endl;
    std::cout << "Auxiliary basis functions: " << aux_basis.get_nao() << std::endl;
    std::cout << "Metric: " << (config.metric == METRIC_TYPE::COULOMB ? "Coulomb" : "Overlap") << std::endl;


    // Compute integrals
    switch (config.metric) {
    case METRIC_TYPE::COULOMB:
        compute2C<Coulomb2C_SPH>(aux_basis, eri2c);
        computeRho<Coulomb3C_SPH>(wavy, wavy_aux, dm, rho, config.asym_atm_list);
        break;
    case METRIC_TYPE::OVERLAP:
        compute2C<Overlap2C_SPH>(aux_basis, eri2c);
        computeRho<Overlap3C_SPH>(normal_basis, aux_basis, dm, rho, config.asym_atm_list);
        break;
    }

    vec expected_populations;
    if (config.restrain_type == RESTRAINT_TYPE::NONE) {
        std::cout << "Solving unrestrained linear system..." << std::endl;
        solve_linear_system(eri2c, aux_basis.get_nao(), rho);
    }
    else {
        //Convert the charge scheme enum to string for display
        std::string charge_scheme;
        switch (config.charge_scheme) {
        case CHARGE_SCHEME::MULLIKEN:
            charge_scheme = "Mulliken";
            break;
        case CHARGE_SCHEME::SANDERSON_ESTIMATE:
            charge_scheme = "Sanderson Estimate";
            break;
        case CHARGE_SCHEME::TFVC:
            charge_scheme = "TFVC";
            break;
        case CHARGE_SCHEME::HIRSHFELD:
            charge_scheme = "Hirshfeld";
            break;
        case CHARGE_SCHEME::NUCLEAR:
            charge_scheme = "Nuclear Charge";
            break;
        case CHARGE_SCHEME::MBIS:
            charge_scheme = "MBIS";
            break;
        case CHARGE_SCHEME::EMBIS:
            charge_scheme = "EMBIS";
            break;
        }

        //Simple does nothing extra here
        if (config.restrain_type == RESTRAINT_TYPE::SIMPLE) {
            std::cout << "Using simple restraints with charges from: " << charge_scheme << std::endl;
        }
        //Hybrid applies Tikhonov regularization
        else if (config.restrain_type == RESTRAINT_TYPE::SIMPLE_AND_TIK) {
            std::cout << "Using simple restraints and Tikhonov regularization with charges from: " << charge_scheme << std::endl;
            // First apply Tikhonov regularization to the original matrix
            size_t n = aux_basis.get_nao();
            for (size_t i = 0; i < n; i++) {
                eri2c[i * n + i] += config.tikhonov_lambda;
            }
        }

        // Calculate expected charges based on chosen scheme
        vec2 targets;
        if (config.multipole_lmax >= 0) {
            //One grid for populations and moments; its l=0 term is the electron count without ECP electrons
            targets = calculate_expected_multipoles(wavy, config.charge_scheme, config.multipole_lmax);
            expected_populations.resize(wavy.get_ncen());
            for (int i = 0; i < wavy.get_ncen(); i++)
                expected_populations[i] = std::sqrt(constants::FOUR_PI) * targets[i][0];
        }
        else {
            expected_populations = calculate_expected_populations(wavy, wavy_aux, config.charge_scheme);
            if (wavy.get_has_ECPs()) {
                //Subtract the ecp electrons from the populations
                for (int i = 0; i < wavy.get_ncen(); i++) {
                    expected_populations[i] -= wavy.get_atom_ECP_electrons(i);
                }
            }
        }
        //The moment rows are meant to pin the partitioned atoms, so they carry the user's strength alone, O(1) against the metric
        const vec weights = config.multipole_lmax >= 0 ? vec(wavy.get_ncen(), config.multipole_strength)
            : restraint_weights(wavy_aux, aux_basis.get_nao(), config.restraint_strength, config.adaptive_restraint);
        // Apply enhanced electron restraints
        add_electron_restraint(eri2c, rho, wavy_aux, weights, expected_populations);
        if (config.multipole_lmax > 0)
            add_multipole_restraint(eri2c, rho, wavy_aux, targets, weights, config.multipole_lmax);
        // Solve the regularized 
        std::cout << "Solving regularized linear system..." << std::endl;
        if (config.multipole_lmax >= 0) {
            //Penalty on the Coulomb-metric objective, (J + R^T R) c = rho + R^T t: a deviation from the targets is paid
            //for in residual self-energy dc^T J dc. Least squares on the stacked [J; R] measures it as |J dc|^2 instead,
            //which costs nothing along the near-dependent diffuse functions and blew |c| up to ~100 on three heavy atoms.
            const int n = (int)aux_basis.get_nao(), m = (int)rho.size();
            vec J(eri2c.begin(), eri2c.begin() + (size_t)n * n), b(rho.begin(), rho.begin() + n);
            for (int k = n; k < m; k++) {
                const double* r = eri2c.data() + (size_t)k * n;
                for (int i = 0; i < n; i++) {
                    if (r[i] == 0.0) continue;
                    b[i] += r[i] * rho[k];
                    for (int j = 0; j < n; j++) J[(size_t)i * n + j] += r[i] * r[j];
                }
            }
            solve_linear_system(J, (size_t)n, b);
            double error = 0.0;
            for (int k = n; k < m; k++) {
                double d = -rho[k];
                for (int i = 0; i < n; i++) d += eri2c[(size_t)k * n + i] * b[i];
                error += d * d;
            }
            std::cout << "Restraint residual: " << std::fixed << std::setprecision(12) << std::sqrt(error) << std::endl;
            rho = b;
        }
        else {
            solve_linear_system(eri2c, rho.size(), aux_basis.get_nao(), rho);
            //dgels leaves the restraint residuals behind the solution
            rho.resize(aux_basis.get_nao());
        }
        if (config.multipole_lmax >= 0) {
            const vec2 fitted = fitted_multipoles(rho, wavy_aux, config.multipole_lmax);
            const double stone = std::sqrt(constants::FOUR_PI);
            std::cout << "\nAtomic multipoles of the fitted density, Racah normalisation, e bohr^l\n"
                << "  Atom  l  m       target       fitted    deviation" << std::endl;
            for (int a = 0; a < wavy_aux.get_ncen(); a++)
                for (int l = 0; l <= config.multipole_lmax; l++)
                    for (int m = -l; m <= l; m++) {
                        const double f = stone / std::sqrt(2.0 * l + 1.0), t = targets[a][l * l + l + m] * f, q = fitted[a][l * l + l + m] * f;
                        std::cout << std::setw(6) << wavy_aux.get_atom_label(a) << std::setw(3) << l << std::setw(3) << m
                            << std::fixed << std::setprecision(6) << std::setw(13) << t << std::setw(13) << q << std::setw(13) << q - t << std::endl;
                    }
        }
    }


    //write out fitted density coefficients
    //std::cout << "\nFitted density coefficients:" << std::endl;
    //for (size_t i = 0; i < rho.size(); i++) {
    //    std::cout << std::fixed << std::setprecision(6) << rho[i] << ", ";
    //}
    //std::cout << "\n";
    //rho = reorder_p(rho, wavy_aux);



    // Analyze quality if requested
    if (config.analyze_quality) {
        analyze_density_fit_quality(rho, wavy_aux, expected_populations);
    }

    std::cout << "==============================================\n" << std::endl;
    return rho;
}

DensityFitting::CONFIG DensityFitting::config_from_options(const options& opt)
{
    CONFIG config;
    config.analyze_quality = opt.debug;
    if (opt.multipole_lmax < 0)
        return config;
    config.restrain_type = RESTRAINT_TYPE::SIMPLE;
    config.multipole_lmax = opt.multipole_lmax;
    config.multipole_strength = opt.multipole_strength;
    switch (opt.multipole_scheme) {
    case PartitionType::TFVC: config.charge_scheme = CHARGE_SCHEME::TFVC; break;
    case PartitionType::MBIS: config.charge_scheme = CHARGE_SCHEME::MBIS; break;
    case PartitionType::EMBIS: config.charge_scheme = CHARGE_SCHEME::EMBIS; break;
    default: config.charge_scheme = CHARGE_SCHEME::HIRSHFELD; break;
    }
    return config;
}


// Enhanced electron restraint with adaptive weighting for s-orbitals only
// Only s-orbitals are restrained as other orbitals don't contribute to 
// spherically averaged electron density used in electron counting
// Base strength scaled by aux basis size (larger basis needs less), atom count (more atoms need more) and Z (heavier atoms more)
vec DensityFitting::restraint_weights(const WFN& wavy_aux, const size_t n_aux, double base_restraint_coef, bool adaptive_weighting)
{
    const size_t n_atoms = wavy_aux.get_ncen();
    double restraint_coef = base_restraint_coef;
    if (adaptive_weighting) {
        restraint_coef *= std::max(0.1, 1.0 - std::log10(n_aux) * 0.1);
        restraint_coef *= std::min(2.0, 1.0 + std::sqrt(n_atoms) * 0.1);
    }
    std::cout << "Setting adaptive restraint coefficient to: " << std::fixed
        << std::showpoint << std::setprecision(6) << restraint_coef << std::endl;
    vec weights(n_atoms, restraint_coef);
    if (adaptive_weighting)
        for (int a = 0; a < n_atoms; a++)
            weights[a] *= std::min(2.0, 1.0 + wavy_aux.get_atom_charge(a) * 0.02);
    return weights;
}

void DensityFitting::add_electron_restraint(vec& eri2c, vec& rho, const WFN& wavy_aux,
    const vec& atom_weights,
    const vec& expected_charges)
{
    double radial;
    basis_set_entry bf;

    size_t original_size = rho.size();
    size_t n_atoms = wavy_aux.get_ncen();

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

        const double atom_weight = atom_weights[atm_idx];
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
            }
            else {
                // Skip non-s orbitals but still increment indices appropriately
                coef_idx += (2 * type + 1);
            }

            prim += current_atom.get_shellcount()[shell];
        }
    }

    std::cout << "Added electron restraints for " << n_atoms << " atoms." << std::endl;
}

double DensityFitting::radial_moment(const double exponent, const double coef, const int l)
{
    primitive p(0, l, exponent, coef);
    return p.normalization_constant() * p.get_coef() * std::tgamma(l + 1.5) / (2.0 * std::pow(exponent, l + 1.5));
}

// One row per atom and (l, m) for 1 <= l <= lmax, non-zero only on the l-shells of that atom, since
// Q_lm of an atom-centred function about its own centre vanishes for every other l.
// Rows are scaled by 1 / r_cov^l so every order enters with the magnitude of the population row.
void DensityFitting::add_multipole_restraint(vec& eri2c, vec& rho, const WFN& wavy_aux,
    const vec2& targets, const vec& atom_weights, const int lmax)
{
    err_checkf(lmax >= 1 && lmax <= 8, "Multipole restraints are implemented for 1 <= l <= 8", std::cout);
    const int n_atoms = wavy_aux.get_ncen(), n_old = (int)rho.size(), n_aux = (int)(eri2c.size() / rho.size());
    const int per_atom = (lmax + 1) * (lmax + 1) - 1, n_rows = n_atoms * per_atom;
    eri2c.resize((size_t)(n_old + n_rows) * n_aux, 0.0);
    rho.resize(n_old + n_rows, 0.0);
    dMatrixRef2 rows(eri2c.data() + (size_t)n_old * n_aux, n_rows, n_aux);
    ivec skipped(lmax + 1, 0);
    int coef_idx = 0;
    for (int a = 0; a < n_atoms; a++) {
        const atom A = wavy_aux.get_atom(a);
        const double r_cov = constants::ang2bohr(constants::covalent_radii[A.get_charge()]);
        const int row0 = a * per_atom;
        bvec has_l(lmax + 1, false);
        int prim = 0;
        for (int shell = 0; shell < (int)A.get_shellcount_size(); shell++) {
            const int l = A.get_basis_set_entry(prim).get_type();
            if (l >= 1 && l <= lmax) {
                double I_l = 0.0;
                for (int e = 0; e < (int)A.get_shellcount(shell); e++) {
                    const basis_set_entry& bf = A.get_basis_set_entry(prim + e);
                    I_l += radial_moment(bf.get_exponent(), bf.get_coefficient(), l);
                }
                const double w = atom_weights[a] / std::pow(r_cov, l);
                for (int m = -l; m <= l; m++)
                    rows(row0 + l * l - 1 + l + m, coef_idx + l + m) += I_l * w;
                has_l[l] = true;
            }
            coef_idx += 2 * l + 1;
            prim += A.get_shellcount(shell);
        }
        for (int l = 1; l <= lmax; l++) {
            if (!has_l[l]) {
                skipped[l]++;
                continue;
            }
            const double w = atom_weights[a] / std::pow(r_cov, l);
            for (int m = -l; m <= l; m++)
                rho[n_old + row0 + l * l - 1 + l + m] = targets[a][l * l + l + m] * w;
        }
    }
    for (int l = 1; l <= lmax; l++)
        if (skipped[l] > 0)
            std::cout << skipped[l] << " atoms carry no l=" << l << " auxiliary functions, their order " << l << " moments are not restrained." << std::endl;
    std::cout << "Added multipole restraints up to l=" << lmax << " for " << n_atoms << " atoms." << std::endl;
}

vec2 DensityFitting::fitted_multipoles(const vec& coefficients, const WFN& wavy_aux, const int lmax)
{
    const int n_atoms = wavy_aux.get_ncen();
    vec2 moments(n_atoms, vec((lmax + 1) * (lmax + 1), 0.0));
    int coef_idx = 0;
    for (int a = 0; a < n_atoms; a++) {
        const atom A = wavy_aux.get_atom(a);
        int prim = 0;
        for (int shell = 0; shell < (int)A.get_shellcount_size(); shell++) {
            const int l = A.get_basis_set_entry(prim).get_type();
            if (l <= lmax) {
                double I_l = 0.0;
                for (int e = 0; e < (int)A.get_shellcount(shell); e++) {
                    const basis_set_entry& bf = A.get_basis_set_entry(prim + e);
                    I_l += radial_moment(bf.get_exponent(), bf.get_coefficient(), l);
                }
                for (int m = -l; m <= l; m++)
                    moments[a][l * l + l + m] += I_l * coefficients[coef_idx + l + m];
            }
            coef_idx += 2 * l + 1;
            prim += A.get_shellcount(shell);
        }
    }
    return moments;
}

double DensityFitting::lower_gamma_half(const int l, const double x)
{
    const double a = l + 1.5;
    if (x < a + 1.0) {
        double term = 1.0 / a, sum = term;
        for (int k = 1; k < 500 && term > 1e-17 * sum; k++) {
            term *= x / (a + k);
            sum += term;
        }
        return std::pow(x, a) * std::exp(-x) * sum;
    }
    double g = std::sqrt(constants::PI) * std::erf(std::sqrt(x)), xa = std::sqrt(x);
    const double ex = std::exp(-x);
    for (int k = 0; k <= l; k++) {
        g = (k + 0.5) * g - xa * ex;
        xa *= x;
    }
    return g;
}

// V(R) = 4pi/(2l+1) N c [R^-l-1 gamma(l+3/2, aR^2)/(2a^(l+3/2)) + R^l exp(-aR^2)/(2a)] Y_lm(R^)
double DensityFitting::aux_potential(const double exponent, const double coef, const int l, const int m, const double* R)
{
    const primitive p(0, l, exponent, coef);
    const double nc = p.normalization_constant() * p.get_coef(), r2 = R[0] * R[0] + R[1] * R[1] + R[2] * R[2], r = std::sqrt(r2);
    if (r < 1e-10) return l == 0 ? constants::FOUR_PI * nc * constants::c_1_4p / (2.0 * exponent) : 0.0;
    const double d[3] = { R[0] / r, R[1] / r, R[2] / r };
    const double rad = std::pow(r, -l - 1) * lower_gamma_half(l, exponent * r2) / (2.0 * std::pow(exponent, l + 1.5)) + std::pow(r, l) * std::exp(-exponent * r2) / (2.0 * exponent);
    return constants::FOUR_PI / (2 * l + 1) * nc * rad * constants::spherical_harmonic(l, m, d);
}

namespace {
    // Atom, rank, m and primitive range of every aux function in coefficient order
    struct aux_index { ivec atom, l, m, prim, nprim; int lmax = 0; };
    aux_index index_aux(const WFN& aux)
    {
        aux_index ix;
        for (int a = 0; a < aux.get_ncen(); a++) {
            const atom A = aux.get_atom(a);
            int prim = 0;
            for (int shell = 0; shell < (int)A.get_shellcount_size(); shell++) {
                const int l = A.get_basis_set_entry(prim).get_type(), n = A.get_shellcount(shell);
                for (int m = -l; m <= l; m++) {
                    ix.atom.push_back(a);
                    ix.l.push_back(l);
                    ix.m.push_back(m);
                    ix.prim.push_back(prim);
                    ix.nprim.push_back(n);
                }
                ix.lmax = std::max(ix.lmax, l);
                prim += n;
            }
        }
        return ix;
    }
    // Potential of aux function i at R, electrons counted positive
    double function_potential(const WFN& aux, const aux_index& ix, const int i, const double* R)
    {
        const atom A = aux.get_atom(ix.atom[i]);
        const double d[3] = { R[0] - A.get_coordinate(0), R[1] - A.get_coordinate(1), R[2] - A.get_coordinate(2) };
        double v = 0.0;
        for (int e = 0; e < ix.nprim[i]; e++) {
            const basis_set_entry& bf = A.get_basis_set_entry(ix.prim[i] + e);
            v += DensityFitting::aux_potential(bf.get_exponent(), bf.get_coefficient(), ix.l[i], ix.m[i], d);
        }
        return v;
    }
    // Field -grad phi at R of the nuclei and fitted density of one molecule, the density part by central differences of the potential
    void molecule_field(const WFN& aux, const aux_index& ix, const vec& coef, const double* R, double* F)
    {
        const double h = 1e-4;
        F[0] = F[1] = F[2] = 0.0;
        for (int a = 0; a < aux.get_ncen(); a++) {
            double d[3], r2 = 0.0;
            for (int x = 0; x < 3; x++) { d[x] = R[x] - aux.get_atom_coordinate(a, x); r2 += d[x] * d[x]; }
            const double f = (aux.get_atom_charge(a) - aux.get_atom_ECP_electrons(a)) / (r2 * std::sqrt(r2));
            for (int x = 0; x < 3; x++) F[x] += f * d[x];
        }
        for (int i = 0; i < (int)coef.size(); i++)
            for (int x = 0; x < 3; x++) {
                double Rp[3] = { R[0], R[1], R[2] }, Rm[3] = { R[0], R[1], R[2] };
                Rp[x] += h, Rm[x] -= h;
                F[x] += coef[i] * (function_potential(aux, ix, i, Rp) - function_potential(aux, ix, i, Rm)) / (2 * h);
            }
    }
    // -1/2 sum alpha_a |F_a|^2 over the atoms of aux in the field of partner, Thakkar polarizabilities as in CrystalExplorer
    double polarization(const WFN& aux, const WFN& partner, const aux_index& ixP, const vec& coef_P)
    {
        occ::IVec Z(aux.get_ncen());
        occ::Mat3N F(3, aux.get_ncen());
        for (int a = 0; a < aux.get_ncen(); a++) {
            const double R[3] = { aux.get_atom_coordinate(a, 0), aux.get_atom_coordinate(a, 1), aux.get_atom_coordinate(a, 2) };
            double f[3];
            molecule_field(partner, ixP, coef_P, R, f);
            Z(a) = aux.get_atom_charge(a);
            for (int x = 0; x < 3; x++) F(x, a) = f[x];
        }
        return occ::interaction::ce_model_polarization_energy(Z, F, aux.get_charge() != 0);
    }
    std::vector<occ::core::Atom> occ_atoms(const WFN& aux)
    {
        std::vector<occ::core::Atom> atoms(aux.get_ncen());
        for (int a = 0; a < aux.get_ncen(); a++)
            atoms[a] = { aux.get_atom_charge(a), aux.get_atom_coordinate(a, 0), aux.get_atom_coordinate(a, 1), aux.get_atom_coordinate(a, 2) };
        return atoms;
    }
    double d4_energy(const std::vector<occ::core::Atom>& atoms, const int charge)
    {
        occ::disp::D4Dispersion d4(atoms);
        d4.set_charge(charge);
        return d4.energy();
    }
    // Gordon-Kim exchange-repulsion of the fitted densities on a Becke grid over the dimer: gk[0] the Thomas-Fermi kinetic
    // energy, gk[1] 1/9 of the von Weizsaecker gradient correction and gk[2] the Dirac exchange of rhoA + rhoB minus the
    // monomers, gk[3], gk[4] the electron counts of A and B on the grid. Gradients by central differences of the aux density
    void gordon_kim(const vec& coef_A, const WFN& aux_A, const vec& coef_B, const WFN& aux_B, const int x_fun, double* gk)
    {
        const int wfn_type[6] = { 1, 2, 5, 11, 21, 36 };
        const WFN* mono[2] = { &aux_A, &aux_B };
        WFN dimer(e_origin::NOT_YET_DEFINED);
        for (int m = 0; m < 2; m++)
            for (int a = 0; a < mono[m]->get_ncen(); a++) {
                const atom& A = mono[m]->get_atom(a);
                dimer.push_back_atom(A);
                for (int b = 0; b < (int)A.get_basis_set_size(); b++)
                    dimer.add_exp(dimer.get_ncen(), wfn_type[std::min(A.get_basis_set_type(b), 5)], A.get_basis_set_exponent(b));
            }
        GridConfiguration config;
        config.partition_type = PartitionType::Becke;
        config.no_density_eval = true;
        GridManager gm(config);
        ivec all(dimer.get_ncen());
        for (int a = 0; a < dimer.get_ncen(); a++) all[a] = a;
        std::ostringstream quiet;
        gm.setup3DGridsForMolecule(dimer, all, bvec(dimer.get_ncen(), true), cell(), false, quiet);
        const std::vector<atom> at_A = aux_A.get_atoms(), at_B = aux_B.get_atoms();
        const GridData& GD = gm.getGridData();
        const double C_TF = 0.3 * std::pow(3 * constants::PI * constants::PI, 2.0 / 3.0), h = 1e-4;
        double tf = 0.0, vw = 0.0, x = 0.0, nA = 0.0, nB = 0.0;
        for (int a = 0; a < dimer.get_ncen(); a++) {
            const vec2& g = GD.atomic_grids[a];
            const int n = GD.num_points_per_atom[a];
#pragma omp parallel for reduction(+:tf, vw, x, nA, nB)
            for (int p = 0; p < n; p++) {
                const double w = g[GridData::GridIndex::BECKE_WEIGHT][p], r[3] = { g[0][p], g[1][p], g[2][p] };
                double rho[3] = { calc_density_ML(r[0], r[1], r[2], coef_A, at_A), calc_density_ML(r[0], r[1], r[2], coef_B, at_B), 0.0 }, g2[3] = { 0.0, 0.0, 0.0 };
                rho[2] = rho[0] + rho[1];
                nA += w * rho[0], nB += w * rho[1];
                for (int c = 0; c < 3; c++) {
                    double rp[3] = { r[0], r[1], r[2] }, rm[3] = { r[0], r[1], r[2] };
                    rp[c] += h, rm[c] -= h;
                    const double dA = (calc_density_ML(rp[0], rp[1], rp[2], coef_A, at_A) - calc_density_ML(rm[0], rm[1], rm[2], coef_A, at_A)) / (2 * h);
                    const double dB = (calc_density_ML(rp[0], rp[1], rp[2], coef_B, at_B) - calc_density_ML(rm[0], rm[1], rm[2], coef_B, at_B)) / (2 * h);
                    g2[0] += dA * dA, g2[1] += dB * dB, g2[2] += (dA + dB) * (dA + dB);
                }
                for (int m = 0; m < 3; m++) {
                    if (rho[m] < 1e-12) continue;
                    const double s = m == 2 ? 1.0 : -1.0;
                    tf += s * w * C_TF * std::pow(rho[m], 5.0 / 3.0);
                    vw += s * w * g2[m] / (72.0 * rho[m]);
                    x += s * w * DensityFitting::exchange_density(rho[m], g2[m], x_fun);
                }
            }
        }
        gk[0] = tf, gk[1] = vw, gk[2] = x, gk[3] = nA, gk[4] = nB;
    }
}

double DensityFitting::exchange_density(const double rho, const double g2, const int x_fun)
{
    const double C_X = -0.75 * std::pow(3 * constants::INV_PI, 1.0 / 3.0), r43 = std::pow(rho, 4.0 / 3.0);
    if (x_fun == 1) {
        const double kappa = 0.804, mu = 0.2195149727645171;
        const double s2 = g2 / (4 * std::pow(3 * constants::PI * constants::PI, 2.0 / 3.0) * r43 * r43);
        return C_X * r43 * (1 + kappa - kappa / (1 + mu * s2 / kappa));
    }
    if (x_fun == 2) {
        const double beta = 0.0042, x = std::cbrt(2.0) * std::sqrt(g2) / r43;
        return C_X * r43 - 2 * beta * std::pow(0.5 * rho, 4.0 / 3.0) * x * x / (1 + 6 * beta * x * std::asinh(x));
    }
    return C_X * r43;
}

DensityFitting::INTERACTION DensityFitting::interaction_energy(const vec& coef_A, const WFN& aux_A, const vec& coef_B, const WFN& aux_B, const double repulsion_K, const int x_fun)
{
    const aux_index ixA = index_aux(aux_A), ixB = index_aux(aux_B);
    err_checkf(coef_A.size() == ixA.atom.size() && coef_B.size() == ixB.atom.size(), "Coefficient count does not match the auxiliary basis", std::cout);
    const int nA = aux_A.get_ncen(), nB = aux_B.get_ncen(), LA = ixA.lmax, LB = ixB.lmax;
    INTERACTION E;
    E.pair.assign(nA, vec(nB, 0.0));
    E.rank.assign(LA + 2, vec(LB + 2, 0.0));
    ivec ZA(nA), ZB(nB);
    for (int a = 0; a < nA; a++) ZA[a] = aux_A.get_atom_charge(a) - aux_A.get_atom_ECP_electrons(a);
    for (int b = 0; b < nB; b++) ZB[b] = aux_B.get_atom_charge(b) - aux_B.get_atom_ECP_electrons(b);
    for (int a = 0; a < nA; a++)
        for (int b = 0; b < nB; b++) {
            double d2 = 0.0;
            for (int x = 0; x < 3; x++) d2 += std::pow(aux_A.get_atom_coordinate(a, x) - aux_B.get_atom_coordinate(b, x), 2);
            const double e = ZA[a] * ZB[b] / std::sqrt(d2);
            E.nuc_nuc += e;
            E.pair[a][b] += e;
        }
    E.rank[0][0] = E.nuc_nuc;
    for (int a = 0; a < nA; a++) {
        const double R[3] = { aux_A.get_atom_coordinate(a, 0), aux_A.get_atom_coordinate(a, 1), aux_A.get_atom_coordinate(a, 2) };
        for (int i = 0; i < (int)coef_B.size(); i++) {
            const double e = -ZA[a] * coef_B[i] * function_potential(aux_B, ixB, i, R);
            E.nucA_rhoB += e;
            E.pair[a][ixB.atom[i]] += e;
            E.rank[0][ixB.l[i] + 1] += e;
        }
    }
    for (int b = 0; b < nB; b++) {
        const double R[3] = { aux_B.get_atom_coordinate(b, 0), aux_B.get_atom_coordinate(b, 1), aux_B.get_atom_coordinate(b, 2) };
        for (int i = 0; i < (int)coef_A.size(); i++) {
            const double e = -ZB[b] * coef_A[i] * function_potential(aux_A, ixA, i, R);
            E.nucB_rhoA += e;
            E.pair[ixA.atom[i]][b] += e;
            E.rank[ixA.l[i] + 1][0] += e;
        }
    }
    Int_Params pA(aux_A), pB(aux_B);
    Int_Params pAB(pA, pB);
    const int na = pA.get_nao(), nao = pAB.get_nao();
    err_checkf(na == (int)coef_A.size() && nao - na == (int)coef_B.size(), "Two-centre integral dimension does not match the coefficients", std::cout);
    vec J;
    compute2C<Coulomb2C_SPH>(pAB, J);
    for (int i = 0; i < na; i++)
        for (int j = 0; j < nao - na; j++) {
            const double e = coef_A[i] * J[i * nao + na + j] * coef_B[j];
            E.rho_rho += e;
            E.pair[ixA.atom[i]][ixB.atom[j]] += e;
            E.rank[ixA.l[i] + 1][ixB.l[j] + 1] += e;
        }
    vec S;
    compute2C<Overlap2C_SPH>(pAB, S);
    for (int i = 0; i < na; i++)
        for (int j = 0; j < nao - na; j++) E.overlap += coef_A[i] * S[i * nao + na + j] * coef_B[j];
    if (repulsion_K > 0.0) E.rep = repulsion_K * E.overlap;
    else {
        double gk[5];
        gordon_kim(coef_A, aux_A, coef_B, aux_B, x_fun, gk);
        E.rep_kin = gk[0], E.rep_vw = gk[1], E.rep_x = gk[2], E.n_A = gk[3], E.n_B = gk[4], E.x_fun = x_fun;
        E.rep = E.rep_kin + E.rep_x;
    }
    E.pol_A = polarization(aux_A, aux_B, ixB, coef_B);
    E.pol_B = polarization(aux_B, aux_A, ixA, coef_A);
    std::vector<occ::core::Atom> atoms = occ_atoms(aux_A), atoms_B = occ_atoms(aux_B);
    E.disp = -d4_energy(atoms, aux_A.get_charge()) - d4_energy(atoms_B, aux_B.get_charge());
    atoms.insert(atoms.end(), atoms_B.begin(), atoms_B.end());
    E.disp += d4_energy(atoms, aux_A.get_charge() + aux_B.get_charge());
    return E;
}

void DensityFitting::print_interaction_energy(const INTERACTION& E, const WFN& aux_A, const WFN& aux_B, std::ostream& file)
{
    const bool KS = E.n_A == 0.0 && E.n_B == 0.0;
    const char* x_names[3] = { "rep. exch. Dirac ", "rep. exch. PBE   ", "rep. exch. B88   " };
    const char* names[13] = { "nucleus-nucleus  ", "nuclei A - rho B ", "nuclei B - rho A ", "rho A - rho B    ", "electrostatic    ",
                              "pol. A in field B", "pol. B in field A", "dispersion D4    ", "rep. kin. TF     ", x_names[E.x_fun], KS ? "repulsion K*S    " : "repulsion GK     ", "total            ", "vW/9 not in total" };
    const double parts[13] = { E.nuc_nuc, E.nucA_rhoB, E.nucB_rhoA, E.rho_rho, E.electrostatic(), E.pol_A, E.pol_B, E.disp, E.rep_kin, E.rep_x, E.rep, E.total(), E.rep_vw };
    file << "\nElectrostatic interaction energy of the fitted densities\n";
    for (int i = 0; i < (KS ? 12 : 13); i++) {
        if (i == 5) {
            file << "\nBeyond electrostatics: Thakkar polarizabilities in the partner's field, D4 with PBE damping, density overlap S = Int rhoA rhoB\n";
            file << "  S                " << std::scientific << std::setprecision(6) << std::setw(14) << E.overlap << " e^2/bohr^3" << std::fixed;
            if (KS) file << "  (repulsion K*S)\n";
            else file << "  (repulsion Gordon-Kim on a Becke grid holding " << std::setprecision(4) << E.n_A << " / " << E.n_B << " e)\n";
        }
        if (KS && (i == 8 || i == 9)) continue;
        file << "  " << names[i] << std::fixed << std::setprecision(6) << std::setw(14) << parts[i] << " Eh" << std::setprecision(4) << std::setw(12) << parts[i] * constants::kcal_mol_per_hartree << " kcal/mol\n";
    }
    file << "\nBy atom pair, kcal/mol, rows A columns B\n      ";
    for (int b = 0; b < aux_B.get_ncen(); b++) file << std::setw(9) << aux_B.get_atom_label(b);
    file << "      sum\n" << std::setprecision(3);
    for (int a = 0; a < aux_A.get_ncen(); a++) {
        double s = 0.0;
        file << std::setw(6) << aux_A.get_atom_label(a);
        for (int b = 0; b < aux_B.get_ncen(); b++) {
            file << std::setw(9) << E.pair[a][b] * constants::kcal_mol_per_hartree;
            s += E.pair[a][b];
        }
        file << std::setw(9) << s * constants::kcal_mol_per_hartree << "\n";
    }
    file << "\nBy rank, kcal/mol, rows A columns B, n = nuclei\n      ";
    for (int j = 0; j < (int)E.rank[0].size(); j++) file << std::setw(11) << (j == 0 ? std::string("n") : "l=" + std::to_string(j - 1));
    file << "      sum\n";
    for (int i = 0; i < (int)E.rank.size(); i++) {
        double s = 0.0;
        file << std::setw(6) << (i == 0 ? std::string("n") : "l=" + std::to_string(i - 1));
        for (int j = 0; j < (int)E.rank[i].size(); j++) {
            file << std::setw(11) << E.rank[i][j] * constants::kcal_mol_per_hartree;
            s += E.rank[i][j];
        }
        file << std::setw(11) << s * constants::kcal_mol_per_hartree << "\n";
    }
    file << std::endl;
}

static PartitionType scheme_partition(const DensityFitting::CHARGE_SCHEME& scheme)
{
    switch (scheme) {
    case DensityFitting::CHARGE_SCHEME::TFVC: return PartitionType::TFVC;
    case DensityFitting::CHARGE_SCHEME::HIRSHFELD: return PartitionType::Hirshfeld;
    case DensityFitting::CHARGE_SCHEME::MBIS: return PartitionType::MBIS;
    case DensityFitting::CHARGE_SCHEME::EMBIS: return PartitionType::EMBIS;
    default: err_not_impl_f("Only TFVC, Hirshfeld, MBIS and EMBIS partition the density on a grid", std::cout);
    }
    return PartitionType::Hirshfeld;
}

//PartitionType and PartitionResults::CHARGE_ORDER are numbered differently
static int charge_order(const PartitionType type)
{
    switch (type) {
    case PartitionType::TFVC: return PartitionResults::CHARGE_ORDER::S_TFVC;
    case PartitionType::MBIS: return PartitionResults::CHARGE_ORDER::S_MBIS;
    case PartitionType::EMBIS: return PartitionResults::CHARGE_ORDER::S_EMBIS;
    default: return PartitionResults::CHARGE_ORDER::S_HIRSH;
    }
}

vec2 DensityFitting::calculate_expected_multipoles(const WFN& wavy, const CHARGE_SCHEME& scheme, const int lmax)
{
    GridConfiguration config;
    config.partition_type = scheme_partition(scheme);
    config.pbc = 0;
    config.debug = false;
    const int ncen = wavy.get_ncen();
    ivec atom_list(ncen);
    for (int a = 0; a < ncen; a++)
        atom_list[a] = a;
    GridManager grid_manager(config);
    WFN temp = wavy;
    temp.delete_unoccupied_MOs();
    grid_manager.setup3DGridsForMolecule(temp, atom_list);
    return grid_manager.calculatePartitionedMultipoles(temp, lmax);
}

// Calculate expected atomic populations based on different partitioning schemes
vec DensityFitting::calculate_expected_populations(const WFN& wavy, const WFN& wavy_aux, const CHARGE_SCHEME& scheme)
{
    vec expected_populations(wavy_aux.get_ncen());

    if (scheme == CHARGE_SCHEME::NUCLEAR) {
        // Simple nuclear populations
        for (int i = 0; i < wavy_aux.get_ncen(); i++) {
            expected_populations[i] = wavy_aux.get_atoms()[i].get_charge();
        }
    }
    // https://pubs.acs.org/doi/10.1021/ed065p227
    else if (scheme == CHARGE_SCHEME::SANDERSON_ESTIMATE) {
        double compound_electronegativity = 1.0;
        for (const auto& atom : wavy.get_atoms()) {
            compound_electronegativity *= constants::allen_electronegativities[atom.get_charge() - 1];
        }
        compound_electronegativity = std::pow(compound_electronegativity, 1.0 / wavy.get_ncen());

        for (int iat = 0; iat < wavy_aux.get_ncen(); iat++) {
            double atom_electronegativity = constants::allen_electronegativities[wavy_aux.get_atoms()[iat].get_charge() - 1];
            expected_populations[iat] = wavy_aux.get_atoms()[iat].get_charge() + (compound_electronegativity - atom_electronegativity) / (1.57 * std::sqrt(atom_electronegativity));
        }
    }
    else if (scheme == CHARGE_SCHEME::MULLIKEN) {
        dMatrix2 dm = wavy.get_dm();
        vec eri2c;
        Int_Params normal_basis(wavy);
        compute2C<Overlap2C_SPH>(normal_basis, eri2c);
        dMatrixRef2 eri2c_ref(eri2c.data(), normal_basis.get_nao(), normal_basis.get_nao());
        const size_t nao = dm.extent(1);

        // optional: symmetric Mulliken operator M = 1/2 (P S + S P)
        // (avoid forming full matrices if memory is tight; just accumulate the diagonal)
        vec mulliken_pop(wavy.get_ncen(), 0.0);

        size_t mu_begin = 0;
        for (unsigned iat = 0; iat < wavy.get_ncen(); iat++) {
            const atom A = wavy.get_atoms()[iat];

            // determine how many AOs belong to this atom (use contracted shells, not primitives!)
            size_t nAO_A = 0;
            int prim = 0;
            for (size_t sh = 0; sh < A.get_shellcount().size(); sh++) {
                int l = A.get_basis_set_entry(prim).get_type() - 1;
                nAO_A += size_t(2 * l + 1);
                prim += A.get_shellcount()[sh];
            }
            size_t mu_end = mu_begin + nAO_A;

            double GA = 0.0;
            for (size_t m = mu_begin; m < mu_end; m++) {
                double diag_PS = 0.0;
                for (size_t n = 0; n < nao; n++)
                    diag_PS += dm(m, n) * eri2c_ref(n, m);
                GA += diag_PS;
            }
            expected_populations[iat] = GA;//A.get_charge() - GA;        // Mulliken charge
            mu_begin = mu_end;
        }
    }
    else if (scheme == CHARGE_SCHEME::TFVC || scheme == CHARGE_SCHEME::HIRSHFELD || scheme == CHARGE_SCHEME::MBIS || scheme == CHARGE_SCHEME::EMBIS) {
        PartitionType type = scheme_partition(scheme);
        GridConfiguration config;
        config.accuracy = 0;
        config.partition_type = type;
        config.pbc = 0;
        config.debug = false;
        const int ncen = wavy.get_ncen();

        ivec asym_atom_list(ncen);
        for (int atom_nr = 0; atom_nr < ncen; atom_nr++) {
            asym_atom_list[atom_nr] = atom_nr;
        }
        svec labels(ncen);
        const auto atoms = wavy.get_atoms();
        for (int i = 0; i < ncen; i++) {
            labels[i] = atoms[i].get_label();
        }

        GridManager grid_manager(config);

        WFN temp = wavy;
        temp.delete_unoccupied_MOs();
        // Setup grids for the molecule
        grid_manager.setup3DGridsForMolecule(temp, asym_atom_list);


        // Calculate partitioned charges
        auto results = grid_manager.calculatePartitionedCharges(temp);
        //results.printChargeTable(labels, temp, std::cout);
        for (int i = 0; i < temp.get_ncen(); i++) {
            expected_populations[i] = results.atom_charges[charge_order(type)][i];
        }
    }
    else {
        std::cerr << "Warning: Unknown charge scheme. Defaulting to nuclear charges.'" << std::endl;
        expected_populations = calculate_expected_populations(wavy, wavy_aux, CHARGE_SCHEME::NUCLEAR);
    }

    return expected_populations;
}

// Analyze the quality of density fitting and detect problematic charges
void DensityFitting::analyze_density_fit_quality(const vec& coefficients, const WFN& wavy_aux,
    const vec& expected_charges)
{
    std::cout << "\n=== Density Fitting Quality Analysis ===" << std::endl;

    vec atomic_populations(wavy_aux.get_ncen(), 0.0);
    int expected_total_electrons = 0;
    int coef_idx = 0;

    // Calculate atomic populations from coefficients using only s-orbitals
    for (int atm_idx = 0; atm_idx < wavy_aux.get_ncen(); atm_idx++) {
        atom current_atom = wavy_aux.get_atoms()[atm_idx];
        expected_total_electrons += current_atom.get_charge();
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
            }
            else {
                // Skip non-s orbitals but still increment coefficient index
                coef_idx += (2 * type + 1);
            }

            prim += current_atom.get_shellcount()[shell];
        }

        atomic_populations[atm_idx] += current_atom.get_ECP_electrons(); // Include ECP electrons if any
        double n_electrons = current_atom.get_charge();
        double computed_charge = n_electrons - atomic_populations[atm_idx];
        double expected_charge = 0.0;

        if (!expected_charges.empty() && atm_idx < expected_charges.size()) {
            expected_charge = n_electrons - (expected_charges[atm_idx] + current_atom.get_ECP_electrons());
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
    double real_total_electrons = std::accumulate(atomic_populations.begin(), atomic_populations.end(), 0.0);
    std::cout << "Expected / Real total electrons: " << expected_total_electrons << " / " << std::fixed << std::setprecision(3) << real_total_electrons << std::endl;
}

#include "SALTED_utilities.h"
//#include "test_functions.h"
// Example usage function demonstrating the enhanced density fitting approaches
void DensityFitting::demonstrate_enhanced_density_fitting(WFN& wavy, const WFN& wavy_aux)
{
    //#include "test_functions.h"
    std::cout << "\n=== Enhanced Density Fitting Demonstration ===" << std::endl;
    CONFIG ri_config;
    ri_config.adaptive_restraint = true;
    ri_config.analyze_quality = true;
    ri_config.charge_scheme = CHARGE_SCHEME::HIRSHFELD;
    ri_config.metric = METRIC_TYPE::COULOMB;
    ri_config.restraint_strength = 2e-4;
    ri_config.tikhonov_lambda = 1e-6;

    _time_point start_time = get_time();
    // Method 0: Unrestrained (original approach)
    std::cout << "\n--- Method 0: Unrestrained (Baseline) ---" << std::endl;
    ri_config.restrain_type = RESTRAINT_TYPE::NONE;
    vec coeff_unrestrained = density_fit(wavy, wavy_aux, ri_config);
    std::cout << "Time for unrestrained fit: "
        << std::chrono::duration<double>(get_time() - start_time).count() << " seconds." << std::endl;

    start_time = get_time();
    // Method 1: Enhanced restraints with adaptive weighting
    std::cout << "\n--- Method 1: Enhanced Adaptive Restraints ---" << std::endl;
    ri_config.restrain_type = RESTRAINT_TYPE::SIMPLE;
    vec coeff_enhanced = density_fit(wavy, wavy_aux, ri_config);     // analyze quality
    std::cout << "Time for enhanced restraint fit: "
        << std::chrono::duration<double>(get_time() - start_time).count() << " seconds." << std::endl;

    start_time = get_time();
    // Method 2: Hybrid approach
    std::cout << "\n--- Method 2: Hybrid Regularization ---" << std::endl;
    ri_config.restrain_type = RESTRAINT_TYPE::SIMPLE_AND_TIK;
    vec coeff_hybrid = density_fit(wavy, wavy_aux, ri_config); // charge scheme
    std::cout << "Time for hybrid fit: "
        << std::chrono::duration<double>(get_time() - start_time).count() << " seconds." << std::endl;

    GridConfiguration config;
    config.accuracy = 2;
    config.partition_type = PartitionType::Hirshfeld;
    config.pbc = 0;
    config.debug = false;
    config.all_charges = true;
    const int weight_index = GridData::HIRSH_WEIGHT;

    ivec asym_atom_list(wavy.get_ncen());
    for (int atom_nr = 0; atom_nr < wavy.get_ncen(); atom_nr++) {
        asym_atom_list[atom_nr] = atom_nr;
    }

    GridManager grid_manager(config);
    WFN temp = wavy;

    temp.delete_unoccupied_MOs();
    // Setup grids for the molecule
    grid_manager.setup3DGridsForMolecule(temp, asym_atom_list);
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
            + ", Sum = " + std::to_string(diff_densities[i][DiffDensityIndex::DIFF_UNRESTRAINED][SumIndex::SUM_NO_DIFF])
            + ", Partitioned Sum = " + std::to_string(partitioned_densities[i]);
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


void DensityFitting::QM_RI_difference_cube(WFN& wavy, const WFN& wavy_aux) {
    //#include "test_functions.h"
    std::cout << "\n=== Enhanced Density Fitting Demonstration ===" << std::endl;
    CONFIG ri_config;
    ri_config.adaptive_restraint = true;
    ri_config.analyze_quality = true;
    ri_config.charge_scheme = CHARGE_SCHEME::HIRSHFELD;
    ri_config.metric = METRIC_TYPE::COULOMB;
    ri_config.restraint_strength = 2e-4;
    ri_config.tikhonov_lambda = 1e-6;

    _time_point start_time = get_time();
    // Method 0: Unrestrained (original approach)
    std::cout << "\n--- Method 0: Unrestrained (Baseline) ---" << std::endl;
    ri_config.restrain_type = RESTRAINT_TYPE::NONE;
    vec coeffs = density_fit(wavy, wavy_aux, ri_config);
    std::cout << "Time for unrestrained fit: "
        << std::chrono::duration<double>(get_time() - start_time).count() << " seconds." << std::endl;
    std::cout << "Time for RI-Fitting: "
        << std::chrono::duration<double>(get_time() - start_time).count() << " seconds." << std::endl;

    WFN dummy = wavy;
    WFN dummy_aux = wavy_aux;

    properties_options props;
    props.radius = 3.0;
    props.resolution = 0.1;

    readxyzMinMax_fromWFN(dummy, props, true);
    dummy.delete_unoccupied_MOs();
    //Cubes: 0=WFN, 1=RI
    cube WFN_cube({ props.NbSteps[0], props.NbSteps[1], props.NbSteps[2] }, dummy.get_ncen(), true);
    cube RI_cube({ props.NbSteps[0], props.NbSteps[1], props.NbSteps[2] }, dummy.get_ncen(), true);
    WFN_cube.give_parent_wfn(dummy);
    RI_cube.give_parent_wfn(dummy_aux);

    for (int i = 0; i < 3; i++) {
        WFN_cube.set_origin(i, props.MinMax[i]);
        WFN_cube.set_vector(i, i, (props.MinMax[i + 3] - props.MinMax[i]) / props.NbSteps[i]);
        RI_cube.set_origin(i, props.MinMax[i]);
        RI_cube.set_vector(i, i, (props.MinMax[i + 3] - props.MinMax[i]) / props.NbSteps[i]);
    }
    WFN_cube.calc_dv();
    RI_cube.calc_dv();


    std::cout << "Starting work..." << std::endl;
    WFN_cube.set_comment1("Calculated density using NoSpherA2 for WFN");
    RI_cube.set_comment1("Calculated density using SALTED RI");
    WFN_cube.set_path((dummy.get_path().parent_path() / dummy.get_path().stem()).string() + "_rho_WFN.cube");
    RI_cube.set_path((dummy.get_path().parent_path() / dummy.get_path().stem()).string() + "_rho_RI.cube");

    std::cout << "Calculating WFN density cube..." << std::endl;
    wavy.calc_rho_cube(WFN_cube);
    WFN_cube.write_file(false, false);
    std::cout << "Calculating RI density cube..." << std::endl;
    calc_cube_ML(coeffs, dummy_aux, RI_cube);
    RI_cube.write_file(false, false);


    //Calculate difference cube
    cube diff = WFN_cube - RI_cube;
    diff.calc_dv(); 
    double rrs = WFN_cube.rrs(RI_cube);
    vec sums = diff.double_sum();
    std::cout << "RRS = " << std::scientific << std::setprecision(6) << rrs
        << ", Sum = " << std::fixed << std::setprecision(6) << sums[0]
        << ", |Sum| = " << std::fixed << std::setprecision(6) << sums[1] << std::endl;
    diff.set_comment1("Difference cube: WFN - RI");
    diff.give_parent_wfn(dummy);
    diff.set_path((dummy.get_path().parent_path() / dummy.get_path().stem()).string() + "_rho_diff_RI.cube");
    diff.write_file(true, false);
}