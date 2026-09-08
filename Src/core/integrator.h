#pragma once
#include "convenience.h"


namespace DensityFitting
{

    enum class RESTRAINT_TYPE {
        NONE,
        SIMPLE,
        SIMPLE_AND_TIK
    };

    enum class METRIC_TYPE {
        COULOMB,
        OVERLAP
    };

    enum class CHARGE_SCHEME {
        NUCLEAR,
        MULLIKEN,
        SANDERSON_ESTIMATE,
        TFVC,
        HIRSHFELD,
        MBIS,
        EMBIS
    };

    struct CONFIG {
        METRIC_TYPE metric = METRIC_TYPE::COULOMB; // Metric to use for density fitting
        bool analyze_quality = false; // Whether to analyze the quality of the density fitting
        RESTRAINT_TYPE restrain_type = RESTRAINT_TYPE::NONE; // Type of electron population restraints to apply

        //Next only neccecary if restraints are used
        double restraint_strength = 5.0e-5; // Base strength of electron population restraints
        double tikhonov_lambda = 1e-6;
        bool adaptive_restraint = true; // Whether to use adaptive weighting for restraints
        CHARGE_SCHEME charge_scheme = CHARGE_SCHEME::TFVC; // Scheme to use for calculating expected electron populations
        int multipole_lmax = -1; // >= 0: restrain the grid moments of charge_scheme's atoms up to this order with weight multipole_strength, replacing the adaptive weights
        double multipole_strength = 1.0;

        std::optional<ivec> asym_atm_list = std::nullopt; //Currently unsued till fixed!// Optional list of atom indices to only compute atoms actually present in the assymetic unit
    };


    vec density_fit(const WFN& wavy, const WFN& wavy_aux, const CONFIG& config);
    // Fit settings from the command line: -multipole_moments switches the restraints on
    CONFIG config_from_options(const options& opt);

    // Helper functions for charge analysis and restraints
    vec calculate_expected_populations(const WFN& wavy, const WFN& wavy_aux, const CHARGE_SCHEME & = CHARGE_SCHEME::NUCLEAR);
    // Grid moments of the partitioned density about each nucleus, [atom][l*l+l+m] for l = 0..lmax, electrons only
    vec2 calculate_expected_multipoles(const WFN& wavy, const CHARGE_SCHEME& scheme, const int lmax);

    void analyze_density_fit_quality(const vec& coefficients, const WFN& wavy_aux, const vec& expected_charges = vec());
    // Per-atom row weight of the restraints
    vec restraint_weights(const WFN& wavy_aux, const size_t n_aux, double base_restraint_coef = 0.00005, bool adaptive_weighting = true);
    void add_electron_restraint(vec& eri2c, vec& rho, const WFN& wavy_aux,
        const vec& atom_weights, const vec& expected_charges = vec());
    // Moment of one aux primitive about its own centre, Int r^l Y_lm chi = N c Gamma(l+3/2) / (2 alpha^(l+3/2))
    double radial_moment(const double exponent, const double coef, const int l);
    void add_multipole_restraint(vec& eri2c, vec& rho, const WFN& wavy_aux,
        const vec2& targets, const vec& atom_weights, const int lmax);
    // Moments of the fitted density, same layout as the targets, from the coefficients alone
    vec2 fitted_multipoles(const vec& coefficients, const WFN& wavy_aux, const int lmax);

    // First-order electrostatics between two fitted densities and their nuclei, Hartree; needs only the
    // coefficients and the aux basis, so RI-fitted and SALTED-predicted coefficients enter alike.
    // pair[a][b] over the atoms of A and B, rank[i][j] with 0 the nuclei and l+1 the aux functions of rank l.
    // Beyond electrostatics: pol_X = -1/2 sum alpha_a F_a^2 over the atoms of X with Thakkar polarizabilities in the
    // partner's field, disp the D4 energy of the dimer minus the monomers, overlap = Int rhoA rhoB and rep = K * overlap
    struct INTERACTION {
        double nuc_nuc = 0.0, nucA_rhoB = 0.0, nucB_rhoA = 0.0, rho_rho = 0.0;
        double pol_A = 0.0, pol_B = 0.0, disp = 0.0, overlap = 0.0, rep = 0.0;
        double electrostatic() const { return nuc_nuc + nucA_rhoB + nucB_rhoA + rho_rho; };
        double total() const { return electrostatic() + pol_A + pol_B + disp + rep; };
        vec2 pair, rank;
    };
    // Lower incomplete gamma function gamma(l+3/2, x)
    double lower_gamma_half(const int l, const double x);
    // Coulomb potential Int chi(r)/|r-R| of one aux primitive centred at the origin, Y_lm(R^) included
    double aux_potential(const double exponent, const double coef, const int l, const int m, const double* R);
    INTERACTION interaction_energy(const vec& coef_A, const WFN& aux_A, const vec& coef_B, const WFN& aux_B, const double repulsion_K = 0.0);
    void print_interaction_energy(const INTERACTION& E, const WFN& aux_A, const WFN& aux_B, std::ostream& file);

    // Demonstration function
    void demonstrate_enhanced_density_fitting(WFN& wavy, const WFN& wavy_aux);
    // Calculate the difference between the QM density and the RI density and write it to a cube file
    void QM_RI_difference_cube(WFN& wavy, const WFN& wavy_aux);

}