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
        HIRSHFELD
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

        std::optional<ivec> asym_atm_list = std::nullopt; //Currently unsued till fixed!// Optional list of atom indices to only compute atoms actually present in the assymetic unit
    };


    vec density_fit(const WFN& wavy, const WFN& wavy_aux, const CONFIG& config);

    // Helper functions for charge analysis and restraints
    vec calculate_expected_populations(const WFN& wavy, const WFN& wavy_aux, const CHARGE_SCHEME & = CHARGE_SCHEME::NUCLEAR);

    void analyze_density_fit_quality(const vec& coefficients, const WFN& wavy_aux, const vec& expected_charges = vec());
    void add_electron_restraint(vec& eri2c, vec& rho, const WFN& wavy_aux,
        double base_restraint_coef = 0.00005, bool adaptive_weighting = true,
        const vec& expected_charges = vec());

    // Demonstration function
    void demonstrate_enhanced_density_fitting(WFN& wavy, const WFN& wavy_aux);
    // Calculate the difference between the QM density and the RI density and write it to a cube file
    void QM_RI_difference_cube(WFN& wavy, const WFN& wavy_aux);

}