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
    };


    vec density_fit(const WFN& wavy, const WFN& wavy_aux, const CONFIG& config);

    //// Enhanced density fitting with adaptive electron restraints
    //vec density_fit_restrain(const WFN& wavy, const WFN& wavy_aux, const char metric,
    //                double restraint_strength = 0.00005, bool adaptive_restraint = true,
    //                const std::string& charge_scheme = "sanderson_estimate", bool analyze_quality = true);
    //
    //// Unrestrained density fitting (original approach)
    //vec density_fit_unrestrained(WFN& wavy, const WFN& wavy_aux,
    //                             const char metric, bool analyze_quality = true);
    //
    //vec density_fit_hybrid(const WFN& wavy, const WFN& wavy_aux,
    //                      const char metric, double restraint_strength = 0.00005,
    //                      double tikhonov_lambda = 1e-6, const std::string& charge_scheme = "mulliken", bool analyze_quality = false);

    // Helper functions for charge analysis and restraints
    vec calculate_expected_populations(const WFN& wavy, const WFN& wavy_aux, const CHARGE_SCHEME & = CHARGE_SCHEME::NUCLEAR);

    void analyze_density_fit_quality(const vec& coefficients, const WFN& wavy_aux, const vec& expected_charges = vec());
    void add_electron_restraint(vec& eri2c, vec& rho, const WFN& wavy_aux,
        double base_restraint_coef = 0.00005, bool adaptive_weighting = true,
        const vec& expected_charges = vec());

    // Demonstration function
    void demonstrate_enhanced_density_fitting(WFN& wavy, const WFN& wavy_aux);

}