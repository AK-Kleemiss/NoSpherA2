#pragma once
#include "convenience.h"

// Enhanced density fitting with adaptive electron restraints
vec density_fit_restrain(const WFN& wavy, const WFN& wavy_aux, const double max_mem, const char metric,
                double restraint_strength = 0.001, bool adaptive_restraint = true, 
                const std::string& charge_scheme = "mulliken_estimate", bool analyze_quality = true);

// Unrestrained density fitting (original approach)
vec density_fit_unrestrained(const WFN& wavy, const WFN& wavy_aux, const double max_mem, 
                             const char metric, bool analyze_quality = true);

vec density_fit_hybrid(const WFN& wavy, const WFN& wavy_aux, const double max_mem, 
                      const char metric, double restraint_strength = 0.001, 
                      double tikhonov_lambda = 1e-8, const std::string& charge_scheme = "nuclear", bool analyze_quality = true);

// Helper functions for charge analysis and restraints
vec calculate_expected_charges(const WFN& wavy, const WFN& wavy_aux, const std::string& scheme = "nuclear");
vec calculate_expected_charges(const dMatrix2& dm, const WFN& wavy);

void analyze_density_fit_quality(const vec& coefficients, const WFN& wavy_aux, const vec& expected_charges = vec());
void add_electron_restraint(vec& eri2c, vec& rho, const WFN& wavy_aux, 
                           double base_restraint_coef = 0.001, bool adaptive_weighting = true,
                           const vec& expected_charges = vec());

// Demonstration function
void demonstrate_enhanced_density_fitting(const WFN& wavy, const WFN& wavy_aux);