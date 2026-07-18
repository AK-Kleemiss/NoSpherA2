#pragma once
#include "convenience.h"

// Distributional (Gaussian) halting criterion for the XCW/XRW lambda scan.
//
// Reference: tests/P1_test/XCW_plan.md ("Distributional (Gaussian) Halting
// Criterion for XCW/XRW"). At each converged lambda step, the standardized
// structure-factor residuals z_h = (|F_obs,h| - |F_calc,h|) / sigma_h are
// tested against N(0,1) using the Anderson-Darling statistic (fully
// specified reference distribution, not parameter-estimated). A minimum in
// A^2(lambda) indicates the point where the fit stops describing systematic
// model deficiency and starts absorbing noise; this is reported as the
// recommended halting lambda* = argmin A^2.
//
// This module only implements generic, dependency-free statistics; the XCW
// class (Src/core/XCW.cpp/.h) owns collecting z_h from F_obs/F_calc/sigma
// and orchestrating the per-lambda evaluation and reporting.

// One entry of the Gaussian halting diagnostics, computed once per
// converged lambda step (see XCW::evaluate_gaussian_halting).
struct GaussianHaltEntry {
    double lambda = 0.0;
    int n_total = 0;      // reflections available before any filtering
    int n_used = 0;        // reflections actually used (after strong cutoff)
    double sigma_scale = 1.0; // global rescale applied so <z^2> ~ 1 before the Gaussian test (see XCW_plan.md 4.1)

    double A2 = 0.0;        // Anderson-Darling statistic against N(0,1)
    bool ad_reject_5pct = false; // true if A2 exceeds the 5% critical value for a fully specified N(0,1) (2.492)

    double pp_slope = 0.0;    // normal probability plot (Abrahams-Keve) slope, ->1 expected
    double pp_intercept = 0.0; // ->0 expected

    double skewness = 0.0;    // ->0 expected
    double excess_kurtosis = 0.0; // ->0 expected
    double jarque_bera = 0.0;

    double resolution_trend_slope = 0.0; // slope of <z^2> vs resolution bin center
    double resolution_trend_r = 0.0;     // Spearman correlation of <z^2> vs resolution bin
    bool resolution_trend_flagged = false;

    double intensity_trend_slope = 0.0; // slope of <z^2> vs |F| decile bin center
    double intensity_trend_r = 0.0;
    bool intensity_trend_flagged = false;
};

// Standard normal CDF (Phi) and its inverse (probit), needed for the
// Anderson-Darling statistic and the normal probability plot respectively.
double std_normal_cdf(const double z);
double std_normal_inv_cdf(const double p);

// Anderson-Darling statistic of `z` against N(0,1) with fully specified
// (not estimated) mean/variance. `z` need not be pre-sorted.
double anderson_darling_statistic(vec z);

// 5% critical value for the Anderson-Darling test against a fully specified
// normal distribution (D'Agostino & Stephens 1986, Table 4.7, case 0).
inline constexpr double ANDERSON_DARLING_CRITICAL_5PCT = 2.492;

struct ProbabilityPlotFit {
    double slope = 1.0;
    double intercept = 0.0;
};
// Least-squares fit of sorted `z` against the expected standard normal
// order statistics Phi^-1((i-0.5)/n) (Abrahams & Keve 1971).
ProbabilityPlotFit normal_probability_plot_fit(vec z);

double sample_skewness(const vec& z);
double sample_excess_kurtosis(const vec& z);
// Combined omnibus test on skewness and kurtosis; large values indicate
// non-normality. No p-value lookup is performed, the raw statistic is
// logged for the user to interpret/compare across lambda.
double jarque_bera_statistic(const vec& z);

struct BinnedTrend {
    double slope = 0.0;
    double spearman_r = 0.0;
    bool flagged = false;
};
// Bins `z` into `n_bins` groups ordered by `key` (e.g. resolution or |F|),
// computes <z^2> per bin, and reports the linear-regression slope and
// Spearman rank correlation of <z^2> against the bin index (a flat, ~1
// trend is expected for i.i.d. residuals; XCW_plan.md 3.3). `flagged` is
// set when |spearman_r| exceeds `flag_threshold`.
BinnedTrend binned_z_squared_trend(const vec& z, const vec& key, int n_bins, double flag_threshold = 0.5);
