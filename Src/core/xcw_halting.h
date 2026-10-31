#pragma once
#include "convenience.h"
#include <limits>

// Gaussian halting criterion for the XCW/XRW lambda scan: at each converged lambda,
// z_h = (|F_obs,h| - |F_calc,h|) / sigma_h is tested against a fully specified N(0,1) by Anderson-Darling.
// lambda* = argmin A^2 is where the fit stops modelling systematic deficiency and starts absorbing noise.

// Diagnostics of one converged lambda step.
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

double std_normal_cdf(const double z);
double std_normal_inv_cdf(const double p);

// Against N(0,1) with fixed, not estimated, mean and variance; z need not be sorted.
double anderson_darling_statistic(vec z);

// D'Agostino & Stephens 1986, Table 4.7, case 0 (fully specified normal).
inline constexpr double ANDERSON_DARLING_CRITICAL_5PCT = 2.492;

struct ProbabilityPlotFit {
	double slope = 1.0;
	double intercept = 0.0;
};
// Sorted z against the expected normal order statistics Phi^-1((i-0.5)/n), Abrahams & Keve 1971.
ProbabilityPlotFit normal_probability_plot_fit(vec z);

double sample_skewness(const vec& z);
double sample_excess_kurtosis(const vec& z);
// Raw statistic, no p-value; only compared across lambda.
double jarque_bera_statistic(const vec& z);

struct BinnedTrend {
	double slope = 0.0;
	double spearman_r = 0.0;
	bool flagged = false;
};
// <z^2> per bin of z ordered by key (resolution or |F|), slope and Spearman r against the bin
// index; i.i.d. residuals give a flat ~1. flagged when |spearman_r| > flag_threshold.
BinnedTrend binned_z_squared_trend(const vec& z, const vec& key, int n_bins, double flag_threshold = 0.5);

struct PolynomialFit {
	bool valid = false;     // false if not enough points for this degree, or the fit is singular
	int degree = 0;
	vec coeffs;              // coeffs[k] is the x^k coefficient; y = Sum_k coeffs[k]*x^k
	double rss = 0.0;         // residual sum of squares, Sum (y_i - yhat_i)^2
	double r_squared = 0.0;   // 1 - RSS/TSS, coefficient of determination
	double aic = std::numeric_limits<double>::infinity(); // Akaike Information Criterion (lower is better)
	bool has_minimum = false; // true if the fitted curve has an interior local minimum
							   // (positive curvature) within the search range used by
							   // choose_best_polynomial_fit
	double vertex_x = 0.0;    // location of that minimum, only meaningful if has_minimum
	double vertex_y = 0.0;    // fitted value at vertex_x
};

// Normal equations, needs degree+3 points (degree+1 would be exactly determined and unstable).
// Minimum by grid search over [min(x), min(x) + 3*(max(x)-min(x))].
PolynomialFit fit_polynomial(const vec& x, const vec& y, int degree);

// Lowest AIC among the degrees with enough points; all_candidates gets every fit, invalid ones
// included, in the order of degrees.
PolynomialFit choose_best_polynomial_fit(const vec& x, const vec& y, const ivec& degrees,
	std::vector<PolynomialFit>* all_candidates = nullptr);

// True while the lowest A^2 of the scan still sits at its last evaluated
// lambda, i.e. the minimum lies beyond the scanned range and the scan would
// end on its own boundary rather than on a found lambda*. Entries with fewer
// than 8 usable reflections are skipped, as in the halting report.
// `estimated_minimum` receives the AIC-chosen polynomial's vertex when that
// lies beyond the last lambda, and 0 otherwise (too few points to fit, or a
// fit whose vertex contradicts the still-falling data) - an extension is
// warranted either way, there is just no target lambda to name.
bool halting_minimum_beyond_scan(const std::vector<GaussianHaltEntry>& history, double& estimated_minimum);
