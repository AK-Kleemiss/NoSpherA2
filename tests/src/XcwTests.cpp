//The XCW lambda-scan halting statistics (Src/core/xcw_halting.cpp) against analytic
//values and goldens derived by hand, and the XCW settings-file parser through the only
//door it has, the XCW constructor: the parse errors it throws need no crystal, the
//successful parses run construct() on the P1 fixture and are checked through XCW.log.
#include "pch.h"
#include <gtest/gtest.h>

#include "core/convenience.h"
#include "core/xcw_halting.h"
#include "core/XCW.h"

namespace
{
	//the ideal N(0,1) sample: the plotting-position quantiles z_i = Phi^-1((i + 0.5) / n),
	//optionally scaled and shifted, in a scrambled order so nothing relies on sorted input
	vec quantile_sample(const int n, const double scale = 1.0, const double shift = 0.0)
	{
		vec z(n);
		for (int i = 0; i < n; i++) {
			const int k = (i * 37) % n; // 37 is coprime to every n used here
			z[i] = scale * std_normal_inv_cdf((k + 0.5) / n) + shift;
		}
		return z;
	}

	//van der Corput base-2 points through the probit: a deterministic, unsorted sample that
	//fills the normal distribution evenly without being its exact quantiles
	vec vdc_sample(const int n, const double shift = 0.0)
	{
		vec z(n);
		for (int i = 1; i <= n; i++) {
			double x = 0.0, f = 0.5;
			for (int m = i; m > 0; m >>= 1, f *= 0.5) x += f * (m & 1);
			z[i - 1] = std_normal_inv_cdf(x) + shift;
		}
		return z;
	}

	//scratch paths carry the running test's name: ctest runs every case in its own
	//process, in parallel, and a shared name lets one case remove another's file
	std::string test_name()
	{
		return testing::UnitTest::GetInstance()->current_test_info()->name();
	}

	std::filesystem::path write_settings(const std::string& text)
	{
		const auto p = std::filesystem::temp_directory_path() / ("nosphera2_xcw_tests_settings_" + test_name() + ".txt");
		std::ofstream(p) << text;
		return p;
	}

	//what loadSettings throws for a settings text; every text handed in here must throw,
	//because a parse that succeeds runs construct() on an empty cif and that exits the process
	std::string parse_error(const std::string& text)
	{
		const auto p = write_settings(text);
		options opt;
		opt.xcw_settings_path = p;
		opt.do_XCW = true;
		std::string what;
		try {
			XCW x(opt);
		}
		catch (const std::runtime_error& e) {
			what = e.what();
		}
		std::filesystem::remove(p);
		return what;
	}

	//construct an XCW on the P1 fixture with this settings text, in a scratch directory
	//because construct() writes XCW.log and log3.txt into the working directory; returns
	//the XCW.log text, empty when the fixture is missing
	std::string construct_on_p1(const std::string& text)
	{
		const auto fixture = nos_test_repo_root() / "tests" / "P1_test";
		if (!std::filesystem::exists(fixture / "P1_test_NA2.cif") || !std::filesystem::exists(fixture / "P1_test.hkl")) {
			return "";
		}
		const auto dir = std::filesystem::temp_directory_path() / ("nosphera2_xcw_tests_p1_" + test_name());
		std::filesystem::create_directories(dir);
		const auto settings = dir / "settings.txt";
		std::ofstream(settings) << text;
		options opt;
		opt.xcw_settings_path = settings;
		opt.cif = std::filesystem::absolute(fixture / "P1_test_NA2.cif");
		opt.hkl = std::filesystem::absolute(fixture / "P1_test.hkl");
		opt.do_XCW = true;
		const auto old_cwd = std::filesystem::current_path();
		std::filesystem::current_path(dir);
		{
			XCW x(opt);
		}
		std::filesystem::current_path(old_cwd);
		std::stringstream log;
		log << std::ifstream(dir / "XCW.log").rdbuf();
		std::filesystem::remove_all(dir);
		std::string text_out = log.str();
		std::transform(text_out.begin(), text_out.end(), text_out.begin(), [](unsigned char c) { return (char)std::tolower(c); });
		return text_out;
	}

	// ------------------------------------------------------------------
	// Normal distribution functions

	// Phi at tabulated points, and Phi(z) + Phi(-z) = 1
	TEST(XcwHaltingTests, StandardNormalCdfMatchesTables)
	{
		EXPECT_NEAR(std_normal_cdf(0.0), 0.5, 1e-15);
		EXPECT_NEAR(std_normal_cdf(1.0), 0.841344746068543, 1e-12);
		EXPECT_NEAR(std_normal_cdf(-1.0), 0.158655253931457, 1e-12);
		EXPECT_NEAR(std_normal_cdf(1.959963984540054), 0.975, 1e-12);
		EXPECT_NEAR(std_normal_cdf(-3.0), 0.001349898031630, 1e-12);
		for (const double z : { -4.5, -0.7, 0.3, 2.2 }) {
			EXPECT_NEAR(std_normal_cdf(z) + std_normal_cdf(-z), 1.0, 1e-15);
		}
	}

	// Acklam's rational approximation in all three of its branches (lower tail, centre,
	// upper tail) against tabulated probits, to its stated 1.15e-9 accuracy
	TEST(XcwHaltingTests, ProbitMatchesTablesInAllThreeBranches)
	{
		EXPECT_NEAR(std_normal_inv_cdf(0.5), 0.0, 1e-9);
		EXPECT_NEAR(std_normal_inv_cdf(0.975), 1.959963984540054, 2e-9);
		EXPECT_NEAR(std_normal_inv_cdf(0.3), -0.524400512708041, 2e-9);
		EXPECT_NEAR(std_normal_inv_cdf(0.01), -2.326347874040841, 2e-9);
		EXPECT_NEAR(std_normal_inv_cdf(0.02), -2.053748910631823, 2e-9);
		EXPECT_NEAR(std_normal_inv_cdf(0.999), 3.090232306167813, 2e-9);
		EXPECT_NEAR(std_normal_inv_cdf(0.999), -std_normal_inv_cdf(0.001), 1e-12);
	}

	// the probit inverts the cdf across the range, including both tail branches
	TEST(XcwHaltingTests, ProbitInvertsTheCdf)
	{
		for (double p = 0.001; p < 0.9995; p += 0.0123) {
			EXPECT_NEAR(std_normal_cdf(std_normal_inv_cdf(p)), p, 2e-9) << "p=" << p;
		}
		for (const double z : { -3.5, -2.0, -0.1, 0.4, 1.7, 3.2 }) {
			EXPECT_NEAR(std_normal_inv_cdf(std_normal_cdf(z)), z, 1e-7) << "z=" << z;
		}
	}

	// ------------------------------------------------------------------
	// Anderson-Darling

	// A^2 of {-1, 1} worked by hand: sum = 2 ln Phi(-1) + 6 ln Phi(1), A^2 = -2 - sum / 2
	TEST(XcwHaltingTests, AndersonDarlingOfTwoPointsMatchesHandCalculation)
	{
		EXPECT_NEAR(anderson_darling_statistic({ 1.0, -1.0 }), 0.359282982, 1e-8);
	}

	// the ideal N(0,1) sample of 100 quantiles gives the golden A^2 (formula evaluated
	// independently in Python), far below the 5 % critical value, in any input order
	TEST(XcwHaltingTests, AndersonDarlingAcceptsTheExactNormalQuantiles)
	{
		const vec z = quantile_sample(100);
		const double a2 = anderson_darling_statistic(z);
		EXPECT_NEAR(a2, 0.0114951327, 1e-8);
		EXPECT_LT(a2, ANDERSON_DARLING_CRITICAL_5PCT);
		vec sorted = z;
		std::sort(sorted.begin(), sorted.end());
		EXPECT_NEAR(anderson_darling_statistic(sorted), a2, 1e-14);
	}

	// the same quantiles shifted by one sigma, or scaled by 2 or 1/2, are rejected: the
	// test is against a fully specified N(0,1), so both location and scale count
	TEST(XcwHaltingTests, AndersonDarlingRejectsShiftedAndScaledSamples)
	{
		EXPECT_NEAR(anderson_darling_statistic(quantile_sample(100, 1.0, 1.0)), 47.8319577, 1e-6);
		EXPECT_NEAR(anderson_darling_statistic(quantile_sample(100, 2.0)), 20.4411291, 1e-6);
		EXPECT_NEAR(anderson_darling_statistic(quantile_sample(100, 0.5)), 8.1148731, 1e-6);
		EXPECT_GT(anderson_darling_statistic(quantile_sample(100, 1.0, 1.0)), ANDERSON_DARLING_CRITICAL_5PCT);
	}

	// a low-discrepancy normal sample that is not the exact quantiles: accepted, and its
	// half-sigma shift rejected, both against the independently evaluated formula
	TEST(XcwHaltingTests, AndersonDarlingOnVanDerCorputSample)
	{
		EXPECT_NEAR(anderson_darling_statistic(vdc_sample(200)), 0.0619843352, 1e-8);
		EXPECT_NEAR(anderson_darling_statistic(vdc_sample(200, 0.5)), 22.1947679, 1e-6);
	}

	// values far out in the tails saturate Phi at 0 or 1; the clamp keeps A^2 finite.
	// Without the clamp log(0) makes A^2 +inf, which is what isfinite catches; the finite
	// value itself is set by the clamp constants (1e-300, 1e-16), not by anything with a
	// statistical meaning, so it is not pinned down
	TEST(XcwHaltingTests, AndersonDarlingStaysFiniteForExtremeOutliers)
	{
		vec z = quantile_sample(50);
		z.push_back(40.0);
		z.push_back(-40.0);
		const double a2 = anderson_darling_statistic(z);
		EXPECT_TRUE(std::isfinite(a2));
		EXPECT_GT(a2, ANDERSON_DARLING_CRITICAL_5PCT);
	}

	// ------------------------------------------------------------------
	// Probability plot, moments

	// sorted z against the expected order statistics: a sample built as a * quantiles + b
	// fits slope a and intercept b exactly, whatever order it comes in
	TEST(XcwHaltingTests, ProbabilityPlotRecoversScaleAndShift)
	{
		const ProbabilityPlotFit ideal = normal_probability_plot_fit(quantile_sample(64));
		EXPECT_NEAR(ideal.slope, 1.0, 1e-12);
		EXPECT_NEAR(ideal.intercept, 0.0, 1e-12);
		const ProbabilityPlotFit fit = normal_probability_plot_fit(quantile_sample(64, 1.5, -0.3));
		EXPECT_NEAR(fit.slope, 1.5, 1e-12);
		EXPECT_NEAR(fit.intercept, -0.3, 1e-12);
	}

	// Bernoulli(1/4) as {0, 0, 0, 1}: skewness (1 - 2p) / sqrt(p(1 - p)) = 2 / sqrt(3) and
	// excess kurtosis (1 - 6p(1 - p)) / (p(1 - p)) = -2/3; {-1, -1, 1, 1} is the flattest
	// possible sample with excess kurtosis -2
	TEST(XcwHaltingTests, SkewnessAndKurtosisOfBernoulliSamples)
	{
		const vec bern = { 0.0, 0.0, 0.0, 1.0 };
		EXPECT_NEAR(sample_skewness(bern), 2.0 / std::sqrt(3.0), 1e-12);
		EXPECT_NEAR(sample_excess_kurtosis(bern), -2.0 / 3.0, 1e-12);
		const vec flat = { -1.0, 1.0, -1.0, 1.0 };
		EXPECT_NEAR(sample_skewness(flat), 0.0, 1e-15);
		EXPECT_NEAR(sample_excess_kurtosis(flat), -2.0, 1e-15);
	}

	// a mirror-symmetric sample has zero skewness; the exact normal quantiles have a
	// slightly negative excess kurtosis (the plotting positions truncate the tails),
	// -0.0990289743 for n = 200 by the same formula evaluated independently in Python
	TEST(XcwHaltingTests, MomentsOfSymmetricSamples)
	{
		const vec q = quantile_sample(200);
		EXPECT_NEAR(sample_skewness(q), 0.0, 1e-12);
		EXPECT_NEAR(sample_excess_kurtosis(q), -0.0990289743, 1e-7);
		EXPECT_NEAR(sample_skewness(vec{ -2.0, -1.0, 0.0, 1.0, 2.0 }), 0.0, 1e-15);
	}

	// a constant sample has no spread: both moments fall back to 0 instead of 0/0
	TEST(XcwHaltingTests, MomentsOfAConstantSampleAreZero)
	{
		const vec c(10, 3.5);
		EXPECT_EQ(sample_skewness(c), 0.0);
		EXPECT_EQ(sample_excess_kurtosis(c), 0.0);
		EXPECT_EQ(jarque_bera_statistic(c), 0.0);
	}

	// JB = n/6 (S^2 + K^2/4) with the Bernoulli moments above: 4/6 (4/3 + 1/9) = 26/27;
	// one 12 among eleven 0s is Bernoulli(1/12): S = 10/sqrt(11), K = 78/11, so
	// JB = 2 (100/11 + 6084/484) = 5242/121, whatever the outlier's size (moments are
	// scale-free); the 100 quantiles give 0.1145946890 by the Python evaluation
	TEST(XcwHaltingTests, JarqueBeraCombinesSkewnessAndKurtosis)
	{
		EXPECT_NEAR(jarque_bera_statistic({ 0.0, 0.0, 0.0, 1.0 }), 26.0 / 27.0, 1e-12);
		EXPECT_NEAR(jarque_bera_statistic(quantile_sample(100)), 0.1145946890, 1e-7);
		EXPECT_NEAR(jarque_bera_statistic({ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0 }), 5242.0 / 121.0, 1e-12);
	}

	// ------------------------------------------------------------------
	// Binned <z^2> trend

	// z^2 = 1 + 0.1 k for k = 0..59 in six bins of ten by key: bin b holds k = 10b..10b+9,
	// its mean is 1.45 + b, so the slope against the bin index is exactly 1 and Spearman +1,
	// flagged; the key handed in scrambled so the ordering by key is what gets tested
	TEST(XcwHaltingTests, BinnedTrendDetectsGrowthWithTheKey)
	{
		const int n = 60;
		vec z(n), key(n);
		for (int i = 0; i < n; i++) {
			const int k = (i * 37) % n;
			key[i] = 0.1 * k;
			z[i] = std::sqrt(1.0 + 0.1 * k);
		}
		const BinnedTrend t = binned_z_squared_trend(z, key, 6);
		EXPECT_NEAR(t.spearman_r, 1.0, 1e-12);
		EXPECT_NEAR(t.slope, 1.0, 1e-12);
		EXPECT_TRUE(t.flagged);
		//the same z against the reversed key: bin b holds k = 59-10b..50-10b, mean 6.45 - b
		vec reversed(n);
		for (int i = 0; i < n; i++) reversed[i] = -key[i];
		const BinnedTrend r = binned_z_squared_trend(z, reversed, 6);
		EXPECT_NEAR(r.spearman_r, -1.0, 1e-12);
		EXPECT_NEAR(r.slope, -1.0, 1e-12);
		EXPECT_TRUE(r.flagged);
	}

	// |z| = 1 everywhere: bin means all 1, slope 0, Spearman 0 by the zero-variance
	// fallback, not flagged
	TEST(XcwHaltingTests, BinnedTrendOfFlatResidualsIsZero)
	{
		vec z(40), key(40);
		for (int i = 0; i < 40; i++) {
			z[i] = (i % 2) ? 1.0 : -1.0;
			key[i] = 0.05 * i;
		}
		const BinnedTrend t = binned_z_squared_trend(z, key, 4);
		EXPECT_NEAR(t.slope, 0.0, 1e-14);
		EXPECT_EQ(t.spearman_r, 0.0);
		EXPECT_FALSE(t.flagged);
	}

	// four bins with means {1, 1, 2, 2} against index {0, 1, 2, 3}: average ranks for the
	// ties give Spearman 4 / sqrt(20) and the least-squares slope 0.4, worked by hand; n_bins
	// is clamped to n so asking for 10 bins of 4 points yields one point per bin
	TEST(XcwHaltingTests, BinnedTrendSpearmanUsesAverageRanksForTies)
	{
		const vec z = { std::sqrt(2.0), 1.0, std::sqrt(2.0), 1.0 };
		const vec key = { 3.0, 0.0, 2.0, 1.0 };
		const BinnedTrend t = binned_z_squared_trend(z, key, 10, 0.95);
		EXPECT_NEAR(t.spearman_r, 4.0 / std::sqrt(20.0), 1e-12);
		EXPECT_NEAR(t.slope, 0.4, 1e-12);
		EXPECT_FALSE(t.flagged) << "0.894 is below the 0.95 threshold";
		EXPECT_TRUE(binned_z_squared_trend(z, key, 10, 0.5).flagged);
	}

	// 5 points into 2 bins: the remainder goes to the first bin (3 + 2), so the bin means
	// are (1 + 4 + 9) / 3 and (16 + 25) / 2 and the slope is their difference
	TEST(XcwHaltingTests, BinnedTrendSpreadsTheRemainderOverTheFirstBins)
	{
		const vec z = { 1.0, 2.0, 3.0, 4.0, 5.0 };
		const vec key = { 1.0, 2.0, 3.0, 4.0, 5.0 };
		const BinnedTrend t = binned_z_squared_trend(z, key, 2);
		EXPECT_NEAR(t.slope, 20.5 - 14.0 / 3.0, 1e-12);
		EXPECT_NEAR(t.spearman_r, 1.0, 1e-12);
		//a single bin is never asked for: n_bins = 1 becomes 2
		const BinnedTrend one = binned_z_squared_trend(z, key, 1);
		EXPECT_NEAR(one.slope, t.slope, 1e-12);
	}

	// ------------------------------------------------------------------
	// Polynomial fits of A^2(lambda)

	// exact quadratic data y = x^2 - 3x + 2 through six points: coefficients recovered,
	// RSS ~ 0 so R^2 = 1 and AIC = -inf, and the vertex at (1.5, -0.25) found by the grid
	// search to its 15/1999 spacing
	TEST(XcwHaltingTests, PolynomialFitRecoversAnExactQuadratic)
	{
		const vec x = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
		vec y(x.size());
		for (size_t i = 0; i < x.size(); i++) y[i] = x[i] * x[i] - 3.0 * x[i] + 2.0;
		const PolynomialFit fit = fit_polynomial(x, y, 2);
		ASSERT_TRUE(fit.valid);
		EXPECT_EQ(fit.degree, 2);
		ASSERT_EQ(fit.coeffs.size(), 3u);
		EXPECT_NEAR(fit.coeffs[0], 2.0, 1e-9);
		EXPECT_NEAR(fit.coeffs[1], -3.0, 1e-9);
		EXPECT_NEAR(fit.coeffs[2], 1.0, 1e-9);
		EXPECT_LT(fit.rss, 1e-18);
		EXPECT_NEAR(fit.r_squared, 1.0, 1e-12);
		EXPECT_TRUE(fit.has_minimum);
		EXPECT_NEAR(fit.vertex_x, 1.5, 0.01);
		EXPECT_NEAR(fit.vertex_y, -0.25, 1e-3);
	}

	// a fit needs degree + 3 points: four points cannot carry a quadratic, five can
	TEST(XcwHaltingTests, PolynomialFitRefusesTooFewPoints)
	{
		const vec x = { 0.0, 1.0, 2.0, 3.0 };
		const vec y = { 1.0, 2.0, 5.0, 10.0 };
		const PolynomialFit none = fit_polynomial(x, y, 2);
		EXPECT_FALSE(none.valid);
		EXPECT_EQ(none.degree, 2);
		EXPECT_TRUE(none.coeffs.empty());
		EXPECT_FALSE(none.has_minimum);
		EXPECT_TRUE(std::isinf(none.aic));
		const PolynomialFit line = fit_polynomial(x, y, 1);
		EXPECT_TRUE(line.valid);
	}

	// a line has no interior minimum in either direction: rising, the grid minimum sits at
	// the left boundary, falling, at the right one; the coefficients are still exact
	TEST(XcwHaltingTests, PolynomialFitOfALineHasNoMinimum)
	{
		const vec x = { 0.0, 0.5, 1.0, 1.5, 2.0 };
		vec up(x.size()), down(x.size());
		for (size_t i = 0; i < x.size(); i++) {
			up[i] = 1.0 + 2.0 * x[i];
			down[i] = 1.0 - 2.0 * x[i];
		}
		const PolynomialFit a = fit_polynomial(x, up, 1);
		ASSERT_TRUE(a.valid);
		EXPECT_NEAR(a.coeffs[0], 1.0, 1e-12);
		EXPECT_NEAR(a.coeffs[1], 2.0, 1e-12);
		EXPECT_FALSE(a.has_minimum);
		const PolynomialFit b = fit_polynomial(x, down, 1);
		ASSERT_TRUE(b.valid);
		EXPECT_NEAR(b.coeffs[1], -2.0, 1e-12);
		EXPECT_FALSE(b.has_minimum);
	}

	// y = (x - 8)^2 sampled on [0, 5]: the minimum is beyond the data but inside the search
	// window of three spans, which is exactly the extrapolation the halting report uses;
	// -(x - 2)^2 curves the wrong way and its grid minimum is a boundary, so no minimum
	TEST(XcwHaltingTests, PolynomialFitExtrapolatesAMinimumBeyondTheData)
	{
		const vec x = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
		vec convex(x.size()), concave(x.size());
		for (size_t i = 0; i < x.size(); i++) {
			convex[i] = (x[i] - 8.0) * (x[i] - 8.0);
			concave[i] = -(x[i] - 2.0) * (x[i] - 2.0);
		}
		const PolynomialFit a = fit_polynomial(x, convex, 2);
		ASSERT_TRUE(a.valid);
		EXPECT_TRUE(a.has_minimum);
		EXPECT_NEAR(a.vertex_x, 8.0, 0.01);
		EXPECT_NEAR(a.vertex_y, 0.0, 1e-3);
		const PolynomialFit b = fit_polynomial(x, concave, 2);
		ASSERT_TRUE(b.valid);
		EXPECT_NEAR(b.coeffs[2], -1.0, 1e-9);
		EXPECT_FALSE(b.has_minimum);
	}

	// noisy quadratic y = (x - 3)^2 + 0.3 (-1)^i on ten points: RSS, R^2 and the AIC
	// n ln(RSS/n) + 2(k + 1) against numpy.polyfit
	TEST(XcwHaltingTests, PolynomialFitStatisticsOfNoisyData)
	{
		vec x(10), y(10);
		double mean_y = 0.0;
		for (int i = 0; i < 10; i++) {
			x[i] = i;
			y[i] = (i - 3.0) * (i - 3.0) + 0.3 * ((i % 2) ? -1.0 : 1.0);
			mean_y += y[i] / 10.0;
		}
		double tss = 0.0;
		for (const double v : y) tss += (v - mean_y) * (v - mean_y);
		const PolynomialFit fit = fit_polynomial(x, y, 2);
		ASSERT_TRUE(fit.valid);
		EXPECT_NEAR(fit.coeffs[0], 9.08181818, 1e-6);
		EXPECT_NEAR(fit.coeffs[1], -6.01818182, 1e-6);
		EXPECT_NEAR(fit.coeffs[2], 1.0, 1e-6);
		EXPECT_NEAR(fit.rss, 0.87272727, 1e-6);
		EXPECT_NEAR(fit.r_squared, 1.0 - 0.87272727 / tss, 1e-6);
		EXPECT_NEAR(fit.aic, -16.38717267, 1e-5);
		EXPECT_TRUE(fit.has_minimum);
		EXPECT_NEAR(fit.vertex_x, 6.01818182 / 2.0, 0.02);
	}

	// with six points only the quadratic has enough data: it is chosen, the quartic is
	// reported as an invalid candidate in the order of the degrees asked for
	TEST(XcwHaltingTests, BestFitSkipsDegreesWithoutEnoughPoints)
	{
		const vec x = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
		vec y(x.size());
		for (size_t i = 0; i < x.size(); i++) y[i] = (x[i] - 2.0) * (x[i] - 2.0) + 1.0;
		std::vector<PolynomialFit> candidates;
		const PolynomialFit best = choose_best_polynomial_fit(x, y, { 2, 4 }, &candidates);
		ASSERT_TRUE(best.valid);
		EXPECT_EQ(best.degree, 2);
		EXPECT_NEAR(best.vertex_x, 2.0, 0.01);
		ASSERT_EQ(candidates.size(), 2u);
		EXPECT_TRUE(candidates[0].valid);
		EXPECT_EQ(candidates[0].degree, 2);
		EXPECT_FALSE(candidates[1].valid);
		EXPECT_EQ(candidates[1].degree, 4);
		//nothing to try, nothing valid
		EXPECT_FALSE(choose_best_polynomial_fit(x, y, {}).valid);
		EXPECT_FALSE(choose_best_polynomial_fit(x, y, { 4 }).valid);
	}

	// AIC arbitrates: exact quartic data on eight points makes the quartic win outright,
	// the noisy quadratic of the test above makes the quadratic win although the quartic
	// has the lower RSS (0.806 vs 0.873), because its two extra parameters cost more
	TEST(XcwHaltingTests, BestFitIsChosenByAicNotByRss)
	{
		vec x(8), y(8);
		for (int i = 0; i < 8; i++) {
			x[i] = i;
			y[i] = std::pow(i, 4.0);
		}
		std::vector<PolynomialFit> quartic_candidates;
		const PolynomialFit quartic = choose_best_polynomial_fit(x, y, { 2, 4 }, &quartic_candidates);
		ASSERT_TRUE(quartic.valid);
		EXPECT_EQ(quartic.degree, 4);
		EXPECT_NEAR(quartic.coeffs[4], 1.0, 1e-6);
		EXPECT_GT(quartic_candidates[0].rss, 1e4);
		EXPECT_LT(quartic_candidates[1].rss, 1e-4);

		vec xn(10), yn(10);
		for (int i = 0; i < 10; i++) {
			xn[i] = i;
			yn[i] = (i - 3.0) * (i - 3.0) + 0.3 * ((i % 2) ? -1.0 : 1.0);
		}
		std::vector<PolynomialFit> noisy_candidates;
		const PolynomialFit quadratic = choose_best_polynomial_fit(xn, yn, { 2, 4 }, &noisy_candidates);
		ASSERT_TRUE(quadratic.valid);
		EXPECT_EQ(quadratic.degree, 2);
		EXPECT_NEAR(noisy_candidates[1].rss, 0.80559441, 1e-6);
		EXPECT_NEAR(noisy_candidates[1].aic, -13.18759975, 1e-5);
		EXPECT_LT(noisy_candidates[1].rss, noisy_candidates[0].rss);
		EXPECT_GT(noisy_candidates[1].aic, noisy_candidates[0].aic);
	}

	// every x equal makes the normal-equations matrix exactly singular; such a fit comes back invalid
	TEST(XcwHaltingTests, PolynomialFitOfASingularSystemIsInvalid)
	{
		const vec x(6, 1.0);
		const vec y = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
		EXPECT_FALSE(fit_polynomial(x, y, 2).valid);
	}

	// ------------------------------------------------------------------
	// XCW settings file, through the constructor

	// a missing settings file is refused before anything else is opened
	//transform_ADPs takes M and applies T'_{ij..} = sum M_pi M_qj .. T_pq.. : for the rank-2 U that is
	//M^T U M, checked against a hand-multiplied three-fold (which tells M from M^T), and for C and D the
	//sign pattern of a two-fold along b (every 0 or 2 index flips the sign)
	TEST(XcwAdpTests, TransformAdpsRotatesUCDAsContravariantTensors)
	{
		vec2 adps = { { 0.01, 0.02, 0.03, 0.004, 0.005, 0.006 },
					  { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 },
					  { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 } };
		const vec2 U = { { adps[0][0], adps[0][3], adps[0][4] },
						 { adps[0][3], adps[0][1], adps[0][5] },
						 { adps[0][4], adps[0][5], adps[0][2] } };
		// x' = -y, y' = x - y, z' = z handed in transposed, as cell stores it
		const vec2 M = { { 0, 1, 0 }, { -1, -1, 0 }, { 0, 0, 1 } };
		vec2 rotated = adps;
		transform_ADPs(rotated, M);
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				double expect = 0.0;
				for (int p = 0; p < 3; p++)
					for (int q = 0; q < 3; q++)
						expect += M[p][i] * M[q][j] * U[p][q];
				const int voigt[3][3] = { { 0, 3, 4 }, { 3, 1, 5 }, { 4, 5, 2 } };
				EXPECT_NEAR(rotated[0][voigt[i][j]], expect, 1e-15) << i << j;
			}
		}
		// U is no longer the input: the three-fold mixes a and b
		EXPECT_NE(rotated[0][0], adps[0][0]);
		// the higher ranks were transformed too (a three-fold moves C111 into a mix of C222 and friends)
		EXPECT_NE(rotated[1][0], adps[1][0]);
		EXPECT_NE(rotated[2][0], adps[2][0]);

		const vec2 twofold_b = { { -1, 0, 0 }, { 0, 1, 0 }, { 0, 0, -1 } };
		rotated = adps;
		transform_ADPs(rotated, twofold_b);
		const int map3[10][3] = { { 0, 0, 0 }, { 0, 0, 1 }, { 0, 0, 2 }, { 0, 1, 1 }, { 0, 1, 2 }, { 0, 2, 2 }, { 1, 1, 1 }, { 1, 1, 2 }, { 1, 2, 2 }, { 2, 2, 2 } };
		const int map4[15][4] = { { 0, 0, 0, 0 }, { 0, 0, 0, 1 }, { 0, 0, 0, 2 }, { 0, 0, 1, 1 }, { 0, 0, 1, 2 }, { 0, 0, 2, 2 }, { 0, 1, 1, 1 }, { 0, 1, 1, 2 }, { 0, 1, 2, 2 }, { 0, 2, 2, 2 }, { 1, 1, 1, 1 }, { 1, 1, 1, 2 }, { 1, 1, 2, 2 }, { 1, 2, 2, 2 }, { 2, 2, 2, 2 } };
		const double sign_u[6] = { 1, 1, 1, -1, 1, -1 };
		for (int i = 0; i < 6; i++)
			EXPECT_DOUBLE_EQ(rotated[0][i], sign_u[i] * adps[0][i]) << i;
		for (int i = 0; i < 10; i++) {
			double sign = 1.0;
			for (int k = 0; k < 3; k++) sign *= twofold_b[map3[i][k]][map3[i][k]];
			EXPECT_DOUBLE_EQ(rotated[1][i], sign * adps[1][i]) << i;
		}
		for (int i = 0; i < 15; i++) {
			double sign = 1.0;
			for (int k = 0; k < 4; k++) sign *= twofold_b[map4[i][k]][map4[i][k]];
			EXPECT_DOUBLE_EQ(rotated[2][i], sign * adps[2][i]) << i;
		}

		// an atom without ADPs is left alone
		vec2 none(3);
		transform_ADPs(none, M);
		EXPECT_TRUE(none[0].empty() && none[1].empty() && none[2].empty());
	}

	TEST(XcwSettingsTests, MissingSettingsFileThrows)
	{
		options opt;
		opt.xcw_settings_path = std::filesystem::temp_directory_path() / "nosphera2_xcw_tests_nowhere.txt";
		opt.do_XCW = true;
		std::string what;
		try {
			XCW x(opt);
		}
		catch (const std::runtime_error& e) {
			what = e.what();
		}
		EXPECT_EQ(what, "Settings file not found! Aborting run!");
	}

	// keywords are matched case-insensitively and an unknown one names itself, lowered
	TEST(XcwSettingsTests, UnknownKeywordIsNamedInTheError)
	{
		EXPECT_EQ(parse_error("Basis_Set sto-3g Bogus 1"), "Unknown keyword 'bogus'");
		EXPECT_EQ(parse_error("basis_set sto-3g -params 1"), "Unknown keyword '-params'");
	}

	// a value keyword at the end of the file, or followed by a non-number, has no value
	TEST(XcwSettingsTests, MissingValuesAreReported)
	{
		EXPECT_EQ(parse_error("basis_set sto-3g conv"), "Expected value after 'conv'");
		EXPECT_EQ(parse_error("basis_set sto-3g params x"), "Expected value after 'params'");
		EXPECT_EQ(parse_error("basis_set sto-3g max_iter"), "Expected value after 'max_iter'");
		EXPECT_EQ(parse_error("basis_set sto-3g charge"), "Expected value after 'charge'");
		EXPECT_EQ(parse_error("basis_set sto-3g mult"), "Expected value after 'mult'");
		EXPECT_EQ(parse_error("basis_set sto-3g damp"), "Expected value after 'damp'");
		EXPECT_EQ(parse_error("basis_set sto-3g shift"), "Expected value after 'shift'");
		EXPECT_EQ(parse_error("basis_set sto-3g start"), "Expected value after 'start'");
		EXPECT_EQ(parse_error("basis_set sto-3g end"), "Expected value after 'end'");
		EXPECT_EQ(parse_error("basis_set sto-3g step_size"), "Expected value after 'step_size'");
		EXPECT_EQ(parse_error("basis_set sto-3g diis_damping"), "Expected value after 'diis_damping'");
		EXPECT_EQ(parse_error("basis_set sto-3g diis_shift"), "Expected value after 'diis_shift'");
		EXPECT_EQ(parse_error("basis_set sto-3g conv_diis"), "Expected value after 'conv_diis'");
		EXPECT_EQ(parse_error("basis_set sto-3g gradient"), "Expected value after 'gradient'");
		EXPECT_EQ(parse_error("basis_set sto-3g maxp_diff"), "Expected value after 'MaxP_diff'");
		EXPECT_EQ(parse_error("basis_set sto-3g rmsp_diff"), "Expected value after 'RMSP_diff'");
		EXPECT_EQ(parse_error("basis_set"), "Expected basis set name");
		EXPECT_EQ(parse_error("basis_set sto-3g df_basis"), "Expected a fitting basis name after 'df_basis'");
		EXPECT_EQ(parse_error("basis_set sto-3g guess_basis"), "Expected a basis name after 'guess_basis'");
		EXPECT_EQ(parse_error("basis_set sto-3g save"), "Expected a path after 'save'");
	}

	// without a basis set the parse completes and then refuses to run
	TEST(XcwSettingsTests, MissingBasisSetIsRefusedAfterParsing)
	{
		EXPECT_EQ(parse_error("normal params 177 charge 0"), "Basis set name not specified in settings file! Aborting run!");
		EXPECT_EQ(parse_error(""), "Basis set name not specified in settings file! Aborting run!");
	}

	// `read` takes an optional path: a keyword after it is put back (so `read basis_set`
	// reaches the basis_set handler, which then wants its name), a bare token is the path
	// (so `bogus` is the unknown keyword, not `tensor.bin`), and `read` at the end of the
	// file is fine (the error is then the missing basis, not a read failure)
	TEST(XcwSettingsTests, ReadKeywordTakesAnOptionalPath)
	{
		EXPECT_EQ(parse_error("read basis_set"), "Expected basis set name");
		EXPECT_EQ(parse_error("read Basis_Set"), "Expected basis set name");
		EXPECT_EQ(parse_error("read tensor.bin bogus"), "Unknown keyword 'bogus'");
		EXPECT_EQ(parse_error("params 1 read"), "Basis set name not specified in settings file! Aborting run!");
	}

	// the successful parses run construct() on the P1 fixture, so the only thing to check
	// after the parser is that the run got as far as opening XCW.log with the basis it
	// was given; this is the "does not fall over" kind of check. The parsed settings are
	// private to the XCW and construct() writes nothing else of them into the log, so a
	// wrong preset value or a swapped keyword cannot be seen from here without an accessor
	// in Src/; what these do prove is that every keyword below is accepted, that `read`
	// before a keyword puts it back (a swallowed `basis_set` would leave `sto-3g` as an
	// unknown keyword and throw), and that the basis token reaches the basis handler
	TEST(XcwSettingsTests, SloppySlowConvergenceWithReadBeforeAKeyword)
	{
		const std::string log = construct_on_p1("sloppy slow_conv params 177 read basis_set sto-3g f charge 0 mult 1 rhf start 0 step_size 0.01 end 0.02");
		if (log.empty()) {
			GTEST_SKIP() << "tests/P1_test fixture not found";
		}
		EXPECT_NE(log.find("xcw orbital basis set: sto-3g"), std::string::npos) << log;
	}

	TEST(XcwSettingsTests, TightFastConvergenceStreamedTensorWithoutEnd)
	{
		const std::string log = construct_on_p1("tight fast_conv i_float i_tensor_mb 100 stream safe f2 params 177 basis_set sto-3g charge 0 mult 1 uhf start 0 step_size 0.1");
		if (log.empty()) {
			GTEST_SKIP() << "tests/P1_test fixture not found";
		}
		EXPECT_NE(log.find("xcw orbital basis set: sto-3g"), std::string::npos) << log;
	}

	TEST(XcwSettingsTests, VeryTightWithSavedTensorAndExplicitCriteria)
	{
		const std::string log = construct_on_p1("very_tight normal_conv read tensor.bin save out.bin i_double weighted load_wfn nbo f params 177 basis_set sto-3g "
			"conv 1e-6 diis_damping 1e-3 diis_shift 1e-2 conv_diis 1e-5 gradient 7e-5 maxp_diff 1e-5 rmsp_diff 1e-6 damp 0.5 shift 0.5 max_iter 50 charge 0 mult 1 end 0.01");
		if (log.empty()) {
			GTEST_SKIP() << "tests/P1_test fixture not found";
		}
		EXPECT_NE(log.find("xcw orbital basis set: sto-3g"), std::string::npos) << log;
	}
}
