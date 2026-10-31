//The crystallographic weighting of the XCW (structure_factors::update_cryst_weights, the SHELX weight and
//the port of cctbx's new_shelx_weighting.optimise_parameters) against hand-computed values and a synthetic
//data set whose optimum is known, and the `weight` settings keyword.
#include "pch.h"
#include <gtest/gtest.h>

#include "core/structure_factors.h"

namespace NoSpherA2UnitTests
{
	namespace
	{
		using ws = options::SF_settings::weight_scheme;

		//structure_factors without a model: three reflections, scale^2 = 2, |Fc|^2 = 40 each
		void fill(structure_factors& sf, options& opt)
		{
			sf.opt = &opt;
			sf.model_data.nr = 3;
			sf.model_data.nr_fit = 3;
			sf.model_data.n_params = 1;
			auto& d = sf.scatter_data;
			d.scale = std::sqrt(2.0);
			//Fo^2 = 100, 10 and -10 (the sign rides on F_obs, F_obs2 is |Fo^2|), sigma(Fo^2) = 5, 4, 4
			d.F_obs = { 10.0, std::sqrt(10.0), -std::sqrt(10.0) };
			d.abs_F_obs = { 10.0, std::sqrt(10.0), std::sqrt(10.0) };
			d.F_obs2 = { 100.0, 10.0, 10.0 };
			d.sigma_obs2 = { 5.0, 4.0, 4.0 };
			d.sigma_obs = { 0.25, 4.0 * 0.5 / std::sqrt(10.0), 4.0 * 0.5 / std::sqrt(10.0) };
			d.F_calc = cvec(3, std::complex<double>(std::sqrt(40.0), 0.0));
			d.hkl_mask = { 1, 1, 1 };
		}

		options::XCW_settings parsed_weight(const std::string& text)
		{
			const auto p = std::filesystem::temp_directory_path() / "nosphera2_xcw_weight_settings.txt";
			std::ofstream(p) << text;
			options opt;
			opt.xcw_settings_path = p;
			opt.loadXCWsettings();
			std::filesystem::remove(p);
			return opt.xcw_settings;
		}
	}

	//w = 1/(sigma^2 + (aP)^2 + b K P), P = (max(Fo^2, 0) + 2 K Fc^2)/3, by hand
	TEST(XcwWeightTests, ShelxWeightByHand)
	{
		//P = (100 + 2*2*40)/3 = 260/3
		EXPECT_DOUBLE_EQ(shelx_weight(100.0, 5.0, 40.0, 2.0, 0.1, 1.0), 1.0 / (25.0 + 0.01 * (260.0 / 3) * (260.0 / 3) + 2.0 * 260.0 / 3));
		//negative Fo^2 enters P as 0: P = 160/3
		const double P = 160.0 / 3.0;
		EXPECT_DOUBLE_EQ(shelx_weight(-10.0, 5.0, 40.0, 2.0, 0.1, 1.0), 1.0 / (25.0 + 0.01 * P * P + 2.0 * P));
		//a = b = 0 is 1/sigma^2
		EXPECT_DOUBLE_EQ(shelx_weight(100.0, 5.0, 40.0, 2.0, 0.0, 0.0), 1.0 / 25.0);
	}

	TEST(XcwWeightTests, UpdateCrystWeightsShelx)
	{
		options opt;
		opt.xcw_settings.weight = ws::shelx;
		opt.xcw_settings.weight_params = { 0.1, 1.0 };
		structure_factors sf;
		fill(sf, opt);
		sf.update_cryst_weights();
		ASSERT_EQ(sf.sig_I_.size(), 3u);
		//sigma(Fo^2) = 1/sqrt(w); against F sigma(F) = sigma(Fo^2)/(2|Fo|)
		const double P0 = 260.0 / 3.0, P1 = 170.0 / 3.0, P2 = 160.0 / 3.0;
		const double s0 = std::sqrt(25.0 + 0.01 * P0 * P0 + 2.0 * P0), s1 = std::sqrt(16.0 + 0.01 * P1 * P1 + 2.0 * P1), s2 = std::sqrt(16.0 + 0.01 * P2 * P2 + 2.0 * P2);
		EXPECT_NEAR(sf.sig_I_[0], s0, 1e-12 * s0);
		EXPECT_NEAR(sf.sig_I_[1], s1, 1e-12 * s1);
		EXPECT_NEAR(sf.sig_I_[2], s2, 1e-12 * s2);
		EXPECT_NEAR(sf.sig_F_[0], s0 / 20.0, 1e-12 * s0);
		EXPECT_NEAR(sf.sig_F_[2], s2 * 0.5 / std::sqrt(10.0), 1e-12 * s2);
		EXPECT_EQ(sf.weight_ab_, (vec{ 0.1, 1.0 }));
	}

	//the default copies the read sigmas, bit for bit
	TEST(XcwWeightTests, UpdateCrystWeightsSigmaIsTheReadSigma)
	{
		options opt;
		structure_factors sf;
		fill(sf, opt);
		sf.update_cryst_weights();
		EXPECT_EQ(sf.sig_I_, sf.scatter_data.sigma_obs2);
		EXPECT_EQ(sf.sig_F_, sf.scatter_data.sigma_obs);
	}

	//sigma^2 chosen as Delta^2 - (a0 P)^2 - b0 K P, so w Delta^2 = 1 for every reflection at (a0, b0): every bin's
	//GooF is 1 (n = N // (N - 2) = 1) and calcres is 0 there, its minimum
	TEST(XcwWeightTests, OptimiserFindsTheKnownOptimum)
	{
		const double a0 = 0.05, b0 = 0.5, K = 1.5;
		const int N = 2000;
		vec fo(N), sig(N), fc(N);
		for (int i = 0; i < N; i++) {
			fc[i] = 20.0 * std::pow(100.0, static_cast<double>(i) / (N - 1));
			fo[i] = K * fc[i] * (i % 2 ? 1.3 : 0.7);
			const double P = (fo[i] + 2.0 * K * fc[i]) / 3.0, d = fo[i] - K * fc[i];
			const double s2 = d * d - a0 * a0 * P * P - b0 * K * P;
			ASSERT_GT(s2, 0.0);
			sig[i] = std::sqrt(s2);
		}
		const vec ab = optimise_shelx_weights(fo, sig, fc, K, 2);
		ASSERT_EQ(ab.size(), 2u);
		EXPECT_NEAR(ab[0], a0, 1e-3);
		EXPECT_NEAR(ab[1], b0, 1e-3);
	}

	TEST(XcwWeightTests, SettingsKeyword)
	{
		const options::XCW_settings s = parsed_weight("basis_set sto-3g params 10 weight shelx 0.1 0");
		EXPECT_EQ(s.weight, ws::shelx);
		EXPECT_EQ(s.weight_params, (vec{ 0.1, 0.0 }));
		const options::XCW_settings o = parsed_weight("basis_set sto-3g params 10 weight shelx_opt 0.05 1.5");
		EXPECT_EQ(o.weight, ws::shelx_opt);
		EXPECT_EQ(o.weight_params, (vec{ 0.05, 1.5 }));
		EXPECT_EQ(parsed_weight("basis_set sto-3g params 10").weight, ws::sigma);
		EXPECT_THROW(parsed_weight("basis_set sto-3g params 10 weight shelx -1 0"), std::runtime_error);
	}

	TEST(XcwWeightTests, SettingsKeywordMoreSchemes)
	{
		const std::string head = "basis_set sto-3g params 10 weight ";
		options::XCW_settings s = parsed_weight(head + "unit");
		EXPECT_EQ(s.weight, ws::unit);
		EXPECT_TRUE(s.weight_params.empty());
		EXPECT_EQ(s.robust_c, 0.0);
		s = parsed_weight(head + "stl 1");
		EXPECT_EQ(s.weight, ws::stl);
		EXPECT_EQ(s.weight_params, (vec{ 1.0 }));
		//a variable count of A_r, stopping at the next keyword
		s = parsed_weight("weight Chebyshev 1 0.5 -0.2 basis_set sto-3g params 10");
		EXPECT_EQ(s.weight, ws::chebyshev);
		EXPECT_EQ(s.weight_params, (vec{ 1.0, 0.5, -0.2 }));
		s = parsed_weight(head + "chebyshev_opt 3");
		EXPECT_EQ(s.weight, ws::chebyshev_opt);
		EXPECT_EQ(s.weight_params, (vec{ 3.0 }));
		//robust: c after the scheme's values, 6 when left out; for chebyshev the last value is c
		s = parsed_weight(head + "robust shelx 0.1 0 4");
		EXPECT_EQ(s.weight, ws::shelx);
		EXPECT_EQ(s.weight_params, (vec{ 0.1, 0.0 }));
		EXPECT_EQ(s.robust_c, 4.0);
		s = parsed_weight(head + "robust sigma");
		EXPECT_EQ(s.weight, ws::sigma);
		EXPECT_EQ(s.robust_c, 6.0);
		s = parsed_weight(head + "robust chebyshev 1 0.5 6");
		EXPECT_EQ(s.weight_params, (vec{ 1.0, 0.5 }));
		EXPECT_EQ(s.robust_c, 6.0);
		for (const char* bad : { "stl", "unit 1", "chebyshev", "chebyshev_opt 0", "chebyshev_opt 1.5", "robust chebyshev 1", "robust sigma -1", "robust sigma 0", "tukey" })
			EXPECT_THROW(parsed_weight(head + bad), std::runtime_error) << bad;
	}

	//w = 1 and w = (sin theta/lambda)^a in the refined metric, sigma(F) and sigma(F^2) linked by 2|Fo|
	TEST(XcwWeightTests, UnitAndStlByHand)
	{
		for (const int against : { 1, 2 }) {
			options opt;
			opt.xcw_settings.weight = ws::unit;
			structure_factors sf;
			fill(sf, opt);
			sf.quality_criteria.refine_against = against;
			sf.update_cryst_weights();
			for (int r = 0; r < 3; r++) {
				const double fo = sf.scatter_data.abs_F_obs[r];
				EXPECT_DOUBLE_EQ(sf.sig_I_[r], against == 2 ? 1.0 : 2.0 * fo);
				EXPECT_DOUBLE_EQ(sf.sig_F_[r], against == 2 ? 0.5 / fo : 1.0);
			}
		}
		options opt;
		opt.xcw_settings.weight = ws::stl;
		opt.xcw_settings.weight_params = { 2.0 };
		structure_factors sf;
		fill(sf, opt);
		sf.quality_criteria.refine_against = 2;
		sf.unit_cell = cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
		//in set order (0 2 0), (1 0 0), (1 1 0): sin theta/lambda = 1/(2d) = 0.1, 0.05, sqrt(2)/20; w = stl^2
		sf.scatter_data.hkl = { i3{ 1, 0, 0 }, i3{ 0, 2, 0 }, i3{ 1, 1, 0 } };
		sf.update_cryst_weights();
		const double stl[3] = { 0.1, 0.05, std::sqrt(2.0) / 20.0 };
		for (int r = 0; r < 3; r++) EXPECT_NEAR(sf.sig_I_[r], 1.0 / stl[r], 1e-12 / stl[r]);
	}

	//w = 1/sum A_r T_r(t), t = 2 Fc/max(Fc) - 1: |Fc| = 2, 4, 8 is t = -0.5, 0, 1
	TEST(XcwWeightTests, ChebyshevByHand)
	{
		options opt;
		opt.xcw_settings.weight = ws::chebyshev;
		opt.xcw_settings.weight_params = { 1.0, 0.5, 0.2 };
		structure_factors sf;
		fill(sf, opt);
		sf.quality_criteria.refine_against = 2;
		sf.scatter_data.F_calc = { 2.0, 4.0, 8.0 };
		sf.update_cryst_weights();
		//1 - 0.25 + 0.2 (2 0.25 - 1), 1 - 0.2, 1 + 0.5 + 0.2
		const double sum[3] = { 0.65, 0.8, 1.7 };
		for (int r = 0; r < 3; r++) EXPECT_NEAR(sf.sig_I_[r], std::sqrt(sum[r]), 1e-12);
		//a sum <= 0 is weight 0, sigma infinite
		EXPECT_EQ(chebyshev_weight(-1.0, { 1.0, 2.0 }), 0.0);
		//t = 0.9: T2 = 0.62, T3 = 0.216
		EXPECT_NEAR(chebyshev_weight(0.9, { 1.0, 2.0, 3.0, 4.0 }), 1.0 / (1.0 + 2.0 * 0.9 + 3.0 * 0.62 + 4.0 * 0.216), 1e-14);
	}

	//Delta^2 planted as sum A_r T_r(t): the least-squares fit gives A back
	TEST(XcwWeightTests, ChebyshevFitRecoversCoefficients)
	{
		const vec A = { 3.0, -1.0, 0.5, 0.25 };
		vec t(200), d2(200);
		for (int i = 0; i < 200; i++) {
			t[i] = -1.0 + 2.0 * i / 199.0;
			const double T[4] = { 1.0, t[i], 2.0 * t[i] * t[i] - 1.0, 4.0 * t[i] * t[i] * t[i] - 3.0 * t[i] };
			d2[i] = 0.0;
			for (int r = 0; r < 4; r++) d2[i] += A[r] * T[r];
		}
		const vec fit = fit_chebyshev(t, d2, 4);
		ASSERT_EQ(fit.size(), 4u);
		for (int r = 0; r < 4; r++) EXPECT_NEAR(fit[r], A[r], 1e-6);

		//through structure_factors: Fo^2 = k^2 Fc^2 + sqrt(sum A_r T_r(t)) on 50 reflections, chebyshev_opt 3 at its first fill
		options opt;
		opt.xcw_settings.weight = ws::chebyshev_opt;
		opt.xcw_settings.weight_params = { 3.0 };
		structure_factors sf;
		sf.opt = &opt;
		sf.quality_criteria.refine_against = 2;
		const int N = 50;
		sf.model_data.nr = sf.model_data.nr_fit = N;
		auto& d = sf.scatter_data;
		d.scale = std::sqrt(2.0);
		const vec B = { 4.0, 1.0, 0.5 };
		for (int r = 0; r < N; r++) {
			const double fc = r + 1.0, tt = 2.0 * fc / N - 1.0;
			const double fo2 = 2.0 * fc * fc + std::sqrt(B[0] + B[1] * tt + B[2] * (2.0 * tt * tt - 1.0));
			d.F_calc.push_back(fc);
			d.F_obs2.push_back(fo2);
			d.F_obs.push_back(std::sqrt(fo2));
			d.abs_F_obs.push_back(std::sqrt(fo2));
			d.sigma_obs2.push_back(1.0);
			d.sigma_obs.push_back(0.5 / std::sqrt(fo2));
			d.hkl_mask.push_back(1);
		}
		sf.update_cryst_weights();
		ASSERT_EQ(sf.weight_ab_.size(), 3u);
		for (int r = 0; r < 3; r++) EXPECT_NEAR(sf.weight_ab_[r], B[r], 1e-6);
	}

	//Tukey's biweight on 1/sigma^2: the planted outliers get w = 0, every other reflection w = f/sigma^2 with
	//f = (1 - (z/(c s))^2)^2, z = Delta/sigma, s = median|z|/0.6745
	TEST(XcwWeightTests, RobustZeroesPlantedOutliers)
	{
		for (const int against : { 1, 2 }) {
			options opt;
			opt.xcw_settings.robust_c = 6.0;
			structure_factors sf;
			sf.opt = &opt;
			sf.quality_criteria.refine_against = against;
			const int N = 41;
			sf.model_data.nr = sf.model_data.nr_fit = N;
			auto& d = sf.scatter_data;
			d.scale = std::sqrt(2.0);
			for (int r = 0; r < N; r++) {
				const double fo2 = 80.0 + (r == 3 || r == 17 ? 500.0 : 0.6 * ((r * 7) % 11 - 5));
				d.F_calc.push_back(std::sqrt(40.0));
				d.F_obs2.push_back(fo2);
				d.F_obs.push_back(std::sqrt(fo2));
				d.abs_F_obs.push_back(std::sqrt(fo2));
				d.sigma_obs2.push_back(2.0);
				d.sigma_obs.push_back(1.0 / std::sqrt(fo2));
				d.hkl_mask.push_back(1);
			}
			sf.update_cryst_weights();
			const vec& s = against == 2 ? d.sigma_obs2 : d.sigma_obs;
			const vec& sig = against == 2 ? sf.sig_I_ : sf.sig_F_;
			vec z(N), a(N);
			for (int r = 0; r < N; r++) {
				z[r] = (against == 2 ? 80.0 - d.F_obs2[r] : std::sqrt(80.0) - d.F_obs[r]) / s[r];
				a[r] = std::abs(z[r]);
			}
			std::sort(a.begin(), a.end());
			const double cs = 6.0 * a[N / 2] / 0.6745;
			for (int r = 0; r < N; r++) {
				const double w = 1.0 / (sig[r] * sig[r]);
				if (r == 3 || r == 17) {
					EXPECT_EQ(w, 0.0) << r;
					continue;
				}
				const double u = z[r] / cs, f = (1.0 - u * u) * (1.0 - u * u);
				ASSERT_LT(std::abs(u), 1.0);
				EXPECT_NEAR(w, f / (s[r] * s[r]), 1e-12 / (s[r] * s[r])) << r;
			}
		}
		//even count: the median is the mean of the middle two; a zero scale leaves every factor 1
		const vec f = tukey_biweight({ 1.0, -2.0, 3.0, 100.0 }, { 1, 1, 1, 1 }, 1.0);
		EXPECT_EQ(f[3], 0.0);
		EXPECT_NEAR(f[0], std::pow(1.0 - std::pow(0.6745 / 2.5, 2), 2), 1e-14);
		EXPECT_EQ(tukey_biweight({ 0.0, 0.0, 5.0 }, { 1, 1, 1 }, 6.0), (vec{ 1.0, 1.0, 1.0 }));
	}
}
