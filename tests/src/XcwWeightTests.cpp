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
}
