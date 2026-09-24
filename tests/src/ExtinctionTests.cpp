//The extinction models of Src/core/extinction.h against their published forms and against
//finite differences, and the settings keyword that switches them on through the only door it
//has, the XCW constructor on the P1 fixture.
#include "pch.h"
#include <gtest/gtest.h>

#include "core/cell.h"
#include "core/extinction.h"
#include "core/XCW.h"

namespace NoSpherA2UnitTests
{
	namespace
	{
		std::string ext_test_name()
		{
			return ::testing::UnitTest::GetInstance()->current_test_info()->name();
		}

		//what the XCW constructor throws for a settings text that cannot be parsed; a text that
		//parses would run construct() on an empty cif and take the process with it
		std::string ext_parse_error(const std::string& text)
		{
			const auto p = std::filesystem::temp_directory_path() / ("nosphera2_ext_settings_" + ext_test_name() + ".txt");
			std::ofstream(p) << text;
			options opt;
			opt.xcw_settings_path = p;
			opt.do_XCW = true;
			std::string what;
			try { XCW x(opt); }
			catch (const std::runtime_error& e) { what = e.what(); }
			std::filesystem::remove(p);
			return what;
		}

		//construct an XCW on the P1 fixture in a scratch directory (construct() writes XCW.log
		//there) and hand back the log, lowercased; empty when the fixture is missing
		std::string ext_construct_on_p1(const std::string& text)
		{
			const auto fixture = nos_test_repo_root() / "tests" / "P1_test";
			if (!std::filesystem::exists(fixture / "P1_test_NA2.cif") || !std::filesystem::exists(fixture / "P1_test.hkl"))
				return "";
			const auto dir = std::filesystem::temp_directory_path() / ("nosphera2_ext_p1_" + ext_test_name());
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
			{ XCW x(opt); }
			std::filesystem::current_path(old_cwd);
			std::stringstream log;
			log << std::ifstream(dir / "XCW.log").rdbuf();
			std::filesystem::remove_all(dir);
			std::string out = log.str();
			std::transform(out.begin(), out.end(), out.begin(), [](unsigned char c) { return (char)std::tolower(c); });
			return out;
		}

		// dy/dt by central differences, for checking the analytic derivative the least
		// squares steps on
		double numeric_dydt(const extinction::model m, const double cos2t, const double t)
		{
			const double h = 1e-6;
			return (extinction::correction(m, cos2t, t + h) - extinction::correction(m, cos2t, t - h)) / (2 * h);
		}
	}

	// SHELX is Zachariasen's (1 + t)^(-1/2) on Fc^2
	TEST(ExtinctionTests, ShelxIsZachariasen)
	{
		for (const double t : {0.0, 0.01, 0.5, 3.0, 50.0})
		{
			EXPECT_NEAR(extinction::correction(extinction::model::shelx, 0.3, t),
				1.0 / std::sqrt(1.0 + t), 1e-12);
		}
	}

	// the same geometry constant SHELX's EXTI multiplies: 0.001 lambda^3 / sin(2 theta)
	TEST(ExtinctionTests, GeometryConstantMatchesShelx)
	{
		const double lambda = 0.71073, stl = 0.35;          // sin(theta)/lambda in A^-1
		const double sin_t = lambda * stl;
		const double sin_2t = 2 * sin_t * std::sqrt(1 - sin_t * sin_t);
		EXPECT_NEAR(extinction::geometry_constant(lambda, stl),
			0.001 * std::pow(lambda, 3) / sin_2t, 1e-15);
		EXPECT_NEAR(extinction::cos_2theta(lambda, stl), 1 - 2 * sin_t * sin_t, 1e-15);
		// a reflection this wavelength cannot reach gets no correction rather than a NaN
		EXPECT_EQ(extinction::geometry_constant(lambda, 2.0), 0.0);
	}

	// Becker-Coppens with its mosaic coefficients zeroed is Zachariasen's leading term, so
	// the two models' refined coefficients differ only by the factor 2 SHELX absorbed
	TEST(ExtinctionTests, BeckerCoppensLeadingTerm)
	{
		// A_L = 0.025 + 0.285 cos2t and B_L vanish together at cos2t = -0.025/0.285
		const double cos2t = -0.025 / 0.285;
		const double t = 1e-4;   // A t^2 is then negligible against 2 t whatever B does
		EXPECT_NEAR(extinction::correction(extinction::model::bc_lorentzian, cos2t, t),
			1.0 / std::sqrt(1.0 + 2 * t), 1e-9);
	}

	// the coefficients of Becker & Coppens (1974), as GSAS-II's SCExtinction writes them
	TEST(ExtinctionTests, BeckerCoppensShape)
	{
		for (const double cos2t : {-0.9, -0.3, 0.0, 0.4, 0.95})
		{
			const double AG = 0.58 + 0.48 * cos2t + 0.24 * cos2t * cos2t;
			const double BG = 0.02 - 0.025 * cos2t;
			const double AL = 0.025 + 0.285 * cos2t;
			const double BL = cos2t < 0.0 ? -0.45 * cos2t : 0.15 - 0.2 * (0.75 - cos2t) * (0.75 - cos2t);
			for (const double t : {0.0, 0.05, 1.5, 20.0})
			{
				EXPECT_NEAR(extinction::correction(extinction::model::bc_gaussian, cos2t, t),
					1.0 / std::sqrt(1 + 2 * t + AG * t * t / (1 + BG * t)), 1e-12);
				EXPECT_NEAR(extinction::correction(extinction::model::bc_lorentzian, cos2t, t),
					1.0 / std::sqrt(1 + 2 * t + AL * t * t / (1 + BL * t)), 1e-12);
			}
		}
	}

	// dy/dt is what the Gauss-Newton step on the coefficients is built from, and it has to be
	// non-zero at t = 0 or a refinement started from zero could never leave it
	TEST(ExtinctionTests, DerivativeMatchesTheShape)
	{
		for (const extinction::model m : {extinction::model::shelx, extinction::model::bc_gaussian,
			extinction::model::bc_lorentzian})
		{
			for (const double cos2t : {-0.5, 0.0, 0.7})
			{
				for (const double t : {0.02, 0.4, 5.0})
				{
					double dydt = 0.0;
					extinction::correction(m, cos2t, t, &dydt);
					EXPECT_NEAR(dydt, numeric_dydt(m, cos2t, t), 1e-7);
				}
				double dydt = 0.0;
				extinction::correction(m, cos2t, 0.0, &dydt);
				EXPECT_NEAR(dydt, m == extinction::model::shelx ? -0.5 : -1.0, 1e-12);
			}
			EXPECT_EQ(extinction::correction(extinction::model::none, 0.0, 3.0), 1.0);
		}
	}

	// the azimuth-averaged anisotropic form reduces exactly to the isotropic one for an
	// isotropic tensor, for every direction - that degeneracy is the whole justification for
	// using it in place of the direction cosines an hkl file does not carry
	TEST(ExtinctionTests, AnisotropicReducesToIsotropic)
	{
		const double x = 3.7e-4;
		const std::array<double, 6> X{ x, x, x, 0.0, 0.0, 0.0 };
		for (const std::array<double, 3>& raw : { std::array<double, 3>{1, 0, 0},
			std::array<double, 3>{0, 1, 0}, std::array<double, 3>{0, 0, 1},
			std::array<double, 3>{1, 1, 1}, std::array<double, 3>{-2, 0.5, 3} })
		{
			const double n = std::sqrt(raw[0] * raw[0] + raw[1] * raw[1] + raw[2] * raw[2]);
			const std::array<double, 3> h_unit{ raw[0] / n, raw[1] / n, raw[2] / n };
			std::array<double, 6> a{};
			extinction::aniso_coefficients(h_unit, a);
			double sum = 0.0;
			for (int p = 0; p < 6; p++) sum += a[p] * X[p];
			EXPECT_NEAR(sum, x, 1e-15);
		}
	}

	// and it is the quadratic form it claims to be: (tr X - h.X.h)/2 for a general tensor
	TEST(ExtinctionTests, AnisotropicIsTheAveragedQuadraticForm)
	{
		const std::array<double, 6> X{ 1.0, 2.0, 3.0, 0.4, -0.7, 0.2 };
		const std::array<double, 3> h_unit{ 0.36, -0.48, 0.8 };   // already a unit vector
		std::array<double, 6> a{};
		extinction::aniso_coefficients(h_unit, a);
		double sum = 0.0;
		for (int p = 0; p < 6; p++) sum += a[p] * X[p];
		const double hXh = X[0] * h_unit[0] * h_unit[0] + X[1] * h_unit[1] * h_unit[1] + X[2] * h_unit[2] * h_unit[2]
			+ 2 * (X[3] * h_unit[0] * h_unit[1] + X[4] * h_unit[0] * h_unit[2] + X[5] * h_unit[1] * h_unit[2]);
		EXPECT_NEAR(sum, 0.5 * (X[0] + X[1] + X[2] - hXh), 1e-14);
	}

	// `extinction <model> [iso|aniso] [fixed] [start]` in any order after the model, and a
	// model name is required - the words are what the user types, so they are worth a check
	TEST(ExtinctionTests, SettingsKeywordNamesTheModel)
	{
		EXPECT_EQ(ext_parse_error("extinction"),
			"Expected a model (shelx, bc_gaussian or bc_lorentzian) after 'extinction'");
		EXPECT_EQ(ext_parse_error("extinction aniso 1e-3"),
			"Expected a model (shelx, bc_gaussian or bc_lorentzian) after 'extinction'");
		EXPECT_EQ(ext_parse_error("extinction shelx sideways"),
			"Could not read 'sideways' after 'extinction'");
		EXPECT_EQ(ext_parse_error("wavelength"), "Expected value after 'wavelength'");
	}

	// the whole path: keyword -> settings -> setup_extinction, with the wavelength coming from
	// _diffrn_radiation_wavelength in the fixture's CIF, which is the only place it can
	TEST(ExtinctionTests, IsotropicModelIsSetUpFromTheCifWavelength)
	{
		const std::string log = ext_construct_on_p1(
			"params 177 basis_set sto-3g charge 0 mult 1 rhf start 0 step_size 0.01 end 0.02 extinction shelx");
		if (log.empty()) GTEST_SKIP() << "tests/P1_test fixture not found";
		EXPECT_NE(log.find("xcw extinction: shelx (zachariasen), isotropic (1 parameter)"), std::string::npos) << log;
		EXPECT_NE(log.find("lambda = 0.71073 a"), std::string::npos) << log;
	}

	// the anisotropic tensor, held fixed, with the wavelength given in the settings file
	// instead - and the words in the other order, which the parser allows
	TEST(ExtinctionTests, AnisotropicFixedModelTakesTheSettingsWavelength)
	{
		const std::string log = ext_construct_on_p1(
			"params 177 basis_set sto-3g charge 0 mult 1 rhf start 0 step_size 0.01 end 0.02 "
			"wavelength 1.54178 extinction 5e-4 aniso fixed bc_lorentzian");
		if (log.empty()) GTEST_SKIP() << "tests/P1_test fixture not found";
		EXPECT_NE(log.find("becker-coppens, lorentzian mosaic, anisotropic (azimuth-averaged, 6 parameters), held fixed"), std::string::npos) << log;
		EXPECT_NE(log.find("lambda = 1.54178 a"), std::string::npos) << log;
		EXPECT_NE(log.find("start value 0.0005"), std::string::npos) << log;
	}

	// setup_extinction builds the scattering vector in Cartesian as sum_j rcm(i,j) h_j and
	// then normalises it; if the index convention were transposed the direction would be
	// wrong for anything but an orthogonal cell, and |h| would stop being 1/d
	TEST(ExtinctionTests, CartesianScatteringVectorHasTheRightLength)
	{
		const cell triclinic(7.21, 9.43, 11.07, 83.4, 95.1, 104.8);
		for (const std::array<int, 3>& hkl : { std::array<int, 3>{1, 0, 0},
			std::array<int, 3>{0, 1, 0}, std::array<int, 3>{0, 0, 1},
			std::array<int, 3>{3, -2, 5}, std::array<int, 3>{-1, 4, 2} })
		{
			double h_cart[3] = { 0.0, 0.0, 0.0 };
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++) h_cart[i] += triclinic.get_rcm_angs(i, j) * hkl[j];
			const double length = std::sqrt(h_cart[0] * h_cart[0] + h_cart[1] * h_cart[1] + h_cart[2] * h_cart[2]);
			EXPECT_NEAR(length, 1.0 / triclinic.get_d_of_hkl(hkl), 1e-10);
			EXPECT_NEAR(length, 2.0 * triclinic.get_stl_of_hkl(hkl), 1e-10);
		}
	}
}
