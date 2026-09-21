
#include "pch.h"

#include "core/convenience.h"
#include "core/constants.h"
#include "core/fchk.h"
#include "core/AtomGrid.h"
#include "core/SALTED_utilities.h"
#include "core/scattering_factors.h"
#include "core/nos_math.h"
#include "core/GridManager.h"
#include "core/atoms.h"
#include "core/tsc_block.h"
#include "core/cell.h"
#include "core/wfn_class.h"
#include "core/properties.h"
#include "core/integrator.h"
#include "core/integration_params.h"
#include "core/basis_set.h"
#include "core/geometry_aid.h"
#include "core/crystal_energies.h"
#include "core/NoSpherA2.h"
#include "core/isosurface.h"
#include "core/npy.h"
#include "core/libCintMain.h"
#include "core/throughput.h"
#include "core/i_tensor_stream.h"
#include "core/tsc_label_converter.h"
#include "core/sphere_lebedev_rule.h"
#include <occ/qm/hf.h>
#include <occ/qm/scf.h>
#include <spdlog/spdlog.h>
#undef I
#ifdef NOSPHERA2_USE_GPU
#include "core/blas_gpu.h"
#include "core/aux_density_gpu.h"
#endif

static constexpr double PI_VAL = 3.14159265358979323846;


namespace NoSpherA2UnitTests
{
	// ---------------------------------------------------------------------
	// find_frontier_orbitals - the HOMO/LUMO detection behind the Fukui
	// functions. There is no stored frontier index in WFN, so this logic is
	// the only thing standing between a correct f+/f- and a silently wrong
	// one, which makes it worth testing directly rather than only through the
	// integration test.
	// ---------------------------------------------------------------------

	// A closed-shell case with orbital energies present: the HOMO is the
	// highest-energy occupied orbital and the LUMO the lowest-energy virtual.
	TEST(FukuiTests, FrontierOrbitalsRestrictedWithEnergies)
	{
		WFN wavy;
		wavy.push_back_MO(0, 2.0, -10.0);
		wavy.push_back_MO(1, 2.0, -5.0);
		wavy.push_back_MO(2, 2.0, -1.0);  // HOMO
		wavy.push_back_MO(3, 0.0, 0.5);   // LUMO
		wavy.push_back_MO(4, 0.0, 1.5);

		int homo = -1, lumo = -1;
		bool unrestricted = true;
		EXPECT_TRUE(find_frontier_orbitals(wavy, homo, lumo, unrestricted));
		EXPECT_EQ(homo, 2);
		EXPECT_EQ(lumo, 3);
		EXPECT_FALSE(unrestricted);
	}

	// Energies are not always stored - a plain .wfn carries none, and every
	// energy then reads back as exactly 0.0. In that case the routine must fall
	// back to orbital ORDER (last occupied / first virtual) instead of picking
	// an arbitrary orbital because all the energies compare equal.
	TEST(FukuiTests, FrontierOrbitalsFallsBackToOrderWithoutEnergies)
	{
		WFN wavy;
		wavy.push_back_MO(0, 2.0, 0.0);
		wavy.push_back_MO(1, 2.0, 0.0);
		wavy.push_back_MO(2, 2.0, 0.0);  // HOMO by position
		wavy.push_back_MO(3, 0.0, 0.0);  // LUMO by position
		wavy.push_back_MO(4, 0.0, 0.0);

		int homo = -1, lumo = -1;
		bool unrestricted = false;
		EXPECT_TRUE(find_frontier_orbitals(wavy, homo, lumo, unrestricted));
		EXPECT_EQ(homo, 2);
		EXPECT_EQ(lumo, 3);
	}

	// Orbitals are not guaranteed to arrive sorted by energy. When energies are
	// available they must win over index order, otherwise the frontier pair is
	// silently wrong for any reader that emits an unsorted set.
	TEST(FukuiTests, FrontierOrbitalsPreferEnergyOverIndexOrder)
	{
		WFN wavy;
		wavy.push_back_MO(0, 2.0, -1.0);  // highest occupied energy -> HOMO
		wavy.push_back_MO(1, 2.0, -10.0);
		wavy.push_back_MO(2, 0.0, 2.0);
		wavy.push_back_MO(3, 0.0, 0.5);   // lowest virtual energy -> LUMO

		int homo = -1, lumo = -1;
		bool unrestricted = false;
		EXPECT_TRUE(find_frontier_orbitals(wavy, homo, lumo, unrestricted));
		EXPECT_EQ(homo, 0);
		EXPECT_EQ(lumo, 3);
	}

	// The failure this guards is the one that actually bites: a wavefunction
	// file storing only the occupied orbitals (the common case for .wfn and for
	// the sucrose.wfx in tests/) has no LUMO at all, so no Fukui function can be
	// formed. That must be reported, not silently turned into an all-zero f+.
	TEST(FukuiTests, FrontierOrbitalsFailWhenNoVirtualsExist)
	{
		WFN wavy;
		wavy.push_back_MO(0, 2.0, -10.0);
		wavy.push_back_MO(1, 2.0, -5.0);
		wavy.push_back_MO(2, 2.0, -1.0);

		int homo = -1, lumo = -1;
		bool unrestricted = false;
		EXPECT_FALSE(find_frontier_orbitals(wavy, homo, lumo, unrestricted));
		EXPECT_EQ(homo, 2);
		EXPECT_EQ(lumo, -1);
	}

	// The mirror case: no occupied orbitals at all.
	TEST(FukuiTests, FrontierOrbitalsFailWhenNoOccupiedExist)
	{
		WFN wavy;
		wavy.push_back_MO(0, 0.0, 1.0);
		wavy.push_back_MO(1, 0.0, 2.0);

		int homo = -1, lumo = -1;
		bool unrestricted = false;
		EXPECT_FALSE(find_frontier_orbitals(wavy, homo, lumo, unrestricted));
		EXPECT_EQ(homo, -1);
		EXPECT_EQ(lumo, 0);
	}

	// An empty wavefunction must fail rather than index into nothing.
	TEST(FukuiTests, FrontierOrbitalsFailOnEmptyWavefunction)
	{
		WFN wavy;
		int homo = -1, lumo = -1;
		bool unrestricted = false;
		EXPECT_FALSE(find_frontier_orbitals(wavy, homo, lumo, unrestricted));
		EXPECT_EQ(homo, -1);
		EXPECT_EQ(lumo, -1);
	}

	// Fractional occupations from correlated methods must count as occupied. A
	// naive "occ == 2.0 means occupied" test would classify a natural orbital
	// with occupation 1.98 as virtual and put the HOMO in the wrong place.
	TEST(FukuiTests, FrontierOrbitalsTreatFractionalOccupationAsOccupied)
	{
		WFN wavy;
		wavy.push_back_MO(0, 2.0, -10.0);
		wavy.push_back_MO(1, 1.98, -2.0);
		wavy.push_back_MO(2, 0.02, -1.0);  // still occupied, so this is the HOMO
		wavy.push_back_MO(3, 0.0, 0.5);    // genuinely empty -> LUMO

		int homo = -1, lumo = -1;
		bool unrestricted = false;
		EXPECT_TRUE(find_frontier_orbitals(wavy, homo, lumo, unrestricted));
		EXPECT_EQ(homo, 2);
		EXPECT_EQ(lumo, 3);
	}

	//Reference values from the per-point primitive-pair implementation (before 18 Sep 2026) on the 0.2 A cube
//orca_vpot (ORCA 6.1.1, same input) gives -0.0782518, 0.3776471, 0.3006136, 0.0397176 at these points
	TEST(EspTests, PairTableMatchesReferenceCube)
	{
		const auto input = nos_test_repo_root() / "tests" / "epoxide_gbw" / "epoxide.gbw";
		if (!std::filesystem::exists(input)) GTEST_SKIP() << "Missing " << input;
		WFN wave(input, false);
		const WFN::ESP_pairs pairs = wave.build_ESP_pairs();
		EXPECT_EQ(pairs.weight.size() + 1, pairs.off.size());
		const std::array<std::pair<d3, double>, 4> reference = { {
			{ { 1.710101, 11.196696, 1.517632 }, -0.0782518 },  // the negative lobe minimum
			{ { -0.557569, 14.220256, 1.517632 }, 0.377647 },
			{ { -2.825239, 14.976146, 2.651467 }, 0.300611 },
			{ { 2.843936, 16.109981, 3.407357 }, 0.039717 } } };
		for (const auto& [pos, esp] : reference)
			EXPECT_NEAR(wave.computeESP(pos, pairs), esp, 1E-5); // the cube header rounds the grid positions to 1E-6 bohr
	}

	//The analytic ELI-D gradient (orbital Hessians) against a central difference of computeRhoELI
	TEST(EliTests, AnalyticGradientMatchesFiniteDifference)
	{
		const auto input = nos_test_repo_root() / "tests" / "epoxide_gbw" / "epoxide.gbw";
		if (!std::filesystem::exists(input)) GTEST_SKIP() << "Missing " << input;
		WFN wave(input, false);
		const d3 points[] = { { 1.710101, 11.196696, 1.517632 }, { -0.557569, 14.220256, 1.517632 },
			{ -2.825239, 14.976146, 2.651467 }, { 2.843936, 16.109981, 3.407357 }, { 0.3, 13.1, 2.2 } };
		const double h = 1E-4;
		for (const d3& p : points)
		{
			double eli, rho, eli_ref;
			d3 grad;
			wave.computeELIGrad(p, eli, grad);
			wave.computeRhoELI(p, rho, eli_ref);
			EXPECT_NEAR(eli, eli_ref, 1E-10 * std::max(1.0, std::abs(eli_ref)));
			for (int k = 0; k < 3; k++)
			{
				d3 a = p, b = p;
				a[k] += h; b[k] -= h;
				double ea, eb;
				wave.computeRhoELI(a, rho, ea);
				wave.computeRhoELI(b, rho, eb);
				const double fd = (ea - eb) / (2 * h);
				EXPECT_NEAR(grad[k], fd, 1E-6 * std::max(1.0, std::abs(fd))) << "axis " << k << " at " << p[0] << " " << p[1] << " " << p[2];
			}
		}
	}

	//The surface must sit around the molecule it belongs to: readxyzMinMax_fromWFN used to guess the unit
	//from the shortest interatomic distance and scaled the grid by 1.89 whenever no bond was shorter
	//than 2 bohr, so the heavy atoms alone (shortest C-O 1.42 A) are the case that failed (18 Sep 2026)
	TEST(IsosurfaceTests, HirshfeldSurfaceCentredOnMolecule)
	{
		const auto dir = nos_test_repo_root() / "tests" / "isosurface";
		if (!std::filesystem::exists(dir / "pack.xyz")) GTEST_SKIP() << "Missing " << dir;
		svec heavy;
		std::ifstream asu(dir / "asu.xyz");
		for (std::string line; std::getline(asu, line);)
			if (heavy.size() < 2 || (!line.empty() && line[0] != 'H')) heavy.push_back(line);
		heavy[0] = std::to_string(heavy.size() - 2);
		const auto noH = std::filesystem::temp_directory_path() / "asu_noH.xyz";
		{
			std::ofstream out(noH);
			for (const std::string& line : heavy) out << line << '\n';
		}
		WFN mol(noH, false), env(dir / "pack.xyz", false);
		std::filesystem::remove(noH);
		EXPECT_EQ(mol.get_ncen(), 23);
		properties_options opts;
		opts.resolution = 0.4;
		std::ostringstream log;
		std::vector<Triangle> triangles = Hirshfeld_surface(mol, env, opts, log);
		ASSERT_GT(triangles.size(), 1000);
		d3 centre{ 0, 0, 0 }, molecule{ 0, 0, 0 };
		for (const Triangle& t : triangles)
			for (int k = 0; k < 3; k++) centre[k] += t.calc_center()[k] / triangles.size();
		for (int a = 0; a < mol.get_ncen(); a++)
			for (int k = 0; k < 3; k++) molecule[k] += mol.get_atom_coordinate(a, k) / mol.get_ncen();
		for (int k = 0; k < 3; k++)
			EXPECT_NEAR(centre[k], molecule[k], 0.5) << "axis " << k;
		//d_norm: negative where the fragments touch closer than the vdW radii, and a point on an atom's vdW sphere with the environment far away scores about -1 + (d_e - r_e) / r_e
		double lo = 1E9, hi = -1E9;
		for (const Triangle& t : triangles)
		{
			const double v = calc_d_norm_term(t.calc_center(), mol) + calc_d_norm_term(t.calc_center(), env);
			lo = std::min(lo, v);
			hi = std::max(hi, v);
		}
		EXPECT_LT(lo, 0.0);
		EXPECT_GT(hi, 0.0);
		{
			std::ofstream out(noH);
			out << "1\n\nO 0 0 0\n";
		}
		WFN one(noH, false);
		std::filesystem::remove(noH);
		const double r = constants::ang2bohr(constants::vdW_radii[8]);
		EXPECT_NEAR(calc_d_norm_term({ 2 * r, 0, 0 }, one), 1.0, 1E-9);
		EXPECT_NEAR(calc_d_norm_term({ 0, r, 0 }, one), 0.0, 1E-9);
		EXPECT_NEAR(calc_d_norm_term({ 0, 0, r / 2 }, one), -0.5, 1E-9);
	}

	//Red at the oxygen, blue over the CH2 groups: the ESP sign is what a chemist expects
	//A sphere of radius R cut from w = 0.5 + (R^2 - r^2) (decreasing outward like the Hirshfeld weight) is convex everywhere:
	//shape index +1 and curvedness 2/pi ln(1/R) with R in Angstrom; the weight field has no cusp at the grid nodes
	TEST(IsosurfaceTests, CurvatureOfASphereIsConvex)
	{
		const double R = 2.0, h = 0.1;
		const int n = 61;
		cube field({ n, n, n }, 0, true);
		for (int k = 0; k < 3; k++) {
			field.set_origin(k, -3.0);
			field.set_vector(k, k, h);
		}
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				for (int k = 0; k < n; k++) {
					const d3 p = field.get_pos(i, j, k);
					field.set_value(i, j, k, 0.5 + R * R - (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]));
				}
		const std::vector<Triangle> triangles = marchingCubes(field, 0.5);
		ASSERT_GT(triangles.size(), 1000);
		vec S, C;
		surface_curvature(triangles, field, S, C);
		const double expected_C = 2.0 / constants::PI * std::log(1.0 / constants::bohr2ang(R));
		for (size_t t = 0; t < triangles.size(); t++) {
			EXPECT_NEAR(S[t], 1.0, 1E-6) << "face " << t;
			EXPECT_NEAR(C[t], expected_C, 0.04) << "face " << t;
		}
	}

	TEST(IsosurfaceTests, EspColourOfRhoIsosurface)
	{
		const auto input = nos_test_repo_root() / "tests" / "epoxide_gbw" / "epoxide.gbw";
		if (!std::filesystem::exists(input)) GTEST_SKIP() << "Missing " << input;
		WFN wave(input, false);
		properties_options opts;
		opts.resolution = 0.2;
		opts.radius = 2.5;
		readxyzMinMax_fromWFN(wave, opts);
		cube rho(opts.NbSteps, wave.get_ncen(), true);
		for (int i = 0; i < 3; i++)
		{
			rho.set_origin(i, opts.MinMax[i]);
			rho.set_vector(i, i, (opts.MinMax[3 + i] - opts.MinMax[i]) / opts.NbSteps[i]);
		}
		std::ostringstream log;
		Calc_Rho(rho, wave, opts.radius, log, false);
		std::vector<Triangle> triangles = marchingCubes(rho, 0.002);
		ASSERT_GT(triangles.size(), 1000);
		colour_by_ESP(triangles, wave, log);
		const d3 O{ wave.get_atom_coordinate(0, 0), wave.get_atom_coordinate(0, 1), wave.get_atom_coordinate(0, 2) };
		double d_min = 1E9;
		RGB at_oxygen{ 0, 0, 0 };
		for (const Triangle& t : triangles)
		{
			const d3 c = t.calc_center();
			const double d = std::hypot(c[0] - O[0], c[1] - O[1], c[2] - O[2]);
			if (d < d_min) d_min = d, at_oxygen = t.get_colour();
		}
		EXPECT_EQ(at_oxygen[0], 255);
		EXPECT_LT(at_oxygen[2], 128) << "the surface above the oxygen must be red";
		EXPECT_NE(log.str().find("ESP on the surface from -0.06"), std::string::npos) << log.str();
	}

} // namespace NoSpherA2UnitTests
