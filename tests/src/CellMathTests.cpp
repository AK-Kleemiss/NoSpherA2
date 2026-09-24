#include "pch.h"

#include "core/constants.h"
#include "core/convenience.h"
#include "core/nos_math.h"
#include "core/atoms.h"
#include "core/tsc_block.h"
#include "core/cell.h"
#include "core/wfn_class.h"
#include "core/cif.h"

namespace NoSpherA2UnitTests
{
	// error_check ends the process with exit(-1); Windows reports the full value,
	// POSIX wait() only the low 8 bits.
#ifdef _WIN32
	constexpr unsigned CELLMATH_ERROR_EXIT_CODE = static_cast<unsigned>(-1);
#else
	constexpr unsigned CELLMATH_ERROR_EXIT_CODE = 255u;
#endif

	namespace
	{
		// every temporary file of this suite lives in the system temp directory
		// under a name unique to the test (ctest runs tests in parallel, and
		// Windows cannot remove a file another process still reads) and is
		// removed by the test that wrote it
		std::filesystem::path cellmath_tmp(const std::string& name)
		{
			const std::string test = ::testing::UnitTest::GetInstance()->current_test_info()->name();
			return std::filesystem::temp_directory_path() / ("nosphera2_cellmath_" + test + "_" + name);
		}

		std::filesystem::path write_text(const std::string& name, const std::string& text)
		{
			const std::filesystem::path p = cellmath_tmp(name);
			std::ofstream f(p);
			f << text;
			return p;
		}

		std::string read_text(const std::filesystem::path& p)
		{
			std::ifstream f(p);
			std::stringstream ss;
			ss << f.rdbuf();
			return ss.str();
		}

		// an orthogonal P-1 cell: 8 x 8 x 6 A, volume 384 A^3
		std::string p1bar_cif(const std::string& extra_before_symmetry = "")
		{
			return "data_test\n"
				"_cell_length_a 8.0\n"
				"_cell_length_b 8.0\n"
				"_cell_length_c 6.0\n"
				"_cell_angle_alpha 90.0\n"
				"_cell_angle_beta 90.0\n"
				"_cell_angle_gamma 90.0\n"
				"_cell_volume 384.0\n"
				+ extra_before_symmetry +
				"loop_\n"
				"_space_group_symop_operation_xyz\n"
				"'x, y, z'\n"
				"'-x, -y, -z'\n";
		}

		asym_atom make_asym(const std::string& label, int type, const d3& frac, bool grown = false)
		{
			asym_atom a;
			a.label = label;
			a.type = type;
			a.pos = { 0.0, 0.0, 0.0 };
			a.frac_pos = frac;
			a.asym_fact = 0.0;
			a.anom = cdouble(0.0, 0.0);
			a.grown = grown;
			return a;
		}

		asym_atom make_xyz(const cell& c, const std::string& label, int type, double fx, double fy, double fz)
		{
			const vec cart = c.get_coords_cartesian(fx, fy, fz, true);
			asym_atom a = make_asym(label, type, { 0.0, 0.0, 0.0 });
			a.pos = { cart[0], cart[1], cart[2] };
			return a;
		}
	}

	// the parameter constructor is never used by the readers, so its metric
	// tensors are checked against the defining relations: cm rows are the
	// lattice vectors (in bohr, like read_CIF), cm * rcm equals two pi times
	// the identity, and the identity is the only symmetry operation
	TEST(CellMathTests, SixParameterCtorMatchesAnalyticMetrics)
	{
		const double a = 5.0, b = 7.0, c = 9.0, al = 80.0, be = 95.0, ga = 110.0;
		const cell cl(a, b, c, al, be, ga);
		const double ca = std::cos(al * constants::PI_180);
		const double cb = std::cos(be * constants::PI_180);
		const double cg = std::cos(ga * constants::PI_180);
		const double V = a * b * c * std::sqrt(1 + 2 * ca * cb * cg - ca * ca - cb * cb - cg * cg);
		EXPECT_NEAR(cl.get_V(), V, 1e-10);
		EXPECT_NEAR(cl.get_ca(), ca, 1e-15);
		EXPECT_NEAR(cl.get_cb(), cb, 1e-15);
		EXPECT_NEAR(cl.get_cg(), cg, 1e-15);
		EXPECT_NEAR(cl.get_sa(), std::sin(al * constants::PI_180), 1e-15);
		EXPECT_NEAR(cl.get_sb(), std::sin(be * constants::PI_180), 1e-15);
		EXPECT_NEAR(cl.get_sg(), std::sin(ga * constants::PI_180), 1e-15);
		EXPECT_EQ(cl.get_crystal_system(), "triclinic");
		ASSERT_EQ(cl.get_sym()[0][0].size(), 1u);
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				EXPECT_EQ(cl.get_sym(i, j, 0), i == j ? 1 : 0) << i << j;
		// row lengths are the cell edges, row dot products encode the angles
		const double lengths[3] = { a, b, c };
		for (int i = 0; i < 3; i++)
		{
			double n2 = 0.0;
			for (int j = 0; j < 3; j++)
				n2 += cl.get_cm_angs(i, j) * cl.get_cm_angs(i, j);
			EXPECT_NEAR(std::sqrt(n2), lengths[i], 1e-10) << "row " << i;
		}
		double ab = 0.0, ac = 0.0, bc = 0.0;
		for (int j = 0; j < 3; j++)
		{
			ab += cl.get_cm_angs(0, j) * cl.get_cm_angs(1, j);
			ac += cl.get_cm_angs(0, j) * cl.get_cm_angs(2, j);
			bc += cl.get_cm_angs(1, j) * cl.get_cm_angs(2, j);
		}
		EXPECT_NEAR(ab, a * b * cg, 1e-10);
		EXPECT_NEAR(ac, a * c * cb, 1e-10);
		EXPECT_NEAR(bc, b * c * ca, 1e-10);
		// cm * rcm = 2 pi * I
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
			{
				double s = 0.0;
				for (int k = 0; k < 3; k++)
					s += cl.get_cm(i, k) * cl.get_rcm(k, j);
				EXPECT_NEAR(s, (i == j) ? constants::TWO_PI : 0.0, 1e-10) << i << " " << j;
			}
	}

	// set_system is only reached through the parameter constructor, so all
	// eight branches (both orthorhombic spellings included) are exercised here;
	// note that b = c with a different a lands in the second orthorhombic
	// branch, the classifier only knows the c-unique tetragonal setting
	TEST(CellMathTests, SetSystemClassifiesEveryBranch)
	{
		EXPECT_EQ(cell(5, 5, 5, 90, 90, 90).get_crystal_system(), "cubic");
		EXPECT_EQ(cell(5, 5, 8, 90, 90, 90).get_crystal_system(), "tetragonal");
		EXPECT_EQ(cell(5, 6, 7, 90, 90, 90).get_crystal_system(), "orthorhombic");
		EXPECT_EQ(cell(5, 5, 8, 90, 90, 120).get_crystal_system(), "hexagonal");
		EXPECT_EQ(cell(5, 5, 5, 80, 80, 80).get_crystal_system(), "trigonal");
		EXPECT_EQ(cell(5, 6, 7, 90, 100, 90).get_crystal_system(), "monoclinic");
		EXPECT_EQ(cell(5, 7, 7, 90, 90, 90).get_crystal_system(), "orthorhombic");
		EXPECT_EQ(cell(5, 6, 7, 80, 95, 110).get_crystal_system(), "triclinic");
		EXPECT_EQ(cell().get_crystal_system(), "triclinic");
	}

	// get_d_of_hkl, get_stl_of_hkl and get_reciprocal_metric must agree with
	// the reciprocal lattice vectors stored in rcm, which is an independent route
	TEST(CellMathTests, DSpacingAgreesWithReciprocalVectorsAndMetric)
	{
		const cell cl(5.0, 7.0, 9.0, 80.0, 95.0, 110.0);
		const std::array<int, 3> hkl = { 1, -2, 3 };
		const std::array<double, 3> hkl_d = { 1.0, -2.0, 3.0 };
		// d* = |h a* + k b* + l c*| with the columns of rcm / 2 pi in 1/A (get_rcm_angs) as a*, b*, c*
		double dstar2 = 0.0;
		for (int j = 0; j < 3; j++)
		{
			double comp = 0.0;
			for (int i = 0; i < 3; i++)
				comp += hkl[i] * cl.get_rcm_angs(j, i);
			dstar2 += comp * comp;
		}
		const double d = cl.get_d_of_hkl(hkl);
		EXPECT_NEAR(d, 1.0 / std::sqrt(dstar2), 1e-10);
		EXPECT_NEAR(cl.get_d_of_hkl(hkl_d), d, 1e-12);
		EXPECT_NEAR(cl.get_stl_of_hkl(hkl), 1.0 / (2.0 * d), 1e-12);
		EXPECT_NEAR(cl.get_stl_of_hkl(hkl_d), 0.5 * std::sqrt(dstar2), 1e-10);
		const std::array<double, 6> G = cl.get_reciprocal_metric();
		const double h = hkl[0], k = hkl[1], l = hkl[2];
		const double from_metric = G[0] * h * h + G[1] * k * k + G[2] * l * l + 2 * G[3] * h * k + 2 * G[4] * h * l + 2 * G[5] * k * l;
		EXPECT_NEAR(from_metric, dstar2, 1e-12);
		// cubic sanity value: d(110) of a 4 A cube is 4 / sqrt(2)
		const cell cube(4.0, 4.0, 4.0, 90.0, 90.0, 90.0);
		EXPECT_NEAR(cube.get_d_of_hkl(std::array<int, 3>{ 1, 1, 0 }), 4.0 / std::sqrt(2.0), 1e-12);
		EXPECT_NEAR(cube.get_stl_of_hkl(std::array<int, 3>{ 2, 0, 0 }), 0.25, 1e-12);
	}

	// fractional to cartesian is the linear map given by the rows of cm (bohr),
	// in angstrom or bohr; both the pointer and the vector overload must agree
	TEST(CellMathTests, CartesianCoordinatesFollowTheCellMatrix)
	{
		const cell cl(5.0, 7.0, 9.0, 80.0, 95.0, 110.0);
		const double f[3] = { 0.1, -0.25, 0.7 };
		const vec ang = cl.get_coords_cartesian(f[0], f[1], f[2], false);
		const vec bohr = cl.get_coords_cartesian(f[0], f[1], f[2], true);
		double ptr_bohr[3] = { 0.0, 0.0, 0.0 };
		cl.make_coords_cartesian(ptr_bohr, f[0], f[1], f[2]);
		for (int j = 0; j < 3; j++)
		{
			double expected = 0.0;
			for (int i = 0; i < 3; i++)
				expected += f[i] * cl.get_cm(i, j);
			EXPECT_NEAR(bohr[j], expected, 1e-10) << j;
			EXPECT_NEAR(ang[j], constants::bohr2ang(expected), 1e-10) << j;
			EXPECT_NEAR(ptr_bohr[j], bohr[j], 1e-14) << j;
		}
	}

	// the indexed accessors return the six parameters and fall back to the
	// documented sentinel values (-400 and -1) with a message for a bad index
	TEST(CellMathTests, AngleAndLengthAccessorsWithSentinels)
	{
		const cell cl(5.0, 6.0, 7.0, 91.0, 92.0, 93.0);
		EXPECT_DOUBLE_EQ(cl.get_length(0), 5.0);
		EXPECT_DOUBLE_EQ(cl.get_length(1), 6.0);
		EXPECT_DOUBLE_EQ(cl.get_length(2), 7.0);
		EXPECT_DOUBLE_EQ(cl.get_angle(0), 91.0);
		EXPECT_DOUBLE_EQ(cl.get_angle(1), 92.0);
		EXPECT_DOUBLE_EQ(cl.get_angle(2), 93.0);
		EXPECT_NEAR(cl.get_angle_rad(2), 93.0 * constants::PI_180, 1e-15);
		testing::internal::CaptureStdout();
		EXPECT_DOUBLE_EQ(cl.get_angle(3), -400.0);
		EXPECT_DOUBLE_EQ(cl.get_angle_rad(-1), -400.0);
		EXPECT_DOUBLE_EQ(cl.get_length(7), -1.0);
		const std::string out = testing::internal::GetCapturedStdout();
		EXPECT_NE(out.find("Wrong angle!"), std::string::npos);
		EXPECT_NE(out.find("Wrong length!"), std::string::npos);
	}

	// the CIF reader and the parameter constructor store the same cm (bohr),
	// rcm (2 pi / bohr) and reciprocal edges
	TEST(CellMathIoTests, FileCtorMatchesParameterCtorUpToUnits)
	{
		const std::filesystem::path p = write_text("units.cif",
			"data_t\n_cell_length_a 5.0\n_cell_length_b 7.0\n_cell_length_c 9.0\n"
			"_cell_angle_alpha 80.0\n_cell_angle_beta 95.0\n_cell_angle_gamma 110.0\n"
			"_cell_volume 291.375(5)\n"
			"loop_\n_symmetry_equiv_pos_as_xyz\n'x, y, z'\n");
		std::ostringstream log;
		const cell from_file(p, log, false, false);
		std::filesystem::remove(p);
		const cell from_params(5.0, 7.0, 9.0, 80.0, 95.0, 110.0);
		EXPECT_NEAR(from_file.get_V(), from_params.get_V(), 1e-10);
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
			{
				EXPECT_NEAR(from_file.get_cm(i, j), from_params.get_cm(i, j), 1e-9) << i << j;
				EXPECT_NEAR(from_file.get_rcm(i, j), from_params.get_rcm(i, j), 1e-9) << i << j;
				EXPECT_NEAR(from_file.get_cm_angs(i, j), constants::bohr2ang(from_params.get_cm(i, j)), 1e-9) << i << j;
				EXPECT_NEAR(from_file.get_rcm_angs(i, j), constants::ang2bohr(from_params.get_rcm(i, j)) / constants::TWO_PI, 1e-9) << i << j;
			}
		// reciprocal edges: a* = b c sin(alpha) / V
		const double sa = from_params.get_sa(), sb = from_params.get_sb(), sg = from_params.get_sg(), V = from_params.get_V();
		EXPECT_NEAR(from_file.get_as(), 7.0 * 9.0 * sa / V, 1e-12);
		EXPECT_NEAR(from_file.get_bs(), 5.0 * 9.0 * sb / V, 1e-12);
		EXPECT_NEAR(from_file.get_cs(), 5.0 * 7.0 * sg / V, 1e-12);
		EXPECT_NEAR(from_params.get_as(), from_file.get_as(), 1e-12);
		// the deprecated tag spelling was accepted and the identity is the only operation
		ASSERT_EQ(from_file.get_sym()[0][0].size(), 1u);
		EXPECT_EQ(from_file.get_sym(0, 0, 0), 1);
		EXPECT_NE(log.str().find("Reading:"), std::string::npos);
	}

	// the debug flag walks every diagnostic branch of the reader; the trace has
	// to name each stage so a broken read can be located from the log alone
	TEST(CellMathIoTests, DebugTraceNamesEveryReaderStage)
	{
		const std::filesystem::path p = write_text("debug.cif", p1bar_cif());
		std::ostringstream log;
		const cell cl(p, log, true, false);
		std::filesystem::remove(p);
		const std::string out = log.str();
		for (const char* stage : { "starting to read cif!", "Starting while !.eof()", "Making cm and rcm",
			"line in loop field definition", "Reading operation!", "Comparing", "This is a new symmetry operation!", "RCM done!" })
			EXPECT_NE(out.find(stage), std::string::npos) << stage;
		// outside XCW mode the inversion is folded into the identity: one rotation, two translations
		EXPECT_EQ(cl.get_sym()[0][0].size(), 1u);
		EXPECT_EQ(cl.get_trans()[0].size(), 2u);
		EXPECT_NEAR(cl.get_V(), 384.0, 1e-10);
	}

	// a loop_ that is not the symmetry block (here atom sites) has to be skipped
	// so that the real symmetry loop further down is still found
	TEST(CellMathIoTests, NonSymmetryLoopIsSkipped)
	{
		const std::string atom_loop =
			"loop_\n_atom_site_label\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n"
			"C1 0.1 0.2 0.3\nN1 0.4 0.5 0.6\n\n";
		const std::filesystem::path p = write_text("atomloop.cif", p1bar_cif(atom_loop));
		std::ostringstream log;
		const cell cl(p, log, true, true);
		std::filesystem::remove(p);
		EXPECT_NE(log.str().find("I don't think this is the symmetry block"), std::string::npos);
		ASSERT_EQ(cl.get_sym()[0][0].size(), 2u);
		EXPECT_EQ(cl.get_sym(0, 0, 1), -1);
		EXPECT_EQ(cl.get_sym(1, 1, 1), -1);
		EXPECT_EQ(cl.get_sym(2, 2, 1), -1);
		EXPECT_EQ(cl.get_sym(0, 1, 1), 0);
		EXPECT_EQ(cl.get_sym(1, 0, 0), 0);
		EXPECT_EQ(cl.get_sym(2, 2, 0), 1);
	}

	// a cell parameter that is not a number, a missing cell keyword and a
	// _cell_volume more than 10 % off all have to leave through error_check
	TEST(CellMathIoDeathTest, MalformedCellHeaderExitsCleanly)
	{
		const std::filesystem::path bad_number = write_text("badnum.cif",
			"data_t\n_cell_length_a abc\n_cell_length_b 8.0\n_cell_length_c 6.0\n"
			"_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n");
		const std::filesystem::path missing = write_text("missing.cif",
			"data_t\n_cell_length_a 8.0\n_cell_length_b 8.0\n"
			"_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n");
		const std::filesystem::path volume = write_text("volume.cif",
			"data_t\n_cell_length_a 8.0\n_cell_length_b 8.0\n_cell_length_c 6.0\n"
			"_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n_cell_volume 100.0\n");
		std::ostringstream log;
		EXPECT_EXIT(cell c(bad_number, log), ::testing::ExitedWithCode(CELLMATH_ERROR_EXIT_CODE), ".*");
		EXPECT_EXIT(cell c(missing, log), ::testing::ExitedWithCode(CELLMATH_ERROR_EXIT_CODE), ".*");
		EXPECT_EXIT(cell c(volume, log), ::testing::ExitedWithCode(CELLMATH_ERROR_EXIT_CODE), ".*");
		EXPECT_EXIT(cell c(cellmath_tmp("does_not_exist.cif"), log), ::testing::ExitedWithCode(CELLMATH_ERROR_EXIT_CODE), ".*");
		std::filesystem::remove(bad_number);
		std::filesystem::remove(missing);
		std::filesystem::remove(volume);
	}

	// a symmetry factor that std::stod cannot parse ("q") must be reported
	// through error_check rather than being silently taken as zero
	TEST(CellMathIoDeathTest, UnparsableSymopFactorExitsCleanly)
	{
		int rot[3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
		double trans[3] = { 0.0, 0.0, 0.0 };
		EXPECT_EXIT(cell::parse_symop("q,y,z", "t.cif", rot, trans, std::cout),
			::testing::ExitedWithCode(CELLMATH_ERROR_EXIT_CODE), ".*");
		EXPECT_EXIT(cell::parse_symop("x+1/q,y,z", "t.cif", rot, trans, std::cout),
			::testing::ExitedWithCode(CELLMATH_ERROR_EXIT_CODE), ".*");
	}

	// an xyz atom that is the inversion image of the asymmetric atom is grown,
	// linked to it through the inversion, and both get half the symmetry weight
	TEST(CellMathIoTests, GrownStructureLinksInversionImage)
	{
		const std::filesystem::path p = write_text("grow.cif", p1bar_cif());
		std::ostringstream log;
		cell cl(p, log, false, true);
		std::filesystem::remove(p);
		ASSERT_EQ(cl.get_sym()[0][0].size(), 2u);

		std::vector<asym_atom> asym = { make_asym("C1", 6, { 0.1, 0.2, 0.3 }) };
		std::vector<asym_atom> xyz = {
			make_xyz(cl, "x0", 6, 0.1, 0.2, 0.3),
			make_xyz(cl, "x1", 6, -0.1, -0.2, -0.3) };
		cl.grow_asym_atoms(asym, xyz);
		ASSERT_EQ(asym.size(), 2u);
		EXPECT_FALSE(asym[0].grown);
		EXPECT_TRUE(asym[1].grown);
		EXPECT_NEAR(asym[1].frac_pos[0], -0.1, 1e-10);
		EXPECT_NEAR(asym[1].frac_pos[1], -0.2, 1e-10);
		EXPECT_NEAR(asym[1].frac_pos[2], -0.3, 1e-10);
		// the xyz atom that was already in the asymmetric unit got its fractions too
		EXPECT_NEAR(xyz[0].frac_pos[2], 0.3, 1e-10);

		ivec3 links;
		cl.eval_symm(asym, 1, links);
		ASSERT_EQ(links.size(), 1u);
		ASSERT_EQ(links[0].size(), 2u);
		ASSERT_EQ(links[0][0].size(), 1u);
		ASSERT_EQ(links[0][1].size(), 1u);
		EXPECT_EQ(links[0][0][0], 0);
		EXPECT_EQ(links[0][1][0], 1);

		const ivec applied = cl.apply_grown(links);
		ASSERT_EQ(applied.size(), 1u);
		EXPECT_EQ(applied[0], 1);

		cl.set_symmetry_factors(asym, links);
		EXPECT_NEAR(asym[0].asym_fact, 0.5, 1e-12);
		EXPECT_NEAR(asym[1].asym_fact, 0.5, 1e-12);
		// the image remembers the operation that made it, the parent has none
		EXPECT_EQ(asym[0].sym_op, -1);
		EXPECT_EQ(asym[1].sym_op, 1);
	}

	// when only one of two asymmetric atoms has its inversion image present the
	// operation is not applied and the improperly grown structure is flagged
	TEST(CellMathIoTests, PartiallyGrownStructureWarns)
	{
		const std::filesystem::path p = write_text("partial.cif", p1bar_cif());
		std::ostringstream log;
		cell cl(p, log, false, true);
		std::filesystem::remove(p);

		std::vector<asym_atom> asym = {
			make_asym("C1", 6, { 0.1, 0.2, 0.3 }),
			make_asym("N1", 7, { 0.35, 0.15, 0.05 }) };
		std::vector<asym_atom> xyz = { make_xyz(cl, "x", 6, 0.9, 0.8, 0.7) };
		cl.grow_asym_atoms(asym, xyz);
		ASSERT_EQ(asym.size(), 3u);
		ivec3 links;
		cl.eval_symm(asym, 2, links);
		EXPECT_EQ(links[0][2].size(), 1u);
		EXPECT_TRUE(links[1][2].empty());
		testing::internal::CaptureStderr();
		const ivec applied = cl.apply_grown(links);
		const std::string err = testing::internal::GetCapturedStderr();
		EXPECT_TRUE(applied.empty());
		EXPECT_NE(err.find("Symmetry operation not fully matched"), std::string::npos);
		cl.set_symmetry_factors(asym, links);
		EXPECT_NEAR(asym[0].asym_fact, 0.5, 1e-12);
		EXPECT_NEAR(asym[1].asym_fact, 1.0, 1e-12);
		EXPECT_NEAR(asym[2].asym_fact, 0.5, 1e-12);
		EXPECT_EQ(asym[2].sym_op, 1);
	}

	namespace
	{
		// P21/c with an orthogonal 8 x 8 x 6 A cell: 0 identity, 1 screw, 2 inversion, 3 glide
		std::string p21c_cif()
		{
			return "data_test\n"
				"_cell_length_a 8.0\n"
				"_cell_length_b 8.0\n"
				"_cell_length_c 6.0\n"
				"_cell_angle_alpha 90.0\n"
				"_cell_angle_beta 90.0\n"
				"_cell_angle_gamma 90.0\n"
				"_cell_volume 384.0\n"
				"loop_\n"
				"_space_group_symop_operation_xyz\n"
				"'x, y, z'\n"
				"'-x, y+1/2, -z+1/2'\n"
				"'-x, -y, -z'\n"
				"'x, -y+1/2, z+1/2'\n";
		}

		cell p21c_cell()
		{
			const std::filesystem::path p = write_text("p21c.cif", p21c_cif());
			std::ostringstream log;
			cell cl(p, log, false, true);
			std::filesystem::remove(p);
			return cl;
		}
	}

	// The six tests below specify cell::compose_ops, grown_subgroup,
	// coset_representatives and set_subgroup_factors - the old grown-structure
	// implementation, which is commented out in cell.h, in cell.cpp and at its
	// XCW.cpp call site as "an old but working implementation ... kept for reference
	// purposes in case something goes wrong with the new implementation". These
	// tests were the only live callers left, so every platform's build stopped on
	// them. They are disabled the same way the code they specify is, and belong in
	// whichever commit brings that implementation back.
#if 0
	// composition is matched modulo lattice translations: the screw squared is
	// (x, y+1, z), i.e. the identity, and screw after inversion is the glide
	TEST(CellMathIoTests, ComposeOpsMatchesModuloLattice)
	{
		cell cl = p21c_cell();
		ASSERT_EQ(cl.get_sym()[0][0].size(), 4u);
		EXPECT_EQ(cl.compose_ops(1, 1), 0);
		EXPECT_EQ(cl.compose_ops(2, 2), 0);
		EXPECT_EQ(cl.compose_ops(1, 2), 3);
		EXPECT_EQ(cl.compose_ops(3, 1), 2);
		EXPECT_EQ(cl.compose_ops(0, 3), 3);
	}

	// a molecule grown across the screw axis is one orbit of H = {1, 2_1}: the
	// structure factor sum needs only the identity and inversion cosets and the
	// atoms keep weight 1 because nothing in the cluster is repeated by a coset
	TEST(CellMathIoTests, GrownScrewImageProjectsOntoTwoCosets)
	{
		cell cl = p21c_cell();
		std::vector<asym_atom> asym = { make_asym("C1", 6, { 0.1, 0.2, 0.3 }) };
		std::vector<asym_atom> xyz = {
			make_xyz(cl, "x0", 6, 0.1, 0.2, 0.3),
			make_xyz(cl, "x1", 6, -0.1, 0.7, 0.2) };
		cl.grow_asym_atoms(asym, xyz);
		ASSERT_EQ(asym.size(), 2u);
		ivec3 links;
		cl.eval_symm(asym, 1, links);
		const ivec H = cl.grown_subgroup(links);
		ASSERT_EQ(H, (ivec{ 0, 1 }));
		const ivec reps = cl.coset_representatives(H);
		ASSERT_EQ(reps.size(), 2u);
		EXPECT_EQ(reps[0], 0);
		EXPECT_EQ(reps[1], 2);
		cl.set_symmetry_factors(asym, links);
		EXPECT_NEAR(asym[0].asym_fact, 0.5, 1e-12);
		cl.set_subgroup_factors(asym, links, H);
		EXPECT_NEAR(asym[0].asym_fact, 1.0, 1e-12);
		EXPECT_NEAR(asym[1].asym_fact, 1.0, 1e-12);
	}

	// an atom on the inversion centre with a ligand and its inversion image: H is
	// {1, -1}, the centre's stabiliser lies inside H so it keeps weight 1 while
	// the full-group scheme would have halved it
	TEST(CellMathIoTests, GrownSiteSymmetricClusterKeepsUnitWeights)
	{
		cell cl = p21c_cell();
		std::vector<asym_atom> asym = {
			make_asym("Fe1", 26, { 0.0, 0.0, 0.0 }),
			make_asym("N1", 7, { 0.1, 0.15, 0.2 }) };
		std::vector<asym_atom> xyz = {
			make_xyz(cl, "x0", 26, 0.0, 0.0, 0.0),
			make_xyz(cl, "x1", 7, 0.1, 0.15, 0.2),
			make_xyz(cl, "x2", 7, -0.1, -0.15, -0.2) };
		cl.grow_asym_atoms(asym, xyz);
		ASSERT_EQ(asym.size(), 3u);
		ivec3 links;
		cl.eval_symm(asym, 2, links);
		ASSERT_EQ(links[0][0].size(), 2u);
		const ivec H = cl.grown_subgroup(links);
		ASSERT_EQ(H, (ivec{ 0, 2 }));
		const ivec reps = cl.coset_representatives(H);
		ASSERT_EQ(reps.size(), 2u);
		EXPECT_EQ(reps[0], 0);
		EXPECT_EQ(reps[1], 1);
		cl.set_symmetry_factors(asym, links);
		EXPECT_NEAR(asym[0].asym_fact, 0.5, 1e-12);
		EXPECT_NEAR(asym[1].asym_fact, 0.5, 1e-12);
		cl.set_subgroup_factors(asym, links, H);
		EXPECT_NEAR(asym[0].asym_fact, 1.0, 1e-12);
		EXPECT_NEAR(asym[1].asym_fact, 1.0, 1e-12);
		EXPECT_NEAR(asym[2].asym_fact, 1.0, 1e-12);
	}

	// a whole unit cell of a general-position atom is one orbit of G itself: a
	// single coset, and the image atoms take their parent's weight
	TEST(CellMathIoTests, GrownUnitCellSumsOneCoset)
	{
		cell cl = p21c_cell();
		std::vector<asym_atom> asym = { make_asym("C1", 6, { 0.1, 0.2, 0.3 }) };
		std::vector<asym_atom> xyz = {
			make_xyz(cl, "x0", 6, 0.1, 0.2, 0.3),
			make_xyz(cl, "x1", 6, -0.1, 0.7, 0.2),
			make_xyz(cl, "x2", 6, -0.1, -0.2, -0.3),
			make_xyz(cl, "x3", 6, 0.1, 0.3, 0.8) };
		cl.grow_asym_atoms(asym, xyz);
		ASSERT_EQ(asym.size(), 4u);
		ivec3 links;
		cl.eval_symm(asym, 1, links);
		const ivec H = cl.grown_subgroup(links);
		ASSERT_EQ(H.size(), 4u);
		EXPECT_EQ(cl.coset_representatives(H), (ivec{ 0 }));
		cl.set_subgroup_factors(asym, links, H);
		for (const asym_atom& a : asym) EXPECT_NEAR(a.asym_fact, 1.0, 1e-12);
	}

	// clusters no operation beyond the identity maps onto themselves keep the full
	// sum with the full-group weights: screw and inversion images without the glide
	// image, an image of only one of two asymmetric atoms, and nothing grown at all
	TEST(CellMathIoTests, ClustersWithoutSymmetryKeepTheFullSum)
	{
		cell cl = p21c_cell();
		{
			std::vector<asym_atom> asym = { make_asym("C1", 6, { 0.1, 0.2, 0.3 }) };
			std::vector<asym_atom> xyz = {
				make_xyz(cl, "x0", 6, 0.1, 0.2, 0.3),
				make_xyz(cl, "x1", 6, -0.1, 0.7, 0.2),
				make_xyz(cl, "x2", 6, -0.1, -0.2, -0.3) };
			cl.grow_asym_atoms(asym, xyz);
			ivec3 links;
			cl.eval_symm(asym, 1, links);
			const ivec H = cl.grown_subgroup(links);
			EXPECT_EQ(H, (ivec{ 0 }));
			EXPECT_EQ(cl.coset_representatives(H).size(), 4u);
			cl.set_subgroup_factors(asym, links, H);
			for (const asym_atom& a : asym) EXPECT_NEAR(a.asym_fact, 1.0 / 3.0, 1e-12);
		}
		{
			std::vector<asym_atom> asym = {
				make_asym("C1", 6, { 0.1, 0.2, 0.3 }),
				make_asym("N1", 7, { 0.35, 0.15, 0.05 }) };
			std::vector<asym_atom> xyz = { make_xyz(cl, "x", 6, -0.1, -0.2, -0.3) };
			cl.grow_asym_atoms(asym, xyz);
			ivec3 links;
			cl.eval_symm(asym, 2, links);
			EXPECT_EQ(cl.grown_subgroup(links), (ivec{ 0 }));
		}
		{
			std::vector<asym_atom> asym = { make_asym("C1", 6, { 0.1, 0.2, 0.3 }) };
			ivec3 links;
			cl.eval_symm(asym, 1, links);
			EXPECT_EQ(cl.grown_subgroup(links), (ivec{ 0 }));
			EXPECT_EQ(cl.coset_representatives({ 0 }).size(), 4u);
		}
	}

	// a cluster that is H-invariant without being one H-orbit per atom: the
	// central atom on the inversion centre with both its screw image and the
	// inversion image of a ligand. H = {1, -1} maps it onto itself, the centre and
	// its screw image are two H-orbits and share the weight |H| / (|stab_G| * copies)
	TEST(CellMathIoTests, GrownTwoOrbitClusterSharesTheWeight)
	{
		cell cl = p21c_cell();
		std::vector<asym_atom> asym = {
			make_asym("Fe1", 26, { 0.0, 0.0, 0.0 }),
			make_asym("N1", 7, { 0.1, 0.15, 0.2 }) };
		std::vector<asym_atom> xyz = {
			make_xyz(cl, "x0", 26, 0.0, 0.0, 0.0),
			make_xyz(cl, "x1", 26, 0.0, 0.5, 0.5),
			make_xyz(cl, "x2", 7, 0.1, 0.15, 0.2),
			make_xyz(cl, "x3", 7, -0.1, -0.15, -0.2) };
		cl.grow_asym_atoms(asym, xyz);
		ASSERT_EQ(asym.size(), 4u);
		ivec3 links;
		cl.eval_symm(asym, 2, links);
		const ivec H = cl.grown_subgroup(links);
		ASSERT_EQ(H, (ivec{ 0, 2 }));
		EXPECT_EQ(cl.coset_representatives(H), (ivec{ 0, 1 }));
		cl.set_subgroup_factors(asym, links, H);
		// asymmetric atoms first (Fe1, N1), then the grown images (Fe1 screw image, N1 inversion image)
		EXPECT_NEAR(asym[0].asym_fact, 0.5, 1e-12);
		EXPECT_NEAR(asym[1].asym_fact, 1.0, 1e-12);
		EXPECT_NEAR(asym[2].asym_fact, 0.5, 1e-12);
		EXPECT_NEAR(asym[3].asym_fact, 1.0, 1e-12);
	}
#endif

	// xyz atoms that coincide with an asymmetric atom, also when shifted by a
	// lattice translation, must not be appended a second time
	TEST(CellMathIoTests, LatticeShiftedXyzAtomsAreNotGrown)
	{
		const std::filesystem::path p = write_text("nogrow.cif", p1bar_cif());
		std::ostringstream log;
		cell cl(p, log, false, true);
		std::filesystem::remove(p);
		std::vector<asym_atom> asym = { make_asym("C1", 6, { 0.1, 0.2, 0.3 }) };
		std::vector<asym_atom> xyz = {
			make_xyz(cl, "x0", 6, 1.1, 0.2, -0.7),
			make_xyz(cl, "x1", 6, 0.1, 0.2, 0.3) };
		cl.grow_asym_atoms(asym, xyz);
		EXPECT_EQ(asym.size(), 1u);
		EXPECT_FALSE(xyz[0].grown);
		EXPECT_FALSE(xyz[1].grown);
	}

	// every argument check of the atomID constructors has to throw the
	// documented exception type instead of encoding garbage
	TEST(CellMathTests, AtomIdCtorsRejectInvalidInput)
	{
		EXPECT_THROW(atomID(0.1, 0.2, 0.3, 0, 0), std::out_of_range);
		EXPECT_THROW(atomID(0.1, 0.2, 0.3, 0, 256), std::out_of_range);
		EXPECT_THROW(atomID(0.1, 0.2, 0.3, 40000, 6), std::out_of_range);
		EXPECT_THROW(atomID(0.1, 0.2, 0.3, 0, 6, 300), std::out_of_range);
		EXPECT_THROW(atomID(0.1, 0.2, 0.3, 0, 6, -1), std::out_of_range);
		// Z = 0 in the binary words
		EXPECT_THROW(atomID(std::uint64_t(0), std::uint64_t(0)), std::runtime_error);
		// INT32_MIN cannot come from the encoder, the loader has to refuse it
		const std::array<std::uint64_t, 2> words = atomID(0.1, 0.2, 0.3, 5, 6).as_uint64();
		const std::uint64_t poisoned = (words[0] & 0xFFFFFFFF00000000ULL) | 0x80000000ULL;
		EXPECT_THROW(atomID(poisoned, words[1]), std::runtime_error);
		EXPECT_NO_THROW(atomID(words[0], words[1]));
		EXPECT_THROW(atomID(std::string_view("abc")), std::invalid_argument);
		EXPECT_THROW(atomID(std::string_view("zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz")), std::invalid_argument);
		EXPECT_THROW(atomID(std::string_view("0123456789abcdef0123456789abcde ")), std::invalid_argument);
	}

	// hex string and binary stream are the two persistence forms; both must
	// round trip bit-exactly and the stream form must report a broken stream
	TEST(CellMathTests, AtomIdHexAndStreamRoundTrips)
	{
		const atomID id(0.125, -0.5, 0.75, -7, 26, 3);
		const std::string hex = id.to_hex_string();
		ASSERT_EQ(hex.size(), 32u);
		// hand-derived bit layout, scale = INT32_MAX / 16 = 134217727.9375:
		//   x  0.125 -> llround( 16777215.99) =  16777216 = 0x01000000
		//   y -0.5   -> llround(-67108863.97) = -67108864 = 0xFC000000
		//   z  0.75  -> llround(100663295.95) = 100663296 = 0x06000000
		//   first  = y << 32 | x                      = fc00000001000000
		//   second = res 03 | Z 1a | data fff9 | z    = 031afff906000000
		EXPECT_EQ(hex, "fc00000001000000031afff906000000");
		const atomID from_hex{ std::string_view(hex) };
		EXPECT_TRUE(from_hex == id);
		EXPECT_NEAR(from_hex.frac_x(), 0.125, 1e-8);
		EXPECT_NEAR(from_hex.frac_y(), -0.5, 1e-8);
		EXPECT_NEAR(from_hex.frac_z(), 0.75, 1e-8);
		EXPECT_EQ(from_hex.data(), -7);
		EXPECT_EQ(from_hex.Z(), 26);
		EXPECT_EQ(from_hex.reserved(), 3);

		std::stringstream binary(std::ios::in | std::ios::out | std::ios::binary);
		id.write_atom_id(binary);
		const atomID from_stream(binary);
		EXPECT_TRUE(from_stream == id);
		EXPECT_FALSE(from_stream == atomID(0.125, -0.5, 0.75, -7, 27, 3));

		std::stringstream empty(std::ios::in | std::ios::binary);
		EXPECT_THROW(atomID reader(empty), std::runtime_error);
		EXPECT_THROW(atomID().write_atom_id(binary), std::runtime_error);
		std::stringstream failed(std::ios::out | std::ios::binary);
		failed.setstate(std::ios::badbit);
		EXPECT_THROW(id.write_atom_id(failed), std::runtime_error);
		EXPECT_FALSE(atomID().is_initialized());
		EXPECT_TRUE(id.is_initialized());
	}

	// basis_set_entry equality compares every field, and the atom setters
	// change exactly the addressed primitive; move assignment keeps the data
	TEST(CellMathTests, BasisSetEntryEqualityAndAtomSetters)
	{
		const basis_set_entry e(0.5, 1.25, 1, 0);
		EXPECT_TRUE(e == basis_set_entry(0.5, 1.25, 1, 0));
		EXPECT_FALSE(e == basis_set_entry(0.6, 1.25, 1, 0));
		EXPECT_FALSE(e == basis_set_entry(0.5, 1.5, 1, 0));
		EXPECT_FALSE(e == basis_set_entry(0.5, 1.25, 2, 0));
		EXPECT_FALSE(e == basis_set_entry(0.5, 1.25, 1, 1));

		atom a("C1", atomID(0.1, 0.2, 0.3, 0, 6), 1, 0.0, 0.0, 0.0, 6);
		ASSERT_TRUE(a.push_back_basis_set(10.0, 0.1, 1, 0));
		ASSERT_TRUE(a.push_back_basis_set(2.0, 0.9, 1, 0));
		ASSERT_EQ(a.get_basis_set_size(), 2u);
		a.set_basis_set_exponent(1, 3.0);
		a.set_basis_set_coefficient(0, 0.25);
		EXPECT_DOUBLE_EQ(a.get_basis_set_exponent(0), 10.0);
		EXPECT_DOUBLE_EQ(a.get_basis_set_exponent(1), 3.0);
		EXPECT_DOUBLE_EQ(a.get_basis_set_coefficient(0), 0.25);
		EXPECT_DOUBLE_EQ(a.get_basis_set_coefficient(1), 0.9);
		//operator== also compares the cached primitive, which the setters leave at the constructor values,
		//so a fresh entry is not equal to the edited one; compare the fields instead
		EXPECT_FALSE(a.get_basis_set()[1] == basis_set_entry(0.9, 3.0, 1, 0));
		EXPECT_EQ(a.get_basis_set()[1].get_type(), 1u);
		EXPECT_EQ(a.get_basis_set()[1].get_shell(), 0u);
		EXPECT_DOUBLE_EQ(a.get_basis_set()[1].get_primitive().get_exp(), 2.0);

		atom target;
		target = std::move(a);
		EXPECT_EQ(target.get_label(), "C1");
		EXPECT_EQ(target.get_charge(), 6);
		ASSERT_EQ(target.get_basis_set_size(), 2u);
		EXPECT_DOUBLE_EQ(target.get_basis_set_exponent(1), 3.0);
		atom moved(std::move(target));
		EXPECT_EQ(moved.get_nr(), 1);
		EXPECT_DOUBLE_EQ(moved.get_basis_set_coefficient(0), 0.25);
	}

	// the single-wavefunction CIF writer with a tsc block: spin labels for
	// alpha, beta and unknown MOs and the form factor loop must all be present
	TEST(CellMathIoTests, WriteWfnCifWithTscBlockAndSpinLabels)
	{
		WFN w;
		w.set_basis_set_name("def2-SVP");
		w.push_back_atom("C1", 0.0, 0.0, 0.0, 6);
		w.push_back_atom("H1", 2.0, 0.0, 0.0, 1);
		w.push_back_MO(0, 2.0, -1.0, 0);
		w.push_back_MO(1, 1.0, -0.5, 1);
		w.push_back_MO(2, 0.0, 0.5, 7);
		const cvec2 ff = { { cdouble(1.5, -0.25), cdouble(2.0, 0.5) } };
		const std::vector<atomID> ids = { atomID(0.1, 0.2, 0.3, 0, 6) };
		const ivec2 hkl = { { 1, 2 }, { 0, -1 }, { 3, 0 } };
		tsc_block<int, cdouble> block(ff, ids, hkl);
		options opt;
		const std::filesystem::path p = cellmath_tmp("single.cif");
		write_wfn_CIF(w, p, block, opt);
		const std::string out = read_text(p);
		std::filesystem::remove(p);
		EXPECT_NE(out.find("1 'def2-SVP' ["), std::string::npos);
		EXPECT_NE(out.find("'spins': [alpha beta unknown]"), std::string::npos);
		EXPECT_NE(out.find("'occupancies': [2 1 0]"), std::string::npos);
		EXPECT_NE(out.find("'atom_site_label': 'H1'"), std::string::npos);
		EXPECT_NE(out.find("_aspheric_ffs_partitioning.software 'NoSpherA2'"), std::string::npos);
		EXPECT_NE(out.find("'energies': [-1 -0.5 0.5]"), std::string::npos);
		EXPECT_NE(out.find("_aspheric_ff.index_h"), std::string::npos);
		// hkl is stored per dimension, so the reflections are (1,0,3) and (2,-1,0)
		EXPECT_NE(out.find("1 0 3 '[1.5 ]' '[-0.25 ]'"), std::string::npos);
		EXPECT_NE(out.find("2 -1 0 '[2 ]' '[0.5 ]'"), std::string::npos);
	}

	// the vector overload numbers basis sets and wavefunctions from 1 and
	// leaves the form factor block out when the tsc block is empty
	TEST(CellMathIoTests, WriteWfnCifVectorOverloadNumbersEntries)
	{
		std::vector<WFN> wavy(2);
		wavy[0].set_basis_set_name("b1");
		wavy[1].set_basis_set_name("b2");
		wavy[0].push_back_atom("O1", 0.0, 0.0, 0.0, 8);
		wavy[1].push_back_atom("N1", 0.0, 0.0, 1.0, 7);
		wavy[0].push_back_MO(0, 2.0, -1.0, 0);
		wavy[1].push_back_MO(0, 1.0, -2.0, 1);
		tsc_block<int, cdouble> empty;
		ASSERT_TRUE(empty.is_empty());
		options opt;
		const std::filesystem::path p = cellmath_tmp("vector.cif");
		write_wfn_CIF(wavy, p, empty, opt);
		const std::string out = read_text(p);
		std::filesystem::remove(p);
		EXPECT_NE(out.find("1 'b1' ["), std::string::npos);
		EXPECT_NE(out.find("2 'b2' ["), std::string::npos);
		EXPECT_NE(out.find("1 'Molecular' 'GTO' 'Cartesian' {"), std::string::npos);
		EXPECT_NE(out.find("2 'Molecular' 'GTO' 'Cartesian' {"), std::string::npos);
		EXPECT_NE(out.find("'spins': [alpha]"), std::string::npos);
		EXPECT_NE(out.find("'spins': [beta]"), std::string::npos);
		EXPECT_EQ(out.find("_aspheric_ff"), std::string::npos);
	}

	// both tsc overloads refuse a wavefunction without a basis set name
	TEST(CellMathIoDeathTest, WriteWfnCifWithoutBasisSetExits)
	{
		WFN w;
		w.push_back_atom("C1", 0.0, 0.0, 0.0, 6);
		tsc_block<int, cdouble> empty;
		options opt;
		const std::filesystem::path p = cellmath_tmp("nobasis.cif");
		EXPECT_EXIT(write_wfn_CIF(w, p, empty, opt), ::testing::ExitedWithCode(CELLMATH_ERROR_EXIT_CODE), ".*");
		std::vector<WFN> wavy(1);
		EXPECT_EXIT(write_wfn_CIF(wavy, p, empty, opt), ::testing::ExitedWithCode(CELLMATH_ERROR_EXIT_CODE), ".*");
		std::filesystem::remove(p);
	}

	// first system: an exact fit x = (1, 0, 1) that the outer loop reaches
	// through columns 0 and 2 alone, the gradient of column 1 is then exactly
	// zero (KKT stop). second system: the third column added drives x1 to
	// -3.44, so the Lawson-Hanson inner loop has to step back and drop it;
	// the final active set {0, 2} gives the normal equations
	//   [[5, -4], [-4, 9]] y = [1, 11]  ->  y = (53, 59) / 29,
	// residual r = (-17, -31, 28, 34) / 29, |r| = sqrt(3190) / 29, and the
	// dropped column's gradient c1 . r = -31 / 29 < 0 confirms optimality
	TEST(CellMathTests, NnlsInnerLoopMovesCoefficientBackToZero)
	{
		dMatrix2 A(4, 3);
		const double a_vals[12] = { -1, -1, 0, 1, 1, -1, -1, 1, 2, 2, 0, -1 };
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 3; j++)
				A(i, j) = a_vals[i * 3 + j];
		dMatrix1 b(4);
		b(0) = -1; b(1) = 0; b(2) = 1; b(3) = 1;
		NNLSResult r = nnls(A, b);
		EXPECT_EQ(r.status, 0);
		ASSERT_EQ(r.x.size(), 3u);
		EXPECT_NEAR(r.x[0], 1.0, 1e-10);
		EXPECT_NEAR(r.x[1], 0.0, 1e-10);
		EXPECT_NEAR(r.x[2], 1.0, 1e-10);
		EXPECT_NEAR(r.rnorm, 0.0, 1e-10);

		dMatrix2 C(4, 3);
		const double c_vals[12] = { 2, 2, -2, 0, 1, 2, 0, 0, 1, 1, 1, 0 };
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 3; j++)
				C(i, j) = c_vals[i * 3 + j];
		dMatrix1 d(4);
		d(0) = -1; d(1) = 3; d(2) = 3; d(3) = 3;
		r = nnls(C, d);
		EXPECT_EQ(r.status, 0);
		EXPECT_NEAR(r.x[0], 53.0 / 29.0, 1e-10);
		EXPECT_NEAR(r.x[1], 0.0, 1e-10);
		EXPECT_NEAR(r.x[2], 59.0 / 29.0, 1e-10);
		EXPECT_NEAR(r.rnorm, std::sqrt(3190.0) / 29.0, 1e-10);
	}

	// an exhausted iteration budget returns status 1 with a zero vector; with a
	// duplicated column the first copy fits b exactly and the second one's
	// gradient is zero up to rounding, so it never enters the active set
	// (whether the "not numerically independent" branch fires on the way
	// depends on the last bit of the Householder reflection, so it is not
	// asserted; the answer is (1, 0) either way)
	TEST(CellMathTests, NnlsIterationLimitAndDependentColumn)
	{
		dMatrix2 Id(2, 2);
		Id(0, 0) = 1; Id(0, 1) = 0; Id(1, 0) = 0; Id(1, 1) = 1;
		dMatrix1 b(2);
		b(0) = 1; b(1) = 2;
		testing::internal::CaptureStderr();
		NNLSResult r = nnls(Id, b, 0);
		const std::string err = testing::internal::GetCapturedStderr();
		EXPECT_EQ(r.status, 1);
		EXPECT_NE(err.find("did not converge"), std::string::npos);
		ASSERT_EQ(r.x.size(), 2u);
		EXPECT_DOUBLE_EQ(r.x[0], 0.0);
		EXPECT_DOUBLE_EQ(r.x[1], 0.0);

		dMatrix2 D(2, 2);
		D(0, 0) = 1; D(0, 1) = 1; D(1, 0) = 1; D(1, 1) = 1;
		dMatrix1 c(2);
		c(0) = 1; c(1) = 1;
		r = nnls(D, c);
		EXPECT_EQ(r.status, 0);
		EXPECT_NEAR(r.x[0], 1.0, 1e-12);
		EXPECT_NEAR(r.x[1], 0.0, 1e-12);
		EXPECT_NEAR(r.rnorm, 0.0, 1e-12);
	}

	// dgesv on an exactly singular matrix reports the zero pivot instead of
	// silently returning garbage, and a regular system is solved exactly
	TEST(CellMathTests, SolveLinearSystemReportsSingularMatrix)
	{
		vec A = { 1, 2, 2, 4 };
		vec b = { 1, 2 };
		testing::internal::CaptureStdout();
		solve_linear_system(A, 2, b);
		const std::string out = testing::internal::GetCapturedStdout();
		EXPECT_NE(out.find("LAPACKE_dgesv returned 2"), std::string::npos);

		const vec2 M = { { 2, 1 }, { 1, 3 } };
		vec rhs = { 5, 10 };
		solve_linear_system(M, rhs);
		EXPECT_NEAR(rhs[0], 1.0, 1e-12);
		EXPECT_NEAR(rhs[1], 3.0, 1e-12);
	}

	// the least-squares overload prints the residual norm; a consistent tall
	// system fits exactly and an all-zero column makes dgels report its info
	TEST(CellMathTests, LeastSquaresOverloadPrintsResidualAndRankFailure)
	{
		vec A = { 1, 0, 0, 1, 1, 1 };
		vec b = { 1, 2, 3 };
		testing::internal::CaptureStdout();
		solve_linear_system(A, 3, 2, b);
		std::string out = testing::internal::GetCapturedStdout();
		EXPECT_NEAR(b[0], 1.0, 1e-12);
		EXPECT_NEAR(b[1], 2.0, 1e-12);
		EXPECT_NE(out.find("Error: 0.000000000000"), std::string::npos);
		EXPECT_EQ(out.find("returned"), std::string::npos);

		vec Z = { 1, 0, 1, 0, 1, 0 };
		vec c = { 1, 2, 3 };
		testing::internal::CaptureStdout();
		solve_linear_system(Z, 3, 2, c);
		out = testing::internal::GetCapturedStdout();
		EXPECT_NE(out.find("LAPACKE_dgels returned 2"), std::string::npos);
	}

	// a matrix-vector product with mismatched shapes throws instead of reading
	// out of bounds, and the dot product of an empty vector is defined as zero
	TEST(CellMathTests, SelfDotShapeCheckAndEmptyDot)
	{
		const vec2 M = { { 1, 2 }, { 3, 4 } };
		EXPECT_THROW(self_dot(M, vec{ 1, 2, 3 }), std::invalid_argument);
		EXPECT_THROW(self_dot(M, vec{ 1 }, true), std::invalid_argument);
		const vec ok = self_dot(M, vec{ 1, 1 });
		ASSERT_EQ(ok.size(), 2u);
		EXPECT_DOUBLE_EQ(ok[0], 3.0);
		EXPECT_DOUBLE_EQ(ok[1], 7.0);
		// M is not symmetric, so the transposed product is distinguishable
		const vec okT = self_dot(M, vec{ 1, 1 }, true);
		ASSERT_EQ(okT.size(), 2u);
		EXPECT_DOUBLE_EQ(okT[0], 4.0);
		EXPECT_DOUBLE_EQ(okT[1], 6.0);
		EXPECT_DOUBLE_EQ(dot(vec{}, vec{ 1, 2 }), 0.0);
		EXPECT_DOUBLE_EQ(dot(vec{ 1, 2 }, vec{}), 0.0);
		EXPECT_EQ(dot(cvec{}, cvec{}), cdouble(0.0, 0.0));
		EXPECT_DOUBLE_EQ(dot(vec{ 1, 2, 3 }, vec{ 4, 5, 6 }), 32.0);
	}

	// flatten is instantiated for int and complex nesting as well; the order
	// has to be row-major over the nested vectors
	TEST(CellMathTests, FlattenIntAndComplexNesting)
	{
		const ivec2 i2 = { { 1, 2 }, { 3 } };
		const ivec fi = flatten<int>(i2);
		EXPECT_EQ(fi, (ivec{ 1, 2, 3 }));
		const cvec2 c2 = { { cdouble(1, 2) }, { cdouble(3, 4), cdouble(5, 6) } };
		const cvec fc = flatten<cdouble>(c2);
		ASSERT_EQ(fc.size(), 3u);
		EXPECT_EQ(fc[2], cdouble(5, 6));
		const ivec3 i3 = { { { 1 }, { 2, 3 } }, { { 4 } } };
		EXPECT_EQ(flatten<int>(i3), (ivec{ 1, 2, 3, 4 }));
		const cvec3 c3 = { { { cdouble(0, 1) } }, { { cdouble(2, 0), cdouble(0, 3) } } };
		const cvec fc3 = flatten<cdouble>(c3);
		ASSERT_EQ(fc3.size(), 3u);
		EXPECT_EQ(fc3[0], cdouble(0, 1));
		EXPECT_EQ(fc3[2], cdouble(0, 3));
		EXPECT_TRUE(flatten<int>(ivec2{}).empty());
	}

	// the SVD pseudoinverse has to drop singular values under the cutoff (rank
	// deficient square case) and handle rectangular input with the transposed shape
	TEST(CellMathTests, PseudoInverseDropsSmallSingularValues)
	{
		dMatrix2 R(2, 2);
		R(0, 0) = 1; R(0, 1) = 1; R(1, 0) = 1; R(1, 1) = 1;
		const dMatrix2 Rp = LAPACKE_invert(R);
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				EXPECT_NEAR(Rp(i, j), 0.25, 1e-12) << i << j;

		dMatrix2 T(3, 2);
		T(0, 0) = 2; T(0, 1) = 0; T(1, 0) = 0; T(1, 1) = 4; T(2, 0) = 0; T(2, 1) = 0;
		const dMatrix2 Tp = LAPACKE_invert(T);
		ASSERT_EQ(Tp.extent(0), 2u);
		ASSERT_EQ(Tp.extent(1), 3u);
		// T+ = [[1/2, 0, 0], [0, 1/4, 0]]; with the cutoff above the singular
		// value 2 the first row is zeroed and only the 1/4 survives
		const dMatrix2 Tc = LAPACKE_invert(T, 3.0);
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 3; j++)
			{
				const double diag = (i == j) ? (i == 0 ? 0.5 : 0.25) : 0.0;
				EXPECT_NEAR(Tp(i, j), diag, 1e-12) << i << j;
				EXPECT_NEAR(Tc(i, j), (i == 1 && j == 1) ? 0.25 : 0.0, 1e-12) << i << j;
			}
	}

	// a symmetric row/column swap on a non-square matrix is a programming
	// error and has to leave through error_check
	TEST(CellMathDeathTest, SwapRowsColsRequiresSquareMatrix)
	{
		dMatrix2 M(2, 3);
		for (unsigned i = 0; i < 6; i++)
			M.data()[i] = static_cast<double>(i);
		EXPECT_EXIT(swap_rows_cols_symm(M, 0, 1), ::testing::ExitedWithCode(CELLMATH_ERROR_EXIT_CODE), ".*");
	}

	// matrix square root through the spectral decomposition: the principal
	// root of [[2,1],[1,2]] (eigenpairs 1 and 3 on (1,-1) and (1,1)) is
	// [[1+r, r-1], [r-1, 1+r]] / 2 with r = sqrt(3), which a sign flip of one
	// eigenvalue's root would not give; eigenvalues below the cutoff are dropped
	TEST(CellMathTests, MatSqrtSquaresBackAndAppliesCutoff)
	{
		vec A = { 2, 1, 1, 2 };
		vec W(2);
		const vec S = mat_sqrt(A, W);
		const double r3 = std::sqrt(3.0);
		EXPECT_NEAR(S[0], (1 + r3) / 2, 1e-12);
		EXPECT_NEAR(S[1], (r3 - 1) / 2, 1e-12);
		EXPECT_NEAR(S[2], (r3 - 1) / 2, 1e-12);
		EXPECT_NEAR(S[3], (1 + r3) / 2, 1e-12);
		EXPECT_NEAR(W[0], 1.0, 1e-12);
		EXPECT_NEAR(W[1], r3, 1e-12);
		// with the cutoff above the smaller eigenvalue only the larger survives
		vec B = { 2, 1, 1, 2 };
		vec Wc(2);
		const vec Sc = mat_sqrt(B, Wc, 1.5);
		EXPECT_DOUBLE_EQ(Wc[0], 0.0);
		EXPECT_NEAR(Wc[1], r3, 1e-12);
		// sqrt(3) times the projector onto (1,1)/sqrt(2): every entry sqrt(3)/2
		for (int i = 0; i < 4; i++)
			EXPECT_NEAR(Sc[i], r3 / 2.0, 1e-12) << i;
	}
}
