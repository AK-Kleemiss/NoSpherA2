//Coverage tests for the parts of Src/core/convenience.cpp that ConvenienceTests.cpp leaves untouched:
//the CIF ADP / U_iso readers, the interactive cube selection, the int Fortran record, the option
//branches that die or finish, the small helpers with uncovered branches and the process helpers.
//Every expectation is derived by hand from the input, never from what the code prints.
#include "pch.h"

#include "core/convenience.h"
#include "core/constants.h"
#include "core/basis_set.h"
#include "core/wfn_class.h"
#include "core/cube.h"
#include "core/cell.h"
#include "core/throughput.h"

#include <array>
#include <cfloat>
#include <cstdint>
#include <cstring>
#include <iomanip>

//defined with external linkage in convenience.cpp, not declared in the header
double get_bessel_ratio(const double nu, const double x);

namespace
{
	//a scratch directory named after the test, removed when the test passed; enter() makes it the cwd
	//and the destructor restores the previous one
	struct ScratchDir
	{
		std::filesystem::path dir;
		std::filesystem::path old_cwd;
		bool entered = false;
		explicit ScratchDir(const std::string& name)
			: dir(std::filesystem::temp_directory_path() / ("nos_convcov_" + name))
		{
			std::error_code ec;
			std::filesystem::remove_all(dir, ec);
			std::filesystem::create_directories(dir);
		}
		void enter()
		{
			old_cwd = std::filesystem::current_path();
			std::filesystem::current_path(dir);
			entered = true;
		}
		std::filesystem::path file(const std::string& name) const { return dir / name; }
		~ScratchDir()
		{
			std::error_code ec;
			if (entered)
				std::filesystem::current_path(old_cwd, ec);
			if (!::testing::Test::HasFailure())
				std::filesystem::remove_all(dir, ec);
		}
	};

	void write_text(const std::filesystem::path& path, const std::string& text)
	{
		std::ofstream f(path);
		f << text;
	}

	struct CoutCapture
	{
		std::ostringstream buf;
		std::streambuf* old;
		CoutCapture() : old(std::cout.rdbuf(buf.rdbuf())) {}
		~CoutCapture() { std::cout.rdbuf(old); }
		std::string str() const { return buf.str(); }
	};

	struct CinFeed
	{
		std::istringstream in;
		std::streambuf* old;
		explicit CinFeed(const std::string& text) : in(text), old(std::cin.rdbuf(in.rdbuf())) {}
		~CinFeed()
		{
			std::cin.rdbuf(old);
			std::cin.clear();
		}
	};

	options parse(const std::vector<std::string>& args, const bool debug = false)
	{
		options opt;
		opt.debug = debug;
		for (const auto& a : args)
			opt.arguments.push_back(a);
		opt.digest_options();
		return opt;
	}

	//options holds a stream reference and cannot be assigned, so the capture wraps the construction
	options parse_quiet(const std::vector<std::string>& args, const bool debug, std::string& out)
	{
		CoutCapture cap;
		options opt = parse(args, debug);
		out = cap.str();
		return opt;
	}

	//C1 at (1,2,3) A, O1 at (5,5,5) A, H1 at (7,7,7) A, the fractions 0.1/0.2/0.3, 0.5, 0.7 of a 10 A cube
	WFN three_atom_wfn()
	{
		WFN w(e_origin::NOT_YET_DEFINED);
		w.push_back_atom("A", constants::ang2bohr(1.0), constants::ang2bohr(2.0), constants::ang2bohr(3.0), 6);
		w.push_back_atom("B", constants::ang2bohr(5.0), constants::ang2bohr(5.0), constants::ang2bohr(5.0), 8);
		w.push_back_atom("C", constants::ang2bohr(7.0), constants::ang2bohr(7.0), constants::ang2bohr(7.0), 1);
		return w;
	}

	//the atom_site loop lists C1, O1 and N9 (no WFN atom at 0.9), the aniso loop has its U columns in the order
	//11,22,33,23,13,12, the Gram-Charlier loops give C1 the values 1..10 and 11..25 in column order
	const char* const cif_text =
		"data_cov\n"
		"_cell_length_a 10.0\n"
		"_cell_length_b 10.0\n"
		"_cell_length_c 10.0\n"
		"\n"
		"loop_\n"
		"_atom_site_label\n"
		"_atom_site_fract_x\n"
		"_atom_site_fract_y\n"
		"_atom_site_fract_z\n"
		"_atom_site_U_iso_or_equiv\n"
		"C1 0.1 0.2 0.3 0.03125\n"
		"O1 0.5 0.5 0.5 ?\n"
		"N9 0.9 0.9 0.9 0.5\n"
		"\n"
		"loop_\n"
		"_atom_site_aniso_label\n"
		"_atom_site_aniso_U_11\n"
		"_atom_site_aniso_U_22\n"
		"_atom_site_aniso_U_33\n"
		"_atom_site_aniso_U_23\n"
		"_atom_site_aniso_U_13\n"
		"_atom_site_aniso_U_12\n"
		"C1 0.0625 0.125 0.1875 0.25 0.3125 0.375\n"
		"O1 0.5 0.5 0.5 0.0 0.0 0.0\n"
		"\n"
		"loop_\n"
		"_atom_site_anharm_GC_C_label\n"
		"_atom_site_anharm_GC_C_111\n"
		"_atom_site_anharm_GC_C_112\n"
		"_atom_site_anharm_GC_C_113\n"
		"_atom_site_anharm_GC_C_122\n"
		"_atom_site_anharm_GC_C_123\n"
		"_atom_site_anharm_GC_C_133\n"
		"_atom_site_anharm_GC_C_222\n"
		"_atom_site_anharm_GC_C_223\n"
		"_atom_site_anharm_GC_C_233\n"
		"_atom_site_anharm_GC_C_333\n"
		"C1 1 2 3 4 5 6 7 8 9 10\n"
		"\n"
		"loop_\n"
		"_atom_site_anharm_GC_D_label\n"
		"_atom_site_anharm_GC_D_1111\n"
		"_atom_site_anharm_GC_D_1112\n"
		"_atom_site_anharm_GC_D_1113\n"
		"_atom_site_anharm_GC_D_1122\n"
		"_atom_site_anharm_GC_D_1123\n"
		"_atom_site_anharm_GC_D_1133\n"
		"_atom_site_anharm_GC_D_1222\n"
		"_atom_site_anharm_GC_D_1223\n"
		"_atom_site_anharm_GC_D_1233\n"
		"_atom_site_anharm_GC_D_1333\n"
		"_atom_site_anharm_GC_D_2222\n"
		"_atom_site_anharm_GC_D_2223\n"
		"_atom_site_anharm_GC_D_2233\n"
		"_atom_site_anharm_GC_D_2333\n"
		"_atom_site_anharm_GC_D_3333\n"
		"C1 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25\n"
		"\n";

	//U_11, U_22, U_33, U_12, U_13, U_23 of C1 in the storage order of the readers
	const vec c1_U = { 0.0625, 0.125, 0.1875, 0.375, 0.3125, 0.25 };

	//the aniso rows of the grown-mode tests: one asymmetric atom C1 with U in 11,22,33,12,13,23 order
	const char* const grown_cif_text =
		"data_grown\n"
		"loop_\n"
		"_atom_site_aniso_label\n"
		"_atom_site_aniso_U_11\n"
		"_atom_site_aniso_U_22\n"
		"_atom_site_aniso_U_33\n"
		"_atom_site_aniso_U_12\n"
		"_atom_site_aniso_U_13\n"
		"_atom_site_aniso_U_23\n"
		"C1 0.0625 0.125 0.1875 0.375 0.3125 0.25\n"
		"\n"
		"loop_\n"
		"_atom_site_anharm_GC_C_label\n"
		"_atom_site_anharm_GC_C_111\n"
		"_atom_site_anharm_GC_C_112\n"
		"_atom_site_anharm_GC_C_113\n"
		"_atom_site_anharm_GC_C_122\n"
		"_atom_site_anharm_GC_C_123\n"
		"_atom_site_anharm_GC_C_133\n"
		"_atom_site_anharm_GC_C_222\n"
		"_atom_site_anharm_GC_C_223\n"
		"_atom_site_anharm_GC_C_233\n"
		"_atom_site_anharm_GC_C_333\n"
		"C1 1 2 3 4 5 6 7 8 9 10\n"
		"\n";

	WFN grown_pair_wfn()
	{
		WFN w(e_origin::NOT_YET_DEFINED);
		w.push_back_atom("C1", 0.0, 0.0, 0.0, 6);
		w.push_back_atom("C1_sym", 2.0, 0.0, 0.0, 6);
		return w;
	}

	std::string hex_digest(const uint8_t hash[32])
	{
		std::ostringstream s;
		s << std::hex << std::setfill('0');
		for (int i = 0; i < 32; i++)
			s << std::setw(2) << static_cast<int>(hash[i]);
		return s.str();
	}

	bool valid_occ_dir(const std::filesystem::path& p)
	{
		return std::filesystem::is_directory(p) && std::filesystem::is_directory(p / "basis") &&
			std::filesystem::is_directory(p / "methods");
	}
}

//---------------------------------------------------------------- CIF readers

TEST(ConvenienceCoverageCifTests, AdpReaderWithCellAssignsLabelsFracsAndUij)
{
	ScratchDir scratch("adp_cell");
	write_text(scratch.file("a.cif"), cif_text);
	std::ofstream log(scratch.file("log.txt"));
	WFN w = three_atom_wfn();
	cell unit_cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);

	EXPECT_TRUE(read_fracs_ADPs_from_CIF(scratch.file("a.cif"), w, unit_cell, log, true));

	EXPECT_EQ(w.get_atom_label(0), "C1");
	EXPECT_EQ(w.get_atom_label(1), "O1");
	EXPECT_EQ(w.get_atom_label(2), "C"); //no CIF atom at (7,7,7) A
	EXPECT_DOUBLE_EQ(w.get_atom(0).get_frac_coordinate(0), 0.1);
	EXPECT_DOUBLE_EQ(w.get_atom(0).get_frac_coordinate(1), 0.2);
	EXPECT_DOUBLE_EQ(w.get_atom(0).get_frac_coordinate(2), 0.3);
	EXPECT_DOUBLE_EQ(w.get_atom(1).get_frac_coordinate(2), 0.5);

	const vec2 c1 = w.get_atom(0).get_ADPs();
	ASSERT_EQ(c1.size(), 3u);
	ASSERT_EQ(c1[0].size(), 6u);
	for (int j = 0; j < 6; j++)
		EXPECT_DOUBLE_EQ(c1[0][j], c1_U[j]) << j;
	ASSERT_EQ(c1[1].size(), 10u);
	ASSERT_EQ(c1[2].size(), 15u);
	for (int j = 0; j < 6; j++)
	{
		EXPECT_DOUBLE_EQ(c1[1][j], 1.0 + j) << j;
		EXPECT_DOUBLE_EQ(c1[2][j], 11.0 + j) << j;
	}

	const vec2 o1 = w.get_atom(1).get_ADPs();
	ASSERT_EQ(o1.size(), 3u);
	ASSERT_EQ(o1[0].size(), 6u);
	EXPECT_DOUBLE_EQ(o1[0][0], 0.5);
	EXPECT_DOUBLE_EQ(o1[0][5], 0.0);
	EXPECT_TRUE(o1[1].empty());
	EXPECT_TRUE(o1[2].empty());

	const vec2 h1 = w.get_atom(2).get_ADPs();
	ASSERT_EQ(h1.size(), 3u);
	EXPECT_TRUE(h1[0].empty());
	EXPECT_TRUE(h1[1].empty());
	EXPECT_TRUE(h1[2].empty());
	log.close();
	std::ifstream logged(scratch.file("log.txt"));
	std::string text((std::istreambuf_iterator<char>(logged)), std::istreambuf_iterator<char>());
	EXPECT_NE(text.find("I DID NOT FIND THIS ATOM"), std::string::npos); //N9
}

TEST(ConvenienceCoverageCifTests, AdpReaderWithCellCopiesAllTenCijk)
{
	ScratchDir scratch("adp_cell_cijk");
	write_text(scratch.file("a.cif"), cif_text);
	std::ofstream log(scratch.file("log.txt"));
	WFN w = three_atom_wfn();
	cell unit_cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
	ASSERT_TRUE(read_fracs_ADPs_from_CIF(scratch.file("a.cif"), w, unit_cell, log, false));
	const vec2 c1 = w.get_atom(0).get_ADPs();
	ASSERT_EQ(c1[1].size(), 10u);
	for (int j = 6; j < 10; j++)
		EXPECT_DOUBLE_EQ(c1[1][j], 1.0 + j) << j;
}

TEST(ConvenienceCoverageCifTests, AdpReaderWithCellCopiesAllFifteenDijkl)
{
	ScratchDir scratch("adp_cell_dijkl");
	write_text(scratch.file("a.cif"), cif_text);
	std::ofstream log(scratch.file("log.txt"));
	WFN w = three_atom_wfn();
	cell unit_cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
	ASSERT_TRUE(read_fracs_ADPs_from_CIF(scratch.file("a.cif"), w, unit_cell, log, false));
	const vec2 c1 = w.get_atom(0).get_ADPs();
	ASSERT_EQ(c1[2].size(), 15u);
	for (int j = 6; j < 15; j++)
		EXPECT_DOUBLE_EQ(c1[2][j], 11.0 + j) << j;
}

TEST(ConvenienceCoverageCifTests, LabelReaderReadsUijCijkDijklByLabel)
{
	ScratchDir scratch("adp_label");
	//the label reader gathers its seven tokens across line breaks, so C1's aniso row is wrapped here
	std::string wrapped = cif_text;
	const std::string one_line = "C1 0.0625 0.125 0.1875 0.25 0.3125 0.375\n";
	const size_t at = wrapped.find(one_line);
	ASSERT_NE(at, std::string::npos);
	wrapped.replace(at, one_line.size(), "C1 0.0625 0.125\n0.1875 0.25 0.3125 0.375\n");
	write_text(scratch.file("a.cif"), wrapped);
	std::ofstream log(scratch.file("log.txt"));
	WFN w(e_origin::NOT_YET_DEFINED);
	w.push_back_atom("C1", 0.0, 0.0, 0.0, 6);
	w.push_back_atom("O1", 3.0, 0.0, 0.0, 8);
	w.push_back_atom("H1", 5.0, 0.0, 0.0, 1);
	const ivec3 no_links;

	EXPECT_TRUE(read_fracs_ADPs_from_CIF(scratch.file("a.cif"), w, log, true, false, no_links));

	const vec2 c1 = w.get_atom(0).get_ADPs();
	ASSERT_EQ(c1.size(), 3u);
	ASSERT_EQ(c1[0].size(), 6u);
	for (int j = 0; j < 6; j++)
		EXPECT_DOUBLE_EQ(c1[0][j], c1_U[j]) << j; //columns 23,13,12 land in slots 5,4,3
	ASSERT_EQ(c1[1].size(), 10u);
	for (int j = 0; j < 10; j++)
		EXPECT_DOUBLE_EQ(c1[1][j], 1.0 + j) << j;
	ASSERT_EQ(c1[2].size(), 15u);
	for (int j = 0; j < 15; j++)
		EXPECT_DOUBLE_EQ(c1[2][j], 11.0 + j) << j;

	const vec2 o1 = w.get_atom(1).get_ADPs();
	ASSERT_EQ(o1.size(), 3u);
	ASSERT_EQ(o1[0].size(), 6u);
	EXPECT_DOUBLE_EQ(o1[0][0], 0.5);
	EXPECT_DOUBLE_EQ(o1[0][2], 0.5);
	EXPECT_DOUBLE_EQ(o1[0][3], 0.0);
	EXPECT_TRUE(o1[1].empty());
	EXPECT_TRUE(o1[2].empty());
	EXPECT_TRUE(w.get_atom(2).get_ADPs().empty());
}

TEST(ConvenienceCoverageCifTests, LabelReaderThrowsForUnknownLabel)
{
	ScratchDir scratch("adp_label_unknown");
	write_text(scratch.file("a.cif"), grown_cif_text);
	std::ofstream log(scratch.file("log.txt"));
	WFN w(e_origin::NOT_YET_DEFINED);
	w.push_back_atom("Zz9", 0.0, 0.0, 0.0, 6);
	const ivec3 no_links;
	EXPECT_THROW(read_fracs_ADPs_from_CIF(scratch.file("a.cif"), w, log, false, false, no_links), std::runtime_error);
}

TEST(ConvenienceCoverageCifTests, LabelReaderGrownCopiesUijToLinkedAtoms)
{
	ScratchDir scratch("adp_grown");
	write_text(scratch.file("a.cif"), grown_cif_text);
	std::ofstream log(scratch.file("log.txt"));
	WFN w = grown_pair_wfn();
	//one asymmetric atom; total atom 0 is itself (no operation), total atom 1 is its image under symop 0
	ivec3 links(1);
	links[0].resize(2);
	links[0][1] = { 0 };

	EXPECT_TRUE(read_fracs_ADPs_from_CIF(scratch.file("a.cif"), w, log, false, true, links));

	for (int a = 0; a < 2; a++)
	{
		const vec2 adps = w.get_atom(a).get_ADPs();
		ASSERT_EQ(adps.size(), 3u) << a;
		ASSERT_EQ(adps[0].size(), 6u) << a;
		for (int j = 0; j < 6; j++)
			EXPECT_DOUBLE_EQ(adps[0][j], c1_U[j]) << a << " " << j;
	}
}

TEST(ConvenienceCoverageCifTests, LabelReaderGrownCopiesCijkToLinkedAtoms)
{
	ScratchDir scratch("adp_grown_cijk");
	write_text(scratch.file("a.cif"), grown_cif_text);
	std::ofstream log(scratch.file("log.txt"));
	WFN w = grown_pair_wfn();
	ivec3 links(1);
	links[0].resize(2);
	links[0][0] = { 0 };
	links[0][1] = { 0 };

	ASSERT_TRUE(read_fracs_ADPs_from_CIF(scratch.file("a.cif"), w, log, false, true, links));

	const vec2 image = w.get_atom(1).get_ADPs();
	ASSERT_EQ(image.size(), 3u);
	ASSERT_EQ(image[1].size(), 10u);
	for (int j = 0; j < 10; j++)
		EXPECT_DOUBLE_EQ(image[1][j], 1.0 + j) << j;
}

TEST(ConvenienceCoverageCifTests, UisoReaderReadsValuesAndTreatsPlaceholderAsZero)
{
	ScratchDir scratch("uiso");
	write_text(scratch.file("a.cif"), cif_text);
	std::ofstream log(scratch.file("log.txt"));
	WFN w = three_atom_wfn();
	cell unit_cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);

	const vec u = read_U_iso_from_CIF(scratch.file("a.cif"), w, unit_cell, log, true);

	ASSERT_EQ(u.size(), 3u);
	EXPECT_DOUBLE_EQ(u[0], 0.03125); //exact in float
	EXPECT_DOUBLE_EQ(u[1], 0.0);     //'?'
	EXPECT_DOUBLE_EQ(u[2], 0.0);     //no CIF atom for H1
	EXPECT_EQ(w.get_atom_label(0), "C1");
	EXPECT_EQ(w.get_atom_label(1), "O1");
	EXPECT_EQ(w.get_atom_label(2), "C");
	EXPECT_DOUBLE_EQ(w.get_atom(1).get_frac_coordinate(0), 0.5);
	log.close();
	std::ifstream logged(scratch.file("log.txt"));
	std::string text((std::istreambuf_iterator<char>(logged)), std::istreambuf_iterator<char>());
	EXPECT_NE(text.find("Found an atom: C1"), std::string::npos);
	EXPECT_NE(text.find("I DID NOT FIND THIS ATOM"), std::string::npos);
}

TEST(ConvenienceCoverageCifTests, UisoReaderWithoutUisoColumnGivesZeros)
{
	ScratchDir scratch("uiso_nocol");
	write_text(scratch.file("a.cif"),
		"data_x\n"
		"loop_\n"
		"_atom_site_label\n"
		"_atom_site_fract_x\n"
		"_atom_site_fract_y\n"
		"_atom_site_fract_z\n"
		"C1 0.1 0.2 0.3\n"
		"\n");
	std::ofstream log(scratch.file("log.txt"));
	WFN w = three_atom_wfn();
	cell unit_cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
	const vec u = read_U_iso_from_CIF(scratch.file("a.cif"), w, unit_cell, log, false);
	ASSERT_EQ(u.size(), 3u);
	EXPECT_DOUBLE_EQ(u[0], 0.0);
	EXPECT_DOUBLE_EQ(u[1], 0.0);
	EXPECT_DOUBLE_EQ(u[2], 0.0);
	EXPECT_EQ(w.get_atom_label(0), "C1");
	EXPECT_EQ(w.get_atom_label(1), "B");
}

//---------------------------------------------------------------- cube selection and unsaved files

TEST(ConvenienceCoverageCubeTests, SelectCubesRejectsBadInputAndStoresPairs)
{
	ScratchDir scratch("select_cubes");
	std::vector<WFN> wavy;
	wavy.emplace_back(e_origin::NOT_YET_DEFINED);
	wavy.emplace_back(e_origin::NOT_YET_DEFINED);
	wavy[0].set_path(scratch.file("first.wfn"));
	wavy[1].set_path(scratch.file("second.wfn"));
	cube saved(std::array<int, 3>{2, 2, 2});
	saved.set_path(scratch.file("saved.cube"));
	write_text(scratch.file("saved.cube"), "x\n");
	cube mem_a(std::array<int, 3>{2, 2, 2});
	mem_a.set_path(scratch.file("mem_a.cube"));
	cube mem_b(std::array<int, 3>{2, 2, 2});
	mem_b.set_path(scratch.file("mem_b.cube"));
	wavy[0].push_back_cube(saved);
	wavy[1].push_back_cube(mem_a);
	wavy[1].push_back_cube(mem_b);

	std::vector<std::vector<unsigned int>> selection(2, std::vector<unsigned int>(2, 99));
	std::string out;
	{
		//no dot, wfn out of range, cube out of range, then two valid picks
		CinFeed feed("x\n5.0\n1.7\n0.0\n1.1\n");
		CoutCapture cap;
		select_cubes(selection, wavy, 2, false, true);
		out = cap.str();
	}
	EXPECT_EQ(selection[0][0], 0u);
	EXPECT_EQ(selection[1][0], 0u);
	EXPECT_EQ(selection[0][1], 1u);
	EXPECT_EQ(selection[1][1], 1u);
	EXPECT_NE(out.find("Need to select 2 files in total."), std::string::npos);
	EXPECT_NE(out.find("(MEM ONLY)"), std::string::npos);
	EXPECT_NE(out.find("no . found in input!"), std::string::npos);
	EXPECT_NE(out.find("Invalid choice!"), std::string::npos);
	EXPECT_NE(out.find("Translated: 1 1"), std::string::npos);
	EXPECT_NE(out.find("Going to return!"), std::string::npos);
}

TEST(ConvenienceCoverageCubeTests, SelectCubesWfnOnlyPushesTheWavefunctionIndex)
{
	std::vector<WFN> wavy;
	wavy.emplace_back(e_origin::NOT_YET_DEFINED);
	wavy.emplace_back(e_origin::NOT_YET_DEFINED);
	std::vector<std::vector<unsigned int>> selection(2);
	std::string out;
	{
		CinFeed feed("7\n1\n");
		CoutCapture cap;
		select_cubes(selection, wavy, 1, true, false);
		out = cap.str();
	}
	ASSERT_EQ(selection[0].size(), 1u);
	EXPECT_EQ(selection[0][0], 1u);
	EXPECT_TRUE(selection[1].empty());
	EXPECT_NE(out.find("Need to select 1 file."), std::string::npos);
	EXPECT_NE(out.find("Ignoring the .!"), std::string::npos);
	EXPECT_NE(out.find("Invalid choice!"), std::string::npos);
}

TEST(ConvenienceCoverageCubeTests, UnsavedFilesReportsCubesWithoutAFile)
{
	ScratchDir scratch("unsaved");
	std::vector<WFN> wavy;
	wavy.emplace_back(e_origin::NOT_YET_DEFINED);
	EXPECT_FALSE(unsaved_files(wavy)); //no cubes at all
	cube saved(std::array<int, 3>{2, 2, 2});
	saved.set_path(scratch.file("saved.cube"));
	write_text(scratch.file("saved.cube"), "x\n");
	wavy[0].push_back_cube(saved);
	EXPECT_FALSE(unsaved_files(wavy));
	cube mem(std::array<int, 3>{2, 2, 2});
	mem.set_path(scratch.file("never_written.cube"));
	wavy[0].push_back_cube(mem);
	EXPECT_TRUE(unsaved_files(wavy));
}

//---------------------------------------------------------------- Fortran record of ints

TEST(ConvenienceCoverageBinaryTests, IntRecordReadsPayloadAndRejectsDamage)
{
	ScratchDir scratch("fortran_int");
	auto write_record = [&](const std::string& name, const int head, const ivec& payload, const int tail)
	{
		std::ofstream f(scratch.file(name), std::ios::binary);
		f.write(reinterpret_cast<const char*>(&head), sizeof(int));
		f.write(reinterpret_cast<const char*>(payload.data()), payload.size() * sizeof(int));
		f.write(reinterpret_cast<const char*>(&tail), sizeof(int));
	};
	write_record("good.bin", 12, { 1, 2, 3 }, 12);
	write_record("tail.bin", 12, { 1, 2, 3 }, 8);
	write_record("neg.bin", -4, {}, -4);

	ivec v;
	{
		std::ifstream f(scratch.file("good.bin"), std::ios::binary);
		EXPECT_TRUE(read_block_from_fortran_binary(f, v));
	}
	ASSERT_EQ(v.size(), 3u);
	EXPECT_EQ(v[0], 1);
	EXPECT_EQ(v[1], 2);
	EXPECT_EQ(v[2], 3);
	{
		std::ifstream f(scratch.file("tail.bin"), std::ios::binary);
		ivec t;
		CoutCapture cap;
		EXPECT_FALSE(read_block_from_fortran_binary(f, t));
		EXPECT_NE(cap.str().find("12 vs. 8"), std::string::npos);
	}
	{
		std::ifstream f(scratch.file("neg.bin"), std::ios::binary);
		ivec n;
		CoutCapture cap;
		EXPECT_FALSE(read_block_from_fortran_binary(f, n));
		EXPECT_NE(cap.str().find("record of -4 bytes"), std::string::npos);
	}
}

//---------------------------------------------------------------- helpers with uncovered branches

TEST(ConvenienceCoverageHelperTests, ShrinkStringRemovesEveryDigit)
{
	std::string s = "N3456789";
	EXPECT_EQ(shrink_string(s), "N");
	std::string t = "Cl(3) 4 5 6 7 8 9 0";
	EXPECT_EQ(shrink_string(t), "Cl");
}

TEST(ConvenienceCoverageHelperTests, CheckBohrDebugPrintsDecision)
{
	WFN close(e_origin::NOT_YET_DEFINED);
	close.push_back_atom("O", 0.0, 0.0, 0.0, 8);
	close.push_back_atom("H", 1.0, 0.0, 0.0, 1);
	close.push_back_atom("H", 0.0, 3.0, 0.0, 1);
	{
		CoutCapture cap;
		EXPECT_FALSE(check_bohr(close, true));
		EXPECT_NE(cap.str().find("Length for: 0;1: 1, min_length: 300"), std::string::npos);
		EXPECT_NE(cap.str().find("Decided it's written in Angstrom"), std::string::npos);
	}
	WFN wide(e_origin::NOT_YET_DEFINED);
	wide.push_back_atom("O", 0.0, 0.0, 0.0, 8);
	wide.push_back_atom("H", 2.5, 0.0, 0.0, 1);
	{
		CoutCapture cap;
		EXPECT_TRUE(check_bohr(wide, true));
		EXPECT_NE(cap.str().find("Decided it's written in Bohr"), std::string::npos);
	}
}

TEST(ConvenienceCoverageHelperTests, Shell2FunctionHighAndNegativeShells)
{
	for (int prim = 0; prim < 4; prim++)
	{
		EXPECT_EQ(shell2function(-5, prim), -32 + prim);
		EXPECT_EQ(shell2function(-4, prim), -21 + prim);
		EXPECT_EQ(shell2function(-3, prim), -12 + prim);
		EXPECT_EQ(shell2function(5, prim), 36 + prim);
	}
	EXPECT_EQ(shell2function(3, 10), 0); //an f shell has ten functions
	EXPECT_EQ(shell2function(3, 11), 0);
}

TEST(ConvenienceCoverageHelperTests, BesselRatioMatchesSphericalClosedForms)
{
	//J_{nu+1}/J_nu at nu = l + 1/2 is j_{l+1}/j_l
	const double x = 1.0;
	const double j0 = std::sin(x) / x;
	const double j1 = (std::sin(x) / x - std::cos(x)) / x;
	const double j2 = ((3.0 / (x * x) - 1.0) * std::sin(x) - 3.0 * std::cos(x) / x) / x;
	EXPECT_NEAR(get_bessel_ratio(0.5, x), j1 / j0, 1e-13);
	EXPECT_NEAR(get_bessel_ratio(1.5, x), j2 / j1, 1e-13);
	const double y = 2.5;
	const double k0 = std::sin(y) / y;
	const double k1 = (std::sin(y) / y - std::cos(y)) / y;
	EXPECT_NEAR(get_bessel_ratio(0.5, y), k1 / k0, 1e-13);
}

TEST(ConvenienceCoverageHelperTests, BesselRatioRescalesTinyContinuedFraction)
{
	//x = DBL_MIN: the first convergent x/3 is subnormal, the rescale by DBL_MIN keeps the ratio x/3;
	//every later partial numerator underflows to zero so the fraction is exact
	const double x = DBL_MIN;
	const double ratio = get_bessel_ratio(0.5, x);
	EXPECT_GT(ratio, 0.0);
	EXPECT_NEAR(ratio, x / 3.0, x / 3.0 * 1e-6);
}

TEST(ConvenienceCoverageHelperTests, Sha256UpdateAccumulatesAcrossCalls)
{
	uint32_t state[8] = { 0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19 };
	uint8_t buffer[64];
	std::memset(buffer, 0, sizeof(buffer));
	uint64_t bitlen = 0;
	sha::sha256_update(state, buffer, reinterpret_cast<const uint8_t*>("a"), 1, bitlen);
	EXPECT_EQ(bitlen, 8u);
	sha::sha256_update(state, buffer, reinterpret_cast<const uint8_t*>("bc"), 2, bitlen);
	EXPECT_EQ(bitlen, 24u);
	uint8_t hash[32];
	sha::sha256_final(state, buffer, bitlen, hash);
	EXPECT_EQ(hex_digest(hash), "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");

	//56 bytes in 5-byte pieces: the block boundary is crossed by the padding in final, not by update
	const std::string msg = "abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq";
	uint32_t state2[8] = { 0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19 };
	uint8_t buffer2[64];
	std::memset(buffer2, 0, sizeof(buffer2));
	uint64_t bitlen2 = 0;
	for (size_t pos = 0; pos < msg.size(); pos += 5)
	{
		const size_t n = std::min<size_t>(5, msg.size() - pos);
		sha::sha256_update(state2, buffer2, reinterpret_cast<const uint8_t*>(msg.data() + pos), n, bitlen2);
	}
	EXPECT_EQ(bitlen2, 448u);
	sha::sha256_final(state2, buffer2, bitlen2, hash);
	EXPECT_EQ(hex_digest(hash), "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1");

	//a message of 64 bytes makes update itself run the transform once
	const std::string block(64, 'a');
	uint32_t state3[8] = { 0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19 };
	uint8_t buffer3[64];
	std::memset(buffer3, 0, sizeof(buffer3));
	uint64_t bitlen3 = 0;
	sha::sha256_update(state3, buffer3, reinterpret_cast<const uint8_t*>(block.data()), block.size(), bitlen3);
	EXPECT_EQ(bitlen3, 512u);
	sha::sha256_final(state3, buffer3, bitlen3, hash);
	EXPECT_EQ(hex_digest(hash), sha::sha256(block));
}

//---------------------------------------------------------------- process and file helpers

TEST(ConvenienceCoverageProcessTests, HomePathComesFromTheEnvironment)
{
	const std::filesystem::path home = get_home_path();
#ifdef _WIN32
	const char* drive = std::getenv("HOMEDRIVE");
	const char* rest = std::getenv("HOMEPATH");
	std::string expected;
	if (drive != nullptr)
		expected += drive;
	if (rest != nullptr)
		expected += rest;
	EXPECT_EQ(home.string(), expected);
#else
	const char* env_home = std::getenv("HOME");
	if (env_home != nullptr && env_home[0] != '\0')
		EXPECT_EQ(home.string(), std::string(env_home));
	else
		EXPECT_EQ(home.string(), "/tmp");
#endif
}

TEST(ConvenienceCoverageProcessTests, CopyFileCopiesBytes)
{
	ScratchDir scratch("copy_file");
	const std::string payload("bytes\x00\x01\xff and text\n", 20);
	{
		std::ofstream f(scratch.file("from.bin"), std::ios::binary);
		f.write(payload.data(), payload.size());
	}
	std::filesystem::path from = scratch.file("from.bin");
	std::filesystem::path to = scratch.file("to.bin");
	copy_file(from, to);
	std::ifstream f(to, std::ios::binary);
	std::string copied((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
	EXPECT_EQ(copied, payload);
	EXPECT_EQ(std::filesystem::file_size(to), payload.size());
}

TEST(ConvenienceCoverageProcessTests, AvailableMemoryIsPositive)
{
	EXPECT_GT(available_memory_bytes(), 0u);
}

TEST(ConvenienceCoverageProcessTests, EnsureOccDataPathAcceptsValidEnvironmentOrSearchesTheExeDirectory)
{
	const char* env = std::getenv("OCC_DATA_PATH");
	if (env != nullptr && valid_occ_dir(env))
	{
		//ctest sets a valid path; the early return needs no executable
		EXPECT_TRUE(ensure_occ_data_path("no_such_program"));
		EXPECT_TRUE(ensure_occ_data_path(nullptr));
		return;
	}
	ScratchDir scratch("occ_data");
	const std::string exe = scratch.file("NoSpherA2.exe").string();
	//nothing valid next to the executable
	EXPECT_FALSE(ensure_occ_data_path(exe.c_str()));
	//exe_dir/occ/share is the second candidate
	std::filesystem::create_directories(scratch.file("occ") / "share" / "basis");
	std::filesystem::create_directories(scratch.file("occ") / "share" / "methods");
	EXPECT_TRUE(ensure_occ_data_path(exe.c_str()));
	const char* set = std::getenv("OCC_DATA_PATH");
	ASSERT_NE(set, nullptr);
	EXPECT_EQ(std::filesystem::path(set), scratch.file("occ") / "share");
	//and now the environment alone satisfies it
	EXPECT_TRUE(ensure_occ_data_path(nullptr));
}

//---------------------------------------------------------------- option digestion

TEST(ConvenienceCoverageOptionTests, FinishedStopsTheParserBeforeTheNextFlag)
{
	const options opt = parse({ "-lahvatest", "-charge", "5" });
	EXPECT_TRUE(opt.finished);
	EXPECT_EQ(opt.charge, 0);
}

TEST(ConvenienceCoverageOptionTests, TscLabelsAloneOnlySetsTheFlag)
{
	const options opt = parse({ "-tsc_labels" });
	EXPECT_TRUE(opt.label_tsc_output);
	EXPECT_FALSE(opt.finished);
	const options next = parse({ "-tsc_labels", "-charge", "3" });
	EXPECT_TRUE(next.label_tsc_output);
	EXPECT_EQ(next.charge, 3);
}

TEST(ConvenienceCoverageOptionTests, AnomDispOccAndSimpleFlags)
{
	ScratchDir scratch("occ_flag");
	write_text(scratch.file("occ.toml"), "x\n");
	const options opt = parse({ "-anom_disp", "disp.dat", "-occ", scratch.file("occ.toml").string(), "-all_charges",
		"-s_rho", "-no_gpu_cublas", "-geometry_aid_metals", "-test" });
	EXPECT_EQ(opt.anom_disp_path, std::filesystem::path("disp.dat"));
	EXPECT_EQ(opt.occ, scratch.file("occ.toml").string());
	EXPECT_TRUE(opt.all_charges);
	EXPECT_TRUE(opt.properties.s_rho);
	EXPECT_FALSE(opt.gpu_cublas);
	EXPECT_TRUE(opt.geometry_aid_metals);
	EXPECT_TRUE(opt.test);
}

TEST(ConvenienceCoverageOptionTests, GflopsEnablesThroughputTracking)
{
	const bool before = throughput::enabled();
	const options opt = parse({ "-gflops" });
	EXPECT_TRUE(opt.track_gflops);
	EXPECT_TRUE(throughput::enabled());
	throughput::set_enabled(before);
}

TEST(ConvenienceCoverageOptionTests, MultipoleMomentsMbisAddsAnAutoAuxBasis)
{
	const options opt = parse({ "-multipole_moments", "MBIS", "2" });
	EXPECT_EQ(opt.multipole_scheme, PartitionType::MBIS);
	EXPECT_EQ(opt.multipole_lmax, 2);
	EXPECT_TRUE(opt.RI_FIT);
	EXPECT_EQ(opt.partition_type, PartitionType::RI);
	EXPECT_EQ(opt.aux_basis.size(), 1u);
	const options embis = parse({ "-multipole_moments", "embis", "0", "-multipole_strength", "0.5" });
	EXPECT_EQ(embis.multipole_scheme, PartitionType::EMBIS);
	EXPECT_EQ(embis.multipole_lmax, 0);
	EXPECT_DOUBLE_EQ(embis.multipole_strength, 0.5);
}

TEST(ConvenienceCoverageOptionTests, TwinLawWithDebugIsEchoed)
{
	std::string out;
	const options opt = parse_quiet({ "-twin", "1", "0", "0", "0", "1", "0", "0", "0", "-1", "-charge", "2" }, true, out);
	ASSERT_EQ(opt.twin_law.size(), 1u);
	ASSERT_EQ(opt.twin_law[0].size(), 9u);
	EXPECT_DOUBLE_EQ(opt.twin_law[0][0], 1.0);
	EXPECT_DOUBLE_EQ(opt.twin_law[0][4], 1.0);
	EXPECT_DOUBLE_EQ(opt.twin_law[0][8], -1.0);
	EXPECT_EQ(opt.charge, 2);
	EXPECT_NE(out.find("twin_law:"), std::string::npos);
}

TEST(ConvenienceCoverageOptionTests, CmtcWithDebugCollectsFilesCifsAndGroups)
{
	std::string out;
	const options opt = parse_quiet({ "-cmtc", "a.wfn", "a.cif", "1,2", "b.wfn", "b.cif", "3", "-charge", "1" }, true, out);
	EXPECT_TRUE(opt.cif_based_combined_tsc_calc);
	ASSERT_EQ(opt.combined_tsc_calc_files.size(), 2u);
	EXPECT_EQ(opt.combined_tsc_calc_files[1], std::filesystem::path("b.wfn"));
	ASSERT_EQ(opt.combined_tsc_calc_cifs.size(), 2u);
	EXPECT_EQ(opt.combined_tsc_calc_cifs[0], std::filesystem::path("a.cif"));
	ASSERT_EQ(opt.groups.size(), 2u);
	EXPECT_EQ(opt.groups[0], (ivec{ 1, 2 }));
	EXPECT_EQ(opt.groups[1], (ivec{ 3 }));
	EXPECT_EQ(opt.charge, 1);
	EXPECT_NE(out.find("--Group: 1,2"), std::string::npos);
	EXPECT_NE(out.find("--Delimiter not found, using ."), std::string::npos);
}

TEST(ConvenienceCoverageOptionTests, MtcWithDebugCollectsFilesAndGroups)
{
	std::string out;
	const options opt = parse_quiet({ "-mtc", "a.wfn", "1,2", "b.wfn", "3" }, true, out);
	EXPECT_TRUE(opt.combined_tsc_calc);
	ASSERT_EQ(opt.combined_tsc_calc_files.size(), 2u);
	EXPECT_EQ(opt.combined_tsc_calc_files[0], std::filesystem::path("a.wfn"));
	ASSERT_EQ(opt.groups.size(), 2u);
	EXPECT_EQ(opt.groups[0], (ivec{ 1, 2 }));
	EXPECT_EQ(opt.groups[1], (ivec{ 3 }));
	EXPECT_NE(out.find("--Group: 3"), std::string::npos);
}

TEST(ConvenienceCoverageOptionTests, NnlsTestFinishes)
{
	std::string out;
	const options opt = parse_quiet({ "-NNLS_TEST" }, false, out);
	EXPECT_TRUE(opt.finished);
	EXPECT_NE(out.find("NNLS solution:"), std::string::npos);
}

TEST(ConvenienceCoverageOptionTests, DrawOrbitsWritesACubeForValidM)
{
	ScratchDir scratch("draw_orbits");
	scratch.enter();
	std::string out;
	const options opt = parse_quiet({ "-draw_orbits", "1,0,0.5,2.0" }, false, out);
	EXPECT_TRUE(opt.finished);
	EXPECT_DOUBLE_EQ(opt.properties.resolution, 0.5);
	EXPECT_DOUBLE_EQ(opt.properties.radius, 2.0);
	EXPECT_TRUE(std::filesystem::exists(scratch.file("Oribital-lam1-m-0.cube")));
	EXPECT_GT(std::filesystem::file_size(scratch.file("Oribital-lam1-m-0.cube")), 0u);
}

TEST(ConvenienceCoverageOptionTests, DrawOrbitsRefusesMOutsideL)
{
	ScratchDir scratch("draw_orbits_bad_m");
	scratch.enter();
	std::string out;
	const options opt = parse_quiet({ "-draw_orbits", "2,3" }, false, out);
	EXPECT_TRUE(opt.finished);
	EXPECT_DOUBLE_EQ(opt.properties.resolution, 0.025);
	EXPECT_DOUBLE_EQ(opt.properties.radius, 3.5);
	EXPECT_NE(out.find("m must be between -l and l"), std::string::npos);
	EXPECT_FALSE(std::filesystem::exists(scratch.file("Oribital-lam2-m-3.cube")));
}

TEST(ConvenienceCoverageOptionTests, SphericalAtomsWritesOneTableForEveryElement)
{
	ScratchDir scratch("spherical_atoms");
	scratch.enter();
	//write_spherical_atoms prints from an omp parallel for; that is only race-free on the synced
	//stdio buffer, so cout is not captured into a stringbuf here
	const options opt = parse({ "-spherical_atoms" });
	EXPECT_TRUE(opt.finished);
	int files = 0;
	for (const auto& entry : std::filesystem::directory_iterator(scratch.dir))
		if (entry.path().filename().string().rfind("spherical_", 0) == 0)
			files++;
	EXPECT_EQ(files, 102);
	ASSERT_TRUE(std::filesystem::exists(scratch.file("spherical_H.txt")));
	ASSERT_TRUE(std::filesystem::exists(scratch.file("spherical_No.txt")));
	std::ifstream h(scratch.file("spherical_H.txt"));
	std::string header;
	std::getline(h, header);
	EXPECT_EQ(header, "1 r Density");
	double r = 0.0, rho = 0.0;
	h >> r >> rho;
	EXPECT_DOUBLE_EQ(r, 1e-7);
	//the table is a contracted hydrogen (rho(0) = zeta^3 / pi with zeta = 1.15, not the free-atom
	//1 / pi), so the invariant checked is the electron count: 4 pi int rho r^2 dr = 1, trapezoid on
	//the geometric grid (2.5 % steps, error ~1e-4)
	double last_r = r, last_rho = rho, electrons = 0.0;
	int rows = 1;
	while (h >> r >> rho)
	{
		EXPECT_NEAR(r, last_r * 1.025, last_r * 1e-9);
		electrons += 0.5 * (last_rho * last_r * last_r + rho * r * r) * (r - last_r);
		last_r = r;
		last_rho = rho;
		rows++;
	}
	EXPECT_NEAR(4.0 * constants::PI * electrons, 1.0, 2e-3);
	EXPECT_LT(last_rho, 1e-15);
	EXPECT_GT(rows, 100);
}

TEST(ConvenienceCoverageOptionTests, ConvertXcwReadsTheJobNameAndFinishes)
{
	ScratchDir scratch("convert_xcw");
	write_text(scratch.file("stdout.txt"), "Tonto header\nName ... job1\nmore\n");
	std::string out;
	const options opt = parse_quiet({ "-convert_XCW", scratch.file("stdout.txt").string(), "0.1" }, false, out);
	EXPECT_TRUE(opt.finished);
	EXPECT_NE(out.find("jobname: job1"), std::string::npos);
	EXPECT_NE(out.find("lambda step 0.100000"), std::string::npos);
	//direct call, no lambda files next to the stdout so nothing is converted
	{
		CoutCapture cap;
		options direct;
		convert_tonto_XCW_lambda_steps(scratch.file("stdout.txt").string(), "0.25", true, direct);
		EXPECT_NE(cap.str().find("lambda step 0.250000 and jobname: job1"), std::string::npos);
	}
}

TEST(ConvenienceCoverageOptionTests, LookForDebugVerboseFlagAndArgumentCopy)
{
	const bool old_report = ProgressBar::report_counts;
	std::vector<std::string> tokens = { "NoSpherA2", "-v", "-charge", "4" };
	std::vector<char*> argv;
	for (auto& t : tokens)
		argv.push_back(t.data());
	int argc = static_cast<int>(argv.size());
	std::ostringstream log;
	{
		CoutCapture cap;
		options opt(argc, argv.data(), log);
		EXPECT_TRUE(opt.debug);
		EXPECT_TRUE(ProgressBar::report_counts);
		ASSERT_EQ(opt.arguments.size(), 4u);
		EXPECT_EQ(opt.arguments[0], "NoSpherA2");
		EXPECT_EQ(opt.arguments[3], "4");
		EXPECT_NE(cap.str().find("Turning on verbose mode!"), std::string::npos);
	}
	ProgressBar::report_counts = old_report;
}

//---------------------------------------------------------------- death tests: options that refuse their input

TEST(ConvenienceCoverageDeathTest, HelpFlagPrintsAndExitsZero)
{
	std::vector<std::string> tokens = { "NoSpherA2", "-h" };
	std::vector<char*> argv;
	for (auto& t : tokens)
		argv.push_back(t.data());
	int argc = static_cast<int>(argv.size());
	std::ostringstream log;
	EXPECT_EXIT({ options opt(argc, argv.data(), log); }, ::testing::ExitedWithCode(0), ".*");
}

TEST(ConvenienceCoverageDeathTest, MissingOrUnparsableValueDies)
{
	EXPECT_EXIT(parse({ "-charge" }), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({ "-charge", "abc" }), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
}

TEST(ConvenienceCoverageDeathTest, PartitionOptionsDie)
{
	EXPECT_EXIT(parse({ "-tscb", "table.txt" }), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({ "-tsc_labels", "table.dat", "a.cif" }), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({ "-occ", "no_such_occ_input.toml" }), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({ "-SALTED_Training" }), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
}

TEST(ConvenienceCoverageDeathTest, PropertyOptionsDie)
{
	EXPECT_EXIT(parse({ "-cube_density", "no_such.cube" }), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({ "-calc_dens_1D", "0", "1", "10", "2.0" }), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({ "-eli_analysis", "a.wfn" }), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({ "-qtaim_eli" }), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({ "-qtaim_eli", "a.wfn", ",,," }), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
}

TEST(ConvenienceCoverageDeathTest, RiOptionsDie)
{
	EXPECT_EXIT(parse({ "-geometry_aid_cutoff" }), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({ "-classify_atoms" }), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({ "-classify_atoms", "no_such_model.bin" }), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({ "-rgbi_basis", "bogus" }), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({ "-multipole_moments", "bogus", "2" }), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({ "-multipole_moments", "hirshfeld", "9" }), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({ "-multipole_strength", "0" }), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({ "-interaction_energy", "a.wfn" }), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({ "-interaction_energies", "no_such.job" }), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({ "-RI_CUBE", "coefs.npy" }), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({ "-test_RI" }), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({ "-RI_WFN_DIFF" }), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
}

TEST(ConvenienceCoverageDeathTest, XcwAndDevOptionsDie)
{
	EXPECT_EXIT(parse({ "-convert_XCW", "stdout.txt" }), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(parse({ "-lukas_test" }), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
}

TEST(ConvenienceCoverageDeathTest, ConvertXcwMissingStdoutDies)
{
	options opt;
	EXPECT_EXIT(convert_tonto_XCW_lambda_steps("no_such_stdout.txt", "0.1", false, opt),
		::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
}

TEST(ConvenienceCoverageDeathTest, BesselFirstKindRefusesNegativeOrderOrArgument)
{
	EXPECT_EXIT(bessel_first_kind(-1, 1.0), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
	EXPECT_EXIT(bessel_first_kind(0, -1.0), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
}

TEST(ConvenienceCoverageDeathTest, AdpReaderDiesWhenTheCifHasMoreAtomsThanTheWavefunction)
{
	ScratchDir scratch("adp_overflow");
	write_text(scratch.file("a.cif"), cif_text);
	EXPECT_EXIT({
		std::ofstream log(scratch.file("log.txt"));
		WFN w(e_origin::NOT_YET_DEFINED);
		w.push_back_atom("A", constants::ang2bohr(1.0), constants::ang2bohr(2.0), constants::ang2bohr(3.0), 6);
		w.push_back_atom("B", constants::ang2bohr(5.0), constants::ang2bohr(5.0), constants::ang2bohr(5.0), 8);
		cell unit_cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
		read_fracs_ADPs_from_CIF(scratch.file("a.cif"), w, unit_cell, log, false);
	}, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
}

TEST(ConvenienceCoverageDeathTest, UisoReaderDiesWhenTheCifHasMoreAtomsThanTheWavefunction)
{
	ScratchDir scratch("uiso_overflow");
	write_text(scratch.file("a.cif"), cif_text);
	EXPECT_EXIT({
		std::ofstream log(scratch.file("log.txt"));
		WFN w(e_origin::NOT_YET_DEFINED);
		w.push_back_atom("A", constants::ang2bohr(1.0), constants::ang2bohr(2.0), constants::ang2bohr(3.0), 6);
		w.push_back_atom("B", constants::ang2bohr(5.0), constants::ang2bohr(5.0), constants::ang2bohr(5.0), 8);
		cell unit_cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
		read_U_iso_from_CIF(scratch.file("a.cif"), w, unit_cell, log, false);
	}, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
}
