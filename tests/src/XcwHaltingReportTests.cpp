//XCW: the Gaussian halting reports (evaluate_gaussian_halting, the progress estimate and
//the summary) and the grown branch of the structure_factors constructor, all driven through
//XCW_wrapper::run() on a 402-reflection subset of the P1 fixture that each test writes into its
//own scratch directory. XcwTests.cpp owns the halting maths and the
//settings parser; nothing here repeats those.
#include "pch.h"
#include <gtest/gtest.h>

#include <array>
#include <iomanip>
#include <set>
#include <string>
#include <vector>

#include "core/constants.h"
#include "core/convenience.h"
#include "core/cell.h"
#include "core/scattering_factors.h"
#include "core/XCW_wrapper.h"

namespace {

	std::string test_name()
	{
		return ::testing::UnitTest::GetInstance()->current_test_info()->name();
	}

	//Everything the run prints to std::cout, restored on scope exit
	struct CoutCapture {
		std::ostringstream buffer;
		std::streambuf* old;
		CoutCapture() : old(std::cout.rdbuf(buffer.rdbuf())) {}
		~CoutCapture() { std::cout.rdbuf(old); }
		std::string str() const { return buffer.str(); }
	};

	//cwd back where it was, whatever the test does in between
	struct CwdGuard {
		std::filesystem::path old = std::filesystem::current_path();
		~CwdGuard() { std::filesystem::current_path(old); }
	};

	std::filesystem::path fixture_dir()
	{
		return nos_test_repo_root() / "tests" / "P1_test";
	}

	bool fixture_present()
	{
		return std::filesystem::exists(fixture_dir() / "P1_test.cif") && std::filesystem::exists(fixture_dir() / "P1_test.hkl");
	}

	std::filesystem::path scratch_dir()
	{
		const auto dir = std::filesystem::temp_directory_path() / ("nosphera2_xcw_halting_report_" + test_name());
		std::filesystem::remove_all(dir);
		std::filesystem::create_directories(dir);
		return dir;
	}

	std::string read_file(const std::filesystem::path& p)
	{
		std::stringstream s;
		s << std::ifstream(p).rdbuf();
		return s.str();
	}

	int count_of(const std::string& text, const std::string& needle)
	{
		int n = 0;
		for (size_t pos = text.find(needle); pos != std::string::npos; pos = text.find(needle, pos + needle.size())) {
			n++;
		}
		return n;
	}

	//Every 8th data line of P1_test.hkl, written as it stands (CR stripped) so read_hkl_full
	//parses it exactly as the fixture. The 0 0 0 terminator and repeated hkl are left out so
	//obs[i] and hkl line up one to one. n_strong mirrors read_hkl_full + evaluate_gaussian_halting
	//by hand: stof on the F2 and sigma(F2) fields, |F| = sqrt(|F2|), sigma(F) = sigma(F2)/(2|F|),
	//excluded when sigma(F) <= 0, strong when |F|/sigma(F) >= cutoff.
	void write_subset_hkl(const std::filesystem::path& out, int& n_written, int& n_strong, const double cutoff)
	{
		n_written = 0;
		n_strong = 0;
		std::ifstream in(fixture_dir() / "P1_test.hkl");
		std::ofstream o(out);
		std::set<std::array<int, 3>> seen;
		const std::regex letters{ R"([abcdefghijklmnopqrstuvwxyz\(\)ABCDEFGHIJKLMNOPQRSTUVW])" };
		std::string line;
		int data_index = 0;
		while (std::getline(in, line)) {
			if (!line.empty() && line.back() == '\r') {
				line.pop_back();
			}
			if (line.size() < 2 || std::regex_search(line, letters)) {
				continue;
			}
			const int idx = data_index++;
			if (idx % 8 != 0) {
				continue;
			}
			std::array<int, 3> hkl{};
			for (int i = 0; i < 3; i++) {
				hkl[i] = std::stoi(line.substr(4 * size_t(i), 4));
			}
			if (hkl == std::array<int, 3>{ 0, 0, 0 } || !seen.insert(hkl).second) {
				continue;
			}
			o << line << "\n";
			n_written++;
			const std::string rest = line.substr(12);
			const size_t dot = rest.find_first_of('.');
			const float F2 = std::stof(rest.substr(0, dot + 3));
			const float sigma2 = std::stof(rest.substr(dot + 3));
			const double abs_F = std::sqrt(static_cast<double>(F2 < 0 ? -F2 : F2));
			const double sigma = sigma2 * 0.5 / abs_F;
			if (sigma <= 0.0) {
				continue;
			}
			if (abs_F / sigma >= cutoff) {
				n_strong++;
			}
		}
	}

	//The P1 molecule's stored two-electron integrals, which XcwScfTests.SharedEriBuild writes under ctest; without
	//them the run builds its own, with the same numbers
	std::string shared_eri(const std::filesystem::path& dir)
	{
		const char* d = std::getenv("NOS_TEST_SHARED_DIR");
		if (d == nullptr || *d == '\0' || !std::filesystem::exists(std::filesystem::path(d) / "xcw_p1_subset_eri.bin")) return "";
		std::filesystem::copy_file(std::filesystem::path(d) / "xcw_p1_subset_eri.bin", dir / "shared_eri.bin", std::filesystem::copy_options::overwrite_existing);
		return " read_eri shared_eri.bin";
	}

	//The options the NoSpherA2 driver would build for -XCW on the fixture cif, the subset hkl
	//and the settings text, CPU only. opt must outlive the solver, which keeps a pointer to it.
	options make_options(const std::filesystem::path& dir, const std::string& settings_text)
	{
		std::ofstream(dir / "settings.txt") << settings_text;
		options opt;
		opt.xcw_settings_path = dir / "settings.txt";
		opt.cif = std::filesystem::absolute(fixture_dir() / "P1_test.cif");
		opt.hkl = dir / "subset.hkl";
		opt.do_XCW = true;
		opt.use_gpu = false;
		opt.no_date = true;
		opt.groups[0].push_back(0);
		return opt;
	}

	//Construct and run in dir; returns what went to std::cout
	std::string run_xcw(const std::filesystem::path& dir, options& opt)
	{
		CwdGuard cwd;
		std::filesystem::current_path(dir);
		CoutCapture out;
		opt.loadXCWsettings();
		{
			XCW_wrapper x(opt);
			x.run();
		}
		return out.str();
	}

	//The cout table row for this lambda (" 0.01000\t..."), split on the runs of tabs
	std::vector<std::string> table_row(const std::string& cout_text, const std::string& lambda)
	{
		const std::string key = "\n " + lambda + "\t";
		const size_t start = cout_text.find(key);
		if (start == std::string::npos) {
			return {};
		}
		const size_t end = cout_text.find('\n', start + 1);
		std::istringstream row(cout_text.substr(start + 1, end - start - 1));
		std::vector<std::string> fields;
		for (std::string f; row >> f;) fields.push_back(f);
		return fields;
	}

	const std::string RUN_BASE = "sloppy normal_conv params 177 basis_set sto-3g max_iter 100 F charge 0 mult 1 rhf start 0 step_size 0.01 ";

}

//Six requested lambda steps with the halting criterion on. A^2 never turns upward on this
//fixture, so the scan extends itself by its cap of six more steps: a per-lambda block for each
//of the twelve, progress updates after the fifth and tenth, the summary table with twelve rows
//and the final recommendation on both streams.
TEST(XcwHaltingReportTests, HaltingReportsAcrossTheExtendedLambdaScan)
{
	if (!fixture_present()) {
		GTEST_SKIP() << "P1 fixture missing under " << fixture_dir();
	}
	const auto dir = scratch_dir();
	int n_written = 0, n_strong = 0;
	write_subset_hkl(dir / "subset.hkl", n_written, n_strong, 3.0);
	ASSERT_EQ(n_written, 402);
	ASSERT_GE(n_strong, 8);

	options opt = make_options(dir, RUN_BASE + "end 0.05 gaussian_halt" + shared_eri(dir));
	opt.cpu_itensor_fp32 = false;
	const std::string out = run_xcw(dir, opt);
	const std::string log = read_file(dir / "SCF.log");

	//halting_minimum_beyond_scan: A^2 falls to the end of the requested range on this
	//fixture, so the scan adds its six permitted steps and then says why it stops anyway
	const std::string extending = "so the minimum is outside the requested range - extending the scan";
	EXPECT_EQ(count_of(out, extending), 6) << out;
	EXPECT_EQ(count_of(log, extending), 6) << log;
	EXPECT_EQ(count_of(out, "the scan has already been extended by its limit of 6 steps"), 1) << out;

	//evaluate_gaussian_halting: every step used the same strong set
	EXPECT_EQ(count_of(log, "Gaussian halting criterion at lambda="), 12) << log;
	const std::string n_used = "n_used=" + std::to_string(n_strong) + "/" + std::to_string(n_written) + " (|F|/sigma >= 3.00000)";
	EXPECT_EQ(count_of(log, n_used), 12) << log;
	EXPECT_EQ(count_of(log, "resolution-binned <z^2> trend: slope="), 12);
	EXPECT_EQ(count_of(log, "|F|-binned <z^2> trend: slope="), 12);

	//report_halting_progress_estimate(false) after steps 5 and 10, on both streams
	for (const std::string& n : { std::string("5"), std::string("10") }) {
		const std::string update = "Gaussian halting criterion: progress update after " + n + " lambda steps";
		EXPECT_EQ(count_of(log, update), 1) << log;
		EXPECT_EQ(count_of(out, update), 1) << out;
	}

	//report_gaussian_halting_summary: header and one row per lambda
	EXPECT_NE(log.find("Gaussian halting criterion summary (tests/P1_test/XCW_plan.md)"), std::string::npos);
	EXPECT_NE(log.find(" Lambda\t\tA^2\treject5%\tpp_slope\tpp_intercept\tskew\tkurt\tres_trend_r\tint_trend_r\tn_used\n"), std::string::npos);
	for (int x = 0; x < 12; x++) {
		const std::string row_key = "\n\t0." + std::string(x < 10 ? "0" : "") + std::to_string(x) + "000\t";
		EXPECT_NE(log.find(row_key), std::string::npos) << "summary row missing for step " << x;
	}
	//the final recommendation, once per progress update and once in the summary
	EXPECT_EQ(count_of(log, "Recommended halting lambda* = 0."), 3) << log;
	EXPECT_EQ(count_of(out, "Recommended halting lambda* = 0."), 3) << out;
	//degree 2 needs 5 points (all three reports), degree 4 needs 7 (only the step-5 report is short)
	EXPECT_EQ(count_of(out, "candidate fit: degree=2 RSS="), 3);
	EXPECT_EQ(count_of(out, "candidate fit: degree=4 -- not enough points yet (need >= degree+3 evaluated lambda steps)"), 1);
	EXPECT_EQ(count_of(log, "candidate fit: degree=2 RSS="), 3);

	//the cout table carries the A^2 column, equal to the summary row's A^2
	EXPECT_NE(out.find("Target quantity\t\tCrit(all)\tR1(all)\t\tA^2 (halt)"), std::string::npos) << out;
	const std::vector<std::string> row0 = table_row(out, "0.00000");
	ASSERT_EQ(row0.size(), 10u) << out;
	const size_t srow = log.find("\n\t0.00000\t");
	ASSERT_NE(srow, std::string::npos);
	const size_t a2_start = srow + std::string("\n\t0.00000\t").size();
	const std::string summary_a2 = log.substr(a2_start, log.find('\t', a2_start) - a2_start);
	EXPECT_EQ(row0[9], summary_a2);

	if (!::testing::Test::HasFailure()) {
		std::filesystem::remove_all(dir);
	}
}

//The halting criterion with an impossible cutoff takes the "too few strong reflections" exit,
//which the summary still lists and the recommendation stays silent about
TEST(XcwHaltingReportTests, TooFewStrongReflectionsSkipTheHalting)
{
	if (!fixture_present()) {
		GTEST_SKIP() << "P1 fixture missing under " << fixture_dir();
	}
	const auto dir = scratch_dir();
	int n_written = 0, n_strong = 0;
	write_subset_hkl(dir / "subset.hkl", n_written, n_strong, 1e9);
	ASSERT_EQ(n_written, 402);
	ASSERT_EQ(n_strong, 0);

	options opt = make_options(dir, RUN_BASE + "end 0 gaussian_halt strong_cutoff 1e9" + shared_eri(dir));
	opt.no_date = false;
	const std::string out = run_xcw(dir, opt);
	const std::string log = read_file(dir / "SCF.log");

	//the skip branch of evaluate_gaussian_halting, then a summary row with n_used 0 and no lambda*
	EXPECT_NE(log.find("Gaussian halting criterion: only 0 strong reflections (|F|/sigma >= "), std::string::npos) << log;
	EXPECT_NE(log.find(", skipping (need >= 8)."), std::string::npos);
	EXPECT_EQ(log.find("Gaussian halting criterion at lambda="), std::string::npos);
	EXPECT_NE(log.find("Gaussian halting criterion summary (tests/P1_test/XCW_plan.md)"), std::string::npos);
	//a skipped entry keeps the GaussianHaltEntry defaults: every statistic 0, n_used 0
	EXPECT_NE(log.find("\n\t0.00000\t0.0000\tno\t\t0.0000\t\t0.0000\t0.0000\t0.0000\t0.0000\t\t0.0000\t0\n"), std::string::npos) << log;
	EXPECT_EQ(log.find("Recommended halting lambda*"), std::string::npos);
	EXPECT_EQ(out.find("Recommended halting lambda*"), std::string::npos);
	const std::vector<std::string> row = table_row(out, "0.00000");
	ASSERT_EQ(row.size(), 10u) << out;
	EXPECT_EQ(row[9], "0.0000");

	//no_date off prints the grid setup's timing
	EXPECT_NE(out.find("Time taken for AO values computation: "), std::string::npos) << out;
	EXPECT_NE(out.find("Finished XCW fitting procedure."), std::string::npos);

	if (!::testing::Test::HasFailure()) {
		std::filesystem::remove_all(dir);
	}
}

//`grown` with an xyz that holds exactly the asymmetric unit: the xyz is read, nothing is
//added, symmetry linking and the grown U_iso / ADP paths run on the 23 atoms and the
//constructors finish as they do without the keyword.
TEST(XcwHaltingReportTests, GrownConstructUsesTheXyzAtoms)
{
	if (!fixture_present()) {
		GTEST_SKIP() << "P1 fixture missing under " << fixture_dir();
	}
	const auto dir = scratch_dir();
	const auto cif = std::filesystem::absolute(fixture_dir() / "P1_test.cif");

	//the asymmetric unit in Angstrom, from the same reader structure_factors uses
	std::vector<asym_atom> atoms;
	{
		cell unit_cell(cif, std::cout, false, true);
		std::ifstream cif_in(cif);
		int ncen = 0;
		vec3 ADPs;
		read_CIF(cif_in, unit_cell, ncen, atoms, ADPs, false);
		ASSERT_EQ(ncen, 23);
		ASSERT_EQ(atoms.size(), 23u);
	}
	{
		std::ofstream xyz(dir / "asym.xyz");
		xyz << atoms.size() << "\nP1 asymmetric unit\n" << std::fixed << std::setprecision(6);
		for (const asym_atom& a : atoms) {
			xyz << constants::atnr2letter(a.type) << " " << constants::bohr2ang(a.pos[0]) << " " << constants::bohr2ang(a.pos[1]) << " " << constants::bohr2ang(a.pos[2]) << "\n";
		}
	}
	int n_written = 0, n_strong = 0;
	write_subset_hkl(dir / "subset.hkl", n_written, n_strong, 3.0);
	ASSERT_EQ(n_written, 402);

	options opt = make_options(dir, "normal params 177 basis_set sto-3g charge 0 mult 1 rhf start 0 step_size 0.01 end 0 grown");
	opt.xyz_file = dir / "asym.xyz";
	std::string out;
	{
		CwdGuard cwd;
		std::filesystem::current_path(dir);
		CoutCapture capture;
		opt.loadXCWsettings();
		{
			XCW_wrapper x(opt);
			EXPECT_TRUE(opt.xcw_settings.grown);
			EXPECT_EQ(x.structure().model_data.ncen, 23);
		}
		out = capture.str();
	}
	EXPECT_EQ(out.find("Grown structures require an xyz file"), std::string::npos) << out;
	EXPECT_NE(out.find("Nr of reflections read from file: 402"), std::string::npos) << out;
	EXPECT_NE(out.find("XCW orbital basis set: sto-3g"), std::string::npos) << out;

	if (!::testing::Test::HasFailure()) {
		std::filesystem::remove_all(dir);
	}
}

//`grown` without -xyz: the structure_factors constructor's err_checkf ends the process with
//the error_check exit code (the message goes to std::cout, so no stderr pattern)
TEST(XcwHaltingReportTests, GrownWithoutXyzFileDies)
{
	if (!fixture_present()) {
		GTEST_SKIP() << "P1 fixture missing under " << fixture_dir();
	}
	const auto dir = scratch_dir();
	int n_written = 0, n_strong = 0;
	write_subset_hkl(dir / "subset.hkl", n_written, n_strong, 3.0);
	options opt = make_options(dir, "normal params 177 basis_set sto-3g charge 0 mult 1 rhf start 0 step_size 0.01 end 0 grown");
	ASSERT_TRUE(opt.xyz_file.empty());
	opt.loadXCWsettings();
	EXPECT_EXIT({
		structure_factors SF(opt);
		}, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), "");

	if (!::testing::Test::HasFailure()) {
		std::filesystem::remove_all(dir);
	}
}
