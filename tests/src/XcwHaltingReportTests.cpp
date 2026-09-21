//XCW: the Gaussian halting reports (evaluate_gaussian_halting, the progress estimate and
//the summary), the disk-backed I tensor (decide_i_storage / i_budget / start_i_save /
//finish_i_save / open_i_stream_for_reading) and the grown branch of construct(), all
//driven through run_XCW_fitting() on a 402-reflection subset of the P1 fixture that each
//test writes into its own scratch directory. XcwTests.cpp owns the halting maths and the
//settings parser; nothing here repeats those.
#include "pch.h"
#include <gtest/gtest.h>

#include <array>
#include <cstdint>
#include <iomanip>
#include <set>
#include <string>
#include <vector>

#include "core/constants.h"
#include "core/convenience.h"
#include "core/cell.h"
#include "core/scattering_factors.h"
#include "core/i_tensor_stream.h"
#include "core/XCW.h"

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
		return std::filesystem::exists(fixture_dir() / "P1_test_NA2.cif") && std::filesystem::exists(fixture_dir() / "P1_test.hkl");
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

	//A short anomalous dispersion table with blank lines, which parse_anom_atoms skips
	void write_anom(const std::filesystem::path& p)
	{
		std::ofstream(p) << "\nC 0.00313 0.00162\n\nCl 0.14908 0.15974\nH 0.0 0.0\nN 0.00611 0.00317\nO 0.01110 0.00600\nS 0.12463 0.12335\n\n";
	}

	//The options the NoSpherA2 driver would build for -XCW on the fixture cif, the subset hkl
	//and the settings text, CPU only. opt must outlive the XCW, which keeps a pointer to it.
	options make_options(const std::filesystem::path& dir, const std::string& settings_text)
	{
		std::ofstream(dir / "settings.txt") << settings_text;
		options opt;
		opt.xcw_settings_path = dir / "settings.txt";
		opt.cif = std::filesystem::absolute(fixture_dir() / "P1_test_NA2.cif");
		opt.hkl = dir / "subset.hkl";
		opt.anom_disp_path = dir / "anom.txt";
		opt.do_XCW = true;
		opt.use_gpu = false;
		opt.no_date = true;
		opt.groups[0].push_back(0);
		return opt;
	}

	//Construct and run in dir; returns what went to std::cout
	std::string run_xcw(const std::filesystem::path& dir, const options& opt)
	{
		CwdGuard cwd;
		std::filesystem::current_path(dir);
		CoutCapture out;
		{
			XCW x(opt);
			x.run_XCW_fitting();
		}
		return out.str();
	}

	//The cout table row for this lambda ("0.01000\t\t..."), split on the double tabs
	std::vector<std::string> table_row(const std::string& cout_text, const std::string& lambda)
	{
		const std::string key = "\n" + lambda + "\t\t";
		const size_t start = cout_text.find(key);
		if (start == std::string::npos) {
			return {};
		}
		const size_t end = cout_text.find('\n', start + 1);
		const std::string row = cout_text.substr(start + 1, end - start - 1);
		std::vector<std::string> fields;
		size_t pos = 0;
		while (true) {
			const size_t next = row.find("\t\t", pos);
			fields.push_back(row.substr(pos, next == std::string::npos ? std::string::npos : next - pos));
			if (next == std::string::npos) {
				break;
			}
			pos = next + 2;
		}
		return fields;
	}

	//"Screened out S unique pairs" -> S; with no_date off eval_I prints its own
	//"Screened out D of C AO-block entries" line first, so the pairs line is located by its tail
	long long screened_pairs(const std::string& cout_text)
	{
		const std::string key = "Screened out ";
		const size_t tail = cout_text.find(" unique pairs of mu, nu");
		if (tail == std::string::npos) {
			return -1;
		}
		const size_t pos = cout_text.rfind(key, tail);
		if (pos == std::string::npos) {
			return -1;
		}
		return std::stoll(cout_text.substr(pos + key.size()));
	}

	//nmo for sto-3g on the P1 molecule: 2 S x 9 + Cl 9 + 7 C x 5 + 4 O x 5 + 3 N x 5 + 6 H x 1
	constexpr int P1_NMO = 103;
	constexpr long long P1_PAIRS = static_cast<long long>(P1_NMO) * (P1_NMO + 1) / 2;
	//i_tensor_file header: 5 int64 fields, then pair_mu and pair_nu as int32 each
	constexpr std::uintmax_t HEADER_FIXED_BYTES = 40;

	const std::string RUN_BASE = "sloppy normal_conv params 177 basis_set sto-3g max_iter 100 F charge 0 mult 1 rhf start 0 step_size 0.01 ";

}

//Six lambda steps with the halting criterion on: a per-lambda block for each, one progress
//update after the fifth step, the summary table with six rows and the final recommendation on
//both streams; the double-precision tensor is held in memory under the i_tensor_mb budget and
//saved on the writer thread, and a second run reads it back instead of rebuilding it.
TEST(XcwHaltingReportTests, HaltingReportsAcrossSixLambdaSteps)
{
	if (!fixture_present()) {
		GTEST_SKIP() << "P1 fixture missing under " << fixture_dir();
	}
	const auto dir = scratch_dir();
	int n_written = 0, n_strong = 0;
	write_subset_hkl(dir / "subset.hkl", n_written, n_strong, 3.0);
	ASSERT_EQ(n_written, 402);
	ASSERT_GE(n_strong, 8);
	write_anom(dir / "anom.txt");

	options opt = make_options(dir, RUN_BASE + "end 0.05 i_tensor_mb 4096 save itensor.bin");
	opt.xcw_gaussian_halt = true;
	opt.cpu_itensor_fp32 = false;
	const std::string out = run_xcw(dir, opt);
	const std::string log = read_file(dir / "XCW.log");

	//evaluate_gaussian_halting: every step used the same strong set
	EXPECT_EQ(count_of(log, "Gaussian halting criterion at lambda="), 6) << log;
	const std::string n_used = "n_used=" + std::to_string(n_strong) + "/" + std::to_string(n_written) + " (|F|/sigma >= 3.00000)";
	EXPECT_EQ(count_of(log, n_used), 6) << log;
	EXPECT_EQ(count_of(log, "resolution-binned <z^2> trend: slope="), 6);
	EXPECT_EQ(count_of(log, "|F|-binned <z^2> trend: slope="), 6);

	//report_halting_progress_estimate(false) once, after step 5 of 6, on both streams
	EXPECT_EQ(count_of(log, "Gaussian halting criterion: progress update after 5 lambda steps"), 1) << log;
	EXPECT_EQ(count_of(out, "Gaussian halting criterion: progress update after 5 lambda steps"), 1) << out;

	//report_gaussian_halting_summary: header and one row per lambda
	EXPECT_NE(log.find("Gaussian halting criterion summary (tests/P1_test/XCW_plan.md)"), std::string::npos);
	EXPECT_NE(log.find(" Lambda\t\tA^2\treject5%\tpp_slope\tpp_intercept\tskew\tkurt\tres_trend_r\tint_trend_r\tn_used\n"), std::string::npos);
	for (int x = 0; x < 6; x++) {
		const std::string row_key = "\n\t0.0" + std::to_string(x) + "000\t";
		EXPECT_NE(log.find(row_key), std::string::npos) << "summary row missing for step " << x;
	}
	//the final recommendation, once in the progress update and once in the summary
	EXPECT_EQ(count_of(log, "Recommended halting lambda* = 0.0"), 2) << log;
	EXPECT_EQ(count_of(out, "Recommended halting lambda* = 0.0"), 2) << out;
	//degree 2 needs 5 points (fits at step 5 and at the end), degree 4 needs 7 (never here)
	EXPECT_EQ(count_of(out, "candidate fit: degree=2 RSS="), 2);
	EXPECT_EQ(count_of(out, "candidate fit: degree=4 -- not enough points yet (need >= degree+3 evaluated lambda steps)"), 2);
	EXPECT_EQ(count_of(log, "candidate fit: degree=2 RSS="), 2);

	//the cout table carries the A^2 column, equal to the summary row's A^2
	EXPECT_NE(out.find("Target quantity \tA^2 (halt)"), std::string::npos) << out;
	const std::vector<std::string> row0 = table_row(out, "0.00000");
	ASSERT_EQ(row0.size(), 7u) << out;
	const size_t srow = log.find("\n\t0.00000\t");
	ASSERT_NE(srow, std::string::npos);
	const size_t a2_start = srow + std::string("\n\t0.00000\t").size();
	const std::string summary_a2 = log.substr(a2_start, log.find('\t', a2_start) - a2_start);
	EXPECT_EQ(row0[6], summary_a2);

	//decide_i_storage: held, double precision, budget named
	EXPECT_NE(out.find("I tensor held in memory: "), std::string::npos) << out;
	EXPECT_NE(out.find(" MB (fits the 4096.00 MB i_tensor_mb budget)"), std::string::npos) << out;
	EXPECT_EQ(out.find("(single precision)"), std::string::npos);
	EXPECT_EQ(out.find("I tensor streamed to disk"), std::string::npos);

	//start_i_save / finish_i_save, the write finished before the run ends
	EXPECT_NE(out.find("Writing the I tensor to itensor.bin in the background; a later run can `read itensor.bin` instead of building it"), std::string::npos) << out;
	const size_t written = out.find("I tensor written to itensor.bin");
	const size_t finished = out.find("Finished XCW fitting procedure.");
	ASSERT_NE(written, std::string::npos) << out;
	ASSERT_NE(finished, std::string::npos) << out;
	EXPECT_LT(written, finished);
	EXPECT_EQ(out.find("Could not write the I tensor"), std::string::npos);

	//the file on disk has the header and nr blocks of kept complex doubles
	const long long screened = screened_pairs(out);
	ASSERT_GT(screened, 0) << out;
	const size_t expected_kept = static_cast<size_t>(P1_PAIRS - screened);
	size_t kept = 0;
	bool single = true;
	ASSERT_TRUE(std::filesystem::exists(dir / "itensor.bin"));
	ASSERT_TRUE(i_tensor_file::matches(dir / "itensor.bin", n_written, P1_NMO, kept, single));
	EXPECT_FALSE(single);
	EXPECT_EQ(kept, expected_kept);
	EXPECT_EQ(std::filesystem::file_size(dir / "itensor.bin"),
		HEADER_FIXED_BYTES + 8 * static_cast<std::uintmax_t>(kept) + static_cast<std::uintmax_t>(n_written) * kept * 16);

	//eval_I_anom_disp read branch: the saved tensor is loaded, not rebuilt, and gives the same rows
	options opt2 = make_options(dir, RUN_BASE + "end 0.01 i_tensor_mb 4096 read itensor.bin");
	opt2.xcw_gaussian_halt = true;
	opt2.cpu_itensor_fp32 = false;
	const std::string out2 = run_xcw(dir, opt2);
	EXPECT_NE(out2.find("I tensor read from itensor.bin ("), std::string::npos) << out2;
	EXPECT_NE(out2.find(" MB), not recomputed, held in memory"), std::string::npos) << out2;
	EXPECT_EQ(out2.find("Screened out"), std::string::npos);
	EXPECT_EQ(out2.find("Writing the I tensor"), std::string::npos);
	EXPECT_EQ(out2.find("NOTE: the tensor on disk"), std::string::npos);
	for (const std::string lambda : { "0.00000", "0.01000" }) {
		const std::vector<std::string> a = table_row(out, lambda);
		const std::vector<std::string> b = table_row(out2, lambda);
		ASSERT_EQ(a.size(), 7u) << out;
		ASSERT_EQ(b.size(), 7u) << out2;
		EXPECT_NEAR(std::stod(a[1]), std::stod(b[1]), 2e-3) << "criterion at " << lambda;
		EXPECT_NEAR(std::stod(a[2]), std::stod(b[2]), 2e-3) << "GooF at " << lambda;
		EXPECT_NEAR(std::stod(a[3]), std::stod(b[3]), 1e-4) << "energy at " << lambda;
		EXPECT_NEAR(std::stod(a[4]), std::stod(b[4]), 2e-3) << "perturbation at " << lambda;
	}

	if (!::testing::Test::HasFailure()) {
		std::filesystem::remove_all(dir);
	}
}

//A -mem budget below one reflection's block forces the single-precision tensor onto disk one
//reflection at a time; the file appears at the default path with the right shape, a second run
//with a bare `read` streams it back through open_i_stream_for_reading and reproduces the row.
//The halting criterion with an impossible cutoff takes the "too few strong reflections" exit,
//which the summary still lists and the recommendation stays silent about.
TEST(XcwHaltingReportTests, StreamedTensorRoundTrip)
{
	if (!fixture_present()) {
		GTEST_SKIP() << "P1 fixture missing under " << fixture_dir();
	}
	const auto dir = scratch_dir();
	int n_written = 0, n_strong = 0;
	write_subset_hkl(dir / "subset.hkl", n_written, n_strong, 1e9);
	ASSERT_EQ(n_written, 402);
	ASSERT_EQ(n_strong, 0);
	write_anom(dir / "anom.txt");

	options opt = make_options(dir, RUN_BASE + "end 0");
	opt.mem_given = true;
	opt.mem = 0.001;
	opt.no_date = false;
	opt.xcw_gaussian_halt = true;
	opt.xcw_strong_cutoff = 1e9;
	const std::string out = run_xcw(dir, opt);
	const std::string log = read_file(dir / "XCW.log");

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
	ASSERT_EQ(row.size(), 7u) << out;
	EXPECT_EQ(row[6], "0.0000");

	//no_date off prints the integral timing
	EXPECT_NE(out.find("Time taken for XCW integrals: "), std::string::npos) << out;

	//decide_i_storage streamed branch with the -mem source and the one-at-a-time note
	const std::regex streamed{ R"(I tensor streamed to disk: [0-9]+\.[0-9]{2} MB \(single precision\) total, 1 of 402 reflections resident \([0-9]+\.[0-9]{2} MB\) to fit -mem \(0\.00 MB\))" };
	EXPECT_TRUE(std::regex_search(out, streamed)) << out;
	const std::regex note{ R"(  NOTE: one reflection alone is [0-9]+\.[0-9]{2} MB, over the budget\. Running one at a time\.)" };
	EXPECT_TRUE(std::regex_search(out, note)) << out;
	EXPECT_EQ(out.find("I tensor held in memory"), std::string::npos);
	EXPECT_EQ(out.find("Writing the I tensor"), std::string::npos);
	EXPECT_NE(out.find("Finished XCW fitting procedure."), std::string::npos);

	//the default file, single precision, one block per reflection
	const long long screened = screened_pairs(out);
	ASSERT_GT(screened, 0) << out;
	const size_t expected_kept = static_cast<size_t>(P1_PAIRS - screened);
	const auto stream_file = dir / "I_tensor_stream.bin";
	ASSERT_TRUE(std::filesystem::exists(stream_file));
	size_t kept = 0;
	bool single = false;
	ASSERT_TRUE(i_tensor_file::matches(stream_file, n_written, P1_NMO, kept, single));
	EXPECT_TRUE(single);
	EXPECT_EQ(kept, expected_kept);
	EXPECT_EQ(std::filesystem::file_size(stream_file),
		HEADER_FIXED_BYTES + 8 * static_cast<std::uintmax_t>(kept) + static_cast<std::uintmax_t>(n_written) * kept * 8);

	//read it back a window at a time under the same budget
	options opt2 = make_options(dir, RUN_BASE + "end 0 read");
	opt2.mem_given = true;
	opt2.mem = 0.001;
	const std::string out2 = run_xcw(dir, opt2);
	EXPECT_NE(out2.find("I tensor read from I_tensor_stream.bin ("), std::string::npos) << out2;
	EXPECT_NE(out2.find(" MB, single precision), not recomputed, read a window at a time"), std::string::npos) << out2;
	EXPECT_EQ(out2.find("Screened out"), std::string::npos);
	EXPECT_EQ(out2.find("I tensor streamed to disk"), std::string::npos);
	const std::vector<std::string> row2 = table_row(out2, "0.00000");
	ASSERT_EQ(row2.size(), 6u) << out2;
	EXPECT_NEAR(std::stod(row[1]), std::stod(row2[1]), 2e-3);
	EXPECT_NEAR(std::stod(row[2]), std::stod(row2[2]), 2e-3);
	EXPECT_NEAR(std::stod(row[3]), std::stod(row2[3]), 1e-4);

	if (!::testing::Test::HasFailure()) {
		std::filesystem::remove_all(dir);
	}
}

//`grown` with an xyz that holds exactly the asymmetric unit: the xyz is read, nothing is
//added, symmetry linking and the grown U_iso / ADP paths run on the 23 atoms and construct
//finishes as it does without the keyword.
TEST(XcwHaltingReportTests, GrownConstructUsesTheXyzAtoms)
{
	if (!fixture_present()) {
		GTEST_SKIP() << "P1 fixture missing under " << fixture_dir();
	}
	const auto dir = scratch_dir();
	const auto cif = std::filesystem::absolute(fixture_dir() / "P1_test_NA2.cif");

	//the asymmetric unit in Angstrom, from the same reader construct uses
	std::vector<asym_atom> atoms;
	{
		cell unit_cell(cif, std::cout, false, true);
		std::ifstream cif_in(cif);
		int ncen = 0;
		bvec needs_grid;
		read_atoms_from_CIF(cif_in, unit_cell, ncen, needs_grid, atoms, false);
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
		{
			XCW x(opt);
		}
		out = capture.str();
	}
	EXPECT_EQ(out.find("I need an xyz file"), std::string::npos) << out;
	EXPECT_NE(out.find("Nr of reflections read from file: 402"), std::string::npos) << out;
	EXPECT_NE(out.find("XCW orbital basis set: sto-3g"), std::string::npos) << out;
	EXPECT_TRUE(std::filesystem::exists(dir / "log3.txt"));
	EXPECT_NE(read_file(dir / "XCW.log").find("XCW orbital basis set: sto-3g"), std::string::npos);

	if (!::testing::Test::HasFailure()) {
		std::filesystem::remove_all(dir);
	}
}

//`grown` without -xyz: construct says so on stderr, then read_xyz's err_checkf on the empty
//path ends the process with the error_check exit code
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
	EXPECT_EXIT({
		XCW x(opt);
		}, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), "I need an xyz file to grow the crystal, but none was provided");

	if (!::testing::Test::HasFailure()) {
		std::filesystem::remove_all(dir);
	}
}
