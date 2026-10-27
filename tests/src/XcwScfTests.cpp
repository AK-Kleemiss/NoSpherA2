//XCW SCF driver: XCW_wrapper::run on the P1 fixture through the settings file, observed through
//SCF.log, stdout and the files a run leaves behind.
//
//Every run uses the sto-3g basis of tests/P1_test at a single lambda on every 8th reflection of
//P1_test.hkl (402 of 3215). The GooF values of the subset are nobody's golden, so the runs check
//determinism (two runs, one answer), physical sanity (finite, bounded, descending) and the one
//number the subset cannot move: at lambda = 0 the perturbation is zero and the converged energy is
//the Hartree-Fock energy of P1_test_XCW.good, -1961.925585145 Eh. One test sets up the grid for H2
//in a cubic cell and compares its structure factors with the analytic Fourier transform of the
//Gaussians.
#include "pch.h"
#include <gtest/gtest.h>

#include "core/convenience.h"
#include "core/XCW_wrapper.h"
#include "core/basis_set.h"
#include "core/tsc_block.h"
#include "core/tuning.h"

#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

	std::string test_name()
	{
		return testing::UnitTest::GetInstance()->current_test_info()->name();
	}

	//Empty when the fixture is not there; the test then skips
	std::filesystem::path p1_fixture()
	{
		const auto fixture = nos_test_repo_root() / "tests" / "P1_test";
		if (!std::filesystem::exists(fixture / "P1_test.cif") || !std::filesystem::exists(fixture / "P1_test.hkl")) {
			return {};
		}
		return fixture;
	}

	std::filesystem::path scratch_dir()
	{
		const auto dir = std::filesystem::temp_directory_path() / ("nosphera2_xcw_scf_" + test_name());
		std::filesystem::remove_all(dir);
		std::filesystem::create_directories(dir);
		return dir;
	}

	//Restores std::cout even when an assertion or an exception leaves the run early
	struct cout_capture {
		std::stringstream buf;
		std::streambuf* old;
		cout_capture() : old(std::cout.rdbuf(buf.rdbuf())) {}
		~cout_capture() { std::cout.rdbuf(old); }
	};

	struct cwd_guard {
		std::filesystem::path old;
		explicit cwd_guard(const std::filesystem::path& dir) : old(std::filesystem::current_path()) { std::filesystem::current_path(dir); }
		~cwd_guard() { std::filesystem::current_path(old); }
	};

	std::string read_text(const std::filesystem::path& p)
	{
		std::stringstream s;
		s << std::ifstream(p).rdbuf();
		return s.str();
	}

	size_t occurrences(const std::string& text, const std::string& needle)
	{
		size_t n = 0;
		for (size_t pos = text.find(needle); pos != std::string::npos; pos = text.find(needle, pos + needle.size())) {
			n++;
		}
		return n;
	}

	//The number that follows key in text; -1 when key is absent or nothing numeric follows,
	//which no expected value below equals
	double number_after(const std::string& text, const std::string& key)
	{
		const size_t pos = text.find(key);
		if (pos == std::string::npos) {
			return -1.0;
		}
		std::istringstream in(text.substr(pos + key.size()));
		double v = -1.0;
		in >> v;
		return in.fail() ? -1.0 : v;
	}

	struct p1_run {
		std::string out;    //what the constructors and XCW_wrapper::run printed
		std::string log;    //SCF.log
		int nr = 0;         //reflections in the hkl file the run was given
	};

	//Every 8th data line of P1_test.hkl as it stands (CR stripped), the 0 0 0 terminator and
	//repeated hkl left out, so the count is what read_hkl_full keeps. 402 lines.
	int write_subset_hkl(const std::filesystem::path& out)
	{
		std::ifstream in(p1_fixture() / "P1_test.hkl");
		std::ofstream o(out);
		std::set<std::string> seen;
		std::string line;
		int n = 0, data_index = 0;
		while (std::getline(in, line)) {
			if (!line.empty() && line.back() == '\r') line.pop_back();
			if (line.size() < 12 || data_index++ % 8 != 0) continue;
			const std::string hkl = line.substr(0, 12);
			if (hkl == "   0   0   0" || !seen.insert(hkl).second) continue;
			o << line << "\n";
			n++;
		}
		return n;
	}

	//The options of an XCW run in dir on cif and hkl with the given settings file, no GPU and no timing
	//lines; the anomalous dispersion comes from the CIF
	options xcw_options(const std::filesystem::path& dir, const std::filesystem::path& cif, const std::filesystem::path& hkl,
		const std::string& settings_text, const bool double_grid = false, const int accuracy = 2)
	{
		const auto settings = dir / "settings.txt";
		std::ofstream(settings) << settings_text;
		options opt;
		opt.xcw_settings_path = settings;
		opt.cif = std::filesystem::absolute(cif);
		opt.hkl = std::filesystem::absolute(hkl);
		opt.do_XCW = true;
		opt.use_gpu = false;
		opt.cpu_itensor_fp32 = !double_grid;
		opt.accuracy = accuracy;
		opt.no_date = true;
		//As the driver does before constructing the solver: the tscb of a converged step is built
		//from CIF disorder group 0 (NoSpherA2.cpp, do_XCW branch)
		opt.groups[0].push_back(0);
		return opt;
	}

	//One XCW run with xcw_options; print_F writes the F_calc table of each converged step
	p1_run run_xcw(const std::filesystem::path& dir, const std::filesystem::path& cif, const std::filesystem::path& hkl,
		const std::string& settings_text, const bool print_F = false)
	{
		options opt = xcw_options(dir, cif, hkl, settings_text);
		opt.print_F = print_F;
		p1_run r;
		{
			cwd_guard cwd(dir);
			cout_capture capture;
			opt.loadXCWsettings();
			{
				XCW_wrapper x(opt);
				x.run();
			}
			r.out = capture.buf.str();
		}
		r.log = read_text(dir / "SCF.log");
		return r;
	}

	//The P1 fixture on the reflection subset, written into dir on the first call
	p1_run run_on_p1(const std::filesystem::path& dir, const std::string& settings_text, const bool print_F = false)
	{
		const auto fixture = p1_fixture();
		const auto hkl = dir / "subset.hkl";
		const int nr = write_subset_hkl(hkl);
		p1_run r = run_xcw(dir, fixture / "P1_test.cif", hkl, settings_text, print_F);
		r.nr = nr;
		return r;
	}

	//The stored two-electron integrals depend on geometry and basis only. Under ctest, SharedEriBuild writes the P1
	//molecule's into NOS_TEST_SHARED_DIR once and the runs that test the SCF rather than these read them; without that
	//directory or file every run builds its own, with the same numbers.
	std::filesystem::path shared_eri()
	{
		const char* d = std::getenv("NOS_TEST_SHARED_DIR");
		if (d == nullptr || *d == '\0') return {};
		return std::filesystem::path(d) / "xcw_p1_subset_eri.bin";
	}

	//A read path is one settings token, so the file is copied into dir rather than named where it lies
	std::string read_shared_eri(const std::filesystem::path& dir)
	{
		const auto p = shared_eri();
		if (p.empty() || !std::filesystem::exists(p)) return "";
		std::filesystem::copy_file(p, dir / "shared_eri.bin", std::filesystem::copy_options::overwrite_existing);
		return " read_eri shared_eri.bin";
	}

	//One line of the lambda table XCW_wrapper::run prints per converged step:
	//lambda(5) criterion(4) GooF2(4) R1(5) energy(9) lambda*criterion(3) quant(9), tab separated.
	//An SCF.log iteration row has the same seven columns behind a leading tab, with the
	//iteration count in the first.
	struct lambda_row {
		std::string lambda, criterion, goof2, r1, energy, penalty, quant;
		double d(const std::string& s) const { return std::stod(s); }
	};

	std::vector<lambda_row> table_rows(const std::string& text, const bool iterations)
	{
		std::vector<lambda_row> rows;
		std::istringstream in(text);
		std::string line;
		while (std::getline(in, line)) {
			//both start with a space or a tab; a lambda carries a decimal point, an iteration count not
			if (line.empty() || (line[0] != ' ' && line[0] != '\t') || line.find('\t') == std::string::npos) continue;
			std::vector<std::string> fields;
			std::istringstream cols(line);
			std::string f;
			while (cols >> f) fields.push_back(f);
			if (fields.size() < 7 || !std::isdigit(static_cast<unsigned char>(fields[0][0])) || (fields[0].find('.') == std::string::npos) == !iterations) {
				continue;
			}
			rows.push_back({ fields[0], fields[1], fields[2], fields[3], fields[4], fields[5], fields[6] });
		}
		return rows;
	}

	std::vector<lambda_row> lambda_rows(const std::string& out) { return table_rows(out, false); }
	std::vector<lambda_row> iteration_rows(const std::string& log) { return table_rows(log, true); }

	//A criterion or GooF that a real run can print: finite, positive, nowhere near a blow-up
	bool physical_criterion(const std::string& s)
	{
		const double v = std::stod(s);
		return std::isfinite(v) && v > 0.0 && v < 1e3;
	}

	//P1_test_XCW.good, lambda = 0: the Hartree-Fock energy of P1_test.cif, which no reflection subset moves
	constexpr double golden_energy = -1961.925585145;
	constexpr int p1_nmo = 103;
	constexpr int p1_pairs = p1_nmo * (p1_nmo + 1) / 2;
	constexpr int p1_atoms = 23;
	const std::string common = "normal normal_conv params 177 basis_set sto-3g max_iter 100 charge 0 mult 1 ";

} //namespace

//------------------------------------------------------------------
//XCW runs on P1

//The FIXTURES_SETUP of xcw_p1_eri (SetIntegrationTestLocks.cmake): one SCF iteration is enough to have the
//stored two-electron integrals written
TEST(XcwScfTests, SharedEriBuild)
{
	if (p1_fixture().empty()) GTEST_SKIP() << "fixture tests/P1_test not found";
	if (shared_eri().empty()) GTEST_SKIP() << "NOS_TEST_SHARED_DIR is not set, so every run builds its own integrals";
	const auto dir = scratch_dir();
	std::filesystem::create_directories(shared_eri().parent_path());
	std::filesystem::remove(shared_eri());
	const p1_run run = run_on_p1(dir, "normal normal_conv params 177 basis_set sto-3g max_iter 1 charge 0 mult 1 f rhf start 0 step_size 0.01 end 0 save_eri eri.bin");
	ASSERT_TRUE(std::filesystem::exists(dir / "eri.bin")) << run.log;
	std::filesystem::copy_file(dir / "eri.bin", shared_eri());
	std::filesystem::remove_all(dir);
}

//lambda = 0 is a plain Hartree-Fock: the energy of P1_test_XCW.good, the .tscb/.wfn/.fchk
//of that step, damping and level shift switched off along the way; a second run lands on the
//same numbers
TEST(XcwScfTests, LambdaZeroIsTheHartreeFockEnergy)
{
	if (p1_fixture().empty()) GTEST_SKIP() << "fixture tests/P1_test not found";
	const auto dir = scratch_dir();
	const std::string settings = common + "f rhf start 0 step_size 0.01 end 0";

	const p1_run first = run_on_p1(dir, settings + read_shared_eri(dir));
	EXPECT_EQ(first.nr, 402);
	EXPECT_NE(first.out.find("XCW orbital basis set: sto-3g"), std::string::npos);
	EXPECT_EQ(first.out.find("I tensor"), std::string::npos) << first.out;
	EXPECT_NE(first.out.find("More detailed output in SCF.log file..."), std::string::npos);
	EXPECT_NE(first.out.find("Finished XCW fitting procedure."), std::string::npos);
	//Distance screening is a property of the structure and the basis, not of the hkl
	const int screened = static_cast<int>(number_after(first.out, "Screened out "));
	EXPECT_GT(screened, 0);
	EXPECT_LT(screened, p1_pairs);

	const auto rows = lambda_rows(first.out);
	ASSERT_EQ(rows.size(), 1u) << first.out;
	EXPECT_EQ(rows[0].lambda, "0.00000");
	EXPECT_NEAR(rows[0].d(rows[0].energy), golden_energy, 1e-6);
	EXPECT_TRUE(physical_criterion(rows[0].criterion)) << rows[0].criterion;
	EXPECT_TRUE(physical_criterion(rows[0].goof2)) << rows[0].goof2;
	EXPECT_EQ(rows[0].penalty, "0.000");
	//No perturbation at lambda = 0: the target quantity is the energy
	EXPECT_EQ(rows[0].quant, rows[0].energy);

	EXPECT_NE(first.log.find("Starting XCW SCF solver with lambda = 0.00000"), std::string::npos);
	EXPECT_NE(first.log.find("***SCF converged in "), std::string::npos);
	EXPECT_EQ(first.log.find("did not converge"), std::string::npos);
	EXPECT_GE(occurrences(first.log, "***Turned off damping***"), 1u);
	EXPECT_EQ(occurrences(first.log, "***Turned off level shift***"), 1u);
	EXPECT_NE(first.log.find("Creating .tscb file from converged SCF calculation..."), std::string::npos);
	EXPECT_NE(first.log.find("of two-electron integrals over"), std::string::npos);

	//The converged step's files: to_string(0.0) minus its dot, padded to seven digits
	EXPECT_TRUE(std::filesystem::exists(dir / "NA2_0000000.wfn"));
	EXPECT_TRUE(std::filesystem::exists(dir / "NA2_0000000.fchk"));
	//the F_calc table only with -print_F
	EXPECT_FALSE(std::filesystem::exists(dir / "NA2_0000000_Fcalc.txt"));
	ASSERT_TRUE(std::filesystem::exists(dir / "NA2_0000000.tscb"));
	{
		tsc_block<int, cdouble> tscb(dir / "NA2_0000000.tscb");
		EXPECT_EQ(tscb.scatterer_size(), static_cast<size_t>(p1_atoms));
		EXPECT_GE(tscb.reflection_size(), static_cast<size_t>(first.nr));
		ASSERT_GT(tscb.scatterer_size(), 0u);
		const cvec& sf = tscb.get_sf_for_scatterer(0);
		EXPECT_EQ(sf.size(), tscb.reflection_size());
		for (const cdouble& v : sf) {
			ASSERT_TRUE(std::isfinite(v.real()) && std::isfinite(v.imag()));
		}
	}

	const p1_run second = run_on_p1(dir, settings + read_shared_eri(dir));
	EXPECT_EQ(static_cast<int>(number_after(second.out, "Screened out ")), screened);
	const auto rows2 = lambda_rows(second.out);
	ASSERT_EQ(rows2.size(), 1u) << second.out;
	EXPECT_EQ(rows2[0].criterion, rows[0].criterion);
	EXPECT_EQ(rows2[0].goof2, rows[0].goof2);
	EXPECT_NEAR(rows2[0].d(rows2[0].energy), rows[0].d(rows[0].energy), 1e-8);
	EXPECT_NE(second.log.find("***SCF converged in "), std::string::npos);

	std::filesystem::remove_all(dir);
}

//At lambda = 0.01 the perturbation raises the energy a little above the Hartree-Fock minimum, and the
//target quantity is the energy plus lambda times the criterion
TEST(XcwScfTests, LambdaPointZeroOneRaisesTheEnergy)
{
	if (p1_fixture().empty()) GTEST_SKIP() << "fixture tests/P1_test not found";
	const auto dir = scratch_dir();
	const p1_run run = run_on_p1(dir, common + "f rhf start 0.01 step_size 0.01 end 0.01" + read_shared_eri(dir));
	const auto rows = lambda_rows(run.out);
	ASSERT_EQ(rows.size(), 1u) << run.out;
	EXPECT_EQ(rows[0].lambda, "0.01000");
	EXPECT_TRUE(std::filesystem::exists(dir / "NA2_0010000.tscb"));
	//quant = energy + lambda * criterion; the printed penalty is that product to 3 decimals
	EXPECT_NEAR(rows[0].d(rows[0].quant), rows[0].d(rows[0].energy) + 0.01 * rows[0].d(rows[0].criterion), 1e-3);
	EXPECT_NEAR(rows[0].d(rows[0].penalty), 0.01 * rows[0].d(rows[0].criterion), 1e-3);
	EXPECT_GT(rows[0].d(rows[0].energy), golden_energy);
	EXPECT_LT(rows[0].d(rows[0].energy), golden_energy + 0.5);
	EXPECT_TRUE(physical_criterion(rows[0].criterion)) << rows[0].criterion;
	EXPECT_TRUE(physical_criterion(rows[0].goof2)) << rows[0].goof2;
	EXPECT_NE(run.log.find("Starting XCW SCF solver with lambda = 0.01000"), std::string::npos);
	EXPECT_NE(run.log.find("***SCF converged in "), std::string::npos);
	std::filesystem::remove_all(dir);
}

//`guess_basis`: the first lambda starts from OCC's own Hartree-Fock in that basis. With the
//orbital basis itself as the guess basis the small SCF is the whole answer, so the log line
//carries the golden energy and the XCW loop has nothing left to do beyond recognising it.
TEST(XcwScfTests, OrbitalJsonOverridesMatchTheLibraryBasis)
{
	if (p1_fixture().empty()) GTEST_SKIP() << "P1 fixture not found";
	const auto dir = scratch_dir();
	const auto basis = BasisSetLibrary::get_basis_set("sto-3g");
	basis->write_occ_json(dir / "override basis.json", { 1, 6, 7, 8, 16, 17 });
	const auto baseline = run_on_p1(dir, common + "start 0 end 0" + read_shared_eri(dir));
	const std::string eri = read_shared_eri(dir);
	const auto overridden = run_on_p1(dir, common + "basis_overrides \"" + (dir / "override basis.json").string() + "\" start 0 end 0" + eri);
	const auto a = lambda_rows(baseline.out), b = lambda_rows(overridden.out);
	ASSERT_EQ(a.size(), 1u);
	ASSERT_EQ(b.size(), 1u);
	EXPECT_NEAR(a[0].d(a[0].energy), b[0].d(b[0].energy), 1e-8);
	EXPECT_NEAR(a[0].d(a[0].criterion), b[0].d(b[0].criterion), 1e-6);
	EXPECT_NE(overridden.log.find("_overrides_override basis.json"), std::string::npos);
	//The library basis's stored integrals are keyed on every exponent, coefficient and centre: read back, the
	//override is that basis to the last bit
	if (!eri.empty()) EXPECT_NE(overridden.log.find("read from shared_eri.bin"), std::string::npos) << overridden.log;
	std::filesystem::remove_all(dir);
}

TEST(XcwScfTests, SmallBasisGuessStartsFromConvergedDensity)
{
	if (p1_fixture().empty()) GTEST_SKIP() << "fixture tests/P1_test not found";
	const auto dir = scratch_dir();
	const p1_run run = run_on_p1(dir, common + "f rhf start 0 step_size 0.01 end 0 guess_basis sto-3g" + read_shared_eri(dir));

	const std::string key = "XCW: initial guess from a sto-3g Hartree-Fock (";
	ASSERT_NE(run.log.find(key), std::string::npos) << run.log;
	EXPECT_EQ(static_cast<int>(number_after(run.log, key)), p1_nmo);
	EXPECT_NEAR(number_after(run.log, "functions), E = "), golden_energy, 1e-5);
	EXPECT_NE(run.log.find(" Eh"), std::string::npos);

	ASSERT_NE(run.log.find("***SCF converged in "), std::string::npos) << run.log;
	//Iteration 1 can never converge (the quantity difference is measured against zero);
	//from a converged density the criteria are met within the next few
	const int iterations = static_cast<int>(number_after(run.log, "***SCF converged in "));
	EXPECT_GE(iterations, 2);
	EXPECT_LE(iterations, 8);

	const auto rows = lambda_rows(run.out);
	ASSERT_EQ(rows.size(), 1u) << run.out;
	EXPECT_NEAR(rows[0].d(rows[0].energy), golden_energy, 1e-6);
	EXPECT_TRUE(physical_criterion(rows[0].goof2)) << rows[0].goof2;
	EXPECT_TRUE(std::filesystem::exists(dir / "NA2_0000000.tscb"));
	std::filesystem::remove_all(dir);
}

//A closed shell in UHF is the RHF answer; the criterion column under `f2` is GooF(F2)
//squared. The unrestricted branches of the effective density, the level shift and the
//orbital gradient all run here, and the RHF run on the same subset is the reference.
TEST(XcwScfTests, UnrestrictedClosedShellAgainstF2MatchesRestricted)
{
	if (p1_fixture().empty()) GTEST_SKIP() << "fixture tests/P1_test not found";
	const auto dir = scratch_dir();
	const p1_run run = run_on_p1(dir, common + "f2 uhf start 0 step_size 0.01 end 0" + read_shared_eri(dir));

	ASSERT_NE(run.log.find("***SCF converged in "), std::string::npos) << run.log;
	EXPECT_EQ(occurrences(run.log, "***Turned off level shift***"), 1u);
	const auto rows = lambda_rows(run.out);
	ASSERT_EQ(rows.size(), 1u) << run.out;
	EXPECT_NEAR(rows[0].d(rows[0].criterion), rows[0].d(rows[0].goof2) * rows[0].d(rows[0].goof2), 1e-4 * rows[0].d(rows[0].goof2) + 1e-4);
	EXPECT_NEAR(rows[0].d(rows[0].energy), golden_energy, 1e-6);
	EXPECT_EQ(rows[0].quant, rows[0].energy);
	EXPECT_TRUE(std::filesystem::exists(dir / "NA2_0000000.tscb"));
	EXPECT_TRUE(std::filesystem::exists(dir / "NA2_0000000.fchk"));

	const p1_run rhf = run_on_p1(dir, common + "f2 rhf start 0 step_size 0.01 end 0" + read_shared_eri(dir));
	const auto rows_rhf = lambda_rows(rhf.out);
	ASSERT_EQ(rows_rhf.size(), 1u) << rhf.out;
	EXPECT_EQ(rows_rhf[0].goof2, rows[0].goof2);
	EXPECT_NEAR(rows_rhf[0].d(rows_rhf[0].energy), rows[0].d(rows[0].energy), 1e-6);
	std::filesystem::remove_all(dir);
}

//`fast_conv` runs DIIS alone: nothing to switch off, so the log never says so, and the
//weighted criterion is announced and used. The preset sets neither damping nor shift;
//SCF_settings defaults both to 0.
TEST(XcwScfTests, FastConvWeightedRunsWithoutDampingOrShift)
{
	if (p1_fixture().empty()) GTEST_SKIP() << "fixture tests/P1_test not found";
	const auto dir = scratch_dir();
	const p1_run run = run_on_p1(dir, "normal fast_conv params 177 basis_set sto-3g max_iter 100 charge 0 mult 1 f rhf start 0 step_size 0.01 end 0 weighted" + read_shared_eri(dir));

	EXPECT_NE(run.out.find("XCW: fitting against the 1/|H|^2-weighted residual self-energy criterion"), std::string::npos) << run.out;
	EXPECT_NE(run.log.find("XCW: fitting against the 1/|H|^2-weighted residual self-energy criterion"), std::string::npos);
	ASSERT_NE(run.log.find("***SCF converged in "), std::string::npos) << run.log;
	EXPECT_EQ(occurrences(run.log, "Turned off"), 0u);
	EXPECT_EQ(occurrences(run.log, "Decreased damping"), 0u);
	const auto rows = lambda_rows(run.out);
	ASSERT_EQ(rows.size(), 1u) << run.out;
	//The weighted criterion is a different number from GooF(F2), not a copy of it
	EXPECT_TRUE(physical_criterion(rows[0].criterion)) << rows[0].criterion;
	EXPECT_TRUE(physical_criterion(rows[0].goof2)) << rows[0].goof2;
	EXPECT_NE(rows[0].criterion, rows[0].goof2);
	EXPECT_NEAR(rows[0].d(rows[0].energy), golden_energy, 1e-6);
	std::filesystem::remove_all(dir);
}

//Five iterations at lambda = 0.01 from the core guess, twice: the same five SCF.log rows to
//the last printed digit (the grid contraction, the Fock build and DIIS are all in
//them), energies within a few Hartree of the Hartree-Fock minimum and lower at the end
//than at the start, criteria positive. No converged result is asked for.
TEST(XcwScfTests, FiveIterationsAreDeterministicAndPhysical)
{
	if (p1_fixture().empty()) GTEST_SKIP() << "fixture tests/P1_test not found";
	const auto dir = scratch_dir();
	const std::string settings = "normal normal_conv params 177 basis_set sto-3g max_iter 5 charge 0 mult 1 f rhf start 0.01 step_size 0.01 end 0.01" + read_shared_eri(dir);

	const p1_run first = run_on_p1(dir, settings);
	const auto rows = iteration_rows(first.log);
	ASSERT_EQ(rows.size(), 5u) << first.log;
	EXPECT_NE(first.log.find("***SCF did not converge***"), std::string::npos);
	EXPECT_TRUE(lambda_rows(first.out).empty()) << first.out;
	for (size_t i = 0; i < rows.size(); i++) {
		EXPECT_EQ(rows[i].lambda, std::to_string(i + 1));
		const double e = rows[i].d(rows[i].energy);
		EXPECT_TRUE(std::isfinite(e)) << rows[i].energy;
		EXPECT_GT(e, golden_energy - 5.0) << rows[i].energy;
		EXPECT_LT(e, golden_energy + 5.0) << rows[i].energy;
		EXPECT_TRUE(physical_criterion(rows[i].criterion)) << rows[i].criterion;
		EXPECT_TRUE(physical_criterion(rows[i].goof2)) << rows[i].goof2;
		//quant = energy + lambda * criterion, both printed to three decimals
		EXPECT_NEAR(rows[i].d(rows[i].penalty), 0.01 * rows[i].d(rows[i].criterion), 1e-3);
		EXPECT_NEAR(rows[i].d(rows[i].quant), e + rows[i].d(rows[i].penalty), 1e-3);
	}
	EXPECT_LT(rows[4].d(rows[4].energy), rows[0].d(rows[0].energy));

	const p1_run second = run_on_p1(dir, settings);
	const auto rows2 = iteration_rows(second.log);
	ASSERT_EQ(rows2.size(), 5u) << second.log;
	for (size_t i = 0; i < rows.size(); i++) {
		EXPECT_EQ(rows2[i].criterion, rows[i].criterion) << "iteration " << i + 1;
		EXPECT_EQ(rows2[i].goof2, rows[i].goof2) << "iteration " << i + 1;
		EXPECT_EQ(rows2[i].energy, rows[i].energy) << "iteration " << i + 1;
		EXPECT_EQ(rows2[i].quant, rows[i].quant) << "iteration " << i + 1;
	}
	std::filesystem::remove_all(dir);
}

//One iteration can never satisfy the perturbed-energy criterion (it is measured against
//zero), so max_iter 1 fails the first lambda: the five-line verdict goes to the log, no
//tscb is written, the scan stops before lambda 0.01 and the run still ends properly.
TEST(XcwScfTests, MaxIterOneStopsScanWithoutTscb)
{
	if (p1_fixture().empty()) GTEST_SKIP() << "fixture tests/P1_test not found";
	const auto dir = scratch_dir();
	const p1_run run = run_on_p1(dir, "normal normal_conv params 177 basis_set sto-3g max_iter 1 charge 0 mult 1 f rhf start 0 step_size 0.01 end 0.01" + read_shared_eri(dir));

	EXPECT_NE(run.log.find("***SCF did not converge***"), std::string::npos) << run.log;
	EXPECT_EQ(run.log.find("***SCF converged in "), std::string::npos);
	EXPECT_NE(run.log.find("NOT CONVERGED for perturbed energy: "), std::string::npos);
	EXPECT_NE(run.log.find(" for DIIS error: "), std::string::npos);
	EXPECT_NE(run.log.find(" for orbital gradient: "), std::string::npos);
	EXPECT_NE(run.log.find(" for maximum difference in density matrix: "), std::string::npos);
	EXPECT_NE(run.log.find(" for RMSD of density matrix: "), std::string::npos);
	EXPECT_EQ(occurrences(run.log, "CONVERGED"), 5u);
	EXPECT_EQ(occurrences(run.log, "Starting XCW SCF solver with lambda = "), 1u);
	EXPECT_EQ(run.log.find("Creating .tscb file"), std::string::npos);
	EXPECT_NE(run.log.find("\n 1\t"), std::string::npos);
	EXPECT_EQ(run.log.find("\n 2\t"), std::string::npos);

	//step 0 has no converged neighbour to halve the step towards, so the message names the
	//iteration cap instead of a continuation step
	const std::string stop = "XCW: unable to converge lambda 0.00000000 in 1 SCF iterations (raise max_iter or loosen the criteria); stopping scan.";
	EXPECT_NE(run.out.find(stop), std::string::npos) << run.out;
	EXPECT_NE(run.log.find(stop), std::string::npos);
	EXPECT_EQ(run.out.find("retrying lambda"), std::string::npos);
	EXPECT_NE(run.out.find("Finished XCW fitting procedure."), std::string::npos);
	EXPECT_TRUE(lambda_rows(run.out).empty()) << run.out;
	EXPECT_FALSE(std::filesystem::exists(dir / "NA2_0000000.tscb"));
	EXPECT_FALSE(std::filesystem::exists(dir / "NA2_0010000.tscb"));
	std::filesystem::remove_all(dir);
}

//------------------------------------------------------------------
//The grid against closed-form numbers: H2 in a cubic P1 cell, STO-3G, U = 0, so the structure factor of
//a density matrix with one element set is that element times
//I_mn(h) = integral phi_m phi_n exp(i q.r) d3r with q = 2 pi/a (h k l), a sum over the
//nine primitive pairs of the Gaussian product theorem, (pi/p)^(3/2) exp(-ab/p |A-B|^2)
//exp(-|q|^2/(4p)) exp(i q.P), p = a + b, P = (aA + bB)/p, with the contraction normalised
//as the basis loader does. What the grid quadrature, the Hirshfeld partition and the phase
//and DW factors produce with the direct-sum transform has to agree with that to the grid's
//accuracy (largest at 0 0 4), F(-h) has to be the conjugate of F(h) and no diagonal element
//can exceed the overlap of a normalised function with itself.
namespace {

	struct h2_analytic {
		//STO-3G H: exponents and coefficients for normalised primitives (basis_data.cpp)
		const double a[3] = { 3.425250914, 0.6239137298, 0.168855404 };
		const double c[3] = { 0.1543289673, 0.5353281423, 0.4446345422 };
		double norm = 1.0;
		h2_analytic()
		{
			double s = 0.0;
			for (int k = 0; k < 3; k++)
				for (int l = 0; l < 3; l++)
					s += c[k] * c[l] * std::pow(2.0 * std::sqrt(a[k] * a[l]) / (a[k] + a[l]), 1.5);
			norm = 1.0 / std::sqrt(s);
		}
		//A and B in bohr, q in 1/bohr
		cdouble I(const vec& A, const vec& B, const vec& q) const
		{
			const double AB2 = (A[0] - B[0]) * (A[0] - B[0]) + (A[1] - B[1]) * (A[1] - B[1]) + (A[2] - B[2]) * (A[2] - B[2]);
			const double q2 = q[0] * q[0] + q[1] * q[1] + q[2] * q[2];
			cdouble sum = 0.0;
			for (int k = 0; k < 3; k++) {
				for (int l = 0; l < 3; l++) {
					const double p = a[k] + a[l];
					const double N = std::pow(2.0 * a[k] / constants::PI, 0.75) * std::pow(2.0 * a[l] / constants::PI, 0.75);
					double qP = 0.0;
					for (int d = 0; d < 3; d++) qP += q[d] * (a[k] * A[d] + a[l] * B[d]) / p;
					sum += c[k] * c[l] * N * std::pow(constants::PI / p, 1.5) * std::exp(-a[k] * a[l] / p * AB2 - q2 / (4.0 * p)) * std::polar(1.0, qP);
				}
			}
			return sum * norm * norm;
		}
	};

} //namespace

TEST(XcwScfTests, GridOfH2MatchesAnalyticGaussianTransform)
{
	const auto dir = scratch_dir();
	constexpr double a_ang = 6.0;
	constexpr double dz = 0.37 / a_ang;   //H-H 0.74 A along z, centred in the cell
	const auto cif = dir / "h2.cif";
	{
		std::ofstream o(cif);
		o << "data_h2\n_cell_length_a 6.0\n_cell_length_b 6.0\n_cell_length_c 6.0\n_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n_cell_volume 216.0\n"
			<< "loop_\n_space_group_symop_operation_xyz\n'x, y, z'\n"
			<< "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n";
		o << std::setprecision(12) << "H1 H 0.5 0.5 " << 0.5 - dz << "\nH2 H 0.5 0.5 " << 0.5 + dz << "\n";
	}
	//Every h,k,l in -1..1 but 000, plus four farther ones and their Friedel mates
	std::set<i3> hkl;
	for (int h = -1; h <= 1; h++)
		for (int k = -1; k <= 1; k++)
			for (int l = -1; l <= 1; l++)
				if (h != 0 || k != 0 || l != 0) hkl.insert({ h, k, l });
	for (const i3& r : { i3{ 2, 0, 0 }, i3{ 0, 3, 1 }, i3{ 2, 2, 2 }, i3{ 0, 0, 4 } }) {
		hkl.insert(r);
		hkl.insert({ -r[0], -r[1], -r[2] });
	}
	const auto hkl_path = dir / "h2.hkl";
	{
		std::ofstream o(hkl_path);
		for (const i3& r : hkl) o << std::setw(4) << r[0] << std::setw(4) << r[1] << std::setw(4) << r[2] << "  100.00    1.00\n";
	}
	struct direct_transform {
		std::string old; bool had = false;
		direct_transform() { if (const char* e = tuning("NOSPHERA2_XCW_GRID_DIRECT")) { old = e; had = true; } set_tuning("NOSPHERA2_XCW_GRID_DIRECT", "1"); }
		~direct_transform() { set_tuning("NOSPHERA2_XCW_GRID_DIRECT", had ? old.c_str() : nullptr); }
	} pin;
	options opt = xcw_options(dir, cif, hkl_path, "normal normal_conv params 1 basis_set sto-3g max_iter 100 charge 0 mult 1 f rhf start 0 step_size 0.01 end 0", true, 3);
	std::vector<cvec> got;
	std::vector<i3> rows;
	{
		cwd_guard cwd(dir);
		cout_capture capture;
		opt.loadXCWsettings();
		structure_factors sf(opt);
		const occ::core::Molecule mol = SCF_solver::setup_SCF_mol(sf.asym_atoms, 0, 1);
		sf.set_DW();
		sf.set_phases();
		sf.set_anom();
		sf.setup_grid(sf.orbital_basis->to_AOBasis(mol.atoms()));
		ASSERT_EQ(sf.model_data.nmo, 2);
		rows.assign(sf.scatter_data.hkl.begin(), sf.scatter_data.hkl.end());
		got.assign(rows.size(), cvec(3));
		//F = anomalous + Sum_{mu<=nu} (2|4) D I: D(0,0) = 1/2 gives I_00, D(0,1) = D(1,0) = 1/4 I_01, D(1,1) = 1/2 I_11
		const double D_entries[3][3] = { { 0.5, 0, 0 }, { 0, 0.25, 0 }, { 0, 0, 0.5 } };
		for (int t = 0; t < 3; t++) {
			dMatrix2 D(2, 2);
			D(0, 0) = D_entries[t][0];
			D(0, 1) = D(1, 0) = D_entries[t][1];
			D(1, 1) = D_entries[t][2];
			sf.calc_F_calc(D);
			for (size_t r = 0; r < rows.size(); r++) got[r][t] = sf.scatter_data.F_calc[r] - sf.scatter_data.anom_correction[r];
		}
	}
	const h2_analytic ref;
	const double a_bohr = constants::ang2bohr(a_ang);
	const vec A = { 0.5 * a_bohr, 0.5 * a_bohr, (0.5 - dz) * a_bohr };
	const vec B = { 0.5 * a_bohr, 0.5 * a_bohr, (0.5 + dz) * a_bohr };
	const vec centres[2] = { A, B };
	const int pair_mu[3] = { 0, 0, 1 }, pair_nu[3] = { 0, 1, 1 };
	double max_dev = 0.0;
	for (size_t r = 0; r < rows.size(); r++) {
		const i3& h = rows[r];
		const vec q = { constants::TWO_PI / a_bohr * h[0], constants::TWO_PI / a_bohr * h[1], constants::TWO_PI / a_bohr * h[2] };
		for (int t = 0; t < 3; t++) {
			const cdouble expected = ref.I(centres[pair_mu[t]], centres[pair_nu[t]], q);
			const double dev = std::abs(got[r][t] - expected);
			max_dev = std::max(max_dev, dev);
			EXPECT_LT(dev, 5e-5) << "hkl " << h[0] << " " << h[1] << " " << h[2] << " pair " << t << ": got " << got[r][t] << " expected " << expected;
		}
		EXPECT_LE(std::abs(got[r][0]), 1.0 + 1e-6);
		EXPECT_LE(std::abs(got[r][2]), 1.0 + 1e-6);
	}
	std::cerr << "H2 grid: max |grid - analytic| = " << max_dev << " over " << rows.size() << " reflections" << std::endl;
	//Friedel mates are complex conjugates: the same grid points and weights enter, only the phase flips
	for (size_t r = 0; r < rows.size(); r++) {
		const size_t mate = std::find(rows.begin(), rows.end(), i3{ -rows[r][0], -rows[r][1], -rows[r][2] }) - rows.begin();
		ASSERT_LT(mate, rows.size());
		for (int t = 0; t < 3; t++) {
			EXPECT_NEAR(got[r][t].real(), got[mate][t].real(), 1e-10);
			EXPECT_NEAR(got[r][t].imag(), -got[mate][t].imag(), 1e-10);
		}
	}
	std::filesystem::remove_all(dir);
}

//`i_sigma <x>` keeps only I/sigma(I) >= x in the fit: the run says how many, the lambda
//table adds Crit(all) and R1(all) behind the seven fit-set columns, and the lambda = 0
//energy does not depend on which reflections are scored.
TEST(XcwScfTests, ISigmaCutoffShrinksFitSetAndReportsAllReflections)
{
	if (p1_fixture().empty()) GTEST_SKIP() << "fixture tests/P1_test not found";
	const auto dir = scratch_dir();
	const p1_run run = run_on_p1(dir, common + "f rhf start 0 step_size 0.01 end 0 i_sigma 20" + read_shared_eri(dir));

	const auto rows = lambda_rows(run.out);
	ASSERT_EQ(rows.size(), 1u) << run.out;
	EXPECT_NEAR(rows[0].d(rows[0].energy), golden_energy, 1e-6);

	const size_t at = run.out.find("XCW: I/sigma(I) >= 20 ");
	ASSERT_NE(at, std::string::npos) << run.out;
	int n_fit = 0, nr = 0;
	std::istringstream(run.out.substr(run.out.find("): ", at) + 3)) >> n_fit;
	std::istringstream(run.out.substr(run.out.find(" of ", at) + 4)) >> nr;
	EXPECT_GT(n_fit, 0);
	EXPECT_LT(n_fit, nr) << run.out;

	//the nine columns of the converged row: fit-set criterion and R1(gt) differ from Crit(all) and R1(all)
	const size_t row = run.out.find("\n 0.00000\t");
	ASSERT_NE(row, std::string::npos) << run.out;
	std::vector<std::string> fields;
	std::istringstream cols(run.out.substr(row + 1, run.out.find('\n', row + 1) - row - 1));
	for (std::string f; cols >> f;) fields.push_back(f);
	ASSERT_EQ(fields.size(), 9u) << run.out;
	EXPECT_NE(fields[7], fields[1]);
	EXPECT_NE(fields[8], fields[3]);
	std::filesystem::remove_all(dir);
}

//`slow_conv` damps every step, the unperturbed first one included, and still reaches the golden lambda = 0 energy
TEST(XcwScfTests, SlowConvDampsFromTheFirstStep)
{
	if (p1_fixture().empty()) GTEST_SKIP() << "fixture tests/P1_test not found";
	const auto dir = scratch_dir();
	const p1_run run = run_on_p1(dir, common + "f rhf start 0 step_size 0.01 end 0.01 slow_conv" + read_shared_eri(dir));
	EXPECT_EQ(run.out.find("normal schedule"), std::string::npos) << run.out;
	const auto rows = lambda_rows(run.out);
	ASSERT_EQ(rows.size(), 2u) << run.out;
	EXPECT_NEAR(rows[0].d(rows[0].energy), golden_energy, 1e-6);
	EXPECT_GT(rows[1].d(rows[1].energy), golden_energy);
	std::filesystem::remove_all(dir);
}

//The scale is the minimiser of the criterion the SCF descends, so at the scale the dump
//prints d(chi^2)/dk = Sum (k|Fc| - |Fo|)|Fc|/sigma^2 vanishes over the fit set (I/sigma(I) >= 2,
//i.e. |F|/sigma(F) >= 4). The unweighted fit it replaced left the ratio at 1.2e-3 on this
//subset; the four-decimal dump columns leave ~1e-5.
TEST(XcwScfTests, ScaleIsStationaryForTheCriterion)
{
	if (p1_fixture().empty()) GTEST_SKIP() << "fixture tests/P1_test not found";
	const auto dir = scratch_dir();
	const p1_run run = run_on_p1(dir, common + "f rhf start 0 step_size 0.01 end 0" + read_shared_eri(dir), true);
	std::istringstream dump(read_text(dir / "NA2_0000000_Fcalc.txt"));
	std::string line;
	double numerator = 0.0, denominator = 0.0;
	int in_fit = 0;
	while (std::getline(dump, line)) {
		if (line.empty() || line[0] == '#') continue;
		std::istringstream cols(line);
		int h, k, l;
		double fo, sig, kfc, phase;
		cols >> h >> k >> l >> fo >> sig >> kfc >> phase;
		if (!(fo > 0.0) || fo / (2.0 * sig) < 2.0) continue;
		in_fit++;
		numerator += (kfc - fo) * kfc / (sig * sig);
		denominator += kfc * kfc / (sig * sig);
	}
	ASSERT_GT(in_fit, 100) << run.out;
	EXPECT_LT(std::abs(numerator / denominator), 1e-4) << "d(chi^2)/dk is not zero at the printed scale: " << numerator / denominator;
	std::filesystem::remove_all(dir);
}

//`check_hessian` hands the step to the second-order solver once the DIIS error is below 1e-2; from
//there the trust-region augmented-Hessian steps on the orbital rotations have to reach the
//same lambda = 0 minimum as the Roothaan/DIIS iterations, with micro-iterations reported
TEST(XcwScfTests, SecondOrderStepsReachTheGoldenEnergy)
{
	if (p1_fixture().empty()) GTEST_SKIP() << "fixture tests/P1_test not found";
	const auto dir = scratch_dir();
	const p1_run run = run_on_p1(dir, common + "f rhf start 0 step_size 0.01 end 0 check_hessian" + read_shared_eri(dir));
	EXPECT_NE(run.log.find("second-order steps on the orbital rotations from here"), std::string::npos) << run.log;
	EXPECT_NE(run.log.find("TRAH: "), std::string::npos) << run.log;
	EXPECT_NE(run.log.find("***SCF converged in"), std::string::npos) << run.log;
	const auto rows = lambda_rows(run.out);
	ASSERT_EQ(rows.size(), 1u) << run.out;
	EXPECT_NEAR(rows[0].d(rows[0].energy), golden_energy, 1e-6);
	std::filesystem::remove_all(dir);
}

//The exact Hessian-vector product behind the second-order steps (Fock response plus the
//response of the perturbation with the scale moving) against a central finite difference of
//the gradient, at a lambda where the perturbation is on, for both spin treatments; the
//grid in double so the difference quotient is not limited by fp32 noise in the gradient
TEST(XcwScfTests, HessianVectorProductMatchesFiniteDifference)
{
	if (p1_fixture().empty()) GTEST_SKIP() << "fixture tests/P1_test not found";
	for (const std::string kind : { "f rhf", "f2 uhf" }) {
		const auto dir = scratch_dir();
		const p1_run run = run_on_p1(dir, common + kind + " start 0.02 step_size 0.01 end 0.02 i_double check_hessian" + read_shared_eri(dir));
		const size_t at = run.log.find("Hessian check: |Hv - FD| / |FD| = ");
		ASSERT_NE(at, std::string::npos) << run.log;
		const double rel = std::stod(run.log.substr(at + std::string("Hessian check: |Hv - FD| / |FD| = ").size()));
		EXPECT_LT(rel, 1e-4) << kind << ": " << run.log.substr(at, 160);
		std::cout << kind << ": " << run.log.substr(at, 130) << std::endl;
		std::filesystem::remove_all(dir);
	}
}
