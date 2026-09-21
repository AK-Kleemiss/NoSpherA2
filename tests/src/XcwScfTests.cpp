//XCW SCF driver and I tensor storage: run_XCW_fitting on the P1 fixture through the
//settings file, observed through XCW.log, stdout and the files a run leaves behind, plus
//the i_tensor_file container on its own.
//
//Every run uses the sto-3g basis of tests/P1_test at a single lambda on every 8th
//reflection of P1_test.hkl (402 of 3215). The I tensor costs reflections x grid points and
//was 10 of the 11 s of a full-hkl run; an SCF iteration is 9 ms. The GooF values of the
//subset are nobody's golden, so the runs check determinism (two runs, one answer), physical
//sanity (finite, bounded, descending) and the one number the subset cannot move: at lambda
//= 0 the perturbation is zero and the converged energy is the Hartree-Fock energy of
//P1_test_XCW.good, -1961.923538820 Eh. The last test builds the tensor for H2 in a cubic
//cell and compares every element with the analytic Fourier transform of the Gaussians.
#include "pch.h"
#include <gtest/gtest.h>

#include "core/convenience.h"
#include "core/XCW.h"
#include "core/i_tensor_stream.h"
#include "core/tsc_block.h"

#include <cmath>
#include <complex>
#include <cstdint>
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
		if (!std::filesystem::exists(fixture / "P1_test_NA2.cif") || !std::filesystem::exists(fixture / "P1_test.hkl")) {
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
		std::string out;    //what run_XCW_fitting and the constructor printed
		std::string log;    //XCW.log
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

	//One XCW run in dir on cif and hkl with the given settings file, no GPU and no timing
	//lines; the anomalous dispersion file is passed when it is given
	p1_run run_xcw(const std::filesystem::path& dir, const std::filesystem::path& cif, const std::filesystem::path& hkl,
		const std::filesystem::path& anom, const std::string& settings_text, const bool double_tensor = false, const int accuracy = 2)
	{
		const auto settings = dir / "settings.txt";
		std::ofstream(settings) << settings_text;
		options opt;
		opt.xcw_settings_path = settings;
		opt.cif = std::filesystem::absolute(cif);
		opt.hkl = std::filesystem::absolute(hkl);
		if (!anom.empty()) {
			opt.anom_disp_path = std::filesystem::absolute(anom);
		}
		opt.do_XCW = true;
		opt.use_gpu = false;
		opt.cpu_itensor_fp32 = !double_tensor;
		opt.accuracy = accuracy;
		opt.no_date = true;
		//As the driver does before constructing XCW: the tscb of a converged step is built
		//from CIF disorder group 0 (NoSpherA2.cpp, do_XCW branch)
		opt.groups[0].push_back(0);
		p1_run r;
		{
			cwd_guard cwd(dir);
			cout_capture capture;
			{
				XCW x(opt);
				x.run_XCW_fitting();
			}
			r.out = capture.buf.str();
		}
		r.log = read_text(dir / "XCW.log");
		return r;
	}

	//The P1 fixture on the reflection subset, written into dir on the first call
	p1_run run_on_p1(const std::filesystem::path& dir, const std::string& settings_text, const bool with_anom)
	{
		const auto fixture = p1_fixture();
		const auto hkl = dir / "subset.hkl";
		const int nr = write_subset_hkl(hkl);
		p1_run r = run_xcw(dir, fixture / "P1_test_NA2.cif", hkl, with_anom ? fixture / "anom_disp.txt" : std::filesystem::path{}, settings_text);
		r.nr = nr;
		return r;
	}

	//One line of the lambda table run_XCW_fitting prints per converged step:
	//lambda(5) criterion(3) GooF2(3) R1(4) energy(9) lambda*criterion(3) quant(9), tab separated.
	//An XCW.log iteration row has the same seven columns behind a leading tab, with the
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
			if (line.empty() || (line[0] == '\t') != iterations) continue;
			const size_t first = iterations ? 1 : 0;
			if (line.size() <= first || !std::isdigit(static_cast<unsigned char>(line[first])) || line.find('\t', first) == std::string::npos) {
				continue;
			}
			std::vector<std::string> fields;
			std::istringstream cols(line);
			std::string f;
			while (std::getline(cols, f, '\t')) {
				if (!f.empty()) fields.push_back(f);
			}
			if (fields.size() < 7) {
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

	//P1_test_XCW.good, lambda = 0: the Hartree-Fock energy, which no reflection subset moves
	constexpr double golden_energy = -1961.923538820;
	constexpr int p1_nmo = 103;
	constexpr int p1_pairs = p1_nmo * (p1_nmo + 1) / 2;
	constexpr int p1_atoms = 23;
	const std::string common = "normal normal_conv params 177 basis_set sto-3g max_iter 100 charge 0 mult 1 ";

	//------------------------------------------------------------------
	//A tiny tensor for the container tests: 3 reflections, 4 of the 10 pairs of a 4-function
	//basis, values that single precision holds exactly

	constexpr int tiny_nr = 3;
	constexpr int tiny_nmo = 4;
	const ivec tiny_mu = { 0, 0, 1, 3 };
	const ivec tiny_nu = { 0, 2, 3, 3 };
	constexpr int64_t magic_compact = 0x4E4132495F544E32LL;
	constexpr int64_t magic_full = 0x4E4132495F54454ELL;

	cdouble tiny_value(const int r, const size_t i)
	{
		return cdouble(10.0 * r + static_cast<double>(i), -(r + 0.5 * static_cast<double>(i)));
	}

	//Written from double blocks into a file of either element type
	void write_tiny(const std::filesystem::path& p, const bool single)
	{
		i_tensor_file f;
		f.create(p, tiny_nr, tiny_nmo, tiny_mu, tiny_nu, single);
		cvec block(tiny_mu.size());
		//Out of order on purpose: the writer seeks by index
		for (const int r : { 2, 0, 1 }) {
			for (size_t i = 0; i < block.size(); i++) block[i] = tiny_value(r, i);
			f.write_block(r, block.data());
		}
		f.finish_write();
		f.close();
	}

	//A file from raw parts: header, pair lists, then payload bytes
	void write_raw(const std::filesystem::path& p, const int64_t h[5], const ivec& mu, const ivec& nu, const size_t payload_bytes)
	{
		std::ofstream out(p, std::ios::binary);
		out.write(reinterpret_cast<const char*>(h), 5 * sizeof(int64_t));
		out.write(reinterpret_cast<const char*>(mu.data()), static_cast<std::streamsize>(mu.size() * sizeof(int)));
		out.write(reinterpret_cast<const char*>(nu.data()), static_cast<std::streamsize>(nu.size() * sizeof(int)));
		const std::vector<char> zeros(payload_bytes, 0);
		out.write(zeros.data(), static_cast<std::streamsize>(zeros.size()));
	}

	std::string open_error(const std::filesystem::path& p)
	{
		i_tensor_file f;
		try {
			f.open(p, 1);
		}
		catch (const std::runtime_error& e) {
			return e.what();
		}
		return "";
	}

} //namespace

//------------------------------------------------------------------
//i_tensor_file

TEST(XcwScfTensorFileTests, SizesFollowElementType)
{
	EXPECT_EQ(i_tensor_file::block_bytes(4, true), 32u);
	EXPECT_EQ(i_tensor_file::block_bytes(4, false), 64u);
	EXPECT_EQ(i_tensor_file::total_bytes(3, 4, true), 96u);
	EXPECT_EQ(i_tensor_file::total_bytes(3, 4, false), 192u);
	//size_t arithmetic: 20000 reflections of 200000 pairs pass 2^31 elements
	EXPECT_EQ(i_tensor_file::total_bytes(20000, 200000, false), 20000ull * 200000ull * 16ull);
}

TEST(XcwScfTensorFileTests, DoubleRoundTripThroughWindow)
{
	const auto dir = scratch_dir();
	const auto p = dir / "tiny.bin";
	write_tiny(p, false);
	//header 40 B + 2*4 ints 32 B + 3 blocks of 4 cdouble 192 B
	EXPECT_EQ(std::filesystem::file_size(p), 40u + 32u + 192u);

	size_t kept = 0;
	bool single = true;
	EXPECT_TRUE(i_tensor_file::matches(p, tiny_nr, tiny_nmo, kept, single));
	EXPECT_EQ(kept, 4u);
	EXPECT_FALSE(single);

	i_tensor_file f;
	f.open(p, 2);
	EXPECT_EQ(f.nr(), tiny_nr);
	EXPECT_EQ(f.nmo(), tiny_nmo);
	EXPECT_EQ(f.kept(), 4u);
	EXPECT_FALSE(f.single());
	EXPECT_EQ(f.window_blocks(), 2u);
	EXPECT_EQ(f.pair_mu(), tiny_mu);
	EXPECT_EQ(f.pair_nu(), tiny_nu);
	EXPECT_EQ(f.path(), p);

	f.load(0, 2);
	for (int r = 0; r < 2; r++) {
		for (size_t i = 0; i < 4; i++) {
			EXPECT_EQ(f.block(r)[i], tiny_value(r, i)) << "r " << r << " i " << i;
		}
	}
	f.load(2, 3);
	for (size_t i = 0; i < 4; i++) {
		EXPECT_EQ(f.block(2)[i], tiny_value(2, i)) << "i " << i;
	}
	//An empty range is allowed, a range past the window or the file is not
	EXPECT_NO_THROW(f.load(1, 1));
	EXPECT_THROW(f.load(0, 3), std::runtime_error);
	EXPECT_THROW(f.load(-1, 0), std::runtime_error);
	EXPECT_THROW(f.load(2, 4), std::runtime_error);
	EXPECT_THROW(f.load(2, 1), std::runtime_error);
	f.close();
	std::filesystem::remove_all(dir);
}

TEST(XcwScfTensorFileTests, SingleFileNarrowsDoubleBlocks)
{
	const auto dir = scratch_dir();
	const auto p = dir / "tiny32.bin";
	write_tiny(p, true);
	EXPECT_EQ(std::filesystem::file_size(p), 40u + 32u + 96u);

	size_t kept = 0;
	bool single = false;
	EXPECT_TRUE(i_tensor_file::matches(p, tiny_nr, tiny_nmo, kept, single));
	EXPECT_EQ(kept, 4u);
	EXPECT_TRUE(single);

	i_tensor_file f;
	f.open(p, 3);
	EXPECT_TRUE(f.single());
	f.load(0, 3);
	for (int r = 0; r < 3; r++) {
		for (size_t i = 0; i < 4; i++) {
			const cdouble v = tiny_value(r, i);
			EXPECT_EQ(f.block32(r)[i], std::complex<float>(static_cast<float>(v.real()), static_cast<float>(v.imag()))) << "r " << r << " i " << i;
		}
	}
	f.close();
	std::filesystem::remove_all(dir);
}

TEST(XcwScfTensorFileTests, DoubleFileWidensSingleBlocks)
{
	const auto dir = scratch_dir();
	const auto p = dir / "widen.bin";
	{
		i_tensor_file f;
		f.create(p, tiny_nr, tiny_nmo, tiny_mu, tiny_nu, false);
		std::vector<std::complex<float>> block(4);
		for (int r = 0; r < tiny_nr; r++) {
			for (size_t i = 0; i < 4; i++) {
				const cdouble v = tiny_value(r, i);
				block[i] = std::complex<float>(static_cast<float>(v.real()), static_cast<float>(v.imag()));
			}
			f.write_block(r, block.data());
		}
		f.finish_write();
	}
	i_tensor_file f;
	f.open(p, 1);
	EXPECT_FALSE(f.single());
	for (int r = 0; r < tiny_nr; r++) {
		f.load(r, r + 1);
		for (size_t i = 0; i < 4; i++) {
			EXPECT_EQ(f.block(r)[i], tiny_value(r, i)) << "r " << r << " i " << i;
		}
	}
	f.close();
	std::filesystem::remove_all(dir);
}

TEST(XcwScfTensorFileTests, WindowClampsToOneAndToNr)
{
	const auto dir = scratch_dir();
	const auto p = dir / "window.bin";
	write_tiny(p, false);
	i_tensor_file f;
	f.open(p, 0);
	EXPECT_EQ(f.window_blocks(), 1u);
	f.set_window(100);
	EXPECT_EQ(f.window_blocks(), 3u);
	EXPECT_NO_THROW(f.load(0, 3));
	f.set_window(1);
	EXPECT_EQ(f.window_blocks(), 1u);
	EXPECT_THROW(f.load(0, 2), std::runtime_error);
	f.close();
	//close() twice and a fresh object are both fine
	f.close();
	std::filesystem::remove_all(dir);
}

TEST(XcwScfTensorFileTests, MatchesRefusesOtherShapesAndLayouts)
{
	const auto dir = scratch_dir();
	const auto good = dir / "good.bin";
	write_tiny(good, false);
	size_t kept = 99;
	bool single = true;
	EXPECT_FALSE(i_tensor_file::matches(good, tiny_nr + 1, tiny_nmo, kept, single));
	EXPECT_EQ(kept, 0u);
	EXPECT_FALSE(single);
	EXPECT_FALSE(i_tensor_file::matches(good, tiny_nr, tiny_nmo + 1, kept, single));
	EXPECT_FALSE(i_tensor_file::matches(dir / "absent.bin", tiny_nr, tiny_nmo, kept, single));

	//The previous layout, every pair in double: refused with a note on stdout
	{
		const int64_t h[5] = { magic_full, tiny_nr, tiny_nmo, 10, 16 };
		write_raw(dir / "full.bin", h, ivec(10, 0), ivec(10, 0), 3 * 10 * 16);
		cout_capture capture;
		EXPECT_FALSE(i_tensor_file::matches(dir / "full.bin", tiny_nr, tiny_nmo, kept, single));
		EXPECT_NE(capture.buf.str().find("stores every pair, the previous layout; it is rebuilt"), std::string::npos) << capture.buf.str();
	}
	//No pairs kept
	{
		const int64_t h[5] = { magic_compact, tiny_nr, tiny_nmo, 0, 16 };
		write_raw(dir / "empty.bin", h, {}, {}, 0);
		EXPECT_FALSE(i_tensor_file::matches(dir / "empty.bin", tiny_nr, tiny_nmo, kept, single));
	}
	//An element size that is neither complex<float> nor complex<double>
	{
		const int64_t h[5] = { magic_compact, tiny_nr, tiny_nmo, 4, 12 };
		write_raw(dir / "elem.bin", h, tiny_mu, tiny_nu, 3 * 4 * 12);
		EXPECT_FALSE(i_tensor_file::matches(dir / "elem.bin", tiny_nr, tiny_nmo, kept, single));
	}
	//Header and pairs but one block short
	{
		const int64_t h[5] = { magic_compact, tiny_nr, tiny_nmo, 4, 8 };
		write_raw(dir / "short.bin", h, tiny_mu, tiny_nu, 2 * 4 * 8);
		EXPECT_FALSE(i_tensor_file::matches(dir / "short.bin", tiny_nr, tiny_nmo, kept, single));
	}
	//Shorter than the header itself
	{
		std::ofstream(dir / "stub.bin", std::ios::binary) << "NA2I";
		EXPECT_FALSE(i_tensor_file::matches(dir / "stub.bin", tiny_nr, tiny_nmo, kept, single));
	}
	//The good file still matches after all that, and reports its shape
	EXPECT_TRUE(i_tensor_file::matches(good, tiny_nr, tiny_nmo, kept, single));
	EXPECT_EQ(kept, 4u);
	EXPECT_FALSE(single);
	std::filesystem::remove_all(dir);
}

TEST(XcwScfTensorFileTests, OpenRejectsDamagedHeaders)
{
	const auto dir = scratch_dir();

	EXPECT_NE(open_error(dir / "absent.bin").find("cannot open"), std::string::npos);
	{
		const int64_t h[5] = { magic_full, tiny_nr, tiny_nmo, 4, 16 };
		write_raw(dir / "magic.bin", h, tiny_mu, tiny_nu, 3 * 4 * 16);
		EXPECT_NE(open_error(dir / "magic.bin").find("is not an I tensor in the compact layout"), std::string::npos);
	}
	{
		const int64_t h[5] = { magic_compact, tiny_nr, tiny_nmo, 4, 12 };
		write_raw(dir / "elem.bin", h, tiny_mu, tiny_nu, 3 * 4 * 12);
		EXPECT_NE(open_error(dir / "elem.bin").find("has an unknown element size"), std::string::npos);
	}
	{
		const int64_t h[5] = { magic_compact, 0, tiny_nmo, 4, 16 };
		write_raw(dir / "nr0.bin", h, tiny_mu, tiny_nu, 0);
		EXPECT_NE(open_error(dir / "nr0.bin").find("has a non-positive dimension"), std::string::npos);
	}
	{
		//11 pairs of a 4-function basis: more than the 10 that exist
		const int64_t h[5] = { magic_compact, tiny_nr, tiny_nmo, 11, 16 };
		write_raw(dir / "pairs.bin", h, ivec(11, 0), ivec(11, 0), 3 * 11 * 16);
		const std::string msg = open_error(dir / "pairs.bin");
		EXPECT_NE(msg.find("claims implausible dimensions"), std::string::npos) << msg;
		EXPECT_NE(msg.find("pairs 11"), std::string::npos) << msg;
	}
	{
		const int64_t h[5] = { magic_compact, tiny_nr, tiny_nmo, 4, 16 };
		write_raw(dir / "trunc.bin", h, { 0, 0 }, {}, 0);
		EXPECT_NE(open_error(dir / "trunc.bin").find("truncated pair list"), std::string::npos);
	}
	{
		//mu > nu at entry 1
		const int64_t h[5] = { magic_compact, tiny_nr, tiny_nmo, 4, 16 };
		write_raw(dir / "swap.bin", h, { 0, 2, 1, 3 }, { 0, 1, 3, 3 }, 3 * 4 * 16);
		EXPECT_NE(open_error(dir / "swap.bin").find("has a corrupt pair list at entry 1"), std::string::npos);
	}
	{
		//nu past the basis at entry 2
		const int64_t h[5] = { magic_compact, tiny_nr, tiny_nmo, 4, 16 };
		write_raw(dir / "range.bin", h, { 0, 0, 1, 3 }, { 0, 2, 4, 3 }, 3 * 4 * 16);
		EXPECT_NE(open_error(dir / "range.bin").find("has a corrupt pair list at entry 2"), std::string::npos);
	}
	{
		//Pairs out of order at entry 3
		const int64_t h[5] = { magic_compact, tiny_nr, tiny_nmo, 4, 16 };
		write_raw(dir / "order.bin", h, { 0, 0, 3, 1 }, { 0, 2, 3, 3 }, 3 * 4 * 16);
		EXPECT_NE(open_error(dir / "order.bin").find("has a corrupt pair list at entry 3"), std::string::npos);
	}
	{
		//Right header, one byte too many
		const int64_t h[5] = { magic_compact, tiny_nr, tiny_nmo, 4, 16 };
		write_raw(dir / "long.bin", h, tiny_mu, tiny_nu, 3 * 4 * 16 + 1);
		const std::string msg = open_error(dir / "long.bin");
		EXPECT_NE(msg.find("is 265 bytes, but its header implies 264"), std::string::npos) << msg;
	}
	//A sound one opens after all the rejections and the rejected handles are closed
	{
		const int64_t h[5] = { magic_compact, tiny_nr, tiny_nmo, 4, 16 };
		write_raw(dir / "sound.bin", h, tiny_mu, tiny_nu, 3 * 4 * 16);
		EXPECT_EQ(open_error(dir / "sound.bin"), "");
	}
	std::filesystem::remove_all(dir);
}

//------------------------------------------------------------------
//XCW runs on P1

//lambda = 0 is a plain Hartree-Fock: the energy of P1_test_XCW.good, the .tscb/.wfn/.fchk
//of that step, damping and level shift switched off along the way. `save` writes the
//tensor on a thread; a second run `read`s it instead of building and lands on the same
//numbers.
TEST(XcwScfTests, SaveThenReadTensorReproducesLambdaZero)
{
	if (p1_fixture().empty()) GTEST_SKIP() << "fixture tests/P1_test not found";
	const auto dir = scratch_dir();
	const std::string settings = common + "f rhf start 0 step_size 0.01 end 0 read tensor.bin save tensor.bin";

	const p1_run first = run_on_p1(dir, settings, true);
	EXPECT_EQ(first.nr, 402);
	EXPECT_NE(first.out.find("XCW orbital basis set: sto-3g-basis"), std::string::npos);
	EXPECT_EQ(first.out.find("Could not open anomalous dispersion file"), std::string::npos);
	EXPECT_EQ(first.out.find("I tensor read from"), std::string::npos);
	EXPECT_NE(first.out.find("Writing the I tensor to tensor.bin in the background"), std::string::npos) << first.out;
	EXPECT_NE(first.out.find("I tensor written to tensor.bin"), std::string::npos) << first.out;
	EXPECT_NE(first.out.find("More detailed output in XCW.log file..."), std::string::npos);
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

	//The saved tensor has this problem's shape and the single precision it was built in
	size_t kept = 0;
	bool single = false;
	ASSERT_TRUE(i_tensor_file::matches(dir / "tensor.bin", first.nr, p1_nmo, kept, single));
	EXPECT_EQ(kept, static_cast<size_t>(p1_pairs - screened));
	EXPECT_TRUE(single);
	EXPECT_EQ(std::filesystem::file_size(dir / "tensor.bin"), 40u + 2u * kept * sizeof(int) + i_tensor_file::total_bytes(first.nr, kept, true));

	const p1_run second = run_on_p1(dir, settings, true);
	EXPECT_NE(second.out.find("I tensor read from tensor.bin"), std::string::npos) << second.out;
	EXPECT_NE(second.out.find(", single precision), not recomputed, held in memory"), std::string::npos) << second.out;
	EXPECT_EQ(second.out.find("Screened out"), std::string::npos);
	EXPECT_EQ(second.out.find("Writing the I tensor"), std::string::npos);
	EXPECT_EQ(second.out.find("NOTE:"), std::string::npos);
	const auto rows2 = lambda_rows(second.out);
	ASSERT_EQ(rows2.size(), 1u) << second.out;
	EXPECT_EQ(rows2[0].criterion, rows[0].criterion);
	EXPECT_EQ(rows2[0].goof2, rows[0].goof2);
	EXPECT_NEAR(rows2[0].d(rows2[0].energy), rows[0].d(rows[0].energy), 1e-8);
	EXPECT_NE(second.log.find("***SCF converged in "), std::string::npos);

	std::filesystem::remove_all(dir);
}

//A one MB budget cannot hold the 11 MB tensor: it goes to I_tensor_stream.bin and the SCF
//walks it a window at a time. The perturbation at lambda = 0.01 then comes from the file,
//and the result is the in-memory one.
TEST(XcwScfTests, StreamedTensorMatchesInMemory)
{
	if (p1_fixture().empty()) GTEST_SKIP() << "fixture tests/P1_test not found";
	const auto dir = scratch_dir();
	const std::string settings = common + "f rhf start 0.01 step_size 0.01 end 0.01";

	const p1_run held = run_on_p1(dir, settings, true);
	EXPECT_EQ(held.out.find("I tensor streamed to disk"), std::string::npos);
	EXPECT_FALSE(std::filesystem::exists(dir / "I_tensor_stream.bin"));
	const auto rows = lambda_rows(held.out);
	ASSERT_EQ(rows.size(), 1u) << held.out;
	EXPECT_EQ(rows[0].lambda, "0.01000");
	EXPECT_TRUE(std::filesystem::exists(dir / "NA2_0010000.tscb"));
	//quant = energy + lambda * criterion; the printed penalty is that product to 3 decimals
	EXPECT_NEAR(rows[0].d(rows[0].quant), rows[0].d(rows[0].energy) + 0.01 * rows[0].d(rows[0].criterion), 1e-3);
	EXPECT_NEAR(rows[0].d(rows[0].penalty), 0.01 * rows[0].d(rows[0].criterion), 1e-3);
	//The perturbation raises the energy above the Hartree-Fock minimum, by little at 0.01
	EXPECT_GT(rows[0].d(rows[0].energy), golden_energy);
	EXPECT_LT(rows[0].d(rows[0].energy), golden_energy + 0.5);
	EXPECT_TRUE(physical_criterion(rows[0].criterion)) << rows[0].criterion;
	EXPECT_TRUE(physical_criterion(rows[0].goof2)) << rows[0].goof2;
	EXPECT_NE(held.log.find("Starting XCW SCF solver with lambda = 0.01000"), std::string::npos);
	EXPECT_NE(held.log.find("***SCF converged in "), std::string::npos);
	const int screened = static_cast<int>(number_after(held.out, "Screened out "));
	EXPECT_GT(screened, 0);
	const size_t kept = static_cast<size_t>(p1_pairs - screened);

	std::filesystem::remove(dir / "NA2_0010000.tscb");
	const p1_run streamed = run_on_p1(dir, settings + " i_tensor_mb 1", true);
	const std::string line_key = "I tensor streamed to disk: ";
	ASSERT_NE(streamed.out.find(line_key), std::string::npos) << streamed.out;
	EXPECT_NE(streamed.out.find("(single precision) total, "), std::string::npos) << streamed.out;
	EXPECT_NE(streamed.out.find("to fit i_tensor_mb (1.00 MB)"), std::string::npos) << streamed.out;
	//The window is what one MiB holds of kept complex<float> blocks
	const size_t window = static_cast<size_t>(number_after(streamed.out, "(single precision) total, "));
	EXPECT_EQ(window, (1024u * 1024u) / i_tensor_file::block_bytes(kept, true));
	EXPECT_NE(streamed.out.find(std::to_string(window) + " of " + std::to_string(held.nr) + " reflections resident"), std::string::npos) << streamed.out;
	EXPECT_NEAR(number_after(streamed.out, line_key), i_tensor_file::total_bytes(held.nr, kept, true) / 1048576.0, 0.01);

	ASSERT_TRUE(std::filesystem::exists(dir / "I_tensor_stream.bin"));
	size_t kept_on_disk = 0;
	bool single = false;
	EXPECT_TRUE(i_tensor_file::matches(dir / "I_tensor_stream.bin", held.nr, p1_nmo, kept_on_disk, single));
	EXPECT_EQ(kept_on_disk, kept);
	EXPECT_TRUE(single);

	const auto rows2 = lambda_rows(streamed.out);
	ASSERT_EQ(rows2.size(), 1u) << streamed.out;
	EXPECT_EQ(rows2[0].criterion, rows[0].criterion);
	EXPECT_EQ(rows2[0].goof2, rows[0].goof2);
	EXPECT_NEAR(rows2[0].d(rows2[0].energy), rows[0].d(rows[0].energy), 1e-8);
	EXPECT_NEAR(rows2[0].d(rows2[0].quant), rows[0].d(rows[0].quant), 1e-8);
	EXPECT_TRUE(std::filesystem::exists(dir / "NA2_0010000.tscb"));
	EXPECT_NE(streamed.log.find("***SCF converged in "), std::string::npos);

	std::filesystem::remove_all(dir);
}

//`guess_basis`: the first lambda starts from OCC's own Hartree-Fock in that basis. With the
//orbital basis itself as the guess basis the small SCF is the whole answer, so the log line
//carries the golden energy and the XCW loop has nothing left to do beyond recognising it.
TEST(XcwScfTests, SmallBasisGuessStartsFromConvergedDensity)
{
	if (p1_fixture().empty()) GTEST_SKIP() << "fixture tests/P1_test not found";
	const auto dir = scratch_dir();
	const p1_run run = run_on_p1(dir, common + "f rhf start 0 step_size 0.01 end 0 guess_basis sto-3g", true);

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
//itself. The unrestricted branches of the effective density, the level shift and the
//orbital gradient all run here, and the RHF run on the same subset is the reference.
TEST(XcwScfTests, UnrestrictedClosedShellAgainstF2MatchesRestricted)
{
	if (p1_fixture().empty()) GTEST_SKIP() << "fixture tests/P1_test not found";
	const auto dir = scratch_dir();
	const p1_run run = run_on_p1(dir, common + "f2 uhf start 0 step_size 0.01 end 0", true);

	ASSERT_NE(run.log.find("***SCF converged in "), std::string::npos) << run.log;
	EXPECT_EQ(occurrences(run.log, "***Turned off level shift***"), 1u);
	const auto rows = lambda_rows(run.out);
	ASSERT_EQ(rows.size(), 1u) << run.out;
	EXPECT_EQ(rows[0].criterion, rows[0].goof2);
	EXPECT_NEAR(rows[0].d(rows[0].energy), golden_energy, 1e-6);
	EXPECT_EQ(rows[0].quant, rows[0].energy);
	EXPECT_TRUE(std::filesystem::exists(dir / "NA2_0000000.tscb"));
	EXPECT_TRUE(std::filesystem::exists(dir / "NA2_0000000.fchk"));

	const p1_run rhf = run_on_p1(dir, common + "f2 rhf start 0 step_size 0.01 end 0", true);
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
	const p1_run run = run_on_p1(dir, "normal fast_conv params 177 basis_set sto-3g max_iter 100 charge 0 mult 1 f rhf start 0 step_size 0.01 end 0 weighted", true);

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

//Five iterations at lambda = 0.01 from the core guess, twice: the same five XCW.log rows to
//the last printed digit (the I tensor contraction, the Fock build and DIIS are all in
//them), energies within a few Hartree of the Hartree-Fock minimum and lower at the end
//than at the start, criteria positive. No converged result is asked for.
TEST(XcwScfTests, FiveIterationsAreDeterministicAndPhysical)
{
	if (p1_fixture().empty()) GTEST_SKIP() << "fixture tests/P1_test not found";
	const auto dir = scratch_dir();
	const std::string settings = "normal normal_conv params 177 basis_set sto-3g max_iter 5 charge 0 mult 1 f rhf start 0.01 step_size 0.01 end 0.01";

	const p1_run first = run_on_p1(dir, settings, true);
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

	const p1_run second = run_on_p1(dir, settings, true);
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
	const p1_run run = run_on_p1(dir, "normal normal_conv params 177 basis_set sto-3g max_iter 1 charge 0 mult 1 f rhf start 0 step_size 0.01 end 0.01", false);

	//No anomalous dispersion file this time: said once, then carried on
	EXPECT_EQ(occurrences(run.out, "Could not open anomalous dispersion file. Continuing without anomalous dispersions."), 1u) << run.out;
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
	EXPECT_NE(run.log.find("\t1\t\t"), std::string::npos);
	EXPECT_EQ(run.log.find("\t2\t\t"), std::string::npos);

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
//The I tensor against closed-form numbers: H2 in a cubic P1 cell, STO-3G, U = 0, so
//I_mn(h) = integral phi_m phi_n exp(i q.r) d3r with q = 2 pi/a (h k l) is a sum over the
//nine primitive pairs of the Gaussian product theorem, (pi/p)^(3/2) exp(-ab/p |A-B|^2)
//exp(-|q|^2/(4p)) exp(i q.P), p = a + b, P = (aA + bB)/p, with the contraction normalised
//as the basis loader does. What the grid quadrature, the Hirshfeld partition, the
//phase and DW factors and the pair packing produce has to agree with that to the grid's
//accuracy (5.6e-6 at accuracy 3 on FLOWOFFICE, 3.3e-4 at the default 2, both largest at
//0 0 4), I(-h) has to be the conjugate of I(h) and no diagonal element can exceed the
//overlap of a normalised function with itself.
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

TEST(XcwScfTests, ITensorOfH2MatchesAnalyticGaussianTransform)
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
	const int nr = static_cast<int>(hkl.size());

	const p1_run run = run_xcw(dir, cif, hkl_path, {}, "normal normal_conv params 1 basis_set sto-3g max_iter 100 charge 0 mult 1 f rhf start 0 step_size 0.01 end 0 save tensor.bin", true, 3);
	EXPECT_NE(run.out.find("I tensor written to tensor.bin"), std::string::npos) << run.out;
	ASSERT_NE(run.log.find("***SCF converged in "), std::string::npos) << run.log;

	i_tensor_file f;
	f.open(dir / "tensor.bin", nr);
	EXPECT_EQ(f.nr(), nr);
	EXPECT_EQ(f.nmo(), 2);
	ASSERT_EQ(f.kept(), 3u);
	EXPECT_FALSE(f.single());
	EXPECT_EQ(f.pair_mu(), (ivec{ 0, 0, 1 }));
	EXPECT_EQ(f.pair_nu(), (ivec{ 0, 1, 1 }));
	f.load(0, nr);

	const h2_analytic ref;
	const double a_bohr = constants::ang2bohr(a_ang);
	const vec A = { 0.5 * a_bohr, 0.5 * a_bohr, (0.5 - dz) * a_bohr };
	const vec B = { 0.5 * a_bohr, 0.5 * a_bohr, (0.5 + dz) * a_bohr };
	const vec centres[2] = { A, B };
	double max_dev = 0.0;
	int r = 0;
	std::vector<cvec> got(nr);
	for (const i3& h : hkl) {
		const vec q = { constants::TWO_PI / a_bohr * h[0], constants::TWO_PI / a_bohr * h[1], constants::TWO_PI / a_bohr * h[2] };
		const cdouble* block = f.block(r);
		got[r].assign(block, block + 3);
		for (size_t t = 0; t < 3; t++) {
			const cdouble expected = ref.I(centres[f.pair_mu()[t]], centres[f.pair_nu()[t]], q);
			const double dev = std::abs(block[t] - expected);
			max_dev = std::max(max_dev, dev);
			EXPECT_LT(dev, 5e-5) << "hkl " << h[0] << " " << h[1] << " " << h[2] << " pair " << t << ": got " << block[t] << " expected " << expected;
		}
		EXPECT_LE(std::abs(block[0]), 1.0 + 1e-6);
		EXPECT_LE(std::abs(block[2]), 1.0 + 1e-6);
		r++;
	}
	std::cerr << "H2 I tensor: max |grid - analytic| = " << max_dev << " over " << nr << " reflections" << std::endl;

	//Friedel mates are complex conjugates: the same grid points and weights enter, only
	//the phase flips, so this holds to rounding
	r = 0;
	for (const i3& h : hkl) {
		const auto mate = std::distance(hkl.begin(), hkl.find({ -h[0], -h[1], -h[2] }));
		for (size_t t = 0; t < 3; t++) {
			EXPECT_NEAR(got[r][t].real(), got[mate][t].real(), 1e-10);
			EXPECT_NEAR(got[r][t].imag(), -got[mate][t].imag(), 1e-10);
		}
		r++;
	}
	f.close();
	std::filesystem::remove_all(dir);
}

//`i_sigma <x>` keeps only I/sigma(I) >= x in the fit: the run says how many, the lambda
//table adds Crit(all) and R1(all) behind the seven fit-set columns, and the lambda = 0
//energy does not depend on which reflections are scored.
TEST(XcwScfTests, ISigmaCutoffShrinksFitSetAndReportsAllReflections)
{
	if (p1_fixture().empty()) GTEST_SKIP() << "fixture tests/P1_test not found";
	const auto dir = scratch_dir();
	const p1_run run = run_on_p1(dir, common + "f rhf start 0 step_size 0.01 end 0 i_sigma 20", true);

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
	const size_t row = run.out.find("\n0.00000\t");
	ASSERT_NE(row, std::string::npos) << run.out;
	std::vector<std::string> fields;
	std::istringstream cols(run.out.substr(row + 1, run.out.find('\n', row + 1) - row - 1));
	for (std::string f; std::getline(cols, f, '\t');) if (!f.empty()) fields.push_back(f);
	ASSERT_EQ(fields.size(), 9u) << run.out;
	EXPECT_NE(fields[7], fields[1]);
	EXPECT_NE(fields[8], fields[3]);
	std::filesystem::remove_all(dir);
}
