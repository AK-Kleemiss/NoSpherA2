//XCW SCF driver and I tensor storage: run_XCW_fitting on the P1 fixture through the
//settings file, observed through XCW.log, stdout and the files a run leaves behind, plus
//the i_tensor_file container on its own.
//
//Every run uses the sto-3g basis of tests/P1_test at a single lambda, the cheapest XCW
//there is, and reads its expectations from tests/P1_test/P1_test_XCW.good (the golden the
//integration test compares against): E(lambda = 0) = -1961.923538820 Eh, GooF(F2) 4.696,
//GooF(F) 4.666 with the anomalous dispersion file, the weighted criterion 9.733
//(P1_test_XCW_h2.good), 1798 of the 103*104/2 = 5356 pairs screened out.
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
#include <iostream>
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
    };

    //One XCW run in dir with the given settings file, no GPU and no timing lines. The
    //anomalous dispersion file of the fixture is passed when with_anom, so that the
    //criteria match P1_test_XCW.good; without it the run says so and continues.
    p1_run run_on_p1(const std::filesystem::path& dir, const std::string& settings_text, const bool with_anom)
    {
        const auto fixture = p1_fixture();
        const auto settings = dir / "settings.txt";
        std::ofstream(settings) << settings_text;
        options opt;
        opt.xcw_settings_path = settings;
        opt.cif = std::filesystem::absolute(fixture / "P1_test_NA2.cif");
        opt.hkl = std::filesystem::absolute(fixture / "P1_test.hkl");
        if (with_anom) {
            opt.anom_disp_path = std::filesystem::absolute(fixture / "anom_disp.txt");
        }
        opt.do_XCW = true;
        opt.use_gpu = false;
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

    //One line of the lambda table run_XCW_fitting prints per converged step:
    //lambda(5) criterion(3) GooF2(3) energy(9) lambda*criterion(3) quant(9), tab separated
    struct lambda_row {
        std::string lambda, criterion, goof2, energy, penalty, quant;
        double d(const std::string& s) const { return std::stod(s); }
    };

    std::vector<lambda_row> lambda_rows(const std::string& out)
    {
        std::vector<lambda_row> rows;
        std::istringstream in(out);
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || !std::isdigit(static_cast<unsigned char>(line[0])) || line.find('\t') == std::string::npos) {
                continue;
            }
            std::vector<std::string> fields;
            std::istringstream cols(line);
            std::string f;
            while (std::getline(cols, f, '\t')) {
                if (!f.empty()) fields.push_back(f);
            }
            if (fields.size() < 6) {
                continue;
            }
            rows.push_back({ fields[0], fields[1], fields[2], fields[3], fields[4], fields[5] });
        }
        return rows;
    }

    //P1_test_XCW.good, lambda = 0
    constexpr double golden_energy = -1961.923538820;
    constexpr double golden_goof2 = 4.696;
    constexpr double golden_goof1 = 4.666;
    constexpr int golden_screened = 1798;
    constexpr int p1_nmo = 103;
    constexpr int p1_nr = 3215;
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
        std::vector<cdouble> block(tiny_mu.size());
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
    EXPECT_NE(first.out.find("XCW orbital basis set: sto-3g-basis"), std::string::npos);
    EXPECT_EQ(first.out.find("Could not open anomalous dispersion file"), std::string::npos);
    EXPECT_EQ(first.out.find("I tensor read from"), std::string::npos);
    EXPECT_NE(first.out.find("Writing the I tensor to tensor.bin in the background"), std::string::npos) << first.out;
    EXPECT_NE(first.out.find("I tensor written to tensor.bin"), std::string::npos) << first.out;
    EXPECT_NE(first.out.find("More detailed output in XCW.log file..."), std::string::npos);
    EXPECT_NE(first.out.find("Finished XCW fitting procedure."), std::string::npos);
    EXPECT_EQ(static_cast<int>(number_after(first.out, "Screened out ")), golden_screened);

    const auto rows = lambda_rows(first.out);
    ASSERT_EQ(rows.size(), 1u) << first.out;
    EXPECT_EQ(rows[0].lambda, "0.00000");
    EXPECT_NEAR(rows[0].d(rows[0].energy), golden_energy, 1e-6);
    EXPECT_NEAR(rows[0].d(rows[0].criterion), golden_goof1, 5e-3);
    EXPECT_NEAR(rows[0].d(rows[0].goof2), golden_goof2, 5e-3);
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
        EXPECT_GE(tscb.reflection_size(), static_cast<size_t>(p1_nr));
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
    ASSERT_TRUE(i_tensor_file::matches(dir / "tensor.bin", p1_nr, p1_nmo, kept, single));
    EXPECT_EQ(kept, static_cast<size_t>(p1_nmo * (p1_nmo + 1) / 2 - golden_screened));
    EXPECT_TRUE(single);
    EXPECT_EQ(std::filesystem::file_size(dir / "tensor.bin"), 40u + 2u * kept * sizeof(int) + i_tensor_file::total_bytes(p1_nr, kept, true));

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

//A one MB budget cannot hold the 87 MB tensor: it goes to I_tensor_stream.bin and the SCF
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
    //The lambda = 0.01 row of P1_test_XCW.good, reached here from the core guess rather than
    //from the lambda = 0 wavefunction: the same minimum
    EXPECT_NEAR(rows[0].d(rows[0].energy), -1961.906493663, 1e-6);
    EXPECT_NEAR(rows[0].d(rows[0].criterion), 4.225, 5e-3);
    EXPECT_NEAR(rows[0].d(rows[0].goof2), 4.241, 5e-3);
    EXPECT_NE(held.log.find("Starting XCW SCF solver with lambda = 0.01000"), std::string::npos);
    EXPECT_NE(held.log.find("***SCF converged in "), std::string::npos);
    const int screened = static_cast<int>(number_after(held.out, "Screened out "));
    EXPECT_EQ(screened, golden_screened);
    const size_t kept = static_cast<size_t>(p1_nmo * (p1_nmo + 1) / 2 - screened);

    std::filesystem::remove(dir / "NA2_0010000.tscb");
    const p1_run streamed = run_on_p1(dir, settings + " i_tensor_mb 1", true);
    const std::string line_key = "I tensor streamed to disk: ";
    ASSERT_NE(streamed.out.find(line_key), std::string::npos) << streamed.out;
    EXPECT_NE(streamed.out.find("(single precision) total, "), std::string::npos) << streamed.out;
    EXPECT_NE(streamed.out.find("to fit i_tensor_mb (1.00 MB)"), std::string::npos) << streamed.out;
    //The window is what one MiB holds of kept complex<float> blocks
    const size_t window = static_cast<size_t>(number_after(streamed.out, "(single precision) total, "));
    EXPECT_EQ(window, (1024u * 1024u) / i_tensor_file::block_bytes(kept, true));
    EXPECT_NE(streamed.out.find(std::to_string(window) + " of " + std::to_string(p1_nr) + " reflections resident"), std::string::npos) << streamed.out;
    EXPECT_NEAR(number_after(streamed.out, line_key), i_tensor_file::total_bytes(p1_nr, kept, true) / 1048576.0, 0.01);

    ASSERT_TRUE(std::filesystem::exists(dir / "I_tensor_stream.bin"));
    size_t kept_on_disk = 0;
    bool single = false;
    EXPECT_TRUE(i_tensor_file::matches(dir / "I_tensor_stream.bin", p1_nr, p1_nmo, kept_on_disk, single));
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
    EXPECT_NEAR(rows[0].d(rows[0].goof2), golden_goof2, 5e-3);
    EXPECT_TRUE(std::filesystem::exists(dir / "NA2_0000000.tscb"));
    std::filesystem::remove_all(dir);
}

//A closed shell in UHF is the RHF answer; the criterion column under `f2` is GooF(F2)
//itself. The unrestricted branches of the effective density, the level shift and the
//orbital gradient all run here.
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
    EXPECT_NEAR(rows[0].d(rows[0].goof2), golden_goof2, 5e-3);
    EXPECT_NEAR(rows[0].d(rows[0].energy), golden_energy, 1e-6);
    EXPECT_EQ(rows[0].quant, rows[0].energy);
    EXPECT_TRUE(std::filesystem::exists(dir / "NA2_0000000.tscb"));
    EXPECT_TRUE(std::filesystem::exists(dir / "NA2_0000000.fchk"));
    std::filesystem::remove_all(dir);
}

//`fast_conv` runs DIIS alone: nothing to switch off, so the log never says so, and the
//weighted criterion is announced and used (9.733 at lambda = 0 in P1_test_XCW_h2.good).
//The preset sets neither damping nor shift; SCF_settings defaults both to 0.
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
    EXPECT_NEAR(rows[0].d(rows[0].criterion), 9.733, 5e-3);
    EXPECT_NEAR(rows[0].d(rows[0].goof2), golden_goof2, 5e-3);
    EXPECT_NEAR(rows[0].d(rows[0].energy), golden_energy, 1e-6);
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

    //step / 128 with step 0.01 = 7.8125e-5, whose double lies just above the tie, so eight
    //decimals round up
    const std::string stop = "XCW: unable to converge lambda 0.00000000 with a continuation step above 0.00007813; stopping scan.";
    EXPECT_NE(run.out.find(stop), std::string::npos) << run.out;
    EXPECT_NE(run.log.find(stop), std::string::npos);
    EXPECT_EQ(run.out.find("retrying lambda"), std::string::npos);
    EXPECT_NE(run.out.find("Finished XCW fitting procedure."), std::string::npos);
    EXPECT_TRUE(lambda_rows(run.out).empty()) << run.out;
    EXPECT_FALSE(std::filesystem::exists(dir / "NA2_0000000.tscb"));
    EXPECT_FALSE(std::filesystem::exists(dir / "NA2_0010000.tscb"));
    std::filesystem::remove_all(dir);
}
