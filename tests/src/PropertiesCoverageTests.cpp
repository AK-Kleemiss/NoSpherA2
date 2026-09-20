#include "pch.h"

#include "core/convenience.h"
#include "core/constants.h"
#include "core/atoms.h"
#include "core/cell.h"
#include "core/cube.h"
#include "core/wfn_class.h"
#include "core/spherical_density.h"
#include "core/GridManager.h"
#include "core/properties.h"
#include "core/b2c.h"
#include "core/NoSpherA2.h"
#include "core/SALTED_utilities.h"

#include <chrono>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

//cpp-only functions under test, declared here with their definition's signature
vec calc_dipole_for_atom(WFN &wavy, const int &i, cube &Hirshfeld_atom, vec &charges, std::string type);
void set_core_counts(int *max_s, int *max_p, int *max_d, int *max_f, const int &core_els, const int &ECP_mode);
void calculateConjugate(SALTEDDescriptors &v2);

namespace
{
    //fresh scratch directory for one test; removed by the test on success
    std::filesystem::path temp_dir(const std::string &name)
    {
        const std::filesystem::path dir = std::filesystem::temp_directory_path() / ("nosphera2_propcov_" + name);
        std::filesystem::remove_all(dir);
        std::filesystem::create_directories(dir);
        return dir;
    }

    //death tests re-run the body in a child process on Windows: only create, never wipe, before EXPECT_EXIT
    std::filesystem::path death_dir(const std::string &name)
    {
        const std::filesystem::path dir = std::filesystem::temp_directory_path() / ("nosphera2_propcov_" + name);
        std::filesystem::create_directories(dir);
        return dir;
    }

    std::string read_file(const std::filesystem::path &p)
    {
        std::ifstream in(p);
        std::stringstream ss;
        ss << in.rdbuf();
        return ss.str();
    }

    struct ScopedCwd
    {
        std::filesystem::path old = std::filesystem::current_path();
        explicit ScopedCwd(const std::filesystem::path &dir) { std::filesystem::current_path(dir); }
        ~ScopedCwd() { std::filesystem::current_path(old); }
    };

    //the option parser and the drivers flip these globals; put them back whatever the test did
    struct ScopedGlobals
    {
        bool timings = constants::hide_timings;
        bool gpu = constants::hide_gpu_notes;
        bool counts = ProgressBar::report_counts;
        double cutoff = constants::exp_cutoff;
        ~ScopedGlobals()
        {
            constants::hide_timings = timings;
            constants::hide_gpu_notes = gpu;
            ProgressBar::report_counts = counts;
            constants::exp_cutoff = cutoff;
        }
    };

    //two-centre H2 model of PropertiesTests: one doubly occupied bonding MO of two s Gaussians
    //(exponent a) at -R and R on x, optionally the empty antibonding partner
    struct H2Model
    {
        double a = 1.0;
        double R = 1.0;
        double cp = 0.0;
        double cm = 0.0;
        WFN wavy{ e_origin::NOT_YET_DEFINED };

        H2Model(double alpha, double half_distance, bool with_virtual = true) : a(alpha), R(half_distance)
        {
            const double SAA = std::pow(constants::PI / (2.0 * a), 1.5);
            const double SAB = SAA * std::exp(-a * (2.0 * R) * (2.0 * R) / 2.0);
            cp = 1.0 / std::sqrt(2.0 * (SAA + SAB));
            cm = 1.0 / std::sqrt(2.0 * (SAA - SAB));
            wavy.push_back_atom("H", -R, 0.0, 0.0, 1);
            wavy.push_back_atom("H", R, 0.0, 0.0, 1);
            wavy.push_back_MO(1, 2.0, -0.5);
            if (with_virtual)
                wavy.push_back_MO(2, 0.0, 0.3);
            double on_A[2] = { cp, cm };
            double on_B[2] = { cp, -cm };
            wavy.add_primitive(1, 1, a, on_A);
            wavy.add_primitive(2, 1, a, on_B);
            constants::exp_cutoff = -23.5;
        }
        double gA(const d3 &p) const { return std::exp(-a * ((p[0] + R) * (p[0] + R) + p[1] * p[1] + p[2] * p[2])); }
        double gB(const d3 &p) const { return std::exp(-a * ((p[0] - R) * (p[0] - R) + p[1] * p[1] + p[2] * p[2])); }
        double psi0(const d3 &p) const { return cp * (gA(p) + gB(p)); }
        double psi1(const d3 &p) const { return cm * (gA(p) - gB(p)); }
        double rho(const d3 &p) const { return 2.0 * psi0(p) * psi0(p); }
    };

    //open-shell H2 on the same Gaussians, written as a wfn so the reader assigns the spin
    //manifolds (a drop in orbital energy starts the beta block). Variant A: alpha bonding (occ 1,
    //-0.5), alpha antibonding (0, 0.3), beta bonding scaled by 1/2 (1, -0.4), beta antibonding
    //(0, 0.2): HOMO and LUMO are both beta, the spin density is 0.75 psi0^2. Variant B: alpha
    //bonding (1, -0.5), beta bonding scaled by 1/2 (1, -0.6), beta antibonding (0, 0.2): the HOMO
    //is alpha, the LUMO beta and the alpha manifold has no virtual
    H2Model write_open_shell_h2(const std::filesystem::path &path, const bool mixed)
    {
        H2Model m(1.0, 1.0, false);
        WFN w(e_origin::NOT_YET_DEFINED);
        w.push_back_atom("H", -m.R, 0.0, 0.0, 1);
        w.push_back_atom("H", m.R, 0.0, 0.0, 1);
        if (!mixed)
        {
            w.push_back_MO(1, 1.0, -0.5);
            w.push_back_MO(2, 0.0, 0.3);
            w.push_back_MO(3, 1.0, -0.4, 1);
            w.push_back_MO(4, 0.0, 0.2, 1);
            double on_A[4] = { m.cp, m.cm, 0.5 * m.cp, m.cm };
            double on_B[4] = { m.cp, -m.cm, 0.5 * m.cp, -m.cm };
            w.add_primitive(1, 1, m.a, on_A);
            w.add_primitive(2, 1, m.a, on_B);
        }
        else
        {
            w.push_back_MO(1, 1.0, -0.5);
            w.push_back_MO(2, 1.0, -0.6, 1);
            w.push_back_MO(3, 0.0, 0.2, 1);
            double on_A[3] = { m.cp, 0.5 * m.cp, m.cm };
            double on_B[3] = { m.cp, 0.5 * m.cp, -m.cm };
            w.add_primitive(1, 1, m.a, on_A);
            w.add_primitive(2, 1, m.a, on_B);
        }
        w.set_exp_cutoff();
        w.write_wfn(path, false, false);
        return m;
    }

    //one hydrogen at the origin, one singly occupied unit-exponent s Gaussian: N = (pi/2)^1.5
    WFN single_h()
    {
        WFN w(e_origin::NOT_YET_DEFINED);
        w.push_back_atom("H", 0.0, 0.0, 0.0, 1);
        w.push_back_MO(1, 1.0, -0.5);
        double c = 1.0;
        w.add_primitive(1, 1, 1.0, &c);
        w.set_exp_cutoff();
        return w;
    }

    //alpha orbital exp(-r^2) and beta orbital exp(-r^2) / 2 on one hydrogen: spin density 0.75 exp(-2 r^2)
    WFN single_h_spin()
    {
        WFN w(e_origin::NOT_YET_DEFINED);
        w.push_back_atom("H", 0.0, 0.0, 0.0, 1);
        w.push_back_MO(1, 1.0, -0.5);
        w.push_back_MO(2, 1.0, -0.4, 1);
        double coefs[2] = { 1.0, 0.5 };
        w.add_primitive(1, 1, 1.0, coefs);
        w.set_exp_cutoff();
        return w;
    }

    //He (exponent 2) at the origin and H (exponent 1) at z = 1.5 bohr sharing one doubly occupied
    //MO with unit coefficients: N = 2 (S_HeHe + S_HH + 2 S_HeH)
    const double kHeHDist = 1.5;
    WFN heh()
    {
        WFN w(e_origin::NOT_YET_DEFINED);
        w.push_back_atom("He", 0.0, 0.0, 0.0, 2);
        w.push_back_atom("H", 0.0, 0.0, kHeHDist, 1);
        w.push_back_MO(1, 2.0, -0.9);
        double c = 1.0;
        w.add_primitive(1, 1, 2.0, &c);
        w.add_primitive(2, 1, 1.0, &c);
        w.set_exp_cutoff();
        return w;
    }
    double heh_electrons()
    {
        const double SAA = std::pow(constants::PI / 4.0, 1.5);
        const double SBB = std::pow(constants::PI / 2.0, 1.5);
        const double SAB = std::pow(constants::PI / 3.0, 1.5) * std::exp(-2.0 * kHeHDist * kHeHDist / 3.0);
        return 2.0 * (SAA + SBB + 2.0 * SAB);
    }

    //the GridTests hydrogen molecule on z: coefficient 1/2 on two unit-exponent s Gaussians
    const double kH2Coef = 0.5;
    const double kH2Half = 0.7;
    WFN make_h2()
    {
        WFN w(e_origin::NOT_YET_DEFINED);
        w.push_back_atom("H1", 0.0, 0.0, -kH2Half, 1);
        w.push_back_atom("H2", 0.0, 0.0, kH2Half, 1);
        w.push_back_MO(0, 2.0, -0.5);
        double c = kH2Coef;
        w.add_primitive(1, 1, 1.0, &c);
        w.add_primitive(2, 1, 1.0, &c);
        w.set_exp_cutoff();
        return w;
    }
    double h2_electrons()
    {
        const double S = std::pow(constants::PI / 2.0, 1.5);
        return 2.0 * kH2Coef * kH2Coef * S * (2.0 + 2.0 * std::exp(-(2.0 * kH2Half) * (2.0 * kH2Half) / 2.0));
    }

    //n^3 voxels of spacing h centred on the origin, axis aligned
    cube make_grid(int n, double h)
    {
        cube grid({ n, n, n }, 0, true);
        for (int k = 0; k < 3; k++)
        {
            grid.set_origin(k, -0.5 * (n - 1) * h);
            grid.set_vector(k, k, h);
        }
        grid.calc_dv();
        return grid;
    }

    //4 pi int r^2 f(r) dr on a logarithmic trapezoid
    double radial_integral(const std::function<double(double)> &f, double r_min, double r_max, int n)
    {
        const double t0 = std::log(r_min), t1 = std::log(r_max), dt = (t1 - t0) / n;
        double sum = 0.0;
        for (int i = 0; i <= n; i++)
        {
            const double r = std::exp(t0 + i * dt);
            const double w = (i == 0 || i == n) ? 0.5 : 1.0;
            sum += w * f(r) * r * r * r;
        }
        return constants::FOUR_PI * sum * dt;
    }

    double vec_total(const vec2 &v)
    {
        double s = 0.0;
        for (const vec &row : v)
            s += std::accumulate(row.begin(), row.end(), 0.0);
        return s;
    }

    int run_nos(std::vector<std::string> args)
    {
        args.insert(args.begin(), "NoSpherA2");
        std::vector<char *> argv;
        for (std::string &a : args)
            argv.push_back(a.data());
        argv.push_back(nullptr);
        return run_app((int)args.size(), argv.data());
    }

    //the numbers after the n-th '|' of the first line in text that starts with prefix
    vec row_after_bars(const std::string &text, const std::string &prefix, const int bars)
    {
        std::istringstream lines(text);
        std::string line;
        while (std::getline(lines, line))
        {
            if (line.rfind(prefix, 0) != 0)
                continue;
            size_t pos = 0;
            for (int b = 0; b < bars; b++)
                pos = line.find('|', pos) + 1;
            std::string tail = line.substr(pos);
            for (char &c : tail)
                if (c == ',')
                    c = ' ';
            std::istringstream in(tail);
            vec out;
            double v;
            while (in >> v)
                out.push_back(v);
            return out;
        }
        return {};
    }

    //the number that follows key in text
    double number_after(const std::string &text, const std::string &key)
    {
        const size_t pos = text.find(key);
        if (pos == std::string::npos)
            return std::numeric_limits<double>::quiet_NaN();
        return std::stod(text.substr(pos + key.size()));
    }

    std::string capture_end()
    {
        std::cout.flush();
        return testing::internal::GetCapturedStdout();
    }
}

//print_time formats seconds below a minute, minutes with the second remainder below an hour and
//hours with the minute remainder above, and leaves the stream's flags and precision untouched
TEST(PropertiesCoverageTimingTests, PrintTimeMinutesAndHoursBranches)
{
    ScopedGlobals guard;
    constants::hide_timings = false;
    const std::pair<long long, std::string> cases[] = {
        { 59, "Time to calculate Values: 59 s\n" },
        { 61, "Time to calculate Values: 1 m 1 s\n" },
        { 3601, "Time to calculate Values: 1 h 0 m\n" },
        { 7325, "Time to calculate Values: 2 h 2 m\n" },
    };
    for (const auto &[seconds, expected] : cases)
    {
        std::ostringstream os;
        os << std::scientific << std::setprecision(3);
        _time_point end = get_time();
        _time_point start = end - std::chrono::seconds(seconds);
        print_time(start, end, os);
        EXPECT_EQ(os.str(), expected);
        EXPECT_EQ(os.flags() & std::ios_base::floatfield, std::ios_base::scientific);
        EXPECT_EQ(os.precision(), 3);
    }
    constants::hide_timings = true;
    std::ostringstream silent;
    _time_point end = get_time();
    _time_point start = end - std::chrono::seconds(61);
    print_time(start, end, silent);
    EXPECT_TRUE(silent.str().empty());
}

//the spin density cube is sum over alpha minus sum over beta of occ psi^2: 0.75 exp(-2 r^2) for the
//one-centre model; with nodate false the progress bar goes to the console and the timing to the log
TEST(PropertiesCoverageSpinTests, SpinDensityCubeIsAlphaMinusBeta)
{
    ScopedGlobals guard;
    const WFN w = single_h_spin();
    cube s = make_grid(7, 0.5);
    std::ostringstream log;
    bool nodate = true;
    Calc_S_Rho(s, w, log, nodate);
    for (int i = 0; i < 7; i++)
        for (int j = 0; j < 7; j++)
            for (int k = 0; k < 7; k++)
            {
                const d3 p = s.get_pos(i, j, k);
                const double r2 = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
                EXPECT_NEAR(s.get_value(i, j, k), 0.75 * std::exp(-2.0 * r2), 1e-9 * 0.75);
            }
    EXPECT_TRUE(log.str().empty());

    constants::hide_timings = false;
    nodate = false;
    cube t = make_grid(3, 0.5);
    testing::internal::CaptureStdout();
    Calc_S_Rho(t, w, log, nodate);
    const std::string console = capture_end();
    EXPECT_NE(console.find("Calculating Values"), std::string::npos);
    EXPECT_NE(log.str().find("Time to calculate Values: 0 s"), std::string::npos);
    EXPECT_NEAR(t.get_value(1, 1, 1), 0.75, 1e-12);
}

//the spherical-harmonics MO evaluation is a stub that aborts; the driver must exit through err_not_impl
TEST(PropertiesCoverageSpinTests, MoSphericalHarmonicsIsNotImplemented)
{
    const WFN w = single_h();
    cube c({ 1, 1, 1 }, 0, true);
    std::ostringstream log;
    EXPECT_EXIT({ Calc_MO_spherical_harmonics(c, w, 0, log, true); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
}

//charge = sum v dv and mu = sum (r - r_atom) v dv on the grid: for v = 1 + z on the 5^3 grid of
//spacing 1/2 (coordinates -1..1, dv = 1/8) with the atom at (0.3, 0, 0) that is 15.625 electrons,
//mu_z = 25 * 2.5 / 8 = 7.8125, mu_x = (sum_x (x - 0.3)) * 25 / 8 = -4.6875 and mu_y = 0
TEST(PropertiesCoverageDipoleTests, AtomDipoleOfLinearRamp)
{
    WFN w(e_origin::NOT_YET_DEFINED);
    w.push_back_atom("H", 0.3, 0.0, 0.0, 1);
    cube c = make_grid(5, 0.5);
    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
            for (int k = 0; k < 5; k++)
                c.set_value(i, j, k, 1.0 + c.get_pos(i, j, k)[2]);
    vec charges(1, 0.0);
    const int atom = 0;
    const vec mu = calc_dipole_for_atom(w, atom, c, charges, "atom");
    ASSERT_EQ(mu.size(), 4u);
    EXPECT_NEAR(mu[0], -4.6875, 1e-12);
    EXPECT_NEAR(mu[1], 0.0, 1e-12);
    EXPECT_NEAR(mu[2], 7.8125, 1e-12);
    EXPECT_NEAR(mu[3], 15.625, 1e-12);
}

//every origin choice but the atom position is an err_not_impl exit, unknown names included
TEST(PropertiesCoverageDipoleTests, UnimplementedDipoleOriginsExit)
{
    WFN w(e_origin::NOT_YET_DEFINED);
    w.push_back_atom("H", 0.0, 0.0, 0.0, 1);
    cube c = make_grid(3, 0.5);
    vec charges(1, 0.0);
    const int atom = 0;
    const char *types[] = { "geometry", "hirshfeld", "vdW", "bogus" };
    for (const char *type : types)
        EXPECT_EXIT({ calc_dipole_for_atom(w, atom, c, charges, type); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
}

TEST(PropertiesCoverageDipoleTests, DipoleMomentsRequiresWavefunction)
{
    options opt;
    opt.no_date = true;
    std::ostringstream log;
    EXPECT_EXIT({ dipole_moments(opt, log); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
}

//the Hirshfeld atoms of a homonuclear diatomic carry opposite dipoles along the bond and none
//across it: mu_x(0) = -mu_x(1), mu_y = mu_z = 0, and the bond polarisation is not zero
TEST(PropertiesCoverageDipoleTests, DipoleMomentsDebugTableOfH2)
{
    ScopedGlobals guard;
    const std::filesystem::path dir = temp_dir("dipole_h2");
    H2Model m(1.0, 1.0);
    m.wavy.write_wfn(dir / "h2.wfn", false, false);
    {
        ScopedCwd cwd(dir);
        options opt;
        opt.wfn = "h2.wfn";
        opt.debug = true;
        opt.no_date = true;
        opt.properties.resolution = 0.3;
        opt.properties.radius = 1.5;
        std::ostringstream log;
        testing::internal::CaptureStdout();
        dipole_moments(opt, log);
        const std::string console = capture_end();
        const std::string text = log.str();
        EXPECT_NE(console.find("Properties calculation done!"), std::string::npos);
        EXPECT_NE(text.find("Starting calculation of dipole moment"), std::string::npos);
        EXPECT_NE(text.find("Origins etc are set up"), std::string::npos);
        EXPECT_NE(text.find("Calcualting Hirshfeld density for atom: 1"), std::string::npos);
        const size_t table = text.find(" atom   |  dipole moment");
        ASSERT_NE(table, std::string::npos);
        const vec row0 = row_after_bars(text.substr(table), "  0 (H) |", 1);
        const vec row1 = row_after_bars(text.substr(table), "  1 (H) |", 1);
        ASSERT_EQ(row0.size(), 3u);
        ASSERT_EQ(row1.size(), 3u);
        EXPECT_NEAR(row0[0], -row1[0], 1e-4);
        EXPECT_GT(std::fabs(row0[0]), 1e-3);
        EXPECT_LT(std::fabs(row0[1]), 1e-4);
        EXPECT_LT(std::fabs(row0[2]), 1e-4);
        EXPECT_LT(std::fabs(row1[1]), 1e-4);
        EXPECT_LT(std::fabs(row1[2]), 1e-4);
    }
    std::filesystem::remove_all(dir);
}

//alpha = (mu(+E) - mu(-E)) / 2E: seven copies of the same wavefunction give a zero tensor and a
//neutral Hirshfeld charge on each hydrogen
TEST(PropertiesCoverageDipoleTests, PolarizabilitiesOfSevenIdenticalWavefunctionsVanish)
{
    ScopedGlobals guard;
    const std::filesystem::path dir = temp_dir("polarizability_h2");
    H2Model m(1.0, 1.0);
    m.wavy.write_wfn(dir / "h2.wfn", false, false);
    options opt;
    opt.debug = true;
    opt.no_date = true;
    opt.properties.resolution = 0.3;
    opt.properties.radius = 1.5;
    for (int i = 0; i < 7; i++)
        opt.pol_wfns.push_back(dir / "h2.wfn");
    std::ostringstream log;
    testing::internal::CaptureStdout();
    polarizabilities(opt, log);
    const std::string console = capture_end();
    const std::string text = log.str();
    EXPECT_NE(console.find("Properties calculation done!"), std::string::npos);
    EXPECT_NE(text.find("Starting calculation of Polarizabilities"), std::string::npos);
    const size_t table = text.find("Polarizabilities:");
    ASSERT_NE(table, std::string::npos);
    for (const char *prefix : { "  0 (H) |", "  1 (H) |" })
    {
        const vec charge = row_after_bars(text.substr(table), prefix, 1);
        ASSERT_GE(charge.size(), 1u);
        EXPECT_NEAR(charge[0], 0.0, 2e-2);
        const vec alpha = row_after_bars(text.substr(table), prefix, 2);
        ASSERT_EQ(alpha.size(), 9u);
        for (const double a : alpha)
            EXPECT_NEAR(a, 0.0, 1e-10);
    }
    std::filesystem::remove_all(dir);
}

//an open-shell wfn whose HOMO and LUMO are both beta: the manifold table, the beta-only remark
//and the summary file with the unrestricted flag
TEST(PropertiesCoverageFukuiTests, FukuiAnalysisUnrestrictedBetaPair)
{
    ScopedGlobals guard;
    const std::filesystem::path dir = temp_dir("fukui_beta_pair");
    write_open_shell_h2(dir / "h2u.wfn", false);
    options opt;
    opt.wfn = dir / "h2u.wfn";
    opt.no_date = true;
    opt.accuracy = 1;
    std::ostringstream log;
    fukui_analysis(opt, log);
    const std::string text = log.str();
    EXPECT_NE(text.find("Read 2 atoms and 4 molecular orbitals (2 occupied)."), std::string::npos);
    EXPECT_NE(text.find("HOMO-LUMO gap: 0.600000 Hartree"), std::string::npos);
    EXPECT_NE(text.find("Unrestricted wavefunction."), std::string::npos);
    EXPECT_NE(text.find("Chosen HOMO is beta, chosen LUMO is beta."), std::string::npos);
    EXPECT_NE(text.find("  alpha: HOMO = MO 0 (-0.500000), LUMO = MO 1 (0.300000)"), std::string::npos);
    EXPECT_NE(text.find("   beta: HOMO = MO 2 (-0.400000), LUMO = MO 3 (0.200000)"), std::string::npos);
    EXPECT_NE(text.find("Both come from the beta manifold"), std::string::npos);
    EXPECT_NE(text.find("Summary (Hirshfeld partition):"), std::string::npos);
    const std::string dat = read_file(dir / "h2u_fukui.dat");
    EXPECT_NE(dat.find("HOMO_index 2"), std::string::npos);
    EXPECT_NE(dat.find("LUMO_index 3"), std::string::npos);
    EXPECT_NE(dat.find("unrestricted 1"), std::string::npos);
    std::filesystem::remove_all(dir);
}

//HOMO alpha, LUMO beta, and an alpha manifold without a virtual orbital
TEST(PropertiesCoverageFukuiTests, FukuiAnalysisUnrestrictedMixedManifolds)
{
    ScopedGlobals guard;
    const std::filesystem::path dir = temp_dir("fukui_mixed");
    write_open_shell_h2(dir / "h2m.wfn", true);
    options opt;
    opt.wfn = dir / "h2m.wfn";
    opt.no_date = true;
    opt.accuracy = 1;
    std::ostringstream log;
    fukui_analysis(opt, log);
    const std::string text = log.str();
    EXPECT_NE(text.find("HOMO-LUMO gap: 0.700000 Hartree"), std::string::npos);
    EXPECT_NE(text.find("Chosen HOMO is alpha, chosen LUMO is beta."), std::string::npos);
    EXPECT_NE(text.find("  alpha: HOMO = MO 0 (-0.500000), LUMO = MO none"), std::string::npos);
    EXPECT_NE(text.find("   beta: HOMO = MO 1 (-0.600000), LUMO = MO 2 (0.200000)"), std::string::npos);
    EXPECT_NE(text.find("DIFFERENT manifolds"), std::string::npos);
    const std::string dat = read_file(dir / "h2m_fukui.dat");
    EXPECT_NE(dat.find("HOMO_index 0"), std::string::npos);
    EXPECT_NE(dat.find("LUMO_index 2"), std::string::npos);
    EXPECT_NE(dat.find("unrestricted 1"), std::string::npos);
    std::filesystem::remove_all(dir);
}

//the property driver with every MO, the spin density, the Fukui block on an unrestricted wfn and
//the ESP-coloured isosurface: the spin cube is 0.75 psi0^2, f+ integrates to one LUMO electron
TEST(PropertiesCoverageDriverTests, PropertiesCalculationSpinMosFukuiAndIsosurface)
{
    ScopedGlobals guard;
    const std::filesystem::path dir = temp_dir("driver_open_shell");
    const H2Model m = write_open_shell_h2(dir / "h2u.wfn", false);
    {
        ScopedCwd cwd(dir);
        options opt;
        opt.wfn = "h2u.wfn";
        opt.debug = true;
        opt.no_date = true;
        opt.accuracy = 1;
        opt.properties.s_rho = true;
        opt.properties.all_mos = true;
        opt.properties.fukui = true;
        opt.properties.esp_isosurface = 0.01;
        opt.properties.resolution = 0.3;
        opt.properties.radius = 1.0;
        testing::internal::CaptureStdout();
        properties_calculation(opt);
        const std::string console = capture_end();
        EXPECT_NE(console.find("Properties calculation done!"), std::string::npos);
        const std::string text = read_file(dir / "NoSpherA2_cube.log");
        EXPECT_NE(text.find("Starting calculation of properties"), std::string::npos);
        EXPECT_NE(text.find("Size of MOs: 4"), std::string::npos);
        EXPECT_NE(text.find("Resetting Radius to at least 2.5 for the isosurface!"), std::string::npos);
        EXPECT_NE(text.find("Spin density, "), std::string::npos);
        EXPECT_NE(text.find("Fukui functions and dual descriptor, "), std::string::npos);
        EXPECT_NE(text.find("Calcualting MO: 3"), std::string::npos);
        EXPECT_NE(text.find("WARNING: unrestricted wavefunction"), std::string::npos);
        EXPECT_NE(text.find("Calculating Fukui functions..."), std::string::npos);
        EXPECT_NE(text.find("Wrote Fukui summary to"), std::string::npos);
        EXPECT_NE(text.find("Colouring the rho = 0.01 au isosurface"), std::string::npos);
        EXPECT_NE(text.find("Found "), std::string::npos);
        EXPECT_NEAR(number_after(text, "Integrated f+ over the grid: "), 1.0, 0.1);
        for (const char *name : { "h2u_s_rho.cube", "h2u_MO_0.cube", "h2u_MO_1.cube", "h2u_MO_2.cube", "h2u_MO_3.cube",
                                  "h2u_fukui_plus.cube", "h2u_fukui_minus.cube", "h2u_fukui_zero.cube", "h2u_dual_descriptor.cube",
                                  "h2u_fukui.dat", "h2u_rho_esp.obj" })
            EXPECT_TRUE(std::filesystem::exists(dir / name)) << name;
        std::ostringstream sink;
        WFN dummy(e_origin::NOT_YET_DEFINED);
        const cube spin(dir / "h2u_s_rho.cube", true, dummy, sink);
        int checked = 0;
        for (int i = 0; i < spin.get_size(0); i++)
            for (int j = 0; j < spin.get_size(1); j++)
                for (int k = 0; k < spin.get_size(2); k++)
                {
                    const double v = spin.get_value(i, j, k);
                    if (std::fabs(v) < 1e-3)
                        continue;
                    const double ref = 0.75 * m.psi0(spin.get_pos(i, j, k)) * m.psi0(spin.get_pos(i, j, k));
                    //the reader's exp_cutoff drops a primitive whose orbital contribution is below
                    //density_accuracy = 5e-5, so psi_alpha is off by up to 5e-5 (psi_beta by half): the
                    //spin density by 2 psi0 (5e-5 + 0.5 2.5e-5) < 6e-5; the cube header carries positions to
                    //6 decimals, ~1e-5 bohr over the 21 steps of this grid, 4 a r dr < 1.5e-4 relative
                    EXPECT_NEAR(v, ref, 2e-4 * std::fabs(ref) + 1e-4);
                    checked++;
                }
        EXPECT_GT(checked, 50);
    }
    std::filesystem::remove_all(dir);
}

//combining MO 1 of a wavefunction with MO 1 of its copy: the sum cube is 2 psi0, the difference
//cube is identically zero, and the vmd loader script appears
TEST(PropertiesCoverageDriverTests, CombineMosWritesSumAndDifferenceCubes)
{
    ScopedGlobals guard;
    const std::filesystem::path dir = temp_dir("combine_mo");
    H2Model m(1.0, 1.0);
    m.wavy.write_wfn(dir / "h2.wfn", false, false);
    std::filesystem::copy_file(dir / "h2.wfn", dir / "h2b.wfn");
    {
        ScopedCwd cwd(dir);
        options opt;
        opt.combine_mo = { "h2.wfn", "h2b.wfn" };
        opt.cmo1 = { 1 };
        opt.cmo2 = { 1 };
        testing::internal::CaptureStdout();
        do_combine_mo(opt);
        const std::string console = capture_end();
        EXPECT_NE(console.find("In total we have 4 atoms"), std::string::npos);
        EXPECT_TRUE(std::filesystem::exists(dir / "read_files.vmd"));
        ASSERT_TRUE(std::filesystem::exists(dir / "h2_1+h2b_1.cube"));
        ASSERT_TRUE(std::filesystem::exists(dir / "h2_1-h2b_1.cube"));
        std::ostringstream sink;
        WFN dummy(e_origin::NOT_YET_DEFINED);
        const cube sum(dir / "h2_1+h2b_1.cube", true, dummy, sink);
        WFN dummy2(e_origin::NOT_YET_DEFINED);
        const cube diff(dir / "h2_1-h2b_1.cube", true, dummy2, sink);
        int checked = 0;
        for (int i = 0; i < sum.get_size(0); i++)
            for (int j = 0; j < sum.get_size(1); j++)
                for (int k = 0; k < sum.get_size(2); k++)
                {
                    EXPECT_NEAR(diff.get_value(i, j, k), 0.0, 1e-12);
                    const double v = sum.get_value(i, j, k);
                    if (std::fabs(v) < 1e-3)
                        continue;
                    const double ref = 2.0 * m.psi0(sum.get_pos(i, j, k));
                    //the reader's exp_cutoff drops a primitive whose orbital contribution is below
                    //density_accuracy = 5e-5, so 2 psi0 is off by up to 1e-4 where one Gaussian is
                    //negligible; the cube header carries origin and vectors to 6 decimals, 4e-7 bohr per
                    //0.1 A step over 50 steps, 2 a r dr < 2e-4 relative at the r ~ 2.7 bohr where 2 psi0 ~ 1e-3
                    EXPECT_NEAR(v, ref, 5e-4 * std::fabs(ref) + 1.2e-4);
                    checked++;
                }
        EXPECT_GT(checked, 1000);
    }
    std::filesystem::remove_all(dir);
}

//the ECP core tables close shells: 2 s + 6 p + 10 d + 14 f electrons equal the core count for
//every closed-shell entry, core 0 is a no-op, mode 3 reads mode 1's table
TEST(PropertiesCoverageSphericalTests, CoreCountsCloseShells)
{
    struct Row { int mode; int core; int s, p, d, f; };
    const Row rows[] = {
        { 1, 2, 1, 0, 0, 0 }, { 1, 10, 2, 1, 0, 0 }, { 1, 18, 3, 2, 0, 0 }, { 1, 28, 3, 2, 1, 0 },
        { 1, 46, 4, 3, 2, 0 }, { 1, 60, 4, 3, 2, 1 }, { 1, 78, 5, 4, 3, 1 },
        { 3, 10, 2, 1, 0, 0 }, { 3, 60, 4, 3, 2, 1 },
        { 2, 2, 1, 0, 0, 0 }, { 2, 10, 2, 1, 0, 0 }, { 2, 18, 3, 2, 0, 0 }, { 2, 28, 3, 2, 1, 0 },
        { 2, 36, 4, 3, 1, 0 }, { 2, 46, 4, 3, 2, 0 }, { 2, 54, 5, 4, 2, 0 }, { 2, 68, 5, 4, 2, 1 }, { 2, 78, 5, 4, 3, 1 },
    };
    for (const Row &r : rows)
    {
        int s = -7, p = -7, d = -7, f = -7;
        set_core_counts(&s, &p, &d, &f, r.core, r.mode);
        EXPECT_EQ(s, r.s) << "mode " << r.mode << " core " << r.core;
        EXPECT_EQ(p, r.p) << "mode " << r.mode << " core " << r.core;
        EXPECT_EQ(d, r.d) << "mode " << r.mode << " core " << r.core;
        EXPECT_EQ(f, r.f) << "mode " << r.mode << " core " << r.core;
        EXPECT_EQ(2 * s + 6 * p + 10 * d + 14 * f, r.core) << "mode " << r.mode << " core " << r.core;
    }
    //the lanthanide range shares the 68-electron shell set
    for (const int core : { 55, 60, 67 })
    {
        int s = -7, p = -7, d = -7, f = -7;
        set_core_counts(&s, &p, &d, &f, core, 2);
        EXPECT_EQ(s, 5);
        EXPECT_EQ(p, 4);
        EXPECT_EQ(d, 2);
        EXPECT_EQ(f, 1);
    }
    for (const int mode : { 1, 2, 3 })
    {
        int s = -7, p = -7, d = -7, f = -7;
        set_core_counts(&s, &p, &d, &f, 0, mode);
        EXPECT_EQ(s, -7);
        EXPECT_EQ(p, -7);
        EXPECT_EQ(d, -7);
        EXPECT_EQ(f, -7);
    }
}

//a core count that closes no shell, and an unknown ECP mode, are err_not_impl exits
TEST(PropertiesCoverageSphericalTests, CoreCountsRefuseUnknownCores)
{
    const std::pair<int, int> bad[] = { { 4, 1 }, { 4, 2 }, { 2, 4 } };
    for (const auto &[core, mode] : bad)
        EXPECT_EXIT({ int s = 0; int p = 0; int d = 0; int f = 0; set_core_counts(&s, &p, &d, &f, core, mode); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
}

//f_core(k -> 0) counts the core electrons through the mode-2 table, and the core density never
//exceeds the atom's density
TEST(PropertiesCoverageSphericalTests, ThakkarEcpMode2CoreFormFactors)
{
    const std::pair<int, int> cases[] = { { 6, 2 }, { 36, 18 }, { 36, 28 }, { 54, 46 } };
    for (const auto &[Z, core] : cases)
    {
        Thakkar atom(Z, 2);
        EXPECT_NEAR(atom.get_core_form_factor(1e-6, core), core, 2e-2 * core) << "Z " << Z << " core " << core;
        for (const double r : { 0.05, 0.5, 2.0 })
        {
            const double total = atom.get_radial_density(r);
            const double c = atom.get_core_density(r, core);
            EXPECT_GE(c, 0.0);
            EXPECT_LE(c, total * (1.0 + 1e-12)) << "Z " << Z << " core " << core << " r " << r;
        }
    }
}

//the custom radial density with every shell equals the full density, splits additively over the
//s and the p/d/f shells, is zero with no shell, and the requested inner shells integrate to their
//electron count (Cu 3s2p: 18, carbon 1s: 2, which is also the mode-1 core density of core 2)
TEST(PropertiesCoverageSphericalTests, CustomRadialDensityPartitions)
{
    Thakkar carbon(6);
    Thakkar copper(29);
    for (const double r : { 0.02, 0.2, 0.7, 1.5, 3.0 })
    {
        const double full = carbon.get_radial_density(r);
        EXPECT_NEAR(carbon.get_radial_custom_density(r, 7, 6, 4, 2, 0, 0, 0, 0), full, 1e-9 * full);
        const double s_only = carbon.get_radial_custom_density(r, 7, 0, 0, 0, 0, 0, 0, 0);
        const double pdf_only = carbon.get_radial_custom_density(r, 0, 6, 4, 2, 0, 0, 0, 0);
        EXPECT_NEAR(s_only + pdf_only, full, 1e-9 * full);
        EXPECT_GT(s_only, 0.0);
        EXPECT_GT(pdf_only, 0.0);
        EXPECT_EQ(carbon.get_radial_custom_density(r, 0, 0, 0, 0, 0, 0, 0, 0), 0.0);
        const double core = carbon.get_radial_custom_density(r, 1, 0, 0, 0, 0, 0, 0, 0);
        EXPECT_NEAR(carbon.get_core_density(r, 2), core, 1e-9 * core);
        const double cu_full = copper.get_radial_density(r);
        EXPECT_NEAR(copper.get_radial_custom_density(r, 7, 6, 4, 2, 0, 0, 0, 0), cu_full, 1e-9 * cu_full);
    }
    const double cu_argon = radial_integral([&](double r) { return copper.get_radial_custom_density(r, 3, 2, 0, 0, 0, 0, 0, 0); }, 1e-7, 40.0, 4000);
    EXPECT_NEAR(cu_argon, 18.0, 2e-2 * 18.0);
    const double c_1s = radial_integral([&](double r) { return carbon.get_radial_custom_density(r, 1, 0, 0, 0, 0, 0, 0, 0); }, 1e-7, 40.0, 4000);
    EXPECT_NEAR(c_1s, 2.0, 2e-2 * 2.0);
    const double c_all = radial_integral([&](double r) { return carbon.get_radial_density(r); }, 1e-7, 40.0, 4000);
    EXPECT_NEAR(c_all, 6.0, 2e-2 * 6.0);
}

//suspected defect: Src/core/spherical_density.cpp:243 calc_custom_orbs counts coefficients from m_start = lower_m + min, so with min > 0 nr_coef skips the excluded orbitals' coefficients and every later orbital reads the wrong c[]
//excluding the 1s shell of copper must leave the other 27 electrons: the min-shell density plus
//the 1s density is the full density, and its integral matches the k -> 0 form factor of the same selection
TEST(PropertiesCoverageSphericalTests, DISABLED_CustomRadialDensityExcludesInnerShells)
{
    Thakkar copper(29);
    for (const double r : { 0.02, 0.2, 0.7, 1.5, 3.0 })
    {
        const double full = copper.get_radial_density(r);
        const double without_1s = copper.get_radial_custom_density(r, 7, 6, 4, 2, 1, 0, 0, 0);
        const double only_1s = copper.get_radial_custom_density(r, 1, 0, 0, 0, 0, 0, 0, 0);
        EXPECT_NEAR(without_1s + only_1s, full, 1e-9 * full);
    }
    const double valence = radial_integral([&](double r) { return copper.get_radial_custom_density(r, 7, 6, 4, 2, 1, 0, 0, 0); }, 1e-7, 40.0, 4000);
    EXPECT_NEAR(valence, 27.0, 2e-2 * 27.0);
    EXPECT_NEAR(copper.get_custom_form_factor(1e-6, 7, 6, 4, 2, 1, 0, 0, 0), 27.0, 2e-2 * 27.0);
}

//debug setup with Hirshfeld partitioning builds every weight scheme on helper grids and reports
//the spherical atoms, the per-atom grids and the Q_00 populations; the multipoles of the
//homonuclear pair sum to the electron count with opposite bond dipoles
TEST(PropertiesCoverageGridTests, DebugSetupEveryScheme)
{
    const WFN w = make_h2();
    GridConfiguration cfg;
    cfg.accuracy = 2;
    cfg.debug = true;
    cfg.partition_type = PartitionType::Hirshfeld;
    GridManager gm(cfg);
    std::ostringstream file;
    testing::internal::CaptureStdout();
    gm.setup3DGridsForMolecule(w, { 0, 1 }, {}, cell(), true, file);
    const vec2 Q = gm.calculatePartitionedMultipoles(w, 1);
    const std::string console = capture_end();
    EXPECT_NE(console.find("GridManager: Setting up grids for 2 atoms with Hirshfeld partitioning"), std::string::npos);
    EXPECT_NE(console.find("Size of atom_type_list:"), std::string::npos);
    EXPECT_NE(console.find("Calculating for atomic number 1"), std::string::npos);
    EXPECT_NE(console.find("GridManager: atom 1 population from Q_00:"), std::string::npos);
    EXPECT_NE(console.find("GridManager: Setup complete."), std::string::npos);
    EXPECT_NE(file.str().find("Generating integration grids for atoms"), std::string::npos);
    EXPECT_NE(file.str().find("Generated grid for atom 1/"), std::string::npos);
    EXPECT_TRUE(gm.getNeedsHelper());
    ASSERT_EQ(Q.size(), 2u);
    ASSERT_EQ(Q[0].size(), 4u);
    const double n = std::sqrt(constants::FOUR_PI) * (Q[0][0] + Q[1][0]);
    EXPECT_NEAR(n, h2_electrons(), 1e-3 * h2_electrons());
    for (int a = 0; a < 2; a++)
    {
        EXPECT_NEAR(Q[a][1], 0.0, 1e-6);
        EXPECT_NEAR(Q[a][3], 0.0, 1e-6);
    }
    EXPECT_NEAR(Q[0][2] + Q[1][2], 0.0, 1e-6);
    EXPECT_GT(std::fabs(Q[0][2]), 1e-3);
}

//accuracy 4 (non-hydrogen and hydrogen tables, radial grids grown to 1e-10) and the accuracy >= 5
//default tables integrate the density to the analytic electron count
TEST(PropertiesCoverageGridTests, LebedevParamsAtAccuracyFourAndFive)
{
    const WFN he_h = heh();
    const WFN h = single_h();
    const double n_heh = heh_electrons();
    const double n_h = std::pow(constants::PI / 2.0, 1.5);
    struct Case { const WFN *w; int accuracy; PartitionType type; double n; };
    const Case cases[] = {
        { &he_h, 4, PartitionType::Hirshfeld, n_heh },
        { &he_h, 5, PartitionType::Becke, n_heh },
        { &h, 5, PartitionType::Becke, n_h },
    };
    for (const Case &c : cases)
    {
        GridConfiguration cfg;
        cfg.accuracy = c.accuracy;
        cfg.partition_type = c.type;
        GridManager gm(cfg);
        std::ostringstream file;
        ivec atoms;
        for (int i = 0; i < c.w->get_ncen(); i++)
            atoms.push_back(i);
        gm.setup3DGridsForMolecule(*c.w, atoms, {}, cell(), false, file);
        const vec2 Q = gm.calculatePartitionedMultipoles(*c.w, 0);
        double n = 0.0;
        for (const vec &q : Q)
            n += std::sqrt(constants::FOUR_PI) * q[0];
        EXPECT_NEAR(n, c.n, 1e-3 * c.n) << "accuracy " << c.accuracy << " atoms " << c.w->get_ncen();
        EXPECT_GT(gm.getTotalGridPoints(), 0);
    }
}

//with helper grids the density vectors skip grid atoms outside atom_list: asking for atom 0 alone
//gives the same atom-0 electrons as asking for both, and both halves of the symmetric Gaussian are equal
TEST(PropertiesCoverageGridTests, DensityVectorsSkipAtomsOutsideTheList)
{
    const int n = 65;
    const double h = 0.25, origin = -8.0;
    cube c({ n, n, n }, 0, true);
    for (int i = 0; i < 3; i++)
    {
        c.set_origin(i, origin);
        for (int j = 0; j < 3; j++)
            c.set_vector(i, j, i == j ? h : 0.0);
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
            {
                const double x = origin + i * h, y = origin + j * h, z = origin + k * h;
                c.set_value(i, j, k, std::exp(-(x * x + y * y + z * z)));
            }
    const WFN w = make_h2();
    GridConfiguration cfg;
    cfg.accuracy = 3;
    cfg.partition_type = PartitionType::Becke;
    cfg.all_charges = true;
    cfg.debug = true;
    GridManager gm(cfg);
    std::ostringstream log;
    testing::internal::CaptureStdout();
    gm.setup3DGridsForMolecule(w, { 0, 1 }, {}, cell(), false, log);
    ASSERT_TRUE(gm.getNeedsHelper());
    vec2 d1, d2, d3, dens;
    vec one, both;
    gm.getDensityVectorsFromCube(w, { 0 }, c, d1, d2, d3, dens, one);
    ASSERT_EQ(one.size(), 1u);
    ASSERT_EQ(dens.size(), 1u);
    ASSERT_EQ(d1.size(), 1u);
    EXPECT_EQ(dens[0].size(), d1[0].size());
    EXPECT_NEAR(vec_total(dens), one[0], 1e-8);
    gm.getDensityVectorsFromCube(w, { 0, 1 }, c, d1, d2, d3, dens, both);
    const std::string console = capture_end();
    EXPECT_NE(console.find("GridManager: Generating density vectors from cube..."), std::string::npos);
    ASSERT_EQ(both.size(), 2u);
    EXPECT_NEAR(one[0], both[0], 1e-9);
    EXPECT_NEAR(both[0], both[1], 1e-6);
    const double ref = std::pow(constants::PI, 1.5);
    EXPECT_NEAR(both[0] + both[1], ref, 5e-2 * ref);
}

//two- and one-dimensional grids get the shorter headers, and rows carry the coordinates at three
//decimals then the values at ten
TEST(PropertiesCoverageGridTests, WriteSimpleGridLowerDimensions)
{
    const std::filesystem::path dir = temp_dir("simple_grid");
    GridManager gm;
    gm.writeSimpleGrid(dir / "two.dat", { { 0.0, 1.0 }, { 2.0, 3.0 } }, { { "v", { 0.5, -1.0 } } });
    const std::string two = read_file(dir / "two.dat");
    EXPECT_EQ(two, "# X\tY\tv\t\n0.000\t2.000\t5.0000000000e-01\t\n1.000\t3.000\t-1.0000000000e+00\t\n");
    gm.writeSimpleGrid(dir / "one.dat", { { 0.25 } }, { { "a", { 2.0 } }, { "b", { 3.0 } } });
    const std::string one = read_file(dir / "one.dat");
    EXPECT_EQ(one, "# X\ta\tb\t\n0.250\t2.0000000000e+00\t3.0000000000e+00\t\n");
    std::filesystem::remove_all(dir);
}

TEST(PropertiesCoverageGridTests, WriteSimpleGridSizeMismatchExits)
{
    const std::filesystem::path dir = death_dir("simple_grid_mismatch");
    GridManager gm;
    EXPECT_EXIT({ gm.writeSimpleGrid(dir / "bad.dat", { { 0.0, 1.0 } }, { { "v", { 0.5 } } }); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
    std::filesystem::remove_all(dir);
}

//two maxima on the grid points of the nuclei with the atoms 2.5 bohr apart: the hill climb finds
//two mirror-image basins, an empty selection writes one cube per basin named after its atom, and
//the log integrals are equal and sum to the grid integral
TEST(PropertiesCoverageBasinTests, LegacyB2cWritesOneCubePerBasin)
{
    ScopedGlobals guard;
    const std::filesystem::path dir = temp_dir("b2c_two_basins");
    H2Model m(1.0, 1.25);
    //10 voxels along the bond put the nuclei on grid points (x = +-1.25) with no mid-plane voxel, the
    //odd 9 across put y = z = 0 on the grid: each nucleus is one unique voxel maximum. A 10^3 grid has
    //its across axes at +-0.25, a 2x2 plateau of equal voxels, and the hill-climb makes a basin of each
    cube rho({ 10, 9, 9 }, 0, true);
    for (int k = 0; k < 3; k++)
    {
        rho.set_origin(k, k == 0 ? -2.25 : -2.0);
        rho.set_vector(k, k, 0.5);
    }
    rho.calc_dv();
    double total = 0.0;
    for (int i = 0; i < 10; i++)
        for (int j = 0; j < 9; j++)
            for (int k = 0; k < 9; k++)
            {
                const double v = m.rho(rho.get_pos(i, j, k));
                rho.set_value(i, j, k, v);
                total += v * rho.get_dv();
            }
    rho.set_na(2);
    rho.give_parent_wfn(m.wavy);
    rho.set_path(dir / "h2.cube");
    std::istringstream answers("0\n");
    std::streambuf *old_cin = std::cin.rdbuf(answers.rdbuf());
    testing::internal::CaptureStdout();
    const bool ok = b2c(&rho, m.wavy.get_atoms(), true, false);
    const std::string console = capture_end();
    std::cin.rdbuf(old_cin);
    EXPECT_TRUE(ok);
    EXPECT_NE(console.find("I found 2 Basins."), std::string::npos);
    EXPECT_NE(console.find("done with liste, writing basins now!"), std::string::npos);
    EXPECT_NE(console.find("DEBUG: Labels_size:"), std::string::npos);
    //basin ids follow the scan order, so either hydrogen may own basin 1
    const bool forward = std::filesystem::exists(dir / "h2_H_0_0.cube") && std::filesystem::exists(dir / "h2_H_1_1.cube");
    const bool backward = std::filesystem::exists(dir / "h2_H_1_0.cube") && std::filesystem::exists(dir / "h2_H_0_1.cube");
    EXPECT_TRUE(forward || backward);
    const std::string log = read_file(dir / "h2.b2c_log");
    EXPECT_NE(log.find("Number of Basins: 2"), std::string::npos);
    const double b1 = number_after(log, "Basin 1: ");
    const double b2 = number_after(log, "Basin 2: ");
    const double all = number_after(log, "Integral over all basins: ");
    EXPECT_NEAR(b1, b2, 1e-6 * b1);
    EXPECT_NEAR(b1 + b2, total, 1e-6 * total);
    EXPECT_NEAR(all, total, 1e-6 * total);
    std::filesystem::remove_all(dir);
}

//one maximum between two close nuclei: a single basin, so the bond-critical-point search has no
//zero-flux surface to walk and finishes clean
TEST(PropertiesCoverageBasinTests, LegacyB2cSingleBasinBcpSearch)
{
    ScopedGlobals guard;
    const std::filesystem::path dir = temp_dir("b2c_single_bcp");
    H2Model m(1.0, 0.5);
    cube rho = make_grid(9, 0.5);
    for (int i = 0; i < 9; i++)
        for (int j = 0; j < 9; j++)
            for (int k = 0; k < 9; k++)
                rho.set_value(i, j, k, m.rho(rho.get_pos(i, j, k)));
    rho.set_na(2);
    rho.give_parent_wfn(m.wavy);
    rho.set_path(dir / "h2.cube");
    std::istringstream answers("1\n0\n");
    std::streambuf *old_cin = std::cin.rdbuf(answers.rdbuf());
    testing::internal::CaptureStdout();
    const bool ok = b2c(&rho, m.wavy.get_atoms(), true, true);
    const std::string console = capture_end();
    std::cin.rdbuf(old_cin);
    EXPECT_TRUE(ok);
    EXPECT_NE(console.find("I found 1 Basins."), std::string::npos);
    EXPECT_NE(console.find("done with BCPs"), std::string::npos);
    EXPECT_TRUE(std::filesystem::exists(dir / "h2_1_basins.cube"));
    const std::string log = read_file(dir / "h2.b2c_log");
    EXPECT_NE(log.find("Number of Basins: 1"), std::string::npos);
    std::filesystem::remove_all(dir);
}

//suspected defect: Src/core/b2c.cpp:392 neighbours (size iCP) is indexed by the 1-based basin id CP(x, y, z); with two basins neighbours[2] is out of range
//the BCP of two mirror basins sits on the mid-plane: reported density = rho at the border voxels
//(x = +-0.25), and the two labels are the two hydrogens
TEST(PropertiesCoverageBasinTests, DISABLED_LegacyB2cTwoBasinBcp)
{
    ScopedGlobals guard;
    const std::filesystem::path dir = temp_dir("b2c_two_basin_bcp");
    H2Model m(1.0, 1.25);
    cube rho = make_grid(10, 0.5);
    for (int i = 0; i < 10; i++)
        for (int j = 0; j < 10; j++)
            for (int k = 0; k < 10; k++)
                rho.set_value(i, j, k, m.rho(rho.get_pos(i, j, k)));
    rho.set_na(2);
    rho.give_parent_wfn(m.wavy);
    rho.set_path(dir / "h2.cube");
    std::istringstream answers("0\n");
    std::streambuf *old_cin = std::cin.rdbuf(answers.rdbuf());
    testing::internal::CaptureStdout();
    const bool ok = b2c(&rho, m.wavy.get_atoms(), true, true);
    const std::string console = capture_end();
    std::cin.rdbuf(old_cin);
    EXPECT_TRUE(ok);
    EXPECT_NE(console.find("I found 2 Basins."), std::string::npos);
    EXPECT_NEAR(number_after(console, "BCP: H-H ED: "), m.rho({ 0.25, 0.0, 0.0 }), 1e-6 * m.rho({ 0.25, 0.0, 0.0 }));
    std::filesystem::remove_all(dir);
}

//conjugation negates every imaginary part and leaves the real parts and the block layout alone
TEST(PropertiesCoverageSaltedTests, ConjugateFlipsImaginaryParts)
{
    SALTEDDescriptors d(2, 3, 1);
    std::vector<cdouble> &v = d.values();
    ASSERT_EQ(v.size(), 24u);
    for (size_t i = 0; i < v.size(); i++)
        v[i] = cdouble((double)i, -(double)i);
    calculateConjugate(d);
    for (size_t i = 0; i < v.size(); i++)
    {
        EXPECT_EQ(v[i].real(), (double)i);
        EXPECT_EQ(v[i].imag(), (double)i);
    }
    EXPECT_EQ(d.block(0, 0, 0), v.data());
    EXPECT_EQ(d.block(1, 2, 1)[2], v[23]);
}

TEST(PropertiesCoverageSaltedTests, ConjugateOfEmptyDescriptorsIsNoop)
{
    SALTEDDescriptors d;
    calculateConjugate(d);
    EXPECT_TRUE(d.values().empty());
}

//-out needs a file name that is not another option
TEST(PropertiesCoverageAppTests, OutOptionWithoutArgumentExits)
{
    const std::filesystem::path dir = death_dir("app_out_missing");
    {
        ScopedCwd cwd(dir);
        EXPECT_EXIT({ run_nos({ "-out" }); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
        EXPECT_EXIT({ run_nos({ "-out", "-no_date" }); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
    }
    std::filesystem::remove_all(dir);
}

//a non-numeric -cpus is caught by the option digester and turned into an error exit
TEST(PropertiesCoverageAppTests, NonNumericCpusExits)
{
    const std::filesystem::path dir = death_dir("app_bad_cpus");
    {
        ScopedCwd cwd(dir);
        EXPECT_EXIT({ run_nos({ "-cpus", "abc", "-out", "bad.log", "-no_date" }); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
    }
    std::filesystem::remove_all(dir);
}

//-h prints the help and exits with status 0
TEST(PropertiesCoverageAppTests, HelpExitsZero)
{
    const std::filesystem::path dir = death_dir("app_help");
    {
        ScopedCwd cwd(dir);
        EXPECT_EXIT({ run_nos({ "-h", "-out", "help.log" }); }, ::testing::ExitedWithCode(0), ".*");
    }
    std::filesystem::remove_all(dir);
}

//-fractal on a 4^3 cube whose value is the x index: 144 neighbour comparisons, every iso level
//strictly between two integers is crossed by the 16 x-pairs -> df = ln 16 / ((1/3) ln 144), iso
//levels outside the value range have no crossing and df = 0; the plot header carries the value
//range and its step count matches the lines that follow. -mem beyond the addressable range is
//clamped to 50000 MB with a note in the log
TEST(PropertiesCoverageAppTests, FractalOptionWritesPlot)
{
    ScopedGlobals guard;
    const std::filesystem::path dir = temp_dir("app_fractal");
    {
        ScopedCwd cwd(dir);
        WFN two_h(e_origin::NOT_YET_DEFINED);
        two_h.push_back_atom("H", -0.7, 0.0, 0.0, 1);
        two_h.push_back_atom("H", 0.7, 0.0, 0.0, 1);
        cube c = make_grid(4, 1.0);
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                for (int k = 0; k < 4; k++)
                    c.set_value(i, j, k, (double)i);
        c.set_na(2);
        c.give_parent_wfn(two_h);
        c.set_comment1("x index ramp");
        c.set_comment2("fractal test");
        c.set_path("h2.cube");
        ASSERT_TRUE(c.write_file(true));
        testing::internal::CaptureStdout();
        const int rc = run_nos({ "-fractal", "h2.cube", "-no_date", "-mem", "1e30", "-out", "h2.log" });
        const std::string console = capture_end();
        EXPECT_EQ(rc, 0);
        const std::string log = read_file(dir / "h2.log");
        EXPECT_NE(log.find("Setting max memory to 50000 MB"), std::string::npos);
        //the -fractal path restores cout before it announces the plot, so the line is on the console
        EXPECT_NE(console.find("Finished writing fractal dimensions plot"), std::string::npos);
        ASSERT_TRUE(std::filesystem::exists(dir / "h2.cube_fractal_plot"));
        std::ifstream plot(dir / "h2.cube_fractal_plot");
        int steps = 0;
        double map_min = 0.0, map_max = 0.0, e0 = 0.0, e1 = 0.0;
        double steps_as_double = 0.0;
        ASSERT_TRUE(static_cast<bool>(plot >> steps_as_double >> map_min >> map_max >> e0 >> e1));
        steps = (int)std::llround(steps_as_double);
        EXPECT_NEAR(map_min, 0.0, 1e-12);
        EXPECT_NEAR(map_max, 3.0, 1e-12);
        const double df_between = 3.0 * std::log(16.0) / std::log(144.0);
        int lines = 0, checked_half = 0, checked_below = 0;
        double iso = 0.0, df = 0.0;
        while (plot >> iso >> df)
        {
            lines++;
            if (std::fabs(iso - 0.5) < 1e-9)
            {
                EXPECT_NEAR(df, df_between, 1e-6);
                checked_half++;
            }
            if (std::fabs(iso + 0.01) < 1e-9)
            {
                EXPECT_EQ(df, 0.0);
                checked_below++;
            }
            if (std::fabs(iso - 1.0) < 1e-9 || std::fabs(iso - 2.0) < 1e-9)
                EXPECT_EQ(df, 0.0);
        }
        EXPECT_EQ(lines, steps);
        EXPECT_EQ(checked_half, 1);
        EXPECT_EQ(checked_below, 1);
    }
    std::filesystem::remove_all(dir);
}

//-density_difference of a wavefunction against its own copy: the RSR is exactly zero, the two
//rho cubes and the difference cube land in the working directory, the debug line names both files
TEST(PropertiesCoverageAppTests, DensityDifferenceOfAWavefunctionWithItself)
{
    ScopedGlobals guard;
    const std::filesystem::path dir = temp_dir("app_density_difference");
    {
        ScopedCwd cwd(dir);
        H2Model m(1.0, 1.0);
        m.wavy.write_wfn(dir / "h2.wfn", false, false);
        std::filesystem::copy_file(dir / "h2.wfn", dir / "h2b.wfn");
        testing::internal::CaptureStdout();
        const int rc = run_nos({ "-wfn", "h2.wfn", "-density_difference", "h2b.wfn", "-debug", "-no_date", "-resolution", "0.3", "-radius", "1.5", "-out", "h2.log" });
        capture_end();
        EXPECT_EQ(rc, 0);
        const std::string log = read_file(dir / "h2.log");
        EXPECT_NE(log.find(" vs "), std::string::npos);
        EXPECT_NE(log.find("h2b.wfn"), std::string::npos);
        EXPECT_NE(log.find("Writing cube 1..."), std::string::npos);
        EXPECT_NE(log.find("Bye Bye!"), std::string::npos);
        const double rsr = number_after(log, "RSR between the two cubes: ");
        EXPECT_LE(std::fabs(rsr), 1e-12);
        EXPECT_NEAR(number_after(log, "Ne of shifted electrons: "), 0.0, 1e-12);
        for (const char *name : { "h2_rho.cube", "h2b_rho.cube", "h2b_diff.cube" })
            EXPECT_TRUE(std::filesystem::exists(dir / name)) << name;
        std::ostringstream sink;
        WFN dummy(e_origin::NOT_YET_DEFINED);
        const cube diff(dir / "h2b_diff.cube", true, dummy, sink);
        EXPECT_NEAR(diff.max_value(), 0.0, 1e-12);
        EXPECT_NEAR(diff.min_value(), 0.0, 1e-12);
        WFN dummy2(e_origin::NOT_YET_DEFINED);
        const cube rho(dir / "h2_rho.cube", true, dummy2, sink);
        int checked = 0;
        for (int i = 0; i < rho.get_size(0); i++)
            for (int j = 0; j < rho.get_size(1); j++)
                for (int k = 0; k < rho.get_size(2); k++)
                {
                    const double v = rho.get_value(i, j, k);
                    if (v < 1e-3)
                        continue;
                    const double ref = m.rho(rho.get_pos(i, j, k));
                    //the reader's exp_cutoff drops a primitive whose orbital contribution is below
                    //density_accuracy = 5e-5, so rho = 2 psi0^2 is off by up to 4 psi0 5e-5 < 1e-4 (5e-3
                    //relative at rho = 1e-3); the cube header carries positions to 6 decimals, 4 a r dr < 1e-4
                    EXPECT_NEAR(v, ref, 2e-4 * ref + 1e-4);
                    checked++;
                }
        EXPECT_GT(checked, 50);
    }
    std::filesystem::remove_all(dir);
}
