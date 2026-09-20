#include "pch.h"
#include "core/convenience.h"
#include "core/constants.h"
#include "core/atoms.h"
#include "core/wfn_class.h"
#include "core/basis_set.h"
#include "core/fchk.h"
#include "core/cube.h"
#include "core/npy.h"
#include "core/nos_math.h"
#include "core/integrator.h"
#include "core/integration_params.h"
#include "core/libCintKernels.h"
#include "core/libCintMain.h"
#include "core/SALTED_utilities.h"
#include <occ/qm/hf.h>
#include <occ/qm/scf.h>
#include <spdlog/spdlog.h>
#undef I
#include <numeric>
#include <array>
#include <complex>

//Non-static helpers of integrator.cpp and basis_set.cpp that have no header declaration.
vec einsum_ijk_ij_p(const dMatrix3& v1, const dMatrix2& v2);
vec reorder_p(vec coefs_in, WFN aux_basis);
std::vector<double> dedup_exponents(std::vector<double> exps, double tol);
std::vector<int> pivoted_cholesky(dMatrix2& A, double threshold);
std::vector<double> prune_element_candidates_for_L(int L, const std::vector<double>& exponents, double threshold);

namespace
{
    //Swaps the buffer of one standard stream for the lifetime of the object.
    struct StreamCapture
    {
        std::ostringstream buffer;
        std::ostream& stream;
        std::streambuf* old;
        explicit StreamCapture(std::ostream& s) : stream(s), old(s.rdbuf(buffer.rdbuf())) {}
        ~StreamCapture() { stream.rdbuf(old); }
        std::string str() const { return buffer.str(); }
    };

    //Feeds a fixed string to std::cin for the lifetime of the object.
    struct CinFeed
    {
        std::istringstream in;
        std::streambuf* old;
        explicit CinFeed(const std::string& text) : in(text), old(std::cin.rdbuf(in.rdbuf())) {}
        ~CinFeed() { std::cin.rdbuf(old); }
    };

    //Unique scratch directory per test; the cwd is restored and the directory removed on success.
    struct Scratch
    {
        std::filesystem::path dir;
        std::filesystem::path old_cwd;
        explicit Scratch(const std::string& name) : old_cwd(std::filesystem::current_path())
        {
            dir = std::filesystem::temp_directory_path() / ("NoSpherA2_FittingIoCoverage_" + name);
            std::filesystem::remove_all(dir);
            std::filesystem::create_directories(dir);
        }
        ~Scratch()
        {
            std::error_code ec;
            std::filesystem::current_path(old_cwd, ec);
            if (!::testing::Test::HasFailure())
                std::filesystem::remove_all(dir, ec);
        }
    };

    void write_text(const std::filesystem::path& p, const std::string& text)
    {
        std::ofstream f(p, std::ios::binary);
        f << text;
    }

    std::string read_text(const std::filesystem::path& p)
    {
        std::ifstream f(p, std::ios::binary);
        std::stringstream ss;
        ss << f.rdbuf();
        return ss.str();
    }

    std::vector<std::string> split_lines(const std::string& text)
    {
        std::vector<std::string> lines;
        std::stringstream ss(text);
        std::string line;
        while (std::getline(ss, line))
        {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            lines.push_back(line);
        }
        return lines;
    }

    std::filesystem::path epoxide_fixture()
    {
        const std::filesystem::path p = nos_test_repo_root() / "tests" / "epoxide_gbw" / "epoxide.gbw";
        return std::filesystem::exists(p) ? p : std::filesystem::path();
    }

    std::vector<std::shared_ptr<BasisSet>> combo_basis()
    {
        return { BasisSetLibrary::get_basis_set("combo_basis_fit") };
    }

    //H2+ doublet from occ, as SaltedFchkTests builds it.
    WFN occ_h2_plus()
    {
        spdlog::set_level(spdlog::level::err);
        const std::vector<occ::core::Atom> atoms{ { 1, 0.0, 0.0, -0.7 }, { 1, 0.0, 0.0, 0.7 } };
        std::vector<occ::gto::Shell> shells;
        for (const auto& at : atoms)
            for (int l = 0; l <= 1; l++)
            {
                shells.emplace_back(l, std::vector<double>{ l == 0 ? 1.2 : 0.9 }, std::vector<vec>{ { 1.0 } }, std::array<double, 3>{ at.x, at.y, at.z });
                shells.back().kind = occ::gto::Shell::Kind::Spherical;
                shells.back().incorporate_shell_norm();
            }
        occ::gto::AOBasis basis(atoms, shells, "sp");
        basis.set_pure(true);
        occ::qm::HartreeFock hf(basis);
        occ::qm::SCF<occ::qm::HartreeFock> scf(hf, occ::qm::SpinorbitalKind::Unrestricted);
        scf.set_charge_multiplicity(1, 2);
        scf.compute_initial_guess();
        scf.compute_scf_energy();
        WFN w(scf.wavefunction(), false);
        w.set_multi(2);
        return w;
    }

    //One population row per atom: the integral of every s aux function of that atom, as
    //integrator.cpp builds them for the charge restraints and the total-electron constraint.
    vec2 population_rows(const WFN& aux)
    {
        const aux_density_table t(aux.get_atoms());
        vec2 rows(t.n_at, vec(t.n_coef, 0.0));
        for (int a = 0; a < t.n_at; a++)
            for (int s = t.sh_start[a]; s < t.sh_start[a + 1]; s++)
                if (t.sh_l[s] == 0)
                    rows[a][t.coef_off[s]] = t.shell_population_integral(s);
        return rows;
    }

    double fitted_electrons(const vec2& rows, const vec& c)
    {
        double total = 0.0;
        for (const vec& row : rows)
            for (size_t i = 0; i < row.size(); i++)
                total += row[i] * c[i];
        return total;
    }

    //max_i |(H c - g)_i| relative to the largest row scale max_i (|g_i| + sum_j |H_ij c_j|).
    double relative_residual(const vec& H, const vec& c, const vec& g)
    {
        const size_t n = c.size();
        double worst = 0.0, scale = 0.0;
        for (size_t i = 0; i < n; i++)
        {
            double hc = 0.0, abs_hc = 0.0;
            for (size_t j = 0; j < n; j++)
            {
                hc += H[i * n + j] * c[j];
                abs_hc += std::fabs(H[i * n + j] * c[j]);
            }
            worst = std::max(worst, std::fabs(hc - g[i]));
            scale = std::max(scale, std::fabs(g[i]) + abs_hc);
        }
        return worst / scale;
    }

    struct CoulombSystem
    {
        vec H;
        vec g;
        size_t n = 0;
    };

    CoulombSystem coulomb_system(const WFN& wave, const WFN& aux)
    {
        CoulombSystem s;
        Int_Params np(wave), ap(aux);
        compute2C<Coulomb2C_SPH>(ap, s.H);
        computeRho<Coulomb3C_SPH>(np, ap, wave.get_dm(), s.g);
        s.n = (size_t)ap.get_nao();
        return s;
    }
}

//----------------------------------------------------------------------------------------------
//integrator.cpp
//----------------------------------------------------------------------------------------------

//v1(i,j,p) = 1 + i + 2j + 4p and v2(i,j) = 1 + 2i + j on a 2x2x2 / 2x2 grid:
//rho[0] = 1*1 + 3*2 + 2*3 + 4*4 = 29, rho[1] = 5*1 + 7*2 + 6*3 + 8*4 = 69
TEST(FittingIoCoverageIntegratorTests, EinsumContractsFirstTwoIndices)
{
    const dMatrix3 v1 = reshape<dMatrix3>(vec{ 1, 5, 3, 7, 2, 6, 4, 8 }, Shape3D(2, 2, 2));
    const dMatrix2 v2 = reshape<dMatrix2>(vec{ 1, 2, 3, 4 }, Shape2D(2, 2));
    const vec rho = einsum_ijk_ij_p(v1, v2);
    ASSERT_EQ(rho.size(), 2u);
    EXPECT_DOUBLE_EQ(rho[0], 29.0);
    EXPECT_DOUBLE_EQ(rho[1], 69.0);
}

//One atom with an s, a p and a d shell (l stored as 0/1/2): only the p triple is rotated
//(x,y,z) -> (y,z,x) in the coefficient vector, the s and the five d entries stay in place.
TEST(FittingIoCoverageIntegratorTests, ReorderPRotatesOnlyPShells)
{
    WFN aux(e_origin::NOT_YET_DEFINED);
    aux.push_back_atom("H", 0.0, 0.0, 0.0, 1);
    ASSERT_TRUE(aux.push_back_atom_basis_set(0, 1.0, 1.0, 0, 0));
    ASSERT_TRUE(aux.push_back_atom_basis_set(0, 0.8, 1.0, 1, 1));
    ASSERT_TRUE(aux.push_back_atom_basis_set(0, 0.6, 1.0, 2, 2));
    ASSERT_EQ(aux.get_atom_shell_count(0), 3);
    const vec in{ 10, 21, 22, 23, 30, 31, 32, 33, 34 };
    const vec out = reorder_p(in, aux);
    const vec expected{ 10, 22, 23, 21, 30, 31, 32, 33, 34 };
    ASSERT_EQ(out.size(), expected.size());
    for (size_t i = 0; i < out.size(); i++)
        EXPECT_DOUBLE_EQ(out[i], expected[i]) << i;
}

//Defaults map to an unrestrained fit with the CONFIG defaults untouched; a non-negative lmax
//switches the charge restraints on, lmax > 0 the multipole restraints, and every partition
//scheme maps to its charge scheme (Becke/RI fall back to Hirshfeld).
TEST(FittingIoCoverageIntegratorTests, ConfigFromOptionsMapsEveryScheme)
{
    options opt;
    DensityFitting::CONFIG cfg = DensityFitting::config_from_options(opt);
    EXPECT_FALSE(cfg.analyze_quality);
    EXPECT_FALSE(cfg.restrain_charges);
    EXPECT_FALSE(cfg.restrain_multipoles);
    EXPECT_FALSE(cfg.use_tikhonov);
    EXPECT_EQ(cfg.multipole_lmax, -1);
    EXPECT_EQ(cfg.charge_scheme, DensityFitting::CHARGE_SCHEME::TFVC);

    opt.debug = true;
    opt.multipole_lmax = 0;
    opt.multipole_strength = 2.5;
    opt.multipole_partition = false;
    opt.multipole_scheme = PartitionType::TFVC;
    cfg = DensityFitting::config_from_options(opt);
    EXPECT_TRUE(cfg.analyze_quality);
    EXPECT_TRUE(cfg.restrain_charges);
    EXPECT_FALSE(cfg.restrain_multipoles);
    EXPECT_FALSE(cfg.partition_restraints);
    EXPECT_EQ(cfg.multipole_lmax, 0);
    EXPECT_DOUBLE_EQ(cfg.multipole_strength, 2.5);
    EXPECT_EQ(cfg.charge_scheme, DensityFitting::CHARGE_SCHEME::TFVC);

    opt.multipole_lmax = 2;
    opt.multipole_partition = true;
    opt.multipole_scheme = PartitionType::MBIS;
    cfg = DensityFitting::config_from_options(opt);
    EXPECT_TRUE(cfg.restrain_multipoles);
    EXPECT_TRUE(cfg.partition_restraints);
    EXPECT_EQ(cfg.multipole_lmax, 2);
    EXPECT_EQ(cfg.charge_scheme, DensityFitting::CHARGE_SCHEME::MBIS);

    opt.multipole_scheme = PartitionType::EMBIS;
    EXPECT_EQ(DensityFitting::config_from_options(opt).charge_scheme, DensityFitting::CHARGE_SCHEME::EMBIS);
    opt.multipole_scheme = PartitionType::Hirshfeld;
    EXPECT_EQ(DensityFitting::config_from_options(opt).charge_scheme, DensityFitting::CHARGE_SCHEME::HIRSHFELD);
    opt.multipole_scheme = PartitionType::Becke;
    EXPECT_EQ(DensityFitting::config_from_options(opt).charge_scheme, DensityFitting::CHARGE_SCHEME::HIRSHFELD);
    opt.multipole_scheme = PartitionType::RI;
    EXPECT_EQ(DensityFitting::config_from_options(opt).charge_scheme, DensityFitting::CHARGE_SCHEME::HIRSHFELD);
}

//Non-adaptive: every atom gets the base coefficient. Adaptive with n_aux = 100 and three atoms
//(H, C, O): base * (1 - 0.1 log10 100) * min(2, 1 + 0.1 sqrt 3) = 0.8 * 1.1732051 = 0.9385641,
//then times min(2, 1 + 0.02 Z) = 1.02 / 1.12 / 1.16.
TEST(FittingIoCoverageIntegratorTests, RestraintWeightsPlainAndAdaptive)
{
    WFN aux(e_origin::NOT_YET_DEFINED);
    aux.push_back_atom("H", 0.0, 0.0, 0.0, 1);
    aux.push_back_atom("C", 1.0, 0.0, 0.0, 6);
    aux.push_back_atom("O", 2.0, 0.0, 0.0, 8);
    {
        StreamCapture out(std::cout);
        const vec w = DensityFitting::restraint_weights(aux, 100, 0.5, false);
        ASSERT_EQ(w.size(), 3u);
        for (double v : w) EXPECT_DOUBLE_EQ(v, 0.5);
        EXPECT_NE(out.str().find("Setting charge-restraint row scale to: 0.500000"), std::string::npos) << out.str();
    }
    {
        StreamCapture out(std::cout);
        const vec w = DensityFitting::restraint_weights(aux, 100, 1.0, true);
        ASSERT_EQ(w.size(), 3u);
        const double coef = 0.8 * (1.0 + std::sqrt(3.0) * 0.1);
        EXPECT_NEAR(w[0], coef * 1.02, 1e-12);
        EXPECT_NEAR(w[1], coef * 1.12, 1e-12);
        EXPECT_NEAR(w[2], coef * 1.16, 1e-12);
    }
}

//Nuclear populations are Z; the Sanderson estimate follows the closed formula from the Allen
//electronegativities; an out-of-range scheme warns on stderr and falls back to nuclear.
TEST(FittingIoCoverageIntegratorTests, ExpectedPopulationsNuclearSandersonAndUnknown)
{
    const std::filesystem::path p = epoxide_fixture();
    if (p.empty()) GTEST_SKIP() << "tests/epoxide_gbw/epoxide.gbw not found";
    WFN wave(p);
    auto basis = combo_basis();
    const WFN aux = generate_aux_wfn(wave, basis);
    ASSERT_EQ(aux.get_ncen(), 7);

    const vec nuclear = DensityFitting::calculate_expected_populations(wave, aux, DensityFitting::CHARGE_SCHEME::NUCLEAR);
    ASSERT_EQ(nuclear.size(), 7u);
    for (int a = 0; a < 7; a++)
        EXPECT_DOUBLE_EQ(nuclear[a], (double)aux.get_atom_charge(a)) << a;
    EXPECT_NEAR(std::accumulate(nuclear.begin(), nuclear.end(), 0.0), 24.0, 1e-12);

    double chi_compound = 1.0;
    for (int a = 0; a < 7; a++)
        chi_compound *= constants::allen_electronegativities[wave.get_atom_charge(a) - 1];
    chi_compound = std::pow(chi_compound, 1.0 / 7.0);
    const vec sanderson = DensityFitting::calculate_expected_populations(wave, aux, DensityFitting::CHARGE_SCHEME::SANDERSON_ESTIMATE);
    ASSERT_EQ(sanderson.size(), 7u);
    for (int a = 0; a < 7; a++)
    {
        const double chi = constants::allen_electronegativities[aux.get_atom_charge(a) - 1];
        EXPECT_NEAR(sanderson[a], aux.get_atom_charge(a) + (chi_compound - chi) / (1.57 * std::sqrt(chi)), 1e-12) << a;
    }

    StreamCapture errcap(std::cerr);
    const vec unknown = DensityFitting::calculate_expected_populations(wave, aux, static_cast<DensityFitting::CHARGE_SCHEME>(99));
    ASSERT_EQ(unknown.size(), 7u);
    for (int a = 0; a < 7; a++)
        EXPECT_DOUBLE_EQ(unknown[a], nuclear[a]) << a;
    EXPECT_NE(errcap.str().find("Warning: Unknown charge scheme"), std::string::npos) << errcap.str();
}

//Mulliken populations sum to Tr(D S) = the electron count of the neutral molecule.
TEST(FittingIoCoverageIntegratorTests, ExpectedPopulationsMullikenSumsToElectronCount)
{
    const std::filesystem::path p = epoxide_fixture();
    if (p.empty()) GTEST_SKIP() << "tests/epoxide_gbw/epoxide.gbw not found";
    WFN wave(p);
    auto basis = combo_basis();
    const WFN aux = generate_aux_wfn(wave, basis);
    const vec pop = DensityFitting::calculate_expected_populations(wave, aux, DensityFitting::CHARGE_SCHEME::MULLIKEN);
    ASSERT_EQ(pop.size(), 7u);
    EXPECT_NEAR(std::accumulate(pop.begin(), pop.end(), 0.0), 24.0, 1e-5);
    for (int a = 0; a < 7; a++)
        EXPECT_GT(pop[a], 0.0) << a;
}

//Hirshfeld populations on the coarse grid: all positive and summing to the electron count.
TEST(FittingIoCoverageIntegratorTests, ExpectedPopulationsHirshfeldSumsToElectronCount)
{
    const std::filesystem::path p = epoxide_fixture();
    if (p.empty()) GTEST_SKIP() << "tests/epoxide_gbw/epoxide.gbw not found";
    WFN wave(p);
    auto basis = combo_basis();
    const WFN aux = generate_aux_wfn(wave, basis);
    StreamCapture out(std::cout);
    const vec pop = DensityFitting::calculate_expected_populations(wave, aux, DensityFitting::CHARGE_SCHEME::HIRSHFELD);
    ASSERT_EQ(pop.size(), 7u);
    EXPECT_NEAR(std::accumulate(pop.begin(), pop.end(), 0.0), 24.0, 0.2);
    for (int a = 0; a < 7; a++)
    {
        EXPECT_GT(pop[a], 0.0) << a;
        EXPECT_LT(pop[a], 10.0) << a;
    }
}

//The unrestrained Coulomb fit is the plain solve of (aux|aux) c = (aux|rho); the fitted
//density carries the 24 electrons of epoxide and the quality report is printed on request.
TEST(FittingIoCoverageIntegratorTests, UnrestrainedCoulombFitMatchesDirectSolve)
{
    const std::filesystem::path p = epoxide_fixture();
    if (p.empty()) GTEST_SKIP() << "tests/epoxide_gbw/epoxide.gbw not found";
    WFN wave(p);
    auto basis = combo_basis();
    const WFN aux = generate_aux_wfn(wave, basis);
    const CoulombSystem sys = coulomb_system(wave, aux);
    vec H = sys.H, c_ref = sys.g;
    ASSERT_EQ(solve_linear_system(H, sys.n, c_ref), 0);

    DensityFitting::CONFIG cfg;
    cfg.analyze_quality = true;
    StreamCapture out(std::cout);
    const vec c = DensityFitting::density_fit(wave, aux, cfg);
    ASSERT_EQ(c.size(), sys.n);
    double cmax = 0.0;
    for (double v : c_ref) cmax = std::max(cmax, std::fabs(v));
    for (size_t i = 0; i < sys.n; i++)
        EXPECT_NEAR(c[i], c_ref[i], 1e-9 * cmax) << i;
    EXPECT_NEAR(fitted_electrons(population_rows(aux), c), 24.0, 0.1);

    const std::string log = out.str();
    EXPECT_NE(log.find("=== Density Fitting ==="), std::string::npos);
    EXPECT_NE(log.find("Normal basis functions: 62"), std::string::npos);
    EXPECT_NE(log.find("Auxiliary basis functions: " + std::to_string(sys.n)), std::string::npos);
    EXPECT_NE(log.find("Metric: Coulomb"), std::string::npos);
    EXPECT_NE(log.find("Solving unrestrained linear system..."), std::string::npos);
    EXPECT_NE(log.find("=== Density Fitting Quality Analysis ==="), std::string::npos);
    EXPECT_NE(log.find("Atom 1 (Z="), std::string::npos);
    EXPECT_NE(log.find("Expected / Real total electrons: 24.000 / "), std::string::npos);
}

//The overlap metric solves S_aux c = <aux|rho> instead.
TEST(FittingIoCoverageIntegratorTests, OverlapMetricFitMatchesDirectSolve)
{
    const std::filesystem::path p = epoxide_fixture();
    if (p.empty()) GTEST_SKIP() << "tests/epoxide_gbw/epoxide.gbw not found";
    WFN wave(p);
    auto basis = combo_basis();
    const WFN aux = generate_aux_wfn(wave, basis);
    Int_Params np(wave), ap(aux);
    vec S, g;
    compute2C<Overlap2C_SPH>(ap, S);
    computeRho<Overlap3C_SPH>(np, ap, wave.get_dm(), g);
    const size_t n = (size_t)ap.get_nao();
    vec S_copy = S, c_ref = g;
    ASSERT_EQ(solve_linear_system(S_copy, n, c_ref), 0);

    DensityFitting::CONFIG cfg;
    cfg.metric = DensityFitting::METRIC_TYPE::OVERLAP;
    StreamCapture out(std::cout);
    const vec c = DensityFitting::density_fit(wave, aux, cfg);
    ASSERT_EQ(c.size(), n);
    double cmax = 0.0;
    for (double v : c_ref) cmax = std::max(cmax, std::fabs(v));
    for (size_t i = 0; i < n; i++)
        EXPECT_NEAR(c[i], c_ref[i], 1e-9 * cmax) << i;
    EXPECT_NE(out.str().find("Metric: Overlap"), std::string::npos);
}

//Tikhonov with lambda = 1 solves (H + I) c = g exactly.
TEST(FittingIoCoverageIntegratorTests, TikhonovShiftsTheDiagonal)
{
    const std::filesystem::path p = epoxide_fixture();
    if (p.empty()) GTEST_SKIP() << "tests/epoxide_gbw/epoxide.gbw not found";
    WFN wave(p);
    auto basis = combo_basis();
    const WFN aux = generate_aux_wfn(wave, basis);
    const CoulombSystem sys = coulomb_system(wave, aux);
    vec H = sys.H;
    for (size_t i = 0; i < sys.n; i++) H[i * sys.n + i] += 1.0;

    DensityFitting::CONFIG cfg;
    cfg.use_tikhonov = true;
    cfg.tikhonov_lambda = 1.0;
    StreamCapture out(std::cout);
    const vec c = DensityFitting::density_fit(wave, aux, cfg);
    ASSERT_EQ(c.size(), sys.n);
    EXPECT_LT(relative_residual(H, c, sys.g), 1e-9);
    vec H_plain = sys.H;
    EXPECT_GT(relative_residual(H_plain, c, sys.g), 1e-6);
    const std::string log = out.str();
    EXPECT_NE(log.find("Fit controls:"), std::string::npos);
    EXPECT_NE(log.find("Tikhonov: on (lambda=1"), std::string::npos);
    EXPECT_NE(log.find("Atomic charge restraints: off"), std::string::npos);
    EXPECT_NE(log.find("Multipole restraints: off"), std::string::npos);
    EXPECT_NE(log.find("Exact total-electron constraint: off"), std::string::npos);
    EXPECT_NE(log.find("Solving density-fitting system..."), std::string::npos);
}

//Atom-centred nuclear-charge restraints with a fixed row scale w = 0.5 add w^2 r_a r_a^T to the
//(Tikhonov-shifted) matrix and w^2 Z_a r_a to the right-hand side, one row per atom.
TEST(FittingIoCoverageIntegratorTests, NuclearChargeRestraintsAddPenaltyRows)
{
    const std::filesystem::path p = epoxide_fixture();
    if (p.empty()) GTEST_SKIP() << "tests/epoxide_gbw/epoxide.gbw not found";
    WFN wave(p);
    auto basis = combo_basis();
    const WFN aux = generate_aux_wfn(wave, basis);
    const CoulombSystem sys = coulomb_system(wave, aux);
    const vec2 rows = population_rows(aux);
    ASSERT_EQ(rows.size(), 7u);
    ASSERT_EQ(rows[0].size(), sys.n);
    vec H = sys.H, g = sys.g;
    for (size_t i = 0; i < sys.n; i++) H[i * sys.n + i] += 1.0;
    for (int a = 0; a < 7; a++)
    {
        const double Z = aux.get_atom_charge(a);
        for (size_t i = 0; i < sys.n; i++)
        {
            g[i] += 0.25 * rows[a][i] * Z;
            for (size_t j = 0; j < sys.n; j++)
                H[i * sys.n + j] += 0.25 * rows[a][i] * rows[a][j];
        }
    }

    DensityFitting::CONFIG cfg;
    cfg.use_tikhonov = true;
    cfg.tikhonov_lambda = 1.0;
    cfg.restrain_charges = true;
    cfg.charge_scheme = DensityFitting::CHARGE_SCHEME::NUCLEAR;
    cfg.adaptive_restraint = false;
    cfg.restraint_strength = 0.5;
    cfg.analyze_quality = true;
    StreamCapture out(std::cout);
    const vec c = DensityFitting::density_fit(wave, aux, cfg);
    ASSERT_EQ(c.size(), sys.n);
    EXPECT_LT(relative_residual(H, c, g), 1e-9);
    const std::string log = out.str();
    EXPECT_NE(log.find("Atomic charge restraints: on (Nuclear Charge)"), std::string::npos);
    EXPECT_NE(log.find("Restraint definition: atom centred"), std::string::npos);
    EXPECT_NE(log.find("Setting charge-restraint row scale to: 0.500000"), std::string::npos);
    EXPECT_NE(log.find("Added charge restraints for 7 atoms."), std::string::npos);
    EXPECT_NE(log.find("=== Density Fitting Quality Analysis ==="), std::string::npos);
    EXPECT_NE(log.find(", Expected = "), std::string::npos);
}

//The exact total-electron constraint pins the fitted electron count to 24.
TEST(FittingIoCoverageIntegratorTests, TotalElectronConstraintIsExact)
{
    const std::filesystem::path p = epoxide_fixture();
    if (p.empty()) GTEST_SKIP() << "tests/epoxide_gbw/epoxide.gbw not found";
    WFN wave(p);
    auto basis = combo_basis();
    const WFN aux = generate_aux_wfn(wave, basis);
    DensityFitting::CONFIG cfg;
    cfg.use_tikhonov = true;
    cfg.tikhonov_lambda = 1.0;
    cfg.constrain_total_electrons = true;
    StreamCapture out(std::cout);
    const vec c = DensityFitting::density_fit(wave, aux, cfg);
    ASSERT_FALSE(c.empty());
    EXPECT_NEAR(fitted_electrons(population_rows(aux), c), 24.0, 1e-6);
    const std::string log = out.str();
    EXPECT_NE(log.find("Exact total-electron constraint: on"), std::string::npos);
    EXPECT_NE(log.find("Exact total-electron constraint: target = 24.0000000000"), std::string::npos);
}

//Atom-centred multipole restraints up to l = 1 print the configuration and the multipole report.
//No Tikhonov term: lambda = 1 on the Coulomb metric damps the diffuse s functions and loses
//about 0.2 e, so only the restraint-only fit keeps the electron count within 0.1.
TEST(FittingIoCoverageIntegratorTests, AtomCentredMultipoleRestraintsReport)
{
    const std::filesystem::path p = epoxide_fixture();
    if (p.empty()) GTEST_SKIP() << "tests/epoxide_gbw/epoxide.gbw not found";
    WFN wave(p);
    auto basis = combo_basis();
    const WFN aux = generate_aux_wfn(wave, basis);
    DensityFitting::CONFIG cfg;
    cfg.restrain_multipoles = true;
    cfg.multipole_lmax = 1;
    cfg.charge_scheme = DensityFitting::CHARGE_SCHEME::HIRSHFELD;
    StreamCapture out(std::cout);
    const vec c = DensityFitting::density_fit(wave, aux, cfg);
    ASSERT_FALSE(c.empty());
    EXPECT_NEAR(fitted_electrons(population_rows(aux), c), 24.0, 0.1);
    const std::string log = out.str();
    EXPECT_NE(log.find("Multipole restraints: on (lmax=1, Hirshfeld)"), std::string::npos);
    EXPECT_NE(log.find("Restraint definition: atom centred"), std::string::npos);
    EXPECT_NE(log.find("Added multipole restraints up to l=1 for 7 atoms."), std::string::npos);
    EXPECT_NE(log.find("Atomic multipoles of the fitted density, Racah normalisation"), std::string::npos);
    EXPECT_NE(log.find("  Atom  l  m       target       fitted    deviation"), std::string::npos);
}

//Grid-partitioned charge restraints (Hirshfeld targets) keep the electron count and are
//labelled as such in the configuration and the report.
TEST(FittingIoCoverageIntegratorTests, GridPartitionedChargeRestraints)
{
    const std::filesystem::path p = epoxide_fixture();
    if (p.empty()) GTEST_SKIP() << "tests/epoxide_gbw/epoxide.gbw not found";
    WFN wave(p);
    auto basis = combo_basis();
    const WFN aux = generate_aux_wfn(wave, basis);
    DensityFitting::CONFIG cfg;
    cfg.restrain_charges = true;
    cfg.partition_restraints = true;
    cfg.charge_scheme = DensityFitting::CHARGE_SCHEME::HIRSHFELD;
    StreamCapture out(std::cout);
    const vec c = DensityFitting::density_fit(wave, aux, cfg);
    ASSERT_FALSE(c.empty());
    EXPECT_NEAR(fitted_electrons(population_rows(aux), c), 24.0, 0.1);
    const std::string log = out.str();
    EXPECT_NE(log.find("Atomic charge restraints: on (Hirshfeld)"), std::string::npos);
    EXPECT_NE(log.find("Restraint definition: grid partitioned"), std::string::npos);
    EXPECT_NE(log.find("Added charge restraints for 7 atoms."), std::string::npos);
    EXPECT_EQ(log.find("Atomic multipoles of the fitted density"), std::string::npos);
}

//The restraint set-up refuses lmax outside 1..8, a non-grid scheme for multipoles or
//partitioned charges, and a negative Tikhonov parameter.
TEST(FittingIoCoverageIntegratorTests, DensityFitRejectsInvalidConfigurations)
{
    const std::filesystem::path p = epoxide_fixture();
    if (p.empty()) GTEST_SKIP() << "tests/epoxide_gbw/epoxide.gbw not found";
    WFN wave(p);
    auto basis = combo_basis();
    const WFN aux = generate_aux_wfn(wave, basis);

    DensityFitting::CONFIG too_high;
    too_high.restrain_multipoles = true;
    too_high.multipole_lmax = 9;
    too_high.charge_scheme = DensityFitting::CHARGE_SCHEME::HIRSHFELD;
    EXPECT_EXIT(DensityFitting::density_fit(wave, aux, too_high), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");

    DensityFitting::CONFIG nuclear_multipoles;
    nuclear_multipoles.restrain_multipoles = true;
    nuclear_multipoles.multipole_lmax = 1;
    nuclear_multipoles.charge_scheme = DensityFitting::CHARGE_SCHEME::NUCLEAR;
    EXPECT_EXIT(DensityFitting::density_fit(wave, aux, nuclear_multipoles), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");

    DensityFitting::CONFIG nuclear_partition;
    nuclear_partition.restrain_charges = true;
    nuclear_partition.partition_restraints = true;
    nuclear_partition.charge_scheme = DensityFitting::CHARGE_SCHEME::NUCLEAR;
    EXPECT_EXIT(DensityFitting::density_fit(wave, aux, nuclear_partition), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");

    DensityFitting::CONFIG negative_lambda;
    negative_lambda.use_tikhonov = true;
    negative_lambda.tikhonov_lambda = -1.0;
    EXPECT_EXIT(DensityFitting::density_fit(wave, aux, negative_lambda), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
}

//A zero coefficient vector fits no electrons: every population is 0.000, every charge is Z,
//the deviation from the nuclear expectation exceeds one electron only for the three heavy atoms
//(the H deviation of exactly 1.000 does not cross the strict > 1 threshold), and the total line
//reads 24.000 expected against 0.000 fitted.
TEST(FittingIoCoverageIntegratorTests, AnalyzeQualityReportsPerAtomPopulations)
{
    const std::filesystem::path p = epoxide_fixture();
    if (p.empty()) GTEST_SKIP() << "tests/epoxide_gbw/epoxide.gbw not found";
    WFN wave(p);
    auto basis = combo_basis();
    const WFN aux = generate_aux_wfn(wave, basis);
    const aux_density_table table(aux.get_atoms());
    const vec zero(table.n_coef, 0.0);
    vec expected(7);
    for (int a = 0; a < 7; a++) expected[a] = aux.get_atom_charge(a);
    StreamCapture out(std::cout);
    DensityFitting::analyze_density_fit_quality(zero, aux, table, expected);
    const std::string log = out.str();
    EXPECT_NE(log.find("=== Density Fitting Quality Analysis ==="), std::string::npos);
    size_t warnings = 0, pos = 0;
    while ((pos = log.find("WARNING: significant deviation", pos)) != std::string::npos) { warnings++; pos++; }
    EXPECT_EQ(warnings, 3u) << log;
    EXPECT_NE(log.find("Population = 0.000, Charge = "), std::string::npos);
    EXPECT_NE(log.find(", Expected = "), std::string::npos);
    EXPECT_NE(log.find(", Deviation = "), std::string::npos);
    EXPECT_NE(log.find("Expected / Real total electrons: 24.000 / 0.000"), std::string::npos) << log;
}

//The demonstration runs the unrestrained, enhanced and hybrid fits and prints the comparison.
TEST(FittingIoCoverageIntegratorTests, DemonstrateEnhancedFittingPrintsComparison)
{
    const std::filesystem::path p = epoxide_fixture();
    if (p.empty()) GTEST_SKIP() << "tests/epoxide_gbw/epoxide.gbw not found";
    WFN wave(p);
    auto basis = combo_basis();
    const WFN aux = generate_aux_wfn(wave, basis);
    StreamCapture out(std::cout);
    DensityFitting::demonstrate_enhanced_density_fitting(wave, aux);
    const std::string log = out.str();
    EXPECT_NE(log.find("=== Enhanced Density Fitting Demonstration ==="), std::string::npos);
    EXPECT_NE(log.find("--- Method 0: Unrestrained (Baseline) ---"), std::string::npos);
    EXPECT_NE(log.find("--- Method 1: Enhanced Adaptive Restraints ---"), std::string::npos);
    EXPECT_NE(log.find("--- Method 2: Hybrid Regularization ---"), std::string::npos);
    EXPECT_NE(log.find("Time for unrestrained fit: "), std::string::npos);
    EXPECT_NE(log.find("(0) | Unrestrained: RRS = "), std::string::npos);
    EXPECT_NE(log.find("=== Method Comparison ==="), std::string::npos);
    EXPECT_NE(log.find("Unrestrained        - Max coeff: "), std::string::npos);
    EXPECT_NE(log.find("=== RMS Coefficient Differences ==="), std::string::npos);
    EXPECT_NE(log.find("Unrestrained vs Enhanced: "), std::string::npos);
    EXPECT_NE(log.find("=== Average Coefficient Magnitudes ==="), std::string::npos);
}

//QM_RI_difference_cube writes the WFN, RI and difference cubes next to the wavefunction on a
//0.1 A grid with a 3 A margin; the difference cube is WFN minus RI point by point.
TEST(FittingIoCoverageIntegratorTests, QmRiDifferenceCubeWritesThreeCubes)
{
    const std::filesystem::path p = epoxide_fixture();
    if (p.empty()) GTEST_SKIP() << "tests/epoxide_gbw/epoxide.gbw not found";
    Scratch s("QmRiDifferenceCube");
    const std::filesystem::path local = s.dir / "epoxide.gbw";
    std::filesystem::copy_file(p, local);
    WFN wave(local);
    auto basis = combo_basis();
    const WFN aux = generate_aux_wfn(wave, basis);

    std::array<int, 3> steps{};
    for (int i = 0; i < 3; i++)
    {
        double lo = wave.get_atom_coordinate(0, i), hi = lo;
        for (int a = 1; a < wave.get_ncen(); a++)
        {
            lo = std::min(lo, wave.get_atom_coordinate(a, i));
            hi = std::max(hi, wave.get_atom_coordinate(a, i));
        }
        steps[i] = (int)std::ceil(constants::bohr2ang(hi - lo + 2.0 * constants::ang2bohr(3.0)) / 0.1);
    }

    StreamCapture out(std::cout);
    DensityFitting::QM_RI_difference_cube(wave, aux);
    const std::string log = out.str();
    EXPECT_NE(log.find("Calculating WFN density cube..."), std::string::npos);
    EXPECT_NE(log.find("Calculating RI density cube..."), std::string::npos);
    EXPECT_NE(log.find("RRS = "), std::string::npos);
    EXPECT_NE(log.find(", |Sum| = "), std::string::npos);

    const std::filesystem::path wfn_cube = s.dir / "epoxide_rho_WFN.cube";
    const std::filesystem::path ri_cube = s.dir / "epoxide_rho_RI.cube";
    const std::filesystem::path diff_cube = s.dir / "epoxide_rho_diff_RI.cube";
    ASSERT_TRUE(std::filesystem::exists(wfn_cube));
    ASSERT_TRUE(std::filesystem::exists(ri_cube));
    ASSERT_TRUE(std::filesystem::exists(diff_cube));

    std::ostringstream log_w, log_r, log_d;
    WFN empty_w(e_origin::NOT_YET_DEFINED), empty_r(e_origin::NOT_YET_DEFINED), empty_d(e_origin::NOT_YET_DEFINED);
    cube w(wfn_cube, true, empty_w, log_w);
    cube r(ri_cube, true, empty_r, log_r);
    cube d(diff_cube, true, empty_d, log_d);
    for (int i = 0; i < 3; i++)
    {
        EXPECT_EQ(w.get_size(i), steps[i]) << i;
        EXPECT_EQ(r.get_size(i), steps[i]) << i;
        EXPECT_EQ(d.get_size(i), steps[i]) << i;
    }
    EXPECT_EQ(empty_w.get_ncen(), 7);
    const int cx = steps[0] / 2, cy = steps[1] / 2, cz = steps[2] / 2;
    for (int dx = -1; dx <= 1; dx++)
    {
        const double vw = w.get_value(cx + dx, cy, cz), vr = r.get_value(cx + dx, cy, cz), vd = d.get_value(cx + dx, cy, cz);
        EXPECT_NEAR(vd, vw - vr, 1e-4 * (std::fabs(vw) + std::fabs(vr)) + 1e-5) << dx;
    }
    EXPECT_GT(w.get_value(cx, cy, cz) + w.get_value(cx - 1, cy, cz) + w.get_value(cx + 1, cy, cz), 0.0);
}

//----------------------------------------------------------------------------------------------
//basis_set.cpp
//----------------------------------------------------------------------------------------------

//Sorted descending, 1.0 is dropped against 1.05 (relative gap 0.048 < 0.1), the rest survive.
TEST(FittingIoCoverageBasisTests, DedupExponentsDropsNearDuplicates)
{
    const std::vector<double> out = dedup_exponents({ 1.0, 1.05, 2.0, 0.5 }, 0.1);
    ASSERT_EQ(out.size(), 3u);
    EXPECT_DOUBLE_EQ(out[0], 2.0);
    EXPECT_DOUBLE_EQ(out[1], 1.05);
    EXPECT_DOUBLE_EQ(out[2], 0.5);
    EXPECT_TRUE(dedup_exponents({}, 0.1).empty());
    const std::vector<double> single = dedup_exponents({ 3.0 }, 0.1);
    ASSERT_EQ(single.size(), 1u);
    EXPECT_DOUBLE_EQ(single[0], 3.0);
}

//Diagonal {1, 4, 0.01}: the pivot order is the largest diagonal first and 0.01 falls under the
//threshold. For [[1, 0.9], [0.9, 1]] the second residual diagonal is 1 - 0.81 = 0.19, kept
//at threshold 0.1 and dropped at 0.5.
TEST(FittingIoCoverageBasisTests, PivotedCholeskyOrdersAndStops)
{
    dMatrix2 D = reshape<dMatrix2>(vec{ 1.0, 0, 0, 0, 4.0, 0, 0, 0, 0.01 }, Shape2D(3, 3));
    const std::vector<int> piv = pivoted_cholesky(D, 0.1);
    ASSERT_EQ(piv.size(), 2u);
    EXPECT_EQ(piv[0], 1);
    EXPECT_EQ(piv[1], 0);

    dMatrix2 C = reshape<dMatrix2>(vec{ 1.0, 0.9, 0.9, 1.0 }, Shape2D(2, 2));
    const std::vector<int> one = pivoted_cholesky(C, 0.5);
    ASSERT_EQ(one.size(), 1u);
    EXPECT_EQ(one[0], 0);
    dMatrix2 C2 = reshape<dMatrix2>(vec{ 1.0, 0.9, 0.9, 1.0 }, Shape2D(2, 2));
    const std::vector<int> two = pivoted_cholesky(C2, 0.1);
    ASSERT_EQ(two.size(), 2u);
    EXPECT_EQ(two[0], 0);
    EXPECT_EQ(two[1], 1);
}

//Two same-centre normalised Gaussians of exponents a, b have the Coulomb overlap
//(2 sqrt(ab) / (a + b))^(l + 1/2): for 10 and 0.1 that is 0.445 (s) and 0.088 (p), so the
//residual diagonal after the first pivot is 0.802 (s) and 0.992 (p). The threshold decides
//whether the second exponent survives; the result is sorted descending.
TEST(FittingIoCoverageBasisTests, PruneCandidatesKeepsWellConditionedExponents)
{
    EXPECT_TRUE(prune_element_candidates_for_L(0, {}, 0.5).empty());
    const std::vector<double> single = prune_element_candidates_for_L(0, { 2.0 }, 0.5);
    ASSERT_EQ(single.size(), 1u);
    EXPECT_DOUBLE_EQ(single[0], 2.0);
    const std::vector<double> near_pair = prune_element_candidates_for_L(0, { 1.0, 1.05 }, 0.5);
    ASSERT_EQ(near_pair.size(), 1u);
    EXPECT_DOUBLE_EQ(near_pair[0], 1.05);

    StreamCapture out(std::cout);
    const std::vector<double> s_both = prune_element_candidates_for_L(0, { 0.1, 10.0 }, 0.5);
    ASSERT_EQ(s_both.size(), 2u);
    EXPECT_DOUBLE_EQ(s_both[0], 10.0);
    EXPECT_DOUBLE_EQ(s_both[1], 0.1);
    EXPECT_EQ(prune_element_candidates_for_L(0, { 0.1, 10.0 }, 0.9).size(), 1u);
    const std::vector<double> p_both = prune_element_candidates_for_L(1, { 0.1, 10.0 }, 0.5);
    ASSERT_EQ(p_both.size(), 2u);
    EXPECT_DOUBLE_EQ(p_both[0], 10.0);
    EXPECT_DOUBLE_EQ(p_both[1], 0.1);
    EXPECT_EQ(prune_element_candidates_for_L(1, { 0.1, 10.0 }, 0.999).size(), 1u);
}

//One hydrogen with a d, f, g and h shell: the coefficient count is 5 + 7 + 9 + 11 = 32 and the
//complete pass emits one primitive per spherical function.
TEST(FittingIoCoverageBasisTests, LoadBasisCompleteEmitsSphericalCountForHighShells)
{
    auto bs = std::make_shared<BasisSet>();
    bs->set_count_for_element(0, 4);
    bs->add_owned_primitive({ 0, 2, 2.0, 1.0, 0 });
    bs->add_owned_primitive({ 0, 3, 1.5, 1.0, 1 });
    bs->add_owned_primitive({ 0, 4, 1.1, 1.0, 2 });
    bs->add_owned_primitive({ 0, 5, 0.7, 1.0, 3 });
    WFN w(e_origin::NOT_YET_DEFINED);
    w.push_back_atom("H", 0.0, 0.0, 0.0, 1);
    EXPECT_EQ(load_basis_into_WFN(w, bs, true, true), 32);
    EXPECT_EQ(w.get_nex(), 32);
    ASSERT_EQ(w.get_atom_basis_set_size(0), 4);
    for (int k = 0; k < 4; k++)
    {
        EXPECT_EQ((int)w.get_atom_basis_set_entry(0, k).get_type(), k + 2) << k;
        EXPECT_EQ(w.get_atom_basis_set_entry(0, k).get_shell(), k) << k;
        EXPECT_DOUBLE_EQ(w.get_atom_basis_set_entry(0, k).get_coefficient(), 1.0) << k;
    }
    for (int i = 0; i < 32; i++)
        EXPECT_EQ(w.get_center(i), 1) << i;
    EXPECT_NEAR(w.get_exponent(0), 2.0, 1e-15);
    EXPECT_NEAR(w.get_exponent(31), 0.7, 1e-15);
}

//The complete pass numbers the primitives in the Cartesian wfn scheme (d = 5..10, f = 11..20,
//g = 21..35, h = 36..56), which has 6/10/15/21 functions per shell, so a d+f+g+h hydrogen
//should carry 52 primitives with the d block being exactly types 5..10.
//suspected defect: Src/core/basis_set.cpp load_basis_into_WFN complete pass emits 2l+1 entries per shell but numbers them in the Cartesian wfn scheme (5..10 for d), so d and higher shells are truncated
TEST(FittingIoCoverageBasisTests, DISABLED_LoadBasisCompleteUsesCartesianNumbering)
{
    auto bs = std::make_shared<BasisSet>();
    bs->set_count_for_element(0, 4);
    bs->add_owned_primitive({ 0, 2, 2.0, 1.0, 0 });
    bs->add_owned_primitive({ 0, 3, 1.5, 1.0, 1 });
    bs->add_owned_primitive({ 0, 4, 1.1, 1.0, 2 });
    bs->add_owned_primitive({ 0, 5, 0.7, 1.0, 3 });
    WFN w(e_origin::NOT_YET_DEFINED);
    w.push_back_atom("H", 0.0, 0.0, 0.0, 1);
    load_basis_into_WFN(w, bs, true, true);
    ASSERT_EQ(w.get_nex(), 52);
    for (int i = 0; i < 6; i++)
        EXPECT_EQ(w.get_type(i), 5 + i) << i;
    EXPECT_EQ(w.get_type(6), 11);
    EXPECT_EQ(w.get_type(16), 21);
    EXPECT_EQ(w.get_type(31), 36);
    EXPECT_EQ(w.get_type(51), 56);
}

//Lithium with s (2.0, 0.5) and p (0.8): l_occ = 1 so the aux set runs to l = 2 with beta 1.8
//throughout. Ranges: l0 a_min 1.0 / a_max 4.0 -> 4 functions, l1 1.3 / 2.8 -> 3, l2 1.6 / 1.6
//-> 1, all as a_min * 1.8^i from the top down.
TEST(FittingIoCoverageBasisTests, AutoAuxLithiumWithSpShells)
{
    atom li("Li", atomID(), 1, 0.0, 0.0, 0.0, 3);
    ASSERT_TRUE(li.push_back_basis_set(2.0, 1.0, 1, 0));
    ASSERT_TRUE(li.push_back_basis_set(0.5, 1.0, 1, 1));
    ASSERT_TRUE(li.push_back_basis_set(0.8, 1.0, 2, 2));
    BasisSet bs;
    StreamCapture out(std::cout);
    bs.gen_auto_aux_for_element(li);
    EXPECT_TRUE(bs.has_element(3));
    ASSERT_EQ(bs.get_owned_primitive_count(), 8u);
    const auto span = bs[2];
    ASSERT_EQ(span.size(), 8u);
    const double b = 1.8;
    const std::vector<double> exps{ 1.0 * b * b * b, 1.0 * b * b, 1.0 * b, 1.0, 1.3 * b * b, 1.3 * b, 1.3, 1.6 };
    const std::vector<int> types{ 0, 0, 0, 0, 1, 1, 1, 2 };
    for (size_t i = 0; i < 8; i++)
    {
        EXPECT_NEAR(span[i].exp, exps[i], 1e-9) << i;
        EXPECT_EQ(span[i].type, types[i]) << i;
        EXPECT_EQ(span[i].shell, (int)i) << i;
        EXPECT_EQ(span[i].center, 0) << i;
        EXPECT_DOUBLE_EQ(span[i].coefficient, 1.0) << i;
    }
    const std::string log = out.str();
    EXPECT_NE(log.find("Beta for Element 3 and l 0 : 1.8 with 4 Functions."), std::string::npos) << log;
    EXPECT_NE(log.find("Beta for Element 3 and l 1 : 1.8 with 3 Functions."), std::string::npos) << log;
    EXPECT_NE(log.find("Beta for Element 3 and l 2 : 1.8 with 1 Functions."), std::string::npos) << log;
}

//A contracted s shell whose second primitive carries a 1e-5 coefficient has an effective
//exponent of ~1e10 for that primitive, so the aux-derived maximum exceeds 1e7 for l = 0 and
//l = 1. l = 0 is capped by the primitive maximum (4.0, beta 1.8, 4 functions); l = 1 takes the
//>1e7 guard and falls back to the primitive maximum 2.8 with beta 2.0 (3 functions: 5.2, 2.6,
//1.3); l = 2 keeps its aux value 2.2635 with beta 2.2 (2 functions: 3.52, 1.6).
TEST(FittingIoCoverageBasisTests, AutoAuxTinyContractionTriggersLargeExponentGuard)
{
    atom h("H", atomID(), 1, 0.0, 0.0, 0.0, 1);
    ASSERT_TRUE(h.push_back_basis_set(2.0, 1.0, 1, 0));
    ASSERT_TRUE(h.push_back_basis_set(0.5, 1e-5, 1, 0));
    ASSERT_TRUE(h.push_back_basis_set(0.8, 1.0, 2, 1));
    BasisSet bs;
    StreamCapture out(std::cout);
    bs.gen_auto_aux_for_element(h);
    ASSERT_EQ(bs.get_owned_primitive_count(), 9u);
    const auto span = bs[0];
    ASSERT_EQ(span.size(), 9u);
    const std::vector<double> exps{ 5.832, 3.24, 1.8, 1.0, 5.2, 2.6, 1.3, 3.52, 1.6 };
    const std::vector<int> types{ 0, 0, 0, 0, 1, 1, 1, 2, 2 };
    for (size_t i = 0; i < 9; i++)
    {
        EXPECT_NEAR(span[i].exp, exps[i], 1e-9) << i;
        EXPECT_EQ(span[i].type, types[i]) << i;
    }
    const std::string log = out.str();
    EXPECT_NE(log.find("Beta for Element 1 and l 0 : 1.8 with 4 Functions."), std::string::npos) << log;
    EXPECT_NE(log.find("Beta for Element 1 and l 1 : 2 with 3 Functions."), std::string::npos) << log;
    EXPECT_NE(log.find("Beta for Element 1 and l 2 : 2.2 with 2 Functions."), std::string::npos) << log;
}

//Potassium (l_occ = 2) with the same s, s, p shells: l_max_aux = min(4, 2) = 2, so the aux set
//is the same eight functions as for lithium.
//suspected defect: Src/core/basis_set.cpp gen_auto_aux_for_element a_max_adjusted loop runs to 2*l_occ_max beyond size l_max_aux+1
TEST(FittingIoCoverageBasisTests, DISABLED_AutoAuxPotassiumKeepsAuxWithinLmax)
{
    atom k("K", atomID(), 1, 0.0, 0.0, 0.0, 19);
    ASSERT_TRUE(k.push_back_basis_set(2.0, 1.0, 1, 0));
    ASSERT_TRUE(k.push_back_basis_set(0.5, 1.0, 1, 1));
    ASSERT_TRUE(k.push_back_basis_set(0.8, 1.0, 2, 2));
    BasisSet bs;
    StreamCapture out(std::cout);
    bs.gen_auto_aux_for_element(k);
    ASSERT_EQ(bs.get_owned_primitive_count(), 8u);
    const auto span = bs[18];
    const double b = 1.8;
    const std::vector<double> exps{ 1.0 * b * b * b, 1.0 * b * b, 1.0 * b, 1.0, 1.3 * b * b, 1.3 * b, 1.3, 1.6 };
    for (size_t i = 0; i < 8; i++)
        EXPECT_NEAR(span[i].exp, exps[i], 1e-9) << i;
}

//An s-only lithium (exponent 2.0): l_max = 0, l_max_aux = 0, one s aux function of exponent 4.0.
//suspected defect: Src/core/basis_set.cpp gen_auto_aux_for_element a_max_adjusted loop runs to 2*l_occ_max beyond size l_max_aux+1
TEST(FittingIoCoverageBasisTests, DISABLED_AutoAuxSOnlyLithium)
{
    atom li("Li", atomID(), 1, 0.0, 0.0, 0.0, 3);
    ASSERT_TRUE(li.push_back_basis_set(2.0, 1.0, 1, 0));
    BasisSet bs;
    StreamCapture out(std::cout);
    bs.gen_auto_aux_for_element(li);
    ASSERT_EQ(bs.get_owned_primitive_count(), 1u);
    EXPECT_NEAR(bs[2][0].exp, 4.0, 1e-9);
    EXPECT_EQ(bs[2][0].type, 0);
}

//Turbomole layout with one shell of each supported type on a bare hydrogen, read with the
//debug trace on: s/p/d/f become types 1..4 in shells 0..3.
TEST(FittingIoCoverageBasisTests, ReadMissingTurbomoleAllShellTypes)
{
    Scratch s("ReadMissingTurbomole");
    write_text(s.dir / "all_shells.basis", "keys= { turbomole= }\nH:\n{\n1 s\n2.0 1.0\n1 p\n0.8 1.0\n1 d\n1.1 1.0\n1 f\n0.7 1.0\n}\n");
    WFN w(e_origin::NOT_YET_DEFINED);
    w.push_back_atom("H", 0.0, 0.0, 0.0, 1);
    w.set_basis_set_name("all_shells.basis");
    StreamCapture out(std::cout);
    ASSERT_TRUE(BasisSetLibrary::read_basis_set_missing(s.dir, w, true));
    ASSERT_EQ(w.get_atom_basis_set_size(0), 4);
    const std::vector<double> exps{ 2.0, 0.8, 1.1, 0.7 };
    for (int k = 0; k < 4; k++)
    {
        EXPECT_EQ((int)w.get_atom_basis_set_entry(0, k).get_type(), k + 1) << k;
        EXPECT_EQ(w.get_atom_basis_set_entry(0, k).get_shell(), k) << k;
        EXPECT_NEAR(w.get_atom_basis_set_entry(0, k).get_exponent(), exps[k], 1e-15) << k;
        EXPECT_NEAR(w.get_atom_basis_set_entry(0, k).get_coefficient(), 1.0, 1e-15) << k;
    }
    EXPECT_EQ(w.get_atom_shell_count(0), 4);
    const std::string log = out.str();
    EXPECT_NE(log.find("basis set is valid, continueing..."), std::string::npos) << log;
    EXPECT_NE(log.find("File of basis set to load: "), std::string::npos);
    EXPECT_NE(log.find("Found keys=!"), std::string::npos);
    EXPECT_NE(log.find("This file is written in turbomole type!"), std::string::npos);
    EXPECT_NE(log.find("I read an additional line!"), std::string::npos);
    EXPECT_NE(log.find("It's a match!"), std::string::npos);
    EXPECT_NE(log.find("I found }!"), std::string::npos);
    EXPECT_NE(log.find("FINISHED WITH READING BASIS SET!"), std::string::npos);
}

//Gamess-US layout: "type count" headers and an index column in front of exponent and coefficient.
TEST(FittingIoCoverageBasisTests, ReadMissingGamessLayout)
{
    Scratch s("ReadMissingGamess");
    write_text(s.dir / "gamess.basis", "keys= { gamess-us= }\nH:\n{\ns 1\n1 2.0 1.0\np 1\n1 0.8 1.0\n}\n");
    WFN w(e_origin::NOT_YET_DEFINED);
    w.push_back_atom("H", 0.0, 0.0, 0.0, 1);
    w.set_basis_set_name("gamess.basis");
    ASSERT_TRUE(BasisSetLibrary::read_basis_set_missing(s.dir, w, false));
    ASSERT_EQ(w.get_atom_basis_set_size(0), 2);
    EXPECT_EQ((int)w.get_atom_basis_set_entry(0, 0).get_type(), 1);
    EXPECT_NEAR(w.get_atom_basis_set_entry(0, 0).get_exponent(), 2.0, 1e-15);
    EXPECT_EQ((int)w.get_atom_basis_set_entry(0, 1).get_type(), 2);
    EXPECT_EQ(w.get_atom_basis_set_entry(0, 1).get_shell(), 1);
    EXPECT_NEAR(w.get_atom_basis_set_entry(0, 1).get_exponent(), 0.8, 1e-15);
}

//Two bare elements in one file, lithium listed after hydrogen but requested first: the reader
//rewinds for every element and only the atoms still without a basis receive the block.
TEST(FittingIoCoverageBasisTests, ReadMissingTwoElementsRewind)
{
    Scratch s("ReadMissingTwoElements");
    write_text(s.dir / "two.basis", "keys= { turbomole= }\nH:\n{\n1 s\n2.0 1.0\n}\nLi:\n{\n1 s\n1.5 1.0\n1 p\n0.9 1.0\n}\n");
    WFN w(e_origin::NOT_YET_DEFINED);
    w.push_back_atom("Li", 0.0, 0.0, 0.0, 3);
    w.push_back_atom("H", 0.0, 0.0, 1.5, 1);
    w.push_back_atom("H", 0.0, 0.0, -1.5, 1);
    w.set_basis_set_name("two.basis");
    ASSERT_TRUE(BasisSetLibrary::read_basis_set_missing(s.dir, w, false));
    EXPECT_EQ(w.get_nr_basis_set_loaded(), 3);
    ASSERT_EQ(w.get_atom_basis_set_size(0), 2);
    EXPECT_NEAR(w.get_atom_basis_set_entry(0, 0).get_exponent(), 1.5, 1e-15);
    EXPECT_EQ((int)w.get_atom_basis_set_entry(0, 1).get_type(), 2);
    EXPECT_NEAR(w.get_atom_basis_set_entry(0, 1).get_exponent(), 0.9, 1e-15);
    for (int a = 1; a < 3; a++)
    {
        ASSERT_EQ(w.get_atom_basis_set_size(a), 1) << a;
        EXPECT_NEAR(w.get_atom_basis_set_entry(a, 0).get_exponent(), 2.0, 1e-15) << a;
        EXPECT_EQ((int)w.get_atom_basis_set_entry(a, 0).get_type(), 1) << a;
    }
}

//Every rejection path returns false: no file, an unknown format key, a g shell, an element
//that is not in the file, and a gamess index that exceeds the shell's primitive count.
TEST(FittingIoCoverageBasisTests, ReadMissingRejections)
{
    Scratch s("ReadMissingRejections");
    StreamCapture out(std::cout);
    WFN w(e_origin::NOT_YET_DEFINED);
    w.push_back_atom("H", 0.0, 0.0, 0.0, 1);
    w.set_basis_set_name("absent.basis");
    EXPECT_FALSE(BasisSetLibrary::read_basis_set_missing(s.dir, w, false));
    EXPECT_NE(out.str().find("could not find this basis set"), std::string::npos) << out.str();

    write_text(s.dir / "nwchem.basis", "keys= { nwchem= }\nH:\n{\n1 s\n2.0 1.0\n}\n");
    WFN w2(e_origin::NOT_YET_DEFINED);
    w2.push_back_atom("H", 0.0, 0.0, 0.0, 1);
    w2.set_basis_set_name("nwchem.basis");
    EXPECT_FALSE(BasisSetLibrary::read_basis_set_missing(s.dir, w2, false));

    write_text(s.dir / "gshell.basis", "keys= { turbomole= }\nH:\n{\n1 g\n2.0 1.0\n}\n");
    WFN w3(e_origin::NOT_YET_DEFINED);
    w3.push_back_atom("H", 0.0, 0.0, 0.0, 1);
    w3.set_basis_set_name("gshell.basis");
    EXPECT_FALSE(BasisSetLibrary::read_basis_set_missing(s.dir, w3, false));

    write_text(s.dir / "honly.basis", "keys= { turbomole= }\nH:\n{\n1 s\n2.0 1.0\n}\n");
    WFN w4(e_origin::NOT_YET_DEFINED);
    w4.push_back_atom("He", 0.0, 0.0, 0.0, 2);
    w4.set_basis_set_name("honly.basis");
    EXPECT_FALSE(BasisSetLibrary::read_basis_set_missing(s.dir, w4, false));

    write_text(s.dir / "badindex.basis", "keys= { gamess-us= }\nH:\n{\ns 1\n5 2.0 1.0\n}\n");
    WFN w5(e_origin::NOT_YET_DEFINED);
    w5.push_back_atom("H", 0.0, 0.0, 0.0, 1);
    w5.set_basis_set_name("badindex.basis");
    EXPECT_FALSE(BasisSetLibrary::read_basis_set_missing(s.dir, w5, false));

    write_text(s.dir / "nokeys.basis", "H:\n{\n1 s\n2.0 1.0\n}\n");
    WFN w6(e_origin::NOT_YET_DEFINED);
    w6.push_back_atom("H", 0.0, 0.0, 0.0, 1);
    w6.set_basis_set_name("nokeys.basis");
    EXPECT_FALSE(BasisSetLibrary::read_basis_set_missing(s.dir, w6, false));
}

//The vanilla reader with the debug trace on the def2-TZVP fixture: same six hydrogen
//primitives in four shells as the sibling test reads silently, plus the trace lines.
TEST(FittingIoCoverageBasisTests, ReadVanillaDebugTrace)
{
    const std::filesystem::path dir = nos_test_repo_root() / "tests" / "NiP3_fchk";
    if (!std::filesystem::exists(dir / "def2-TZVP")) GTEST_SKIP() << "fixture tests/NiP3_fchk/def2-TZVP not found";
    WFN w(e_origin::NOT_YET_DEFINED);
    w.push_back_atom("H", 0.0, 0.0, 0.0, 1);
    w.set_basis_set_name("def2-TZVP");
    StreamCapture out(std::cout);
    ASSERT_TRUE(BasisSetLibrary::read_basis_set_vanilla(dir, w, true));
    EXPECT_EQ(w.get_atom_basis_set_size(0), 6);
    EXPECT_EQ(w.get_atom(0).get_shellcount_size(), 4u);
    const std::string log = out.str();
    EXPECT_NE(log.find("basis set is valid, continuing..."), std::string::npos) << log;
    EXPECT_NE(log.find("File of basis set to load: "), std::string::npos);
    EXPECT_NE(log.find("Found keys=!"), std::string::npos);
    EXPECT_NE(log.find("This file is written in turbomole type!"), std::string::npos);
    EXPECT_NE(log.find("I found }: "), std::string::npos);
    EXPECT_NE(log.find("FINISHED WITH READING BASIS SET!"), std::string::npos);
}

//----------------------------------------------------------------------------------------------
//fchk.cpp
//----------------------------------------------------------------------------------------------

//H2+ doublet (charge 1, multiplicity 2): the debug run reports "alpha, beta, elcount: 1 0 1"
//in width-5 columns, drops the four .debug dumps into the cwd (the spin density one only
//because the wavefunction is unrestricted) and still writes the fchk.
TEST(FittingIoCoverageFchkTests, FreeFchkDebugWritesTraceFiles)
{
    Scratch s("FreeFchkDebug");
    std::filesystem::current_path(s.dir);
    WFN w = occ_h2_plus();
    w.set_origin(e_origin::wfn);
    w.set_method("rhf");
    ASSERT_TRUE(w.get_is_unrestricted());
    ASSERT_EQ(w.get_charge(), 1);
    const std::filesystem::path tmp = s.dir / "h2p.fchk";
    std::ostringstream log;
    StreamCapture out(std::cout);
    ASSERT_TRUE(free_fchk(log, tmp, "", w, true, true));
    EXPECT_TRUE(std::filesystem::exists(tmp));
    EXPECT_TRUE(std::filesystem::exists(s.dir / "norm_prim.debug"));
    EXPECT_TRUE(std::filesystem::exists(s.dir / "cmo.debug"));
    EXPECT_TRUE(std::filesystem::exists(s.dir / "dm.debug"));
    EXPECT_TRUE(std::filesystem::exists(s.dir / "sdm.debug"));
    const std::string text = log.str();
    EXPECT_NE(text.find("alpha, beta, elcount:     1    0    1"), std::string::npos) << text;
    EXPECT_NE(text.find("Origin: "), std::string::npos);
    EXPECT_NE(text.find("DM is in dm.debug"), std::string::npos);
    EXPECT_NE(text.find("SDM is in sdm.debug"), std::string::npos);
    EXPECT_NE(text.find("Starting to write fchk now!"), std::string::npos);
    EXPECT_NE(text.find("Writing "), std::string::npos);
    EXPECT_NE(out.str().find("Done with DM!"), std::string::npos) << out.str();
    std::ifstream in(tmp);
    EXPECT_NEAR(read_fchk_double(in, "Number of electrons"), 1.0, 1e-12);
    EXPECT_NEAR(read_fchk_double(in, "Charge"), 1.0, 1e-12);
    EXPECT_NEAR(read_fchk_double(in, "Multiplicity"), 2.0, 1e-12);
}

//With the basis of the second atom cleared the writer tries to reload the missing atoms from
//the basis-set directory; there is no such file, so err_checkf exits.
TEST(FittingIoCoverageFchkTests, FreeFchkExitsWhenMissingBasisCannotBeRead)
{
    Scratch s("FreeFchkMissingBasis");
    WFN w = occ_h2_plus();
    w.set_origin(e_origin::wfn);
    w.set_method("rhf");
    w.clear_atom_basis_set(1);
    w.set_basis_set_name("no_such_basis");
    ASSERT_EQ(w.get_nr_basis_set_loaded(), 1);
    std::ostringstream log;
    EXPECT_EXIT(free_fchk(log, s.dir / "h2p.fchk", s.dir, w, false, true), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
}

//The stream overload finds the heading, rewinds first, and parses the value from column 49;
//a heading that is absent, or a string shorter than 50 characters, ends the program.
TEST(FittingIoCoverageFchkTests, ReadFchkDoubleFromStreamAndErrors)
{
    Scratch s("ReadFchkDouble");
    const std::string energy = "Total Energy" + std::string(31, ' ') + "R     -1.0500000000E+00";
    ASSERT_EQ(energy.size(), 66u);
    write_text(s.dir / "small.fchk", energy + "\nVirial Ratio" + std::string(31, ' ') + "R      2.0000000000E+00\n");
    {
        std::ifstream in(s.dir / "small.fchk");
        ASSERT_TRUE(in.good());
        EXPECT_NEAR(read_fchk_double(in, "Virial Ratio"), 2.0, 1e-12);
        EXPECT_NEAR(read_fchk_double(in, "Total Energy"), -1.05, 1e-12);
    }
    EXPECT_NEAR(read_fchk_double(energy), -1.05, 1e-12);
    //The death-test child re-runs this body on Windows; its Scratch cannot remove a file the
    //parent still holds open, so no stream stays open across the EXPECT_EXITs
    EXPECT_EXIT({ std::ifstream in(s.dir / "small.fchk"); read_fchk_double(in, "Nonexistent Key"); }, ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
    EXPECT_EXIT(read_fchk_double(std::string("Total Energy R 1.0")), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
}

//----------------------------------------------------------------------------------------------
//npy.h
//----------------------------------------------------------------------------------------------

//The version-1 header of a 2x3 unsigned-short array is the dict padded with spaces to a
//newline-terminated length that makes magic + version + length + header a multiple of 16.
TEST(FittingIoCoverageNpyTests, ReadHeaderReturnsPaddedDict)
{
    npy::npy_data<unsigned short> d;
    d.data = { 1, 2, 3, 4, 5, 6 };
    d.shape = { 2, 3 };
    d.fortran_order = false;
    std::stringstream stream;
    npy::write_npy<unsigned short>(stream, d);
    const std::string header = npy::read_header(stream);
    const std::string dict = std::string("{'descr': '") + npy::host_endian_char + "u2', 'fortran_order': False, 'shape': (2, 3), }";
    ASSERT_GE(header.size(), dict.size());
    EXPECT_EQ(header.substr(0, dict.size()), dict) << header;
    EXPECT_EQ(header.back(), '\n');
    EXPECT_EQ((10 + header.size()) % 16, 0u);
    for (size_t i = dict.size(); i + 1 < header.size(); i++)
        EXPECT_EQ(header[i], ' ') << i;
    std::stringstream again;
    npy::write_npy<unsigned short>(again, d);
    const npy::npy_data<unsigned short> back = npy::read_npy<unsigned short>(again);
    ASSERT_EQ(back.data.size(), 6u);
    EXPECT_EQ(back.data[5], 6);
    ASSERT_EQ(back.shape.size(), 2u);
    EXPECT_EQ(back.shape[1], 3u);
}

//The four-argument loader appends to the data vector it is handed and reports the Fortran flag
//and the shape; a shape-less (0-d) array holds exactly one value.
TEST(FittingIoCoverageNpyTests, LoadArrayFromNumpyLongAppendsWithFortranFlag)
{
    Scratch s("LoadArrayFromNumpyLong");
    npy::npy_data<long> d;
    d.data = { -3, 4, 5, 6, 7, 8 };
    d.shape = { 2, 3 };
    d.fortran_order = true;
    const std::filesystem::path p = s.dir / "long.npy";
    npy::write_npy<long>(p.string(), d);
    std::vector<unsigned long> shape;
    bool fortran = false;
    std::vector<long> data{ 99 };
    npy::LoadArrayFromNumpy<long>(p, shape, fortran, data);
    EXPECT_TRUE(fortran);
    ASSERT_EQ(shape.size(), 2u);
    EXPECT_EQ(shape[0], 2u);
    EXPECT_EQ(shape[1], 3u);
    ASSERT_EQ(data.size(), 7u);
    EXPECT_EQ(data[0], 99);
    EXPECT_EQ(data[1], -3);
    EXPECT_EQ(data[6], 8);

    npy::npy_data<long> scalar;
    scalar.data = { 7 };
    scalar.shape = {};
    const std::filesystem::path q = s.dir / "scalar.npy";
    npy::write_npy<long>(q.string(), scalar);
    std::vector<unsigned long> shape0;
    bool fortran0 = true;
    std::vector<long> data0;
    npy::LoadArrayFromNumpy<long>(q, shape0, fortran0, data0);
    EXPECT_TRUE(shape0.empty());
    EXPECT_FALSE(fortran0);
    ASSERT_EQ(data0.size(), 1u);
    EXPECT_EQ(data0[0], 7);
}

//----------------------------------------------------------------------------------------------
//cube.cpp
//----------------------------------------------------------------------------------------------

//A 2x1x8 cube stored six values per line: the second line carries the last two z values of the
//first row and the first four of the next x row, the third line the remaining four.
TEST(FittingIoCoverageCubeTests, ReadValuesSpillsAcrossRows)
{
    Scratch s("ReadValuesSpill");
    const std::filesystem::path p = s.dir / "spill.cube";
    write_text(p, "h\nh\nh\nh\nh\nh\n0 1 2 3 4 5\n6 7 8 9 10 11\n12 13 14 15\n");
    cube c({ 2, 1, 8 }, 0, true);
    std::ifstream f(p);
    ASSERT_TRUE(f.good());
    ASSERT_TRUE(c.read_values(f));
    for (int k = 0; k < 8; k++)
    {
        EXPECT_DOUBLE_EQ(c.get_value(0, 0, k), (double)k) << k;
        EXPECT_DOUBLE_EQ(c.get_value(1, 0, k), 8.0 + k) << k;
    }
}

//Six values for a 1x1x4 cube overflow the last row with nowhere to spill; a 1x1x8 cube whose
//only line ends the file leaves two values unread. Both report and return false.
TEST(FittingIoCoverageCubeTests, ReadValuesRejectsOverflowAndEof)
{
    Scratch s("ReadValuesReject");
    const std::filesystem::path over = s.dir / "over.cube";
    write_text(over, "h\nh\nh\nh\nh\nh\n0 1 2 3 4 5\n");
    {
        cube c({ 1, 1, 4 }, 0, true);
        std::ifstream f(over);
        StreamCapture out(std::cout);
        EXPECT_FALSE(c.read_values(f));
        EXPECT_NE(out.str().find("This should not happen! Read a value outside of range! Run_x: 0 Run_y: 0 rest2: 2 run_z: 6"), std::string::npos) << out.str();
        EXPECT_DOUBLE_EQ(c.get_value(0, 0, 3), 3.0);
    }
    const std::filesystem::path eof = s.dir / "eof.cube";
    write_text(eof, "h\nh\nh\nh\nh\nh\n0 1 2 3 4 5");
    {
        cube c({ 1, 1, 8 }, 0, true);
        std::ifstream f(eof);
        StreamCapture out(std::cout);
        EXPECT_FALSE(c.read_values(f));
        const std::string log = out.str();
        EXPECT_NE(log.find("This file ended before i read all expected values!"), std::string::npos) << log;
        EXPECT_NE(log.find("ENCOUNTERED EOF!"), std::string::npos);
        EXPECT_NE(log.find("x,y,reads1,z: 1 1 6,6"), std::string::npos);
        EXPECT_DOUBLE_EQ(c.get_value(0, 0, 5), 5.0);
    }
}

//The interactive super cube reads the three multipliers from std::cin, rejecting 0 and 25 for
//x before accepting 2; the 1x1x2 cube with one hydrogen doubles along x: two atoms shifted by
//size*vector = 0.5 bohr, two value lines of the same two values, and the file reads back.
TEST(FittingIoCoverageCubeTests, SuperCubeFromStdin)
{
    Scratch s("SuperCubeStdin");
    WFN parent(e_origin::NOT_YET_DEFINED);
    parent.push_back_atom("H", 0.1, 0.2, 0.3, 1);
    cube c({ 1, 1, 2 }, 1, true);
    c.give_parent_wfn(parent);
    for (int i = 0; i < 3; i++)
    {
        c.set_origin(i, 0.0);
        for (int j = 0; j < 3; j++)
            c.set_vector(i, j, i == j ? 0.5 : 0.0);
    }
    c.calc_dv();
    c.set_value(0, 0, 0, 1.0);
    c.set_value(0, 0, 1, 2.0);
    c.set_comment1("c1");
    c.set_comment2("c2");
    c.set_path(s.dir / "sc.cube");

    std::filesystem::path out_path;
    std::string prompts;
    {
        CinFeed feed("0 25 2\n1\n1\n");
        StreamCapture out(std::cout);
        out_path = c.super_cube();
        prompts = out.str();
    }
    EXPECT_EQ(prompts, "How many times in X-direction? This is unreasonable, try again! (between 1-20): This is unreasonable, try again! (between 1-20): Y-direction? Z-direction? ");
    EXPECT_EQ(out_path.extension().string(), ".cube_super");
    EXPECT_EQ(out_path.stem().string(), "sc");
    ASSERT_TRUE(std::filesystem::exists(out_path));

    const std::vector<std::string> lines = split_lines(read_text(out_path));
    ASSERT_EQ(lines.size(), 10u);
    EXPECT_EQ(lines[0], "c1 SUPER CUBE");
    EXPECT_EQ(lines[1], "c2");
    EXPECT_EQ(lines[2], "   2 0 0 0");
    EXPECT_EQ(lines[3], "     2    0.500000    0.000000    0.000000");
    EXPECT_EQ(lines[4], "     1    0.000000    0.500000    0.000000");
    EXPECT_EQ(lines[5], "     2    0.000000    0.000000    0.500000");
    EXPECT_EQ(lines[6], "    1    1.000000    0.100000    0.200000    0.300000");
    EXPECT_EQ(lines[7], "    1    1.000000    0.600000    0.200000    0.300000");
    EXPECT_EQ(lines[8], "  1.00000E+00  2.00000E+00");
    EXPECT_EQ(lines[9], "  1.00000E+00  2.00000E+00");

    WFN empty(e_origin::NOT_YET_DEFINED);
    std::ostringstream log;
    cube back(out_path, true, empty, log);
    EXPECT_EQ(back.get_size(0), 2);
    EXPECT_EQ(back.get_size(1), 1);
    EXPECT_EQ(back.get_size(2), 2);
    EXPECT_EQ(empty.get_ncen(), 2);
    EXPECT_NEAR(empty.get_atom_coordinate(1, 0), 0.6, 1e-6);
    EXPECT_DOUBLE_EQ(back.get_value(0, 0, 0), 1.0);
    EXPECT_DOUBLE_EQ(back.get_value(0, 0, 1), 2.0);
    EXPECT_DOUBLE_EQ(back.get_value(1, 0, 0), 1.0);
    EXPECT_DOUBLE_EQ(back.get_value(1, 0, 1), 2.0);
}

//----------------------------------------------------------------------------------------------
//mat_nos_math.cpp
//----------------------------------------------------------------------------------------------

//Complex values survive a 3-D reshape, a flatten, a 4-D reshape and a 2-D reshape from the
//4-D matrix, all in row-major order.
TEST(FittingIoCoverageMathTests, ComplexReshapeFlattenThreeAndFourD)
{
    const cvec values{ { 1, -1 }, { 2, -2 }, { 3, -3 }, { 4, -4 } };
    const cMatrix3 m3 = reshape<cMatrix3>(values, Shape3D(1, 2, 2));
    EXPECT_EQ(m3(0, 0, 0), values[0]);
    EXPECT_EQ(m3(0, 0, 1), values[1]);
    EXPECT_EQ(m3(0, 1, 0), values[2]);
    EXPECT_EQ(m3(0, 1, 1), values[3]);
    const cMatrix1 flat3 = flatten<cMatrix1>(m3);
    ASSERT_EQ(flat3.size(), 4u);
    for (size_t i = 0; i < 4; i++)
        EXPECT_EQ(flat3(i), values[i]) << i;

    const cMatrix4 m4 = reshape<cMatrix4>(values, Shape4D(1, 1, 2, 2));
    EXPECT_EQ(m4(0, 0, 1, 0), values[2]);
    EXPECT_EQ(m4(0, 0, 1, 1), values[3]);
    const cMatrix1 flat4 = flatten<cMatrix1>(m4);
    ASSERT_EQ(flat4.size(), 4u);
    for (size_t i = 0; i < 4; i++)
        EXPECT_EQ(flat4(i), values[i]) << i;

    const cMatrix2 m2 = reshape<cMatrix2>(m4, Shape2D(2, 2));
    EXPECT_EQ(m2(0, 1), values[1]);
    EXPECT_EQ(m2(1, 0), values[2]);
    const cMatrix3 m3b = reshape<cMatrix3>(m2, Shape3D(2, 1, 2));
    EXPECT_EQ(m3b(1, 0, 0), values[2]);
    EXPECT_EQ(m3b(1, 0, 1), values[3]);
}

//Integer 3-D and 4-D matrices flatten back to the row-major source order.
TEST(FittingIoCoverageMathTests, IntegerFlattenThreeAndFourD)
{
    const ivec values{ 5, -6, 7, -8, 9, -10 };
    const iMatrix3 m3 = reshape<iMatrix3>(values, Shape3D(1, 2, 3));
    EXPECT_EQ(m3(0, 1, 2), -10);
    const iMatrix1 flat3 = flatten<iMatrix1>(m3);
    ASSERT_EQ(flat3.size(), 6u);
    for (size_t i = 0; i < 6; i++)
        EXPECT_EQ(flat3(i), values[i]) << i;

    const iMatrix4 m4 = reshape<iMatrix4>(values, Shape4D(1, 1, 3, 2));
    EXPECT_EQ(m4(0, 0, 2, 0), 9);
    EXPECT_EQ(m4(0, 0, 0, 1), -6);
    const iMatrix1 flat4 = flatten<iMatrix1>(m4);
    ASSERT_EQ(flat4.size(), 6u);
    for (size_t i = 0; i < 6; i++)
        EXPECT_EQ(flat4(i), values[i]) << i;
}
