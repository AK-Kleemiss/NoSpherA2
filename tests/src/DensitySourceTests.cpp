#include "pch.h"

#include "core/density_source.h"
#include "core/gaussian_atom.h"
#include "core/spherical_density.h"
#include "core/basis_set.h"
#include "core/SALTED_utilities.h"
#include "core/properties.h"

#include <numeric>

vec2 sigma_to_alpha(const std::pair<vec, vec>& sigma);

//err_checkf exits with -1
constexpr unsigned ERROR_CHECK_EXIT_CODE =
#ifdef _WIN32
    static_cast<unsigned>(-1);
#else
    255u;
#endif

//density_source.h: one calculate_density / calculate_laplacian / calculate_rdg / calculate_eli for every density
//NoSpherA2 knows. Each source must agree with the point function it wraps, the derivatives with central
//differences of rho, and the batch with the point call; a Gaussian_Molecule must partition itself.
namespace
{
    //cartesian type index with exponents (a, b, c), as in WfnOpsTests
    int type_of(const int a, const int b, const int c)
    {
        int v[3];
        for (int t = 1; t <= 286; t++) {
            constants::type2vector(t, v);
            if (v[0] == a && v[1] == b && v[2] == c) return t;
        }
        return -1;
    }

    struct prim_spec
    {
        int cent, l[3];
        double e, c[3];
    };

    //He + H, three MOs, primitives up to i: the WFN of WfnOpsTests
    WFN make_wfn()
    {
        WFN w(e_origin::NOT_YET_DEFINED);
        w.push_back_atom("He", 0.0, 0.0, 0.0, 2);
        w.push_back_atom("H", 1.4, 0.3, -0.2, 1);
        w.push_back_MO(1, 2.0, -0.9);
        w.push_back_MO(2, 1.0, -0.3);
        w.push_back_MO(3, 0.0, 0.5);
        prim_spec specs[] = {
            {1, {0, 0, 0}, 1.20, {0.70, -0.20, 0.10}},
            {1, {1, 0, 0}, 0.80, {0.30, 0.40, -0.20}},
            {1, {0, 2, 0}, 0.90, {-0.20, 0.30, 0.50}},
            {1, {4, 0, 0}, 0.70, {0.15, -0.25, 0.20}},
            {2, {0, 0, 0}, 1.00, {0.40, 0.50, -0.10}},
            {2, {0, 0, 1}, 0.60, {0.25, -0.35, 0.30}},
            {2, {1, 1, 1}, 0.55, {0.20, 0.10, -0.40}},
        };
        for (prim_spec& p : specs)
            w.add_primitive(p.cent, type_of(p.l[0], p.l[1], p.l[2]), p.e, p.c);
        w.set_exp_cutoff();
        return w;
    }

    //O + H in combo_basis_fit with a positive s part, the fitted density of CrystalEnergyTests
    Gaussian_Molecule make_molecule(const bool s_only = false)
    {
        WFN oh(e_origin::NOT_YET_DEFINED);
        oh.push_back_atom("O", 0.0, 0.0, 0.0, 8);
        oh.push_back_atom("H", 1.8, 0.0, 0.0, 1);
        std::vector<std::shared_ptr<BasisSet>> basis{ BasisSetLibrary::get_basis_set("combo_basis_fit") };
        WFN aux = generate_aux_wfn(oh, basis);
        const aux_density_table t(aux.get_atoms());
        vec c(t.n_coef);
        for (int i = 0; i < t.n_coef; i++) c[i] = 0.3 * std::sin(1.0 + i);
        for (int s = 0; s < t.n_sh; s++) if (t.sh_l[s] == 0) c[t.coef_off[s]] += 1.0;
        //s_only: the positive s part alone, a density that is positive everywhere
        if (s_only) for (int s = 0; s < t.n_sh; s++) if (t.sh_l[s] != 0) std::fill(c.begin() + t.coef_off[s], c.begin() + t.coef_off[s] + 2 * t.sh_l[s] + 1, 0.0);
        return Gaussian_Molecule(std::move(aux), std::move(c));
    }

    const d3 probe{ 0.6, 0.2, -0.1 };
    const d3 centre{ 0.5, -0.4, 0.3 };

    //gradient and Laplacian of any source by central differences of its rho, the reference for the analytic ones
    template<class S> double fd_values(const S& s, const d3& p, d3& g, double& lap, const double h = 1e-4)
    {
        const double rho = calculate_density(s, p);
        lap = 0.0;
        for (int k = 0; k < 3; k++) {
            d3 pp = p, pm = p;
            pp[k] += h;
            pm[k] -= h;
            const double rp = calculate_density(s, pp), rm = calculate_density(s, pm);
            g[k] = (rp - rm) / (2 * h);
            lap += (rp + rm - 2 * rho) / (h * h);
        }
        return rho;
    }

    //The generic checks every source passes: derivatives against finite differences, the derived properties
    //against their definitions, the batch against the point call
    template<class S> void check_source(const S& s, const std::string& name)
    {
        d3 g, g_fd;
        double lap, lap_fd;
        const double rho = calculate_density(s, probe, g, lap);
        const double rho_fd = fd_values(s, probe, g_fd, lap_fd);
        EXPECT_NEAR(rho, calculate_density(s, probe), 1e-12 * rho) << name;
        EXPECT_NEAR(rho, rho_fd, 1e-12 * rho) << name;
        for (int k = 0; k < 3; k++) EXPECT_NEAR(g[k], g_fd[k], 1e-5 * std::abs(g_fd[k]) + 1e-7) << name << " k " << k;
        EXPECT_NEAR(lap, lap_fd, 1e-4 * std::abs(lap_fd) + 1e-6) << name;
        EXPECT_NEAR(calculate_laplacian(s, probe), lap, 1e-10 * std::abs(lap)) << name;
        EXPECT_NEAR(calculate_rdg(s, probe), constants::alpha_coef * std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]) / std::pow(rho, constants::c_43), 1e-12) << name;

        const int n = 5;
        const double x[n] = { 0.0, 0.6, -0.3, 1.1, 0.2 }, y[n] = { 0.0, 0.2, 0.4, -0.5, 0.9 }, z[n] = { 0.0, -0.1, 0.7, 0.3, -0.8 };
        double batch[n];
        calculate_density(s, n, x, y, z, batch);
        for (int p = 0; p < n; p++) EXPECT_DOUBLE_EQ(batch[p], calculate_density(s, d3{ x[p], y[p], z[p] })) << name << " p " << p;
        DensityBatch cb = density_batch(s);
        double via_cb[n];
        cb(n, x, y, z, via_cb);
        for (int p = 0; p < n; p++) EXPECT_DOUBLE_EQ(via_cb[p], batch[p]) << name << " p " << p;
    }
}

//The WFN goes through its orbitals: compute_dens, the computeValues kernel, computeLap and computeELI
TEST(DensitySourceTests, WfnMatchesItsOrbitalFunctions)
{
    const WFN w = make_wfn();
    EXPECT_DOUBLE_EQ(calculate_density(w, probe), w.compute_dens(probe));
    check_source(w, "WFN");
    EXPECT_NEAR(calculate_laplacian(w, probe), w.computeLap(probe), 1e-12);
    EXPECT_DOUBLE_EQ(calculate_eli(w, probe), w.computeELI(probe));
    //the old computeValues normGrad is the reduced density gradient
    double rho, ng, hess[9], elf, eli, lap;
    w.computeValues(probe, rho, ng, hess, elf, eli, lap);
    EXPECT_NEAR(calculate_rdg(w, probe), ng, 1e-12);
    EXPECT_NEAR(calculate_laplacian(w, probe), lap, 1e-12);
    EXPECT_NEAR(calculate_eli(w, probe), eli, 1e-12);
    //the gradient of the kernel is computeGrad
    d3 g, g_ref;
    double tau;
    w.computeValues(probe, rho, g, hess, tau);
    w.computeGrad(probe, g_ref);
    for (int k = 0; k < 3; k++) EXPECT_NEAR(g[k], g_ref[k], 1e-12);
    EXPECT_GT(tau, 0.0);
}

//The fitted density goes through its aux_density_table; the atoms of a molecule add up to it, ELI-D is the PC07 estimate
TEST(DensitySourceTests, GaussianMoleculeAndAtomsMatchTheTable)
{
    const Gaussian_Molecule M = make_molecule();
    EXPECT_DOUBLE_EQ(calculate_density(M, probe), M.rho(probe));
    check_source(M, "Gaussian_Molecule");
    check_source(M.atom(0), "Gaussian_Atom O");
    check_source(M.atom(1), "Gaussian_Atom H");
    EXPECT_NEAR(calculate_density(M.atom(0), probe) + calculate_density(M.atom(1), probe), calculate_density(M, probe), 1e-12);
    EXPECT_NEAR(calculate_laplacian(M, probe), M.lap(probe), 1e-10);
    EXPECT_NEAR(calculate_eli(M, probe), M.eli(probe), 1e-10);
    //a WFN and a Gaussian_Molecule through the same template
    const WFN w = make_wfn();
    EXPECT_GT(calculate_eli(w, probe), 0.0);
    EXPECT_GT(calculate_eli(M, probe), 0.0);
}

//A spherical model at a nucleus: rho(p) = rho(|p - c|), the radial differences agree with the 3D ones
TEST(DensitySourceTests, CentredRadialModelsMatchTheirRadialDensity)
{
    const Thakkar O(8);
    const Centred<Thakkar> cO{ O, centre };
    EXPECT_DOUBLE_EQ(calculate_density(cO, probe), O.get_radial_density(array_length(probe, centre)));
    check_source(cO, "Centred<Thakkar>");
    //at the nucleus the gradient vanishes and the Laplacian is 3 rho''
    d3 g;
    double lap;
    calculate_density(cO, centre, g, lap);
    for (int k = 0; k < 3; k++) EXPECT_EQ(g[k], 0.0);
    EXPECT_LT(lap, 0.0);

    const MBIS_Atom m(8, { 0.2, 0.5 }, { 2.0, 6.0 });
    const Centred<MBIS_Atom> cm{ m, centre };
    EXPECT_DOUBLE_EQ(calculate_density(cm, probe), m.get_radial_density(array_length(probe, centre)));
    check_source(cm, "Centred<MBIS_Atom>");
    //the analytic Slater derivatives: rho = sum p/(8 pi s^3) exp(-r/s)
    const double r = array_length(probe, centre);
    double d1 = 0.0, d2 = 0.0;
    for (int i = 0; i < 2; i++) {
        const double s = i == 0 ? 0.2 : 0.5, p = i == 0 ? 2.0 : 6.0;
        const double f = p * constants::INV_EIGHT_PI / (s * s * s) * std::exp(-r / s);
        d1 -= f / s;
        d2 += f / (s * s);
    }
    calculate_density(cm, probe, g, lap);
    EXPECT_NEAR(lap, d2 + 2 * d1 / r, 1e-6 * std::abs(d2));
    for (int k = 0; k < 3; k++) EXPECT_NEAR(g[k], d1 * (probe[k] - centre[k]) / r, 1e-6 * std::abs(d1));
}

//An EMBIS model at a nucleus: rho(p) = get_density(p - c); with diagonal alpha it is the MBIS atom
TEST(DensitySourceTests, CentredEmbisMatchesGetDensity)
{
    const vec sig = { 0.3, 0.8 }, pop = { 2.0, 6.0 };
    const EMBIS_Atom e(8, sigma_to_alpha({ sig, pop }), pop);
    const Centred<EMBIS_Atom> ce{ e, centre };
    EXPECT_DOUBLE_EQ(calculate_density(ce, probe), e.get_density({ probe[0] - centre[0], probe[1] - centre[1], probe[2] - centre[2] }));
    check_source(ce, "Centred<EMBIS_Atom>");
    const MBIS_Atom m(8, sig, pop);
    const Centred<MBIS_Atom> cm{ m, centre };
    EXPECT_NEAR(calculate_density(ce, probe), calculate_density(cm, probe), 1e-12 * calculate_density(cm, probe));
    EXPECT_NEAR(calculate_laplacian(ce, probe), calculate_laplacian(cm, probe), 1e-5 * std::abs(calculate_laplacian(cm, probe)));
    EXPECT_NEAR(calculate_rdg(ce, probe), calculate_rdg(cm, probe), 1e-6 * calculate_rdg(cm, probe));
    EXPECT_NEAR(calculate_eli(ce, probe), calculate_eli(cm, probe), 1e-5 * calculate_eli(cm, probe));
}

//TFVC, MBIS and EMBIS partition the fitted density itself on the grid: the populations of every scheme sum to the
//electron count and stay positive, and the density-based schemes differ from Hirshfeld
TEST(DensitySourceTests, GaussianMoleculePartitionsItself)
{
    using DensityFitting::CHARGE_SCHEME;
    const Gaussian_Molecule M = make_molecule(true);
    const double N = M.electrons();
    vec hirsh;
    for (const CHARGE_SCHEME scheme : { CHARGE_SCHEME::HIRSHFELD, CHARGE_SCHEME::TFVC, CHARGE_SCHEME::MBIS, CHARGE_SCHEME::EMBIS }) {
        const vec pop = M.populations(scheme);
        ASSERT_EQ(pop.size(), 2u) << (int)scheme;
        EXPECT_GT(pop[0], 0.0) << (int)scheme;
        EXPECT_GT(pop[1], 0.0) << (int)scheme;
        //grid quadrature; the MBIS pro-atoms stop short of the diffuse aux tail, 2.4e-3 there
        EXPECT_NEAR(pop[0] + pop[1], N, 5e-3) << (int)scheme;
        if (scheme == CHARGE_SCHEME::HIRSHFELD) hirsh = pop;
        else EXPECT_GT(std::abs(pop[0] - hirsh[0]), 1e-4) << (int)scheme;
    }
    //without the density on the grid the aux-only basis cannot be partitioned by a density-based scheme
    EXPECT_EXIT(DensityFitting::partition_multipole_rows(M.basis(), M.table(), CHARGE_SCHEME::MBIS, 0), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
    EXPECT_EQ(DensityFitting::partition_multipole_rows(M.basis(), M.table(), CHARGE_SCHEME::HIRSHFELD, 0).size(), 2u);
}

//The cube functions of properties.h on a fitted density: every loaded cube equals the point functions, the ESP the
//potential, and ELF (orbitals) is refused
TEST(DensitySourceTests, CubeFunctionsTakeAFittedDensity)
{
    const Gaussian_Molecule M = make_molecule(true);
    const int n = 5;
    const double h = 0.8;
    auto grid = [&](bool on) {
        cube c({ n, n, n }, 2, on);
        for (int k = 0; k < 3; k++) {
            c.set_origin(k, -0.5 * (n - 1) * h + (k == 0 ? 0.9 : 0.0));
            c.set_vector(k, k, h);
        }
        c.calc_dv();
        return c;
    };
    std::ostringstream quiet;
    const double radius = 5.0; //Angstrom, covers the box
    std::vector<cube> cubes;
    for (int t = 0; t < 16; t++) cubes.push_back(grid(t == cube_type::Rho || t == cube_type::RDG || t == cube_type::Lap || t == cube_type::Eli));
    Calc_Rho(cubes[cube_type::Rho], M, radius, quiet, false);
    Calc_Prop(cubes, M, radius, quiet, true, false);
    cube esp = grid(true), eli = grid(true);
    Calc_ESP(esp, M, radius, true, quiet, false);
    Calc_Eli(eli, M, radius, quiet);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++) {
                const d3 p = esp.get_pos(i, j, k);
                const double rho = M.rho(p);
                double H[9];
                calculate_hessian(M, p, H);
                EXPECT_NEAR(cubes[cube_type::Rho].get_value(i, j, k), get_lambda_1(H) < 0 ? -rho : rho, 1e-12);
                EXPECT_NEAR(cubes[cube_type::RDG].get_value(i, j, k), calculate_rdg(M, p), 1e-9);
                EXPECT_NEAR(cubes[cube_type::Lap].get_value(i, j, k), M.lap(p), 1e-10);
                EXPECT_NEAR(cubes[cube_type::Eli].get_value(i, j, k), M.eli(p), 1e-9);
                EXPECT_NEAR(eli.get_value(i, j, k), M.eli(p), 1e-9);
                EXPECT_NEAR(esp.get_value(i, j, k), M.esp(p), 1e-10);
            }
    //outside the radius the cube stays 0
    cube outside = grid(true);
    Calc_Rho(outside, M, 0.01, quiet, false);
    EXPECT_EQ(outside.get_value(0, 0, 0), 0.0);
    std::vector<cube> with_elf;
    for (int t = 0; t < 16; t++) with_elf.push_back(grid(t == cube_type::Rho || t == cube_type::Elf));
    EXPECT_EXIT(Calc_Prop(with_elf, M, radius, quiet, true, false), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
}
