#include "pch.h"

#include "core/convenience.h"
#include "core/constants.h"
#include "core/AtomGrid.h"
#include "core/GridManager.h"
#include "core/cell.h"
#include "core/cube.h"
#include "core/wfn_class.h"
#include "core/spherical_density.h"
#include "core/integrator.h"
#include "core/sphere_lebedev_rule.h"

#include <iomanip>
#include <numeric>
#include <set>

//The MBIS shell tables live in AtomGrid.cpp without a header
std::pair<vec, vec> get_shsig_shpop(const int &atom_type);
std::pair<vec2, vec> get_shalpha_shpop(const int &atom_type);
vec2 sigma_to_alpha(const std::pair<vec, vec> &sigma);

namespace
{
    //two hydrogens on z, one doubly occupied MO of two unit-exponent s primitives with
    //coefficient 1/2 each: rho = (phi1 + phi2)^2 / 2, so the electron count is
    //2 c^2 (pi/2)^(3/2) (2 + 2 exp(-R^2/2)) with R = 1.4 bohr
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
        const double R = 2.0 * kH2Half;
        return 2.0 * kH2Coef * kH2Coef * std::pow(constants::PI / 2.0, 1.5) * (2.0 + 2.0 * std::exp(-R * R / 2.0));
    }

    //4 pi Int r^2 f(r) g(r) dr on a logarithmic trapezoid grid
    double radial_integral(const std::function<double(double)> &f, const std::function<double(double)> &g,
        const double r_min, const double r_max, const int n)
    {
        const double dl = std::log(r_max / r_min) / (n - 1);
        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            const double r = r_min * std::exp(i * dl);
            const double t = r * r * r * f(r) * g(r);
            sum += (i == 0 || i == n - 1) ? 0.5 * t : t;
        }
        return constants::FOUR_PI * sum * dl;
    }

    double vec_total(const vec2 &v)
    {
        double s = 0.0;
        for (const vec &row : v) s += std::accumulate(row.begin(), row.end(), 0.0);
        return s;
    }

    //Becke JCP 88, 2547 (1988) eq. 19-21: the step s(nu) = (1 - p^k(nu)) / 2 with
    //p(x) = 3x/2 - x^3/2 iterated k times; for two centres the normalised weight of
    //centre a is s(nu_ab) itself
    double becke_step(double nu, const int iterations)
    {
        for (int i = 0; i < iterations; i++) nu = 1.5 * nu - 0.5 * nu * nu * nu;
        return 0.5 * (1.0 - nu);
    }
}

//every tabulated Lebedev order must have unit total weight and reproduce the exact
//sphere averages of x, x^2 and x^2 y^2 (degree 4, so only from order 14 on)
TEST(GridLebedevTests, EveryOrderIntegratesLowMoments)
{
    lebedev_sphere sphere;
    for (int order : constants::lebedev_table) {
        if (order == 0) continue; //the table is declared [33] with 32 entries, the last is zero padding
        vec x(order), y(order), z(order), w(order);
        sphere.ld_by_order(order, x.data(), y.data(), z.data(), w.data());
        double sw = 0.0, sx = 0.0, sxx = 0.0, sxxyy = 0.0;
        for (int i = 0; i < order; i++) {
            EXPECT_NEAR(x[i] * x[i] + y[i] * y[i] + z[i] * z[i], 1.0, 1e-12) << "order " << order;
            sw += w[i];
            sx += w[i] * x[i];
            sxx += w[i] * x[i] * x[i];
            sxxyy += w[i] * x[i] * x[i] * y[i] * y[i];
        }
        EXPECT_NEAR(sw, 1.0, 1e-12) << "order " << order;
        EXPECT_NEAR(sx, 0.0, 1e-12) << "order " << order;
        EXPECT_NEAR(sxx, 1.0 / 3.0, 1e-12) << "order " << order;
        if (order >= 14)
            EXPECT_NEAR(sxxyy, 1.0 / 15.0, 1e-12) << "order " << order;
    }
}

//the radial grid of a one-exponent atom integrates r^2 exp(-r^2) to sqrt(pi)/4 and the
//plain and omp accessors return the same table
TEST(GridAtomGridTests, RadialGridIntegratesGaussian)
{
    const double alpha_min[1] = { 0.5 };
    std::ostringstream log;
    AtomGrid g(1e-12, 50, 50, 1, 1.0, 0, alpha_min, log);
    const int n = g.get_num_radial_grid_points();
    ASSERT_GT(n, 10);
    vec r(n), w(n), r2(n), w2(n), r3(n), r4(n);
    g.get_radial_grid(r.data(), w.data());
    g.get_radial_grid_omp(r2.data(), w2.data());
    g.get_radial_distances(r3.data());
    g.get_radial_distances_omp(r4.data());
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        EXPECT_GT(w[i], 0.0);
        if (i > 0) EXPECT_GT(r[i], r[i - 1]);
        EXPECT_EQ(r[i], r2[i]);
        EXPECT_EQ(w[i], w2[i]);
        EXPECT_EQ(r[i], r3[i]);
        EXPECT_EQ(r[i], r4[i]);
        sum += w[i] * std::exp(-r[i] * r[i]);
    }
    EXPECT_NEAR(sum, std::sqrt(constants::PI) / 4.0, 1e-7);
    EXPECT_EQ(g.get_num_grid_points(), n * 50);
}

//the full atomic grid of a single centre carries 4 pi in its weights: Int exp(-r^2) = pi^(3/2)
TEST(GridAtomGridTests, FullGridIntegratesGaussianWithFourPi)
{
    const double alpha_min[1] = { 0.5 };
    std::ostringstream log;
    AtomGrid g(1e-12, 110, 110, 1, 1.0, 0, alpha_min, log);
    const int n = g.get_num_grid_points();
    const double cx = 0.3, cy = -0.2, cz = 0.1;
    const int Z = 1;
    vec x(n), y(n), z(n), aw(n), bw(n), tw(n);
    g.get_grid(1, 0, &cx, &cy, &cz, &Z, x.data(), y.data(), z.data(), aw.data(), bw.data(), tw.data(), vec());
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        const double dx = x[i] - cx, dy = y[i] - cy, dz = z[i] - cz;
        sum += aw[i] * std::exp(-(dx * dx + dy * dy + dz * dz));
        EXPECT_EQ(bw[i], aw[i]);
        EXPECT_EQ(tw[i], aw[i]);
        EXPECT_NEAR(g.get_gridx(i), x[i] - cx, 1e-12);
    }
    const double ref = std::pow(constants::PI, 1.5);
    EXPECT_NEAR(sum, ref, 1e-7 * ref);
}

//Becke weights of two equal centres are the three-fold iterated step of
//mu = (r_a - r_b) / R_ab (hand formula below), 1/2 on the mirror plane, and the
//internally computed pair table matches an explicit one; without a chi table the
//TFVC weight is 1/2 per centre by construction (pa_tv stays 1), so only the sum is checked
TEST(GridAtomGridTests, BeckeWeightsOfTwoCentresSumToOne)
{
    const int Z[2] = { 1, 1 };
    const double x[2] = { 0.0, 0.0 }, y[2] = { 0.0, 0.0 }, z[2] = { -kH2Half, kH2Half };
    const double table[4] = { 0.0, 2.0 * kH2Half, 2.0 * kH2Half, 0.0 };
    vec pa_b(2), pa_tv(2);
    const vec chi;
    const double pts[3][3] = { { 0.4, -0.3, 0.0 }, { 0.1, 0.2, 0.5 }, { -1.0, 0.7, -1.3 } };
    for (const auto &p : pts) {
        const auto w0 = get_integration_weights(2, Z, x, y, z, 0, p[0], p[1], p[2], pa_b, pa_tv, chi);
        const auto w1 = get_integration_weights(2, Z, x, y, z, 1, p[0], p[1], p[2], pa_b, pa_tv, chi);
        const auto w0t = get_integration_weights(2, Z, x, y, z, 0, p[0], p[1], p[2], pa_b, pa_tv, chi, table);
        EXPECT_NEAR(w0[0] + w1[0], 1.0, 1e-12);
        EXPECT_NEAR(w0[1] + w1[1], 1.0, 1e-12);
        EXPECT_EQ(w0[0], w0t[0]);
        EXPECT_EQ(w0[1], w0t[1]);
        const double ra = std::sqrt(p[0] * p[0] + p[1] * p[1] + (p[2] - z[0]) * (p[2] - z[0]));
        const double rb = std::sqrt(p[0] * p[0] + p[1] * p[1] + (p[2] - z[1]) * (p[2] - z[1]));
        EXPECT_NEAR(w0[0], becke_step((ra - rb) / (2.0 * kH2Half), 3), 1e-12);
        EXPECT_NEAR(w0[1], 0.5, 1e-12);
    }
    const auto mid = get_integration_weights(2, Z, x, y, z, 0, 0.4, -0.3, 0.0, pa_b, pa_tv, chi);
    EXPECT_NEAR(mid[0], 0.5, 1e-12);
    EXPECT_NEAR(mid[1], 0.5, 1e-12);

    //the two atomic grids together integrate a Gaussian at the origin to pi^(3/2)
    const double alpha_min[1] = { 0.5 };
    std::ostringstream log;
    AtomGrid g(1e-12, 110, 110, 1, 1.0, 0, alpha_min, log);
    const int n = g.get_num_grid_points();
    vec gx(n), gy(n), gz(n), aw(n), bw(n), tw(n);
    double becke = 0.0, tfvc = 0.0;
    for (int c = 0; c < 2; c++) {
        g.get_grid(2, c, x, y, z, Z, gx.data(), gy.data(), gz.data(), aw.data(), bw.data(), tw.data(), chi);
        for (int i = 0; i < n; i++) {
            const double f = std::exp(-(gx[i] * gx[i] + gy[i] * gy[i] + gz[i] * gz[i]));
            becke += bw[i] * f;
            tfvc += tw[i] * f;
        }
    }
    const double ref = std::pow(constants::PI, 1.5);
    EXPECT_NEAR(becke, ref, 1e-6 * ref);
    EXPECT_NEAR(tfvc, ref, 1e-6 * ref);
}

//with a chi table the TFVC branch runs: unit chi makes JCP 139, 071103 (2013) eq. 7
//collapse to nu = mu, so the TFVC weight is the four-fold iterated Becke step of mu
//(three-fold for Becke); for H next to C the Becke size adjustment A3 with
//chi = R_H / R_C = 0.35 / 0.70 gives a_HC = 3/8, so at the midpoint (mu = 0) the
//C weight is 1 - s(3/8)
TEST(GridAtomGridTests, TfvcWeightsWithUnitChi)
{
    const int Z[2] = { 1, 1 };
    const double x[2] = { 0.0, 0.0 }, y[2] = { 0.0, 0.0 }, z[2] = { -kH2Half, kH2Half };
    vec pa_b(2), pa_tv(2);
    const vec chi(4, 1.0);
    const auto mid0 = get_integration_weights(2, Z, x, y, z, 0, 0.2, 0.1, 0.0, pa_b, pa_tv, chi);
    const auto mid1 = get_integration_weights(2, Z, x, y, z, 1, 0.2, 0.1, 0.0, pa_b, pa_tv, chi);
    EXPECT_NEAR(mid0[0], 0.5, 1e-12);
    EXPECT_NEAR(mid0[1], 0.5, 1e-12);
    EXPECT_NEAR(mid0[0] + mid1[0], 1.0, 1e-12);
    EXPECT_NEAR(mid0[1] + mid1[1], 1.0, 1e-12);
    const auto off0 = get_integration_weights(2, Z, x, y, z, 0, 0.0, 0.0, -0.3, pa_b, pa_tv, chi);
    const auto off1 = get_integration_weights(2, Z, x, y, z, 1, 0.0, 0.0, -0.3, pa_b, pa_tv, chi);
    const double mu = (0.4 - 1.0) / (2.0 * kH2Half);
    EXPECT_NEAR(off0[0], becke_step(mu, 3), 1e-12);
    EXPECT_NEAR(off0[1], becke_step(mu, 4), 1e-12);
    EXPECT_NE(off0[0], off0[1]);
    EXPECT_NEAR(off0[0] + off1[0], 1.0, 1e-12);
    EXPECT_NEAR(off0[1] + off1[1], 1.0, 1e-12);

    const int ZHC[2] = { 1, 6 };
    const double zhc[2] = { -1.0, 1.0 };
    const auto hc = get_integration_weights(2, ZHC, x, y, zhc, 1, 0.0, 0.0, 0.0, pa_b, pa_tv, vec());
    EXPECT_NEAR(hc[0], 1.0 - becke_step(3.0 / 8.0, 3), 1e-12);
}

//the MBIS shell tables hold Z electrons in MBIS_function[Z] shells and the diagonal
//alpha tensors are the inverse squared sigmas
TEST(GridAtomGridTests, ShellTablesPopulationsSumToZ)
{
    for (int Z : { 1, 6, 11, 20, 40, 60 }) {
        const auto sp = get_shsig_shpop(Z);
        const auto ap = get_shalpha_shpop(Z);
        const size_t shells = constants::MBIS_function[Z];
        EXPECT_EQ(sp.first.size(), shells) << "Z " << Z;
        EXPECT_EQ(sp.second.size(), shells) << "Z " << Z;
        EXPECT_EQ(ap.first.size(), shells) << "Z " << Z;
        EXPECT_NEAR(std::accumulate(sp.second.begin(), sp.second.end(), 0.0), Z, 1e-9) << "Z " << Z;
        EXPECT_NEAR(std::accumulate(ap.second.begin(), ap.second.end(), 0.0), Z, 1e-9) << "Z " << Z;
        const vec2 alpha = sigma_to_alpha(sp);
        for (size_t m = 0; m < shells; m++) {
            EXPECT_GT(sp.first[m], 0.0);
            ASSERT_EQ(alpha[m].size(), 6u);
            for (int k : { 0, 3, 5 })
                EXPECT_NEAR(alpha[m][k], 1.0 / (sp.first[m] * sp.first[m]), 1e-9 * alpha[m][k]) << "Z " << Z;
            for (int k : { 1, 2, 4 })
                EXPECT_EQ(alpha[m][k], 0.0);
            for (int k = 0; k < 6; k++)
                EXPECT_NEAR(ap.first[m][k], alpha[m][k], 1e-9 * std::abs(alpha[m][k]) + 1e-12) << "Z " << Z;
        }
    }
}

//the radial grid parameter solvers return the roots of the TCA 106, 178 error estimates
TEST(GridAtomGridTests, RadialParameterSolversHitTheirErrorTargets)
{
    const double err = 1e-8;
    const double r_inner = get_r_inner(err, 2.0);
    EXPECT_NEAR(r_inner, std::sqrt(std::exp((1.9 - std::log(1.0 / err)) * 2.0 / 3.0) / 2.0), 1e-14);
    EXPECT_LT(get_r_inner(err, 4.0), r_inner);

    const double alpha = 1.0;
    const double r_outer = get_r_outer(err, alpha, 0, 4.0);
    const double f_outer = std::tgamma(1.5) * std::sqrt(alpha * r_outer * r_outer) * std::exp(-alpha * r_outer * r_outer);
    EXPECT_NEAR(f_outer, err, 1e-6 * err);
    EXPECT_GT(get_r_outer(1e-12, alpha, 0, 4.0), r_outer);
    EXPECT_LT(get_r_outer(err, 2.0 * alpha, 0, 4.0), r_outer);

    const double h = get_h(err, 0, 0.3);
    const double f_h = constants::C0 / h * std::exp(-constants::PI2 / (2.0 * h));
    EXPECT_NEAR(f_h, err, 1e-6 * err);
    EXPECT_LT(get_h(1e-12, 0, 0.3), h);
}

//the Thakkar core form factor at k -> 0 counts the core electrons of every closed
//core the mode-1 table knows, and is continuous towards small k
TEST(GridSphericalTests, ThakkarCoreFormFactorAtZeroCountsCoreElectrons)
{
    const std::pair<int, int> cases[] = { { 26, 10 }, { 36, 18 }, { 47, 28 }, { 54, 46 }, { 80, 60 }, { 86, 78 } };
    for (const auto &[Z, core] : cases) {
        Thakkar t(Z);
        const double f0 = t.get_core_form_factor(1e-6, core);
        EXPECT_NEAR(f0, core, 2e-2) << "Z " << Z;
        EXPECT_NEAR(t.get_core_form_factor(1e-5, core), f0, 1e-6 * f0) << "Z " << Z;
        EXPECT_LT(t.get_core_form_factor(2.0, core), f0) << "Z " << Z;
        EXPECT_NEAR(t.get_form_factor(1e-6), Z, 2e-2) << "Z " << Z;
    }
}

//the exact k == 0 branch of the Thakkar form factor is the k -> 0 limit of the k > 0 branch,
//and the hydrogen default counts one electron there
TEST(GridSphericalTests, ThakkarFormFactorExactlyAtZeroMatchesTheLimit)
{
    Thakkar t(26);
    EXPECT_NEAR(t.get_core_form_factor(0.0, 10), t.get_core_form_factor(1e-6, 10), 1e-3);
    EXPECT_NEAR(t.get_form_factor(0.0), t.get_form_factor(1e-6), 1e-3);
    EXPECT_NEAR(Thakkar().get_form_factor(0.0), 1.0, 1e-6);
}

//ECP mode 2 (xTB) knows the 36-, 54- and 68-electron cores as well
TEST(GridSphericalTests, ThakkarCoreFormFactorEcpMode2)
{
    const std::pair<int, int> cases[] = { { 26, 10 }, { 36, 36 }, { 54, 54 }, { 80, 68 }, { 86, 78 } };
    for (const auto &[Z, core] : cases) {
        Thakkar t(Z, 2);
        EXPECT_NEAR(t.get_core_form_factor(1e-6, core), core, 2e-2) << "Z " << Z;
    }
}

//the real-space core density of Ag integrates to its 28 core electrons and to the
//analytic k -> 0 form factor of the same orbitals
TEST(GridSphericalTests, ThakkarCoreDensityIntegratesToCoreCount)
{
    Thakkar t(47);
    const auto core = [&t](const double r) { return t.get_core_density(r, 28); };
    const auto one = [](const double) { return 1.0; };
    const double n_core = radial_integral(core, one, 1e-6, 40.0, 20000);
    EXPECT_NEAR(n_core, 28.0, 2e-2);
    EXPECT_NEAR(n_core, t.get_core_form_factor(1e-6, 28), 1e-3);
    EXPECT_LT(t.get_core_density(3.0, 28), t.get_radial_density(3.0));
}

//the default constructor is the hydrogen table: same density as Thakkar(1), one electron
TEST(GridSphericalTests, ThakkarDefaultCtorMatchesHydrogen)
{
    Thakkar h_default;
    Thakkar h(1);
    for (double r : { 0.1, 0.5, 1.0, 2.0 })
        EXPECT_NEAR(h_default.get_radial_density(r), h.get_radial_density(r), 1e-14 * h.get_radial_density(r));
    const auto rho = [&h_default](const double r) { return h_default.get_radial_density(r); };
    const auto one = [](const double) { return 1.0; };
    EXPECT_NEAR(radial_integral(rho, one, 1e-6, 40.0, 20000), 1.0, 1e-3);
    EXPECT_NEAR(h_default.get_form_factor(1e-6), 1.0, 1e-3);
    EXPECT_EQ(h_default.get_atomic_number(), 1);
}

//an MBIS atom is a sum of normalised Slater shells: hand value at one radius, the
//population comes back from the radial integral, and the table interpolator is close
TEST(GridSphericalTests, MbisAtomDensityIsSlaterSum)
{
    const vec sig = { 0.2, 0.5 };
    const vec pop = { 2.0, 4.0 };
    MBIS_Atom a(6, sig, pop);
    const double r = 0.3;
    double ref = 0.0;
    for (int m = 0; m < 2; m++)
        ref += pop[m] / (8.0 * constants::PI * std::pow(sig[m], 3)) * std::exp(-r / sig[m]);
    EXPECT_NEAR(a.get_radial_density(r), ref, 1e-14 * ref);
    const auto rho = [&a](const double x) { return a.get_radial_density(x); };
    const auto one = [](const double) { return 1.0; };
    EXPECT_NEAR(radial_integral(rho, one, 1e-7, 40.0, 40000), 6.0, 1e-6);
    a.make_interpolator(1.0005, 1e-7);
    for (double x : { 0.3, 1.0, 2.5 })
        EXPECT_NEAR(a.get_interpolated_density(x), a.get_radial_density(x), 1e-5 * a.get_radial_density(x));
    EXPECT_EQ(a.get_interpolated_density(1e3), 0.0);
}

//the linear table interpolators interpolate on [x[nr], x[nr+1]], the interval
//log_spline_index brackets dist with
TEST(GridSphericalTests, LinearInterpolatorAnchorsOnTheLowerNode)
{
    const vec x = { 1.0, 2.0, 4.0, 8.0 };
    const vec y = { 8.0, 4.0, 2.0, 1.0 };
    EXPECT_NEAR(linear_interpolate_spherical_density(y, x, 3.0, std::log(2.0), 1.0), 3.0, 1e-12);

    Thakkar c(6);
    c.make_interpolator(1.015, 1e-7);
    for (double r : { 1.0, 2.0, 3.0 })
        EXPECT_NEAR(c.get_interpolated_density(r), c.get_radial_density(r), 2e-3 * c.get_radial_density(r));

    MBIS_Atom a(6, { 0.2, 0.5 }, { 2.0, 4.0 });
    a.make_interpolator(1.015, 1e-7);
    for (double r : { 1.0, 2.0, 3.0 })
        EXPECT_NEAR(a.get_interpolated_density(r), a.get_radial_density(r), 2e-3 * a.get_radial_density(r));
}

//an EMBIS atom with the diagonal alpha of sigma_to_alpha is the isotropic MBIS atom
TEST(GridSphericalTests, EmbisDiagonalAlphaMatchesMbis)
{
    const vec sig = { 0.3, 0.8 };
    const vec pop = { 2.0, 4.0 };
    MBIS_Atom iso(6, sig, pop);
    EMBIS_Atom aniso(6, sigma_to_alpha({ sig, pop }), pop);
    const d3 pts[3] = { { 0.1, 0.2, 0.3 }, { -0.5, 0.0, 0.4 }, { 1.0, -1.0, 0.5 } };
    for (const d3 &p : pts) {
        const double r = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
        const double ref = iso.get_radial_density(r);
        EXPECT_NEAR(aniso.get_density(p), ref, 1e-12 * ref);
    }
    //a default EMBIS_Atom is hydrogen with empty tables, get_density on it reads alpha[0] out of bounds
}

//the Gaussian-expanded ECP core corrections transform to their form factors:
//4 pi Int r^2 rho sin(kr)/(kr) for every mode that has a table for the element; the tables are
//difference densities (they integrate to ~0), and the xtb table (mode 2) is empty for every element
TEST(GridSphericalTests, SphericalGaussianDensityFormFactorMatchesRadialIntegral)
{
    const std::pair<int, int> cases[] = { { 1, 37 }, { 3, 6 } };
    for (const auto &[mode, Z] : cases) {
        Spherical_Gaussian_Density d(Z, mode);
        EXPECT_NE(d.get_radial_density(0.5), 0.0) << "mode " << mode;
        const auto rho = [&d](const double r) { return d.get_radial_density(r); };
        const auto one = [](const double) { return 1.0; };
        const double k = 0.5;
        const auto sinc = [k](const double r) { return std::sin(k * r) / (k * r); };
        const double n = radial_integral(rho, one, 1e-6, 30.0, 40000);
        EXPECT_NEAR(d.get_form_factor(1e-6), n, 1e-6 * std::abs(n) + 1e-9) << "mode " << mode;
        const double fk = radial_integral(rho, sinc, 1e-6, 30.0, 40000);
        EXPECT_NEAR(d.get_form_factor(k), fk, 1e-6 * std::abs(fk) + 1e-9) << "mode " << mode;
    }
    for (const int mode : { 1, 2 }) {
        Spherical_Gaussian_Density none(6, mode);
        EXPECT_EQ(none.get_radial_density(0.5), 0.0) << "mode " << mode;
        EXPECT_EQ(none.get_form_factor(0.5), 0.0) << "mode " << mode;
    }
}

//the k = 0 branch must be the limit of the k > 0 branch, c (pi/z)^(3/2) per Gaussian
TEST(GridSphericalTests, SphericalGaussianDensityFormFactorContinuousAtZero)
{
    Spherical_Gaussian_Density d(37, 1);
    const double f_small = d.get_form_factor(1e-8);
    EXPECT_NEAR(d.get_form_factor(0.0), f_small, 1e-6 * f_small);
}

//the header spline helpers reproduce a smooth function on a logarithmic table, the index
//lookup brackets dist, and the Thakkar spline interpolator matches the analytic density
TEST(GridSphericalTests, CubicSplineHelpersReproduceSmoothFunction)
{
    const double start = 1e-3, incr = 1.01, lincr = std::log(incr);
    vec x, y;
    for (double r = start; r < 20.0; r *= incr) {
        x.push_back(r);
        y.push_back(std::exp(-r));
    }
    const vec y2 = natural_cubic_spline_second_derivatives(x, y);
    ASSERT_EQ(y2.size(), x.size());
    EXPECT_EQ(y2.front(), 0.0);
    EXPECT_EQ(y2.back(), 0.0);
    for (double r : { 0.37, 1.3, 2.9, 7.7 }) {
        const int nr = log_spline_index(x, r, lincr, start);
        EXPECT_LE(x[nr], r);
        EXPECT_GT(x[nr + 1], r);
        EXPECT_NEAR(cubic_spline_interpolate_spherical_density(y, x, y2, r, lincr, start), std::exp(-r), 1e-6 * std::exp(-r));
    }
    EXPECT_EQ(cubic_spline_interpolate_spherical_density(y, x, y2, 50.0, lincr, start), 0.0);
    EXPECT_EQ(cubic_spline_interpolate_spherical_density(y, x, y2, 1e-5, lincr, start), y[0]);
    EXPECT_EQ(natural_cubic_spline_second_derivatives({ 1.0, 2.0 }, { 1.0, 1.0 }).size(), 2u);

    Thakkar c(6);
    c.make_interpolator(1.005, 1e-7);
    EXPECT_NEAR(c.get_lincr(), std::log(1.005), 1e-15);
    EXPECT_EQ(c.get_start(), 1e-7);
    EXPECT_EQ(c.get_radial_dist().size(), c.get_radial_density_table().size());
    for (double r : { 0.5, 1.0, 2.0 })
        EXPECT_NEAR(c.get_interpolated_density_spline(r), c.get_radial_density(r), 1e-4 * c.get_radial_density(r));
}

//the configuration helpers: cutoff by accuracy and the scheme names
TEST(GridManagerTests, ConfigurationCutoffAndPartitionName)
{
    GridConfiguration cfg;
    EXPECT_EQ(cfg.accuracy, 2);
    EXPECT_EQ(cfg.getCutoff(), 1e-10);
    cfg.accuracy = 3;
    EXPECT_EQ(cfg.getCutoff(), 1e-14);
    cfg.accuracy = 4;
    EXPECT_EQ(cfg.getCutoff(), 1e-30);
    const std::pair<PartitionType, std::string> names[] = {
        { PartitionType::Hirshfeld, "Hirshfeld" }, { PartitionType::Becke, "Becke" }, { PartitionType::TFVC, "TFVC" },
        { PartitionType::MBIS, "MBIS" }, { PartitionType::EMBIS, "EMBIS" }, { PartitionType::RI, "RI" } };
    for (const auto &[type, name] : names) {
        cfg.partition_type = type;
        EXPECT_EQ(cfg.getPartitionName(), name);
    }
}

//GridData sizing and the two static atom helpers skip the charge-119 dummies
TEST(GridManagerTests, GridDataResizeAndStatics)
{
    GridData gd;
    gd.resizeForAtoms(3);
    ASSERT_EQ(gd.atomic_grids.size(), 3u);
    EXPECT_EQ(gd.atomic_grids[2].size(), 11u);
    EXPECT_EQ(gd.num_points_per_atom.size(), 3u);
    gd.resizeForAtoms(2, true);
    ASSERT_EQ(gd.helper_grids.size(), 2u);
    EXPECT_EQ(gd.helper_grids[1].size(), 11u);
    EXPECT_EQ(gd.helper_num_points_per_atom.size(), 2u);
    gd.clear();
    EXPECT_EQ(gd.atomic_grids.size(), 0u);
    EXPECT_EQ(gd.total_points, 0);

    WFN w(e_origin::NOT_YET_DEFINED);
    w.push_back_atom("C1", 0.0, 0.0, 0.0, 6);
    w.push_back_atom("H1", 1.0, 0.0, 0.0, 1);
    w.push_back_atom("Q1", 2.0, 0.0, 0.0, 119);
    w.push_back_atom("H2", 3.0, 0.0, 0.0, 1);
    EXPECT_EQ(GridManager::identifyAtomTypes(w), (ivec{ 1, 6 }));
    const bvec need = GridManager::determineAtomsNeedingGrids(w, { 0, 2, 7, -1 });
    EXPECT_EQ(need, (bvec{ true, false, false, false }));
}

//with all_charges every scheme runs on the helper grids: Becke, TFVC, Hirshfeld, MBIS
//and EMBIS all integrate the H2 density to its analytic electron count, split evenly
TEST(GridManagerTests, EveryPartitionIntegratesTheDensity)
{
    const WFN w = make_h2();
    GridConfiguration cfg;
    cfg.accuracy = 3;
    cfg.all_charges = true;
    GridManager gm(cfg);
    std::ostringstream log;
    gm.setup3DGridsForMolecule(w, { 0, 1 }, {}, cell(), false, log);
    EXPECT_TRUE(gm.getNeedsHelper());
    EXPECT_GT(gm.getTotalGridPoints(), 1000);
    const PartitionResults res = gm.calculatePartitionedCharges(w);
    const double N = h2_electrons();
    ASSERT_EQ(res.overall_charges.size(), 5u);
    for (int s = 0; s < 5; s++) {
        EXPECT_NEAR(res.overall_charges[s], N, 5e-3) << "scheme " << s;
        ASSERT_EQ(res.atom_charges[s].size(), 2u);
        EXPECT_NEAR(res.atom_charges[s][0], res.atom_charges[s][1], 1e-6) << "scheme " << s;
    }
    EXPECT_NEAR(res.overall_charges[PartitionResults::S_BECKE], N, 1e-4);
    EXPECT_NEAR(res.overall_charges[PartitionResults::S_HIRSH], N, 1e-4);
    std::ostringstream table;
    gm.printChargeTable({ "H1", "H2" }, w, { 0, 1 }, table, res);
    EXPECT_NE(table.str().find("Becke"), std::string::npos);
    EXPECT_NE(table.str().find("EMBIS"), std::string::npos);
    EXPECT_NE(table.str().find("H2"), std::string::npos);
}

//the Becke monopole equals the Becke population, the z dipole of atom 0 is the plain
//Sum_p rho w (z - z_0) sqrt(3 / 4 pi) over its grid columns (independent of the
//spherical_harmonic table and of the l*l+l+m indexing), the two z dipoles of the
//symmetric pair cancel and the x and y dipoles vanish by symmetry
TEST(GridManagerTests, BeckeMultipolesMonopoleMatchesCharge)
{
    const WFN w = make_h2();
    GridConfiguration cfg;
    cfg.accuracy = 2;
    cfg.partition_type = PartitionType::Becke;
    GridManager gm(cfg);
    std::ostringstream log;
    gm.setup3DGridsForMolecule(w, { 0, 1 }, {}, cell(), false, log);
    EXPECT_FALSE(gm.getNeedsHelper());
    const PartitionResults res = gm.calculatePartitionedCharges(w);
    const vec2 Q = gm.calculatePartitionedMultipoles(w, 1);
    ASSERT_EQ(Q.size(), 2u);
    ASSERT_EQ(Q[0].size(), 4u);
    for (int a = 0; a < 2; a++) {
        EXPECT_NEAR(std::sqrt(constants::FOUR_PI) * Q[a][0], res.atom_charges[PartitionResults::S_BECKE][a], 1e-10);
        EXPECT_NEAR(Q[a][1], 0.0, 1e-8);
        EXPECT_NEAR(Q[a][3], 0.0, 1e-8);
    }
    EXPECT_GT(std::abs(Q[0][2]), 1e-3);
    EXPECT_NEAR(Q[0][2] + Q[1][2], 0.0, 1e-8);
    EXPECT_NEAR(res.overall_charges[PartitionResults::S_BECKE], h2_electrons(), 1e-3);
    const vec2 &g0 = gm.getGridData().atomic_grids[0];
    double q10 = 0.0;
    for (int p = 0; p < gm.getNumPointsForAtom(0); p++)
        q10 += g0[GridData::WFN_DENSITY][p] * g0[GridData::BECKE_WEIGHT][p] * (g0[GridData::Z][p] + kH2Half);
    q10 *= std::sqrt(3.0 / constants::FOUR_PI);
    EXPECT_NEAR(Q[0][2], q10, 1e-10);
}

//a second identical setup reuses the grid (only densities are re-evaluated) and the
//g(r) column is Sum_MO phi^2 without coefficients: twice rho for this wavefunction
TEST(GridManagerTests, RepeatedSetupReusesGridAndGDensityIsTwiceRho)
{
    const WFN w = make_h2();
    GridConfiguration cfg;
    cfg.accuracy = 2;
    cfg.partition_type = PartitionType::Hirshfeld;
    GridManager gm(cfg);
    std::ostringstream log;
    gm.setup3DGridsForMolecule(w, { 0, 1 }, {}, cell(), false, log);
    const int points = gm.getTotalGridPoints();
    vec2 d1, d2, d3, rho, g;
    gm.getDensityVectors(w, { 0, 1 }, d1, d2, d3, rho, false);
    gm.setup3DGridsForMolecule(w, { 0, 1 }, {}, cell(), true, log);
    EXPECT_EQ(gm.getTotalGridPoints(), points);
    EXPECT_EQ(gm.getNumPointsForAtom(0) + gm.getNumPointsForAtom(1), points);
    gm.getDensityVectors(w, { 0, 1 }, d1, d2, d3, g, true);
    ASSERT_EQ(rho.size(), 2u);
    ASSERT_EQ(g.size(), 2u);
    const double N = h2_electrons();
    EXPECT_NEAR(vec_total(rho), N, 1e-3);
    EXPECT_NEAR(vec_total(g), 2.0 * N, 2e-3);
    EXPECT_NEAR(vec_total(g), 2.0 * vec_total(rho), 1e-8 * N);
    //the displacement columns are relative to the owning nucleus: the first moments of
    //rho w about it equal the same sums over the grid columns with z_0 = -kH2Half taken
    //off (an absolute z column would be off by z_0 N_A); getDensityVectors drops points
    //below the cutoff, so the sums are compared rather than single indices
    ASSERT_EQ(d3[0].size(), rho[0].size());
    const vec2 &g0 = gm.getGridData().atomic_grids[0];
    double m[3] = { 0.0, 0.0, 0.0 }, ref[3] = { 0.0, 0.0, 0.0 };
    for (size_t p = 0; p < rho[0].size(); p++) {
        m[0] += rho[0][p] * d1[0][p];
        m[1] += rho[0][p] * d2[0][p];
        m[2] += rho[0][p] * d3[0][p];
    }
    for (int p = 0; p < gm.getNumPointsForAtom(0); p++) {
        const double f = g0[GridData::WFN_DENSITY][p] * g0[GridData::HIRSH_WEIGHT][p];
        ref[0] += f * g0[GridData::X][p];
        ref[1] += f * g0[GridData::Y][p];
        ref[2] += f * (g0[GridData::Z][p] + kH2Half);
    }
    for (int i = 0; i < 3; i++)
        EXPECT_NEAR(m[i], ref[i], 1e-8) << "axis " << i;
    EXPECT_GT(std::abs(m[2]), 1e-3);
}

//the 1D grid runs from atom 1 to atom 2 padded on both ends with unit weights, its
//Hirshfeld weight is the Thakkar hydrogen density ratio rho_H(r_1) / (rho_H(r_1) + rho_H(r_2))
//(the manager interpolates that table on a cubic spline, hence 1e-6), its Becke weight
//the hand Becke step, and its density column is compute_dens
TEST(GridManagerTests, OneDimensionalGridWeightsSumToOne)
{
    const WFN w = make_h2();
    GridManager gm;
    const int n = 21;
    gm.setup1DGridsForMolecule(w, 0, 1, n, 0.5);
    const GridData &gd = gm.getGridData();
    ASSERT_EQ(gd.atomic_grids.size(), 2u);
    const vec2 &g0 = gd.atomic_grids[0];
    const vec2 &g1 = gd.atomic_grids[1];
    ASSERT_EQ(g0[GridData::Z].size(), (size_t)n);
    EXPECT_NEAR(g0[GridData::Z][0], -kH2Half - 0.5, 1e-12);
    EXPECT_NEAR(g0[GridData::Z][n - 1], kH2Half + 0.5, 1e-12);
    Thakkar h(1);
    for (int p = 0; p < n; p++) {
        EXPECT_EQ(g0[GridData::X][p], 0.0);
        EXPECT_EQ(g0[GridData::Y][p], 0.0);
        EXPECT_EQ(g0[GridData::WEIGHT][p], 1.0);
        EXPECT_EQ(g1[GridData::Z][p], g0[GridData::Z][p]);
        EXPECT_NEAR(g0[GridData::HIRSH_WEIGHT][p] + g1[GridData::HIRSH_WEIGHT][p], 1.0, 1e-10);
        EXPECT_NEAR(g0[GridData::BECKE_WEIGHT][p] + g1[GridData::BECKE_WEIGHT][p], 1.0, 1e-10);
        EXPECT_NEAR(g0[GridData::TFVC_WEIGHT][p] + g1[GridData::TFVC_WEIGHT][p], 1.0, 1e-10);
        const double zp = g0[GridData::Z][p];
        const double r1 = std::abs(zp + kH2Half), r2 = std::abs(zp - kH2Half);
        const double rho1 = h.get_radial_density(r1), rho2 = h.get_radial_density(r2);
        EXPECT_NEAR(g0[GridData::HIRSH_WEIGHT][p], rho1 / (rho1 + rho2), 1e-6) << "point " << p;
        EXPECT_NEAR(g0[GridData::BECKE_WEIGHT][p], becke_step((r1 - r2) / (2.0 * kH2Half), 3), 1e-12) << "point " << p;
        const double rho = w.compute_dens({ 0.0, 0.0, zp });
        EXPECT_NEAR(g0[GridData::WFN_DENSITY][p], rho, 1e-12 * rho);
    }
    EXPECT_GT(g0[GridData::HIRSH_WEIGHT][0], 0.9);
    EXPECT_GT(g1[GridData::HIRSH_WEIGHT][n - 1], 0.9);
    EXPECT_NEAR(g0[GridData::HIRSH_WEIGHT][n / 2], 0.5, 1e-10);
    EXPECT_NEAR(g0[GridData::BECKE_WEIGHT][n / 2], 0.5, 1e-10);
}

//evaluateFunctionOnGrid returns f times the integration weight per point
TEST(GridManagerTests, EvaluateFunctionOnGridMultipliesWeights)
{
    GridManager gm;
    vec2 grid(4, vec(3));
    grid[GridData::X] = { 0.0, 1.0, 2.0 };
    grid[GridData::Y] = { 1.0, 0.0, 3.0 };
    grid[GridData::Z] = { 2.0, 1.0, 0.0 };
    grid[GridData::WEIGHT] = { 0.5, 2.0, 3.0 };
    const vec r = gm.evaluateFunctionOnGrid(grid, [](double x, double y, double z) { return x + 2.0 * y + 3.0 * z; });
    ASSERT_EQ(r.size(), 3u);
    EXPECT_NEAR(r[0], 0.5 * 8.0, 1e-14);
    EXPECT_NEAR(r[1], 2.0 * 4.0, 1e-14);
    EXPECT_NEAR(r[2], 3.0 * 8.0, 1e-14);
}

//writeSimpleGrid writes a header plus one line per point that parses back to the input
TEST(GridIoTests, WriteSimpleGridRoundTrip)
{
    GridManager gm;
    vec2 grid(3, vec(2));
    grid[0] = { 0.25, -1.5 };
    grid[1] = { 1.0, 2.0 };
    grid[2] = { -0.125, 3.0 };
    const vec f = { 1.25e-3, -7.5e2 };
    const std::filesystem::path path = std::filesystem::temp_directory_path() / "GridIoTests_WriteSimpleGridRoundTrip.dat";
    gm.writeSimpleGrid(path, grid, { { "f", f } });
    std::ifstream in(path);
    ASSERT_TRUE(in.is_open());
    std::string line;
    ASSERT_TRUE(std::getline(in, line));
    EXPECT_EQ(line, "# X\tY\tZ\tf\t");
    for (int p = 0; p < 2; p++) {
        ASSERT_TRUE(std::getline(in, line));
        std::istringstream row(line);
        double x, y, z, v;
        row >> x >> y >> z >> v;
        EXPECT_NEAR(x, grid[0][p], 1e-3);
        EXPECT_NEAR(y, grid[1][p], 1e-3);
        EXPECT_NEAR(z, grid[2][p], 1e-3);
        EXPECT_NEAR(v, f[p], 1e-9 * std::abs(f[p]));
    }
    EXPECT_FALSE(std::getline(in, line));
    in.close();
    std::filesystem::remove(path);
}

//a cube holding exp(-r^2) integrated over the Becke grids of H2 gives pi^(3/2), with
//the per-atom electron sums split evenly by symmetry
TEST(GridManagerTests, DensityVectorsFromCubeIntegrateGaussian)
{
    const int n = 65;
    const double h = 0.25, origin = -8.0;
    cube c({ n, n, n }, 0, true);
    for (int i = 0; i < 3; i++) {
        c.set_origin(i, origin);
        for (int j = 0; j < 3; j++)
            c.set_vector(i, j, i == j ? h : 0.0);
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++) {
                const double x = origin + i * h, y = origin + j * h, z = origin + k * h;
                c.set_value(i, j, k, std::exp(-(x * x + y * y + z * z)));
            }
    const WFN w = make_h2();
    GridConfiguration cfg;
    cfg.accuracy = 3;
    cfg.partition_type = PartitionType::Becke;
    GridManager gm(cfg);
    std::ostringstream log;
    gm.setup3DGridsForMolecule(w, { 0, 1 }, {}, cell(), false, log);
    vec2 d1, d2, d3, dens;
    vec electrons;
    gm.getDensityVectorsFromCube(w, { 0, 1 }, c, d1, d2, d3, dens, electrons);
    ASSERT_EQ(electrons.size(), 2u);
    const double ref = std::pow(constants::PI, 1.5);
    EXPECT_NEAR(electrons[0] + electrons[1], ref, 5e-2 * ref);
    EXPECT_NEAR(electrons[0], electrons[1], 1e-6);
    EXPECT_NEAR(vec_total(dens), electrons[0] + electrons[1], 1e-8);
    EXPECT_EQ(gm.getTotalGridPoints(), (int)(dens[0].size() + dens[1].size()));
}

//periodic images of a 100 A cell are too far away to change the Hirshfeld partition
TEST(GridManagerTests, PbcLargeCellMatchesNoPbc)
{
    const WFN w = make_h2();
    GridConfiguration cfg;
    cfg.accuracy = 1;
    GridManager plain(cfg);
    cfg.pbc = 1;
    GridManager periodic(cfg);
    std::ostringstream log;
    const cell big(100.0, 100.0, 100.0, 90.0, 90.0, 90.0);
    plain.setup3DGridsForMolecule(w, { 0, 1 }, {}, cell(), false, log);
    periodic.setup3DGridsForMolecule(w, { 0, 1 }, {}, big, false, log);
    const PartitionResults a = plain.calculatePartitionedCharges(w);
    const PartitionResults b = periodic.calculatePartitionedCharges(w, big);
    EXPECT_EQ(plain.getTotalGridPoints(), periodic.getTotalGridPoints());
    for (int i = 0; i < 2; i++)
        EXPECT_NEAR(a.atom_charges[PartitionResults::S_HIRSH][i], b.atom_charges[PartitionResults::S_HIRSH][i], 1e-8);
    EXPECT_NEAR(a.overall_charges[PartitionResults::S_HIRSH], h2_electrons(), 5e-3);
}

//print_interaction_energy lists the energy terms, the per-pair table with row sums and
//the per-rank table; the KS variant drops the kinetic and exchange rows
TEST(GridIntegratorTests, PrintInteractionEnergyListsPairsAndRanks)
{
    WFN A(e_origin::NOT_YET_DEFINED), B(e_origin::NOT_YET_DEFINED);
    A.push_back_atom("C1", 0.0, 0.0, 0.0, 6);
    A.push_back_atom("O1", 2.0, 0.0, 0.0, 8);
    B.push_back_atom("N1", 0.0, 5.0, 0.0, 7);
    B.push_back_atom("H1", 0.0, 7.0, 0.0, 1);
    DensityFitting::Interaction_Energy E;
    E.nuc_nuc = 1.0;
    E.nucA_rhoB = -0.5;
    E.nucB_rhoA = -0.25;
    E.rho_rho = 0.125;
    E.pol_A = -0.01;
    E.pol_B = -0.02;
    E.disp = -0.03;
    E.rep = 0.04;
    E.rep_kin = 0.05;
    E.rep_x = -0.01;
    E.overlap = 1e-3;
    E.n_A = 14.0;
    E.n_B = 8.0;
    E.x_fun = 1;
    E.pair = { { 0.1, 0.2 }, { 0.3, 0.4 } };
    E.rank = { { 1.0, -0.5 }, { 0.25, 0.125 } };
    EXPECT_NEAR(E.electrostatic(), 0.375, 1e-15);
    EXPECT_NEAR(E.total(), 0.375 - 0.01 - 0.02 - 0.03 + 0.04, 1e-15);

    std::ostringstream out;
    DensityFitting::print_interaction_energy(E, A, B, out);
    const std::string s = out.str();
    EXPECT_NE(s.find("By atom pair"), std::string::npos);
    EXPECT_NE(s.find("By rank"), std::string::npos);
    EXPECT_NE(s.find("rep. kin. TF"), std::string::npos);
    EXPECT_NE(s.find("rep. exch. PBE"), std::string::npos);
    EXPECT_NE(s.find("Gordon-Kim"), std::string::npos);
    const double k = constants::kcal_mol_per_hartree;
    std::ostringstream row;
    row << std::fixed << std::setprecision(3) << std::setw(6) << "C1" << std::setw(9) << 0.1 * k << std::setw(9) << 0.2 * k << std::setw(9) << 0.3 * k << "\n";
    EXPECT_NE(s.find(row.str()), std::string::npos) << s;
    std::ostringstream head;
    head << std::setw(9) << "N1" << std::setw(9) << "H1" << "      sum\n";
    EXPECT_NE(s.find(head.str()), std::string::npos) << s;
    std::ostringstream rank;
    rank << std::fixed << std::setprecision(3) << std::setw(6) << "l=0" << std::setw(11) << 0.25 * k << std::setw(11) << 0.125 * k << std::setw(11) << 0.375 * k << "\n";
    EXPECT_NE(s.find(rank.str()), std::string::npos) << s;
    std::ostringstream total;
    total << "total            " << std::fixed << std::setprecision(6) << std::setw(14) << E.total() << " Eh";
    EXPECT_NE(s.find(total.str()), std::string::npos) << s;

    E.n_A = 0.0;
    E.n_B = 0.0;
    std::ostringstream ks;
    DensityFitting::print_interaction_energy(E, A, B, ks);
    EXPECT_EQ(ks.str().find("rep. kin. TF"), std::string::npos);
    EXPECT_EQ(ks.str().find("rep. exch."), std::string::npos);
    EXPECT_NE(ks.str().find("repulsion K*S"), std::string::npos);
}
