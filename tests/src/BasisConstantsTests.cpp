#include "pch.h"

#include "core/convenience.h"
#include "core/constants.h"
#include "core/atoms.h"
#include "core/wfn_class.h"
#include "core/basis_set.h"
#include "core/integration_params.h"
#include <occ/gto/io/json_basis.h>
#include <occ/gto/shell.h>
#include <occ/core/atom.h>
#undef I

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <set>
#include <string>
#include <vector>

namespace
{
    //independent reference for the associated Legendre polynomial P_l^m(x), m >= 0, without the
    //Condon-Shortley phase: P_m^m = (2m-1)!! (1-x^2)^(m/2), then the upward recurrence in l
    double legendre_ref(const int l, const int m, const double x)
    {
        double pmm = 1.0;
        const double s = std::sqrt(1.0 - x * x);
        for (int k = 1; k <= m; k++) pmm *= (2 * k - 1) * s;
        if (l == m) return pmm;
        double pm1 = x * (2 * m + 1) * pmm;
        if (l == m + 1) return pm1;
        double p = 0.0;
        for (int ll = m + 2; ll <= l; ll++)
        {
            p = ((2 * ll - 1) * x * pm1 - (ll + m - 1) * pmm) / (ll - m);
            pmm = pm1;
            pm1 = p;
        }
        return p;
    }

    double fact(const int n)
    {
        double f = 1.0;
        for (int i = 2; i <= n; i++) f *= i;
        return f;
    }

    //real spherical harmonic N_l|m| P_l^|m|(cos theta) {cos, sin}(|m| phi) of a unit vector,
    //the convention of constants::real_spherical, written out independently
    double ylm_ref(const int l, const int m, const double x, const double y, const double z)
    {
        const int am = std::abs(m);
        const double theta = std::acos(z), phi = std::atan2(y, x);
        const double pi = 3.14159265358979323846;
        const double norm = am == 0 ? std::sqrt((2 * l + 1) / (4 * pi)) : std::sqrt((2 * l + 1) * fact(l - am) / (2 * pi * fact(l + am)));
        return norm * legendre_ref(l, am, std::cos(theta)) * (m >= 0 ? std::cos(am * phi) : std::sin(am * phi));
    }

    //three unit vectors: generic, all-negative octant, close to the z pole
    const std::vector<std::array<double, 3>>& unit_points()
    {
        static const std::vector<std::array<double, 3>> pts = [] {
            std::vector<std::array<double, 3>> raw{ { 0.3, -0.5, 0.8 }, { -0.7, -0.2, -0.4 }, { 0.05, 0.02, 0.9 } };
            for (auto& p : raw)
            {
                const double r = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
                for (double& c : p) c /= r;
            }
            return raw;
        }();
        return pts;
    }

    //STO-3G hydrogen as the library stores it
    const double sto3g_h_exp[3] = { 3.4252509139999998, 0.6239137298, 0.168855404 };
    const double sto3g_h_coef[3] = { 0.1543289673, 0.5353281423, 0.4446345422 };

    std::filesystem::path temp_file(const std::string& name, const std::string& content)
    {
        const std::filesystem::path p = std::filesystem::temp_directory_path() / name;
        std::ofstream(p) << content;
        return p;
    }
}

//the l <= 8 template with unit coefficients against an independent Legendre recurrence
TEST(BasisConstantsTests, TemplateSphericalHarmonicMatchesReference)
{
    for (const auto& p : unit_points())
        for (int l = 0; l <= 8; l++)
            for (int m = -l; m <= l; m++)
            {
                double coefs[17] = { 0.0 };
                coefs[m + l] = 1.0;
                EXPECT_NEAR(constants::spherical_harmonic<double>(l, p[0], p[1], p[2], coefs), ylm_ref(l, m, p[0], p[1], p[2]), 1e-10) << "l " << l << " m " << m;
            }
}

//the template is linear in the coefficients and returns 0 above l = 8
TEST(BasisConstantsTests, TemplateSphericalHarmonicIsLinearInCoefficients)
{
    const auto& p = unit_points()[0];
    for (int l = 0; l <= 8; l++)
    {
        double coefs[17] = { 0.0 };
        double expected = 0.0;
        for (int m = -l; m <= l; m++)
        {
            coefs[m + l] = 0.1 * (m + 3) * (l % 2 ? -1.0 : 1.0);
            expected += coefs[m + l] * ylm_ref(l, m, p[0], p[1], p[2]);
        }
        EXPECT_NEAR(constants::spherical_harmonic<double>(l, p[0], p[1], p[2], coefs), expected, 1e-10) << "l " << l;
    }
    double ones[19];
    std::fill(ones, ones + 19, 1.0);
    EXPECT_EQ(constants::spherical_harmonic<double>(9, p[0], p[1], p[2], ones), 0.0);
}

//the hand-written (l, m, d) switch up to l = 9 against the reference, so both paths of
//spherical_harmonic(l, d, coefs) share one convention
TEST(BasisConstantsTests, SwitchSphericalHarmonicMatchesReferenceUpToL9)
{
    for (const auto& p : unit_points())
    {
        const double d[3] = { p[0], p[1], p[2] };
        for (int l = 0; l <= 9; l++)
            for (int m = -l; m <= l; m++)
                EXPECT_NEAR(constants::spherical_harmonic(l, m, d), ylm_ref(l, m, p[0], p[1], p[2]), 1e-10) << "l " << l << " m " << m;
    }
}

//spherical_harmonic(l, d, coefs) forwards to the template for l <= 8 and sums the switch above,
//so a contraction at l = 9 equals the coefficient-weighted (9, m, d) values
TEST(BasisConstantsTests, CollapsedSphericalHarmonicSumsAboveL8)
{
    const auto& p = unit_points()[1];
    const double d[3] = { p[0], p[1], p[2] };
    double coefs[19];
    for (int i = 0; i < 19; i++) coefs[i] = std::sin(0.7 * i + 0.1);
    double expected = 0.0;
    for (int m = -9; m <= 9; m++) expected += coefs[m + 9] * ylm_ref(9, m, p[0], p[1], p[2]);
    EXPECT_NEAR(constants::spherical_harmonic(9, d, coefs), expected, 1e-10);
    EXPECT_NEAR(constants::spherical_harmonic(4, d, coefs), constants::spherical_harmonic<double>(4, p[0], p[1], p[2], coefs), 1e-14);
}

//the default branch of the switch (l > 9) goes through cartesian_to_spherical and real_spherical;
//the Apple build has no Legendre polynomials above l = 8 on that path
TEST(BasisConstantsTests, SphericalHarmonicDefaultBranchUsesRealSpherical)
{
#ifdef __APPLE__
    GTEST_SKIP() << "associated_legendre_polynomial above l = 8 is not implemented on Apple";
#else
    for (const auto& p : unit_points())
    {
        const double d[3] = { p[0], p[1], p[2] };
        for (int l = 10; l <= 12; l++)
            for (int m = -l; m <= l; m++)
            {
                EXPECT_NEAR(constants::spherical_harmonic(l, m, d), ylm_ref(l, m, p[0], p[1], p[2]), 1e-10) << "l " << l << " m " << m;
                EXPECT_NEAR(constants::real_spherical(l, m, std::acos(p[2]), std::atan2(p[1], p[0])), ylm_ref(l, m, p[0], p[1], p[2]), 1e-10) << "l " << l << " m " << m;
            }
    }
#endif
}

//addition theorem: sum_m Y_lm^2 = (2l+1)/(4 pi) at every point, the normalisation of each l
TEST(BasisConstantsTests, SphericalHarmonicsSatisfyAdditionTheorem)
{
    for (const auto& p : unit_points())
        for (int l = 0; l <= 8; l++)
        {
            double sum = 0.0;
            for (int m = -l; m <= l; m++)
            {
                double coefs[17] = { 0.0 };
                coefs[m + l] = 1.0;
                const double y = constants::spherical_harmonic<double>(l, p[0], p[1], p[2], coefs);
                sum += y * y;
            }
            EXPECT_NEAR(sum, (2 * l + 1) / (4 * constants::PI), 1e-10) << "l " << l;
        }
}

//the manual Legendre polynomials for m >= 0 against the recurrence at three arguments
TEST(BasisConstantsTests, AssociatedLegendreNonNegativeMMatchesRecurrence)
{
    for (const double x : { -0.7, 0.13, 0.9 })
        for (int l = 0; l <= 8; l++)
            for (int m = 0; m <= l; m++)
                EXPECT_NEAR(constants::associated_legendre_polynomial(l, m, x), legendre_ref(l, m, x), 1e-9 * std::max(1.0, std::abs(legendre_ref(l, m, x)))) << "l " << l << " m " << m << " x " << x;
}

//the manual negative-m branches all follow P_l^-m = -(l-m)!/(l+m)! P_l^m, the convention the
//l <= 8 cases were written in (the sign differs from the textbook (-1)^m for even m); no caller
//reaches m < 0 (real_spherical passes |m|), so this pins the dead branch, not a physics claim
TEST(BasisConstantsTests, AssociatedLegendreNegativeMFollowsFactorialRatio)
{
    for (const double x : { -0.7, 0.13, 0.9 })
        for (int l = 1; l <= 8; l++)
            for (int m = 1; m <= l; m++)
            {
                const double expected = -fact(l - m) / fact(l + m) * legendre_ref(l, m, x);
                EXPECT_NEAR(constants::associated_legendre_polynomial(l, -m, x), expected, 1e-9 * std::max(1.0, std::abs(expected))) << "l " << l << " m " << -m << " x " << x;
            }
}

//above l = 8 the function falls back to std::assoc_legendre, which must agree with the recurrence
TEST(BasisConstantsTests, AssociatedLegendreDefaultBranchAboveL8)
{
#ifdef __APPLE__
    GTEST_SKIP() << "associated_legendre_polynomial above l = 8 is not implemented on Apple";
#else
    for (const double x : { -0.7, 0.13, 0.9 })
        for (int l = 9; l <= 12; l++)
            for (int m = 0; m <= l; m++)
                EXPECT_NEAR(constants::associated_legendre_polynomial(l, m, x), legendre_ref(l, m, x), 1e-9 * std::max(1.0, std::abs(legendre_ref(l, m, x)))) << "l " << l << " m " << m;
#endif
}

//the constexpr norm table is sqrt((2l+1)(l-m)!/(2 pi (l+m)!)) at column l + m with the m = 0 column at sqrt((2l+1)/4pi);
//only the m >= 0 columns are read (spherical_harmonic indexes l + |m|), the m < 0 columns hold the inverse ratio
TEST(BasisConstantsTests, SphericalNormsTableMatchesFormula)
{
    EXPECT_NEAR(constants::spherical_norms[0][0], std::sqrt(1.0 / (4 * constants::PI)), 1e-15);
    for (int l = 1; l <= constants::ASSOCIATED_LEGENDRE_MAX_L; l++)
    {
        EXPECT_NEAR(constants::spherical_norms[l][l], std::sqrt((2 * l + 1) / (4 * constants::PI)), 1e-14) << "l " << l;
        for (int m = 1; m <= l; m++)
        {
            const double expected = std::sqrt((2 * l + 1) * fact(l - m) / (2 * constants::PI * fact(l + m)));
            EXPECT_NEAR(constants::spherical_norms[l][l + m], expected, 1e-14 * expected) << "l " << l << " m " << m;
        }
    }
}

//cartesian -> spherical -> cartesian round trip, the origin special case, and the normalised
//variant agreeing with the full one on unit vectors
TEST(BasisConstantsTests, CartesianToSphericalRoundTrip)
{
    const std::vector<std::array<double, 3>> pts{ { 1.0, -2.0, 0.5 }, { -0.3, 0.1, -4.0 }, { 0.0, 0.0, 2.0 }, { 0.0, 3.0, 0.0 } };
    for (const auto& p : pts)
    {
        const vec s = constants::cartesian_to_spherical(p[0], p[1], p[2]);
        ASSERT_EQ(s.size(), 3u);
        const double r = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
        EXPECT_NEAR(s[0], r, 1e-14);
        EXPECT_NEAR(r * std::sin(s[1]) * std::cos(s[2]), p[0], 1e-13);
        EXPECT_NEAR(r * std::sin(s[1]) * std::sin(s[2]), p[1], 1e-13);
        EXPECT_NEAR(r * std::cos(s[1]), p[2], 1e-13);
        const std::pair<double, double> n = constants::norm_cartesian_to_spherical(p[0] / r, p[1] / r, p[2] / r);
        EXPECT_NEAR(n.first, s[1], 1e-14);
        EXPECT_NEAR(n.second, s[2], 1e-14);
    }
    const vec origin = constants::cartesian_to_spherical(0.0, 0.0, 0.0);
    EXPECT_EQ(origin[0], 0.0);
    EXPECT_EQ(origin[1], 0.0);
    EXPECT_EQ(origin[2], 0.0);
}

//constexpr sqrt against std::sqrt over eight decades, and NaN for the two documented rejects
TEST(BasisConstantsTests, ConstexprSqrtMatchesStdSqrt)
{
    for (const double x : { 1e-8, 0.5, 2.0, 12345.678, 1e6, 3.0e15 })
        EXPECT_NEAR(constants::sqrt(x), std::sqrt(x), 1e-14 * std::sqrt(x)) << x;
    static_assert(constants::sqrt(16.0) == 4.0, "sqrt of a perfect square is exact");
#if !defined(__FAST_MATH__)
    EXPECT_TRUE(std::isnan(constants::sqrt(-2.0)));
    EXPECT_TRUE(std::isnan(constants::sqrt(std::numeric_limits<double>::infinity())));
#endif
}

//the integer, double and double-factorial tables agree with each other and with ft_fun
TEST(BasisConstantsTests, FactorialTablesAreConsistent)
{
    for (int n = 0; n <= 20; n++)
    {
        EXPECT_EQ(static_cast<double>(constants::ft[n]), constants::ftd[n]) << n;
        EXPECT_EQ(static_cast<double>(constants::ft_fun(n)), constants::ftd[n]) << n;
    }
    for (int n = 1; n < constants::MAX_FACTORIAL; n++)
        EXPECT_EQ(constants::ftd[n], constants::ftd[n - 1] * n) << n;
    EXPECT_NEAR(constants::ftd[29], 8.841761993739701954543616e30, 1e17);
    for (int n = 2; n < 25; n++)
        EXPECT_EQ(constants::double_ft[n], n * constants::double_ft[n - 2]) << n;
    //(2n)! = 2^n n! (2n-1)!! ties the two tables together
    for (int n = 1; n <= 10; n++)
        EXPECT_EQ(constants::ftd[2 * n], std::pow(2.0, n) * constants::ftd[n] * static_cast<double>(constants::double_ft[2 * n - 1])) << n;
}

//the ORCA (0, +1, -1, +2, -2, ...) to PySCF (-l..l) map is idx -> l + m for every l up to 10,
//nullopt beyond
TEST(BasisConstantsTests, OrcaToPyscfMapIsShiftedMForEveryL)
{
    for (unsigned int l = 0; l <= 10; l++)
    {
        std::set<std::size_t> seen;
        for (int idx = 0; idx < 2 * static_cast<int>(l) + 1; idx++)
        {
            const int m = idx == 0 ? 0 : (idx % 2 ? (idx + 1) / 2 : -idx / 2);
            const std::optional<std::size_t> r = constants::orca_2_pySCF(l, idx);
            ASSERT_TRUE(r.has_value()) << "l " << l << " idx " << idx;
            EXPECT_EQ(static_cast<int>(*r), static_cast<int>(l) + m) << "l " << l << " idx " << idx;
            seen.insert(*r);
        }
        EXPECT_EQ(seen.size(), 2 * l + 1) << "l " << l;
    }
    EXPECT_FALSE(constants::orca_2_pySCF(11, 0).has_value());
}

//type2vector enumerates the (l+1)(l+2)/2 Cartesian components of each l in order, 286 types up
//to l = 10, with the -1 triple outside that range; types 1..20 in the order the WFN format
//prescribes (S, PX PY PZ, DXX DYY DZZ DXY DXZ DYZ, FXXX FYYY FZZZ FXXY FXXZ FYYZ FXYY FXZZ FYZZ FXYZ)
TEST(BasisConstantsTests, TypeToVectorEnumeratesCartesianComponents)
{
    const int wfn_order[20][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 },
        { 2, 0, 0 }, { 0, 2, 0 }, { 0, 0, 2 }, { 1, 1, 0 }, { 1, 0, 1 }, { 0, 1, 1 },
        { 3, 0, 0 }, { 0, 3, 0 }, { 0, 0, 3 }, { 2, 1, 0 }, { 2, 0, 1 }, { 0, 2, 1 }, { 1, 2, 0 }, { 1, 0, 2 }, { 0, 1, 2 }, { 1, 1, 1 } };
    for (int t = 1; t <= 20; t++)
    {
        int v[3];
        constants::type2vector(t, v);
        EXPECT_EQ(v[0], wfn_order[t - 1][0]) << t;
        EXPECT_EQ(v[1], wfn_order[t - 1][1]) << t;
        EXPECT_EQ(v[2], wfn_order[t - 1][2]) << t;
    }
    int index = 1;
    for (int l = 0; l <= 10; l++)
    {
        std::set<int> seen;
        for (int k = 0; k < (l + 1) * (l + 2) / 2; k++, index++)
        {
            int v[3];
            constants::type2vector(index, v);
            EXPECT_GE(std::min({ v[0], v[1], v[2] }), 0) << index;
            EXPECT_EQ(v[0] + v[1] + v[2], l) << index;
            seen.insert(v[0] * 121 + v[1] * 11 + v[2]);
        }
        EXPECT_EQ(static_cast<int>(seen.size()), (l + 1) * (l + 2) / 2) << "l " << l;
    }
    EXPECT_EQ(index, 287);
    for (const int bad : { 0, -5, 287, 1000 })
    {
        int v[3] = { 7, 7, 7 };
        constants::type2vector(bad, v);
        EXPECT_EQ(v[0], -1) << bad;
        EXPECT_EQ(v[1], -1) << bad;
        EXPECT_EQ(v[2], -1) << bad;
    }
}

//get_Z_from_label(Labels[Z]) == Z - 1 from hydrogen to polonium, deuterium and tritium fold onto
//hydrogen, anything unknown is -1
TEST(BasisConstantsTests, ZFromLabelRoundTripsThroughLabelsUpToPo)
{
    for (int z = 1; z <= 84; z++)
        EXPECT_EQ(constants::get_Z_from_label(constants::Labels[z]), z - 1) << constants::Labels[z];
    EXPECT_EQ(constants::get_Z_from_label("D"), 0);
    EXPECT_EQ(constants::get_Z_from_label("T"), 0);
    EXPECT_EQ(constants::get_Z_from_label("Xx"), -1);
    EXPECT_EQ(constants::get_Z_from_label(""), -1);
    EXPECT_EQ(constants::get_Z_from_label("Hg1"), -1);
    static_assert(constants::get_Z_from_label("Pd") == 45, "palladium is Z = 46");
}

//get_Z_from_label returns Z - 1 through the actinides and up to Rg
TEST(BasisConstantsTests, ZFromLabelAtToRgIsOffByOne)
{
    for (int z = 85; z <= 103; z++)
        EXPECT_EQ(constants::get_Z_from_label(constants::Labels[z]), z - 1) << constants::Labels[z];
    EXPECT_EQ(constants::get_Z_from_label("Rf"), 103);
    EXPECT_EQ(constants::get_Z_from_label("Rg"), 110);
}

//atnr2letter is the inverse table for 1..103 plus the two Q-peak conventions
TEST(BasisConstantsTests, AtnrToLetterMatchesLabels)
{
    for (int z = 1; z <= 103; z++)
        EXPECT_EQ(std::string(constants::atnr2letter(z)), std::string(constants::Labels[z])) << z;
    EXPECT_EQ(std::string(constants::atnr2letter(0)), "Q");
    EXPECT_EQ(std::string(constants::atnr2letter(119)), "Q");
    EXPECT_EQ(std::string(constants::atnr2letter(104)), "PROBLEM");
}

//get_closest_num_angular is the smallest Lebedev size >= n and get_angular_order its index
TEST(BasisConstantsTests, LebedevLookupsRoundTrip)
{
    const int last = 31; //lebedev_table[32] is the zero-filled tail of the 33-entry array
    for (int i = 0; i <= last; i++)
        EXPECT_EQ(constants::get_angular_order(constants::lebedev_table[i]), i) << i;
    for (const int n : { 1, 6, 7, 50, 51, 302, 303, 5810 })
    {
        const int* it = std::lower_bound(constants::lebedev_table, constants::lebedev_table + last + 1, n);
        EXPECT_EQ(constants::get_closest_num_angular(n), *it) << n;
    }
    EXPECT_EQ(constants::get_closest_num_angular(5811), -1);
    EXPECT_EQ(constants::get_angular_order(7), -1);
    EXPECT_EQ(constants::get_angular_order(5811), -1);
}

//the power variants of the unit conversions are the plain ones applied p times and invert each other
TEST(BasisConstantsTests, UnitConversionPowersRoundTrip)
{
    const double x = 2.75;
    EXPECT_NEAR(constants::bohr2ang_p(x, 3), constants::cubic_bohr2ang(x), 1e-14);
    EXPECT_NEAR(constants::ang2bohr_p(x, 3), constants::cubic_ang2bohr(x), 1e-14);
    EXPECT_NEAR(constants::bohr2ang_p(x, 2), constants::bohr2ang(constants::bohr2ang(x)), 1e-15);
    EXPECT_NEAR(constants::ang2bohr_p(constants::bohr2ang_p(x, 4), 4), x, 1e-13);
    EXPECT_NEAR(constants::bohr2ang(constants::ang2bohr(x)), x, 1e-14);
}

//name lookup: case and underscores are normalised, an orbital set is found by its bare name, a
//unique suffix match is accepted, several substring matches are refused rather than guessed
TEST(BasisConstantsLibraryTests, BasisSetNamesResolveByRule)
{
    EXPECT_EQ(BasisSetLibrary::get_basis_set("STO_3G")->get_name(), "sto-3g-basis");
    EXPECT_EQ(BasisSetLibrary::get_basis_set("def2-TZVP")->get_name(), "def2-tzvp-basis");
    EXPECT_EQ(BasisSetLibrary::get_basis_set("def2-svpd")->get_name(), "def2-svpd-rifit");
    EXPECT_EQ(BasisSetLibrary::get_basis_set("x2c-jfit")->get_name(), "x2c-jfit");
    EXPECT_TRUE(BasisSetLibrary::check_basis_set_exists("CC_PVTZ"));
    EXPECT_TRUE(BasisSetLibrary::check_basis_set_exists("def2-universal-jkfit"));
    EXPECT_FALSE(BasisSetLibrary::check_basis_set_exists("def2-tzv"));
    EXPECT_FALSE(BasisSetLibrary::check_basis_set_exists("no-such-basis"));
}

//STO-3G hydrogen and carbon as stored: exponents, contraction coefficients, 0-based l and shells
TEST(BasisConstantsLibraryTests, Sto3gPrimitivesAndElementRange)
{
    std::shared_ptr<BasisSet> b = BasisSetLibrary::get_basis_set("sto-3g");
    const std::span<const SimplePrimitive> h = (*b)[0];
    ASSERT_EQ(h.size(), 3u);
    for (int i = 0; i < 3; i++)
    {
        EXPECT_NEAR(h[i].exp, sto3g_h_exp[i], 1e-12);
        EXPECT_NEAR(h[i].coefficient, sto3g_h_coef[i], 1e-12);
        EXPECT_EQ(h[i].type, 0);
        EXPECT_EQ(h[i].shell, 0);
    }
    const std::span<const SimplePrimitive> c = (*b)[5];
    ASSERT_EQ(c.size(), 9u);
    EXPECT_NEAR(c[0].exp, 71.61683735, 1e-8);
    EXPECT_NEAR(c[3].exp, 2.941249355, 1e-9);
    EXPECT_NEAR(c[3].coefficient, -0.09996722919, 1e-11);
    for (int i = 0; i < 9; i++)
    {
        EXPECT_EQ(c[i].shell, i / 3) << i;
        EXPECT_EQ(c[i].type, i < 6 ? 0 : 1) << i;
    }
    EXPECT_NEAR(c[6].exp, c[3].exp, 1e-15); //the 2s and 2p of STO-3G share their exponents
    EXPECT_TRUE(b->has_element(1));
    EXPECT_TRUE(b->has_element(54));
    EXPECT_FALSE(b->has_element(55));
    EXPECT_FALSE(b->has_element(118));
    //3 per s shell: H-He 3, Li-Ne 9 (1s 2sp), Na-Ar 15, K-Ca 21, Sc-Kr 24 (+3d), Rb-Sr 30, Y-Xe 33
    EXPECT_EQ(b->get_primitive_count(), 2u * 3 + 8 * 9 + 8 * 15 + 2 * 21 + 16 * 24 + 2 * 30 + 16 * 33);
}

//get_data converts every element to primitive objects, empty for elements the set lacks, and
//hands out the same converted table on the second call
TEST(BasisConstantsLibraryTests, GetDataConvertsToPrimitives)
{
    std::shared_ptr<BasisSet> b = BasisSetLibrary::get_basis_set("sto-3g");
    std::shared_ptr<std::array<std::vector<primitive>, 118>> data = b->get_data();
    ASSERT_EQ((*data)[0].size(), 3u);
    EXPECT_NEAR((*data)[0][1].get_exp(), sto3g_h_exp[1], 1e-12);
    EXPECT_NEAR((*data)[0][1].get_coef(), sto3g_h_coef[1], 1e-12);
    EXPECT_EQ((*data)[0][1].get_type(), 0);
    EXPECT_EQ((*data)[5].size(), 9u);
    EXPECT_EQ((*data)[5][8].get_type(), 1);
    EXPECT_EQ((*data)[117].size(), 0u);
    std::shared_ptr<std::array<std::vector<primitive>, 118>> again = b->get_data();
    EXPECT_EQ((*again)[5].size(), 9u);
    EXPECT_TRUE((*again)[0][2] == (*data)[0][2]);
}

//operator+= keeps this set's elements, takes the other's for elements this one lacks, and
//rebuilds the owned storage with fresh offsets
TEST(BasisConstantsLibraryTests, OwnedBasisSetsMergeWithPlusEquals)
{
    BasisSet a;
    a.set_name("a");
    a.set_count_for_element(0, 2);
    a.add_owned_primitive({ 0, 0, 2.0, 1.0, 0 });
    a.add_owned_primitive({ 0, 1, 0.5, 1.0, 1 });
    BasisSet b;
    b.set_name("b");
    b.set_count_for_element(0, 1);
    b.add_owned_primitive({ 0, 0, 9.0, 1.0, 0 });
    b.set_count_for_element(1, 1);
    b.add_owned_primitive({ 0, 2, 0.7, 1.0, 0 });
    EXPECT_EQ(b[1].size(), 1u);
    EXPECT_NEAR(b[1][0].exp, 0.7, 1e-15);

    a += b;
    EXPECT_EQ(a.get_name(), "a_plus_b");
    EXPECT_EQ(a.get_primitive_count(), 3u);
    EXPECT_EQ(a.get_owned_primitive_count(), 3u);
    ASSERT_EQ(a[0].size(), 2u);
    EXPECT_NEAR(a[0][0].exp, 2.0, 1e-15);
    EXPECT_NEAR(a[0][1].exp, 0.5, 1e-15);
    ASSERT_EQ(a[1].size(), 1u);
    EXPECT_NEAR(a[1][0].exp, 0.7, 1e-15);
    EXPECT_EQ(a[1][0].type, 2);
    EXPECT_TRUE(a.has_element(1));
    EXPECT_TRUE(a.has_element(2));
    EXPECT_FALSE(a.has_element(3));
}

//write_occ_json is read back by OCC's own JSON basis reader with the same shells, exponents and
//coefficients per element
TEST(BasisConstantsLibraryTests, WriteOccJsonRoundTripsThroughOccReader)
{
    std::shared_ptr<BasisSet> b = BasisSetLibrary::get_basis_set("sto-3g");
    const std::filesystem::path path = std::filesystem::temp_directory_path() / "BasisConstants_sto3g.json";
    b->write_occ_json(path, { 1, 6 });
    occ::gto::io::JsonBasisReader reader(path.string());
    std::filesystem::remove(path);
    EXPECT_EQ(reader.element_map().size(), 2u);
    for (const int z : { 1, 6 })
    {
        const occ::gto::io::ElementBasis& eb = reader.element_basis(z);
        const std::span<const SimplePrimitive> prims = (*b)[z - 1];
        ASSERT_EQ(static_cast<int>(eb.electron_shells.size()), prims.back().shell + 1) << z;
        int p = 0;
        for (const occ::gto::io::ElectronShell& sh : eb.electron_shells)
        {
            ASSERT_EQ(sh.angular_momentum.size(), 1u);
            EXPECT_EQ(sh.angular_momentum[0], prims[p].type);
            EXPECT_EQ(sh.function_type, "gto");
            ASSERT_EQ(sh.coefficients.size(), 1u);
            ASSERT_EQ(sh.exponents.size(), sh.coefficients[0].size());
            for (size_t i = 0; i < sh.exponents.size(); i++, p++)
            {
                EXPECT_NEAR(sh.exponents[i], prims[p].exp, 1e-14 * prims[p].exp);
                EXPECT_NEAR(sh.coefficients[0][i], prims[p].coefficient, 1e-14);
            }
        }
        EXPECT_EQ(p, static_cast<int>(prims.size())) << z;
    }
}

//to_AOBasis builds one spherical shell per library shell: H2 is two s functions, C adds 2s and 2p
TEST(BasisConstantsLibraryTests, ToAOBasisCountsShellsAndFunctions)
{
    std::shared_ptr<BasisSet> b = BasisSetLibrary::get_basis_set("sto-3g");
    const std::vector<occ::core::Atom> h2{ { 1, 0.0, 0.0, 0.0 }, { 1, 0.0, 0.0, 1.4 } };
    const occ::qm::AOBasis ao_h2 = b->to_AOBasis(h2);
    EXPECT_EQ(ao_h2.nbf(), 2u);
    EXPECT_EQ(ao_h2.size(), 2u);
    const std::vector<occ::core::Atom> ch{ { 6, 0.0, 0.0, 0.0 }, { 1, 0.0, 0.0, 2.0 } };
    const occ::qm::AOBasis ao_ch = b->to_AOBasis(ch);
    EXPECT_EQ(ao_ch.nbf(), 6u);
    EXPECT_EQ(ao_ch.size(), 4u);
    EXPECT_EQ(static_cast<int>(ao_ch.l_max()), 1);
    EXPECT_EQ(ao_ch.shells()[2].l, 1);
    EXPECT_EQ(ao_ch.shells()[2].num_primitives(), 3u);
    EXPECT_NEAR(ao_ch.shells()[2].exponents(0), (*b)[5][6].exp, 1e-15);
    EXPECT_EQ(ao_ch.shells()[3].l, 0);
    EXPECT_NEAR(ao_ch.shells()[3].exponents(2), sto3g_h_exp[2], 1e-12);
    EXPECT_NEAR(ao_ch.shells()[3].origin[2], 2.0, 1e-15);
}

//load_basis_into_WFN: decontracted every primitive is its own shell with coefficient 1, contracted
//the library shells and coefficients survive, and the returned count is the number of spherical functions
TEST(BasisConstantsLibraryTests, LoadBasisDecontractedVersusContracted)
{
    std::shared_ptr<BasisSet> b = BasisSetLibrary::get_basis_set("sto-3g");
    WFN dec(e_origin::NOT_YET_DEFINED);
    dec.push_back_atom("H", 0.0, 0.0, 0.0, 1);
    dec.push_back_atom("C", 0.0, 0.0, 1.5, 6);
    EXPECT_EQ(load_basis_into_WFN(dec, b), 3 + 6 + 9);
    EXPECT_EQ(dec.get_atom_basis_set_size(0), 3);
    EXPECT_EQ(dec.get_atom_basis_set_size(1), 9);
    EXPECT_EQ(dec.get_atom(1).get_shellcount_size(), 9u);
    for (int i = 0; i < 9; i++)
    {
        const basis_set_entry e = dec.get_atom_basis_set_entry(1, i);
        EXPECT_EQ(e.get_coefficient(), 1.0) << i;
        EXPECT_EQ(e.get_shell(), i) << i;
        EXPECT_EQ(e.get_type(), i < 6 ? 0 : 1) << i;
    }
    EXPECT_NEAR(dec.get_atom_basis_set_entry(0, 2).get_exponent(), sto3g_h_exp[2], 1e-12);

    WFN con(e_origin::NOT_YET_DEFINED);
    con.push_back_atom("H", 0.0, 0.0, 0.0, 1);
    con.push_back_atom("C", 0.0, 0.0, 1.5, 6);
    EXPECT_EQ(load_basis_into_WFN(con, b, false), 1 + 1 + 1 + 3);
    EXPECT_EQ(con.get_atom_basis_set_size(1), 9);
    EXPECT_EQ(con.get_atom(1).get_shellcount_size(), 3u);
    EXPECT_EQ(con.get_atom(1).get_shellcount(2), 3u);
    EXPECT_NEAR(con.get_atom_basis_set_entry(0, 1).get_coefficient(), sto3g_h_coef[1], 1e-12);
    EXPECT_EQ(con.get_atom_basis_set_entry(1, 8).get_shell(), 2);
    EXPECT_EQ(con.get_nr_basis_set_loaded(), 2);
}

//complete = true also fills the WFN primitive arrays: one entry per Cartesian component, 1-based
//types (s = 1, p = 2..4) and 1-based centres
TEST(BasisConstantsLibraryTests, LoadBasisCompleteFillsPrimitiveArrays)
{
    std::shared_ptr<BasisSet> b = BasisSetLibrary::get_basis_set("sto-3g");
    WFN w(e_origin::NOT_YET_DEFINED);
    w.push_back_atom("H", 0.0, 0.0, 0.0, 1);
    w.push_back_atom("C", 0.0, 0.0, 1.5, 6);
    load_basis_into_WFN(w, b, false, true);
    ASSERT_EQ(w.get_nex(), 3 + 6 + 3 * 3);
    for (int i = 0; i < 9; i++)
    {
        EXPECT_EQ(w.get_type(i), 1) << i;
        EXPECT_EQ(w.get_center(i), i < 3 ? 1 : 2) << i;
    }
    const std::span<const SimplePrimitive> c = (*b)[5];
    for (int p = 0; p < 3; p++)
        for (int k = 0; k < 3; k++)
        {
            EXPECT_EQ(w.get_type(9 + 3 * p + k), 2 + k) << p << " " << k;
            EXPECT_EQ(w.get_center(9 + 3 * p + k), 2);
            EXPECT_NEAR(w.get_exponent(9 + 3 * p + k), c[6 + p].exp, 1e-15);
        }
    EXPECT_NEAR(w.get_exponent(0), sto3g_h_exp[0], 1e-12);
    EXPECT_NEAR(w.get_exponent(3), c[0].exp, 1e-15);
}

//AutoAux on an s-only hydrogen: one even-tempered s series with beta 1.8 starting at twice the
//smallest orbital exponent, and a second call leaves the generated set alone.
//Series length by hand: the top is min(20 a_eff, 2 a_max). a_eff = 2/(pi r_i^2) with
//r_i = c_i sum_j c_j G(2, e_i + e_j) over the normalised contraction, so sum_i r_i = c^T G c = 1
//and the smallest r_i <= 1/3 (all three STO-3G coefficients are positive), which puts
//20 a_eff >= 20 * 18/pi = 115 above 2 a_max = 6.85. Hence the cap is 2 a_max, the ratio to
//2 a_min is 20.29, log_1.8 of that is 5.12, ceil + 1 = 7 functions with the top at 2 a_min 1.8^6
TEST(BasisConstantsLibraryTests, AutoAuxForHydrogenIsEvenTemperedSeries)
{
    WFN orb(e_origin::NOT_YET_DEFINED);
    orb.push_back_atom("H", 0.0, 0.0, 0.0, 1);
    for (int i = 0; i < 3; i++) orb.push_back_atom_basis_set(0, sto3g_h_exp[i], sto3g_h_coef[i], 1, 0);

    BasisSet aux;
    aux.gen_auto_aux(orb);
    EXPECT_EQ(aux.get_name(), "auto-aux-basis");
    EXPECT_TRUE(aux.has_element(1));
    const std::span<const SimplePrimitive> s = aux[0];
    ASSERT_GE(s.size(), 2u);
    EXPECT_EQ(s.size(), aux.get_owned_primitive_count());
    for (size_t i = 0; i < s.size(); i++)
    {
        EXPECT_EQ(s[i].type, 0) << i;
        EXPECT_EQ(s[i].shell, static_cast<int>(i)) << i;
        EXPECT_EQ(s[i].coefficient, 1.0) << i;
        if (i > 0) EXPECT_NEAR(s[i - 1].exp / s[i].exp, 1.8, 1e-9) << i;
    }
    EXPECT_EQ(s.size(), 7u);
    EXPECT_NEAR(s.back().exp, 2 * sto3g_h_exp[2], 1e-12);
    EXPECT_NEAR(s.front().exp, 2 * sto3g_h_exp[2] * std::pow(1.8, 6), 1e-9);

    const uint32_t count = aux.get_primitive_count();
    aux.gen_auto_aux(orb);
    EXPECT_EQ(aux.get_primitive_count(), count);
}

//with a p shell on hydrogen AutoAux reaches l = 2: the lowest exponent of each aux l is the smallest
//orbital pair sum that couples to it, and the ratios are 1.8 / 2.0 / 2.2
TEST(BasisConstantsLibraryTests, AutoAuxWithPShellCoversLUpToTwo)
{
    WFN orb(e_origin::NOT_YET_DEFINED);
    orb.push_back_atom("H", 0.0, 0.0, 0.0, 1);
    for (int i = 0; i < 3; i++) orb.push_back_atom_basis_set(0, sto3g_h_exp[i], sto3g_h_coef[i], 1, 0);
    orb.push_back_atom_basis_set(0, 0.8, 1.0, 2, 1);

    BasisSet aux;
    aux.gen_auto_aux_for_element(orb.get_atom(0));
    const std::span<const SimplePrimitive> s = aux[0];
    const double beta[3] = { 1.8, 2.0, 2.2 };
    const double a_min[3] = { 2 * sto3g_h_exp[2], sto3g_h_exp[2] + 0.8, 1.6 };
    std::vector<std::vector<double>> per_l(3);
    for (const SimplePrimitive& p : s)
    {
        ASSERT_GE(p.type, 0);
        ASSERT_LE(p.type, 2);
        per_l[p.type].push_back(p.exp);
    }
    EXPECT_EQ(per_l[0].size(), 7u); //the s series of the test above, the p shell does not lower its cap
    for (int l = 0; l < 3; l++)
    {
        ASSERT_FALSE(per_l[l].empty()) << "l " << l;
        EXPECT_NEAR(per_l[l].back(), a_min[l], 1e-12) << "l " << l;
        for (size_t i = 1; i < per_l[l].size(); i++)
            EXPECT_NEAR(per_l[l][i - 1] / per_l[l][i], beta[l], 1e-9) << "l " << l << " i " << i;
    }
    EXPECT_EQ(s.size(), per_l[0].size() + per_l[1].size() + per_l[2].size());
}

//generate_aux_wfn with an empty set restricted to hydrogen: the restriction leaves helium out of
//AutoAux, the completeness pass then generates it, and both atoms end up with 0-based aux shells
TEST(BasisConstantsLibraryTests, GenerateAuxWfnFillsMissingElementsWithAutoAux)
{
    WFN orb(e_origin::NOT_YET_DEFINED);
    orb.push_back_atom("H", 0.0, 0.0, 0.0, 1);
    orb.push_back_atom("He", 0.0, 0.0, 2.0, 2);
    const double he_exp[3] = { 6.362421394, 1.158922999, 0.3136497915 };
    for (int i = 0; i < 3; i++)
    {
        orb.push_back_atom_basis_set(0, sto3g_h_exp[i], sto3g_h_coef[i], 1, 0);
        orb.push_back_atom_basis_set(1, he_exp[i], sto3g_h_coef[i], 1, 0);
    }
    std::shared_ptr<BasisSet> restricted = std::make_shared<BasisSet>();
    restricted->set_auto_aux_elements({ 1 });
    std::vector<std::shared_ptr<BasisSet>> sets{ restricted };
    const WFN aux = generate_aux_wfn(orb, sets);

    EXPECT_EQ(restricted->get_name(), "auto-aux-basis_plus_");
    EXPECT_TRUE(restricted->has_element(1));
    EXPECT_TRUE(restricted->has_element(2));
    ASSERT_EQ(aux.get_ncen(), 2);
    EXPECT_EQ(aux.get_nr_basis_set_loaded(), 2);
    const atom he = aux.get_atom(1);
    ASSERT_EQ(he.get_basis_set_size(), 7u); //STO-3G He is scaled H, same 20.29 exponent ratio, same 7-term series
    EXPECT_NEAR(he.get_basis_set_exponent(he.get_basis_set_size() - 1), 2 * he_exp[2], 1e-12);
    for (unsigned int i = 0; i < he.get_basis_set_size(); i++)
    {
        EXPECT_EQ(he.get_basis_set_type(i), 0) << i;
        EXPECT_EQ(he.get_basis_set_coefficient(i), 1.0) << i;
    }
    EXPECT_EQ(aux.get_atom(0).get_basis_set_size(), (*restricted)[0].size());
    //the orbital wavefunction keeps its own 1-based shells
    EXPECT_EQ(orb.get_atom_basis_set_entry(0, 0).get_type(), 1);
}

//two library-style owned sets for different elements are combined element-wise into the aux
//wavefunction, decontracted with the library's 0-based l
TEST(BasisConstantsLibraryTests, GenerateAuxWfnCombinesSeveralSets)
{
    std::shared_ptr<BasisSet> a = std::make_shared<BasisSet>();
    a->set_name("a");
    a->set_count_for_element(0, 2);
    a->add_owned_primitive({ 0, 0, 2.0, 0.3, 0 });
    a->add_owned_primitive({ 0, 1, 0.5, 0.4, 1 });
    std::shared_ptr<BasisSet> b = std::make_shared<BasisSet>();
    b->set_name("b");
    b->set_count_for_element(1, 1);
    b->add_owned_primitive({ 0, 2, 0.7, 0.5, 0 });

    WFN orb(e_origin::NOT_YET_DEFINED);
    orb.push_back_atom("H", 0.0, 0.0, 0.0, 1);
    orb.push_back_atom("He", 0.0, 0.0, 2.0, 2);
    std::vector<std::shared_ptr<BasisSet>> sets{ a, b };
    const WFN aux = generate_aux_wfn(orb, sets);
    EXPECT_EQ(a->get_name(), "a_plus_b");
    ASSERT_EQ(aux.get_atom_basis_set_size(0), 2);
    ASSERT_EQ(aux.get_atom_basis_set_size(1), 1);
    EXPECT_NEAR(aux.get_atom_basis_set_entry(0, 0).get_exponent(), 2.0, 1e-15);
    EXPECT_EQ(aux.get_atom_basis_set_entry(0, 1).get_type(), 1);
    EXPECT_EQ(aux.get_atom_basis_set_entry(0, 1).get_shell(), 1);
    EXPECT_EQ(aux.get_atom_basis_set_entry(0, 1).get_coefficient(), 1.0);
    EXPECT_NEAR(aux.get_atom_basis_set_entry(1, 0).get_exponent(), 0.7, 1e-15);
    EXPECT_EQ(aux.get_atom_basis_set_entry(1, 0).get_type(), 2);
    EXPECT_EQ(aux.get_atom_label(1), "He");
}

//the tonto-style turbomole fixture shipped with the NiP3 test: every atom of an element gets the
//element's block, shells numbered per block, types s..f as 1..4
TEST(BasisConstantsIoTests, ReadVanillaTurbomoleFixture)
{
    const std::filesystem::path dir = nos_test_repo_root() / "tests" / "NiP3_fchk";
    if (!std::filesystem::exists(dir / "def2-TZVP")) GTEST_SKIP() << "fixture tests/NiP3_fchk/def2-TZVP not found";
    WFN w(e_origin::NOT_YET_DEFINED);
    w.push_back_atom("H", 0.0, 0.0, 0.0, 1);
    w.push_back_atom("H", 0.0, 0.0, 1.4, 1);
    w.push_back_atom("C", 0.0, 2.0, 0.0, 6);
    w.set_basis_set_name("def2-TZVP");
    ASSERT_TRUE(BasisSetLibrary::read_basis_set_vanilla(dir, w, false));
    EXPECT_EQ(w.get_basis_set_name(), "def2-TZVP");
    EXPECT_FALSE(w.get_d_f_switch());
    EXPECT_EQ(w.get_nr_basis_set_loaded(), 3);
    for (int a = 0; a < 2; a++)
    {
        ASSERT_EQ(w.get_atom_basis_set_size(a), 6) << a;
        EXPECT_EQ(w.get_atom(a).get_shellcount_size(), 4u) << a;
        EXPECT_NEAR(w.get_atom_basis_set_entry(a, 0).get_exponent(), 34.061341, 1e-9);
        EXPECT_NEAR(w.get_atom_basis_set_entry(a, 0).get_coefficient(), 0.0060251978, 1e-12);
        EXPECT_EQ(w.get_atom_basis_set_entry(a, 2).get_shell(), 0);
        EXPECT_EQ(w.get_atom_basis_set_entry(a, 3).get_shell(), 1);
        EXPECT_EQ(w.get_atom_basis_set_entry(a, 5).get_type(), 2);
        EXPECT_NEAR(w.get_atom_basis_set_entry(a, 5).get_exponent(), 0.8, 1e-12);
    }
    ASSERT_EQ(w.get_atom_basis_set_size(2), 20);
    EXPECT_EQ(w.get_atom(2).get_shellcount_size(), 11u);
    EXPECT_EQ(w.get_atom(2).get_shellcount(0), 6u);
    EXPECT_EQ(w.get_atom(2).get_shellcount(5), 4u);
    EXPECT_NEAR(w.get_atom_basis_set_entry(2, 0).get_exponent(), 13575.349682, 1e-6);
    EXPECT_EQ(w.get_atom_basis_set_entry(2, 11).get_type(), 2);
    EXPECT_EQ(w.get_atom_basis_set_entry(2, 17).get_type(), 3);
    EXPECT_EQ(w.get_atom_basis_set_entry(2, 19).get_type(), 4);
    EXPECT_EQ(w.get_atom_basis_set_entry(2, 19).get_shell(), 10);
    EXPECT_NEAR(w.get_atom_basis_set_entry(2, 19).get_exponent(), 0.761, 1e-12);
}

//the gamess-us layout (type count, then index exponent coefficient) with the brace on the label
//line, and the gaussian layout (type count, then exponent coefficient) with the brace on its own line
TEST(BasisConstantsIoTests, ReadVanillaGamessAndGaussianLayouts)
{
    const std::filesystem::path dir = std::filesystem::temp_directory_path();
    const std::filesystem::path gamess = temp_file("BasisConstants_gamess.basis",
        "keys= { gamess-us= }\ndata= {\nH: {\ns 2\n1 3.0 0.4\n2 0.5 0.7\np 1\n1 0.8 1.0\n}\n}\n");
    const std::filesystem::path gaussian = temp_file("BasisConstants_gaussian.basis",
        "! comment\nkeys= { gaussian= }\ndata= {\nHe:\n{\nS 1\n9.0 1.0\nD 2\n2.0 0.6\n1.0 0.5\n}\n}\n");

    WFN w(e_origin::NOT_YET_DEFINED);
    w.push_back_atom("H", 0.0, 0.0, 0.0, 1);
    w.push_back_atom("H", 0.0, 0.0, 1.4, 1);
    w.set_basis_set_name("BasisConstants_gamess.basis");
    EXPECT_TRUE(BasisSetLibrary::read_basis_set_vanilla(dir, w, false));
    std::filesystem::remove(gamess);
    for (int a = 0; a < 2; a++)
    {
        ASSERT_EQ(w.get_atom_basis_set_size(a), 3) << a;
        EXPECT_NEAR(w.get_atom_basis_set_entry(a, 0).get_exponent(), 3.0, 1e-15);
        EXPECT_NEAR(w.get_atom_basis_set_entry(a, 1).get_coefficient(), 0.7, 1e-15);
        EXPECT_EQ(w.get_atom_basis_set_entry(a, 1).get_type(), 1);
        EXPECT_EQ(w.get_atom_basis_set_entry(a, 1).get_shell(), 0);
        EXPECT_EQ(w.get_atom_basis_set_entry(a, 2).get_type(), 2);
        EXPECT_EQ(w.get_atom_basis_set_entry(a, 2).get_shell(), 1);
        EXPECT_NEAR(w.get_atom_basis_set_entry(a, 2).get_exponent(), 0.8, 1e-15);
    }

    WFN v(e_origin::NOT_YET_DEFINED);
    v.push_back_atom("He", 0.0, 0.0, 0.0, 2);
    v.set_basis_set_name("BasisConstants_gaussian.basis");
    EXPECT_TRUE(BasisSetLibrary::read_basis_set_vanilla(dir, v, false));
    std::filesystem::remove(gaussian);
    ASSERT_EQ(v.get_atom_basis_set_size(0), 3);
    EXPECT_EQ(v.get_atom(0).get_shellcount_size(), 2u);
    EXPECT_NEAR(v.get_atom_basis_set_entry(0, 0).get_exponent(), 9.0, 1e-15);
    EXPECT_EQ(v.get_atom_basis_set_entry(0, 0).get_type(), 1);
    EXPECT_EQ(v.get_atom_basis_set_entry(0, 2).get_type(), 3);
    EXPECT_EQ(v.get_atom_basis_set_entry(0, 2).get_shell(), 1);
    EXPECT_NEAR(v.get_atom_basis_set_entry(0, 2).get_coefficient(), 0.5, 1e-15);
}

//a CRYSTAL keys line reads like turbomole (count type) and flips the d/f ordering switch
TEST(BasisConstantsIoTests, ReadVanillaCrystalSetsDfSwitch)
{
    const std::filesystem::path dir = std::filesystem::temp_directory_path();
    const std::filesystem::path crystal = temp_file("BasisConstants_crystal.basis",
        "keys= { CRYSTAL= }\nLi: {\n2 s\n3.0 0.4\n0.5 0.7\n1 f\n1.5 1.0\n}\n");
    WFN w(e_origin::NOT_YET_DEFINED);
    w.push_back_atom("Li", 0.0, 0.0, 0.0, 3);
    w.set_basis_set_name("BasisConstants_crystal.basis");
    EXPECT_FALSE(w.get_d_f_switch());
    EXPECT_TRUE(BasisSetLibrary::read_basis_set_vanilla(dir, w, false));
    std::filesystem::remove(crystal);
    EXPECT_TRUE(w.get_d_f_switch());
    ASSERT_EQ(w.get_atom_basis_set_size(0), 3);
    EXPECT_NEAR(w.get_atom_basis_set_entry(0, 1).get_exponent(), 0.5, 1e-15);
    EXPECT_EQ(w.get_atom_basis_set_entry(0, 2).get_type(), 4);
    EXPECT_EQ(w.get_atom_basis_set_entry(0, 2).get_shell(), 1);
}

//every refusal path returns false: no file, a format tonto does not name, a file that never
//reaches the element, and a shell above f
TEST(BasisConstantsIoTests, ReadVanillaRejectsBadInputs)
{
    const std::filesystem::path dir = std::filesystem::temp_directory_path();
    auto make_wfn = [](const std::string& name) {
        WFN w(e_origin::NOT_YET_DEFINED);
        w.push_back_atom("H", 0.0, 0.0, 0.0, 1);
        w.set_basis_set_name(name);
        return w;
    };
    WFN missing = make_wfn("BasisConstants_does_not_exist.basis");
    EXPECT_FALSE(BasisSetLibrary::read_basis_set_vanilla(dir, missing, false));

    const std::filesystem::path unsupported = temp_file("BasisConstants_unsupported.basis", "keys= { molpro= }\nH: {\n1 s\n1.0 1.0\n}\n");
    WFN w_unsupported = make_wfn("BasisConstants_unsupported.basis");
    EXPECT_FALSE(BasisSetLibrary::read_basis_set_vanilla(dir, w_unsupported, false));
    EXPECT_EQ(w_unsupported.get_atom_basis_set_size(0), 0);

    const std::filesystem::path no_element = temp_file("BasisConstants_no_element.basis", "keys= { gaussian= }\nHe: {\ns 1\n1.0 1.0\n}\n");
    WFN w_no_element = make_wfn("BasisConstants_no_element.basis");
    EXPECT_FALSE(BasisSetLibrary::read_basis_set_vanilla(dir, w_no_element, false));
    EXPECT_EQ(w_no_element.get_atom_basis_set_size(0), 0);

    const std::filesystem::path g_shell = temp_file("BasisConstants_g_shell.basis", "keys= { turbomole= }\nH: {\n1 g\n1.0 1.0\n}\n");
    WFN w_g = make_wfn("BasisConstants_g_shell.basis");
    EXPECT_FALSE(BasisSetLibrary::read_basis_set_vanilla(dir, w_g, false));
    EXPECT_EQ(w_g.get_atom_basis_set_size(0), 0);

    std::filesystem::remove(unsupported);
    std::filesystem::remove(no_element);
    std::filesystem::remove(g_shell);
}

//read_basis_set_missing only fills atoms without a basis set and gives them the whole block
TEST(BasisConstantsIoTests, ReadMissingSkipsAtomsThatHaveABasis)
{
    const std::filesystem::path dir = std::filesystem::temp_directory_path();
    const std::filesystem::path file = temp_file("BasisConstants_missing.basis",
        "keys= { gaussian= }\ndata= {\nH:\n{\ns 2\n3.0 0.4\n0.5 0.7\np 1\n0.8 1.0\n}\n}\n");
    WFN w(e_origin::NOT_YET_DEFINED);
    w.push_back_atom("H ", 0.0, 0.0, 0.0, 1);
    w.push_back_atom("H ", 0.0, 0.0, 1.4, 1);
    w.push_back_atom_basis_set(0, 42.0, 1.0, 1, 0);
    w.set_basis_set_name("BasisConstants_missing.basis");
    EXPECT_TRUE(BasisSetLibrary::read_basis_set_missing(dir, w, false));
    std::filesystem::remove(file);
    ASSERT_EQ(w.get_atom_basis_set_size(0), 1);
    EXPECT_NEAR(w.get_atom_basis_set_entry(0, 0).get_exponent(), 42.0, 1e-15);
    ASSERT_EQ(w.get_atom_basis_set_size(1), 3);
    EXPECT_NEAR(w.get_atom_basis_set_entry(1, 0).get_exponent(), 3.0, 1e-15);
    EXPECT_NEAR(w.get_atom_basis_set_entry(1, 1).get_coefficient(), 0.7, 1e-15);
    EXPECT_EQ(w.get_atom_basis_set_entry(1, 1).get_shell(), 0);
    EXPECT_EQ(w.get_atom_basis_set_entry(1, 2).get_type(), 2);
    EXPECT_EQ(w.get_atom_basis_set_entry(1, 2).get_shell(), 1);
    EXPECT_EQ(w.get_nr_basis_set_loaded(), 2);
}

//a label without a blank is read like a padded one
TEST(BasisConstantsIoTests, ReadMissingAcceptsPlainLabels)
{
    const std::filesystem::path dir = std::filesystem::temp_directory_path();
    const std::filesystem::path file = temp_file("BasisConstants_missing_plain.basis",
        "keys= { gaussian= }\nH:\n{\ns 1\n3.0 1.0\n}\n");
    WFN w(e_origin::NOT_YET_DEFINED);
    w.push_back_atom("H", 0.0, 0.0, 0.0, 1);
    w.set_basis_set_name("BasisConstants_missing_plain.basis");
    bool ok = false;
    EXPECT_NO_THROW(ok = BasisSetLibrary::read_basis_set_missing(dir, w, false));
    std::filesystem::remove(file);
    EXPECT_TRUE(ok);
    EXPECT_EQ(w.get_atom_basis_set_size(0), 1);
}
