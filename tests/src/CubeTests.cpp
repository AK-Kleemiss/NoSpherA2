#include "pch.h"

#include "core/convenience.h"
#include "core/constants.h"
#include "core/wfn_class.h"
#include "core/cube.h"
#include "core/isosurface.h"
#include "core/npy.h"

#include <atomic>
#include <complex>
#include <map>
#include <set>

//unused helper in isosurface.cpp with external linkage, declared here so the test can call it
std::vector<d3> subdivideCube(const d3& p1, const d3& p2, int level);

namespace
{
    std::filesystem::path tmp(const std::string& name)
    {
        return std::filesystem::temp_directory_path() / ("nos_cube_" + name);
    }

    //orthogonal grid n^3 with spacing h at origin o, values f(x, y, z) of the node indices
    template <class F>
    cube grid(const int n, const double h, const double o, F f)
    {
        cube c({ n, n, n }, 0, true);
        for (int k = 0; k < 3; k++)
        {
            c.set_origin(k, o);
            c.set_vector(k, k, h);
        }
        c.calc_dv();
        for (int x = 0; x < n; x++)
            for (int y = 0; y < n; y++)
                for (int z = 0; z < n; z++)
                    c.set_value(x, y, z, f(x, y, z));
        return c;
    }

    //2x2x2 cube with values 1..8 in x-major order, spacing 0.5 (dv = 0.125), path stem given
    cube counting_cube(const std::string& stem)
    {
        cube c = grid(2, 0.5, 0.0, [](int x, int y, int z) { return 1.0 + 4 * x + 2 * y + z; });
        c.set_path(tmp(stem + ".cube"));
        return c;
    }

    //one hydrogen at the origin read from a temporary xyz file, the parent every atom-carrying cube needs
    WFN one_hydrogen(const std::string& tag)
    {
        const auto xyz = tmp(tag + ".xyz");
        {
            std::ofstream out(xyz);
            out << "1\n\nH 0 0 0\n";
        }
        WFN w(xyz, false);
        std::filesystem::remove(xyz);
        return w;
    }

    std::vector<std::string> lines_of(const std::filesystem::path& p)
    {
        std::ifstream in(p);
        std::vector<std::string> l;
        for (std::string s; std::getline(in, s);)
            l.push_back(s);
        return l;
    }

    //the number printed after key in a captured log, NaN when the key is missing
    double logged(const std::string& out, const std::string& key)
    {
        const size_t at = out.find(key);
        return at == std::string::npos ? std::nan("") : std::stod(out.substr(at + key.size()));
    }

    template <class F>
    std::string what_throws(F f)
    {
        try
        {
            f();
        }
        catch (const std::exception& e)
        {
            return e.what();
        }
        return "";
    }
}

//the vector-of-vectors constructor copies size, origin, vectors and values and computes dv = |det|
TEST(CubeTests, VectorConstructorCopiesEverything)
{
    const ivec size{ 2, 1, 1 };
    const vec origin{ 1.0, 2.0, 3.0 };
    const vec2 vectors{ { 0.5, 0.0, 0.0 }, { 0.0, 0.5, 0.0 }, { 0.0, 0.0, 2.0 } };
    const vec3 values{ { { 7.0 } }, { { -3.0 } } };
    testing::internal::CaptureStdout();
    cube c(0, size, origin, vectors, values);
    testing::internal::GetCapturedStdout();
    EXPECT_TRUE(c.get_loaded());
    EXPECT_EQ(c.get_sizes(), (i3{ 2, 1, 1 }));
    EXPECT_EQ(c.get_origin(2), 3.0);
    EXPECT_EQ(c.get_vector(2, 2), 2.0);
    EXPECT_EQ(c.get_value(0, 0, 0), 7.0);
    EXPECT_EQ(c.get_value(1, 0, 0), -3.0);
    EXPECT_NEAR(c.get_dv(), 0.5, 1e-15);
    EXPECT_EQ(c.get_size(3), -1);
}

//resize allocates the value array of an empty cube; get/set_value refuse indices outside it
TEST(CubeTests, ResizeAndValueBounds)
{
    cube c;
    EXPECT_EQ(c.get_size(0), 0);
    EXPECT_EQ(c.sum(), -1.0);
    EXPECT_EQ(c.double_sum(), (vec{ 0.0, 0.0 }));
    EXPECT_FALSE(c.thresh(1.0));
    c.resize({ 2, 3, 4 });
    EXPECT_EQ(c.get_sizes(), (i3{ 2, 3, 4 }));
    EXPECT_TRUE(c.set_value(1, 2, 3, 5.0));
    EXPECT_EQ(c.get_value(1, 2, 3), 5.0);
    EXPECT_FALSE(c.set_value(2, 0, 0, 1.0));
    EXPECT_FALSE(c.set_value(0, -1, 0, 1.0));
    EXPECT_EQ(c.get_value(0, 3, 0), -1.0);
    EXPECT_EQ(c.get_value(0, 0, 4), -1.0);
    EXPECT_FALSE(c.set_origin(3, 1.0));
    EXPECT_FALSE(c.set_vector(3, 0, 1.0));
    EXPECT_EQ(c.get_vector(0, 3), -1.0);
}

//set_vectors and set_vector both recompute dv as |det| of a non-diagonal matrix, set_dv overrides
TEST(CubeTests, VolumeElementFollowsTheVectors)
{
    cube c({ 1, 1, 1 }, 0, true);
    c.set_vectors({ { { 1.0, 2.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 3.0 } } });
    EXPECT_NEAR(c.get_dv(), 3.0, 1e-15);
    EXPECT_TRUE(c.set_vector(0, 0, 2.0));
    EXPECT_NEAR(c.get_dv(), 6.0, 1e-15);
    c.set_dv(7.0);
    EXPECT_EQ(c.get_dv(), 7.0);
    EXPECT_EQ(c.get_vectors()[0][1], 2.0);
}

//+ - * / act element-wise and name the result <a><op><b>.cube next to the left operand ("_" for the quotient)
TEST(CubeTests, ArithmeticOperatorsAndResultPaths)
{
    const cube a = counting_cube("a");
    cube b = counting_cube("b");
    for (int x = 0; x < 2; x++)
        for (int y = 0; y < 2; y++)
            for (int z = 0; z < 2; z++)
                b.set_value(x, y, z, 2.0);
    const cube sum = a + b, diff = a - b, prod = a * b, quot = a / b;
    EXPECT_EQ(sum.get_path().filename().string(), "nos_cube_a+nos_cube_b.cube");
    EXPECT_EQ(diff.get_path().filename().string(), "nos_cube_a-nos_cube_b.cube");
    EXPECT_EQ(prod.get_path().filename().string(), "nos_cube_a*nos_cube_b.cube");
    EXPECT_EQ(quot.get_path().filename().string(), "nos_cube_a_nos_cube_b.cube");
    EXPECT_EQ(sum.get_path().parent_path(), a.get_path().parent_path());
    for (int x = 0; x < 2; x++)
        for (int y = 0; y < 2; y++)
            for (int z = 0; z < 2; z++)
            {
                const double v = a.get_value(x, y, z);
                EXPECT_EQ(sum.get_value(x, y, z), v + 2.0);
                EXPECT_EQ(diff.get_value(x, y, z), v - 2.0);
                EXPECT_EQ(prod.get_value(x, y, z), 2.0 * v);
                EXPECT_EQ(quot.get_value(x, y, z), v / 2.0);
            }
    EXPECT_EQ(a.get_value(0, 0, 0), 1.0) << "operands are untouched";
}

//a zero numerator over a zero denominator gives 0 rather than NaN
TEST(CubeTests, ZeroOverZeroIsZero)
{
    cube a = counting_cube("a");
    cube b = counting_cube("b");
    a.set_value(0, 0, 0, 0.0);
    b.set_value(0, 0, 0, 0.0);
    const cube q = a / b;
    EXPECT_EQ(q.get_value(0, 0, 0), 0.0);
    EXPECT_EQ(q.get_value(1, 1, 1), 1.0);
}

//x/0 with x != 0 is the 1E100 guard value
TEST(CubeTests, NonzeroOverZeroIsTheGuardValue)
{
    cube a = counting_cube("a");
    cube b = counting_cube("b");
    b.set_value(0, 0, 0, 0.0);
    const cube q = a / b;
    EXPECT_EQ(q.get_value(0, 0, 0), 1E100);
}

//a size mismatch makes the binary operators return an empty cube and the compound ones return false
TEST(CubeTests, SizeMismatchIsRejected)
{
    cube a = counting_cube("a");
    const cube other({ 2, 2, 3 }, 0, true);
    EXPECT_EQ((a + other).get_size(0), 0);
    EXPECT_EQ((a - other).get_size(1), 0);
    EXPECT_EQ((a * other).get_size(2), 0);
    EXPECT_EQ((a / other).get_size(0), 0);
    EXPECT_FALSE(a += other);
    EXPECT_FALSE(a -= other);
    EXPECT_FALSE(a *= other);
    EXPECT_FALSE(a /= other);
    EXPECT_FALSE(a.mask(other));
    EXPECT_FALSE(a.thresh(other));
    EXPECT_FALSE(a.negative_mask(other));
    EXPECT_EQ(a.rrs(other), -1.0);
    EXPECT_EQ(a.jaccard(other), -1.0);
    EXPECT_EQ(a.get_value(0, 0, 0), 1.0) << "a rejected operation leaves the values alone";
}

//compound operators modify in place, one after the other: ((a += b) -= b) *= b then /= b returns to a
TEST(CubeTests, CompoundOperatorsInPlace)
{
    cube a = counting_cube("a");
    cube b = counting_cube("b");
    EXPECT_TRUE(a += b);
    EXPECT_EQ(a.get_value(1, 1, 1), 16.0);
    EXPECT_TRUE(a -= b);
    EXPECT_EQ(a.get_value(1, 1, 1), 8.0);
    EXPECT_TRUE(a *= b);
    EXPECT_EQ(a.get_value(1, 0, 1), 36.0);
    EXPECT_TRUE(a /= b);
    for (int x = 0; x < 2; x++)
        for (int y = 0; y < 2; y++)
            for (int z = 0; z < 2; z++)
                EXPECT_EQ(a.get_value(x, y, z), 1.0 + 4 * x + 2 * y + z);
}

//mask keeps where the mask is nonzero, negative_mask where it is zero, thresh(cube, t) where the mask is >= t, thresh(t) where the value is >= t
TEST(CubeTests, MaskAndThreshold)
{
    cube m = grid(2, 0.5, 0.0, [](int x, int y, int z) { return double((4 * x + 2 * y + z) % 2); });
    cube a = counting_cube("a");
    EXPECT_TRUE(a.mask(m));
    EXPECT_EQ(a.get_value(0, 0, 0), 0.0);
    EXPECT_EQ(a.get_value(0, 0, 1), 2.0);
    EXPECT_EQ(a.get_value(1, 1, 0), 0.0);
    EXPECT_EQ(a.get_value(1, 1, 1), 8.0);

    cube n = counting_cube("a");
    EXPECT_TRUE(n.negative_mask(m));
    EXPECT_EQ(n.get_value(0, 0, 0), 1.0);
    EXPECT_EQ(n.get_value(0, 0, 1), 0.0);
    EXPECT_EQ(n.get_value(1, 1, 1), 0.0);

    cube t = counting_cube("a");
    EXPECT_TRUE(t.thresh(m, 0.5));
    EXPECT_EQ(t.get_value(0, 0, 0), 0.0);
    EXPECT_EQ(t.get_value(0, 0, 1), 2.0);
    cube d = counting_cube("a");
    EXPECT_TRUE(d.thresh(m));
    EXPECT_EQ(d.get_value(0, 0, 0), 1.0) << "the default threshold of -1234 zeroes nothing";

    cube s = counting_cube("a");
    EXPECT_TRUE(s.thresh(5.0));
    EXPECT_EQ(s.get_value(0, 1, 1), 0.0);
    EXPECT_EQ(s.get_value(1, 0, 0), 5.0);
    EXPECT_EQ(s.get_value(1, 1, 1), 8.0);
}

//sums use dv, rrs is sum|a-b| / sum|a+b| and jaccard sum min / sum max, on a = 1..8 and b = 2a
TEST(CubeTests, SumsAndSimilarityMeasures)
{
    const cube a = counting_cube("a");
    const cube b = a * a * (a / a) - a;
    EXPECT_NEAR(a.get_dv(), 0.125, 1e-15);
    EXPECT_NEAR(a.sum(), 36.0 * 0.125, 1e-12);
    EXPECT_NEAR(a.diff_sum(), 18.0 * 0.125, 1e-12);
    const vec ds = a.double_sum();
    ASSERT_EQ(ds.size(), 2u);
    EXPECT_NEAR(ds[0], 18.0 * 0.125, 1e-12);
    EXPECT_NEAR(ds[1], 36.0 * 0.125, 1e-12);
    EXPECT_EQ(a.min_value(), 1.0);
    EXPECT_EQ(a.max_value(), 8.0);
    for (int x = 0; x < 2; x++)
        for (int y = 0; y < 2; y++)
            for (int z = 0; z < 2; z++)
                EXPECT_EQ(b.get_value(x, y, z), a.get_value(x, y, z) * (a.get_value(x, y, z) - 1.0));
    const cube two_a = a + a;
    EXPECT_NEAR(a.rrs(two_a), 1.0 / 3.0, 1e-12);
    EXPECT_NEAR(a.jaccard(two_a), 0.5, 1e-12);
    EXPECT_NEAR(a.rrs(a), 0.0, 1e-15);
    EXPECT_NEAR(a.jaccard(a), 1.0, 1e-15);
}

//trilinear interpolation reproduces a linear field exactly inside the grid and is 0 outside the box
TEST(CubeTests, InterpolationOfALinearField)
{
    const cube c = grid(4, 0.5, 0.0, [](int x, int y, int z) { return 1.0 + 2 * (0.5 * x) + 3 * (0.5 * y) + 4 * (0.5 * z); });
    EXPECT_NEAR(c.get_interpolated_value(0.7, 0.3, 0.9), 1.0 + 1.4 + 0.9 + 3.6, 1e-12);
    EXPECT_NEAR(c.get_interpolated_value(0.0, 0.0, 0.0), 1.0, 1e-12);
    EXPECT_NEAR(c.get_interpolated_value(1.0, 1.0, 1.0), 10.0, 1e-12);
    EXPECT_EQ(c.get_interpolated_value(-0.1, 0.5, 0.5), 0.0);
    EXPECT_EQ(c.get_interpolated_value(0.5, 2.1, 0.5), 0.0);
    EXPECT_EQ(c.get_interpolated_value(0.5, 0.5, 5.0), 0.0);
}

//the last node is exact and points beyond it are 0 like every other point outside the data (the old expectation
//of 4.5 at 1.75 was extrapolation, which the interpolator does not do)
TEST(CubeTests, InterpolationInTheLastCell)
{
    const cube c = grid(4, 0.5, 0.0, [](int x, int, int) { return 1.0 + 2 * (0.5 * x); });
    EXPECT_NEAR(c.get_interpolated_value(1.5, 0.0, 0.0), 4.0, 1e-12);
    EXPECT_EQ(c.get_interpolated_value(1.75, 0.0, 0.0), 0.0);
}

//evaluate_on_grid stores func(node position); with wrap every node collects the 27 periodic images of a constant
TEST(CubeTests, EvaluateOnGridPlainAndWrapped)
{
    cube c = grid(3, 1.0, -1.0, [](int, int, int) { return 0.0; });
    EXPECT_TRUE(c.evaluate_on_grid([](const d3& p) { return p[0] + 10 * p[1] + 100 * p[2]; }));
    EXPECT_TRUE(c.get_loaded());
    for (int x = 0; x < 3; x++)
        for (int y = 0; y < 3; y++)
            for (int z = 0; z < 3; z++)
                EXPECT_NEAR(c.get_value(x, y, z), (x - 1.0) + 10 * (y - 1.0) + 100 * (z - 1.0), 1e-12);

    std::atomic<int> bad_mapping{ 0 };
    std::atomic<int> calls{ 0 };
    EXPECT_TRUE(c.evaluate_on_grid([&](const d3&, const i3& raw, const i3& mapped) {
        calls++;
        for (int d = 0; d < 3; d++)
            if (mapped[d] != ((raw[d] % 3) + 3) % 3)
                bad_mapping++;
        return 1.0;
        }, true));
    EXPECT_EQ(calls.load(), 9 * 9 * 9);
    EXPECT_EQ(bad_mapping.load(), 0);
    for (int x = 0; x < 3; x++)
        for (int y = 0; y < 3; y++)
            for (int z = 0; z < 3; z++)
                EXPECT_EQ(c.get_value(x, y, z), 27.0);

    cube empty;
    EXPECT_FALSE(empty.evaluate_on_grid([](const d3&) { return 1.0; }));
}

//find_value_bounds boxes the non-ignored values with padding, skips NaN, and shrink_to_values crops the grid to it
TEST(CubeTests, ShrinkToTheNonzeroBlock)
{
    auto block = [](int x, int y, int z) { return (x >= 2 && x <= 3 && y >= 2 && y <= 3 && z >= 2 && z <= 3) ? 1.0 : 0.0; };
    cube c = grid(6, 0.5, -1.0, block);
    c.set_value(0, 0, 0, std::nan(""));
    i3 lower, upper;
    ASSERT_TRUE(c.find_value_bounds(lower, upper, 0.0));
    EXPECT_EQ(lower, (i3{ 1, 1, 1 }));
    EXPECT_EQ(upper, (i3{ 4, 4, 4 }));
    ASSERT_TRUE(c.find_value_bounds(lower, upper, 0.0, 1e-12, 0));
    EXPECT_EQ(lower, (i3{ 2, 2, 2 }));
    EXPECT_EQ(upper, (i3{ 3, 3, 3 }));
    ASSERT_TRUE(c.find_value_bounds(lower, upper, 0.0, 1e-12, 5));
    EXPECT_EQ(lower, (i3{ 0, 0, 0 }));
    EXPECT_EQ(upper, (i3{ 5, 5, 5 })) << "padding is clamped to the grid";

    EXPECT_FALSE(c.find_value_bounds(lower, upper, 1.0, 1.5)) << "everything within tolerance of the ignore value";
    EXPECT_FALSE(c.shrink_to_bounds({ 2, 2, 2 }, { 1, 3, 3 }));
    EXPECT_FALSE(c.shrink_to_bounds({ 0, 0, 0 }, { 6, 5, 5 }));
    EXPECT_FALSE(c.shrink_to_bounds({ -1, 0, 0 }, { 5, 5, 5 }));

    const d3 expected_origin = c.get_pos(1, 1, 1);
    ASSERT_TRUE(c.shrink_to_values(0.0));
    EXPECT_EQ(c.get_sizes(), (i3{ 4, 4, 4 }));
    for (int k = 0; k < 3; k++)
        EXPECT_NEAR(c.get_origin(k), expected_origin[k], 1e-15);
    EXPECT_EQ(c.get_value(1, 1, 1), 1.0);
    EXPECT_EQ(c.get_value(0, 0, 0), 0.0);
    EXPECT_EQ(c.get_value(3, 3, 3), 0.0);
    EXPECT_NEAR(c.sum(), 8.0 * 0.125, 1e-12);

    cube unloaded({ 2, 2, 2 }, 0, false);
    EXPECT_FALSE(unloaded.find_value_bounds(lower, upper, 0.0));
    EXPECT_FALSE(unloaded.shrink_to_bounds({ 0, 0, 0 }, { 1, 1, 1 }));
    EXPECT_FALSE(unloaded.shrink_to_values(0.0));
}

//cube_t<int> built from a cube rounds half away from zero and copies the geometry
TEST(CubeTests, IntegerCubeRoundsValues)
{
    cube c({ 1, 1, 3 }, 2, true);
    c.set_vector(2, 2, 0.5);
    c.set_origin(1, 4.0);
    c.set_path(tmp("int.cube"));
    c.set_value(0, 0, 0, 2.5);
    c.set_value(0, 0, 1, -2.5);
    c.set_value(0, 0, 2, 1.4);
    const cubei i(c);
    EXPECT_TRUE(i.get_loaded());
    EXPECT_EQ(i.get_na(), 2);
    EXPECT_EQ(i.get_sizes(), (i3{ 1, 1, 3 }));
    EXPECT_EQ(i.get_origin(1), 4.0);
    EXPECT_EQ(i.get_vector(2, 2), 0.5);
    EXPECT_EQ(i.get_path(), c.get_path());
    EXPECT_EQ(i.get_t_value(0, 0, 0), 3);
    EXPECT_EQ(i.get_t_value(0, 0, 1), -3);
    EXPECT_EQ(i.get_t_value(0, 0, 2), 1);
    EXPECT_EQ(i.get_value(0, 0, 3), -1);
    cubed d({ 1, 1, 1 }, 0, true);
    EXPECT_TRUE(d.set_t_value(0, 0, 0, 0.25));
    EXPECT_EQ(d.get_t_value(0, 0, 0), 0.25);
}

//the Ewald energy is quadratic in the charges: scaling the density by -2 scales the energy by 4, and no charge means no energy;
//the self term is -alpha / sqrt(pi) * sum q^2 * dv^2 with alpha = 2 sqrt(pi) / L = sqrt(pi) / 2 on a 4 bohr cell, so -1 for
//two unit charges and -4 for two double ones (the quadratic check alone would pass any wrong prefactor)
TEST(CubeTests, EwaldEnergyIsQuadraticInTheCharges)
{
    auto dipole = [](int x, int y, int z) { return (x == 1 && y == 1 && z == 1) ? 1.0 : (x == 2 && y == 2 && z == 2) ? -1.0 : 0.0; };
    cube a = grid(4, 1.0, 0.0, dipole);
    cube b = grid(4, 1.0, 0.0, [&](int x, int y, int z) { return -2.0 * dipole(x, y, z); });
    testing::internal::CaptureStdout();
    const double e_a = a.ewald_sum(6);
    const std::string out_a = testing::internal::GetCapturedStdout();
    testing::internal::CaptureStdout();
    const double e_b = b.ewald_sum(6);
    const std::string out = testing::internal::GetCapturedStdout();
    EXPECT_TRUE(std::isfinite(e_a));
    EXPECT_NE(e_a, 0.0);
    EXPECT_NEAR(e_b, 4.0 * e_a, 1e-9 * std::abs(e_a));
    EXPECT_NEAR(logged(out_a, "alpha: "), constants::sqr_pi / 2, 1e-5);
    EXPECT_NEAR(logged(out_a, "Self-energy: "), -1.0, 1e-9);
    EXPECT_NEAR(logged(out, "Self-energy: "), -4.0, 1e-9);
    EXPECT_NEAR(logged(out, "Total energy: "), e_b, 1e-5 * std::abs(e_b));

    cube zero = grid(2, 1.0, 0.0, [](int, int, int) { return 0.0; });
    testing::internal::CaptureStdout();
    const double e_zero = zero.ewald_sum(2);
    const std::string zero_out = testing::internal::GetCapturedStdout();
    EXPECT_EQ(e_zero, 0.0);
    //the 0/0 relative criterion is NaN; under /fp:fast (NoSpherA2Optimizations.cmake) the MSVC build prints no
    //verdict line at all for it, so only the energy and the summary lines are checked
    EXPECT_NE(zero_out.find("Total energy: 0"), std::string::npos) << zero_out;
}

//the real-space energy of +1 at (1,1,1) and -1 at (2,2,2) is 0.5 dv^2 * 2 * (-1) erfc(alpha sqrt3) / sqrt3
TEST(CubeTests, EwaldRealSpaceTermOfADipole)
{
    cube a = grid(4, 1.0, 0.0, [](int x, int y, int z) { return (x == 1 && y == 1 && z == 1) ? 1.0 : (x == 2 && y == 2 && z == 2) ? -1.0 : 0.0; });
    testing::internal::CaptureStdout();
    a.ewald_sum(6);
    const std::string out = testing::internal::GetCapturedStdout();
    const double alpha = constants::sqr_pi / 2, r = std::sqrt(3.0);
    EXPECT_NEAR(logged(out, "Real-space energy: "), -std::erfc(alpha * r) / r, 1e-6);
}

//super_cube(x, y, z) tiles the values x, y, z times, puts the image atoms into the result's own WFN and renames the path
TEST(CubeTests, SuperCubeTilesTheValues)
{
    const auto xyz = tmp("super.xyz");
    {
        std::ofstream out(xyz);
        out << "1\n\nH 0 0 0\n";
    }
    WFN parent(xyz, false);
    cube c = grid(2, 0.5, -0.25, [](int x, int y, int z) { return 1.0 + 4 * x + 2 * y + z; });
    c.set_na(1);
    c.give_parent_wfn(parent);
    c.set_path(tmp("super.cube"));
    c.set_comment1("c1");
    c.set_comment2("c2");
    std::filesystem::remove(xyz);
    const cube s = c.super_cube(2, 1, 3);
    EXPECT_EQ(s.get_sizes(), (i3{ 4, 2, 6 }));
    EXPECT_EQ(s.get_na(), 6);
    ASSERT_EQ(s.get_parent_wfn_atoms().size(), 6u);
    EXPECT_NEAR(s.get_parent_wfn_atoms()[5].get_coordinate(0), 1.0, 1e-12) << "second x image is one cell (2 * 0.5) over";
    EXPECT_NEAR(s.get_parent_wfn_atoms()[5].get_coordinate(2), 2.0, 1e-12) << "third z image is two cells over";
    EXPECT_EQ(s.get_comment1(), "c1 SUPER CUBE");
    EXPECT_EQ(s.get_comment2(), "c2");
    EXPECT_EQ(s.get_path().extension().string(), ".cube_super");
    EXPECT_EQ(s.get_origin(0), -0.25);
    EXPECT_EQ(s.get_vector(1, 1), 0.5);
    for (int x = 0; x < 4; x++)
        for (int y = 0; y < 2; y++)
            for (int z = 0; z < 6; z++)
                EXPECT_EQ(s.get_value(x, y, z), 1.0 + 4 * (x % 2) + 2 * y + z % 2) << x << " " << y << " " << z;
}

//adaptive_refine prints the initial integral as sum(); values within [0, 1] have |inhomogeneity| < 1 so a target of 1 leaves every point fine and converges at once
TEST(CubeTests, AdaptiveRefineConvergesOnAnEasyTarget)
{
    cube c = grid(4, 0.5, -1.0, [](int, int, int) { return 0.0; });
    auto gauss = [](const d3 p) -> const double { return std::exp(-(p[0] * p[0] + p[1] * p[1] + p[2] * p[2])); };
    c.evaluate_on_grid([&](const d3& p) { return gauss(p); });
    const double before = c.sum();
    testing::internal::CaptureStdout();
    c.adaptive_refine(gauss, 1.0, 2);
    const std::string out = testing::internal::GetCapturedStdout();
    std::ostringstream expected;
    expected << "Initial Integral: " << before;
    EXPECT_NE(out.find(expected.str()), std::string::npos) << out;
    EXPECT_NE(out.find("Out of 64 initially fine: 64"), std::string::npos) << out;
    EXPECT_NE(out.find("Refinement Level 1"), std::string::npos) << out;
    EXPECT_NE(out.find("Converged!"), std::string::npos) << out;
    EXPECT_NEAR(c.sum(), before, 1e-15) << "the coarse grid is not modified";
}

//a tight target on a smooth field sends points into refinement and reports what was computed
//(last resort: only the log lines are checked, the refined integral is not exposed)
TEST(CubeTests, AdaptiveRefineRefinesOnATightTarget)
{
    cube c = grid(4, 0.5, -1.0, [](int, int, int) { return 0.0; });
    auto gauss = [](const d3 p) -> const double { return std::exp(-(p[0] * p[0] + p[1] * p[1] + p[2] * p[2])); };
    c.evaluate_on_grid([&](const d3& p) { return gauss(p); });
    testing::internal::CaptureStdout();
    c.adaptive_refine(gauss, 1e-12, 1);
    const std::string out = testing::internal::GetCapturedStdout();
    EXPECT_NE(out.find("Out of 64 initially fine: 0"), std::string::npos) << out;
    EXPECT_NE(out.find("Refinement Level 1"), std::string::npos) << out;
    EXPECT_EQ(out.find("Converged!"), std::string::npos) << out;
}

//the 32-point inhomogeneity of a constant field is 0, so the 8 interior points of a 6^3 grid are initially fine
TEST(CubeTests, AdaptiveRefineConstantFieldIsHomogeneous)
{
    cube c = grid(6, 0.5, -1.0, [](int, int, int) { return 1.0; });
    testing::internal::CaptureStdout();
    c.adaptive_refine([](const d3) -> const double { return 1.0; }, 1e-3, 1);
    const std::string out = testing::internal::GetCapturedStdout();
    EXPECT_NE(out.find("Out of 216 initially fine: 8"), std::string::npos) << out;
}

//write_file(force) then the reading constructor: header, atoms (into the fresh WFN) and values survive the round trip;
//an existing file is refused without force; absolute writes |v|
TEST(CubeIoTests, WriteReadRoundTrip)
{
    WFN atoms(e_origin::NOT_YET_DEFINED);
    atoms.push_back_atom("O", 0.5, -1.0, 2.0, 8);
    atoms.push_back_atom("H", 1.5, 0.0, 0.0, 1);
    cube c = grid(2, 0.5, -0.25, [](int x, int y, int z) { return 0.25 * (4 * x + 2 * y + z) - 0.5; });
    c.set_na(2);
    c.give_parent_wfn(atoms);
    c.set_comment1("first comment");
    c.set_comment2("second comment");
    const auto path = tmp("roundtrip.cube");
    c.set_path(path);
    std::filesystem::remove(path);
    ASSERT_TRUE(c.write_file(true));
    ASSERT_TRUE(std::filesystem::exists(path));

    testing::internal::CaptureStdout();
    EXPECT_FALSE(c.write_file(false));
    EXPECT_NE(testing::internal::GetCapturedStdout().find("File already exists, aborting!"), std::string::npos);

    WFN fresh(e_origin::NOT_YET_DEFINED);
    std::ostringstream log;
    const cube r(path, true, fresh, log);
    EXPECT_TRUE(r.get_loaded());
    EXPECT_EQ(r.get_comment1(), "first comment");
    EXPECT_EQ(r.get_comment2(), "second comment");
    EXPECT_EQ(r.get_na(), 2);
    EXPECT_EQ(r.get_sizes(), (i3{ 2, 2, 2 }));
    EXPECT_NEAR(r.get_origin(0), -0.25, 1e-6);
    EXPECT_NEAR(r.get_vector(1, 1), 0.5, 1e-6);
    EXPECT_NEAR(r.get_vector(0, 1), 0.0, 1e-6);
    EXPECT_NEAR(r.get_dv(), 0.125, 1e-9);
    ASSERT_EQ(fresh.get_ncen(), 2);
    EXPECT_EQ(fresh.get_atom_charge(0), 8);
    EXPECT_EQ(fresh.get_atom_charge(1), 1);
    EXPECT_NEAR(fresh.get_atom_coordinate(0, 1), -1.0, 1e-6);
    EXPECT_NEAR(fresh.get_atom_coordinate(1, 0), 1.5, 1e-6);
    for (int x = 0; x < 2; x++)
        for (int y = 0; y < 2; y++)
            for (int z = 0; z < 2; z++)
                EXPECT_EQ(r.get_value(x, y, z), c.get_value(x, y, z)) << "values k/4 are exact in %.5E";

    ASSERT_TRUE(c.write_file(true, true));
    WFN again(e_origin::NOT_YET_DEFINED);
    const cube abs_cube(path, true, again, log);
    EXPECT_EQ(abs_cube.get_value(0, 0, 0), 0.5);
    EXPECT_EQ(abs_cube.get_value(0, 0, 1), 0.25);
    EXPECT_EQ(abs_cube.get_value(1, 1, 1), 1.25);
    std::filesystem::remove(path);
}

//a header-only cube copies the raw value lines when written to a new path, and write_file(path, debug) reports progress
TEST(CubeIoTests, HeaderOnlyCubeCopiesValuesToANewPath)
{
    WFN atoms(e_origin::NOT_YET_DEFINED);
    atoms.push_back_atom("C", 0.0, 0.0, 0.0, 6);
    cube c = grid(2, 1.0, 0.0, [](int x, int y, int z) { return 0.25 * (4 * x + 2 * y + z); });
    c.set_na(1);
    c.give_parent_wfn(atoms);
    const auto src = tmp("headeronly.cube"), dst = tmp("headeronly_copy.cube"), dbg = tmp("headeronly_debug.cube");
    c.set_path(src);
    std::filesystem::remove(src);
    ASSERT_TRUE(c.write_file(true));

    WFN w(e_origin::NOT_YET_DEFINED);
    std::ostringstream log;
    cube h(src, false, w, log);
    EXPECT_FALSE(h.get_loaded());
    EXPECT_EQ(h.get_na(), 1);
    EXPECT_EQ(w.get_ncen(), 1);
    ASSERT_TRUE(h.write_file(dst));
    EXPECT_EQ(h.get_path(), dst);

    WFN w2(e_origin::NOT_YET_DEFINED);
    const cube copy(dst, true, w2, log);
    EXPECT_EQ(copy.get_na(), 1);
    for (int x = 0; x < 2; x++)
        for (int y = 0; y < 2; y++)
            for (int z = 0; z < 2; z++)
                EXPECT_EQ(copy.get_value(x, y, z), c.get_value(x, y, z));

    testing::internal::CaptureStdout();
    ASSERT_TRUE(c.write_file(dbg, true));
    const std::string out = testing::internal::GetCapturedStdout();
    EXPECT_NE(out.find("Finished atoms!"), std::string::npos);
    EXPECT_NE(out.find("Write Z-line!"), std::string::npos);
    EXPECT_EQ(lines_of(dbg).size(), lines_of(src).size());
    std::filesystem::remove(src);
    std::filesystem::remove(dst);
    std::filesystem::remove(dbg);
}

//write_file(path) writes all three origin components in their own 12-wide field, so a copy's header parses
TEST(CubeIoTests, HeaderOnlyCopyKeepsTheOrigin)
{
    WFN atoms(e_origin::NOT_YET_DEFINED);
    cube c = grid(2, 1.0, 0.0, [](int x, int, int) { return double(x); });
    c.set_origin(2, 1.5);
    c.give_parent_wfn(atoms);
    const auto src = tmp("origin.cube"), dst = tmp("origin_copy.cube");
    c.set_path(src);
    std::filesystem::remove(src);
    ASSERT_TRUE(c.write_file(true));
    WFN w(e_origin::NOT_YET_DEFINED);
    std::ostringstream log;
    cube h(src, false, w, log);
    ASSERT_TRUE(h.write_file(dst));
    WFN w2(e_origin::NOT_YET_DEFINED);
    const cube copy(dst, true, w2, log);
    EXPECT_NEAR(copy.get_origin(1), 0.0, 1e-6);
    EXPECT_NEAR(copy.get_origin(2), 1.5, 1e-6);
    std::filesystem::remove(src);
    std::filesystem::remove(dst);
}

//write_xdgraph reads the values of a header-only cube and writes them x-outer, z, then y fastest after "! Values"
TEST(CubeIoTests, XdGraphOrderIsXZY)
{
    WFN atoms(e_origin::NOT_YET_DEFINED);
    atoms.push_back_atom("N", 1.0, 2.0, 3.0, 7);
    cube c = grid(2, 1.0, 0.0, [](int x, int y, int z) { return 0.25 * (4 * x + 2 * y + z); });
    c.set_na(1);
    c.give_parent_wfn(atoms);
    const auto src = tmp("xd.cube"), dst = tmp("xd.grd");
    c.set_path(src);
    std::filesystem::remove(src);
    ASSERT_TRUE(c.write_file(true));

    WFN w(e_origin::NOT_YET_DEFINED);
    std::ostringstream log;
    cube h(src, false, w, log);
    testing::internal::CaptureStdout();
    ASSERT_TRUE(h.write_xdgraph(dst, true));
    EXPECT_NE(testing::internal::GetCapturedStdout().find("Write Y-line!"), std::string::npos);
    EXPECT_TRUE(h.get_loaded());
    EXPECT_EQ(h.get_path(), dst);

    const std::vector<std::string> l = lines_of(dst);
    ASSERT_GT(l.size(), 10u);
    EXPECT_EQ(l[0], "2DGRDFIL  0");
    EXPECT_EQ(l[3], "! Gridpoints, Origin, Physical Dimensions");
    EXPECT_EQ(l[7], "! Objects");
    EXPECT_NE(l[9].find(atoms.get_atom_label(0)), std::string::npos);
    EXPECT_NE(l[9].find(" ATOM"), std::string::npos);
    size_t values_at = 0;
    for (size_t i = 0; i < l.size(); i++)
        if (l[i] == "! Values") values_at = i;
    ASSERT_GT(values_at, 0u);
    vec numbers;
    for (size_t i = values_at + 1; i < l.size(); i++)
    {
        std::istringstream iss(l[i]);
        for (double v; iss >> v;) numbers.push_back(v);
    }
    ASSERT_EQ(numbers.size(), 8u);
    size_t n = 0;
    for (int x = 0; x < 2; x++)
        for (int z = 0; z < 2; z++)
            for (int y = 0; y < 2; y++)
                EXPECT_EQ(numbers[n++], c.get_value(x, y, z)) << x << " " << y << " " << z;
    std::filesystem::remove(src);
    std::filesystem::remove(dst);
}

//a file that ends before all values are read fails read_file with an EOF message instead of a partial cube
TEST(CubeIoTests, TruncatedFileIsRejected)
{
    WFN atoms(e_origin::NOT_YET_DEFINED);
    cube c = grid(2, 1.0, 0.0, [](int x, int y, int z) { return 1.0 * (4 * x + 2 * y + z); });
    c.give_parent_wfn(atoms);
    const auto src = tmp("truncated.cube");
    c.set_path(src);
    std::filesystem::remove(src);
    ASSERT_TRUE(c.write_file(true));
    std::vector<std::string> l = lines_of(src);
    ASSERT_EQ(l.size(), 10u) << "6 header lines, no atoms, four rows of two values";
    {
        std::ofstream out(src);
        for (size_t i = 0; i < 7; i++) out << l[i] << "\n";
    }
    WFN w(e_origin::NOT_YET_DEFINED);
    std::ostringstream log;
    cube h(src, false, w, log);
    testing::internal::CaptureStdout();
    EXPECT_FALSE(h.read_file(true, false));
    const std::string out = testing::internal::GetCapturedStdout();
    EXPECT_NE(out.find("ENCOUNTERED EOF!"), std::string::npos) << out;
    EXPECT_FALSE(h.get_loaded());
    std::filesystem::remove(src);
}

//a value stream that runs on across a z-row break (six per line, the row ends mid-line) is read; the reader
//supports a spill of less than one z-row, so the grid is 1x2x8 (the old 2x2x2 layout spilled a whole row and
//was outside what read_values handles)
TEST(CubeIoTests, ContinuousValueStreamIsRead)
{
    const auto src = tmp("stream.cube");
    {
        std::ofstream out(src);
        out << "c1\nc2\n    0    0.000000    0.000000    0.000000\n"
            << "    1    1.000000    0.000000    0.000000\n"
            << "    2    0.000000    1.000000    0.000000\n"
            << "    8    0.000000    0.000000    1.000000\n"
            << "  0.0E+00  1.0E+00  2.0E+00  3.0E+00  4.0E+00  5.0E+00\n"
            << "  6.0E+00  7.0E+00  8.0E+00  9.0E+00  1.0E+01  1.1E+01\n"
            << "  1.2E+01  1.3E+01  1.4E+01  1.5E+01\n";
    }
    WFN w(e_origin::NOT_YET_DEFINED);
    std::ostringstream log;
    cube h(src, false, w, log);
    ASSERT_TRUE(h.read_file(true, false));
    EXPECT_EQ(h.get_value(0, 0, 7), 7.0);
    EXPECT_EQ(h.get_value(0, 1, 0), 8.0);
    EXPECT_EQ(h.get_value(0, 1, 3), 11.0);
    EXPECT_EQ(h.get_value(0, 1, 4), 12.0);
    EXPECT_EQ(h.get_value(0, 1, 7), 15.0);
    std::filesystem::remove(src);
}

//fractal_dimension on values = x index: every iso level between two integer values is crossed on all 16 x-pairs
//and nowhere else, df = ln(16) / ln(comparisons^(1/3)); the header carries the value range and the two sums
TEST(CubeIoTests, FractalDimensionOfAStep)
{
    cube c = grid(4, 1.0, 0.0, [](int x, int, int) { return double(x); });
    const auto src = tmp("fractal.cube");
    const auto plot = tmp("fractal.cube_fractal_plot");
    c.set_path(src);
    EXPECT_FALSE(c.fractal_dimension(0.0));
    ASSERT_TRUE(c.fractal_dimension(0.5));
    ASSERT_TRUE(std::filesystem::exists(plot));
    const std::vector<std::string> l = lines_of(plot);
    ASSERT_EQ(l.size(), 13u);
    std::istringstream head(l[0]);
    int steps;
    double map_min, map_max, e0, e1;
    head >> steps >> map_min >> map_max >> e0 >> e1;
    EXPECT_EQ(steps, 12);
    EXPECT_EQ(map_min, 0.0);
    EXPECT_EQ(map_max, 3.0);
    const double a3 = std::pow(0.529177249, 3);
    EXPECT_NEAR(e0, 48.0 * a3, 1e-6);
    EXPECT_NEAR(e1, 96.0 * a3, 1e-6);
    const double epsilon = std::log(144.0) / 3.0;
    for (int i = 0; i < 12; i++)
    {
        std::istringstream iss(l[i + 1]);
        double iso, df;
        iss >> iso >> df;
        EXPECT_NEAR(iso, -1.0 + 0.5 * i, 1e-6) << i;
        const bool crossing = iso == 0.5 || iso == 1.5 || iso == 2.5;
        EXPECT_NEAR(df, crossing ? std::log(16.0) / epsilon : 0.0, 1e-6) << "iso " << iso;
    }
    std::filesystem::remove(plot);
    cube flat({ 1, 4, 4 }, 0, true);
    EXPECT_FALSE(flat.fractal_dimension(0.5)) << "fewer than two nodes in a direction";
}

//box_cube spans the atoms plus opts.radius on both sides with ceil(extent / resolution) steps
TEST(CubeIsoTests, BoxCubeAroundOneAtom)
{
    WFN h = one_hydrogen("box");
    properties_options opts;
    opts.radius = 2.0;
    opts.resolution = 0.5;
    const cube box = box_cube(h, opts);
    const double r = constants::ang2bohr(2.0);
    EXPECT_EQ(box.get_sizes(), (i3{ 8, 8, 8 }));
    EXPECT_EQ(box.get_na(), 1);
    EXPECT_TRUE(box.get_loaded());
    for (int k = 0; k < 3; k++)
    {
        EXPECT_NEAR(box.get_origin(k), -r, 1e-12);
        EXPECT_NEAR(box.get_vector(k, k), 2 * r / 8, 1e-12);
    }
    EXPECT_EQ(box.get_vector(0, 1), 0.0);
    EXPECT_EQ(box.get_value(7, 7, 7), 0.0);
}

//marching cubes on R^2 - r^2 at iso 0 gives a closed sphere: vertices at |r| = R, area 4 pi R^2 and volume 4/3 pi R^3,
//and every directed edge is matched by its reverse exactly once (closed and consistently oriented, so a dropped table
//entry or a flipped winding fails here where the 8 % area tolerance would not); the origin -2.53 keeps every node off
//the sphere so no zero-length edges appear
TEST(CubeIsoTests, SphereAreaAndVolume)
{
    const double R = 2.0, h = 0.25, o = -2.53;
    const cube field = grid(21, h, o, [&](int x, int y, int z) {
        const double px = o + h * x, py = o + h * y, pz = o + h * z;
        return R * R - (px * px + py * py + pz * pz);
        });
    const std::vector<Triangle> tri = marchingCubes(field, 0.0);
    ASSERT_GT(tri.size(), 500u);
    double area = 0.0, volume = 0.0;
    using node = std::array<long long, 3>;
    auto key = [](const d3& p) { return node{ std::llround(p[0] * 1e8), std::llround(p[1] * 1e8), std::llround(p[2] * 1e8) }; };
    std::map<std::pair<node, node>, int> edges;
    for (const Triangle& t : tri)
    {
        area += t.calc_area();
        volume += t.calc_inner_volume();
        for (int v = 1; v <= 3; v++)
        {
            EXPECT_NEAR(array_length(t.get_v(v)), R, 0.02);
            edges[{ key(t.get_v(v)), key(t.get_v(v % 3 + 1)) }]++;
        }
        const d3 c = t.calc_center();
        EXPECT_NEAR(array_length(c), R, 0.03);
    }
    int unmatched = 0;
    for (const auto& [e, n] : edges)
        if (n != 1 || edges.count({ e.second, e.first }) != 1)
            unmatched++;
    EXPECT_EQ(unmatched, 0) << "of " << edges.size() << " directed edges";
    EXPECT_NEAR(area, 4 * constants::PI * R * R, 0.08 * 4 * constants::PI * R * R);
    EXPECT_NEAR(std::abs(volume), 4.0 / 3.0 * constants::PI * R * R * R, 0.03 * 4.0 / 3.0 * constants::PI * R * R * R);
    EXPECT_EQ(Triangle({ 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }).calc_area(), 0.5);
    EXPECT_TRUE(marchingCubes(cube({ 1, 3, 3 }, 0, true), 0.0).empty());
    EXPECT_TRUE(marchingCubes(cube({ 3, 3, 1 }, 0, true), 0.0).empty());
}

//an edge whose corner values differ by less than 1e-12 is cut at its midpoint: one corner at iso, the rest 1e-13 below
TEST(CubeIsoTests, DegenerateEdgeIsCutAtTheMidpoint)
{
    cube c = grid(2, 1.0, 0.0, [](int, int, int) { return 0.5 - 1e-13; });
    c.set_value(0, 0, 0, 0.5);
    const std::vector<Triangle> tri = marchingCubes(c, 0.5);
    ASSERT_EQ(tri.size(), 1u);
    std::set<std::array<double, 3>> got, expected{ { 0.5, 0.0, 0.0 }, { 0.0, 0.5, 0.0 }, { 0.0, 0.0, 0.5 } };
    for (int v = 1; v <= 3; v++)
        got.insert(tri[0].get_v(v));
    EXPECT_EQ(got, expected);
    EXPECT_EQ(tri[0].get_v(4), (d3{ 0, 0, 0 }));
}

//mix_colour: below low is colour 0, above high colour 2, linear to colour 1 at the mid point, then clamped to 0..255
TEST(CubeIsoTests, MixColourRamp)
{
    const std::array<std::array<int, 3>, 3> code{ { { 0, 0, 0 }, { 100, 100, 100 }, { 200, 200, 200 } } };
    EXPECT_EQ(mix_colour(-1.0, code, 0.0, 2.0), (RGB{ 0, 0, 0 }));
    EXPECT_EQ(mix_colour(3.0, code, 0.0, 2.0), (RGB{ 200, 200, 200 }));
    EXPECT_EQ(mix_colour(0.5, code, 0.0, 2.0), (RGB{ 50, 50, 50 }));
    EXPECT_EQ(mix_colour(1.5, code, 0.0, 2.0), (RGB{ 150, 150, 150 }));
    EXPECT_EQ(mix_colour(1.0, code, 0.0, 2.0), (RGB{ 100, 100, 100 }));
    const std::array<std::array<int, 3>, 3> wild{ { { -100, 0, 0 }, { 100, 100, 100 }, { 300, -5, 255 } } };
    EXPECT_EQ(mix_colour(3.0, wild, 0.0, 2.0), (RGB{ 255, 0, 255 }));
    EXPECT_EQ(mix_colour(0.25, wild, 0.0, 2.0), (RGB{ 0, 25, 25 }));
    EXPECT_EQ(mtl_name({ 1, 22, 255 }), "FaceMaterial_1_22_255");
}

//get_colour from a cube interpolates at the face centre and mixes towards colour 1 by val/low or val/high
TEST(CubeIsoTests, ColourFromAField)
{
    const cube field = grid(3, 1.0, 0.0, [](int x, int, int) { return x - 1.0; });
    const std::array<std::array<int, 3>, 3> code{ { { 0, 0, 0 }, { 100, 100, 100 }, { 200, 200, 200 } } };
    auto at = [](double x) { return Triangle({ x - 0.5, 0.5, 0.5 }, { x + 0.5, 0.5, 0.5 }, { x, 0.5, 0.5 }); };
    Triangle t = at(0.25);
    get_colour(t, field, code, -1.0, 1.0);
    EXPECT_EQ(t.get_colour(), (RGB{ 25, 25, 25 })) << "val -0.75: 0.75 of colour 0 and 0.25 of colour 1";
    t = at(1.75);
    get_colour(t, field, code, -1.0, 1.0);
    EXPECT_EQ(t.get_colour(), (RGB{ 175, 175, 175 })) << "val 0.75: 0.75 of colour 2 and 0.25 of colour 1";
    t = at(1.0);
    get_colour(t, field, code, -1.0, 1.0);
    EXPECT_EQ(t.get_colour(), (RGB{ 100, 100, 100 })) << "val 0 is the middle colour";
    t = at(0.25);
    get_colour(t, field, code, -0.5, 0.5);
    EXPECT_EQ(t.get_colour(), (RGB{ 0, 0, 0 })) << "below low";
    t = at(1.75);
    get_colour(t, field, code, -0.5, 0.5);
    EXPECT_EQ(t.get_colour(), (RGB{ 200, 200, 200 })) << "above high";
    t = at(-5.0);
    get_colour(t, field, code, -1.0, 1.0);
    EXPECT_EQ(t.get_colour(), (RGB{ 100, 100, 100 })) << "outside the box the field reads 0";
}

//get_colour from a point function: calc_d_i is the distance to the nearest atom, ramped through mix_colour
TEST(CubeIsoTests, ColourFromTheNearestAtomDistance)
{
    WFN h = one_hydrogen("colour");
    const std::array<std::array<int, 3>, 3> code{ { { 0, 0, 0 }, { 100, 100, 100 }, { 200, 200, 200 } } };
    auto at = [](double x) { return Triangle({ x, 1.0, 0.0 }, { x, -1.0, 0.0 }, { x, 0.0, 0.0 }); };
    EXPECT_NEAR(calc_d_i({ 3.0, 4.0, 0.0 }, h), 5.0, 1e-12);
    Triangle t = at(1.0);
    get_colour(t, calc_d_i, h, code, 0.0, 2.0);
    EXPECT_EQ(t.get_colour(), (RGB{ 100, 100, 100 }));
    t = at(0.5);
    get_colour(t, calc_d_i, h, code, 0.0, 2.0);
    EXPECT_EQ(t.get_colour(), (RGB{ 50, 50, 50 }));
    t = at(1.5);
    get_colour(t, calc_d_i, h, code, 0.0, 2.0);
    EXPECT_EQ(t.get_colour(), (RGB{ 150, 150, 150 }));
    t.set_colour_index(size_t(4));
    EXPECT_EQ(t.get_colour_index(), 4);
}

//colour_by_ESP from a given ESP vector: the most negative face is red, the most positive blue, zero white, and the log names the range
TEST(CubeIsoTests, ColourByEspVector)
{
    std::vector<Triangle> tri(3, Triangle({ 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }));
    const vec esp = surface_ESP(tri, [](const d3& p) { return p[0] + p[1]; });
    ASSERT_EQ(esp.size(), 3u);
    EXPECT_NEAR(esp[0], 2.0 / 3.0, 1e-12);
    std::ostringstream log;
    colour_by_ESP(tri, vec{ -0.04, 0.0, 0.02 }, log);
    EXPECT_EQ(tri[0].get_colour(), (RGB{ 255, 0, 0 }));
    EXPECT_EQ(tri[1].get_colour(), (RGB{ 255, 255, 255 }));
    EXPECT_EQ(tri[2].get_colour(), (RGB{ 127, 127, 255 }));
    EXPECT_NE(log.str().find("ESP on the surface from -0.04 to 0.02 au, coloured red (-0.04) white (0) blue (+0.04)"), std::string::npos) << log.str();
}

//subdivideCube samples (level + 1)^3 points of the box p1..p2 with p1 first and p2 last
TEST(CubeIsoTests, SubdivideCubeSamplesTheBox)
{
    const std::vector<d3> pts = subdivideCube({ 0.0, 0.0, 0.0 }, { 1.0, 2.0, 4.0 }, 2);
    ASSERT_EQ(pts.size(), 27u);
    EXPECT_EQ(pts.front(), (d3{ 0.0, 0.0, 0.0 }));
    EXPECT_EQ(pts.back(), (d3{ 1.0, 2.0, 4.0 }));
    EXPECT_EQ(pts[13], (d3{ 0.5, 1.0, 2.0 }));
    EXPECT_EQ(pts[1], (d3{ 0.0, 0.0, 2.0 })) << "z runs fastest";
}

//writeObj lists three v lines per triangle and one f line per triangle with 1-based consecutive indices; an unwritable path fails
TEST(CubeIoTests, ObjWriter)
{
    const std::vector<Triangle> tri{ Triangle({ 0, 0, 0 }, { 1, 0, 0 }, { 0, 0.5, 0 }), Triangle({ 2, 2, 2 }, { 3, 2, 2 }, { 2, 3, 2 }) };
    const auto obj = tmp("mesh.obj");
    testing::internal::CaptureStdout();
    ASSERT_TRUE(writeObj(obj, tri));
    EXPECT_NE(testing::internal::GetCapturedStdout().find("OBJ file written to"), std::string::npos);
    const std::vector<std::string> l = lines_of(obj);
    ASSERT_EQ(l.size(), 8u);
    EXPECT_EQ(l[0], "v 0 0 0");
    EXPECT_EQ(l[2], "v 0 0.5 0");
    EXPECT_EQ(l[5], "v 2 3 2");
    EXPECT_EQ(l[6], "f 1 2 3");
    EXPECT_EQ(l[7], "f 4 5 6");
    std::filesystem::remove(obj);
    EXPECT_FALSE(writeObj(tmp("no_such_dir") / "mesh.obj", tri));
}

//writeColourObj adds a mtllib line, a usemtl line per face and one newmtl block per distinct colour in the mtl next to the obj
TEST(CubeIoTests, ColourObjWriterAndMtl)
{
    std::vector<Triangle> tri{ Triangle({ 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 255, 0, 0 }),
                               Triangle({ 0, 0, 1 }, { 1, 0, 1 }, { 0, 1, 1 }, { 0, 0, 255 }),
                               Triangle({ 0, 0, 2 }, { 1, 0, 2 }, { 0, 1, 2 }, { 255, 0, 0 }) };
    const auto dir = tmp("colourobj_dir");
    std::filesystem::create_directories(dir);
    const auto obj = dir / "mesh.obj", mtl = dir / "mesh.mtl";
    testing::internal::CaptureStdout();
    ASSERT_TRUE(writeColourObj(obj, tri));
    const std::string out = testing::internal::GetCapturedStdout();
    EXPECT_NE(out.find("MTL file written to " + mtl.string()), std::string::npos) << out;
    EXPECT_NE(out.find("OBJ file written to"), std::string::npos);
    const std::vector<std::string> l = lines_of(obj);
    ASSERT_EQ(l.size(), 16u);
    EXPECT_EQ(l[0], "mtllib mesh.mtl");
    EXPECT_EQ(l[10], "usemtl FaceMaterial_255_0_0");
    EXPECT_EQ(l[11], "f 1 2 3");
    EXPECT_EQ(l[12], "usemtl FaceMaterial_0_0_255");
    EXPECT_EQ(l[14], "usemtl FaceMaterial_255_0_0");
    EXPECT_EQ(l[15], "f 7 8 9");
    ASSERT_TRUE(std::filesystem::exists(mtl));
    const std::vector<std::string> m = lines_of(mtl);
    int newmtl = 0;
    for (const std::string& s : m)
        if (s.rfind("newmtl ", 0) == 0) newmtl++;
    EXPECT_EQ(newmtl, 2) << "the repeated red is one material";
    EXPECT_NE(std::find(m.begin(), m.end(), "newmtl FaceMaterial_255_0_0"), m.end());
    EXPECT_NE(std::find(m.begin(), m.end(), "Kd 1 0 0"), m.end());
    EXPECT_NE(std::find(m.begin(), m.end(), "Kd 0 0 1"), m.end());
    std::filesystem::remove_all(dir);
    EXPECT_FALSE(writeMTL((tmp("no_such_dir") / "x.mtl").string(), tri));
}

namespace
{
    template <typename T>
    void npy_round_trip(const std::string& tag)
    {
        npy::npy_data<T> d;
        d.shape = { 2, 3 };
        d.fortran_order = false;
        for (int k = 0; k < 6; k++)
            d.data.push_back(T(k + 1));
        const auto file = tmp("npy_" + tag + ".npy");
        npy::write_npy<T>(file.string(), d);
        const npy::npy_data<T> r = npy::read_npy<T>(file);
        EXPECT_EQ(r.shape, d.shape) << tag;
        EXPECT_FALSE(r.fortran_order) << tag;
        EXPECT_TRUE(r.data == d.data) << tag;
        std::filesystem::remove(file);
    }
}

//every dtype in dtype_map survives a write_npy / read_npy round trip through a file
TEST(CubeNpyTests, EveryDtypeRoundTrips)
{
    npy_round_trip<float>("float");
    npy_round_trip<double>("double");
    npy_round_trip<long double>("longdouble");
    npy_round_trip<char>("char");
    npy_round_trip<signed char>("schar");
    npy_round_trip<short>("short");
    npy_round_trip<int>("int");
    npy_round_trip<long>("long");
    npy_round_trip<long long>("longlong");
    npy_round_trip<unsigned char>("uchar");
    npy_round_trip<unsigned short>("ushort");
    npy_round_trip<unsigned int>("uint");
    npy_round_trip<unsigned long>("ulong");
    npy_round_trip<unsigned long long>("ulonglong");
    npy_round_trip<std::complex<float>>("cfloat");
    npy_round_trip<std::complex<double>>("cdouble");
    npy_round_trip<std::complex<long double>>("clongdouble");
}

//the old SaveArrayAsNumpy / LoadArrayFromNumpy interface keeps shape, fortran order and data (pointer and vector overloads)
TEST(CubeNpyTests, OldInterfaceRoundTrip)
{
    const unsigned long shape[2] = { 2, 3 };
    const std::vector<double> data{ 1.5, 2.5, 3.5, 4.5, 5.5, 6.5 };
    const auto file = tmp("npy_old.npy");
    npy::SaveArrayAsNumpy<double>(file, true, 2, shape, data);
    std::vector<unsigned long> got_shape;
    bool fortran = false;
    std::vector<double> got{ -1.0 };
    npy::LoadArrayFromNumpy<double>(file, got_shape, fortran, got);
    EXPECT_EQ(got_shape, (std::vector<unsigned long>{ 2, 3 }));
    EXPECT_TRUE(fortran);
    EXPECT_EQ(got, (std::vector<double>{ -1.0, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5 })) << "loading appends";

    const int ints[3] = { 7, 8, 9 };
    const unsigned long one[1] = { 3 };
    npy::SaveArrayAsNumpy<int>(file, false, 1, one, ints);
    std::vector<int> got_ints;
    npy::LoadArrayFromNumpy<int>(file, got_shape, got_ints);
    EXPECT_EQ(got_shape, (std::vector<unsigned long>{ 3 }));
    EXPECT_EQ(got_ints, (std::vector<int>{ 7, 8, 9 }));
    std::filesystem::remove(file);
}

//the header is magic, version 1.0, a little-endian length and a python dict padded to a multiple of 16 bytes
TEST(CubeNpyTests, HeaderBytes)
{
    npy::npy_data<double> d;
    d.shape = { 2, 3 };
    d.data.assign(6, 0.0);
    d.data[0] = 1.0;
    std::ostringstream out(std::ios::binary);
    npy::write_npy<double>(out, d);
    const std::string s = out.str();
    ASSERT_GT(s.size(), 10u);
    EXPECT_EQ(s.substr(0, 6), std::string("\x93NUMPY", 6));
    EXPECT_EQ(s[6], 1);
    EXPECT_EQ(s[7], 0);
    const size_t header_len = (unsigned char)s[8] | ((unsigned char)s[9] << 8);
    EXPECT_EQ((10 + header_len) % 16, 0u);
    const std::string dict = s.substr(10, header_len);
    EXPECT_EQ(dict.back(), '\n');
    const std::string descr = npy::big_endian ? ">f8" : "<f8";
    EXPECT_EQ(dict.substr(0, dict.find_last_not_of(" \n") + 1), "{'descr': '" + descr + "', 'fortran_order': False, 'shape': (2, 3), }");
    EXPECT_EQ(s.size(), 10 + header_len + 6 * sizeof(double));
    const std::string one = npy::big_endian ? "\x3F\xF0" + std::string(6, '\0') : std::string(6, '\0') + "\xF0\x3F";
    EXPECT_EQ(s.substr(10 + header_len, 8), one) << "1.0 is the IEEE bytes 3FF0000000000000 in the declared byte order";
    EXPECT_EQ(npy::pyparse::write_tuple(std::vector<unsigned long>{ 5 }), "(5,)");
    EXPECT_EQ(npy::pyparse::write_tuple(std::vector<unsigned long>{}), "()");
    EXPECT_EQ(npy::pyparse::write_boolean(true), "True");
    EXPECT_EQ(npy::dtype_map.at(std::type_index(typeid(unsigned char))).str(), "|u1");
    EXPECT_EQ(npy::dtype_map.at(std::type_index(typeid(std::complex<double>))).str(), descr.substr(0, 1) + "c16");
}

//a header dict of 65025 bytes or more switches to format 2.0 with a four-byte length, and reads back
TEST(CubeNpyTests, LongHeaderUsesVersionTwo)
{
    npy::npy_data<float> d;
    d.shape.assign(22000, 1);
    d.data = { 3.5f };
    std::stringstream io(std::ios::in | std::ios::out | std::ios::binary);
    npy::write_npy<float>(io, d);
    const std::string s = io.str();
    EXPECT_EQ(s[6], 2);
    EXPECT_EQ(s[7], 0);
    const size_t header_len = (unsigned char)s[8] | ((unsigned char)s[9] << 8) | ((unsigned char)s[10] << 16) | ((unsigned char)s[11] << 24);
    EXPECT_EQ((12 + header_len) % 16, 0u);
    EXPECT_GE(header_len, 65000u);
    const npy::npy_data<float> r = npy::read_npy<float>(io);
    EXPECT_EQ(r.shape.size(), 22000u);
    ASSERT_EQ(r.data.size(), 1u);
    EXPECT_EQ(r.data[0], 3.5f);
}

//every io failure throws a runtime_error with the message that names it
TEST(CubeNpyTests, IoErrorsThrow)
{
    const auto missing = tmp("no_such_dir") / "x.npy";
    EXPECT_EQ(what_throws([&] { npy::read_npy<double>(missing); }), "io error: failed to open a file.");
    npy::npy_data<double> d;
    d.shape = { 1 };
    d.data = { 1.0 };
    EXPECT_EQ(what_throws([&] { npy::write_npy<double>(missing.string(), d); }), "io error: failed to open a file.");
    const npy::npy_data_ptr<double> p{ d.data.data(), { 1 }, false };
    EXPECT_EQ(what_throws([&] { npy::write_npy<double>(missing, p); }), "io error: failed to open a file.");

    std::istringstream empty;
    EXPECT_EQ(what_throws([&] { npy::read_npy<double>(empty); }), "io error: failed reading file");
    std::istringstream junk("XXXXXXXXXXXXXXXX");
    EXPECT_EQ(what_throws([&] { npy::read_npy<double>(junk); }), "this file does not have a valid npy format.");
    std::istringstream v3(std::string("\x93NUMPY\x03\x00\x10\x00", 10) + std::string(16, ' '));
    EXPECT_EQ(what_throws([&] { npy::read_npy<double>(v3); }), "unsupported file format version");

    std::stringstream io(std::ios::in | std::ios::out | std::ios::binary);
    npy::write_npy<double>(io, d);
    testing::internal::CaptureStdout();
    EXPECT_EQ(what_throws([&] { npy::read_npy<float>(io); }), "formatting error: typestrings not matching");
    const std::string out = testing::internal::GetCapturedStdout();
    EXPECT_NE(out.find("header.dtype: "), std::string::npos);
    EXPECT_NE(out.find("f8"), std::string::npos);
    EXPECT_NE(out.find("f4"), std::string::npos);
}

//the python literal parsers accept the numpy header forms and reject everything else with a message
TEST(CubeNpyTests, PythonLiteralParsers)
{
    const npy::dtype_t f8 = npy::parse_descr("<f8");
    EXPECT_EQ(f8.byteorder, '<');
    EXPECT_EQ(f8.kind, 'f');
    EXPECT_EQ(f8.itemsize, 8u);
    EXPECT_EQ(npy::parse_descr("|u1").str(), "|u1");
    EXPECT_EQ(npy::parse_descr(">c16").itemsize, 16u);
    EXPECT_EQ(what_throws([] { npy::parse_descr("f8"); }), "invalid typestring (length)");
    EXPECT_EQ(what_throws([] { npy::parse_descr("xf8"); }), "invalid typestring (byteorder)");
    EXPECT_EQ(what_throws([] { npy::parse_descr("<x8"); }), "invalid typestring (kind)");
    EXPECT_EQ(what_throws([] { npy::parse_descr("<f8a"); }), "invalid typestring (itemsize)");

    const auto dict = npy::pyparse::parse_dict("  {'descr': '<f8', 'fortran_order': False, 'shape': (2, 3), } ", { "descr", "fortran_order", "shape" });
    ASSERT_EQ(dict.size(), 3u);
    EXPECT_EQ(dict.at("descr"), "'<f8'");
    EXPECT_EQ(dict.at("fortran_order"), "False");
    EXPECT_EQ(dict.at("shape"), "(2, 3)");
    EXPECT_TRUE(npy::pyparse::parse_dict("{'a': 1}", {}).empty());
    EXPECT_EQ(what_throws([] { npy::pyparse::parse_dict("'a': 1", { "a" }); }), "Not a Python dictionary.");
    EXPECT_EQ(what_throws([] { npy::pyparse::parse_dict("{'a': 1}", { "b" }); }), "Missing 'b' key.");

    EXPECT_TRUE(npy::pyparse::parse_bool("True"));
    EXPECT_FALSE(npy::pyparse::parse_bool("False"));
    EXPECT_EQ(what_throws([] { npy::pyparse::parse_bool("maybe"); }), "Invalid python boolan.");
    EXPECT_EQ(npy::pyparse::parse_str("'abc'"), "abc");
    EXPECT_EQ(what_throws([] { npy::pyparse::parse_str("abc"); }), "Invalid python string.");
    EXPECT_EQ(npy::pyparse::parse_tuple(" (1, 2) "), (std::vector<std::string>{ "1", " 2" }));
    EXPECT_EQ(npy::pyparse::parse_tuple("(7,)"), (std::vector<std::string>{ "7" }));
    EXPECT_EQ(what_throws([] { npy::pyparse::parse_tuple("1, 2"); }), "Invalid Python tuple.");

    const npy::header_t h = npy::parse_header("{'descr': '<i4', 'fortran_order': True, 'shape': (4,), }      \n");
    EXPECT_EQ(h.dtype.str(), "<i4");
    EXPECT_TRUE(h.fortran_order);
    EXPECT_EQ(h.shape, (npy::shape_t{ 4 }));
    EXPECT_EQ(what_throws([] { npy::parse_header("{'descr': '<i4', 'fortran_order': True, 'shape': (4,), }"); }), "invalid header");
}
