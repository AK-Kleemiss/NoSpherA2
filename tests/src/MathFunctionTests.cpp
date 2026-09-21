
#include "pch.h"

#include "core/convenience.h"
#include "core/constants.h"
#include "core/fchk.h"
#include "core/AtomGrid.h"
#include "core/SALTED_utilities.h"
#include "core/scattering_factors.h"
#include "core/nos_math.h"
#include "core/GridManager.h"
#include "core/atoms.h"
#include "core/tsc_block.h"
#include "core/cell.h"
#include "core/wfn_class.h"
#include "core/properties.h"
#include "core/integrator.h"
#include "core/integration_params.h"
#include "core/basis_set.h"
#include "core/geometry_aid.h"
#include "core/crystal_energies.h"
#include "core/NoSpherA2.h"
#include "core/isosurface.h"
#include "core/npy.h"
#include "core/libCintMain.h"
#include "core/throughput.h"
#include "core/i_tensor_stream.h"
#include "core/tsc_label_converter.h"
#include "core/sphere_lebedev_rule.h"
#include <occ/qm/hf.h>
#include <occ/qm/scf.h>
#include <spdlog/spdlog.h>
#undef I
#ifdef NOSPHERA2_USE_GPU
#include "core/blas_gpu.h"
#include "core/aux_density_gpu.h"
#endif

static constexpr double PI_VAL = 3.14159265358979323846;

namespace {
    template <typename MatType1, typename MatType2>
    void compare_matrices(const MatType1& mat, const MatType2& vecMat) {
        const double tol = 1e-9;
        for (size_t i = 0; i < mat.extent(0); i++) {
            for (size_t j = 0; j < mat.extent(1); j++) {
                //Compile different versions based on the type of the matrix dMatrix and dMatrix2 are accessed with mat(i, j) and vecMat[i][j] respectively
                if constexpr (std::is_same_v<MatType1, dMatrix2> && std::is_same_v<MatType2, vec2>) {
                    EXPECT_NEAR(mat(i, j), vecMat[i][j], tol);
                }
                else if constexpr (std::is_same_v<MatType1, cMatrix2> && std::is_same_v<MatType2, cvec2>) {
                    EXPECT_NEAR(mat(i, j).real(), vecMat[i][j].real(), tol);
                    EXPECT_NEAR(mat(i, j).imag(), vecMat[i][j].imag(), tol);
                }
            }
        }
    }

    void test_solve_linear_equations() {
        // Small, well-conditioned 3x3 test (precomputed)
        // A * x_expected = b
        const vec2 A = {
            {3.0,  2.0, -1.0},
            {2.0, -2.0,  4.0},
            {-1.0, 0.5, -1.0}
        };

        // Precomputed true solution
        const vec x_expected = { 1.0, -2.0, -2.0 };

        // Build right-hand side b = A * x_expected
        vec b(3, 0.0);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                b[i] += A[i][j] * x_expected[j];

        // Keep original RHS for logging if the test fails
        vec b_orig = b;

        // Call the solver -- it replaces b with the solution x (in-place)
        solve_linear_system(A, b);

        // Validate result
        const double tol = 1e-9;
        bool ok = true;
        for (int i = 0; i < 3; ++i) {
            EXPECT_NEAR(b[i], x_expected[i], tol);
        }

        if (ok) {
            std::cout << "test_solve_linear_equations: PASSED\n";
        }
        else {
            std::cout << "test_solve_linear_equations: FAILED\n";
            std::cout << "Matrix A:\n";
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j)
                    std::cout << std::setw(12) << A[i][j];
                std::cout << std::endl;
            }
            std::cout << "Target   b: " << b_orig[0] << " " << b_orig[1] << " " << b_orig[2] << std::endl;
            std::cout << "Original x: ";
            for (double v : x_expected) std::cout << v << " ";
            std::cout << std::endl;
            std::cout << "Returned x: ";
            for (double v : b) std::cout << v << " ";
            std::cout << std::endl;
        }
    }

    void test_openblas()
    {
        ivec dims = { 10, 10 };
        // Init Mat A with some values as a 3x3 matrix
        vec2 A(dims[0], vec(dims[1]));
        vec2 B(dims[0], vec(dims[1]));
        // Init A and B with random values between -100 and 100
        for (int i = 0; i < dims[0]; i++)
        {
            for (int j = 0; j < dims[1]; j++)
            {
                A[i][j] = rand() % 200 - 100;
                B[i][j] = rand() % 200 - 100;
            }
        }

        vec fA = flatten<double>(A);
        vec fB = flatten<double>(B);

        dMatrix2 matA(dims[0], dims[1]);
        std::copy(fA.data(), fA.data() + fA.size(), matA.data());
        dMatrix2 matB(dims[0], dims[1]);
        std::copy(fB.data(), fB.data() + fB.size(), matB.data());
        // Init Mat A and Mat B as 3x3 matrices

        std::cout << "Testing matrices directly\n";
        compare_matrices(matA, A);
        compare_matrices(matB, B);

        std::cout << "Testing untransposed matrices\n";
        // First test regular dot-product
        compare_matrices(dot(matA, matB, false, false), self_dot(A, B));

        std::cout << "Testing transpose A\n";
        ////Second compare first transpose
        compare_matrices(dot(matA, matB, true, false), self_dot(transpose(A), B));

        std::cout << "Testing transpose B\n";
        ////Third comparte second transpose
        compare_matrices(dot(matA, matB, false, true), self_dot(A, transpose(B)));

        std::cout << "Testing transpose A and B\n";
        ////Fourth compare both transposed
        compare_matrices(dot(matA, matB, true, true), self_dot(transpose(A), transpose(B)));

        // Init Complex matrices
        cvec2 C(dims[0], cvec(dims[1])), D(dims[0], cvec(dims[1]));
        for (int i = 0; i < dims[0]; i++)
        {
            for (int j = 0; j < dims[1]; j++)
            {
                C[i][j] = cdouble(rand() % 200 - 100, rand() % 200 - 100);
                D[i][j] = cdouble(rand() % 200 - 100, rand() % 200 - 100);
            }
        }

        cvec fC = flatten<cdouble>(C);
        cvec fD = flatten<cdouble>(D);
        cMatrix2 matC(dims[0], dims[1]);
        std::copy(fC.data(), fC.data() + fC.size(), matC.data());
        cMatrix2 matD(dims[0], dims[1]);
        std::copy(fD.data(), fD.data() + fD.size(), matD.data());

        std::cout << "Testing C-matrices directly\n";
        compare_matrices(matC, C);
        compare_matrices(matD, D);

        std::cout << "Testing untransposed C-matrices\n";
        // First test regular dot-product
        compare_matrices(dot(matC, matD, false, false), self_dot(C, D));

        std::cout << "Testing transpose C\n";
        ////Second compare first transpose
        compare_matrices(dot(matC, matD, true, false), self_dot(transpose(C), D));

        std::cout << "Testing transpose D\n";
        ////Third comparte second transpose
        compare_matrices(dot(matC, matD, false, true), self_dot(C, transpose(D)));

        std::cout << "Testing transpose C and D\n";
        ////Fourth compare both transposed
        compare_matrices(dot(matC, matD, true, true), self_dot(transpose(C), transpose(D)));

        // Test 2D x 1D matrix multiplication
        dims[0] = 12;
        vec E(dims[1]);
        vec2 F(dims[0], vec(dims[1]));
        for (int i = 0; i < dims[1]; i++)
        {
            E[i] = rand() % 200 - 100;
            for (int j = 0; j < dims[0]; j++)
            {
                F[j][i] = rand() % 200 - 100;
            }
        }
        vec fE = flatten<double>(F);
        dMatrix1 matE(dims[1]);
        dMatrix2 matF(dims[0], dims[1]);
        std::copy(E.data(), E.data() + E.size(), matE.data());
        std::copy(fE.data(), fE.data() + fE.size(), matF.data());

        // For matrix just reuse matA
        std::cout << "Testing 2D x 1D matrix multiplication\n";
        compare_matrices(dot(matF, matE, false), self_dot(F, E));

        std::cout << "All BLAS tests passed!\n";
    }
}

namespace NoSpherA2UnitTests
{
    // -----------------------------------------------------------------------

    TEST(GeometryTests, ArrayDistance_UnitCubeDiagonal)
    {
        double d = array_length(d3{ 0.0, 0.0, 0.0 }, d3{ 1.0, 1.0, 1.0 });
        EXPECT_NEAR(std::sqrt(3.0), d, 1e-12);
    }

    TEST(GeometryTests, ArrayDistance_SamePoint)
    {
        EXPECT_NEAR(0.0, array_length(d3{ 1.5, 2.3, -4.7 }, d3{ 1.5, 2.3, -4.7 }), 1e-12);
    }

    TEST(GeometryTests, VecDiff_Basic)
    {
        d3 r = vec_diff({ 5.0, 10.0, 15.0 }, { 1.0, 2.0, 3.0 });
        EXPECT_NEAR(4.0, r[0], 1e-12);
        EXPECT_NEAR(8.0, r[1], 1e-12);
        EXPECT_NEAR(12.0, r[2], 1e-12);
    }

    TEST(GeometryTests, VecDiff_ZeroResult)
    {
        d3 r = vec_diff({ 1.0, 2.0, 3.0 }, { 1.0, 2.0, 3.0 });
        EXPECT_NEAR(0.0, r[0], 1e-12);
        EXPECT_NEAR(0.0, r[1], 1e-12);
        EXPECT_NEAR(0.0, r[2], 1e-12);
    }

    TEST(GeometryTests, VecCross_BasisVectors)
    {
        // i x j = k
        d3 r = vec_cross({ 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 });
        EXPECT_NEAR(0.0, r[0], 1e-12);
        EXPECT_NEAR(0.0, r[1], 1e-12);
        EXPECT_NEAR(1.0, r[2], 1e-12);
    }

    TEST(GeometryTests, VecCross_AntiCommutative)
    {
        // a x b = -(b x a)
        d3 ab = vec_cross({ 1.0, 2.0, 3.0 }, { 4.0, 5.0, 6.0 });
        d3 ba = vec_cross({ 4.0, 5.0, 6.0 }, { 1.0, 2.0, 3.0 });
        for (int i = 0; i < 3; ++i) {
            EXPECT_NEAR(-ab[i], ba[i], 1e-12);
        }
    }

    TEST(GeometryTests, VecCross_ParallelVectors)
    {
        // Parallel vectors → zero cross product
        d3 r = vec_cross({ 2.0, 4.0, 6.0 }, { 1.0, 2.0, 3.0 });
        EXPECT_NEAR(0.0, r[0], 1e-12);
        EXPECT_NEAR(0.0, r[1], 1e-12);
        EXPECT_NEAR(0.0, r[2], 1e-12);
    }

    TEST(GeometryTests, VecDot_Perpendicular)
    {
        EXPECT_NEAR(0.0, vec_dot({ 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }), 1e-12);
    }

    TEST(GeometryTests, VecDot_Parallel)
    {
        EXPECT_NEAR(3.0, vec_dot({ 1.0, 1.0, 1.0 }, { 1.0, 1.0, 1.0 }), 1e-12);
    }

    TEST(GeometryTests, VecDot_KnownValue)
    {
        EXPECT_NEAR(32.0, vec_dot({ 1.0, 2.0, 3.0 }, { 4.0, 5.0, 6.0 }), 1e-12);
    }

    // -----------------------------------------------------------------------

    TEST(NumericTests, IsSimilarRel_WithinTolerance)
    {
        // 1.0 vs 1.005: relative diff ≈ 0.5 %, within 1 %
        EXPECT_TRUE(is_similar_rel(1.0, 1.005, 0.01));
    }

    TEST(NumericTests, IsSimilarRel_OutsideTolerance)
    {
        // 1.0 vs 1.05: relative diff = 5 %, outside 1 %
        EXPECT_FALSE(is_similar_rel(1.0, 1.05, 0.01));
    }

    TEST(NumericTests, IsSimilarRel_EqualValues)
    {
        EXPECT_TRUE(is_similar_rel(42.0, 42.0, 1e-6));
    }

    TEST(NumericTests, IsSimilarAbs_WithinTolerance)
    {
        EXPECT_TRUE(is_similar_abs(1.0, 1.0009, 0.001));
    }

    TEST(NumericTests, IsSimilarAbs_OutsideTolerance)
    {
        EXPECT_FALSE(is_similar_abs(1.0, 1.002, 0.001));
    }

    TEST(NumericTests, FastExpNeg_AtZero)
    {
        // fast_exp_neg uses std::exp for x > -ln(2) ≈ -0.693
        double val = fast_exp_neg(0.0);
        EXPECT_NEAR(1.0, val, 1e-12);
    }

    TEST(NumericTests, FastExpNeg_AtMinusOne)
    {
        double approx = fast_exp_neg(-1.0);
        double exact = std::exp(-1.0);
        // Approximation tolerance: 0.5 %
        EXPECT_NEAR(exact, approx, exact * 0.005);
    }

    TEST(NumericTests, FastExpNeg_AtMinusTen)
    {
        // (1 + x/1024)^1024 underestimates exp(x) by ~x²/(2N) ≈ 4.9% at x=-10.
        double approx = fast_exp_neg(-10.0);
        double exact = std::exp(-10.0);
        EXPECT_NEAR(exact, approx, exact * 0.06);
    }

    TEST(NumericTests, FastExpNeg_BeyondCutoff)
    {
        // Values below -42 must return 0 exactly
        EXPECT_NEAR(0.0, fast_exp_neg(-100.0), 1e-30);
    }

    TEST(NumericTests, BLAS_tests)
    {
        test_openblas();
    }

    TEST(NumericTests, SolveLinearEquations)
    {
        test_solve_linear_equations();
    }

    // -----------------------------------------------------------------------

    TEST(ConstantsTests, ConstAbs_Positive) { EXPECT_EQ(5, constants::const_abs(5)); }
    TEST(ConstantsTests, ConstAbs_Negative) { EXPECT_EQ(5, constants::const_abs(-5)); }
    TEST(ConstantsTests, ConstAbs_Zero) { EXPECT_EQ(0, constants::const_abs(0)); }

    TEST(ConstantsTests, ConstexprPow_Integer)
    {
        EXPECT_NEAR(1024.0, constants::constexpr_pow(2.0, 10), 1e-12);
        EXPECT_NEAR(1.0, constants::constexpr_pow(10.0, 0), 1e-12);
        EXPECT_NEAR(8.0, constants::constexpr_pow(2.0, 3), 1e-12);
    }

    TEST(ConstantsTests, ConstantsSqrt_KnownValues)
    {
        EXPECT_NEAR(2.0, constants::sqrt(4.0), 1e-10);
        EXPECT_NEAR(std::sqrt(2.0), constants::sqrt(2.0), 1e-10);
        EXPECT_NEAR(0.0, constants::sqrt(0.0), 1e-12);
    }

    TEST(ConstantsTests, ConstantsSqrt_NegativeIsNaN)
    {
#if defined(__FAST_MATH__)
    GTEST_SKIP() << "Skipping NaN test under -ffast-math";
#else
    double r = constants::sqrt(-1.0);
    EXPECT_TRUE(std::isnan(r));
#endif
    }

    TEST(ConstantsTests, ExpApprox_AtZero)
    {
        EXPECT_NEAR(1.0, constants::exp_approx(0.0, 25), 1e-12);
    }

    TEST(ConstantsTests, ExpApprox_AtOne)
    {
        // 25-term Taylor series matches std::exp to machine precision
        EXPECT_NEAR(std::exp(1.0), constants::exp_approx(1.0, 25), 1e-12);
    }

    TEST(ConstantsTests, ExpApprox_AtMinusOne)
    {
        EXPECT_NEAR(std::exp(-1.0), constants::exp_approx(-1.0, 25), 1e-10);
    }

    TEST(ConstantsTests, LogApprox_AtOne)
    {
        EXPECT_NEAR(0.0, constants::log_approx(1.0, 25), 1e-12);
    }

    TEST(ConstantsTests, LogApprox_AtE)
    {
        // Arctanh series converges; 25 iterations accurate to 1e-8 for x=e
        EXPECT_NEAR(1.0, constants::log_approx(std::exp(1.0), 25), 1e-6);
    }

    TEST(ConstantsTests, LogApprox_NonPositiveReturnsSentinel)
    {
        EXPECT_NEAR(-1.0, constants::log_approx(0.0, 25), 1e-12);
        EXPECT_NEAR(-1.0, constants::log_approx(-5.0, 25), 1e-12);
    }

    TEST(ConstantsTests, Bohr2Ang_OneBohr)
    {
        // 1 Bohr = a₀ Å = 0.529177210903 Å
        EXPECT_NEAR(0.529177210903, constants::bohr2ang(1.0), 1e-10);
    }

    TEST(ConstantsTests, Bohr2Ang_Zero)
    {
        EXPECT_NEAR(0.0, constants::bohr2ang(0.0), 1e-12);
    }

    TEST(ConstantsTests, Ang2Bohr_OneAngstrom)
    {
        EXPECT_NEAR(1.0 / 0.529177210903, constants::ang2bohr(1.0), 1e-8);
    }

    TEST(ConstantsTests, Ang2Bohr_Roundtrip)
    {
        // ang2bohr(bohr2ang(x)) ≈ x
        double x = 2.5;
        EXPECT_NEAR(x, constants::ang2bohr(constants::bohr2ang(x)), 1e-10);
    }

    TEST(ConstantsTests, CubicBohr2Ang_OneBohr3)
    {
        // 1 Bohr³ = a₀³ Å³
        double expected = 0.529177210903 * 0.529177210903 * 0.529177210903;
        EXPECT_NEAR(expected, constants::cubic_bohr2ang(1.0), 1e-10);
    }

    TEST(ConstantsTests, CubicAng2Bohr_Roundtrip)
    {
        double x = 3.0;
        EXPECT_NEAR(x, constants::cubic_ang2bohr(constants::cubic_bohr2ang(x)), 1e-8);
    }

    TEST(ConstantsTests, Factorial_SmallValues)
    {
        EXPECT_EQ(1LL, static_cast<long long>(constants::ft_fun(0)));
        EXPECT_EQ(1LL, static_cast<long long>(constants::ft_fun(1)));
        EXPECT_EQ(120LL, static_cast<long long>(constants::ft_fun(5)));
        EXPECT_EQ(3628800LL, static_cast<long long>(constants::ft_fun(10)));
    }

    // -----------------------------------------------------------------------

    TEST(BesselTests, BesselJ0_AtZero)
    {
        // j₀(0) = 1 (limit of sin(x)/x as x→0)
        EXPECT_NEAR(1.0, bessel_first_kind(0, 0.0), 1e-12);
    }

    TEST(BesselTests, BesselJ1_AtZero)
    {
        // j_l(0) = 0 for l > 0
        EXPECT_NEAR(0.0, bessel_first_kind(1, 0.0), 1e-12);
        EXPECT_NEAR(0.0, bessel_first_kind(5, 0.0), 1e-12);
    }

    TEST(BesselTests, BesselJ0_AtOne)
    {
        // j₀(1) = sin(1)/1
        EXPECT_NEAR(std::sin(1.0), bessel_first_kind(0, 1.0), 1e-12);
    }

    TEST(BesselTests, BesselJ1_AtOne)
    {
        // j₁(1) = (sin(1) - cos(1)) / 1
        double expected = std::sin(1.0) - std::cos(1.0);
        EXPECT_NEAR(expected, bessel_first_kind(1, 1.0), 1e-12);
    }

    TEST(BesselTests, BesselJ2_AtOne)
    {
        // j₂(1) = (2·sin(1) - 3·cos(1)) / 1
        double expected = 2.0 * std::sin(1.0) - 3.0 * std::cos(1.0);
        EXPECT_NEAR(expected, bessel_first_kind(2, 1.0), 1e-12);
    }

    TEST(BesselTests, BesselJ0_AtPi)
    {
        // j₀(π) = sin(π)/π ≈ 0
        EXPECT_NEAR(std::sin(PI_VAL) / PI_VAL, bessel_first_kind(0, PI_VAL), 1e-12);
    }

    TEST(BesselTests, BesselJ_HigherOrder_PositiveAndFinite)
    {
        // l=7 exercises the continued-fraction fallback path
        double r = bessel_first_kind(7, 2.0);
        EXPECT_TRUE(std::isfinite(r));
        EXPECT_TRUE(r > 0.0);
    }

    TEST(BesselTests, BesselJ_RecurrenceCheck)
    {
        // Recurrence: j_{l-1}(x) + j_{l+1}(x) = (2l+1)/x · j_l(x)
        double x = 3.0;
        int l = 3;
        double jlm1 = bessel_first_kind(l - 1, x);
        double jl = bessel_first_kind(l, x);
        double jlp1 = bessel_first_kind(l + 1, x);
        double lhs = jlm1 + jlp1;
        double rhs = (2.0 * l + 1.0) / x * jl;
        EXPECT_NEAR(lhs, rhs, 1e-10);
    }

    // -----------------------------------------------------------------------
    // HypergeometricTests — 2F1(a,b;c;x)
    // -----------------------------------------------------------------------
    TEST(HypergeometricTests, Identity_2F1_Zero_IsOne)
    {
        // 2F1(a,b;c;0) = 1 for any a,b,c
        double r = hypergeometric(1.0, 2.0, 3.0, 0.0);
        EXPECT_NEAR(1.0, r, 1e-12);
    }

    TEST(HypergeometricTests, KnownValue_2F1_1_1_2_Half)
    {
        // 2F1(1,1;2;0.5) = -2*ln(0.5) = 2*ln(2) ≈ 1.386294...
        double expected = -2.0 * std::log(0.5);
        double r = hypergeometric(1.0, 1.0, 2.0, 0.5);
        EXPECT_NEAR(expected, r, 1e-6);
    }

    TEST(HypergeometricTests, KnownValue_2F1_Half_Half_ThreeHalves_Half)
    {
        // 2F1(0.5,0.5;1.5;0.5) — finite, positive
        double r = hypergeometric(0.5, 0.5, 1.5, 0.5);
        EXPECT_TRUE(std::isfinite(r));
        EXPECT_TRUE(r > 1.0);
    }

    TEST(HypergeometricTests, Symmetry_ab_equals_ba)
    {
        // 2F1(a,b;c;x) = 2F1(b,a;c;x)
        double r1 = hypergeometric(2.0, 3.0, 5.0, 0.3);
        double r2 = hypergeometric(3.0, 2.0, 5.0, 0.3);
        EXPECT_NEAR(r1, r2, 1e-10);
    }

    TEST(HypergeometricTests, NegativeX_ReturnsFinite)
    {
        double r = hypergeometric(1.0, 2.0, 3.0, -0.5);
        EXPECT_TRUE(std::isfinite(r));
    }

    // -----------------------------------------------------------------------
    // NormGaussTests — constants::normgauss(type, exponent)
    // -----------------------------------------------------------------------
    TEST(NormGaussTests, SShell_IsPositive)
    {
        double n = constants::normgauss(1, 1.0);
        EXPECT_TRUE(n > 0.0);
        EXPECT_TRUE(std::isfinite(n));
    }

    TEST(NormGaussTests, SShell_ScalesWithExponent)
    {
        // Normalization grows with exponent for s-type
        double n1 = constants::normgauss(1, 1.0);
        double n2 = constants::normgauss(1, 4.0);
        EXPECT_GT(n2, n1);
    }

    TEST(NormGaussTests, PShell_IsPositive)
    {
        double n = constants::normgauss(2, 1.0);
        EXPECT_TRUE(n > 0.0);
    }

    TEST(NormGaussTests, DShell_IsPositive)
    {
        double n = constants::normgauss(5, 1.0);
        EXPECT_TRUE(n > 0.0);
    }

    TEST(NormGaussTests, SShell_KnownValue)
    {
        // normgauss for s-type (0,0,0), exponent α:
        // N = (2α/π)^(9/4) * sqrt(1/1) = (2α/π)^(9/4)
        // For α=1: N = (2/π)^(9/4)
        double alpha = 1.0;
        double expected = std::pow(2.0 * alpha / PI_VAL, 9.0 / 4.0);
        double result = constants::normgauss(1, alpha);
        EXPECT_NEAR(expected, result, 1e-10);
    }

    // -----------------------------------------------------------------------
    // AssocLegendreTests — constants::associated_legendre_polynomial(l, m, x)
    // -----------------------------------------------------------------------
    TEST(AssocLegendreTests, P00_Is_One)
    {
        EXPECT_NEAR(1.0, constants::associated_legendre_polynomial(0, 0, 0.5), 1e-12);
    }

    TEST(AssocLegendreTests, P10_At_Half_Is_Half)
    {
        // P_1^0(x) = x
        EXPECT_NEAR(0.5, constants::associated_legendre_polynomial(1, 0, 0.5), 1e-12);
    }

    TEST(AssocLegendreTests, P11_AtOne_Is_Zero)
    {
        // P_1^1(x) = sqrt(1-x²), at x=1: 0
        EXPECT_NEAR(0.0, constants::associated_legendre_polynomial(1, 1, 1.0), 1e-12);
    }

    TEST(AssocLegendreTests, P20_Is_Legendre_Polynomial)
    {
        // P_2^0(x) = (3x²-1)/2
        double x = 0.6;
        double expected = 0.5 * (3 * x * x - 1);
        EXPECT_NEAR(expected, constants::associated_legendre_polynomial(2, 0, x), 1e-12);
    }

    TEST(AssocLegendreTests, P21_At_Zero)
    {
        // P_2^1(x) = 3x*sqrt(1-x²), at x=0: 0
        EXPECT_NEAR(0.0, constants::associated_legendre_polynomial(2, 1, 0.0), 1e-12);
    }

    TEST(AssocLegendreTests, P22_At_Zero)
    {
        // P_2^2(x) = -3(x²-1) = 3(1-x²), at x=0: 3
        EXPECT_NEAR(-3.0 * (0.0 - 1.0), constants::associated_legendre_polynomial(2, 2, 0.0), 1e-12);
    }

    TEST(AssocLegendreTests, NegativeM_P1m1)
    {
        // P_1^{-1}(x) = -0.5*sqrt(1-x²)
        double x = 0.5;
        double expected = -0.5 * std::sqrt(1 - x * x);
        EXPECT_NEAR(expected, constants::associated_legendre_polynomial(1, -1, x), 1e-12);
    }

    // -----------------------------------------------------------------------
    // CartesianToSphericalTests — constants::cartesian_to_spherical(x,y,z)
    // -----------------------------------------------------------------------
    TEST(CartesianToSphericalTests, Origin_Has_Zero_Radius)
    {
        vec r = constants::cartesian_to_spherical(0.0, 0.0, 0.0);
        EXPECT_NEAR(0.0, r[0], 1e-12); // r
    }

    TEST(CartesianToSphericalTests, UnitX_Gives_Correct_Angles)
    {
        vec r = constants::cartesian_to_spherical(1.0, 0.0, 0.0);
        EXPECT_NEAR(1.0, r[0], 1e-12); // r=1
        EXPECT_NEAR(PI_VAL / 2.0, r[1], 1e-12); // theta=pi/2
        EXPECT_NEAR(0.0, r[2], 1e-12); // phi=0
    }

    TEST(CartesianToSphericalTests, UnitZ_Has_Zero_Theta)
    {
        vec r = constants::cartesian_to_spherical(0.0, 0.0, 1.0);
        EXPECT_NEAR(1.0, r[0], 1e-12); // r=1
        EXPECT_NEAR(0.0, r[1], 1e-12); // theta=0
    }

    TEST(CartesianToSphericalTests, UnitY_Gives_PhiHalfPi)
    {
        vec r = constants::cartesian_to_spherical(0.0, 1.0, 0.0);
        EXPECT_NEAR(1.0, r[0], 1e-12); // r=1
        EXPECT_NEAR(PI_VAL / 2.0, r[2], 1e-12); // phi=pi/2
    }

    TEST(CartesianToSphericalTests, Radius_Is_Euclidean_Norm)
    {
        vec r = constants::cartesian_to_spherical(3.0, 4.0, 0.0);
        EXPECT_NEAR(5.0, r[0], 1e-12); // r=5
    }

    TEST(CartesianToSphericalTests, Inverse_Recover_Cartesian)
    {
        // Convert (1,1,1) → spherical → back to Cartesian
        vec r_sph = constants::cartesian_to_spherical(1.0, 1.0, 1.0);
        double r = r_sph[0], theta = r_sph[1], phi = r_sph[2];
        double x = r * std::sin(theta) * std::cos(phi);
        double y = r * std::sin(theta) * std::sin(phi);
        double z = r * std::cos(theta);
        EXPECT_NEAR(1.0, x, 1e-12);
        EXPECT_NEAR(1.0, y, 1e-12);
        EXPECT_NEAR(1.0, z, 1e-12);
    }

    // -----------------------------------------------------------------------
    // GetLambda1Tests — median eigenvalue of a 3x3 symmetric matrix
    // -----------------------------------------------------------------------
    TEST(GetLambda1Tests, DiagonalMatrix_ReturnsMidEigenvalue)
    {
        // Eigenvalues are 1, 3, 5 → median = 3
        const double a[9] = {
            1, 0, 0,
            0, 3, 0,
            0, 0, 5
        };
        double tmp[9];
        for (int i = 0; i < 9; ++i) tmp[i] = a[i];
        double lam = get_lambda_1(tmp);
        EXPECT_NEAR(3.0, lam, 1e-10);
    }

    TEST(GetLambda1Tests, DiagonalMatrix_AllEqual_ReturnsThat)
    {
        const double a[9] = {
            2, 0, 0,
            0, 2, 0,
            0, 0, 2
        };
        double tmp[9];
        for (int i = 0; i < 9; ++i) tmp[i] = a[i];
        double lam = get_lambda_1(tmp);
        // all eigenvalues equal 2; any is "middle"
        EXPECT_NEAR(2.0, lam, 1e-10);
    }

    TEST(GetLambda1Tests, DiagonalDescending_ReturnsMid)
    {
        // Eigenvalues 10, 5, 1 → median = 5
        const double a[9] = {
            10, 0, 0,
            0,  5, 0,
            0,  0, 1
        };
        double tmp[9];
        for (int i = 0; i < 9; ++i) tmp[i] = a[i];
        double lam = get_lambda_1(tmp);
        EXPECT_NEAR(5.0, lam, 1e-10);
    }

    TEST(GetLambda1Tests, SymmetricMatrix_ReturnsFiniteValue)
    {
        // 3x3 symmetric: known eigenvalues can be verified with characteristic polynomial
        // A = [[4,2,0],[2,3,1],[0,1,2]] → trace=9, det=14
        const double a[9] = {
            4, 2, 0,
            2, 3, 1,
            0, 1, 2
        };
        double tmp[9];
        for (int i = 0; i < 9; ++i) tmp[i] = a[i];
        double lam = get_lambda_1(tmp);
        EXPECT_TRUE(std::isfinite(lam));
        // Eigenvalues must be between min(diag)=2 and max(diag)=4 for a diag-dominant case
        EXPECT_TRUE(lam >= 1.0 && lam <= 5.0);
    }

    TEST(GetLambda1Tests, SymmetricMatrix_MedianIsInBounds)
    {
        // For any symmetric matrix, median eigenvalue is bounded by extremes
        const double a[9] = {
            3, 1, 0,
            1, 4, 2,
            0, 2, 5
        };
        double tmp[9];
        for (int i = 0; i < 9; ++i) tmp[i] = a[i];
        double lam = get_lambda_1(tmp);
        EXPECT_TRUE(std::isfinite(lam));
        EXPECT_TRUE(lam > 0.0); // positive definite matrix
    }

    TEST(Sph2CartTests, ReturnsEmbeddedMatrices)
    {
        EXPECT_DOUBLE_EQ(constants::sph2cart(0)[0], 1.0);
        EXPECT_DOUBLE_EQ(constants::sph2cart(1)[2 * 3 + 0], 1.0);
        EXPECT_DOUBLE_EQ(constants::sph2cart(2)[2 * 5 + 0], 0.57735026918962576);
        EXPECT_DOUBLE_EQ(constants::sph2cart(3)[9 * 7 + 4], 1.0);
        EXPECT_DOUBLE_EQ(constants::sph2cart(4)[0 * 9 + 0], 0.1690308509457033);
        EXPECT_DOUBLE_EQ(constants::sph2cart(4)[14 * 9 + 0], 0.06338656910463875);
        EXPECT_DOUBLE_EQ(constants::sph2cart(5)[20 * 11 + 1], 0.48412291827592713);
        EXPECT_DOUBLE_EQ(constants::sph2cart(6)[27 * 13 + 11], 0.67169328938139627);
        EXPECT_DOUBLE_EQ(constants::sph2cart(10)[65 * 21 + 19], 0.59362791713657326);
    }

    //Lukas Seifert's derivation: libcint's c2s with ORCA's phase for |m| = 3, 4, 7, 8 over ORCA's angular norm, which drops sqrt((2l-1)(2l-3)) from l = 5 on
    TEST(Sph2CartTests, MatchesLibcintWithOrcaPhase)
    {
        for (int l = 2; l <= 10; l++)
        {
            const int nc = constants::n_cart(l), nsph = constants::n_spher(l);
            vec identity(nc * nc, 0.0), sph(nsph * nc, 0.0);
            for (int i = 0; i < nc; i++) identity[i + nc * i] = 1.0;
            libcint::CINTc2s_bra_sph(sph.data(), nc, identity.data(), l);
            const double norm = l == 2 ? 0.5 * std::sqrt(15.0 / constants::PI) : l == 3 ? 0.5 * std::sqrt(105.0 / constants::PI) : l == 4 ? 1.5 * std::sqrt(35.0 / constants::PI) : 0.5 * std::sqrt((2.0 * l + 1.0) / constants::PI);
            for (int cart = 0; cart < nc; cart++)
            {
                int e[3];
                constants::type2vector(constants::first_type[l] + cart, e);
                int lc = 0;
                for (int lx = l; lx > e[0]; lx--) lc += l - lx + 1;
                lc += l - e[0] - e[1];
                for (int m = -l; m <= l; m++)
                {
                    const double phase = std::abs(m) % 4 == 3 || std::abs(m) % 4 == 0 && m != 0 ? -1.0 : 1.0;
                    EXPECT_NEAR(constants::sph2cart(l)[cart * nsph + (m == 0 ? 0 : m > 0 ? 2 * m - 1 : -2 * m)], phase / norm * sph[(m + l) + nsph * lc], 1e-12) << "l " << l << " cart " << cart << " m " << m;
                }
            }
        }
    }

    //a spherical-coefficient wavefunction built through push_back_spherical_shell against the same MOs summed
    //directly as sum_m c_m R_lm(d) exp(-a r^2) with R_lm = phase/K_l * r^l * spherical_harmonic(d/r) (the ORCA convention
    //of MatchesLibcintWithOrcaPhase); gradient and Laplacian against central differences of that reference
    TEST(Sph2CartTests, SphericalCoefficientsRoundTripThroughTheCartesianWavefunction)
    {
        const std::array<std::array<double, 3>, 2> centres{ { { 0.1, -0.2, 0.3 }, { 1.7, 0.4, -0.9 } } };
        struct sh { int atom, l; vec exps; };
        const std::vector<sh> shells{ { 0, 0, { 3.1, 0.7 } }, { 0, 1, { 1.3 } }, { 0, 2, { 1.9, 0.5 } }, { 0, 3, { 0.9 } }, { 0, 4, { 0.8 } }, { 0, 5, { 0.7 } },
                                      { 1, 0, { 2.2 } }, { 1, 1, { 1.1, 0.4 } }, { 1, 2, { 0.6 } } };
        const vec occ{ 2.0, 1.0 };
        //coef[mo][shell][m][s]: deterministic pseudo-random, m in ORCA order 0,+1,-1,+2,-2,...
        std::vector<vec3> coef(occ.size(), vec3(shells.size()));
        int seed = 0;
        for (int mo = 0; mo < (int)occ.size(); mo++)
            for (int i = 0; i < (int)shells.size(); i++)
            {
                coef[mo][i].assign(constants::n_spher(shells[i].l), vec(shells[i].exps.size()));
                for (auto& row : coef[mo][i]) for (double& c : row) c = std::sin(1.3 * ++seed);
            }

        WFN wavy(e_origin::NOT_YET_DEFINED);
        wavy.push_back_atom("O", centres[0][0], centres[0][1], centres[0][2], 8);
        wavy.push_back_atom("C", centres[1][0], centres[1][1], centres[1][2], 6);
        for (int mo = 0; mo < (int)occ.size(); mo++) wavy.push_back_MO(mo + 1, occ[mo], -1.0 - mo);
        std::vector<primitive> prims;
        ivec start;
        for (const sh& s : shells)
        {
            start.push_back((int)prims.size());
            for (const double e : s.exps) prims.emplace_back(s.atom + 1, s.l + 1, e, 1.0);
        }
        for (int mo = 0; mo < (int)occ.size(); mo++)
            for (int i = 0; i < (int)shells.size(); i++)
                wavy.push_back_spherical_shell(mo, shells[i].l, coef[mo][i], prims, start[i], (int)shells[i].exps.size());
        wavy.set_exp_cutoff();
        ASSERT_EQ(wavy.get_nex(), 76);

        auto R_lm = [](const int l, const int m, const double* d) {
            const double K = l == 2 ? 0.5 * std::sqrt(15.0 / constants::PI) : l == 3 ? 0.5 * std::sqrt(105.0 / constants::PI) : l == 4 ? 1.5 * std::sqrt(35.0 / constants::PI) : 0.5 * std::sqrt((2.0 * l + 1.0) / constants::PI);
            const double phase = std::abs(m) % 4 == 3 || std::abs(m) % 4 == 0 && m != 0 ? -1.0 : 1.0;
            const double r = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
            if (r == 0.0) return l == 0 ? phase / K * constants::c_1_4p : 0.0;
            const double u[3] = { d[0] / r, d[1] / r, d[2] / r };
            return phase / K * std::pow(r, l) * constants::spherical_harmonic(l, m, u);
        };
        auto mo_ref = [&](const d3& P, const int mo) {
            double v = 0;
            for (int i = 0; i < (int)shells.size(); i++)
            {
                const double d[3] = { P[0] - centres[shells[i].atom][0], P[1] - centres[shells[i].atom][1], P[2] - centres[shells[i].atom][2] };
                const double r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                for (int s = 0; s < (int)shells[i].exps.size(); s++)
                    for (int m = -shells[i].l; m <= shells[i].l; m++)
                        v += coef[mo][i][m == 0 ? 0 : m > 0 ? 2 * m - 1 : -2 * m][s] * R_lm(shells[i].l, m, d) * std::exp(-shells[i].exps[s] * r2);
            }
            return v;
        };
        auto dens_ref = [&](const d3& P) { double v = 0; for (int mo = 0; mo < (int)occ.size(); mo++) v += occ[mo] * std::pow(mo_ref(P, mo), 2); return v; };

        const std::vector<d3> points{ { 0.0, 0.0, 0.0 }, { 0.5, -0.3, 0.8 }, { 1.7, 0.4, -0.9 }, { 1.0, 0.1, -0.3 }, { -0.7, 1.1, 0.2 }, { 2.4, -0.8, -1.5 } };
        const double h = 1e-4; // central differences: truncation h^2 rho''', rounding eps rho / h^2, both below 1e-6 relative
        for (const d3& P : points)
        {
            for (int mo = 0; mo < (int)occ.size(); mo++) EXPECT_NEAR(wavy.computeMO(P, mo), mo_ref(P, mo), 1e-10) << "MO " << mo << " at " << P[0] << " " << P[1] << " " << P[2];
            const double rho = dens_ref(P);
            EXPECT_NEAR(wavy.compute_dens(P), rho, 1e-10) << "density at " << P[0] << " " << P[1] << " " << P[2];
            d3 grad;
            wavy.computeGrad(P, grad);
            double lap = 0;
            for (int k = 0; k < 3; k++)
            {
                d3 Pp = P, Pm = P;
                Pp[k] += h; Pm[k] -= h;
                const double rp = dens_ref(Pp), rm = dens_ref(Pm);
                EXPECT_NEAR(grad[k], (rp - rm) / (2 * h), 1e-5 * (1 + std::abs(grad[k]))) << "gradient " << k << " at " << P[0] << " " << P[1] << " " << P[2];
                lap += (rp - 2 * rho + rm) / (h * h);
            }
            EXPECT_NEAR(wavy.computeLap(P), lap, 1e-5 * (1 + std::abs(lap))) << "Laplacian at " << P[0] << " " << P[1] << " " << P[2];
        }
    }

    // ------------------------------------------------------------------
    // Coverage of the small numerics in nos_math / mat_nos_math / vec_nos_math
    // (7 %, 32 % and 45 % of lines before these, 19 Sep 2026)

    TEST(MathTests, ReshapeFlattenTransposeRoundTrip)
    {
        const vec flat = { 1, 2, 3, 4, 5, 6 };
        const dMatrix2 m = reshape<dMatrix2>(flat, Shape2D(2, 3));
        EXPECT_EQ(m(1, 0), 4.0);
        const dMatrix2 t = transpose(m);
        EXPECT_EQ(t.extent(0), 3u);
        EXPECT_EQ(t(0, 1), 4.0);
        const dMatrix1 back = flatten<dMatrix1>(m);
        for (int i = 0; i < 6; i++) EXPECT_EQ(back(i), flat[i]);
        const dMatrix3 m3 = reshape<dMatrix3>(back, Shape3D(1, 2, 3));
        EXPECT_EQ(m3(0, 1, 2), 6.0);
        const dMatrix4 m4 = reshape<dMatrix4>(m3, Shape4D(1, 1, 2, 3));
        EXPECT_EQ(flatten<dMatrix1>(m4)(5), 6.0);
        EXPECT_EQ(reshape<dMatrix2>(m4, Shape2D(3, 2))(2, 1), 6.0);
        const cMatrix2 cm = reshape<cMatrix2>(cvec{ {1, 1}, {2, 2} }, Shape2D(1, 2));
        EXPECT_EQ(cm(0, 1).imag(), 2.0);
        EXPECT_EQ(transpose(cm)(1, 0).real(), 2.0);
        const iMatrix2 im = reshape<iMatrix2>(ivec{ 1, 2, 3, 4 }, Shape2D(2, 2));
        EXPECT_EQ(transpose(im)(0, 1), 3);
        EXPECT_EQ(flatten<iMatrix1>(im)(3), 4);
        // std::vector flavours
        const vec2 v2 = { { 1, 2, 3 }, { 4, 5, 6 } };
        EXPECT_EQ(flatten<double>(v2), flat);
        EXPECT_EQ(transpose(v2)[2][1], 6.0);
        EXPECT_EQ(transpose(flat, 2, 3)[1], 4.0);
        EXPECT_EQ(transpose(flat)[0].size(), 6u);
        const vec3 v3 = { { { 1, 2 }, { 3, 4 } }, { { 5, 6 }, { 7, 8 } } };
        EXPECT_EQ(flatten<double>(v3), (vec{ 1, 2, 3, 4, 5, 6, 7, 8 }));
        // both put the old third index first: [k][i][j] = v3[i][j][k]
        EXPECT_EQ(transpose(v3)[1][0][1], 4.0);
        EXPECT_EQ(reorder3D(v3)[1][0][1], 4.0);
        EXPECT_TRUE(transpose(vec3{}).empty());
        EXPECT_TRUE(reorder3D(vec3{}).empty());
        const ivec2 iv2 = { { 1, 2 }, { 3, 4 } };
        EXPECT_EQ(transpose(iv2)[0][1], 3);
        EXPECT_EQ(flatten<int>(iv2)[2], 3);
        EXPECT_EQ(transpose(ivec{ 1, 2, 3, 4 }, 2, 2)[1], 3);
        const cvec2 cv2 = { { { 1, 0 }, { 0, 1 } } };
        EXPECT_EQ(transpose(cv2)[1][0].imag(), 1.0);
        EXPECT_EQ(flatten<cdouble>(cv2)[1].imag(), 1.0);
        EXPECT_EQ(transpose(cvec{ { 1, 2 } })[0][0].imag(), 2.0);
        EXPECT_EQ(transpose(cvec{ { 1, 2 }, { 3, 4 } }, 1, 2)[1].real(), 3.0);
        EXPECT_EQ(transpose(ivec{ 7 })[0][0], 7);
    }

    TEST(MathTests, DotProductsAgreeWithSelfDot)
    {
        dMatrix2 A(2, 3), B(3, 2);
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 3; j++) {
                A(i, j) = i + j + 1.0;
                B(j, i) = i * j + 2.0;
            }
        const dMatrix2 C = dot(A, B);
        ASSERT_EQ(C.extent(0), 2u);
        ASSERT_EQ(C.extent(1), 2u);
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++) {
                double s = 0;
                for (int k = 0; k < 3; k++) s += A(i, k) * B(k, j);
                EXPECT_DOUBLE_EQ(C(i, j), s);
            }
        const dMatrix2 AtA = dot(A, A, true, false);
        ASSERT_EQ(AtA.extent(0), 3u);
        EXPECT_DOUBLE_EQ(AtA(0, 0), A(0, 0) * A(0, 0) + A(1, 0) * A(1, 0));
        const dMatrix2 ABt = dot(A, transpose(B), false, true);
        EXPECT_DOUBLE_EQ(ABt(1, 1), C(1, 1));
        EXPECT_EQ(dot(dMatrix2{}, B).size(), 0u);
        // the vec2 route: self_dot against the BLAS product of the flattened matrices
        const vec2 a = { { 1, 2, 3 }, { 4, 5, 6 } }, b = { { 1, 0 }, { 0, 1 }, { 2, 2 } };
        const vec2 ab = self_dot(a, b);
        EXPECT_EQ(ab[1][0], 16.0);
        EXPECT_EQ(ab[0][1], 8.0);
        const vec2 atbt = self_dot(transpose(a), transpose(b), true, true);
        EXPECT_EQ(atbt[1][0], ab[1][0]);
        EXPECT_TRUE(self_dot(vec2{}, b).empty());
        const dMatrix2 ab_blas = dot_BLAS(flatten<double>(a), flatten<double>(b), 2, 3, 3, 2);
        EXPECT_DOUBLE_EQ(ab_blas(1, 0), ab[1][0]);
        // matrix-vector
        const vec x = { 1, 1, 1 };
        EXPECT_EQ(dot(a, x)[1], 15.0);
        EXPECT_EQ(self_dot(a, x)[1], 15.0);
        EXPECT_EQ(self_dot(a, vec{ 1, 1 }, true)[2], 9.0);
        EXPECT_EQ(dot_BLAS(flatten<double>(a), x, 2, 3)[0], 6.0);
        EXPECT_EQ(dot_BLAS(flatten<double>(a), vec{ 1, 1 }, 2, 3, true)[2], 9.0);
        dMatrix1 xm(3);
        xm(0) = xm(1) = xm(2) = 1.0;
        EXPECT_DOUBLE_EQ(dot(A, xm)(0), 6.0);
        dMatrix1 ym(2);
        ym(0) = ym(1) = 1.0;
        EXPECT_DOUBLE_EQ(dot(A, ym, true)(2), A(0, 2) + A(1, 2));
        // vector-vector, real and complex, with and without conjugation
        EXPECT_DOUBLE_EQ(dot(vec{ 1, 2, 3 }, vec{ 4, 5, 6 }), 32.0);
        EXPECT_DOUBLE_EQ(self_dot(vec{ 1, 2, 3 }, vec{ 4, 5, 6 }), 32.0);
        const cvec u = { { 1, 2 }, { 3, 4 } }, v = { { 5, 6 }, { 7, 8 } };
        EXPECT_EQ(dot(u, v, false), self_dot(u, v, false));
        EXPECT_EQ(dot(u, v, true), self_dot(u, v, true));
        EXPECT_EQ(self_dot(u, v, true), std::conj(u[0]) * v[0] + std::conj(u[1]) * v[1]);
        // complex GEMM and mat-vec
        cMatrix2 U(2, 2), I(2, 2);
        U(0, 0) = { 1, 1 }; U(0, 1) = { 2, 0 }; U(1, 0) = { 0, 3 }; U(1, 1) = { 4, 4 };
        I(0, 0) = 1.0; I(0, 1) = 0.0; I(1, 0) = 0.0; I(1, 1) = 1.0;
        const cMatrix2 UI = dot(U, I);
        EXPECT_EQ(UI(1, 1), U(1, 1));
        EXPECT_EQ(dot(U, I, true, true)(0, 1), U(1, 0));
        cMatrix1 one(2);
        one(0) = one(1) = 1.0;
        EXPECT_EQ(dot(U, one)(0), U(0, 0) + U(0, 1));
        EXPECT_EQ(dot(U, one, true)(1), U(0, 1) + U(1, 1));
        const cvec2 cu = { { { 1, 1 }, { 2, 0 } }, { { 0, 3 }, { 4, 4 } } };
        EXPECT_EQ(self_dot(cu, cu)[0][0], cu[0][0] * cu[0][0] + cu[0][1] * cu[1][0]);
        EXPECT_EQ(dot(cu, cvec{ 1.0, 1.0 })[1], cu[1][0] + cu[1][1]);
        EXPECT_EQ(self_dot(cu, cvec{ 1.0, 1.0 })[1], cu[1][0] + cu[1][1]);
        EXPECT_EQ(dot_BLAS(flatten<cdouble>(cu), cvec{ 1.0, 1.0 }, 2, 2)[1], cu[1][0] + cu[1][1]);
        EXPECT_EQ(dot_BLAS(flatten<cdouble>(cu), flatten<cdouble>(cu), 2, 2, 2, 2)(0, 0), self_dot(cu, cu)[0][0]);
        // diagonal scaling
        const dMatrix2 D = diag_dot(A, vec{ 1, 2, 3 });
        EXPECT_DOUBLE_EQ(D(1, 2), A(1, 2) * 3);
        const dMatrix2 Dt = diag_dot(A, vec{ 1, 2 }, true);
        EXPECT_DOUBLE_EQ(Dt(2, 1), A(1, 2) * 2);
        EXPECT_EQ(diag_dot(U, cvec{ 1.0, { 0, 1 } })(1, 1), U(1, 1) * cdouble(0, 1));
        // trace products: tr(S T) and tr(S^T T)
        dMatrix2 S(2, 2), T(2, 2);
        S(0, 0) = 1; S(0, 1) = 2; S(1, 0) = 3; S(1, 1) = 4;
        T(0, 0) = 5; T(0, 1) = 6; T(1, 0) = 7; T(1, 1) = 8;
        EXPECT_DOUBLE_EQ(trace_product<double>(S, T), 1 * 5 + 2 * 7 + 3 * 6 + 4 * 8);
        EXPECT_DOUBLE_EQ(trace_product<double>(S, T, true), 1 * 5 + 2 * 6 + 3 * 7 + 4 * 8);
        iMatrix2 Si(2, 2), Ti(2, 2);
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++) {
                Si(i, j) = (int)S(i, j);
                Ti(i, j) = (int)T(i, j);
            }
        EXPECT_EQ(trace_product<int>(Si, Ti), 69);
        EXPECT_EQ(trace_product<int>(Si, Ti, true), 70);
        EXPECT_EQ(trace_product<cdouble>(U, I), U(0, 0) + U(1, 1));
        EXPECT_EQ(trace_product<cdouble>(U, I, true), U(0, 0) + U(1, 1));
    }

    TEST(MathTests, SubmatricesRectanglesAndSymmetricSwaps)
    {
        dMatrix2 M(3, 3);
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) M(i, j) = 10 * i + j;
        const dMatrix2 R = get_rectangle(M, ivec{ 2, 0 });
        ASSERT_EQ(R.extent(0), 2u);
        EXPECT_EQ(R(0, 1), 21.0);
        EXPECT_EQ(R(1, 2), 2.0);
        vec sub(4);
        get_submatrix(M, sub, ivec{ 0, 2 });
        EXPECT_EQ(sub, (vec{ 0, 2, 20, 22 }));
        vec sub2(2);
        get_submatrix(M, sub2, ivec{ 1 }, ivec{ 0, 2 });
        EXPECT_EQ(sub2, (vec{ 10, 12 }));
        vec Dsub(1), Ssub(1);
        get_submatrices(M, M, Dsub, Ssub, ivec{ 1 });
        EXPECT_EQ(Dsub[0], 11.0);
        EXPECT_EQ(Ssub[0], 11.0);
        iMatrix2 Mi(2, 2);
        Mi(0, 0) = 1; Mi(0, 1) = 2; Mi(1, 0) = 3; Mi(1, 1) = 4;
        ivec isub(1), isub2(2);
        get_submatrix(Mi, isub, ivec{ 1 });
        get_submatrix(Mi, isub2, ivec{ 0 }, ivec{ 0, 1 });
        EXPECT_EQ(isub[0], 4);
        EXPECT_EQ(isub2[1], 2);
        EXPECT_EQ(get_rectangle(Mi, ivec{ 1 })(0, 0), 3);
        cMatrix2 Mc(2, 2);
        Mc(0, 0) = { 1, 1 }; Mc(0, 1) = { 2, 2 }; Mc(1, 0) = { 3, 3 }; Mc(1, 1) = { 4, 4 };
        cvec csub(1), csub2(2);
        get_submatrix(Mc, csub, ivec{ 1 });
        get_submatrix(Mc, csub2, ivec{ 1 }, ivec{ 0, 1 });
        EXPECT_EQ(csub[0].imag(), 4.0);
        EXPECT_EQ(csub2[0].real(), 3.0);
        EXPECT_EQ(get_rectangle(Mc, ivec{ 1, 0 })(1, 1).real(), 2.0);
        // a symmetric swap keeps the matrix symmetric and moves the diagonal
        dMatrix2 Sy(3, 3);
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) Sy(i, j) = i * i + j * j + (i == j ? 100.0 * i : 0.0);
        swap_rows_cols_symm(Sy, 0, 2);
        EXPECT_EQ(Sy(0, 0), 8 + 200.0);
        EXPECT_EQ(Sy(2, 2), 0.0);
        EXPECT_EQ(Sy(0, 1), Sy(1, 0));
        EXPECT_EQ(Sy(0, 1), 5.0);
        cMatrix2 Sc(2, 2);
        Sc(0, 0) = 1.0; Sc(0, 1) = { 0, 2 }; Sc(1, 0) = { 0, 2 }; Sc(1, 1) = 3.0;
        swap_rows_cols_symm(Sc, 0, 1);
        EXPECT_EQ(Sc(0, 0).real(), 3.0);
        EXPECT_EQ(Sc(1, 0).imag(), 2.0);
        // element-wise powers, both containers
        const vec2 sq = elementWiseExponentiation(vec2{ { 2, 3 }, { 4, 5 } }, 2.0);
        EXPECT_EQ(sq[1][1], 25.0);
        dMatrix2 base(2, 2);
        base(0, 0) = 4; base(0, 1) = 9; base(1, 0) = 16; base(1, 1) = 25;
        EXPECT_DOUBLE_EQ(elementWiseExponentiation(base, 0.5)(1, 0), 4.0);
    }

    TEST(MathTests, LapackWrappersOnSmallSystems)
    {
        // eigenpairs of [[2,1],[1,2]]: 1 and 3, eigenvectors as columns of A
        vec A = { 2, 1, 1, 2 }, W(2);
        make_Eigenvalues(A, W);
        EXPECT_NEAR(W[0], 1.0, 1e-12);
        EXPECT_NEAR(W[1], 3.0, 1e-12);
        EXPECT_NEAR(std::abs(A[0 * 2 + 1]), std::sqrt(0.5), 1e-12);
        EXPECT_TRUE(isSymmetricViaEigenvalues(vec{ 2, 1, 1, 2 }, 2));
        // a rotation has complex eigenvalues
        EXPECT_FALSE(isSymmetricViaEigenvalues(vec{ 0, -1, 1, 0 }, 2));
        // (sqrt M)^2 = M
        vec M = { 2, 1, 1, 2 }, Wm(2);
        const vec S = mat_sqrt(M, Wm);
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++) {
                double s = 0;
                for (int k = 0; k < 2; k++) s += S[i * 2 + k] * S[k * 2 + j];
                EXPECT_NEAR(s, (i == j) ? 2.0 : 1.0, 1e-12);
            }
        // pseudo-inverse of a 2x3 matrix of full row rank: A A+ = I
        dMatrix2 P(2, 3);
        P(0, 0) = 1; P(0, 1) = 0; P(0, 2) = 1;
        P(1, 0) = 0; P(1, 1) = 1; P(1, 2) = 1;
        const dMatrix2 Pinv = LAPACKE_invert(P);
        ASSERT_EQ(Pinv.extent(0), 3u);
        ASSERT_EQ(Pinv.extent(1), 2u);
        const dMatrix2 PP = dot(P, Pinv);
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++) EXPECT_NEAR(PP(i, j), i == j ? 1.0 : 0.0, 1e-12);
        // square solve, three entry points, one system: [[2,1],[1,3]] x = [3,4] -> x = [1,1]
        vec b = { 3, 4 };
        solve_linear_system(vec2{ { 2, 1 }, { 1, 3 } }, b);
        EXPECT_NEAR(b[0], 1.0, 1e-12);
        EXPECT_NEAR(b[1], 1.0, 1e-12);
        vec flat = { 2, 1, 1, 3 };
        b = { 3, 4 };
        solve_linear_system(flat, 2, b);
        EXPECT_NEAR(b[1], 1.0, 1e-12);
        // least squares, consistent 3x2 system -> exact solution [1,2]
        vec Als = { 1, 0, 0, 1, 1, 1 }, bls = { 1, 2, 3 };
        solve_linear_system(Als, 3, 2, bls);
        EXPECT_NEAR(bls[0], 1.0, 1e-12);
        EXPECT_NEAR(bls[1], 2.0, 1e-12);
    }

    TEST(MathTests, NnlsClampsAtZero)
    {
        // consistent system with a non-negative solution: recovered exactly
        dMatrix2 A(3, 2);
        A(0, 0) = 1; A(0, 1) = 0; A(1, 0) = 0; A(1, 1) = 1; A(2, 0) = 1; A(2, 1) = 1;
        dMatrix1 b(3);
        b(0) = 2; b(1) = 3; b(2) = 5;
        NNLSResult r = nnls(A, b);
        EXPECT_EQ(r.status, 0);
        ASSERT_EQ(r.x.size(), 2u);
        EXPECT_NEAR(r.x[0], 2.0, 1e-12);
        EXPECT_NEAR(r.x[1], 3.0, 1e-12);
        EXPECT_NEAR(r.rnorm, 0.0, 1e-12);
        // unconstrained optimum has a negative entry: it is clamped and the residual is what remains
        dMatrix2 I(2, 2);
        I(0, 0) = 1; I(0, 1) = 0; I(1, 0) = 0; I(1, 1) = 1;
        dMatrix1 c(2);
        c(0) = 1; c(1) = -1;
        r = nnls(I, c);
        EXPECT_EQ(r.status, 0);
        EXPECT_NEAR(r.x[0], 1.0, 1e-12);
        EXPECT_NEAR(r.x[1], 0.0, 1e-12);
        EXPECT_NEAR(r.rnorm, 1.0, 1e-12);
    }

    // ------------------------------------------------------------------
    // Real spherical harmonics: every l branch of constants::spherical_harmonic
    // (l = 0..9 explicit, l >= 10 through real_spherical) is orthonormal on a
    // Lebedev sphere, and the collapsed (l, d, coefs) form is the m-sum

    TEST(SphericalHarmonicTests, OrthonormalOnLebedevSphereThroughL11)
    {
        const int order = 1202; // exact to degree 59, products up to l = 11 need 22
        vec x(order), y(order), z(order), w(order);
        lebedev_sphere().ld_by_order(order, x.data(), y.data(), z.data(), w.data());
        for (int l = 0; l <= 11; l++)
            for (int m = -l; m <= l; m++)
                for (int m2 = -l; m2 <= l; m2++) {
                    double s = 0.0;
                    for (int p = 0; p < order; p++) {
                        const double d[3] = { x[p], y[p], z[p] };
                        s += w[p] * constants::spherical_harmonic(l, m, d) * constants::spherical_harmonic(l, m2, d);
                    }
                    EXPECT_NEAR(constants::FOUR_PI * s, m == m2 ? 1.0 : 0.0, 1e-9) << "l " << l << " m " << m << " m' " << m2;
                }
        // different l are orthogonal too (one pair per l suffices to catch a wrong prefactor)
        for (int l = 1; l <= 11; l++) {
            double s = 0.0;
            for (int p = 0; p < order; p++) {
                const double d[3] = { x[p], y[p], z[p] };
                s += w[p] * constants::spherical_harmonic(l, 0, d) * constants::spherical_harmonic(l - 1, 0, d);
            }
            EXPECT_NEAR(s, 0.0, 1e-10) << l;
        }
    }

    // scipy.special.lpmv reference (Condon-Shortley phase removed) at theta 0.7, phi 1.9: the same
    // sign convention through the explicit formulas (l 2, 9) and the l >= 10 fallback, which used
    // to hand a negative m to std::assoc_legendre and gave 0 with libstdc++
    TEST(SphericalHarmonicTests, MatchesScipyReferenceValuesAcrossTheFallback)
    {
        const double d[3] = { -0.208268857072881, 0.6096232539228104, 0.7648421872844885 };
        const std::vector<std::tuple<int, int, double>> ref = {
            { 2, 1, -0.1740351075891621 }, { 9, -2, 0.34609970521259786 }, { 9, 3, -0.06954492290853755 },
            { 10, -3, 0.22976462110859847 }, { 10, 4, 0.0550889727830941 },
            { 11, -11, 0.005521174076469921 }, { 11, 0, 0.22438475059582277 } };
        for (const auto& [l, m, y] : ref)
            EXPECT_NEAR(constants::spherical_harmonic(l, m, d), y, 1e-12) << "l " << l << " m " << m;
    }

    TEST(SphericalHarmonicTests, CollapsedFormIsTheCoefficientSum)
    {
        const double d[3] = { 0.3, -0.5, 0.812 };
        for (int l = 0; l <= 10; l++) {
            vec coefs(2 * l + 1);
            for (int m = -l; m <= l; m++) coefs[m + l] = 0.1 * (m + 1) + 0.01 * l;
            double expect = 0.0;
            for (int m = -l; m <= l; m++) expect += coefs[m + l] * constants::spherical_harmonic(l, m, d);
            EXPECT_NEAR(constants::spherical_harmonic(l, d, coefs.data()), expect, 1e-12) << l;
        }
    }

} // namespace NoSpherA2UnitTests
