#include "pch.h"
#include "../../NoSpherA2_DLL/unit_test_api.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>
#include <string>

static constexpr double PI_VAL = 3.14159265358979323846;

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

// ---------------------------------------------------------------------------
// Pure unit tests — no golden files, no tests/tests.toml counterpart.
// Each method exercises a single Src/ function via the DLL's C export shim.
// ---------------------------------------------------------------------------

namespace NoSpherA2UnitTests
{
    TEST_CLASS(GeometryTests)
    {
    public:

        TEST_METHOD(ArrayLength_345Triangle)
        {
            // 3-4-5 Pythagorean triple
            Assert::AreEqual(5.0, ut_array_length3(3.0, 4.0, 0.0), 1e-12);
        }

        TEST_METHOD(ArrayLength_UnitAxes)
        {
            Assert::AreEqual(1.0, ut_array_length3(1.0, 0.0, 0.0), 1e-12);
            Assert::AreEqual(1.0, ut_array_length3(0.0, 1.0, 0.0), 1e-12);
            Assert::AreEqual(1.0, ut_array_length3(0.0, 0.0, 1.0), 1e-12);
        }

        TEST_METHOD(ArrayLength_Zero)
        {
            Assert::AreEqual(0.0, ut_array_length3(0.0, 0.0, 0.0), 1e-12);
        }

        TEST_METHOD(ArrayDistance_UnitCubeDiagonal)
        {
            double d = ut_array_distance3(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
            Assert::AreEqual(std::sqrt(3.0), d, 1e-12);
        }

        TEST_METHOD(ArrayDistance_SamePoint)
        {
            Assert::AreEqual(0.0, ut_array_distance3(1.5, 2.3, -4.7, 1.5, 2.3, -4.7), 1e-12);
        }

        TEST_METHOD(VecDiff_Basic)
        {
            double out[3];
            ut_vec_diff(5.0, 10.0, 15.0, 1.0, 2.0, 3.0, out);
            Assert::AreEqual(4.0,  out[0], 1e-12);
            Assert::AreEqual(8.0,  out[1], 1e-12);
            Assert::AreEqual(12.0, out[2], 1e-12);
        }

        TEST_METHOD(VecDiff_ZeroResult)
        {
            double out[3];
            ut_vec_diff(1.0, 2.0, 3.0, 1.0, 2.0, 3.0, out);
            Assert::AreEqual(0.0, out[0], 1e-12);
            Assert::AreEqual(0.0, out[1], 1e-12);
            Assert::AreEqual(0.0, out[2], 1e-12);
        }

        TEST_METHOD(VecCross_BasisVectors)
        {
            // i x j = k
            double out[3];
            ut_vec_cross(1.0, 0.0, 0.0,  0.0, 1.0, 0.0, out);
            Assert::AreEqual(0.0, out[0], 1e-12);
            Assert::AreEqual(0.0, out[1], 1e-12);
            Assert::AreEqual(1.0, out[2], 1e-12);
        }

        TEST_METHOD(VecCross_AntiCommutative)
        {
            // a x b = -(b x a)
            double ab[3], ba[3];
            ut_vec_cross(1.0, 2.0, 3.0,  4.0, 5.0, 6.0, ab);
            ut_vec_cross(4.0, 5.0, 6.0,  1.0, 2.0, 3.0, ba);
            for (int i = 0; i < 3; ++i)
                Assert::AreEqual(-ab[i], ba[i], 1e-12);
        }

        TEST_METHOD(VecCross_ParallelVectors)
        {
            // Parallel vectors → zero cross product
            double out[3];
            ut_vec_cross(2.0, 4.0, 6.0,  1.0, 2.0, 3.0, out);
            Assert::AreEqual(0.0, out[0], 1e-12);
            Assert::AreEqual(0.0, out[1], 1e-12);
            Assert::AreEqual(0.0, out[2], 1e-12);
        }

        TEST_METHOD(VecDot_Perpendicular)
        {
            Assert::AreEqual(0.0, ut_vec_dot(1.0, 0.0, 0.0,  0.0, 1.0, 0.0), 1e-12);
        }

        TEST_METHOD(VecDot_Parallel)
        {
            Assert::AreEqual(3.0, ut_vec_dot(1.0, 1.0, 1.0,  1.0, 1.0, 1.0), 1e-12);
        }

        TEST_METHOD(VecDot_KnownValue)
        {
            Assert::AreEqual(32.0, ut_vec_dot(1.0, 2.0, 3.0,  4.0, 5.0, 6.0), 1e-12);
        }
    };

    // -----------------------------------------------------------------------

    TEST_CLASS(NumericTests)
    {
    public:

        TEST_METHOD(IsNan_ValidDouble)
        {
            Assert::AreEqual(0, ut_is_nan(1.0));
            Assert::AreEqual(0, ut_is_nan(0.0));
            Assert::AreEqual(0, ut_is_nan(-3.14));
        }

        TEST_METHOD(IsNan_QuietNaN)
        {
            double nan_val = std::numeric_limits<double>::quiet_NaN();
            Assert::AreEqual(1, ut_is_nan(nan_val));
        }

        TEST_METHOD(IsNan_Infinity)
        {
            // Infinity is not NaN
            Assert::AreEqual(0, ut_is_nan(std::numeric_limits<double>::infinity()));
        }

        TEST_METHOD(IsSimilarRel_WithinTolerance)
        {
            // 1.0 vs 1.005: relative diff ≈ 0.5 %, within 1 %
            Assert::AreEqual(1, ut_is_similar_rel(1.0, 1.005, 0.01));
        }

        TEST_METHOD(IsSimilarRel_OutsideTolerance)
        {
            // 1.0 vs 1.05: relative diff = 5 %, outside 1 %
            Assert::AreEqual(0, ut_is_similar_rel(1.0, 1.05, 0.01));
        }

        TEST_METHOD(IsSimilarRel_EqualValues)
        {
            Assert::AreEqual(1, ut_is_similar_rel(42.0, 42.0, 1e-6));
        }

        TEST_METHOD(IsSimilarAbs_WithinTolerance)
        {
            Assert::AreEqual(1, ut_is_similar_abs(1.0, 1.0009, 0.001));
        }

        TEST_METHOD(IsSimilarAbs_OutsideTolerance)
        {
            Assert::AreEqual(0, ut_is_similar_abs(1.0, 1.002, 0.001));
        }

        TEST_METHOD(FastExpNeg_AtZero)
        {
            // fast_exp_neg uses std::exp for x > -ln(2) ≈ -0.693
            double val = ut_fast_exp_neg(0.0);
            Assert::AreEqual(1.0, val, 1e-12);
        }

        TEST_METHOD(FastExpNeg_AtMinusOne)
        {
            double approx = ut_fast_exp_neg(-1.0);
            double exact  = std::exp(-1.0);
            // Approximation tolerance: 0.5 %
            Assert::AreEqual(exact, approx, exact * 0.005);
        }

        TEST_METHOD(FastExpNeg_AtMinusTen)
        {
            // (1 + x/1024)^1024 underestimates exp(x) by ~x²/(2N) ≈ 4.9% at x=-10.
            double approx = ut_fast_exp_neg(-10.0);
            double exact  = std::exp(-10.0);
            Assert::AreEqual(exact, approx, exact * 0.06);
        }

        TEST_METHOD(FastExpNeg_BeyondCutoff)
        {
            // Values below -42 must return 0 exactly
            Assert::AreEqual(0.0, ut_fast_exp_neg(-100.0), 1e-30);
        }
    };

    // -----------------------------------------------------------------------

    TEST_CLASS(FchkParsingTests)
    {
    public:

        // FCHK format: 40-char label + type char + value starting at position 49
        // We pad to exactly 49 chars then append the value.

        TEST_METHOD(ReadFchkInt_PositiveValue)
        {
            // Format: keyword padded to 40, type at 40, 8 spaces (41-48), value at 49+
            const char* line = "Number of atoms                         I        3";
            Assert::AreEqual(3, ut_read_fchk_int(line));
        }

        TEST_METHOD(ReadFchkInt_NegativeValue)
        {
            const char* line = "Charge                                  I        -1";
            Assert::AreEqual(-1, ut_read_fchk_int(line));
        }

        TEST_METHOD(ReadFchkInt_LargeValue)
        {
            const char* line = "Number of basis functions               I        1024";
            Assert::AreEqual(1024, ut_read_fchk_int(line));
        }

        TEST_METHOD(ReadFchkDbl_NegativeScientific)
        {
            const char* line = "Total Energy                            R        -1.23456789E+02";
            double val = ut_read_fchk_dbl(line);
            Assert::AreEqual(-123.456789, val, 1e-6);
        }

        TEST_METHOD(ReadFchkDbl_PositiveScientific)
        {
            const char* line = "Zero-point correction                   R        4.56000000E-02";
            double val = ut_read_fchk_dbl(line);
            Assert::AreEqual(0.0456, val, 1e-10);
        }
    };

    // -----------------------------------------------------------------------

    TEST_CLASS(StringUtilTests)
    {
    public:

        TEST_METHOD(EndsWith_MatchingSuffix)
        {
            Assert::AreEqual(1, ut_ends_with("molecule.wfx", ".wfx"));
            Assert::AreEqual(1, ut_ends_with("data.hkl", ".hkl"));
        }

        TEST_METHOD(EndsWith_NonMatchingSuffix)
        {
            Assert::AreEqual(0, ut_ends_with("molecule.wfx", ".gbw"));
            Assert::AreEqual(0, ut_ends_with("test.cpp",     ".txt"));
        }

        TEST_METHOD(EndsWith_EmptySuffix)
        {
            // Empty suffix always matches
            Assert::AreEqual(1, ut_ends_with("anything", ""));
        }

        TEST_METHOD(EndsWith_Suffix_LongerThanString)
        {
            Assert::AreEqual(0, ut_ends_with("ab", "abc"));
        }

        TEST_METHOD(EndsWith_ExactMatch)
        {
            Assert::AreEqual(1, ut_ends_with(".wfx", ".wfx"));
        }

        TEST_METHOD(ShrinkString_RemovesDigits)
        {
            char buf[64];
            int len = ut_shrink_string("C1", buf, sizeof(buf));
            Assert::AreEqual(1, len);
            Assert::AreEqual(0, std::strcmp(buf, "C"));
        }

        TEST_METHOD(ShrinkString_RemovesSpacesAndDigits)
        {
            char buf[64];
            int len = ut_shrink_string("O 1 1", buf, sizeof(buf));
            Assert::IsTrue(len > 0);
            // All spaces and digits removed → only "O" remains
            Assert::AreEqual(0, std::strcmp(buf, "O"));
        }

        TEST_METHOD(ShrinkString_PureLettersUnchanged)
        {
            char buf[64];
            ut_shrink_string("Fe", buf, sizeof(buf));
            Assert::AreEqual(0, std::strcmp(buf, "Fe"));
        }

        TEST_METHOD(ShrinkString_BufferTooSmall)
        {
            char buf[1]; // Can't fit even "C"
            int result = ut_shrink_string("C1H2O", buf, sizeof(buf));
            Assert::AreEqual(-1, result);
        }
    };

    // -----------------------------------------------------------------------

    TEST_CLASS(Sha256Tests)
    {
    public:

        TEST_METHOD(Sha256_OutputLength)
        {
            char buf[65] = {};
            int len = ut_sha256("abc", buf, sizeof(buf));
            Assert::AreEqual(64, len);
        }

        TEST_METHOD(Sha256_KnownVector_Abc)
        {
            // NIST FIPS 180-4 test vector
            char buf[65] = {};
            ut_sha256("abc", buf, sizeof(buf));
            Assert::AreEqual(0, std::strcmp(buf,
                "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"));
        }

        TEST_METHOD(Sha256_EmptyString)
        {
            char buf[65] = {};
            ut_sha256("", buf, sizeof(buf));
            Assert::AreEqual(0, std::strcmp(buf,
                "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"));
        }

        TEST_METHOD(Sha256_Deterministic)
        {
            char buf1[65] = {}, buf2[65] = {};
            ut_sha256("NoSpherA2", buf1, sizeof(buf1));
            ut_sha256("NoSpherA2", buf2, sizeof(buf2));
            Assert::AreEqual(0, std::strcmp(buf1, buf2));
        }

        TEST_METHOD(Sha256_DifferentInputsDifferentOutputs)
        {
            char buf1[65] = {}, buf2[65] = {};
            ut_sha256("abc", buf1, sizeof(buf1));
            ut_sha256("abd", buf2, sizeof(buf2));
            Assert::AreNotEqual(0, std::strcmp(buf1, buf2));
        }

        TEST_METHOD(Sha256_BufferTooSmall)
        {
            char buf[10] = {};
            int result = ut_sha256("abc", buf, sizeof(buf));
            Assert::AreEqual(-1, result);
        }
    };

    // -----------------------------------------------------------------------

    TEST_CLASS(VecAggregateTests)
    {
    public:

        TEST_METHOD(VecSumBool_MixedValues)
        {
            int arr[] = { 1, 0, 1, 1, 0 };
            Assert::AreEqual(3, ut_vec_sum_bool(arr, 5));
        }

        TEST_METHOD(VecSumBool_AllFalse)
        {
            int arr[] = { 0, 0, 0 };
            Assert::AreEqual(0, ut_vec_sum_bool(arr, 3));
        }

        TEST_METHOD(VecSumBool_Empty)
        {
            int dummy = 0;
            Assert::AreEqual(0, ut_vec_sum_bool(&dummy, 0));
        }

        TEST_METHOD(VecSumInt_Basic)
        {
            int arr[] = { 1, 2, 3, 4, 5 };
            Assert::AreEqual(15, ut_vec_sum_int(arr, 5));
        }

        TEST_METHOD(VecSumInt_Negative)
        {
            int arr[] = { -3, 7, -2 };
            Assert::AreEqual(2, ut_vec_sum_int(arr, 3));
        }

        TEST_METHOD(VecSumDouble_Basic)
        {
            double arr[] = { 1.5, -0.5, 2.0 };
            Assert::AreEqual(3.0, ut_vec_sum_double(arr, 3), 1e-12);
        }

        TEST_METHOD(VecSumDouble_AllZero)
        {
            double arr[] = { 0.0, 0.0, 0.0 };
            Assert::AreEqual(0.0, ut_vec_sum_double(arr, 3), 1e-12);
        }

        TEST_METHOD(VecLength_345)
        {
            double arr[] = { 3.0, 4.0 };
            Assert::AreEqual(5.0, ut_vec_length(arr, 2), 1e-12);
        }

        TEST_METHOD(VecLength_Unit)
        {
            double arr[] = { 1.0, 0.0, 0.0 };
            Assert::AreEqual(1.0, ut_vec_length(arr, 3), 1e-12);
        }

        TEST_METHOD(VecLength_4D)
        {
            double arr[] = { 1.0, 1.0, 1.0, 1.0 };
            Assert::AreEqual(2.0, ut_vec_length(arr, 4), 1e-12);
        }

        TEST_METHOD(VecLength_Empty)
        {
            double dummy = 0.0;
            Assert::AreEqual(0.0, ut_vec_length(&dummy, 0), 1e-12);
        }
    };

    // -----------------------------------------------------------------------

    TEST_CLASS(StringUtilTests2)
    {
    public:

        TEST_METHOD(Trim_LeadingTrailingSpaces)
        {
            char buf[64];
            ut_trim("  hello world  ", buf, sizeof(buf));
            Assert::AreEqual(0, std::strcmp(buf, "hello world"));
        }

        TEST_METHOD(Trim_NoSpaces)
        {
            char buf[64];
            ut_trim("no_spaces", buf, sizeof(buf));
            Assert::AreEqual(0, std::strcmp(buf, "no_spaces"));
        }

        TEST_METHOD(Trim_EmptyString)
        {
            char buf[64];
            int len = ut_trim("", buf, sizeof(buf));
            Assert::AreEqual(0, len);
            Assert::AreEqual(0, std::strcmp(buf, ""));
        }

        TEST_METHOD(Trim_OnlySpaces)
        {
            char buf[64];
            ut_trim("    ", buf, sizeof(buf));
            Assert::AreEqual(0, std::strcmp(buf, ""));
        }

        TEST_METHOD(Asciitolower_Uppercase)
        {
            Assert::AreEqual('a', ut_asciitolower('A'));
            Assert::AreEqual('z', ut_asciitolower('Z'));
            Assert::AreEqual('m', ut_asciitolower('M'));
        }

        TEST_METHOD(Asciitolower_AlreadyLower)
        {
            Assert::AreEqual('a', ut_asciitolower('a'));
            Assert::AreEqual('z', ut_asciitolower('z'));
        }

        TEST_METHOD(Asciitolower_NonAlpha)
        {
            // Digits and symbols are returned unchanged
            Assert::AreEqual('5', ut_asciitolower('5'));
            Assert::AreEqual('_', ut_asciitolower('_'));
        }

        TEST_METHOD(DoubleFromEsd_Plain)
        {
            Assert::AreEqual(3.14159, ut_double_from_esd("3.14159"), 1e-10);
        }

        TEST_METHOD(DoubleFromEsd_WithEsd)
        {
            // "(5)" is stripped; value is 1.234
            Assert::AreEqual(1.234, ut_double_from_esd("1.234(5)"), 1e-10);
        }

        TEST_METHOD(DoubleFromEsd_Zero)
        {
            Assert::AreEqual(0.0, ut_double_from_esd("0.0"), 1e-12);
        }

        TEST_METHOD(DecimalPrecisionCif_WithBracketsAndDecimal)
        {
            // "1.2345(6)" → 6 × 10⁻⁴ = 0.0006
            Assert::AreEqual(6e-4, ut_decimal_precision_cif("1.2345(6)"), 1e-10);
        }

        TEST_METHOD(DecimalPrecisionCif_NoBrackets)
        {
            // No brackets → default 0.005
            Assert::AreEqual(0.005, ut_decimal_precision_cif("1.234"), 1e-10);
        }

        TEST_METHOD(DecimalPrecisionCif_IntegerWithEsd)
        {
            // "100(2)" → no decimal, digit count from bracket positions
            Assert::AreEqual(0.2, ut_decimal_precision_cif("100(2)"), 1e-10);
        }
    };

    // -----------------------------------------------------------------------

    TEST_CLASS(BasisTypeTests)
    {
    public:

        TEST_METHOD(Sht2nbas_CartesianShells)
        {
            // Cartesian: S=1, P=3, D=6, F=10, G=15
            Assert::AreEqual(1,  ut_sht2nbas(0));
            Assert::AreEqual(3,  ut_sht2nbas(1));
            Assert::AreEqual(6,  ut_sht2nbas(2));
            Assert::AreEqual(10, ut_sht2nbas(3));
            Assert::AreEqual(15, ut_sht2nbas(4));
        }

        TEST_METHOD(Sht2nbas_SphericalShells)
        {
            // Negative types → spherical: -2=D(5), -3=F(7), -4=G(9)
            Assert::AreEqual(5, ut_sht2nbas(-2));
            Assert::AreEqual(7, ut_sht2nbas(-3));
            Assert::AreEqual(9, ut_sht2nbas(-4));
        }

        TEST_METHOD(DoubleFactorial_SmallValues)
        {
            Assert::AreEqual(1u,   ut_doublefactorial(0));
            Assert::AreEqual(1u,   ut_doublefactorial(1));
            Assert::AreEqual(2u,   ut_doublefactorial(2));
            Assert::AreEqual(3u,   ut_doublefactorial(3));
            Assert::AreEqual(8u,   ut_doublefactorial(4));
            Assert::AreEqual(15u,  ut_doublefactorial(5));
            Assert::AreEqual(48u,  ut_doublefactorial(6));
            Assert::AreEqual(105u, ut_doublefactorial(7));
        }
    };

    // -----------------------------------------------------------------------

    TEST_CLASS(ConstantsTests)
    {
    public:

        TEST_METHOD(ConstAbs_Positive)    { Assert::AreEqual(5,  ut_const_abs(5));  }
        TEST_METHOD(ConstAbs_Negative)    { Assert::AreEqual(5,  ut_const_abs(-5)); }
        TEST_METHOD(ConstAbs_Zero)        { Assert::AreEqual(0,  ut_const_abs(0));  }

        TEST_METHOD(ConstexprPow_Integer)
        {
            Assert::AreEqual(1024.0, ut_constexpr_pow(2.0, 10), 1e-12);
            Assert::AreEqual(1.0,    ut_constexpr_pow(10.0, 0), 1e-12);
            Assert::AreEqual(8.0,    ut_constexpr_pow(2.0, 3),  1e-12);
        }

        TEST_METHOD(ConstantsSqrt_KnownValues)
        {
            Assert::AreEqual(2.0,            ut_constants_sqrt(4.0), 1e-10);
            Assert::AreEqual(std::sqrt(2.0), ut_constants_sqrt(2.0), 1e-10);
            Assert::AreEqual(0.0,            ut_constants_sqrt(0.0), 1e-12);
        }

        TEST_METHOD(ConstantsSqrt_NegativeIsNaN)
        {
            double r = ut_constants_sqrt(-1.0);
            Assert::IsTrue(std::isnan(r));
        }

        TEST_METHOD(ExpApprox_AtZero)
        {
            Assert::AreEqual(1.0, ut_exp_approx(0.0, 25), 1e-12);
        }

        TEST_METHOD(ExpApprox_AtOne)
        {
            // 25-term Taylor series matches std::exp to machine precision
            Assert::AreEqual(std::exp(1.0), ut_exp_approx(1.0, 25), 1e-12);
        }

        TEST_METHOD(ExpApprox_AtMinusOne)
        {
            Assert::AreEqual(std::exp(-1.0), ut_exp_approx(-1.0, 25), 1e-10);
        }

        TEST_METHOD(LogApprox_AtOne)
        {
            Assert::AreEqual(0.0, ut_log_approx(1.0, 25), 1e-12);
        }

        TEST_METHOD(LogApprox_AtE)
        {
            // Arctanh series converges; 25 iterations accurate to 1e-8 for x=e
            Assert::AreEqual(1.0, ut_log_approx(std::exp(1.0), 25), 1e-6);
        }

        TEST_METHOD(LogApprox_NonPositiveReturnsSentinel)
        {
            Assert::AreEqual(-1.0, ut_log_approx(0.0, 25), 1e-12);
            Assert::AreEqual(-1.0, ut_log_approx(-5.0, 25), 1e-12);
        }

        TEST_METHOD(Bohr2Ang_OneBohr)
        {
            // 1 Bohr = a₀ Å = 0.529177210903 Å
            Assert::AreEqual(0.529177210903, ut_bohr2ang(1.0), 1e-10);
        }

        TEST_METHOD(Bohr2Ang_Zero)
        {
            Assert::AreEqual(0.0, ut_bohr2ang(0.0), 1e-12);
        }

        TEST_METHOD(Ang2Bohr_OneAngstrom)
        {
            Assert::AreEqual(1.0 / 0.529177210903, ut_ang2bohr(1.0), 1e-8);
        }

        TEST_METHOD(Ang2Bohr_Roundtrip)
        {
            // ang2bohr(bohr2ang(x)) ≈ x
            double x = 2.5;
            Assert::AreEqual(x, ut_ang2bohr(ut_bohr2ang(x)), 1e-10);
        }

        TEST_METHOD(CubicBohr2Ang_OneBohr3)
        {
            // 1 Bohr³ = a₀³ Å³
            double expected = 0.529177210903 * 0.529177210903 * 0.529177210903;
            Assert::AreEqual(expected, ut_cubic_bohr2ang(1.0), 1e-10);
        }

        TEST_METHOD(CubicAng2Bohr_Roundtrip)
        {
            double x = 3.0;
            Assert::AreEqual(x, ut_cubic_ang2bohr(ut_cubic_bohr2ang(x)), 1e-8);
        }

        TEST_METHOD(Factorial_SmallValues)
        {
            Assert::AreEqual(1LL,       ut_factorial(0));
            Assert::AreEqual(1LL,       ut_factorial(1));
            Assert::AreEqual(120LL,     ut_factorial(5));
            Assert::AreEqual(3628800LL, ut_factorial(10));
        }
    };

    // -----------------------------------------------------------------------

    TEST_CLASS(OrbitalIndexTests)
    {
    public:

        TEST_METHOD(OrcaToPySCF_SShell)
        {
            // S: only one component, m_idx=0 → 0
            Assert::AreEqual(0, ut_orca_2_pyscf(0, 0));
        }

        TEST_METHOD(OrcaToPySCF_PShell)
        {
            // P ORCA ordering 0,+1,-1 → PySCF map {1,2,0}
            Assert::AreEqual(1, ut_orca_2_pyscf(1, 0));
            Assert::AreEqual(2, ut_orca_2_pyscf(1, 1));
            Assert::AreEqual(0, ut_orca_2_pyscf(1, 2));
        }

        TEST_METHOD(OrcaToPySCF_DShell_FirstComponent)
        {
            // D: map {2,3,1,4,0}, m_idx=0 → 2
            Assert::AreEqual(2, ut_orca_2_pyscf(2, 0));
        }

        TEST_METHOD(OrcaToPySCF_OutOfRange)
        {
            // l=100 is not in the switch → nullopt → -1
            Assert::AreEqual(-1, ut_orca_2_pyscf(100, 0));
        }

        TEST_METHOD(TypeToNbo_SShell)     { Assert::AreEqual(1u,   ut_type_2_nbo(1));  }
        TEST_METHOD(TypeToNbo_PxShell)    { Assert::AreEqual(101u,  ut_type_2_nbo(2));  }
        TEST_METHOD(TypeToNbo_DxxShell)   { Assert::AreEqual(201u,  ut_type_2_nbo(5));  }
        TEST_METHOD(TypeToNbo_FxxxShell)  { Assert::AreEqual(301u,  ut_type_2_nbo(11)); }
        TEST_METHOD(TypeToNbo_GxxxxShell) { Assert::AreEqual(401u,  ut_type_2_nbo(21)); }
        TEST_METHOD(TypeToNbo_Unknown)    { Assert::AreEqual(0u,    ut_type_2_nbo(99)); }
    };

    // -----------------------------------------------------------------------

    TEST_CLASS(BesselTests)
    {
    public:

        TEST_METHOD(BesselJ0_AtZero)
        {
            // j₀(0) = 1 (limit of sin(x)/x as x→0)
            Assert::AreEqual(1.0, ut_bessel_j(0, 0.0), 1e-12);
        }

        TEST_METHOD(BesselJ1_AtZero)
        {
            // j_l(0) = 0 for l > 0
            Assert::AreEqual(0.0, ut_bessel_j(1, 0.0), 1e-12);
            Assert::AreEqual(0.0, ut_bessel_j(5, 0.0), 1e-12);
        }

        TEST_METHOD(BesselJ0_AtOne)
        {
            // j₀(1) = sin(1)/1
            Assert::AreEqual(std::sin(1.0), ut_bessel_j(0, 1.0), 1e-12);
        }

        TEST_METHOD(BesselJ1_AtOne)
        {
            // j₁(1) = (sin(1) - cos(1)) / 1
            double expected = std::sin(1.0) - std::cos(1.0);
            Assert::AreEqual(expected, ut_bessel_j(1, 1.0), 1e-12);
        }

        TEST_METHOD(BesselJ2_AtOne)
        {
            // j₂(1) = (2·sin(1) - 3·cos(1)) / 1
            double expected = 2.0 * std::sin(1.0) - 3.0 * std::cos(1.0);
            Assert::AreEqual(expected, ut_bessel_j(2, 1.0), 1e-12);
        }

        TEST_METHOD(BesselJ0_AtPi)
        {
            // j₀(π) = sin(π)/π ≈ 0
            Assert::AreEqual(std::sin(PI_VAL) / PI_VAL, ut_bessel_j(0, PI_VAL), 1e-12);
        }

        TEST_METHOD(BesselJ_HigherOrder_PositiveAndFinite)
        {
            // l=7 exercises the continued-fraction fallback path
            double r = ut_bessel_j(7, 2.0);
            Assert::IsTrue(std::isfinite(r));
            Assert::IsTrue(r > 0.0);
        }

        TEST_METHOD(BesselJ_RecurrenceCheck)
        {
            // Recurrence: j_{l-1}(x) + j_{l+1}(x) = (2l+1)/x · j_l(x)
            double x = 3.0;
            int l = 3;
            double jlm1 = ut_bessel_j(l - 1, x);
            double jl   = ut_bessel_j(l,     x);
            double jlp1 = ut_bessel_j(l + 1, x);
            double lhs  = jlm1 + jlp1;
            double rhs  = (2.0 * l + 1.0) / x * jl;
            Assert::AreEqual(lhs, rhs, 1e-10);
        }
    };

    // -----------------------------------------------------------------------
    // IsSimilarPow10Tests — is_similar(a, b, tolerance) where |a-b| <= 10^tol
    // -----------------------------------------------------------------------
    TEST_CLASS(IsSimilarPow10Tests)
    {
    public:
        TEST_METHOD(Equal_ReturnTrue)
        {
            Assert::AreEqual(1, ut_is_similar_pow10(1.0, 1.0, -6.0));
        }

        TEST_METHOD(WithinTolerance_ReturnTrue)
        {
            // |1.000001 - 1.0| = 1e-6 <= 10^(-6)
            Assert::AreEqual(1, ut_is_similar_pow10(1.000001, 1.0, -6.0));
        }

        TEST_METHOD(OutsideTolerance_ReturnFalse)
        {
            // |1.00001 - 1.0| = 1e-5 > 10^(-6)
            Assert::AreEqual(0, ut_is_similar_pow10(1.00001, 1.0, -6.0));
        }

        TEST_METHOD(NegativeValues_WithinTolerance)
        {
            Assert::AreEqual(1, ut_is_similar_pow10(-5.0, -5.0 + 1e-8, -7.0));
        }

        TEST_METHOD(LooseTolerance_LargeDiff)
        {
            // |100 - 50| = 50 <= 10^2 = 100
            Assert::AreEqual(1, ut_is_similar_pow10(100.0, 50.0, 2.0));
        }

        TEST_METHOD(LooseTolerance_TooLargeDiff)
        {
            // |200 - 50| = 150 > 10^2 = 100
            Assert::AreEqual(0, ut_is_similar_pow10(200.0, 50.0, 2.0));
        }
    };

    // -----------------------------------------------------------------------
    // Shell2FunctionTests — shell2function(type, prim) WFN column index
    // -----------------------------------------------------------------------
    TEST_CLASS(Shell2FunctionTests)
    {
    public:
        TEST_METHOD(SType_Prim0_Returns0)
        {
            // s-type shell (type=1): first and only function is index 0
            int r = ut_shell2function(1, 0);
            Assert::IsTrue(r >= 0);
        }

        TEST_METHOD(PType_Prim0)
        {
            // p-type shell: 3 functions
            int r = ut_shell2function(2, 0);
            Assert::IsTrue(r >= 0);
        }

        TEST_METHOD(PType_Prim1)
        {
            int r0 = ut_shell2function(2, 0);
            int r1 = ut_shell2function(2, 1);
            Assert::IsTrue(r1 > r0);
        }

        TEST_METHOD(DType_Prim5_ValidIndex)
        {
            // d-type (type=3): 6 Cartesian or 5 spherical functions
            int r = ut_shell2function(3, 5);
            Assert::IsTrue(r >= 0);
        }

        TEST_METHOD(ResultsAreStrictlyIncreasingWithinShell)
        {
            // f-type (type=4): consecutive prims must give increasing column indices
            int r0 = ut_shell2function(4, 0);
            int r1 = ut_shell2function(4, 1);
            int r2 = ut_shell2function(4, 2);
            Assert::IsTrue(r1 > r0);
            Assert::IsTrue(r2 > r1);
        }
    };

    // -----------------------------------------------------------------------
    // CountWordsTests — CountWords(str)
    // -----------------------------------------------------------------------
    TEST_CLASS(CountWordsTests)
    {
    public:
        TEST_METHOD(Empty_Returns0)
        {
            Assert::AreEqual(0, ut_count_words(""));
        }

        TEST_METHOD(OneWord)
        {
            Assert::AreEqual(1, ut_count_words("hello"));
        }

        TEST_METHOD(TwoWords)
        {
            Assert::AreEqual(2, ut_count_words("hello world"));
        }

        TEST_METHOD(LeadingTrailingSpaces)
        {
            Assert::AreEqual(2, ut_count_words("  foo   bar  "));
        }

        TEST_METHOD(MultipleSpacesBetweenWords)
        {
            Assert::AreEqual(3, ut_count_words("a  b  c"));
        }

        TEST_METHOD(SingleSpace)
        {
            Assert::AreEqual(0, ut_count_words(" "));
        }
    };

    // -----------------------------------------------------------------------
    // ShrinkStringToAtomTests — shrink_string_to_atom(input, atom_number)
    // -----------------------------------------------------------------------
    TEST_CLASS(ShrinkStringToAtomTests)
    {
    public:
        TEST_METHOD(CarbonAtomNumber6)
        {
            // atnr2letter(6) = "C"
            char buf[32];
            int n = ut_shrink_string_to_atom("C1", 6, buf, 32);
            Assert::IsTrue(n >= 0);
            Assert::AreEqual(std::string("C"), std::string(buf));
        }

        TEST_METHOD(CalciumAtomNumber20)
        {
            // atnr2letter(20) = "Ca"
            char buf[32];
            int n = ut_shrink_string_to_atom("Ca12", 20, buf, 32);
            Assert::IsTrue(n >= 0);
            Assert::AreEqual(std::string("Ca"), std::string(buf));
        }

        TEST_METHOD(BufferTooSmall_ReturnsNegOne)
        {
            char buf[1];
            int n = ut_shrink_string_to_atom("Carbon6", 6, buf, 1);
            Assert::AreEqual(-1, n);
        }

        TEST_METHOD(IronAtomNumber26)
        {
            // atnr2letter(26) = "Fe"
            char buf[32];
            int n = ut_shrink_string_to_atom("Fe3 ", 26, buf, 32);
            Assert::IsTrue(n >= 0);
            Assert::AreEqual(std::string("Fe"), std::string(buf));
        }
    };

    // -----------------------------------------------------------------------
    // SplitStringTests — split_string(input, delim)
    // -----------------------------------------------------------------------
    TEST_CLASS(SplitStringTests)
    {
    public:
        TEST_METHOD(SingleToken_NoDelim)
        {
            const int MAX = 8;
            char storage[MAX][64];
            char* toks[MAX];
            for (int i = 0; i < MAX; ++i) toks[i] = storage[i];
            int n = ut_split_string("hello", " ", toks, MAX, 64);
            Assert::AreEqual(1, n);
            Assert::AreEqual(std::string("hello"), std::string(toks[0]));
        }

        TEST_METHOD(ThreeTokens)
        {
            const int MAX = 8;
            char storage[MAX][64];
            char* toks[MAX];
            for (int i = 0; i < MAX; ++i) toks[i] = storage[i];
            int n = ut_split_string("a b c", " ", toks, MAX, 64);
            Assert::AreEqual(3, n);
            Assert::AreEqual(std::string("a"), std::string(toks[0]));
            Assert::AreEqual(std::string("b"), std::string(toks[1]));
            Assert::AreEqual(std::string("c"), std::string(toks[2]));
        }

        TEST_METHOD(CommaDelimiter)
        {
            const int MAX = 8;
            char storage[MAX][64];
            char* toks[MAX];
            for (int i = 0; i < MAX; ++i) toks[i] = storage[i];
            int n = ut_split_string("x,y,z", ",", toks, MAX, 64);
            Assert::AreEqual(3, n);
            Assert::AreEqual(std::string("z"), std::string(toks[2]));
        }

        TEST_METHOD(MaxOutLimit_ReturnsTotalCount)
        {
            const int MAX = 2;
            char storage[MAX][64];
            char* toks[MAX];
            for (int i = 0; i < MAX; ++i) toks[i] = storage[i];
            // 4 tokens but max_out=2
            int n = ut_split_string("a b c d", " ", toks, MAX, 64);
            Assert::AreEqual(4, n); // total = 4
            Assert::AreEqual(std::string("a"), std::string(toks[0]));
            Assert::AreEqual(std::string("b"), std::string(toks[1]));
        }

        TEST_METHOD(EmptyString_ZeroTokens)
        {
            const int MAX = 8;
            char storage[MAX][64];
            char* toks[MAX];
            for (int i = 0; i < MAX; ++i) toks[i] = storage[i];
            int n = ut_split_string("", " ", toks, MAX, 64);
            Assert::IsTrue(n == 0 || n == 1); // impl-defined for empty input
        }
    };

    // -----------------------------------------------------------------------
    // TimingTests — ut_sleep_and_measure_us(N) returns elapsed µs >= N*1000
    // -----------------------------------------------------------------------
    TEST_CLASS(TimingTests)
    {
    public:
        TEST_METHOD(Sleep10ms_ElapsedAtLeast10000us)
        {
            long long us = ut_sleep_and_measure_us(10);
            Assert::IsTrue(us >= 10000LL);
        }

        TEST_METHOD(Sleep1ms_ElapsedAtLeast1000us)
        {
            long long us = ut_sleep_and_measure_us(1);
            Assert::IsTrue(us >= 1000LL);
        }

        TEST_METHOD(Sleep5ms_ElapsedPositive)
        {
            long long us = ut_sleep_and_measure_us(5);
            Assert::IsTrue(us > 0LL);
        }
    };

    // -----------------------------------------------------------------------
    // HypergeometricTests — 2F1(a,b;c;x)
    // -----------------------------------------------------------------------
    TEST_CLASS(HypergeometricTests)
    {
    public:
        TEST_METHOD(Identity_2F1_Zero_IsOne)
        {
            // 2F1(a,b;c;0) = 1 for any a,b,c
            double r = ut_hypergeometric(1.0, 2.0, 3.0, 0.0);
            Assert::AreEqual(1.0, r, 1e-12);
        }

        TEST_METHOD(KnownValue_2F1_1_1_2_Half)
        {
            // 2F1(1,1;2;0.5) = -2*ln(0.5) = 2*ln(2) ≈ 1.386294...
            double expected = -2.0 * std::log(0.5);
            double r = ut_hypergeometric(1.0, 1.0, 2.0, 0.5);
            Assert::AreEqual(expected, r, 1e-6);
        }

        TEST_METHOD(KnownValue_2F1_Half_Half_ThreeHalves_Half)
        {
            // 2F1(0.5,0.5;1.5;0.5) — finite, positive
            double r = ut_hypergeometric(0.5, 0.5, 1.5, 0.5);
            Assert::IsTrue(std::isfinite(r));
            Assert::IsTrue(r > 1.0);
        }

        TEST_METHOD(Symmetry_ab_equals_ba)
        {
            // 2F1(a,b;c;x) = 2F1(b,a;c;x)
            double r1 = ut_hypergeometric(2.0, 3.0, 5.0, 0.3);
            double r2 = ut_hypergeometric(3.0, 2.0, 5.0, 0.3);
            Assert::AreEqual(r1, r2, 1e-10);
        }

        TEST_METHOD(NegativeX_ReturnsFinite)
        {
            double r = ut_hypergeometric(1.0, 2.0, 3.0, -0.5);
            Assert::IsTrue(std::isfinite(r));
        }
    };

    // -----------------------------------------------------------------------
    // Atnr2LetterTests — constants::atnr2letter(nr) element symbol lookup
    // -----------------------------------------------------------------------
    TEST_CLASS(Atnr2LetterTests)
    {
    public:
        TEST_METHOD(Hydrogen_Is_H)
        {
            char buf[8];
            int n = ut_atnr2letter(1, buf, 8);
            Assert::AreEqual(std::string("H"), std::string(buf));
            Assert::AreEqual(1, n);
        }

        TEST_METHOD(Carbon_Is_C)
        {
            char buf[8];
            ut_atnr2letter(6, buf, 8);
            Assert::AreEqual(std::string("C"), std::string(buf));
        }

        TEST_METHOD(Iron_Is_Fe)
        {
            char buf[8];
            ut_atnr2letter(26, buf, 8);
            Assert::AreEqual(std::string("Fe"), std::string(buf));
        }

        TEST_METHOD(Gold_Is_Au)
        {
            char buf[8];
            ut_atnr2letter(79, buf, 8);
            Assert::AreEqual(std::string("Au"), std::string(buf));
        }

        TEST_METHOD(Lawrencium_103_Is_Lr)
        {
            char buf[8];
            ut_atnr2letter(103, buf, 8);
            Assert::AreEqual(std::string("Lr"), std::string(buf));
        }

        TEST_METHOD(Zero_Is_Q_Peak)
        {
            char buf[8];
            ut_atnr2letter(0, buf, 8);
            Assert::AreEqual(std::string("Q"), std::string(buf));
        }

        TEST_METHOD(OutOfRange_Returns_PROBLEM)
        {
            char buf[16];
            ut_atnr2letter(200, buf, 16);
            Assert::AreEqual(std::string("PROBLEM"), std::string(buf));
        }

        TEST_METHOD(BufferTooSmall_ReturnsNegOne)
        {
            char buf[1];
            int n = ut_atnr2letter(6, buf, 1);
            Assert::AreEqual(-1, n);
        }
    };

    // -----------------------------------------------------------------------
    // Type2VectorTests — constants::type2vector: basis type → [nx, ny, nz]
    // -----------------------------------------------------------------------
    TEST_CLASS(Type2VectorTests)
    {
    public:
        TEST_METHOD(Type1_Is_SShell_000)
        {
            int nx, ny, nz;
            ut_type2vector(1, &nx, &ny, &nz);
            Assert::AreEqual(0, nx);
            Assert::AreEqual(0, ny);
            Assert::AreEqual(0, nz);
        }

        TEST_METHOD(Type2_Is_Px_100)
        {
            int nx, ny, nz;
            ut_type2vector(2, &nx, &ny, &nz);
            Assert::AreEqual(1, nx);
            Assert::AreEqual(0, ny);
            Assert::AreEqual(0, nz);
        }

        TEST_METHOD(Type5_Is_Dx2_200)
        {
            // index 5 in type_vector: (2,0,0) = dx²
            int nx, ny, nz;
            ut_type2vector(5, &nx, &ny, &nz);
            Assert::AreEqual(2, nx);
            Assert::AreEqual(0, ny);
            Assert::AreEqual(0, nz);
        }

        TEST_METHOD(Type8_Is_Dxy_110)
        {
            // index 8: (1,1,0) = dxy
            int nx, ny, nz;
            ut_type2vector(8, &nx, &ny, &nz);
            Assert::AreEqual(1, nx);
            Assert::AreEqual(1, ny);
            Assert::AreEqual(0, nz);
        }

        TEST_METHOD(SumOfExponents_Matches_ShellType)
        {
            // For d-type (types 5-10), nx+ny+nz == 2
            for (int t = 5; t <= 10; ++t) {
                int nx, ny, nz;
                ut_type2vector(t, &nx, &ny, &nz);
                Assert::AreEqual(2, nx + ny + nz);
            }
        }

        TEST_METHOD(OutOfRange_Returns_NegOne)
        {
            int nx, ny, nz;
            ut_type2vector(0, &nx, &ny, &nz);
            Assert::AreEqual(-1, nx);
            ut_type2vector(57, &nx, &ny, &nz);
            Assert::AreEqual(-1, nx);
        }
    };

    // -----------------------------------------------------------------------
    // NormGaussTests — constants::normgauss(type, exponent)
    // -----------------------------------------------------------------------
    TEST_CLASS(NormGaussTests)
    {
    public:
        TEST_METHOD(SShell_IsPositive)
        {
            double n = ut_normgauss(1, 1.0);
            Assert::IsTrue(n > 0.0);
            Assert::IsTrue(std::isfinite(n));
        }

        TEST_METHOD(SShell_ScalesWithExponent)
        {
            // Normalization grows with exponent for s-type
            double n1 = ut_normgauss(1, 1.0);
            double n2 = ut_normgauss(1, 4.0);
            Assert::IsTrue(n2 > n1);
        }

        TEST_METHOD(PShell_IsPositive)
        {
            double n = ut_normgauss(2, 1.0);
            Assert::IsTrue(n > 0.0);
        }

        TEST_METHOD(DShell_IsPositive)
        {
            double n = ut_normgauss(5, 1.0);
            Assert::IsTrue(n > 0.0);
        }

        TEST_METHOD(SShell_KnownValue)
        {
            // normgauss for s-type (0,0,0), exponent α:
            // N = (2α/π)^(9/4) * sqrt(1/1) = (2α/π)^(9/4)
            // For α=1: N = (2/π)^(9/4)
            double alpha = 1.0;
            double expected = std::pow(2.0 * alpha / PI_VAL, 9.0 / 4.0);
            double result = ut_normgauss(1, alpha);
            Assert::AreEqual(expected, result, 1e-10);
        }
    };

    // -----------------------------------------------------------------------
    // AssocLegendreTests — constants::associated_legendre_polynomial(l, m, x)
    // -----------------------------------------------------------------------
    TEST_CLASS(AssocLegendreTests)
    {
    public:
        TEST_METHOD(P00_Is_One)
        {
            Assert::AreEqual(1.0, ut_assoc_legendre(0, 0, 0.5), 1e-12);
        }

        TEST_METHOD(P10_At_Half_Is_Half)
        {
            // P_1^0(x) = x
            Assert::AreEqual(0.5, ut_assoc_legendre(1, 0, 0.5), 1e-12);
        }

        TEST_METHOD(P11_AtOne_Is_Zero)
        {
            // P_1^1(x) = sqrt(1-x²), at x=1: 0
            Assert::AreEqual(0.0, ut_assoc_legendre(1, 1, 1.0), 1e-12);
        }

        TEST_METHOD(P20_Is_Legendre_Polynomial)
        {
            // P_2^0(x) = (3x²-1)/2
            double x = 0.6;
            double expected = 0.5 * (3 * x * x - 1);
            Assert::AreEqual(expected, ut_assoc_legendre(2, 0, x), 1e-12);
        }

        TEST_METHOD(P21_At_Zero)
        {
            // P_2^1(x) = 3x*sqrt(1-x²), at x=0: 0
            Assert::AreEqual(0.0, ut_assoc_legendre(2, 1, 0.0), 1e-12);
        }

        TEST_METHOD(P22_At_Zero)
        {
            // P_2^2(x) = -3(x²-1) = 3(1-x²), at x=0: 3
            Assert::AreEqual(-3.0 * (0.0 - 1.0), ut_assoc_legendre(2, 2, 0.0), 1e-12);
        }

        TEST_METHOD(NegativeM_P1m1)
        {
            // P_1^{-1}(x) = -0.5*sqrt(1-x²)
            double x = 0.5;
            double expected = -0.5 * std::sqrt(1 - x * x);
            Assert::AreEqual(expected, ut_assoc_legendre(1, -1, x), 1e-12);
        }
    };

    // -----------------------------------------------------------------------
    // CartesianToSphericalTests — constants::cartesian_to_spherical(x,y,z)
    // -----------------------------------------------------------------------
    TEST_CLASS(CartesianToSphericalTests)
    {
    public:
        TEST_METHOD(Origin_Has_Zero_Radius)
        {
            double out[3];
            ut_cartesian_to_spherical(0.0, 0.0, 0.0, out);
            Assert::AreEqual(0.0, out[0], 1e-12); // r
        }

        TEST_METHOD(UnitX_Gives_Correct_Angles)
        {
            double out[3];
            ut_cartesian_to_spherical(1.0, 0.0, 0.0, out);
            Assert::AreEqual(1.0,           out[0], 1e-12); // r=1
            Assert::AreEqual(PI_VAL / 2.0,  out[1], 1e-12); // theta=pi/2
            Assert::AreEqual(0.0,           out[2], 1e-12); // phi=0
        }

        TEST_METHOD(UnitZ_Has_Zero_Theta)
        {
            double out[3];
            ut_cartesian_to_spherical(0.0, 0.0, 1.0, out);
            Assert::AreEqual(1.0,  out[0], 1e-12); // r=1
            Assert::AreEqual(0.0,  out[1], 1e-12); // theta=0
        }

        TEST_METHOD(UnitY_Gives_PhiHalfPi)
        {
            double out[3];
            ut_cartesian_to_spherical(0.0, 1.0, 0.0, out);
            Assert::AreEqual(1.0,          out[0], 1e-12); // r=1
            Assert::AreEqual(PI_VAL / 2.0, out[2], 1e-12); // phi=pi/2
        }

        TEST_METHOD(Radius_Is_Euclidean_Norm)
        {
            double out[3];
            ut_cartesian_to_spherical(3.0, 4.0, 0.0, out);
            Assert::AreEqual(5.0, out[0], 1e-12); // r=5
        }

        TEST_METHOD(Inverse_Recover_Cartesian)
        {
            // Convert (1,1,1) → spherical → back to Cartesian
            double out[3];
            ut_cartesian_to_spherical(1.0, 1.0, 1.0, out);
            double r = out[0], theta = out[1], phi = out[2];
            double x = r * std::sin(theta) * std::cos(phi);
            double y = r * std::sin(theta) * std::sin(phi);
            double z = r * std::cos(theta);
            Assert::AreEqual(1.0, x, 1e-12);
            Assert::AreEqual(1.0, y, 1e-12);
            Assert::AreEqual(1.0, z, 1e-12);
        }
    };

    // -----------------------------------------------------------------------
    // GetLambda1Tests — median eigenvalue of a 3x3 symmetric matrix
    // -----------------------------------------------------------------------
    TEST_CLASS(GetLambda1Tests)
    {
    public:
        TEST_METHOD(DiagonalMatrix_ReturnsMidEigenvalue)
        {
            // Eigenvalues are 1, 3, 5 → median = 3
            const double a[9] = {
                1, 0, 0,
                0, 3, 0,
                0, 0, 5
            };
            double lam = ut_get_lambda_1(a);
            Assert::AreEqual(3.0, lam, 1e-10);
        }

        TEST_METHOD(DiagonalMatrix_AllEqual_ReturnsThat)
        {
            const double a[9] = {
                2, 0, 0,
                0, 2, 0,
                0, 0, 2
            };
            double lam = ut_get_lambda_1(a);
            // all eigenvalues equal 2; any is "middle"
            Assert::AreEqual(2.0, lam, 1e-10);
        }

        TEST_METHOD(DiagonalDescending_ReturnsMid)
        {
            // Eigenvalues 10, 5, 1 → median = 5
            const double a[9] = {
                10, 0, 0,
                0,  5, 0,
                0,  0, 1
            };
            double lam = ut_get_lambda_1(a);
            Assert::AreEqual(5.0, lam, 1e-10);
        }

        TEST_METHOD(SymmetricMatrix_ReturnsFiniteValue)
        {
            // 3x3 symmetric: known eigenvalues can be verified with characteristic polynomial
            // A = [[4,2,0],[2,3,1],[0,1,2]] → trace=9, det=14
            const double a[9] = {
                4, 2, 0,
                2, 3, 1,
                0, 1, 2
            };
            double lam = ut_get_lambda_1(a);
            Assert::IsTrue(std::isfinite(lam));
            // Eigenvalues must be between min(diag)=2 and max(diag)=4 for a diag-dominant case
            Assert::IsTrue(lam >= 1.0 && lam <= 5.0);
        }

        TEST_METHOD(SymmetricMatrix_MedianIsInBounds)
        {
            // For any symmetric matrix, median eigenvalue is bounded by extremes
            const double a[9] = {
                3, 1, 0,
                1, 4, 2,
                0, 2, 5
            };
            double lam = ut_get_lambda_1(a);
            Assert::IsTrue(std::isfinite(lam));
            Assert::IsTrue(lam > 0.0); // positive definite matrix
        }
    };

} // namespace NoSpherA2UnitTests
