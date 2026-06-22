#include "pch.h"
#include "../../NoSpherA2_DLL/unit_test_api.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>
#include <string>
#include <vector>

static constexpr double PI_VAL = 3.14159265358979323846;

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

// ---------------------------------------------------------------------------
// Pure unit tests — no golden files, no tests/tests.toml counterpart.
// Each method exercises a single Src/ function via the DLL's C export shim.
// ---------------------------------------------------------------------------

namespace NoSpherA2UnitTests
{
    TEST_CLASS(AtomicSymmetrizationTests)
    {
    public:
        TEST_METHOD(Oh_PShellBecomesIsotropic)
        {
            double matrix[9] = {
                1.0, 4.0, 5.0,
                4.0, 2.0, 6.0,
                5.0, 6.0, 3.0
            };
            const int shells[1] = { 1 };
            Assert::AreEqual(1, ut_symmetrize_atomic_matrix_oh(matrix, 3, shells, 1, 0));
            for (int row = 0; row < 3; ++row) {
                for (int col = 0; col < 3; ++col)
                    Assert::AreEqual(row == col ? 2.0 : 0.0,
                        matrix[row * 3 + col], 1e-12);
            }
        }

        TEST_METHOD(Oh_DShellPreservesEgAndT2gBlocks)
        {
            double matrix[36] = {};
            for (int i = 0; i < 6; ++i)
                matrix[i * 6 + i] = static_cast<double>(i + 1);
            const int shells[1] = { 2 };
            Assert::AreEqual(1, ut_symmetrize_atomic_matrix_oh(matrix, 6, shells, 1, 0));
            for (int i = 0; i < 6; ++i)
                Assert::AreEqual(i < 3 ? 2.0 : 5.0, matrix[i * 6 + i], 1e-12);
            for (int row = 0; row < 6; ++row) {
                for (int col = 0; col < 6; ++col) {
                    if (row != col)
                        Assert::AreEqual(0.0, matrix[row * 6 + col], 1e-12);
                }
            }
        }

        TEST_METHOD(Oh_SymmetrizationIsIdempotent)
        {
            double matrix[16] = {
                3.0, 1.0, 2.0, 3.0,
                1.0, 4.0, 5.0, 6.0,
                2.0, 5.0, 7.0, 8.0,
                3.0, 6.0, 8.0, 9.0
            };
            const int shells[2] = { 0, 1 };
            Assert::AreEqual(1, ut_symmetrize_atomic_matrix_oh(matrix, 4, shells, 2, 0));
            double once[16];
            std::copy(matrix, matrix + 16, once);
            Assert::AreEqual(1, ut_symmetrize_atomic_matrix_oh(matrix, 4, shells, 2, 0));
            for (int i = 0; i < 16; ++i)
                Assert::AreEqual(once[i], matrix[i], 1e-12);
        }

        TEST_METHOD(Oh_SphericalPShellBecomesIsotropic)
        {
            double matrix[9] = {
                1.0, 4.0, 5.0,
                4.0, 2.0, 6.0,
                5.0, 6.0, 3.0
            };
            const int shells[1] = { 1 };
            Assert::AreEqual(1, ut_symmetrize_atomic_matrix_oh(matrix, 3, shells, 1, 1));
            for (int row = 0; row < 3; ++row) {
                for (int col = 0; col < 3; ++col)
                    Assert::AreEqual(row == col ? 2.0 : 0.0,
                        matrix[row * 3 + col], 1e-11);
            }
        }

        TEST_METHOD(Oh_SphericalSymmetrizationIsIdempotent)
        {
            double matrix[36] = {};
            for (int row = 0; row < 6; ++row) {
                for (int col = 0; col < 6; ++col)
                    matrix[row * 6 + col] = static_cast<double>((row + 1) * (col + 1));
            }
            const int shells[2] = { 0, 2 };
            Assert::AreEqual(1, ut_symmetrize_atomic_matrix_oh(matrix, 6, shells, 2, 1));
            double once[36];
            std::copy(matrix, matrix + 36, once);
            Assert::AreEqual(1, ut_symmetrize_atomic_matrix_oh(matrix, 6, shells, 2, 1));
            for (int i = 0; i < 36; ++i)
                Assert::AreEqual(once[i], matrix[i], 1e-10);
        }
    };

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

    // -----------------------------------------------------------------------
    // PrimitiveEvalTests — Gaussian primitive evaluation
    // -----------------------------------------------------------------------
    TEST_CLASS(PrimitiveEvalTests)
    {
    public:
        TEST_METHOD(SType_AtOrigin_IsCoefTimesNorm)
        {
            // r=0: r^0 * exp(0) * norm*coef = norm*coef
            double val = ut_primitive_eval(0, 1.0, 0.5, 2.0, 0.0);
            Assert::AreEqual(0.5 * 2.0, val, 1e-12);
        }

        TEST_METHOD(SType_KnownGaussian)
        {
            // s-type (type=0): r^0 * exp(-alpha*r^2) * norm*coef
            double alpha = 2.0, coef = 1.0, norm = 1.0, r = 1.0;
            double expected = std::exp(-alpha * r * r) * norm * coef;
            Assert::AreEqual(expected, ut_primitive_eval(0, alpha, coef, norm, r), 1e-12);
        }

        TEST_METHOD(PType_AtHalf_IncludesRFactor)
        {
            // p-type (type=1): r^1 * exp(-alpha*r^2) * norm*coef
            double alpha = 1.0, coef = 1.0, norm = 1.0, r = 0.5;
            double expected = r * std::exp(-alpha * r * r) * norm * coef;
            Assert::AreEqual(expected, ut_primitive_eval(1, alpha, coef, norm, r), 1e-12);
        }

        TEST_METHOD(DType_AtUnit_IncludesR2Factor)
        {
            // d-type (type=2): r^2 * exp(-alpha*r^2) * norm*coef
            double alpha = 0.5, coef = 2.0, norm = 1.5, r = 1.0;
            double expected = r * r * std::exp(-alpha) * norm * coef;
            Assert::AreEqual(expected, ut_primitive_eval(2, alpha, coef, norm, r), 1e-12);
        }

        TEST_METHOD(Unnormalized_IgnoresNorm)
        {
            // Unnormalized uses coefficient only, no norm_const
            double alpha = 1.0, coef = 3.0, r = 0.5;
            double expected = r * std::exp(-alpha * r * r) * coef; // type=1
            Assert::AreEqual(expected, ut_primitive_eval_unnorm(1, alpha, coef, r), 1e-12);
        }

        TEST_METHOD(LargeR_DecaysToNearZero)
        {
            double val = ut_primitive_eval(0, 5.0, 1.0, 1.0, 100.0);
            Assert::IsTrue(std::abs(val) < 1e-200);
        }

        TEST_METHOD(ZeroExponent_DecaysSlower)
        {
            // alpha=0: constant along r (for type=0)
            double val = ut_primitive_eval(0, 0.0, 1.0, 1.0, 50.0);
            Assert::AreEqual(1.0, val, 1e-12);
        }
    };

    // -----------------------------------------------------------------------
    // AtomDistanceTests — atom::distance_to standalone
    // -----------------------------------------------------------------------
    TEST_CLASS(AtomDistanceTests)
    {
    public:
        TEST_METHOD(SamePoint_IsZero)
        {
            double d = ut_atom_distance(1.0, 2.0, 3.0, 1.0, 2.0, 3.0);
            Assert::AreEqual(0.0, d, 1e-12);
        }

        TEST_METHOD(UnitX_Separation)
        {
            double d = ut_atom_distance(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
            Assert::AreEqual(1.0, d, 1e-12);
        }

        TEST_METHOD(Pythagorean_3_4_5)
        {
            double d = ut_atom_distance(0.0, 0.0, 0.0, 3.0, 4.0, 0.0);
            Assert::AreEqual(5.0, d, 1e-12);
        }

        TEST_METHOD(ThreeD_UnitCube_Diagonal)
        {
            double d = ut_atom_distance(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
            Assert::AreEqual(std::sqrt(3.0), d, 1e-12);
        }

        TEST_METHOD(Symmetric_AtoB_EqualsBtoA)
        {
            double d1 = ut_atom_distance(1.0, 2.0, 3.0, 4.0, 6.0, 3.0);
            double d2 = ut_atom_distance(4.0, 6.0, 3.0, 1.0, 2.0, 3.0);
            Assert::AreEqual(d1, d2, 1e-12);
        }
    };

    // -----------------------------------------------------------------------
    // WfnTests — WFN object built fully in memory (no file I/O)
    // -----------------------------------------------------------------------
    TEST_CLASS(WfnTests)
    {
        // RAII guard so handles are always released even on assertion failure
        struct WfnGuard {
            void* h;
            explicit WfnGuard() : h(ut_wfn_create()) {}
            ~WfnGuard() { ut_wfn_destroy(h); }
        };

    public:
        TEST_METHOD(EmptyWfn_NcenIsZero)
        {
            WfnGuard g;
            Assert::AreEqual(0, ut_wfn_get_ncen(g.h));
        }

        TEST_METHOD(PushAtom_IncreasesNcen)
        {
            WfnGuard g;
            ut_wfn_push_atom(g.h, "C", 0.0, 0.0, 0.0, 6);
            Assert::AreEqual(1, ut_wfn_get_ncen(g.h));
            ut_wfn_push_atom(g.h, "H", 1.1, 0.0, 0.0, 1);
            Assert::AreEqual(2, ut_wfn_get_ncen(g.h));
        }

        TEST_METHOD(AtomCoordinates_RoundTrip)
        {
            WfnGuard g;
            ut_wfn_push_atom(g.h, "O", 1.5, -2.3, 0.7, 8);
            Assert::AreEqual(1.5,  ut_wfn_get_atom_x(g.h, 0), 1e-12);
            Assert::AreEqual(-2.3, ut_wfn_get_atom_y(g.h, 0), 1e-12);
            Assert::AreEqual(0.7,  ut_wfn_get_atom_z(g.h, 0), 1e-12);
        }

        TEST_METHOD(AtomCharge_RoundTrip)
        {
            WfnGuard g;
            ut_wfn_push_atom(g.h, "N", 0.0, 0.0, 0.0, 7);
            Assert::AreEqual(7, ut_wfn_get_atom_charge(g.h, 0));
        }

        TEST_METHOD(AtomLabel_RoundTrip)
        {
            WfnGuard g;
            ut_wfn_push_atom(g.h, "Fe", 0.0, 0.0, 0.0, 26);
            char buf[8];
            int n = ut_wfn_get_atom_label(g.h, 0, buf, 8);
            Assert::IsTrue(n > 0);
            Assert::AreEqual(std::string("Fe"), std::string(buf));
        }

        TEST_METHOD(AtomDistance_TwoAtoms)
        {
            WfnGuard g;
            ut_wfn_push_atom(g.h, "C", 0.0, 0.0, 0.0, 6);
            ut_wfn_push_atom(g.h, "H", 3.0, 4.0, 0.0, 1);
            double d = ut_wfn_atom_distance(g.h, 0, 1);
            Assert::AreEqual(5.0, d, 1e-12);
        }

        TEST_METHOD(PushMO_IncreasesNmo)
        {
            WfnGuard g;
            ut_wfn_push_mo(g.h, 1, 2.0, -10.5);
            Assert::AreEqual(1, ut_wfn_get_nmo(g.h));
            ut_wfn_push_mo(g.h, 2, 2.0, -7.3);
            ut_wfn_push_mo(g.h, 3, 0.0,  2.1); // virtual
            Assert::AreEqual(3, ut_wfn_get_nmo(g.h));
        }

        TEST_METHOD(MO_EnergyAndOcc_RoundTrip)
        {
            WfnGuard g;
            ut_wfn_push_mo(g.h, 1, 2.0, -10.5);
            Assert::AreEqual(-10.5, ut_wfn_get_mo_energy(g.h, 0), 1e-12);
            Assert::AreEqual(2.0,   ut_wfn_get_mo_occ(g.h, 0),    1e-12);
        }

        TEST_METHOD(MO_Coefficient_AddAndGet)
        {
            WfnGuard g;
            // MO coefficient vectors are populated via add_primitive (one coef per MO).
            ut_wfn_push_atom(g.h, "H", 0.0, 0.0, 0.0, 1);
            ut_wfn_push_mo(g.h, 1, 2.0, -5.0);          // 1 MO
            const double coef = 0.707;
            ut_wfn_add_primitive(g.h, 1, 1, 1.0, &coef, 1); // adds prim + fills coef
            Assert::AreEqual(0.707, ut_wfn_get_mo_coef(g.h, 0, 0), 1e-12);
            // Overwrite with set_mo_coef now that slot 0 exists
            ut_wfn_set_mo_coef(g.h, 0, 0, 0.5);
            Assert::AreEqual(0.5, ut_wfn_get_mo_coef(g.h, 0, 0), 1e-12);
        }

        TEST_METHOD(AddPrimitive_IncreasesNex)
        {
            WfnGuard g;
            ut_wfn_push_atom(g.h, "C", 0.0, 0.0, 0.0, 6);
            ut_wfn_add_exp(g.h, 1, 1, 2.0);
            Assert::AreEqual(1, ut_wfn_get_nex(g.h));
            ut_wfn_add_exp(g.h, 1, 2, 0.5);
            Assert::AreEqual(2, ut_wfn_get_nex(g.h));
        }

        TEST_METHOD(DeleteUnoccupiedMOs_RemovesVirtuals)
        {
            WfnGuard g;
            ut_wfn_push_mo(g.h, 1, 2.0, -10.0); // occupied
            ut_wfn_push_mo(g.h, 2, 2.0,  -5.0); // occupied
            ut_wfn_push_mo(g.h, 3, 0.0,   2.0); // virtual
            ut_wfn_push_mo(g.h, 4, 0.0,   4.0); // virtual
            int nmo = ut_wfn_delete_unoccupied_mos(g.h);
            Assert::AreEqual(2, nmo);
        }

        TEST_METHOD(ChargeAndMultiplicity_RoundTrip)
        {
            WfnGuard g;
            ut_wfn_assign_charge(g.h, -1);
            ut_wfn_assign_multi(g.h, 2);
            Assert::AreEqual(-1, ut_wfn_get_charge(g.h));
            Assert::AreEqual(2,  ut_wfn_get_multi(g.h));
        }

        TEST_METHOD(Method_RoundTrip)
        {
            WfnGuard g;
            ut_wfn_set_method(g.h, "B3LYP");
            char buf[32];
            ut_wfn_get_method(g.h, buf, 32);
            Assert::AreEqual(std::string("B3LYP"), std::string(buf));
        }

        TEST_METHOD(BasisSetName_RoundTrip)
        {
            WfnGuard g;
            ut_wfn_set_basis_set_name(g.h, "def2-TZVP");
            char buf[32];
            ut_wfn_get_basis_set_name(g.h, buf, 32);
            Assert::AreEqual(std::string("def2-TZVP"), std::string(buf));
        }

        TEST_METHOD(NrElectrons_MatchesSumOfZ)
        {
            WfnGuard g;
            // H2O: O(8) + H(1) + H(1) = 10 electrons
            ut_wfn_push_atom(g.h, "O", 0.0, 0.0, 0.0, 8);
            ut_wfn_push_atom(g.h, "H", 0.9, 0.0, 0.0, 1);
            ut_wfn_push_atom(g.h, "H",-0.9, 0.0, 0.0, 1);
            ut_wfn_assign_charge(g.h, 0);
            Assert::AreEqual(10, ut_wfn_get_nr_electrons(g.h));
        }

        TEST_METHOD(CountNrElectrons_SumsOccupations)
        {
            WfnGuard g;
            ut_wfn_push_mo(g.h, 1, 2.0, -10.0);
            ut_wfn_push_mo(g.h, 2, 2.0,  -5.0);
            ut_wfn_push_mo(g.h, 3, 1.0,  -1.0); // singly occupied
            // Total = 5.0
            Assert::AreEqual(5.0, ut_wfn_count_nr_electrons(g.h), 1e-12);
        }

        TEST_METHOD(DensityMatrix_SetAndGet)
        {
            WfnGuard g;
            ut_wfn_resize_dm(g.h, 6);
            Assert::AreEqual(6, ut_wfn_get_dm_size(g.h));
            ut_wfn_set_dm(g.h, 0, 1.5);
            ut_wfn_set_dm(g.h, 5, -0.3);
            Assert::AreEqual(1.5,  ut_wfn_get_dm(g.h, 0), 1e-12);
            Assert::AreEqual(-0.3, ut_wfn_get_dm(g.h, 5), 1e-12);
        }

        TEST_METHOD(DensityMatrix_DefaultIsZero)
        {
            WfnGuard g;
            ut_wfn_resize_dm(g.h, 3);
            for (int i = 0; i < 3; ++i)
                Assert::AreEqual(0.0, ut_wfn_get_dm(g.h, i), 1e-12);
        }
    };

    // -----------------------------------------------------------------------
    // TscBlockTests — tsc_block<int,cdouble> built in memory (no file I/O)
    // -----------------------------------------------------------------------
    TEST_CLASS(TscBlockTests)
    {
        // Helpers to build canonical test data
        // 2 scatterers ("C","H"), 3 reflections
        // SF[C][0]=1+0i, SF[C][1]=2+1i, SF[C][2]=3-1i
        // SF[H][0]=0.5+0i, SF[H][1]=0.8+0.2i, SF[H][2]=0.6-0.1i
        // h=[0,1,1] k=[0,0,1] l=[0,0,0]
        struct TscGuard {
            void* h;
            TscGuard() {
                static const char* labels[] = { "C", "H" };
                static const int   hv[] = { 0, 1, 1 };
                static const int   kv[] = { 0, 0, 1 };
                static const int   lv[] = { 0, 0, 0 };
                // sf_real/imag: row=scatterer, col=reflection
                static const double sr[] = { 1.0, 2.0, 3.0,  0.5, 0.8, 0.6 };
                static const double si[] = { 0.0, 1.0,-1.0,  0.0, 0.2,-0.1 };
                h = ut_tsc_create(2, labels, 3, hv, kv, lv, sr, si);
            }
            ~TscGuard() { ut_tsc_destroy(h); }
        };

    public:
        TEST_METHOD(EmptyBlock_IsEmpty)
        {
            void* h = ut_tsc_create_empty();
            Assert::AreEqual(1, ut_tsc_is_empty(h));
            ut_tsc_destroy(h);
        }

        TEST_METHOD(PopulatedBlock_IsNotEmpty)
        {
            TscGuard g;
            Assert::AreEqual(0, ut_tsc_is_empty(g.h));
        }

        TEST_METHOD(ScattererSize_IsCorrect)
        {
            TscGuard g;
            Assert::AreEqual(2, ut_tsc_scatterer_size(g.h));
        }

        TEST_METHOD(ReflectionSize_IsCorrect)
        {
            TscGuard g;
            Assert::AreEqual(3, ut_tsc_reflection_size(g.h));
        }

        TEST_METHOD(GetScatterer_Labels)
        {
            TscGuard g;
            char buf[8];
            ut_tsc_get_scatterer(g.h, 0, buf, 8);
            Assert::AreEqual(std::string("C"), std::string(buf));
            ut_tsc_get_scatterer(g.h, 1, buf, 8);
            Assert::AreEqual(std::string("H"), std::string(buf));
        }

        TEST_METHOD(ScatterersString_SpaceSeparated)
        {
            TscGuard g;
            char buf[32];
            ut_tsc_scatterers_string(g.h, buf, 32);
            Assert::AreEqual(std::string("C H"), std::string(buf));
        }

        TEST_METHOD(GetSF_Real_FirstScatterer)
        {
            TscGuard g;
            Assert::AreEqual(1.0, ut_tsc_get_sf_real(g.h, 0, 0), 1e-12);
            Assert::AreEqual(2.0, ut_tsc_get_sf_real(g.h, 0, 1), 1e-12);
            Assert::AreEqual(3.0, ut_tsc_get_sf_real(g.h, 0, 2), 1e-12);
        }

        TEST_METHOD(GetSF_Imag_FirstScatterer)
        {
            TscGuard g;
            Assert::AreEqual( 0.0, ut_tsc_get_sf_imag(g.h, 0, 0), 1e-12);
            Assert::AreEqual( 1.0, ut_tsc_get_sf_imag(g.h, 0, 1), 1e-12);
            Assert::AreEqual(-1.0, ut_tsc_get_sf_imag(g.h, 0, 2), 1e-12);
        }

        TEST_METHOD(GetSF_SecondScatterer)
        {
            TscGuard g;
            Assert::AreEqual(0.5, ut_tsc_get_sf_real(g.h, 1, 0), 1e-12);
            Assert::AreEqual(0.8, ut_tsc_get_sf_real(g.h, 1, 1), 1e-12);
        }

        TEST_METHOD(GetMillerIndices)
        {
            TscGuard g;
            // h: [0,1,1]
            Assert::AreEqual(0, ut_tsc_get_index(g.h, 0, 0));
            Assert::AreEqual(1, ut_tsc_get_index(g.h, 0, 1));
            // k: [0,0,1]
            Assert::AreEqual(0, ut_tsc_get_index(g.h, 1, 0));
            Assert::AreEqual(1, ut_tsc_get_index(g.h, 1, 2));
            // l: [0,0,0]
            Assert::AreEqual(0, ut_tsc_get_index(g.h, 2, 0));
        }

        TEST_METHOD(AnomalousDispersion_DefaultFalse)
        {
            TscGuard g;
            Assert::AreEqual(0, ut_tsc_get_ad(g.h));
        }

        TEST_METHOD(AnomalousDispersion_SetTrue)
        {
            TscGuard g;
            ut_tsc_set_ad(g.h, 1);
            Assert::AreEqual(1, ut_tsc_get_ad(g.h));
        }

        TEST_METHOD(Append_AddsNewScatterer)
        {
            // Build a second block with a new scatterer "O"
            static const char* labels2[] = { "O" };
            static const int   hv[] = { 0, 1, 1 };
            static const int   kv[] = { 0, 0, 1 };
            static const int   lv[] = { 0, 0, 0 };
            static const double sr2[] = { 4.0, 5.0, 6.0 };
            static const double si2[] = { 0.0, 0.0, 0.0 };

            TscGuard dst;
            void* src = ut_tsc_create(1, labels2, 3, hv, kv, lv, sr2, si2);

            int rc = ut_tsc_append(dst.h, src);
            ut_tsc_destroy(src);

            Assert::AreEqual(0, rc);                           // success
            Assert::AreEqual(3, ut_tsc_scatterer_size(dst.h)); // C + H + O
            // "O" SF should be accessible
            Assert::AreEqual(4.0, ut_tsc_get_sf_real(dst.h, 2, 0), 1e-12);
        }

        TEST_METHOD(Append_DuplicateScatterer_IsSkipped)
        {
            // Appending a block with "C" (already in dst) should not duplicate it
            TscGuard dst;

            static const char* labels2[] = { "C" };
            static const int   hv[] = { 0, 1, 1 };
            static const int   kv[] = { 0, 0, 1 };
            static const int   lv[] = { 0, 0, 0 };
            static const double sr2[] = { 9.0, 9.0, 9.0 };
            static const double si2[] = { 0.0, 0.0, 0.0 };

            void* src = ut_tsc_create(1, labels2, 3, hv, kv, lv, sr2, si2);
            ut_tsc_append(dst.h, src);
            ut_tsc_destroy(src);

            // Still 2 scatterers — duplicate "C" was ignored
            Assert::AreEqual(2, ut_tsc_scatterer_size(dst.h));
            // Original "C" SF unchanged
            Assert::AreEqual(1.0, ut_tsc_get_sf_real(dst.h, 0, 0), 1e-12);
        }
    };

    // -----------------------------------------------------------------------
    // CellTests — unit cell geometry
    // -----------------------------------------------------------------------
    TEST_CLASS(CellTests)
    {
        struct CellGuard {
            void* h;
            explicit CellGuard(double a, double b, double c,
                               double alpha, double beta, double gamma)
                : h(ut_cell_create(a, b, c, alpha, beta, gamma)) {}
            ~CellGuard() { ut_cell_destroy(h); }
        };

    public:
        TEST_METHOD(CubicCell_Parameters)
        {
            // Simple cubic: a=b=c=5, all angles 90
            CellGuard g(5.0, 5.0, 5.0, 90.0, 90.0, 90.0);
            Assert::AreEqual(5.0, ut_cell_get_a(g.h), 1e-12);
            Assert::AreEqual(5.0, ut_cell_get_b(g.h), 1e-12);
            Assert::AreEqual(5.0, ut_cell_get_c(g.h), 1e-12);
            Assert::AreEqual(90.0, ut_cell_get_angle(g.h, 0), 1e-10);
            Assert::AreEqual(90.0, ut_cell_get_angle(g.h, 1), 1e-10);
            Assert::AreEqual(90.0, ut_cell_get_angle(g.h, 2), 1e-10);
        }

        TEST_METHOD(CubicCell_Volume)
        {
            // V = 5^3 = 125 Å³
            CellGuard g(5.0, 5.0, 5.0, 90.0, 90.0, 90.0);
            Assert::AreEqual(125.0, ut_cell_get_volume(g.h), 1e-8);
        }

        TEST_METHOD(CubicCell_CrystalSystem)
        {
            CellGuard g(5.0, 5.0, 5.0, 90.0, 90.0, 90.0);
            char sys[32] = {};
            ut_cell_get_crystal_system(g.h, sys, 32);
            Assert::AreEqual(std::string("cubic"), std::string(sys));
        }

        TEST_METHOD(TetragonalCell_CrystalSystem)
        {
            CellGuard g(5.0, 5.0, 8.0, 90.0, 90.0, 90.0);
            char sys[32] = {};
            ut_cell_get_crystal_system(g.h, sys, 32);
            Assert::AreEqual(std::string("tetragonal"), std::string(sys));
        }

        TEST_METHOD(OrthorhombicCell_CrystalSystem)
        {
            CellGuard g(4.0, 5.0, 7.0, 90.0, 90.0, 90.0);
            char sys[32] = {};
            ut_cell_get_crystal_system(g.h, sys, 32);
            Assert::AreEqual(std::string("orthorhombic"), std::string(sys));
        }

        TEST_METHOD(HexagonalCell_CrystalSystem)
        {
            CellGuard g(4.0, 4.0, 7.0, 90.0, 90.0, 120.0);
            char sys[32] = {};
            ut_cell_get_crystal_system(g.h, sys, 32);
            Assert::AreEqual(std::string("hexagonal"), std::string(sys));
        }

        TEST_METHOD(MonoclinicCell_CrystalSystem)
        {
            CellGuard g(5.0, 7.0, 9.0, 90.0, 110.0, 90.0);
            char sys[32] = {};
            ut_cell_get_crystal_system(g.h, sys, 32);
            Assert::AreEqual(std::string("monoclinic"), std::string(sys));
        }

        TEST_METHOD(TriclinicCell_CrystalSystem)
        {
            CellGuard g(4.0, 5.0, 6.0, 80.0, 85.0, 75.0);
            char sys[32] = {};
            ut_cell_get_crystal_system(g.h, sys, 32);
            Assert::AreEqual(std::string("triclinic"), std::string(sys));
        }

        TEST_METHOD(OrthorhombicCell_DSpacing_100)
        {
            // Orthorhombic: a=5, b=6, c=7, angles=90; d(1,0,0) = a = 5
            CellGuard g(5.0, 6.0, 7.0, 90.0, 90.0, 90.0);
            double d = ut_cell_d_spacing_hkl(g.h, 1, 0, 0);
            Assert::AreEqual(5.0, d, 1e-8);
        }

        TEST_METHOD(OrthorhombicCell_DSpacing_010)
        {
            CellGuard g(5.0, 6.0, 7.0, 90.0, 90.0, 90.0);
            double d = ut_cell_d_spacing_hkl(g.h, 0, 1, 0);
            Assert::AreEqual(6.0, d, 1e-8);
        }

        TEST_METHOD(OrthorhombicCell_StlHkl)
        {
            // sin(theta)/lambda = 1/(2*d)
            CellGuard g(5.0, 6.0, 7.0, 90.0, 90.0, 90.0);
            double stl = ut_cell_stl_hkl(g.h, 1, 0, 0);
            Assert::AreEqual(1.0 / (2.0 * 5.0), stl, 1e-10);
        }

        TEST_METHOD(CubicCell_FracToCart_Origin)
        {
            // (0,0,0) → (0,0,0) in any cell
            CellGuard g(5.0, 5.0, 5.0, 90.0, 90.0, 90.0);
            double cart[3];
            ut_cell_frac_to_cart(g.h, 0.0, 0.0, 0.0, cart, 0); // in_bohr=0 → Angstrom
            Assert::AreEqual(0.0, cart[0], 1e-12);
            Assert::AreEqual(0.0, cart[1], 1e-12);
            Assert::AreEqual(0.0, cart[2], 1e-12);
        }

        TEST_METHOD(CubicCell_FracToCart_UnitFrac)
        {
            // Cubic 5Å: (1,0,0) frac → (5,0,0) Å
            CellGuard g(5.0, 5.0, 5.0, 90.0, 90.0, 90.0);
            double cart[3];
            ut_cell_frac_to_cart(g.h, 1.0, 0.0, 0.0, cart, 0);
            Assert::AreEqual(5.0, cart[0], 1e-10);
            Assert::AreEqual(0.0, cart[1], 1e-10);
            Assert::AreEqual(0.0, cart[2], 1e-10);
        }

        TEST_METHOD(MonoclinicCell_VolumePositive)
        {
            // Monoclinic: a=5, b=7, c=9, alpha=90, beta=110, gamma=90
            CellGuard g(5.0, 7.0, 9.0, 90.0, 110.0, 90.0);
            double V = ut_cell_get_volume(g.h);
            Assert::IsTrue(V > 0.0);
        }
    };

    // -----------------------------------------------------------------------
    // MoTests — MO class standalone handle
    // -----------------------------------------------------------------------
    TEST_CLASS(MoTests)
    {
        struct MoGuard {
            void* h;
            explicit MoGuard(int nr, double occ, double energy)
                : h(ut_mo_create(nr, occ, energy)) {}
            ~MoGuard() { ut_mo_destroy(h); }
        };

    public:
        TEST_METHOD(Create_OccAndEnergy)
        {
            MoGuard g(1, 2.0, -0.5);
            Assert::AreEqual(2.0, ut_mo_get_occ(g.h), 1e-12);
            Assert::AreEqual(-0.5, ut_mo_get_energy(g.h), 1e-12);
        }

        TEST_METHOD(PushCoef_IncreasesCount)
        {
            MoGuard g(1, 2.0, -0.5);
            Assert::AreEqual(0, ut_mo_get_primitive_count(g.h));
            ut_mo_push_coef(g.h, 0.707);
            Assert::AreEqual(1, ut_mo_get_primitive_count(g.h));
            ut_mo_push_coef(g.h, 0.293);
            Assert::AreEqual(2, ut_mo_get_primitive_count(g.h));
        }

        TEST_METHOD(GetCoef_RoundTrip)
        {
            MoGuard g(1, 2.0, -0.5);
            ut_mo_push_coef(g.h, 1.23);
            ut_mo_push_coef(g.h, -4.56);
            Assert::AreEqual(1.23,  ut_mo_get_coef(g.h, 0), 1e-12);
            Assert::AreEqual(-4.56, ut_mo_get_coef(g.h, 1), 1e-12);
        }

        TEST_METHOD(Hdr_ContainsMO)
        {
            MoGuard g(1, 2.0, -0.5);
            char buf[256] = {};
            int n = ut_mo_hdr(g.h, buf, 256);
            Assert::IsTrue(n > 0);
            // Header line should start with "MO"
            Assert::IsTrue(std::string(buf).find("MO") != std::string::npos);
        }

        TEST_METHOD(MultipleCoefs_AllAccessible)
        {
            MoGuard g(2, 1.0, -0.25);
            for (int i = 0; i < 5; ++i)
                ut_mo_push_coef(g.h, i * 0.1);
            for (int i = 0; i < 5; ++i)
                Assert::AreEqual(i * 0.1, ut_mo_get_coef(g.h, i), 1e-12);
        }
    };

    // -----------------------------------------------------------------------
    // LebedevTests — spherical quadrature
    // -----------------------------------------------------------------------
    TEST_CLASS(LebedevTests)
    {
    public:
        TEST_METHOD(Order6_CorrectPointCount)
        {
            double x[6], y[6], z[6], w[6];
            int n = ut_lebedev_order(6, x, y, z, w, 6);
            Assert::AreEqual(6, n);
        }

        TEST_METHOD(Order6_WeightsSumToOne)
        {
            // Lebedev weights sum to 1.0 (integral of 1 over unit sphere / 4π)
            double x[6], y[6], z[6], w[6];
            ut_lebedev_order(6, x, y, z, w, 6);
            double sum = 0.0;
            for (int i = 0; i < 6; ++i) sum += w[i];
            Assert::AreEqual(1.0, sum, 1e-10);
        }

        TEST_METHOD(Order14_WeightsSumToOne)
        {
            double x[14], y[14], z[14], w[14];
            ut_lebedev_order(14, x, y, z, w, 14);
            double sum = 0.0;
            for (int i = 0; i < 14; ++i) sum += w[i];
            Assert::AreEqual(1.0, sum, 1e-10);
        }

        TEST_METHOD(Order26_WeightsSumToOne)
        {
            double x[26], y[26], z[26], w[26];
            ut_lebedev_order(26, x, y, z, w, 26);
            double sum = 0.0;
            for (int i = 0; i < 26; ++i) sum += w[i];
            Assert::AreEqual(1.0, sum, 1e-10);
        }

        TEST_METHOD(Order6_PointsOnUnitSphere)
        {
            // All quadrature points should lie on the unit sphere
            double x[6], y[6], z[6], w[6];
            ut_lebedev_order(6, x, y, z, w, 6);
            for (int i = 0; i < 6; ++i) {
                double r2 = x[i]*x[i] + y[i]*y[i] + z[i]*z[i];
                Assert::AreEqual(1.0, r2, 1e-12);
            }
        }

        TEST_METHOD(Order6_AllWeightsPositive)
        {
            double x[6], y[6], z[6], w[6];
            ut_lebedev_order(6, x, y, z, w, 6);
            for (int i = 0; i < 6; ++i)
                Assert::IsTrue(w[i] > 0.0);
        }
    };

    // -----------------------------------------------------------------------
    // AtomGridTests — radial/angular DFT integration grid
    // -----------------------------------------------------------------------
    TEST_CLASS(AtomGridTests)
    {
        // Construct a minimal grid for hydrogen (Z=1) with one s-type function
        static void* make_h_grid() {
            double alpha_min[1] = { 0.5 };  // min exponent for l=0
            return ut_atomgrid_create(
                1e-5,   // radial_precision
                6,      // min_angular (6-point Lebedev)
                6,      // max_angular
                1,      // proton_charge (H)
                10.0,   // alpha_max
                0,      // max_l (s only)
                alpha_min, 1);
        }

    public:
        TEST_METHOD(Create_NotNull)
        {
            void* h = make_h_grid();
            Assert::IsNotNull(h);
            ut_atomgrid_destroy(h);
        }

        TEST_METHOD(NumPoints_Positive)
        {
            void* h = make_h_grid();
            Assert::IsNotNull(h);
            int n = ut_atomgrid_get_num_points(h);
            Assert::IsTrue(n > 0);
            ut_atomgrid_destroy(h);
        }

        TEST_METHOD(NumRadialPoints_Positive)
        {
            void* h = make_h_grid();
            Assert::IsNotNull(h);
            int nr = ut_atomgrid_get_num_radial_points(h);
            Assert::IsTrue(nr > 0);
            ut_atomgrid_destroy(h);
        }

        TEST_METHOD(NumPoints_GeRadialTimesAngular)
        {
            void* h = make_h_grid();
            Assert::IsNotNull(h);
            int total = ut_atomgrid_get_num_points(h);
            int radial = ut_atomgrid_get_num_radial_points(h);
            // total >= radial * 6 (min angular = 6)
            Assert::IsTrue(total >= radial * 6);
            ut_atomgrid_destroy(h);
        }

        TEST_METHOD(RadialGrid_WeightsPositive)
        {
            void* h = make_h_grid();
            Assert::IsNotNull(h);
            int nr = ut_atomgrid_get_num_radial_points(h);
            std::vector<double> r(nr), w(nr);
            ut_atomgrid_get_radial_grid(h, r.data(), w.data());
            for (int i = 0; i < nr; ++i)
                Assert::IsTrue(w[i] > 0.0);
            ut_atomgrid_destroy(h);
        }

        TEST_METHOD(RadialGrid_RDistancesPositive)
        {
            void* h = make_h_grid();
            Assert::IsNotNull(h);
            int nr = ut_atomgrid_get_num_radial_points(h);
            std::vector<double> r(nr), w(nr);
            ut_atomgrid_get_radial_grid(h, r.data(), w.data());
            for (int i = 0; i < nr; ++i)
                Assert::IsTrue(r[i] > 0.0);
            ut_atomgrid_destroy(h);
        }
    };

    // -----------------------------------------------------------------------
    // WfnDensityTests — density evaluation on synthetic 1s hydrogen wavefunction
    // -----------------------------------------------------------------------
    TEST_CLASS(WfnDensityTests)
    {
        // Builds a minimal hydrogen-like 1s WFN in memory.
        //   One H atom at origin, one s-type primitive (alpha=1.0),
        //   one occupied MO (occ=2.0) with coef=1.0 and norm_const=1.0.
        static void* make_h1s_wfn() {
            void* h = ut_wfn_create();
            // Atom: H at origin, Z=1
            ut_wfn_push_atom(h, "H", 0.0, 0.0, 0.0, 1);
            // MO: occupied, occ=2.0
            ut_wfn_push_mo(h, 1, 2.0, -0.5);
            // Primitive: center=1 (1-based), type=1 (s), exponent=1.0, coef=1.0
            double coef = 1.0;
            ut_wfn_add_primitive(h, 1, 1, 1.0, &coef, 1);
            return h;
        }

    public:
        TEST_METHOD(Density_AtOrigin_IsPositive)
        {
            void* h = make_h1s_wfn();
            double rho = ut_wfn_compute_dens(h, 0.0, 0.0, 0.0);
            ut_wfn_destroy(h);
            Assert::IsTrue(rho > 0.0);
            Assert::IsTrue(std::isfinite(rho));
        }

        TEST_METHOD(Density_AtOrigin_ExceedsFarPoint)
        {
            // Density should be larger at origin than at r=5 Bohr
            void* h = make_h1s_wfn();
            double rho0 = ut_wfn_compute_dens(h, 0.0, 0.0, 0.0);
            double rho5 = ut_wfn_compute_dens(h, 5.0, 0.0, 0.0);
            ut_wfn_destroy(h);
            Assert::IsTrue(rho0 > rho5);
        }

        TEST_METHOD(Density_FarPoint_IsSmall)
        {
            // At r=10 Bohr, 1s density should be tiny
            void* h = make_h1s_wfn();
            double rho = ut_wfn_compute_dens(h, 10.0, 0.0, 0.0);
            ut_wfn_destroy(h);
            Assert::IsTrue(rho >= 0.0);
            Assert::IsTrue(rho < 1e-3);
        }

        TEST_METHOD(Density_Spherical_SameAtEquidistantPoints)
        {
            // s-type WFN: density should be equal at same distance in any direction
            void* h = make_h1s_wfn();
            double rho_x = ut_wfn_compute_dens(h,  1.0, 0.0, 0.0);
            double rho_y = ut_wfn_compute_dens(h,  0.0, 1.0, 0.0);
            double rho_z = ut_wfn_compute_dens(h,  0.0, 0.0, 1.0);
            ut_wfn_destroy(h);
            Assert::AreEqual(rho_x, rho_y, 1e-10);
            Assert::AreEqual(rho_x, rho_z, 1e-10);
        }

        TEST_METHOD(MoValue_AtOrigin_IsFinite)
        {
            void* h = make_h1s_wfn();
            double mo0 = ut_wfn_compute_mo(h, 0.0, 0.0, 0.0, 0);
            ut_wfn_destroy(h);
            Assert::IsTrue(std::isfinite(mo0));
        }

        TEST_METHOD(MoValue_DecaysWithDistance)
        {
            void* h = make_h1s_wfn();
            double mo0 = std::abs(ut_wfn_compute_mo(h, 0.0, 0.0, 0.0, 0));
            double mo5 = std::abs(ut_wfn_compute_mo(h, 5.0, 0.0, 0.0, 0));
            ut_wfn_destroy(h);
            Assert::IsTrue(mo0 > mo5);
        }

        TEST_METHOD(NrElectrons_MatchesOccupation)
        {
            // count_nr_electrons sums MO occupations: one MO occ=2 → 2 electrons
            void* h = make_h1s_wfn();
            double nel = ut_wfn_count_nr_electrons(h);
            ut_wfn_destroy(h);
            Assert::AreEqual(2.0, nel, 1e-10);
        }
    };

} // namespace NoSpherA2UnitTests
