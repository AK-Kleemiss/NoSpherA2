#include "pch.h"
#include "../../NoSpherA2_DLL/unit_test_api.h"

#include <cmath>
#include <limits>
#include <string>

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

} // namespace NoSpherA2UnitTests
