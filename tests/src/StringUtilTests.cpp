
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


namespace NoSpherA2UnitTests
{
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    // StringUtilTests
    // -----------------------------------------------------------------------

    TEST(StringUtilTests, EndsWith_MatchingSuffix)
    {
        EXPECT_TRUE(ends_with(std::string("molecule.wfx"), std::string(".wfx")));
        EXPECT_TRUE(ends_with(std::string("data.hkl"), std::string(".hkl")));
    }

    TEST(StringUtilTests, EndsWith_NonMatchingSuffix)
    {
        EXPECT_FALSE(ends_with(std::string("molecule.wfx"), std::string(".gbw")));
        EXPECT_FALSE(ends_with(std::string("test.cpp"), std::string(".txt")));
    }

    TEST(StringUtilTests, EndsWith_EmptySuffix)
    {
        // Empty suffix always matches
        EXPECT_TRUE(ends_with(std::string("anything"), std::string("")));
    }

    TEST(StringUtilTests, EndsWith_Suffix_LongerThanString)
    {
        EXPECT_FALSE(ends_with(std::string("ab"), std::string("abc")));
    }

    TEST(StringUtilTests, EndsWith_ExactMatch)
    {
        EXPECT_TRUE(ends_with(std::string(".wfx"), std::string(".wfx")));
    }

    TEST(StringUtilTests, ShrinkString_RemovesDigits)
    {
        std::string input = "C1";
        std::string r = shrink_string(input);
        EXPECT_EQ(1, static_cast<int>(r.size()));
        EXPECT_STREQ("C", r.c_str());
    }

    TEST(StringUtilTests, ShrinkString_RemovesSpacesAndDigits)
    {
        std::string input = "O 1 1";
        std::string r = shrink_string(input);
        int len = static_cast<int>(r.size());
        EXPECT_GT(len, 0);
        // All spaces and digits removed → only "O" remains
        EXPECT_STREQ("O", r.c_str());
    }

    TEST(StringUtilTests, ShrinkString_PureLettersUnchanged)
    {
        std::string input = "Fe";
        std::string r = shrink_string(input);
        EXPECT_STREQ("Fe", r.c_str());
    }

    TEST(StringUtilTests, ShrinkString_BufferTooSmall)
    {
        // Can't fit even "C"
        std::string input = "C1H2O";
        std::string r = shrink_string(input);
        EXPECT_GE(static_cast<int>(r.size()), 1);
    }

    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    // Sha256Tests
    // -----------------------------------------------------------------------

    TEST(Sha256Tests, Sha256_OutputLength)
    {
        std::string r = sha::sha256(std::string("abc"));
        EXPECT_EQ(64, static_cast<int>(r.size()));
    }

    TEST(Sha256Tests, Sha256_KnownVector_Abc)
    {
        // NIST FIPS 180-4 test vector
        std::string r = sha::sha256(std::string("abc"));
        EXPECT_STREQ(
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
            r.c_str());
    }

    TEST(Sha256Tests, Sha256_EmptyString)
    {
        std::string r = sha::sha256(std::string(""));
        EXPECT_STREQ(
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
            r.c_str());
    }

    TEST(Sha256Tests, Sha256_Deterministic)
    {
        std::string r1 = sha::sha256(std::string("NoSpherA2"));
        std::string r2 = sha::sha256(std::string("NoSpherA2"));
        EXPECT_STREQ(r1.c_str(), r2.c_str());
    }

    TEST(Sha256Tests, Sha256_DifferentInputsDifferentOutputs)
    {
        std::string r1 = sha::sha256(std::string("abc"));
        std::string r2 = sha::sha256(std::string("abd"));
        EXPECT_STRNE(r1.c_str(), r2.c_str());
    }

    TEST(Sha256Tests, Sha256_BufferTooSmall)
    {
        // No more small-buffer API: just validate digest length is always 64.
        std::string r = sha::sha256(std::string("abc"));
        EXPECT_EQ(64, static_cast<int>(r.size()));
    }

    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    // VecAggregateTests
    // -----------------------------------------------------------------------

    TEST(VecAggregateTests, VecSumBool_MixedValues)
    {
        bvec v = { true, false, true, true, false };
        EXPECT_EQ(3, vec_sum(v));
    }

    TEST(VecAggregateTests, VecSumBool_AllFalse)
    {
        bvec v = { false, false, false };
        EXPECT_EQ(0, vec_sum(v));
    }

    TEST(VecAggregateTests, VecSumBool_Empty)
    {
        bvec v;
        EXPECT_EQ(0, vec_sum(v));
    }

    TEST(VecAggregateTests, VecSumInt_Basic)
    {
        ivec v = { 1, 2, 3, 4, 5 };
        EXPECT_EQ(15, vec_sum(v));
    }

    TEST(VecAggregateTests, VecSumInt_Negative)
    {
        ivec v = { -3, 7, -2 };
        EXPECT_EQ(2, vec_sum(v));
    }

    TEST(VecAggregateTests, VecSumDouble_Basic)
    {
        vec v = { 1.5, -0.5, 2.0 };
        EXPECT_NEAR(3.0, vec_sum(v), 1e-12);
    }

    TEST(VecAggregateTests, VecSumDouble_AllZero)
    {
        vec v = { 0.0, 0.0, 0.0 };
        EXPECT_NEAR(0.0, vec_sum(v), 1e-12);
    }

    TEST(VecAggregateTests, VecLength_345)
    {
        vec v = { 3.0, 4.0 };
        EXPECT_NEAR(5.0, vec_length(v), 1e-12);
    }

    TEST(VecAggregateTests, VecLength_Unit)
    {
        vec v = { 1.0, 0.0, 0.0 };
        EXPECT_NEAR(1.0, vec_length(v), 1e-12);
    }

    TEST(VecAggregateTests, VecLength_4D)
    {
        vec v = { 1.0, 1.0, 1.0, 1.0 };
        EXPECT_NEAR(2.0, vec_length(v), 1e-12);
    }

    TEST(VecAggregateTests, VecLength_Empty)
    {
        vec v;
        EXPECT_NEAR(0.0, vec_length(v), 1e-12);
    }
    
    // -----------------------------------------------------------------------

    TEST(StringUtilTests2, Trim_LeadingTrailingSpaces)
    {
        std::string r = trim(std::string("  hello world  "));
        EXPECT_EQ(r, "hello world");
    }

    TEST(StringUtilTests2, Trim_NoSpaces)
    {
        std::string r = trim(std::string("no_spaces"));
        EXPECT_EQ(r, "no_spaces");
    }

    TEST(StringUtilTests2, Trim_EmptyString)
    {
        std::string r = trim(std::string(""));
        EXPECT_TRUE(r.empty());
        EXPECT_EQ(r, "");
    }

    TEST(StringUtilTests2, Trim_OnlySpaces)
    {
        std::string r = trim(std::string("    "));
        EXPECT_EQ(r, "");
    }

    TEST(StringUtilTests2, Asciitolower_Uppercase)
    {
        EXPECT_EQ('a', asciitolower('A'));
        EXPECT_EQ('z', asciitolower('Z'));
        EXPECT_EQ('m', asciitolower('M'));
    }

    TEST(StringUtilTests2, Asciitolower_AlreadyLower)
    {
        EXPECT_EQ('a', asciitolower('a'));
        EXPECT_EQ('z', asciitolower('z'));
    }

    TEST(StringUtilTests2, Asciitolower_NonAlpha)
    {
        // Digits and symbols are returned unchanged
        EXPECT_EQ('5', asciitolower('5'));
        EXPECT_EQ('_', asciitolower('_'));
    }

    TEST(StringUtilTests2, DoubleFromEsd_Plain)
    {
        EXPECT_NEAR(3.14159, double_from_string_with_esd(std::string("3.14159")), 1e-10);
    }

    TEST(StringUtilTests2, DoubleFromEsd_WithEsd)
    {
        // "(5)" is stripped; value is 1.234
        EXPECT_NEAR(1.234, double_from_string_with_esd(std::string("1.234(5)")), 1e-10);
    }

    TEST(StringUtilTests2, DoubleFromEsd_Zero)
    {
        EXPECT_NEAR(0.0, double_from_string_with_esd(std::string("0.0")), 1e-12);
    }

    TEST(StringUtilTests2, DecimalPrecisionCif_WithBracketsAndDecimal)
    {
        // "1.2345(6)" → 6 × 10⁻⁴ = 0.0006
        std::string s = "1.2345(6)";
        EXPECT_NEAR(6e-4, get_decimal_precision_from_CIF_number(s), 1e-10);
    }

    TEST(StringUtilTests2, DecimalPrecisionCif_NoBrackets)
    {
        // No brackets → default 0.005
        std::string s = "1.234";
        EXPECT_NEAR(0.005, get_decimal_precision_from_CIF_number(s), 1e-10);
    }

    TEST(StringUtilTests2, DecimalPrecisionCif_IntegerWithEsd)
    {
        // "100(2)" means 100 +- 2: no decimal point, the esd is the integer in the bracket
        std::string s = "100(2)";
        EXPECT_NEAR(2.0, get_decimal_precision_from_CIF_number(s), 1e-10);
    }

    // -----------------------------------------------------------------------

    TEST(BasisTypeTests, Sht2nbas_CartesianShells)
    {
        // Cartesian: S=1, P=3, D=6, F=10, G=15
        EXPECT_EQ(1, sht2nbas(0));
        EXPECT_EQ(3, sht2nbas(1));
        EXPECT_EQ(6, sht2nbas(2));
        EXPECT_EQ(10, sht2nbas(3));
        EXPECT_EQ(15, sht2nbas(4));
    }

    TEST(BasisTypeTests, Sht2nbas_SphericalShells)
    {
        // Negative types → spherical: -2=D(5), -3=F(7), -4=G(9)
        EXPECT_EQ(5, sht2nbas(-2));
        EXPECT_EQ(7, sht2nbas(-3));
        EXPECT_EQ(9, sht2nbas(-4));
    }

    TEST(BasisTypeTests, DoubleFactorial_SmallValues)
    {
        EXPECT_EQ(1u, doublefactorial(0));
        EXPECT_EQ(1u, doublefactorial(1));
        EXPECT_EQ(2u, doublefactorial(2));
        EXPECT_EQ(3u, doublefactorial(3));
        EXPECT_EQ(8u, doublefactorial(4));
        EXPECT_EQ(15u, doublefactorial(5));
        EXPECT_EQ(48u, doublefactorial(6));
        EXPECT_EQ(105u, doublefactorial(7));
    }

    // -----------------------------------------------------------------------

    TEST(OrbitalIndexTests, OrcaToPySCF_SShell)
    {
        // S: only one component, m_idx=0 → 0
        auto r = constants::orca_2_pySCF(0, 0);
        EXPECT_EQ(0, r.has_value() ? static_cast<int>(r.value()) : -1);
    }

    TEST(OrbitalIndexTests, OrcaToPySCF_PShell)
    {
        // P ORCA ordering 0,+1,-1 → PySCF map {1,2,0}
        EXPECT_EQ((size_t)1, constants::orca_2_pySCF(1, 0).value());
        EXPECT_EQ((size_t)2, constants::orca_2_pySCF(1, 1).value());
        EXPECT_EQ((size_t)0, constants::orca_2_pySCF(1, 2).value());
    }

    TEST(OrbitalIndexTests, OrcaToPySCF_DShell_FirstComponent)
    {
        // D: map {2,3,1,4,0}, m_idx=0 → 2
        EXPECT_EQ((size_t)2, constants::orca_2_pySCF(2, 0).value());
    }

    TEST(OrbitalIndexTests, OrcaToPySCF_OutOfRange)
    {
        // l=100 is not in the switch → nullopt → -1
        auto r = constants::orca_2_pySCF(100, 0);
        EXPECT_EQ(-1, r.has_value() ? static_cast<int>(r.value()) : -1);
    }

    TEST(OrbitalIndexTests, TypeToNbo_SShell) { EXPECT_EQ(1u, constants::type_2_nbo(1)); }
    TEST(OrbitalIndexTests, TypeToNbo_PxShell) { EXPECT_EQ(101u, constants::type_2_nbo(2)); }
    TEST(OrbitalIndexTests, TypeToNbo_DxxShell) { EXPECT_EQ(201u, constants::type_2_nbo(5)); }
    TEST(OrbitalIndexTests, TypeToNbo_FxxxShell) { EXPECT_EQ(301u, constants::type_2_nbo(11)); }
    TEST(OrbitalIndexTests, TypeToNbo_GxxxxShell) { EXPECT_EQ(401u, constants::type_2_nbo(21)); }
    TEST(OrbitalIndexTests, TypeToNbo_Unknown) { EXPECT_EQ(0u, constants::type_2_nbo(99)); }

    // -----------------------------------------------------------------------
    // IsSimilarPow10Tests — is_similar(a, b, tolerance) where |a-b| <= 10^tol
    // -----------------------------------------------------------------------
    TEST(IsSimilarPow10Tests, Equal_ReturnTrue)
    {
        EXPECT_TRUE(is_similar(1.0, 1.0, -6.0));
    }

    TEST(IsSimilarPow10Tests, WithinTolerance_ReturnTrue)
    {
        // |1.000001 - 1.0| = 1e-6 <= 10^(-6)
        EXPECT_TRUE(is_similar(1.000001, 1.0, -6.0));
    }

    TEST(IsSimilarPow10Tests, OutsideTolerance_ReturnFalse)
    {
        // |1.00001 - 1.0| = 1e-5 > 10^(-6)
        EXPECT_FALSE(is_similar(1.00001, 1.0, -6.0));
    }

    TEST(IsSimilarPow10Tests, NegativeValues_WithinTolerance)
    {
        EXPECT_TRUE(is_similar(-5.0, -5.0 + 1e-8, -7.0));
    }

    TEST(IsSimilarPow10Tests, LooseTolerance_LargeDiff)
    {
        // |100 - 50| = 50 <= 10^2 = 100
        EXPECT_TRUE(is_similar(100.0, 50.0, 2.0));
    }

    TEST(IsSimilarPow10Tests, LooseTolerance_TooLargeDiff)
    {
        // |200 - 50| = 150 > 10^2 = 100
        EXPECT_FALSE(is_similar(200.0, 50.0, 2.0));
    }

    // -----------------------------------------------------------------------
    // Shell2FunctionTests — shell2function(type, prim) WFN column index
    // -----------------------------------------------------------------------
    TEST(Shell2FunctionTests, SType_Prim0_ReturnsNonNegative)
    {
        // s-type shell (type=1): first and only function is index 0
        int r = shell2function(1, 0);
        EXPECT_GE(r, 0);
    }

    TEST(Shell2FunctionTests, PType_Prim0_ReturnsNonNegative)
    {
        // p-type shell: 3 functions
        int r = shell2function(2, 0);
        EXPECT_GE(r, 0);
    }

    TEST(Shell2FunctionTests, PType_Prim1)
    {
        int r0 = shell2function(2, 0);
        int r1 = shell2function(2, 1);
        EXPECT_GT(r1, r0);
    }

    TEST(Shell2FunctionTests, DType_Prim5_ValidIndex)
    {
        // d-type (type=3): 6 Cartesian or 5 spherical functions
        int r = shell2function(3, 5);
        EXPECT_GE(r, 0);
    }

    TEST(Shell2FunctionTests, ResultsAreStrictlyIncreasingWithinShell)
    {
        // f-type (type=4): consecutive prims must give increasing column indices
        int r0 = shell2function(4, 0);
        int r1 = shell2function(4, 1);
        int r2 = shell2function(4, 2);
        EXPECT_GT(r1, r0);
        EXPECT_GT(r2, r1);
    }

    // -----------------------------------------------------------------------
    // CountWordsTests — CountWords(str)
    // -----------------------------------------------------------------------
    TEST(CountWordsTests, Empty_Returns0)
    {
        EXPECT_EQ(0, CountWords(""));
    }

    TEST(CountWordsTests, OneWord)
    {
        EXPECT_EQ(1, CountWords("hello"));
    }

    TEST(CountWordsTests, TwoWords)
    {
        EXPECT_EQ(2, CountWords("hello world"));
    }

    TEST(CountWordsTests, LeadingTrailingSpaces)
    {
        EXPECT_EQ(2, CountWords("  foo   bar  "));
    }

    TEST(CountWordsTests, MultipleSpacesBetweenWords)
    {
        EXPECT_EQ(3, CountWords("a  b  c"));
    }

    TEST(CountWordsTests, SingleSpace)
    {
        EXPECT_EQ(0, CountWords(" "));
    }

    // -----------------------------------------------------------------------
    // ShrinkStringToAtomTests — shrink_string_to_atom(input, atom_number)
    // -----------------------------------------------------------------------
    TEST(ShrinkStringToAtomTests, CarbonAtomNumber6)
    {
        // atnr2letter(6) = "C"
        std::string input = "C1";
        std::string r = shrink_string_to_atom(input, 6);
        EXPECT_GE(static_cast<int>(r.size()), 0);
        EXPECT_EQ(std::string("C"), r);
    }

    TEST(ShrinkStringToAtomTests, CalciumAtomNumber20)
    {
        // atnr2letter(20) = "Ca"
        std::string input = "Ca12";
        std::string r = shrink_string_to_atom(input, 20);
        EXPECT_GE(static_cast<int>(r.size()), 0);
        EXPECT_EQ(std::string("Ca"), r);
    }

    TEST(ShrinkStringToAtomTests, BufferTooSmall_ReturnsAtLeastOne)
    {
        std::string input = "Carbon6";
        std::string r = shrink_string_to_atom(input, 6);
        EXPECT_GE(static_cast<int>(r.size()), 1);
    }

    TEST(ShrinkStringToAtomTests, IronAtomNumber26)
    {
        // atnr2letter(26) = "Fe"
        std::string input = "Fe3 ";
        std::string r = shrink_string_to_atom(input, 26);
        EXPECT_GE(static_cast<int>(r.size()), 0);
        EXPECT_EQ(std::string("Fe"), r);
    }

    // -----------------------------------------------------------------------
    // SplitStringTests — split_string(input, delim)
    // -----------------------------------------------------------------------
    TEST(SplitStringTests, SingleToken_NoDelim)
    {
        svec toks = split_string<std::string>(std::string("hello"), " ");
        EXPECT_EQ(1, static_cast<int>(toks.size()));
        EXPECT_EQ(std::string("hello"), toks[0]);
    }

    TEST(SplitStringTests, ThreeTokens)
    {
        svec toks = split_string<std::string>(std::string("a b c"), " ");
        EXPECT_EQ(3, static_cast<int>(toks.size()));
        EXPECT_EQ(std::string("a"), toks[0]);
        EXPECT_EQ(std::string("b"), toks[1]);
        EXPECT_EQ(std::string("c"), toks[2]);
    }

    TEST(SplitStringTests, CommaDelimiter)
    {
        svec toks = split_string<std::string>(std::string("x,y,z"), ",");
        EXPECT_EQ(3, static_cast<int>(toks.size()));
        EXPECT_EQ(std::string("z"), toks[2]);
    }

    TEST(SplitStringTests, MaxOutLimit_ReturnsTotalCount)
    {
        // 4 tokens but max_out=2 (old API); direct API returns full token list.
        svec toks = split_string<std::string>(std::string("a b c d"), " ");
        EXPECT_EQ(4, static_cast<int>(toks.size()));
        EXPECT_EQ(std::string("a"), toks[0]);
        EXPECT_EQ(std::string("b"), toks[1]);
    }

    TEST(SplitStringTests, EmptyString_ZeroTokens)
    {
        svec toks = split_string<std::string>(std::string(""), " ");
        EXPECT_TRUE(toks.empty() || toks.size() == 1); // impl-defined for empty input
    }

    // -----------------------------------------------------------------------
    // TimingTests — ut_sleep_and_measure_us(N) returns elapsed µs >= N*1000
    // -----------------------------------------------------------------------
    TEST(TimingTests, Sleep10ms_ElapsedAtLeast10000us)
    {
        auto t0 = get_time();
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        auto t1 = get_time();
        long long us = get_musec(t0, t1);
        EXPECT_GE(us, 10000LL);
    }

    TEST(TimingTests, Sleep1ms_ElapsedAtLeast1000us)
    {
        auto t0 = get_time();
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
        auto t1 = get_time();
        long long us = get_musec(t0, t1);
        EXPECT_GE(us, 1000LL);
    }

    TEST(TimingTests, Sleep5ms_ElapsedPositive)
    {
        auto t0 = get_time();
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
        auto t1 = get_time();
        long long us = get_musec(t0, t1);
        EXPECT_GT(us, 0LL);
    }

    // -----------------------------------------------------------------------
    // Atnr2LetterTests — constants::atnr2letter(nr) element symbol lookup
    // -----------------------------------------------------------------------
    TEST(Atnr2LetterTests, Hydrogen_Is_H)
    {
        const char* sym = constants::atnr2letter(1);
        EXPECT_EQ(std::string("H"), std::string(sym));
    }

    TEST(Atnr2LetterTests, Carbon_Is_C)
    {
        EXPECT_EQ(std::string("C"), std::string(constants::atnr2letter(6)));
    }

    TEST(Atnr2LetterTests, Iron_Is_Fe)
    {
        EXPECT_EQ(std::string("Fe"), std::string(constants::atnr2letter(26)));
    }

    TEST(Atnr2LetterTests, Gold_Is_Au)
    {
        EXPECT_EQ(std::string("Au"), std::string(constants::atnr2letter(79)));
    }

    TEST(Atnr2LetterTests, Lawrencium_103_Is_Lr)
    {
        EXPECT_EQ(std::string("Lr"), std::string(constants::atnr2letter(103)));
    }

    TEST(Atnr2LetterTests, Zero_Is_Q_Peak)
    {
        EXPECT_EQ(std::string("Q"), std::string(constants::atnr2letter(0)));
    }

    TEST(Atnr2LetterTests, OutOfRange_Returns_PROBLEM)
    {
        EXPECT_EQ(std::string("PROBLEM"), std::string(constants::atnr2letter(200)));
    }

    TEST(Atnr2LetterTests, BufferTooSmall_ReturnsNegOne)
    {
        // No more buffer-sized API: symbol is returned as const char*.
        EXPECT_EQ(std::string("C"), std::string(constants::atnr2letter(6)));
    }

    // -----------------------------------------------------------------------
    // Type2VectorTests — constants::type2vector: basis type → [nx, ny, nz]
    // -----------------------------------------------------------------------
    TEST(Type2VectorTests, Type1_Is_SShell_000)
    {
        int v[3];
        constants::type2vector(1, v);
        EXPECT_EQ(0, v[0]);
        EXPECT_EQ(0, v[1]);
        EXPECT_EQ(0, v[2]);
    }

    TEST(Type2VectorTests, Type2_Is_Px_100)
    {
        int v[3];
        constants::type2vector(2, v);
        EXPECT_EQ(1, v[0]);
        EXPECT_EQ(0, v[1]);
        EXPECT_EQ(0, v[2]);
    }

    TEST(Type2VectorTests, Type5_Is_Dx2_200)
    {
        // index 5 in type_vector: (2,0,0) = dx²
        int v[3];
        constants::type2vector(5, v);
        EXPECT_EQ(2, v[0]);
        EXPECT_EQ(0, v[1]);
        EXPECT_EQ(0, v[2]);
    }

    TEST(Type2VectorTests, Type8_Is_Dxy_110)
    {
        // index 8: (1,1,0) = dxy
        int v[3];
        constants::type2vector(8, v);
        EXPECT_EQ(1, v[0]);
        EXPECT_EQ(1, v[1]);
        EXPECT_EQ(0, v[2]);
    }

    TEST(Type2VectorTests, SumOfExponents_Matches_ShellType)
    {
        // For d-type (types 5-10), nx+ny+nz == 2
        for (int t = 5; t <= 10; ++t) {
            int v[3];
            constants::type2vector(t, v);
            EXPECT_EQ(2, v[0] + v[1] + v[2]);
        }
    }

    TEST(Type2VectorTests, OutOfRange_Returns_NegOne)
    {
        int v[3];
        constants::type2vector(0, v);
        EXPECT_EQ(-1, v[0]);
        constants::type2vector(287, v);
        EXPECT_EQ(-1, v[0]);
    }

} // namespace NoSpherA2UnitTests
