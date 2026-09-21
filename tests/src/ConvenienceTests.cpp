#include "pch.h"

#include "core/convenience.h"
#include "core/constants.h"
#include "core/basis_set.h"
#include "core/wfn_class.h"

#include <complex>
#include <cstdio>
#include <cstring>
#include <numeric>

namespace
{
    // a temp file that is removed when the test ends, so nothing is left in the temp dir or the repo
    struct TempFile
    {
        std::filesystem::path path;
        explicit TempFile(const std::string& stem, const std::string& ext = ".txt")
        {
            static int counter = 0;
            const auto tick = std::chrono::steady_clock::now().time_since_epoch().count();
            path = std::filesystem::temp_directory_path() /
                ("nos_conv_" + stem + "_" + std::to_string(tick) + "_" + std::to_string(counter++) + ext);
        }
        void write_text(const std::string& text) const
        {
            std::ofstream f(path);
            f << text;
        }
        ~TempFile()
        {
            std::error_code ec;
            std::filesystem::remove(path, ec);
        }
    };

    // spherical bessel j_l(x) from its power series, independent of the closed forms and recursion in Src
    double bessel_series(const int l, const double x)
    {
        double prefactor = 1.0;
        for (int m = 1; m <= l; m++)
            prefactor *= x / (2.0 * m + 1.0);
        double term = 1.0, sum = 1.0;
        for (int k = 1; k < 200; k++)
        {
            term *= (-x * x / 2.0) / (k * (2.0 * l + 2.0 * k + 1.0));
            sum += term;
            if (std::abs(term) < 1e-18 * std::abs(sum))
                break;
        }
        return prefactor * sum;
    }

    // Laplacian of x^a y^b z^c at a point
    double monomial_laplacian(const int a, const int b, const int c, const double x, const double y, const double z)
    {
        double res = 0.0;
        if (a >= 2) res += a * (a - 1) * pow(x, a - 2) * pow(y, b) * pow(z, c);
        if (b >= 2) res += b * (b - 1) * pow(x, a) * pow(y, b - 2) * pow(z, c);
        if (c >= 2) res += c * (c - 1) * pow(x, a) * pow(y, b) * pow(z, c - 2);
        return res;
    }

    // n!! in double; the repo's doublefactorial is unsigned int and 21!! (needed for the h Gram matrix) overflows it
    double dfact(int n)
    {
        double r = 1.0;
        for (; n > 1; n -= 2)
            r *= n;
        return r;
    }

    // the cart2sph matrices carry the Gaussian normalisation ratio sqrt(prod (2e-1)!! / (2l-1)!!) of each
    // cartesian function; dividing it out gives the plain polynomial coefficient
    double cartesian_norm_ratio(const int a, const int b, const int c)
    {
        return sqrt(dfact(2 * a - 1) * dfact(2 * b - 1) * dfact(2 * c - 1) / dfact(2 * (a + b + c) - 1));
    }

    // integral of x^a y^b z^c over the unit sphere: 4 pi (a-1)!! (b-1)!! (c-1)!! / (a+b+c+1)!! for even exponents, 0 otherwise
    double sphere_monomial_integral(const int a, const int b, const int c)
    {
        if (a % 2 || b % 2 || c % 2)
            return 0.0;
        return 4.0 * constants::PI * dfact(a - 1) * dfact(b - 1) * dfact(c - 1) / dfact(a + b + c + 1);
    }

    // every column is a harmonic polynomial (a solid harmonic), and the columns are orthogonal over the unit sphere with the
    // common norm 4 pi / (2l+1) of the Racah-normalised real solid harmonics (z^l has coefficient 1 in the m = 0 column),
    // which pins the relative and absolute scale of every block, not only its shape
    void expect_columns_harmonic(const vec2& mat, const std::vector<std::array<int, 3>>& exps, const char* what)
    {
        ASSERT_EQ(mat.size(), exps.size()) << what;
        const size_t n_sph = mat[0].size(), n_cart = mat.size();
        const int l = exps[0][0] + exps[0][1] + exps[0][2];
        vec2 coef(n_cart, vec(n_sph));
        for (size_t row = 0; row < n_cart; row++)
            for (size_t col = 0; col < n_sph; col++)
                coef[row][col] = mat[row][col] / cartesian_norm_ratio(exps[row][0], exps[row][1], exps[row][2]);
        const double pts[2][3] = { {1.3, -0.7, 0.9}, {-0.4, 1.1, 0.6} };
        for (size_t col = 0; col < n_sph; col++)
            for (const double* p : pts)
            {
                double lap = 0.0;
                for (size_t row = 0; row < n_cart; row++)
                    lap += coef[row][col] * monomial_laplacian(exps[row][0], exps[row][1], exps[row][2], p[0], p[1], p[2]);
                EXPECT_NEAR(lap, 0.0, 1e-9) << what << " column " << col;
            }
        for (size_t m = 0; m < n_sph; m++)
            for (size_t n = m; n < n_sph; n++)
            {
                double gram = 0.0;
                for (size_t i = 0; i < n_cart; i++)
                    for (size_t j = 0; j < n_cart; j++)
                        gram += coef[i][m] * coef[j][n] * sphere_monomial_integral(exps[i][0] + exps[j][0], exps[i][1] + exps[j][1], exps[i][2] + exps[j][2]);
                const double expected = m == n ? 4.0 * constants::PI / (2 * l + 1) : 0.0;
                EXPECT_NEAR(gram, expected, 1e-9) << what << " columns " << m << " and " << n;
            }
    }

    options parse(const std::vector<std::string>& args)
    {
        options opt;
        for (const auto& a : args)
            opt.arguments.push_back(a);
        opt.digest_options();
        return opt;
    }

    WFN two_atom_wfn(const double d)
    {
        WFN w(e_origin::NOT_YET_DEFINED);
        w.push_back_atom("O", 0.0, 0.0, 0.0, 8);
        w.push_back_atom("H", d, 0.0, 0.0, 1);
        return w;
    }
}

// the budget helper keeps everything when it fits, sizes a window otherwise and never returns 0 items for an oversize item
TEST(ConvenienceTests, ItemsWithinBudgetKeepsAllWhenTheyFitAndOneWhenNothingFits)
{
    EXPECT_EQ(items_within_budget(0, 8, 100), 0u);
    EXPECT_EQ(items_within_budget(10, 0, 100), 0u);
    EXPECT_EQ(items_within_budget(10, 8, 0), 0u);
    EXPECT_EQ(items_within_budget(10, 8, 100), 0u);
    EXPECT_EQ(items_within_budget(100, 8, 100), 12u);
    EXPECT_EQ(items_within_budget(1, 1000, 100), 1u);
}

// an explicit -tsc_block wins, then a -mem budget of 3 complex copies per scatterer, then the default
TEST(ConvenienceTests, TscBlockForFollowsExplicitBlockThenMemBudget)
{
    options opt;
    EXPECT_EQ(opt.tsc_block_for(100000, 50), 1000u);
    opt.mem_given = true;
    opt.mem = 1.0;
    const size_t item = 3 * 100 * sizeof(std::complex<double>);
    EXPECT_EQ(opt.tsc_block_for(10000, 100), (1024u * 1024u) / item);
    EXPECT_EQ(opt.tsc_block_for(100, 100), 0u);
    EXPECT_EQ(opt.tsc_block_for(1000000, 0), (1024u * 1024u) / (3 * sizeof(std::complex<double>)));
    opt.tsc_block_given = true;
    opt.tsc_block_size = 250;
    EXPECT_EQ(opt.tsc_block_for(10000, 100), 250u);
}

// the bounds-checked string vectors throw instead of reading past the end, and the option flavour is an out_of_range too
TEST(ConvenienceTests, SvecAndCheckedSvecThrowPastTheEnd)
{
    std::vector<std::string> plain{ "a", "b" };
    svec s(plain);
    EXPECT_EQ(s[1], "b");
    EXPECT_THROW(s[2], std::out_of_range);
    const svec& cs = s;
    EXPECT_EQ(cs[0], "a");
    EXPECT_THROW(cs[5], std::out_of_range);
    svec moved(std::vector<std::string>{ "x" });
    EXPECT_EQ(moved.size(), 1u);

    options::checked_svec args;
    args.push_back("-acc");
    EXPECT_EQ(args[0], "-acc");
    EXPECT_THROW(args[1], options::missing_argument);
    EXPECT_THROW(args[1], std::out_of_range);
}

// the primitive constructor caches exp^(l+3/2), the shell norm and the normalised coefficient; the s case has a hand value and
// every norm satisfies the radial identity N^2 int_0^inf r^(2l+2) exp(-2 a r^2) dr = N^2 Gamma(l+3/2) / (2 (2a)^(l+3/2)) = 1
TEST(ConvenienceTests, PrimitiveConstructorCachesNormalisationAndCompares)
{
    for (int l = 0; l <= 4; l++)
        for (const double a : { 0.3, 1.0, 5.5 })
        {
            const primitive p(0, l, a, 1.0);
            const double radial = std::tgamma(l + 1.5) / (2.0 * pow(2.0 * a, l + 1.5));
            EXPECT_NEAR(p.normalization_constant() * p.normalization_constant() * radial, 1.0, 1e-12) << "l = " << l << " a = " << a;
        }
    primitive s(0, 0, 1.0, 0.5);
    EXPECT_NEAR(s.normalization_constant(), pow(128.0 / constants::PI, 0.25), 1e-12);
    EXPECT_NEAR(s.get_exp_l_plus_3_2(), 1.0, 1e-15);
    EXPECT_NEAR(s.get_normalized_coefficient(), 0.5 * pow(128.0 / constants::PI, 0.25), 1e-12);

    primitive d(3, 2, 1.7, 0.4);
    EXPECT_EQ(d.get_center(), 3);
    EXPECT_EQ(d.get_type(), 2);
    EXPECT_NEAR(d.get_exp(), 1.7, 1e-15);
    EXPECT_NEAR(d.get_coef(), 0.4, 1e-15);
    EXPECT_NEAR(d.get_exp_l_plus_3_2(), pow(1.7, 3.5), 1e-12);
    const double expected_norm = pow(pow(2.0, 15) * pow(1.7, 7) / constants::PI / (15.0 * 15.0), 0.25);
    EXPECT_NEAR(d.normalization_constant(), expected_norm, 1e-12);
    EXPECT_NEAR(d.get_normalized_coefficient(), 0.4 * expected_norm, 1e-12);
    EXPECT_NEAR(d.eval_gaussian(0.5), pow(0.5, 2) * std::exp(-1.7 * 0.25) * 0.4 * expected_norm, 1e-12);

    SimplePrimitive sp{ 3, 2, 1.7, 0.4, 0 };
    primitive from_simple(sp);
    EXPECT_TRUE(from_simple == d);
    primitive other = d;
    other.set_exp(1.8);
    EXPECT_FALSE(other == d);

    primitive tweak = d;
    tweak.set_norm_const(2.0);
    EXPECT_NEAR(tweak.get_normalized_coefficient(), 0.8, 1e-15);

    ECP_primitive ecp(1, 1, 2.0, 0.3, 2);
    EXPECT_EQ(ecp.n, 2);
    EXPECT_EQ(ecp.get_type(), 1);
    EXPECT_NEAR(ecp.get_exp_l_plus_3_2(), pow(2.0, 2.5), 1e-12);
}

// is_similar takes a log10 exponent, is_similar_rel a relative fraction, is_similar_abs a plain difference
TEST(ConvenienceTests, IsSimilarFamilyInterpretsToleranceDifferently)
{
    EXPECT_TRUE(is_similar(1.0, 1.0005, -3));
    EXPECT_FALSE(is_similar(1.0, 1.002, -3));
    EXPECT_TRUE(is_similar_abs(1.0, 1.4, 0.5));
    EXPECT_FALSE(is_similar_abs(1.0, 1.6, 0.5));
    EXPECT_TRUE(is_similar_rel(100.0, 101.0, 0.05));
    EXPECT_FALSE(is_similar_rel(100.0, 110.0, 0.05));
}

// the small vector helpers against hand values
TEST(ConvenienceTests, VectorHelpersMatchHandValues)
{
    const d3 x{ 1.0, 0.0, 0.0 }, y{ 0.0, 1.0, 0.0 };
    const d3 z = vec_cross(x, y);
    EXPECT_NEAR(z[0], 0.0, 1e-15);
    EXPECT_NEAR(z[1], 0.0, 1e-15);
    EXPECT_NEAR(z[2], 1.0, 1e-15);
    const d3 a{ 1.0, 2.0, 3.0 }, b{ 4.0, -5.0, 6.0 };
    EXPECT_NEAR(vec_dot(a, b), 12.0, 1e-15);
    const d3 diff = vec_diff(a, b);
    EXPECT_NEAR(diff[0], -3.0, 1e-15);
    EXPECT_NEAR(diff[1], 7.0, 1e-15);
    EXPECT_NEAR(diff[2], -3.0, 1e-15);
    EXPECT_NEAR(vec_length(vec{ 3.0, 4.0 }), 5.0, 1e-15);
    EXPECT_NEAR(array_length(d3{ 1.0, 1.0, 1.0 }, d3{ 2.0, 3.0, 3.0 }), 3.0, 1e-15);
    EXPECT_EQ(vec_sum(bvec{ true, false, true }), 2);
    EXPECT_EQ(vec_sum(ivec{ 1, 2, 3 }), 6);
    EXPECT_NEAR(vec_sum(vec{ 0.5, 0.25 }), 0.75, 1e-15);
    const cdouble cs = vec_sum(cvec{ cdouble(1.0, 2.0), cdouble(3.0, -1.0) });
    EXPECT_NEAR(cs.real(), 4.0, 1e-15);
    EXPECT_NEAR(cs.imag(), 1.0, 1e-15);
}

// fast_exp_neg is exact near zero, zero below -42 and the (1+x/n)^n approximation between; from n ln(1+x/n) = x - x^2/(2n) + O(x^3/n^2)
// its relative error is x^2/(2n) to within a few percent for |x| <= 10, so the band [0.9, 1.05] x^2/2048 pins n = 1024:
// one squaring fewer or more doubles or halves the error and a plain exp gives none
TEST(ConvenienceTests, FastExpNegTracksExpWithinItsLeadingError)
{
    EXPECT_DOUBLE_EQ(fast_exp_neg(-0.1), exp(-0.1));
    EXPECT_DOUBLE_EQ(fast_exp_neg(0.0), 1.0);
    EXPECT_EQ(fast_exp_neg(-50.0), 0.0);
    for (const double x : { -1.0, -3.0, -5.0, -10.0 })
    {
        const double ref = exp(x);
        const double rel_err = std::abs(fast_exp_neg(x) - ref) / ref;
        const double leading = x * x / 2048.0;
        EXPECT_GT(rel_err, 0.9 * leading) << "x = " << x;
        EXPECT_LT(rel_err, 1.05 * leading) << "x = " << x;
        EXPECT_LT(fast_exp_neg(x), ref) << "(1+x/n)^n approaches exp(x) from below for x < 0";
    }
}

// sha256 against three published digests, the third one crossing the 64 byte block boundary
TEST(ConvenienceTests, Sha256MatchesPublishedDigests)
{
    EXPECT_EQ(sha::sha256(""), "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
    EXPECT_EQ(sha::sha256("abc"), "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    EXPECT_EQ(sha::sha256("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"),
        "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1");
}

// the string helpers used by the readers: trim, split, empty removal, word count, suffix test, lower case
TEST(ConvenienceTests, StringHelpersTrimSplitCountAndCompare)
{
    EXPECT_EQ(trim("  a b \t"), "a b");
    EXPECT_EQ(trim(""), "");
    EXPECT_EQ(trim("x"), "x");
    const ivec ints = split_string<int>("1,2,3", ",");
    ASSERT_EQ(ints.size(), 3u);
    EXPECT_EQ(ints[2], 3);
    const vec dbls = split_string<double>("0.5 1.5", " ");
    ASSERT_EQ(dbls.size(), 2u);
    EXPECT_NEAR(dbls[1], 1.5, 1e-15);
    svec words = split_string<std::string>("a  b", " ");
    EXPECT_EQ(words.size(), 3u);
    remove_empty_elements(words);
    ASSERT_EQ(words.size(), 2u);
    EXPECT_EQ(words[1], "b");
    svec custom{ "x", "-", "y", "" };
    remove_empty_elements(custom, "-");
    EXPECT_EQ(custom.size(), 2u);
    EXPECT_EQ(CountWords("  one two\tthree\n"), 3);
    EXPECT_EQ(CountWords(""), 0);
    EXPECT_EQ(CountWords(nullptr), -1);
    EXPECT_TRUE(ends_with("file.tscb", ".tscb"));
    EXPECT_FALSE(ends_with("file.tsc", ".tscb"));
    EXPECT_FALSE(ends_with("b", "tscb"));
    EXPECT_EQ(asciitolower('Q'), 'q');
    EXPECT_EQ(asciitolower('q'), 'q');
    EXPECT_EQ(asciitolower('1'), '1');
    EXPECT_NEAR(double_from_string_with_esd("1.234(5)"), 1.234, 1e-15);
    EXPECT_NEAR(double_from_string_with_esd("-0.5"), -0.5, 1e-15);
}

// shrink_string drops spaces, digits and brackets in place and returns the result
TEST(ConvenienceTests, ShrinkStringStripsDigitsSpacesAndBrackets)
{
    std::string label = "C12 (a) 0";
    EXPECT_EQ(shrink_string(label), "Ca");
    EXPECT_EQ(label, "Ca");
    std::string clean = "Fe";
    EXPECT_EQ(shrink_string(clean), "Fe");
}

// shrink_string_to_atom returns the element symbol of the given Z once the label is stripped
TEST(ConvenienceTests, ShrinkStringToAtomReturnsTheElementSymbol)
{
    std::string c = "C12(a)";
    EXPECT_EQ(shrink_string_to_atom(c, 6), "C");
    std::string cl = "Cl1";
    EXPECT_EQ(shrink_string_to_atom(cl, 17), "Cl");
    std::string h = "H3A";
    EXPECT_EQ(shrink_string_to_atom(h, 1), "H");
}

// a label that does not start with the symbol of its atom number gives the symbol itself
TEST(ConvenienceTests, ShrinkStringToAtomSymbolAtIndexOne)
{
    std::string label = "XCl";
    EXPECT_EQ(shrink_string_to_atom(label, 17), "Cl");
}

// the shell tables: cartesian counts (l+1)(l+2)/2 for l >= 0, spherical 2|l|+1 (sp = 4) below, and the Gaussian f ordering in shell2function
TEST(ConvenienceTests, ShellTablesMatchTheirCounts)
{
    for (int l = 0; l <= 5; l++)
        EXPECT_EQ(sht2nbas(l), (l + 1) * (l + 2) / 2) << "l = " << l;
    EXPECT_EQ(sht2nbas(-1), 4);
    EXPECT_EQ(sht2nbas(-2), 5);
    EXPECT_EQ(sht2nbas(-3), 7);
    EXPECT_EQ(sht2nbas(-4), 9);
    EXPECT_EQ(sht2nbas(-5), 11);
    EXPECT_EQ(shell2function(0, 4), 1);
    EXPECT_EQ(shell2function(1, 2), 4);
    EXPECT_EQ(shell2function(2, 0), 5);
    EXPECT_EQ(shell2function(3, 3), 17);
    EXPECT_EQ(shell2function(3, 9), 20);
    EXPECT_EQ(shell2function(4, 1), 22);
    EXPECT_EQ(shell2function(-1, 2), 3);
    EXPECT_EQ(shell2function(-2, 1), -4);
    EXPECT_EQ(shell2function(7, 0), 0);
    // the f permutation covers 11..20 exactly once
    ivec f;
    for (int p = 0; p < 10; p++)
        f.push_back(shell2function(3, p));
    std::sort(f.begin(), f.end());
    for (int p = 0; p < 10; p++)
        EXPECT_EQ(f[p], 11 + p);
}

// the constexpr double factorial and the remove-all helper
TEST(ConvenienceTests, DoublefactorialAndRemoveElement)
{
    static_assert(doublefactorial(7) == 105u, "7!! = 105");
    static_assert(doublefactorial(0) == 1u, "0!! = 1");
    static_assert(doublefactorial(-1) == 1u, "(-1)!! = 1");
    EXPECT_EQ(doublefactorial(9), 945u);
    ivec v{ 1, 2, 1, 3, 1 };
    removeElement(v, 1);
    ASSERT_EQ(v.size(), 2u);
    EXPECT_EQ(v[0], 2);
    EXPECT_EQ(v[1], 3);
}

// the CIF esd parser: digits after the point scale the bracket, no point gives 0.001, no bracket 0.005
TEST(ConvenienceTests, DecimalPrecisionFromCifNumberCoversTheThreeBranches)
{
    std::string a = "1.234(5)";
    EXPECT_NEAR(get_decimal_precision_from_CIF_number(a), 0.005, 1e-15);
    std::string b = "0.12(12)";
    EXPECT_NEAR(get_decimal_precision_from_CIF_number(b), 0.12, 1e-15);
    std::string c = "1.(5)";
    EXPECT_NEAR(get_decimal_precision_from_CIF_number(c), 0.001, 1e-15);
    std::string d = "1.5";
    EXPECT_NEAR(get_decimal_precision_from_CIF_number(d), 0.005, 1e-15);
    // without a decimal point the bracket is the esd of the integer
    std::string e = "12(34)";
    EXPECT_NEAR(get_decimal_precision_from_CIF_number(e), 34.0, 1e-15);
}

// the esd of an integer is the bracket value itself
TEST(ConvenienceTests, DecimalPrecisionOfIntegerEsd)
{
    std::string c = "12(3)";
    EXPECT_NEAR(get_decimal_precision_from_CIF_number(c), 3.0, 1e-15);
}

// the centred text helpers pad to the bar width, the bracketed one loses one column to the closing bracket
TEST(ConvenienceTests, PrintCenteredTextAndMessagePadToWidth)
{
    std::ostringstream bar, msg;
    int width = 10;
    print_centered_text("ab", width, bar);
    EXPECT_EQ(bar.str(), "[    ab   ]\n");
    print_centered_message("ab", 10, msg);
    EXPECT_EQ(msg.str(), "    ab   \n");
}

// one update per item and a redraw only when a percent boundary is crossed: 1000 items give 100 serialised writes
TEST(ConvenienceTests, ProgressBarWritesOncePerPercent)
{
    std::ostringstream out;
    const bool old_report = ProgressBar::report_counts;
    ProgressBar::report_counts = true;
    {
        ProgressBar bar(1000, 20, "#", " ", "work", out);
        EXPECT_EQ(out.str().substr(0, 1), "[");
        for (int i = 0; i < 1000; i++)
            bar.update();
        EXPECT_EQ(bar.update_calls(), 1000u);
        EXPECT_EQ(bar.bar_writes(), 100u);
    }
    ProgressBar::report_counts = old_report;
    const std::string text = out.str();
    EXPECT_NE(text.find("100%"), std::string::npos);
    EXPECT_NE(text.find("[progress] work: 1000 updates, 100 serialised writes"), std::string::npos);
}

// a batched update crossing several boundaries at once still writes only once
TEST(ConvenienceTests, ProgressBarBatchedUpdateWritesOnce)
{
    std::ostringstream out;
    {
        ProgressBar bar(1000, 10, "#", " ", "", out);
        bar.update(3);
        EXPECT_EQ(bar.bar_writes(), 0u) << "3 of 1000 does not reach the first percent";
        bar.update(7);
        EXPECT_EQ(bar.bar_writes(), 1u);
        bar.update(50);
        EXPECT_EQ(bar.bar_writes(), 2u) << "five percent boundaries in one call, one write";
        bar.update(1);
        EXPECT_EQ(bar.bar_writes(), 2u);
        EXPECT_EQ(bar.update_calls(), 4u);
    }
    EXPECT_NE(out.str().find("100%"), std::string::npos);
}

// the contributor block is the part -no_date suppresses; the banner stays
TEST(ConvenienceTests, MessageOmitsContributorsWithNoDate)
{
    const std::string full = NoSpherA2_message(false);
    const std::string bare = NoSpherA2_message(true);
    EXPECT_NE(full.find("Florian Kleemiss,"), std::string::npos);
    EXPECT_EQ(bare.find("Florian Kleemiss,"), std::string::npos);
    EXPECT_EQ(full.substr(0, bare.size()), bare);
    EXPECT_NE(bare.find("BSD-2"), std::string::npos);
}

// every column of the cartesian-to-spherical matrices is a harmonic polynomial once the cartesian norms are divided out
TEST(ConvenienceMathTests, Cart2SphColumnsAreHarmonicPolynomials)
{
    vec2 d, f, g, h;
    ASSERT_TRUE(generate_cart2sph_mat(d, f, g, h));
    ASSERT_EQ(d.size(), 6u);
    ASSERT_EQ(d[0].size(), 5u);
    ASSERT_EQ(f.size(), 10u);
    ASSERT_EQ(f[0].size(), 7u);
    ASSERT_EQ(g.size(), 15u);
    ASSERT_EQ(g[0].size(), 9u);
    ASSERT_EQ(h.size(), 21u);
    ASSERT_EQ(h[0].size(), 11u);
    expect_columns_harmonic(d, { {2,0,0}, {0,2,0}, {0,0,2}, {1,1,0}, {1,0,1}, {0,1,1} }, "d");
    expect_columns_harmonic(f, { {3,0,0}, {0,3,0}, {0,0,3}, {1,2,0}, {2,1,0}, {2,0,1}, {1,0,2}, {0,1,2}, {0,2,1}, {1,1,1} }, "f");
    expect_columns_harmonic(g, { {0,0,4}, {0,1,3}, {0,2,2}, {0,3,1}, {0,4,0}, {1,0,3}, {1,1,2}, {1,2,1}, {1,3,0}, {2,0,2}, {2,1,1}, {2,2,0}, {3,0,1}, {3,1,0}, {4,0,0} }, "g");
    expect_columns_harmonic(h, { {0,0,5}, {0,1,4}, {0,2,3}, {0,3,2}, {0,4,1}, {0,5,0}, {1,0,4}, {1,1,3}, {1,2,2}, {1,3,1}, {1,4,0},
        {2,0,3}, {2,1,2}, {2,2,1}, {2,3,0}, {3,0,2}, {3,1,1}, {3,2,0}, {4,0,1}, {4,1,0}, {5,0,0} }, "h");
}

// the d block spelled out: D0 = zz - (xx + yy)/2, D+2 = sqrt(3)/2 (xx - yy), the m = +-1 and -2 columns are single products
TEST(ConvenienceMathTests, Cart2SphDBlockMatchesTheRealSolidHarmonics)
{
    vec2 d, f, g, h;
    generate_cart2sph_mat(d, f, g, h);
    EXPECT_NEAR(d[0][0], -0.5, 1e-15);
    EXPECT_NEAR(d[1][0], -0.5, 1e-15);
    EXPECT_NEAR(d[2][0], 1.0, 1e-15);
    EXPECT_NEAR(d[4][1], 1.0, 1e-15);
    EXPECT_NEAR(d[5][2], 1.0, 1e-15);
    EXPECT_NEAR(d[0][3], sqrt(3.0) / 2.0, 1e-15);
    EXPECT_NEAR(d[1][3], -sqrt(3.0) / 2.0, 1e-15);
    EXPECT_NEAR(d[3][4], 1.0, 1e-15);
    double total = 0.0;
    for (const auto& row : d)
        for (const double v : row)
            total += std::abs(v);
    EXPECT_NEAR(total, 5.0 + sqrt(3.0), 1e-14);
}

// the median eigenvalue: diagonal shortcut, a 2x2 block with known eigenvalues 1 and 3, and a rotated diag(1,2,3)
TEST(ConvenienceMathTests, GetLambda1ReturnsTheMedianEigenvalue)
{
    double diag[9] = { 5.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 3.0 };
    EXPECT_NEAR(get_lambda_1(diag), 3.0, 1e-15);
    double block[9] = { 2.0, 1.0, 0.0, 1.0, 2.0, 0.0, 0.0, 0.0, 5.0 };
    EXPECT_NEAR(get_lambda_1(block), 3.0, 1e-12);
    double rank1[9] = { 2.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0 };
    EXPECT_NEAR(get_lambda_1(rank1), 1.0, 1e-12);

    const double th = 0.7, ph = 1.1;
    const double R[3][3] = {
        { cos(th), -sin(th) * cos(ph), sin(th) * sin(ph) },
        { sin(th), cos(th) * cos(ph), -cos(th) * sin(ph) },
        { 0.0, sin(ph), cos(ph) } };
    const double lam[3] = { 1.0, 2.0, 3.0 };
    double rotated[9];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            double s = 0.0;
            for (int k = 0; k < 3; k++)
                s += R[i][k] * lam[k] * R[j][k];
            rotated[3 * i + j] = s;
        }
    EXPECT_NEAR(get_lambda_1(rotated), 2.0, 1e-12);
}

// the closed forms l = 0..6 and the recursion for l >= 7 against the power series; x = 0 is a branch of its own
TEST(ConvenienceMathTests, BesselFirstKindMatchesThePowerSeries)
{
    EXPECT_DOUBLE_EQ(bessel_first_kind(0, 0.0), 1.0);
    EXPECT_DOUBLE_EQ(bessel_first_kind(3, 0.0), 0.0);
    for (int l = 0; l <= 9; l++)
        for (const double x : { 1.5, 4.0, 9.0 })
        {
            const double ref = bessel_series(l, x);
            EXPECT_NEAR(bessel_first_kind(l, x), ref, 1e-8 * std::abs(ref)) << "l = " << l << " x = " << x;
        }
    // small x: the closed forms l >= 3 cancel catastrophically there (not tested), the recursion does not
    for (const int l : { 0, 1, 2, 7, 8, 9 })
    {
        const double ref = bessel_series(l, 0.3);
        EXPECT_NEAR(bessel_first_kind(l, 0.3), ref, 1e-8 * std::abs(ref)) << "l = " << l;
    }
    // at x = pi j0 vanishes, which sends the recursion through its j1 normalisation branch
    EXPECT_NEAR(bessel_first_kind(7, constants::PI), bessel_series(7, constants::PI), 1e-12);
    EXPECT_NEAR(bessel_first_kind(1, 2.0), (sin(2.0) / 2.0 - cos(2.0)) / 2.0, 1e-15);
}

// 2F1(1,1;2;x) = -ln(1-x)/x and 2F1(a,b;b;x) = (1-x)^-a, in double and complex
TEST(ConvenienceMathTests, HypergeometricMatchesClosedForms)
{
    for (const double x : { 0.3, -0.5, 0.05 })
        EXPECT_NEAR(hypergeometric(1.0, 1.0, 2.0, x), -std::log(1.0 - x) / x, 1e-9) << "x = " << x;
    EXPECT_NEAR(hypergeometric(0.5, 1.0, 1.0, 0.4), 1.0 / sqrt(0.6), 1e-9);
    const cdouble z(0.2, 0.3);
    const cdouble got = hypergeometric(1.0, 1.0, 2.0, z);
    const cdouble ref = -std::log(cdouble(1.0) - z) / z;
    EXPECT_NEAR(got.real(), ref.real(), 1e-9);
    EXPECT_NEAR(got.imag(), ref.imag(), 1e-9);
    const cdouble pw = hypergeometric(2.0, 1.0, 1.0, z);
    const cdouble pw_ref = 1.0 / ((cdouble(1.0) - z) * (cdouble(1.0) - z));
    EXPECT_NEAR(pw.real(), pw_ref.real(), 1e-9);
    EXPECT_NEAR(pw.imag(), pw_ref.imag(), 1e-9);
}

// swap_sort puts the values in ascending key order; the multi variant does the same along every row of a table.
// The values are deliberately not monotone in the key, so sorting the values themselves or gathering v[order[i]] both fail
TEST(ConvenienceMathTests, SwapSortOrdersValuesByKey)
{
    ivec order{ 2, 0, 3, 1 };
    cvec v{ cdouble(5, 1), cdouble(7, 2), cdouble(1, 3), cdouble(9, 4) };
    swap_sort(order, v);
    const double expected_re[4] = { 7.0, 9.0, 5.0, 1.0 }, expected_im[4] = { 2.0, 4.0, 1.0, 3.0 };
    for (int i = 0; i < 4; i++)
    {
        EXPECT_NEAR(v[i].real(), expected_re[i], 1e-15) << "position " << i;
        EXPECT_NEAR(v[i].imag(), expected_im[i], 1e-15) << "position " << i;
    }

    ivec order3{ 1, 2, 0 };
    ivec2 table{ { 7, 5, 9 }, { 17, 15, 19 }, { 27, 25, 29 } };
    swap_sort_multi(order3, table);
    for (int row = 0; row < 3; row++)
    {
        EXPECT_EQ(table[row][0], 9 + 10 * row);
        EXPECT_EQ(table[row][1], 7 + 10 * row);
        EXPECT_EQ(table[row][2], 5 + 10 * row);
    }
}

// the cube box is the atom extent padded by the radius (Angstrom to bohr) and stepped at the resolution in Angstrom
TEST(ConvenienceMathTests, ReadxyzMinMaxFromWFNPadsAndSteps)
{
    WFN w(e_origin::NOT_YET_DEFINED);
    w.push_back_atom("O", 0.0, 0.0, 0.0, 8);
    w.push_back_atom("H", 2.0, 1.0, -0.5, 1);
    w.push_back_atom("H", -1.0, 3.0, 0.5, 1);
    properties_options opts;
    opts.radius = 2.0;
    opts.resolution = 0.1;
    readxyzMinMax_fromWFN(w, opts);
    const double pad = constants::ang2bohr(2.0);
    EXPECT_NEAR(opts.MinMax[0], -1.0 - pad, 1e-12);
    EXPECT_NEAR(opts.MinMax[3], 2.0 + pad, 1e-12);
    EXPECT_NEAR(opts.MinMax[1], 0.0 - pad, 1e-12);
    EXPECT_NEAR(opts.MinMax[4], 3.0 + pad, 1e-12);
    EXPECT_NEAR(opts.MinMax[2], -0.5 - pad, 1e-12);
    EXPECT_NEAR(opts.MinMax[5], 0.5 + pad, 1e-12);
    EXPECT_EQ(opts.NbSteps[0], (int)ceil(constants::bohr2ang(3.0 + 2.0 * pad) / 0.1));
    EXPECT_EQ(opts.NbSteps[1], (int)ceil(constants::bohr2ang(3.0 + 2.0 * pad) / 0.1));
    EXPECT_EQ(opts.NbSteps[2], (int)ceil(constants::bohr2ang(1.0 + 2.0 * pad) / 0.1));
    EXPECT_EQ(opts.n_grid_points(), size_t(opts.NbSteps[0]) * opts.NbSteps[1] * opts.NbSteps[2]);
    EXPECT_FALSE(opts.calc());
    opts.rho = true;
    EXPECT_TRUE(opts.calc());
}

// a shortest interatomic distance below 2 reads as Angstrom, above as bohr
TEST(ConvenienceMathTests, CheckBohrDecidesOnTheShortestDistance)
{
    EXPECT_FALSE(check_bohr(two_atom_wfn(1.5), false));
    EXPECT_TRUE(check_bohr(two_atom_wfn(2.5), false));
    WFN three = two_atom_wfn(5.0);
    three.push_back_atom("H", 5.0, 1.0, 0.0, 1);
    EXPECT_FALSE(check_bohr(three, false));
}

// go_get_string returns the first line holding the tag, "" at the end, and continues without rewind
TEST(ConvenienceIoTests, GoGetStringFindsLinesAndReportsEof)
{
    TempFile tf("gogets");
    tf.write_text("alpha line\r\nbeta line\ngamma line\n");
    std::ifstream f(tf.path);
    ASSERT_TRUE(f.is_open());
    EXPECT_EQ(go_get_string(f, "beta"), "beta line");
    EXPECT_EQ(go_get_string(f, "alpha"), "alpha line");
    EXPECT_EQ(go_get_string(f, "gamma", false), "gamma line");
    EXPECT_EQ(go_get_string(f, "zeta"), "");
}

// fortran records round-trip through both readers, and a damaged trailer or an oversize record is refused
TEST(ConvenienceIoTests, FortranBinaryRecordsRoundTripAndRejectDamage)
{
    const double payload[3] = { 1.5, -2.25, 1e-3 };
    const int len = 3 * sizeof(double);
    TempFile good("fbin_good", ".bin"), bad("fbin_bad", ".bin");
    {
        std::ofstream o(good.path, std::ios::binary);
        o.write(reinterpret_cast<const char*>(&len), sizeof(int));
        o.write(reinterpret_cast<const char*>(payload), len);
        o.write(reinterpret_cast<const char*>(&len), sizeof(int));
        o.write(reinterpret_cast<const char*>(&len), sizeof(int));
        o.write(reinterpret_cast<const char*>(payload), len);
        o.write(reinterpret_cast<const char*>(&len), sizeof(int));
    }
    {
        const int wrong = len + 4;
        std::ofstream o(bad.path, std::ios::binary);
        o.write(reinterpret_cast<const char*>(&len), sizeof(int));
        o.write(reinterpret_cast<const char*>(payload), len);
        o.write(reinterpret_cast<const char*>(&wrong), sizeof(int));
    }
    std::ifstream f(good.path, std::ios::binary);
    vec v;
    ASSERT_TRUE(read_block_from_fortran_binary(f, v));
    ASSERT_EQ(v.size(), 3u);
    EXPECT_NEAR(v[1], -2.25, 1e-15);
    double raw[3] = { 0.0, 0.0, 0.0 };
    ASSERT_TRUE(read_block_from_fortran_binary(f, raw, sizeof(raw)));
    EXPECT_NEAR(raw[2], 1e-3, 1e-15);
    EXPECT_FALSE(read_block_from_fortran_binary(f, v)) << "nothing left to read";

    std::ifstream b(bad.path, std::ios::binary);
    vec w;
    EXPECT_FALSE(read_block_from_fortran_binary(b, w));
    std::ifstream b2(bad.path, std::ios::binary);
    double tiny[1] = { 0.0 };
    EXPECT_FALSE(read_block_from_fortran_binary(b2, tiny, sizeof(tiny))) << "record larger than capacity";
    std::ifstream b3(bad.path, std::ios::binary);
    double room[3] = { 0.0, 0.0, 0.0 };
    EXPECT_FALSE(read_block_from_fortran_binary(b3, room, sizeof(room))) << "trailer mismatch";
}

// the line readers: getline_universal strips \r, seek_line lands on the tag, append_numbers parses a whole line
TEST(ConvenienceIoTests, LineReadersStripCarriageReturnAndSeek)
{
    std::istringstream is("first\r\nsecond\nTAG here\n1 2.5 -3\n");
    std::string line;
    getline_universal(is, line);
    EXPECT_EQ(line, "first");
    seek_line(is, line, "TAG", std::cout);
    EXPECT_EQ(line, "TAG here");
    read_line_or_fail(is, line, "numbers", std::cout);
    vec numbers;
    append_numbers(line, numbers, "numbers", std::cout);
    ASSERT_EQ(numbers.size(), 3u);
    EXPECT_NEAR(numbers[1], 2.5, 1e-15);
    EXPECT_NEAR(numbers[2], -3.0, 1e-15);
    ivec ints;
    append_numbers(std::string("4 5"), ints, "ints", std::cout);
    ASSERT_EQ(ints.size(), 2u);
    EXPECT_EQ(ints[1], 5);
}

// the timing table names every interval and the total, and an empty input is reported instead of dereferenced
TEST(ConvenienceIoTests, WriteTimingToFileListsIntervalsAndTotal)
{
    const _time_point t0 = get_time();
    const _time_point t1 = t0 + std::chrono::milliseconds(1500);
    const _time_point t2 = t1 + std::chrono::milliseconds(500);
    EXPECT_EQ(get_msec(t0, t1), 1500);
    EXPECT_EQ(get_musec(t1, t2), 500000);
    EXPECT_EQ(get_sec(t0, t2), 2);
    std::ostringstream out;
    write_timing_to_file(out, { t0, t1, t2 }, { "setup", "run" });
    const std::string text = out.str();
    EXPECT_NE(text.find("Time Breakdown"), std::string::npos);
    EXPECT_NE(text.find("... for setup"), std::string::npos);
    EXPECT_NE(text.find("00:01:500  (75.00%)"), std::string::npos);
    EXPECT_NE(text.find("... for run"), std::string::npos);
    EXPECT_NE(text.find("00:00:500  (25.00%)"), std::string::npos);
    EXPECT_NE(text.find("Total Time"), std::string::npos);
    EXPECT_NE(text.find("00:02:000"), std::string::npos);
    std::ostringstream empty;
    write_timing_to_file(empty, {}, {});
    EXPECT_NE(empty.str().find("Empty vector"), std::string::npos);
}

// look_for_debug keeps every argv token and turns -debug into the debug flag and the progress-bar counters
TEST(ConvenienceOptionsTests, LookForDebugStoresArgvAndSetsDebug)
{
    const bool old_report = ProgressBar::report_counts;
    std::vector<std::string> tokens{ "NoSpherA2", "-debug", "-acc", "3" };
    std::vector<char*> argv;
    for (auto& t : tokens)
        argv.push_back(t.data());
    int argc = (int)argv.size();
    std::ostringstream log;
    options opt(argc, argv.data(), log);
    EXPECT_TRUE(opt.debug);
    EXPECT_TRUE(ProgressBar::report_counts);
    ASSERT_EQ(opt.arguments.size(), 4u);
    EXPECT_EQ(opt.arguments[3], "3");
    ProgressBar::report_counts = old_report;

    std::vector<std::string> quiet{ "NoSpherA2", "-acc", "3" };
    std::vector<char*> argv2;
    for (auto& t : quiet)
        argv2.push_back(t.data());
    int argc2 = (int)argv2.size();
    options opt2(argc2, argv2.data(), log);
    EXPECT_FALSE(opt2.debug);
    opt2.digest_options();
    EXPECT_EQ(opt2.accuracy, 3);
}

// bare values and unknown flags are skipped, the flags after them still parse
TEST(ConvenienceOptionsTests, DigestSkipsBareValuesAndUnknownFlags)
{
    const options opt = parse({ "12", "-unknown_flag_xyz", "-acc", "5", "-mult", "3" });
    EXPECT_EQ(opt.accuracy, 5);
    EXPECT_EQ(opt.mult, 3u);
}

// the scalar run options land in their fields; -mem beyond the addressable size is clamped to 50000 MB
TEST(ConvenienceOptionsTests, RunOptionsParseScalars)
{
    const options opt = parse({ "-acc", "3", "-charge", "-2", "-mult", "2", "-dmin", "0.5", "-mem", "512",
        "-method", "PBE", "-pbc", "1", "-radius", "3", "-resolution", "0.2", "-e_field", "0.01",
        "-basin_grid", "0", "-refine", "0.05", "-skpts", "-rkpts", "-profiling", "myroot", "-b", "def2-SVP" });
    EXPECT_EQ(opt.accuracy, 3);
    EXPECT_EQ(opt.charge, -2);
    EXPECT_EQ(opt.mult, 2u);
    EXPECT_NEAR(opt.dmin, 0.5, 1e-15);
    EXPECT_NEAR(opt.mem, 512.0, 1e-15);
    EXPECT_TRUE(opt.mem_given);
    EXPECT_EQ(opt.method, "PBE");
    EXPECT_EQ(opt.pbc, 1);
    EXPECT_NEAR(opt.properties.radius, 3.0, 1e-15);
    EXPECT_NEAR(opt.properties.resolution, 0.2, 1e-15);
    EXPECT_NEAR(opt.efield, 0.01, 1e-15);
    EXPECT_EQ(opt.basin_grid, 1) << "-basin_grid is clamped to at least 1";
    EXPECT_NEAR(opt.properties.integral_accuracy, 0.05, 1e-15);
    EXPECT_TRUE(opt.save_k_pts);
    EXPECT_TRUE(opt.read_k_pts);
    EXPECT_TRUE(opt.profiling);
    EXPECT_EQ(opt.profiling_tests_root, std::filesystem::path("myroot"));
    EXPECT_EQ(opt.basis_set, "def2-SVP");

    const options huge = parse({ "-mem", "1e300" });
    EXPECT_NEAR(huge.mem, 50000.0, 1e-15);
    const options plain = parse({ "-profiling" });
    EXPECT_TRUE(plain.profiling);
    EXPECT_EQ(plain.profiling_tests_root, std::filesystem::path("tests"));
}

// -no_date hides the timings and the gpu notes, -no_date_but_gpu keeps the gpu notes; both globals are restored afterwards
TEST(ConvenienceOptionsTests, NoDateVariantsSetTheHideGlobals)
{
    const bool old_gpu = constants::hide_gpu_notes, old_timings = constants::hide_timings;
    const options a = parse({ "-no_date" });
    EXPECT_TRUE(a.no_date);
    EXPECT_TRUE(constants::hide_gpu_notes);
    EXPECT_TRUE(constants::hide_timings);
    const options b = parse({ "-no_date_but_gpu" });
    EXPECT_TRUE(b.no_date);
    EXPECT_FALSE(constants::hide_gpu_notes);
    EXPECT_TRUE(constants::hide_timings);
    constants::hide_gpu_notes = old_gpu;
    constants::hide_timings = old_timings;
}

// -group takes numbers until the next flag, a + prefix negates; -Anion/-Cation split their space separated lists
TEST(ConvenienceOptionsTests, GroupAndIonListsStopAtTheNextFlag)
{
    const options opt = parse({ "-group", "1", "+2", "3", "-acc", "4", "-Anion", "Cl Br", "I", "-Cation", "Na" });
    ASSERT_EQ(opt.groups.size(), 1u);
    ASSERT_EQ(opt.groups[0].size(), 3u);
    EXPECT_EQ(opt.groups[0][0], 1);
    EXPECT_EQ(opt.groups[0][1], -2);
    EXPECT_EQ(opt.groups[0][2], 3);
    EXPECT_EQ(opt.accuracy, 4);
    ASSERT_EQ(opt.Anions.size(), 3u);
    EXPECT_EQ(opt.Anions[0], "Cl");
    EXPECT_EQ(opt.Anions[1], "Br");
    EXPECT_EQ(opt.Anions[2], "I");
    ASSERT_EQ(opt.Cations.size(), 1u);
    EXPECT_EQ(opt.Cations[0], "Na");
}

// -hkl_min_max eats six ints, -twin nine doubles per law, -ECP an optional mode
TEST(ConvenienceOptionsTests, FixedCountOptionsConsumeTheirValues)
{
    const options opt = parse({ "-hkl_min_max", "-1", "1", "-2", "2", "-3", "3",
        "-twin", "-1", "0", "0", "0", "-1", "0", "0", "0", "1",
        "-twin", "0", "1", "0", "1", "0", "0", "0", "0", "-1",
        "-ECP", "2", "-acc", "1" });
    EXPECT_EQ(opt.hkl_min_max[0][0], -1);
    EXPECT_EQ(opt.hkl_min_max[0][1], 1);
    EXPECT_EQ(opt.hkl_min_max[1][0], -2);
    EXPECT_EQ(opt.hkl_min_max[2][1], 3);
    ASSERT_EQ(opt.twin_law.size(), 2u);
    ASSERT_EQ(opt.twin_law[0].size(), 9u);
    EXPECT_NEAR(opt.twin_law[0][0], -1.0, 1e-15);
    EXPECT_NEAR(opt.twin_law[0][8], 1.0, 1e-15);
    EXPECT_NEAR(opt.twin_law[1][1], 1.0, 1e-15);
    EXPECT_NEAR(opt.twin_law[1][8], -1.0, 1e-15);
    EXPECT_TRUE(opt.ECP);
    EXPECT_EQ(opt.ECP_mode, 2);
    EXPECT_EQ(opt.accuracy, 1);
    const options bare = parse({ "-ecp", "-acc", "1" });
    EXPECT_TRUE(bare.ECP);
    EXPECT_EQ(bare.ECP_mode, 0);
}

// the partition switches, the tsc block size, the combined-tsc charge list with its n prefix and the SALTED model dir
TEST(ConvenienceOptionsTests, PartitionOptionsSetSchemeBlocksAndLists)
{
    EXPECT_EQ(parse({ "-tfvc" }).partition_type, PartitionType::TFVC);
    EXPECT_EQ(parse({ "-MBIS" }).partition_type, PartitionType::MBIS);
    EXPECT_EQ(parse({ "-embis" }).partition_type, PartitionType::EMBIS);
    EXPECT_EQ(parse({ "-becke" }).partition_type, PartitionType::Becke);
    const options opt = parse({ "-hirsh", "3", "-tsc_block", "250", "-mtc_charge", "n2", "3", "-mtc_ECP", "1", "2",
        "-IAM", "-ED", "-def", "-HDEF", "-old_tsc", "-SALTED", "model_dir", "-xyz", "structure.xyz" });
    EXPECT_TRUE(opt.properties.hirsh);
    EXPECT_EQ(opt.properties.hirsh_number, 3);
    EXPECT_EQ(opt.tsc_block_size, 250u);
    EXPECT_TRUE(opt.tsc_block_given);
    ASSERT_EQ(opt.combined_tsc_calc_charge.size(), 2u);
    EXPECT_EQ(opt.combined_tsc_calc_charge[0], -2);
    EXPECT_EQ(opt.combined_tsc_calc_charge[1], 3);
    ASSERT_EQ(opt.combined_tsc_calc_ECP.size(), 2u);
    EXPECT_EQ(opt.combined_tsc_calc_ECP[1], 2);
    EXPECT_TRUE(opt.iam_switch);
    EXPECT_TRUE(opt.electron_diffraction);
    EXPECT_TRUE(opt.properties.def);
    EXPECT_TRUE(opt.properties.hdef);
    EXPECT_TRUE(opt.old_tsc);
    EXPECT_TRUE(opt.SALTED);
    EXPECT_EQ(opt.salted_model_dir, std::filesystem::path("model_dir"));
    EXPECT_EQ(opt.xyz_file, std::filesystem::path("structure.xyz"));
    EXPECT_EQ(opt.wfn, std::filesystem::path("structure.xyz")) << "-xyz stands in for -wfn with SALTED";
    const options both = parse({ "-SALTED", "m", "-xyz", "s.xyz", "-fchk", "f.fchk" });
    EXPECT_EQ(both.fchk, std::filesystem::path("f.fchk"));
}

// the property flags and the -MO list, with all as the everything switch
TEST(ConvenienceOptionsTests, PropertyFlagsAndMoList)
{
    const options opt = parse({ "-MO", "3", "-MO", "7", "-rho", "-lap", "-eli", "-elf", "-esp", "-rdg", "-fukui",
        "-QCT", "-get_g", "-fractal", "frac.cube", "-hirshfeld_surface", "a.wfn", "b.wfn", "-cmos1", "1", "2", "-cmos2", "3",
        "-gbw2wfn", "-wfn_cif", "-d", "basis_dir", "-test", "-promol_nci_single_thread" });
    ASSERT_EQ(opt.properties.MO_numbers.size(), 2u);
    EXPECT_EQ(opt.properties.MO_numbers[1], 7);
    EXPECT_FALSE(opt.properties.all_mos);
    EXPECT_TRUE(opt.properties.rho);
    EXPECT_TRUE(opt.properties.lap);
    EXPECT_TRUE(opt.properties.eli);
    EXPECT_TRUE(opt.properties.elf);
    EXPECT_TRUE(opt.properties.esp);
    EXPECT_TRUE(opt.properties.rdg);
    EXPECT_TRUE(opt.properties.fukui);
    EXPECT_TRUE(opt.properties.calc());
    EXPECT_TRUE(opt.qct);
    EXPECT_TRUE(opt.get_g);
    EXPECT_TRUE(opt.fract);
    EXPECT_EQ(opt.fract_name, std::filesystem::path("frac.cube"));
    EXPECT_EQ(opt.hirshfeld_surface, std::filesystem::path("a.wfn"));
    EXPECT_EQ(opt.hirshfeld_surface2, std::filesystem::path("b.wfn"));
    ASSERT_EQ(opt.cmo1.size(), 2u);
    EXPECT_EQ(opt.cmo1[1], 2);
    ASSERT_EQ(opt.cmo2.size(), 1u);
    EXPECT_EQ(opt.cmo2[0], 3);
    EXPECT_TRUE(opt.gbw2wfn);
    EXPECT_TRUE(opt.write_CIF);
    EXPECT_EQ(opt.basis_set_path, std::filesystem::path("basis_dir"));
    EXPECT_TRUE(opt.test);
    EXPECT_TRUE(opt.properties.promol_nci_single_threaded);
    const options all = parse({ "-MO", "all" });
    EXPECT_TRUE(all.properties.all_mos);
    EXPECT_TRUE(all.properties.MO_numbers.empty());
}

// -esp_isosurface defaults to 0.002 au and reads an explicit value; -eli_analysis and -fukui_analysis take an inline wavefunction
TEST(ConvenienceOptionsTests, IsosurfaceAndAnalysisOptionsReadOptionalValues)
{
    const options def = parse({ "-esp_isosurface", "-rho" });
    EXPECT_NEAR(def.properties.esp_isosurface, 0.002, 1e-15);
    EXPECT_TRUE(def.properties.rho);
    const options val = parse({ "-esp_isosurface", "0.01" });
    EXPECT_NEAR(val.properties.esp_isosurface, 0.01, 1e-15);
    const options eli = parse({ "-eli_analysis", "mol.wfn", "0.05", "3.5", "-acc", "4" });
    EXPECT_TRUE(eli.eli_analysis_run);
    EXPECT_EQ(eli.wfn, std::filesystem::path("mol.wfn"));
    EXPECT_NEAR(eli.properties.resolution, 0.05, 1e-15);
    EXPECT_NEAR(eli.properties.radius, 3.5, 1e-15);
    EXPECT_EQ(eli.accuracy, 4);
    const options fukui = parse({ "-fukui_analysis", "mol.gbw", "-acc", "1" });
    EXPECT_TRUE(fukui.fukui_analysis_run);
    EXPECT_EQ(fukui.wfn, std::filesystem::path("mol.gbw"));
    EXPECT_EQ(fukui.accuracy, 1);
    const options fukui_bare = parse({ "-fukui_analysis", "-acc", "1" });
    EXPECT_TRUE(fukui_bare.fukui_analysis_run);
    EXPECT_TRUE(fukui_bare.wfn.empty());
    const options pol = parse({ "-polarizabilities", "a", "b", "c", "d", "e", "f", "g" });
    ASSERT_EQ(pol.pol_wfns.size(), 7u);
    EXPECT_EQ(pol.pol_wfns[6], std::filesystem::path("g"));
}

// the gpu and xcw toggles flip their fields, the numeric ones read a value
TEST(ConvenienceOptionsTests, GpuAndXcwTogglesFlipTheirFields)
{
    const options opt = parse({ "-no_gpu", "-gpu_fp64", "-gpu_fp32", "-no_gpu_itensor", "-gpu_itensor_tensor", "-no_gpu_cublas",
        "-no_gpu_salted", "-salted_charge_constraint", "-no_gpu_grid", "-no_gpu_density", "-gpu_blas", "-no_cpu_itensor_fp32",
        "-itensor_hybrid", "-no_xcw_extrapolate", "-xcw_incremental", "-xcw_int_precision", "1e-12",
        "-do_XCW", "-calc_F", "-xcw_gaussian_halt", "-xcw_strong_cutoff", "2.5", "-XCW_settings", "settings.toml" });
    EXPECT_FALSE(opt.use_gpu);
    EXPECT_TRUE(opt.gpu_fp64);
    EXPECT_TRUE(opt.gpu_fp32);
    EXPECT_FALSE(opt.gpu_itensor);
    EXPECT_TRUE(opt.gpu_itensor_tensor);
    EXPECT_FALSE(opt.gpu_cublas);
    EXPECT_FALSE(opt.gpu_salted);
    EXPECT_TRUE(opt.salted_charge_constraint);
    EXPECT_FALSE(opt.gpu_grid);
    EXPECT_FALSE(opt.gpu_density);
    EXPECT_TRUE(opt.gpu_blas);
    EXPECT_FALSE(opt.cpu_itensor_fp32);
    EXPECT_TRUE(opt.itensor_hybrid);
    EXPECT_FALSE(opt.xcw_extrapolate);
    EXPECT_TRUE(opt.xcw_incremental);
    EXPECT_NEAR(opt.xcw_int_precision, 1e-12, 1e-25);
    EXPECT_TRUE(opt.do_XCW);
    EXPECT_TRUE(opt.calc_F_calc);
    EXPECT_TRUE(opt.xcw_gaussian_halt);
    EXPECT_NEAR(opt.xcw_strong_cutoff, 2.5, 1e-15);
    EXPECT_EQ(opt.xcw_settings_path, std::filesystem::path("settings.toml"));
    const options back = parse({ "-no_gpu_itensor", "-gpu_itensor", "-no_itensor_hybrid", "-xcw_extrapolate", "-no_xcw_incremental", "-gpu_cublas", "-gpu_grid", "-gpu_density", "-gpu_salted", "-cpu_itensor_fp32" });
    EXPECT_TRUE(back.gpu_itensor);
    EXPECT_FALSE(back.itensor_hybrid);
    EXPECT_TRUE(back.xcw_extrapolate);
    EXPECT_FALSE(back.xcw_incremental);
    EXPECT_TRUE(back.gpu_cublas);
    EXPECT_TRUE(back.cpu_itensor_fp32);
}

// -promol_nci collects every existing .xyz that follows and then up to five numeric cutoffs, stopping at the next flag
TEST(ConvenienceOptionsTests, PromolNciCollectsFragmentsAndCutoffs)
{
    TempFile a("fragA", ".xyz"), b("fragB", ".xyz"), c("fragC", ".xyz");
    a.write_text("1\n\nH 0 0 0\n");
    b.write_text("1\n\nH 1 0 0\n");
    c.write_text("1\n\nH 2 0 0\n");
    const options two = parse({ "-promol_nci", a.path.string(), b.path.string(), "0.9", "0.7", "-acc", "1" });
    EXPECT_TRUE(two.promol_nci);
    ASSERT_EQ(two.promol_nci_xyz.size(), 2u);
    EXPECT_EQ(two.promol_nci_xyz[1], b.path);
    EXPECT_NEAR(two.properties.promol_nci_rcut1, 0.9, 1e-15);
    EXPECT_NEAR(two.properties.promol_nci_rcut2, 0.7, 1e-15);
    EXPECT_NEAR(two.properties.promol_nci_rho_abs_max, 0.5, 1e-15) << "third cutoff keeps its default";
    EXPECT_EQ(two.accuracy, 1);
    const options three = parse({ "-promol_nci", a.path.string(), b.path.string(), c.path.string(), "-acc", "3" });
    ASSERT_EQ(three.promol_nci_xyz.size(), 3u);
    EXPECT_EQ(three.promol_nci_xyz[2], c.path);
    EXPECT_NEAR(three.properties.promol_nci_rcut1, 0.95, 1e-15);
    EXPECT_EQ(three.accuracy, 3);
    const options five = parse({ "-promol_nci", a.path.string(), b.path.string(), "0.9", "0.7", "0.4", "0.8", "0.02" });
    EXPECT_NEAR(five.properties.promol_nci_rho_abs_max, 0.4, 1e-15);
    EXPECT_NEAR(five.properties.promol_nci_rdg_max, 0.8, 1e-15);
    EXPECT_NEAR(five.properties.promol_nci_colour_max, 0.02, 1e-15);
}

// -multipole_moments picks the scheme by name, bounds the order, turns on the RI fit and gets an auto_aux basis when none was named
TEST(ConvenienceOptionsTests, MultipoleAndRepulsionOptions)
{
    const options opt = parse({ "-multipole_moments", "TFVC", "3", "-multipole_strength", "2.5", "-multipole_centre",
        "-repulsion_overlap", "0.3", "-repulsion_exchange", "b88", "-geometry_aid_cutoff", "3.0" });
    EXPECT_EQ(opt.multipole_scheme, PartitionType::TFVC);
    EXPECT_EQ(opt.multipole_lmax, 3);
    EXPECT_TRUE(opt.RI_FIT);
    EXPECT_EQ(opt.partition_type, PartitionType::RI);
    EXPECT_EQ(opt.aux_basis.size(), 1u) << "auto_aux pushed after parsing";
    EXPECT_NEAR(opt.multipole_strength, 2.5, 1e-15);
    EXPECT_FALSE(opt.multipole_partition);
    EXPECT_NEAR(opt.repulsion_overlap, 0.3, 1e-15);
    EXPECT_EQ(opt.repulsion_exchange, 2);
    EXPECT_NEAR(opt.geometry_aid_cutoff, 3.0, 1e-15);
    EXPECT_EQ(parse({ "-multipole_moments", "hirsh", "0" }).multipole_scheme, PartitionType::Hirshfeld);
    EXPECT_EQ(parse({ "-multipole_moments", "mbis", "1" }).multipole_scheme, PartitionType::MBIS);
    EXPECT_EQ(parse({ "-multipole_moments", "EMBIS", "8" }).multipole_scheme, PartitionType::EMBIS);
    EXPECT_EQ(parse({ "-repulsion_exchange", "dirac" }).repulsion_exchange, 0);
    EXPECT_EQ(parse({ "-repulsion_exchange", "pbe" }).repulsion_exchange, 1);
    EXPECT_EQ(parse({ "-repulsion_exchange", "r2scan" }).repulsion_exchange, 3);
    const options ri = parse({ "-ri_fit", "-acc", "2" });
    EXPECT_TRUE(ri.RI_FIT);
    EXPECT_EQ(ri.partition_type, PartitionType::RI);
    ASSERT_EQ(ri.aux_basis.size(), 1u) << "no name falls back to one auto_aux basis";
    EXPECT_EQ(ri.aux_basis[0]->get_primitive_count(), 0u) << "the fallback is an empty set gen_auto_aux fills later";
    EXPECT_TRUE(ri.aux_basis[0]->_auto_aux_elements.empty()) << "and not restricted to any element";
    EXPECT_EQ(ri.multipole_lmax, -1);
    const options ri_auto = parse({ "-ri_fit", "auto_aux", "H", "C", "-acc", "2" });
    ASSERT_EQ(ri_auto.aux_basis.size(), 1u);
    EXPECT_EQ(ri_auto.aux_basis[0]->_auto_aux_elements, (ivec{ 1, 6 })) << "element symbols after auto_aux are stored as Z";
    EXPECT_EQ(ri_auto.accuracy, 2);
}

// -rgbi-groups reads comma lists with ranges in both directions; every rgbi flag turns rgbi on
TEST(ConvenienceOptionsTests, RgbiGroupsParseRangesAndBasis)
{
    const options opt = parse({ "-rgbi-groups", "1-3,5", "4", "-rgbi-groups", "3-1", "-rgbi_basis", "NAO", "-rgbi_EVs" });
    EXPECT_TRUE(opt.rgbi);
    ASSERT_EQ(opt.rgbi_group_sets.size(), 2u);
    ASSERT_EQ(opt.rgbi_group_sets[0].size(), 2u);
    const ivec expected{ 1, 2, 3, 5 };
    EXPECT_EQ(opt.rgbi_group_sets[0][0], expected);
    EXPECT_EQ(opt.rgbi_group_sets[0][1], ivec{ 4 });
    const ivec descending{ 3, 2, 1 };
    ASSERT_EQ(opt.rgbi_group_sets[1].size(), 1u);
    EXPECT_EQ(opt.rgbi_group_sets[1][0], descending);
    EXPECT_EQ(opt.rgbi_orbital_basis, RGBIOrbitalBasis::NAO);
    EXPECT_TRUE(opt.rgbi_EVs);
    const options ano = parse({ "-rgbi_basis", "ano" });
    EXPECT_EQ(ano.rgbi_orbital_basis, RGBIOrbitalBasis::ANO);
    const options theta = parse({ "-rgbi_theta" });
    EXPECT_TRUE(theta.rgbi);
    EXPECT_TRUE(theta.rgbi_theta);
    const options nosym = parse({ "-rgbi_no_sym" });
    EXPECT_TRUE(nosym.rgbi);
    EXPECT_TRUE(nosym.rgbi_no_sym);
    EXPECT_EQ(parse({ "-rgbi-groups" }).rgbi_group_sets.size(), 0u) << "no group after the flag adds no set";
}

// a one-shot option sets finished and stops the parse, so nothing after it is read
TEST(ConvenienceOptionsTests, FinishedStopsTheParse)
{
    const options opt = parse({ "-charge", "1", "-lahvatest", "-acc", "9" });
    EXPECT_TRUE(opt.finished);
    EXPECT_EQ(opt.charge, 1);
    EXPECT_EQ(opt.accuracy, 2) << "-acc after the one-shot option must not be read";
}

// the file options that check existence accept a temp file that does exist
TEST(ConvenienceOptionsTests, FileOptionsAcceptExistingFiles)
{
    TempFile cif("opt", ".cif"), hkl("opt", ".hkl"), wfn("opt", ".wfn"), job("opt", ".toml"), cube("opt", ".cube");
    cif.write_text("data_x\n");
    hkl.write_text("   0   0   1  1.0  0.1\n");
    wfn.write_text("x\n");
    job.write_text("x\n");
    cube.write_text("x\n");
    const options opt = parse({ "-cif", cif.path.string(), "-hkl", hkl.path.string(), "-wfn", wfn.path.string(),
        "-interaction_energies", job.path.string(), "-cube_density", cube.path.string() });
    EXPECT_EQ(opt.cif, cif.path);
    EXPECT_EQ(opt.hkl, hkl.path);
    EXPECT_EQ(opt.wfn, wfn.path);
    EXPECT_EQ(opt.interaction_energies_job, job.path);
    EXPECT_EQ(opt.cube_density, cube.path);
}
