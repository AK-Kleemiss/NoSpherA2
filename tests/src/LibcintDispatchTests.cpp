#include "pch.h"

#include "core/convenience.h"
#include "core/constants.h"
#include "core/atoms.h"
#include "core/wfn_class.h"
#include "core/integration_params.h"
#include "core/libCintKernels.h"
#include "core/libCintMain.h"
#include "core/NoSpherA2.h"

// free functions of libCintMain.cpp / integration_params.cpp that have no header declaration
void calc_screend_functions_and_max_ij(const std::vector<atom> &atoms, const ivec &aoloc, const ivec &bas_orbital_indices, bvec2 &screened, int &max_ij);
ivec generate_bas_indices_per_atom(const Int_Params &params);
double NOS_CINTcommon_fac_sp(int l);

namespace LibcintDispatchTestHelpers
{
    // one primitive per shell: (l, exponent), contraction coefficient 1 before normalisation
    typedef std::vector<std::pair<int, double>> shell_list;

    atom make_atom(const std::string &label, int Z, double x, double y, double z, const shell_list &shells, int type_offset = 0)
    {
        atom a(label, {}, 1, x, y, z, Z);
        for (int s = 0; s < (int)shells.size(); s++)
            a.push_back_basis_set(shells[s].second, 1.0, shells[s].first + type_offset, s);
        return a;
    }

    WFN make_wfn(const std::vector<atom> &atoms, e_origin origin = e_origin::NOT_YET_DEFINED)
    {
        WFN w(origin);
        for (const atom &a : atoms) w.push_back_atom(a);
        w.set_origin(origin);
        return w;
    }

    // int_0^inf r^n exp(-e r^2) dr, the same expression Int_Params::gaussian_int uses
    double gauss_int(int n, double e)
    {
        const double n1 = (n + 1) * 0.5;
        return std::tgamma(n1) / (2.0 * std::pow(e, n1));
    }

    // radial normalisation of a single primitive r^l exp(-a r^2) times a unit-normalised real solid harmonic
    double radial_norm(int l, double a)
    {
        return 1.0 / std::sqrt(gauss_int(2 * l + 2, 2.0 * a));
    }

    double s_norm(double a)
    {
        return std::pow(2.0 * a / constants::PI, 0.75);
    }

    double boys0(double T)
    {
        return T < 1e-14 ? 1.0 : 0.5 * std::sqrt(constants::PI / T) * std::erf(std::sqrt(T));
    }

    // overlap of two unit-normalised s primitives a distance sqrt(R2) apart
    double s_overlap(double a, double b, double R2)
    {
        return std::pow(2.0 * std::sqrt(a * b) / (a + b), 1.5) * std::exp(-a * b * R2 / (a + b));
    }

    // two-centre Coulomb integral (a|b) of two unit-normalised s primitives
    double s_coulomb2c(double a, double b, double R2)
    {
        return s_norm(a) * s_norm(b) * 2.0 * std::pow(constants::PI, 2.5) / (a * b * std::sqrt(a + b)) * boys0(a * b * R2 / (a + b));
    }

    // block-diagonal cart -> sph matrix of a 0-based-type wavefunction, one cart2sph(l, true) block per shell
    dMatrix2 block_cart2sph(const WFN &w0)
    {
        ivec ls;
        for (const atom &a : w0.get_atoms())
            for (int shell = 0, prim = 0; shell < a.get_shellcount_size(); prim += a.get_shellcount(shell++))
                ls.push_back(a.get_basis_set_type(prim));
        int nc = 0, ns = 0;
        for (const int l : ls) nc += (l + 1) * (l + 2) / 2, ns += 2 * l + 1;
        dMatrix2 C(nc, ns);
        int ic = 0, is = 0;
        for (const int l : ls)
        {
            const dMatrix2 b = cart2sph(l, true);
            for (int i = 0; i < (int)b.extent(0); i++)
                for (int j = 0; j < (int)b.extent(1); j++)
                    C(ic + i, is + j) = b(i, j);
            ic += (int)b.extent(0);
            is += (int)b.extent(1);
        }
        return C;
    }

    // C^T M C for a row-major n_cart x n_cart matrix M and a n_cart x n_sph transformation C
    vec transform(const vec &M, const dMatrix2 &C)
    {
        const int nc = (int)C.extent(0), ns = (int)C.extent(1);
        vec out((size_t)ns * ns, 0.0);
        for (int p = 0; p < ns; p++)
            for (int q = 0; q < ns; q++)
            {
                double acc = 0.0;
                for (int i = 0; i < nc; i++)
                    for (int j = 0; j < nc; j++)
                        acc += C(i, p) * M[(size_t)i * nc + j] * C(j, q);
                out[(size_t)p * ns + q] = acc;
            }
        return out;
    }

    // rho_k = sum_ij (ij|k) dm_ij from a C-ordered eri3c (i * naoj * naok + j * naok + k)
    vec contract_eri(const vec &eri, const dMatrix2 &dm, int nao, int naok)
    {
        vec rho(naok, 0.0);
        for (int i = 0; i < nao; i++)
            for (int j = 0; j < nao; j++)
                for (int k = 0; k < naok; k++)
                    rho[k] += eri[((size_t)i * nao + j) * naok + k] * dm(i, j);
        return rho;
    }

    dMatrix2 symmetric_dm(int n)
    {
        dMatrix2 dm(n, n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                dm(i, j) = 0.05 * (i + j + 1) + (i == j ? 1.0 : 0.0);
        return dm;
    }

    std::string read_file(const std::filesystem::path &p)
    {
        std::ifstream f(p);
        std::stringstream ss;
        ss << f.rdbuf();
        return ss.str();
    }
}

using namespace LibcintDispatchTestHelpers;

// the default constructor is the empty state computeEri3c and computeRho start their combined object from
TEST(LibcintDispatchTests, DefaultConstructedParamsAreEmpty)
{
    Int_Params p;
    EXPECT_EQ(p.get_natoms(), 0u);
    EXPECT_EQ(p.get_nbas(), 0u);
    EXPECT_EQ(p.get_nao(), 0u);
    EXPECT_TRUE(p.get_atm().empty());
    EXPECT_TRUE(p.get_bas().empty());
    EXPECT_TRUE(p.get_env().empty());
    EXPECT_TRUE(p.get_basis_sets().empty());
}

// the libcint atm/bas/env layout for two s atoms, hand-derived from the libcint slot conventions
TEST(LibcintDispatchTests, TwoAtomLayoutOfAtmBasEnv)
{
    WFN w = make_wfn({ make_atom("H", 1, 0.0, 0.0, 0.0, { {0, 1.2} }), make_atom("He", 2, 0.0, 0.0, 1.5, { {0, 0.7} }) });
    Int_Params p(w);
    EXPECT_EQ(p.get_natoms(), 2u);
    EXPECT_EQ(p.get_nbas(), 2u);
    EXPECT_EQ(p.get_nao(), 2u);

    const ivec expected_atm = { 1, 20, 1, 23, 0, 0, 2, 24, 1, 27, 0, 0 };
    EXPECT_EQ(p.get_atm(), expected_atm);
    const ivec expected_bas = { 0, 0, 1, 1, 0, 28, 29, 0, 1, 0, 1, 1, 0, 30, 31, 0 };
    EXPECT_EQ(p.get_bas(), expected_bas);

    const vec env = p.get_env();
    ASSERT_EQ(env.size(), 32u);
    for (int i = 0; i < 20; i++) EXPECT_EQ(env[i], 0.0);
    EXPECT_EQ(env[22], 0.0);
    EXPECT_EQ(env[26], 1.5);
    EXPECT_EQ(env[28], 1.2);
    EXPECT_EQ(env[30], 0.7);
    // a single primitive is normalised to 1 / sqrt(int r^2 exp(-2 a r^2) dr)
    EXPECT_NEAR(env[29], 1.0 / std::sqrt(gauss_int(2, 2.4)), 1e-12);
    EXPECT_NEAR(env[31], 1.0 / std::sqrt(gauss_int(2, 1.4)), 1e-12);
    EXPECT_NEAR(Int_Params::normalize_gto({ 1.0 }, { 1.2 }, 0)[0], env[29], 1e-12);

    std::map<int, LibCintBasis> sets = p.get_basis_sets();
    ASSERT_EQ(sets.size(), 2u);
    EXPECT_EQ(sets[1].env_idx, 28);
    EXPECT_EQ(sets[2].env_idx, 30);
    EXPECT_EQ(sets[2].shelltypes, (ivec{ 0 }));
}

// the combined object shifts the second basis' atom index and env pointers behind the first one
TEST(LibcintDispatchTests, CombinedParamsOffsetTheSecondBasis)
{
    WFN a = make_wfn({ make_atom("H", 1, 0.0, 0.0, 0.0, { {0, 1.2}, {1, 0.8} }) });
    WFN b = make_wfn({ make_atom("He", 2, 0.0, 0.0, 1.5, { {0, 0.7} }) });
    Int_Params pa(a), pb(b);
    Int_Params c(pa, pb);
    EXPECT_EQ(c.get_natoms(), 2u);
    EXPECT_EQ(c.get_nbas(), 3u);
    EXPECT_EQ(c.get_nao(), 5u);
    const int off = (int)pa.get_env().size();
    ivec atm = c.get_atm(), bas = c.get_bas();
    ASSERT_EQ(atm.size(), 12u);
    ASSERT_EQ(bas.size(), 24u);
    // the first basis is copied verbatim
    EXPECT_EQ(ivec(atm.begin(), atm.begin() + 6), pa.get_atm());
    EXPECT_EQ(ivec(bas.begin(), bas.begin() + 16), pa.get_bas());
    EXPECT_EQ(atm[6 + 0], 2);
    EXPECT_EQ(atm[6 + 1], 20 + off);
    EXPECT_EQ(atm[6 + 3], 23 + off);
    EXPECT_EQ(bas[16 + 0], 1);
    EXPECT_EQ(bas[16 + 5], pb.get_bas()[5] + off);
    EXPECT_EQ(bas[16 + 6], pb.get_bas()[6] + off);
    const vec env = c.get_env();
    EXPECT_EQ(vec(env.begin(), env.begin() + off), pa.get_env());
    EXPECT_EQ(env[bas[16 + 5]], 0.7);
    EXPECT_EQ(env[20 + off + 2], 1.5);
}

// spherical and cartesian overlap of two normalised s primitives against the closed-form Gaussian overlap
TEST(LibcintDispatchTests, OverlapSphAndCartMatchTheAnalyticValue)
{
    const double a = 1.2, b = 0.7, R = 1.5;
    WFN w = make_wfn({ make_atom("H", 1, 0.0, 0.0, 0.0, { {0, a} }), make_atom("He", 2, 0.0, 0.0, R, { {0, b} }) });
    Int_Params p(w);
    vec sph, crt;
    compute2C<Overlap2C_SPH>(p, sph);
    compute2C<Overlap2C_CRT>(p, crt);
    ASSERT_EQ(sph.size(), 4u);
    ASSERT_EQ(crt.size(), 4u);
    const double S = s_overlap(a, b, R * R);
    EXPECT_NEAR(sph[0], 1.0, 1e-12);
    EXPECT_NEAR(sph[3], 1.0, 1e-12);
    EXPECT_NEAR(sph[1], S, 1e-12);
    EXPECT_NEAR(sph[2], S, 1e-12);
    for (int i = 0; i < 4; i++) EXPECT_NEAR(crt[i], sph[i], 1e-14);
}

// the two-centre Coulomb metric of two normalised s primitives against the Boys-function formula
TEST(LibcintDispatchTests, Coulomb2cSphAndCartMatchTheAnalyticValue)
{
    const double a = 1.2, b = 0.7, R = 1.5;
    WFN w = make_wfn({ make_atom("H", 1, 0.0, 0.0, 0.0, { {0, a} }), make_atom("He", 2, 0.0, 0.0, R, { {0, b} }) });
    Int_Params p(w);
    vec sph, crt;
    compute2C<Coulomb2C_SPH>(p, sph);
    compute2C<Coulomb2C_CRT>(p, crt);
    ASSERT_EQ(sph.size(), 4u);
    EXPECT_NEAR(sph[0], s_coulomb2c(a, a, 0.0), 1e-10);
    EXPECT_NEAR(sph[3], s_coulomb2c(b, b, 0.0), 1e-10);
    EXPECT_NEAR(sph[1], s_coulomb2c(a, b, R * R), 1e-10);
    EXPECT_NEAR(sph[2], sph[1], 1e-12);
    for (int i = 0; i < 4; i++) EXPECT_NEAR(crt[i], sph[i], 1e-12);
}

// a two-primitive contraction is renormalised as a whole, so its self overlap is 1 and its overlap
// with the second atom is the coefficient-weighted sum of primitive overlaps
TEST(LibcintDispatchTests, ContractedShellIsNormalisedAsAWhole)
{
    atom h("H", {}, 1, 0.0, 0.0, 0.0, 1);
    h.push_back_basis_set(2.0, 0.4, 0, 0);
    h.push_back_basis_set(0.5, 0.7, 0, 0);
    WFN w = make_wfn({ h, make_atom("He", 2, 0.0, 0.0, 1.1, { {0, 0.9} }) });
    Int_Params p(w);
    EXPECT_EQ(p.get_nbas(), 2u);
    vec S;
    compute2C<Overlap2C_SPH>(p, S);
    ASSERT_EQ(S.size(), 4u);
    EXPECT_NEAR(S[0], 1.0, 1e-12);
    EXPECT_NEAR(S[3], 1.0, 1e-12);

    const vec c = Int_Params::normalize_gto({ 0.4, 0.7 }, { 2.0, 0.5 }, 0);
    const vec env = p.get_env();
    EXPECT_NEAR(env[30], c[0], 1e-12);
    EXPECT_NEAR(env[31], c[1], 1e-12);
    // primitive-normalised coefficients relate to the libcint ones by 1 / sqrt(gaussian_int)
    const double d0 = c[0] * std::sqrt(gauss_int(2, 4.0)), d1 = c[1] * std::sqrt(gauss_int(2, 1.0));
    EXPECT_NEAR(d0 * d0 + d1 * d1 + 2.0 * d0 * d1 * s_overlap(2.0, 0.5, 0.0), 1.0, 1e-12);
    EXPECT_NEAR(S[1], d0 * s_overlap(2.0, 0.9, 1.21) + d1 * s_overlap(0.5, 0.9, 1.21), 1e-12);
}

// make_loc counts 2l+1 spherical and (l+1)(l+2)/2 cartesian functions per shell
TEST(LibcintDispatchTests, MakeLocCountsSphericalAndCartesianFunctions)
{
    WFN w = make_wfn({ make_atom("H", 1, 0.0, 0.0, 0.0, { {0, 1.0}, {1, 0.8}, {2, 0.6}, {3, 0.5} }) });
    Int_Params p(w);
    ivec bas = p.get_bas();
    const int nbas = (int)p.get_nbas();
    ASSERT_EQ(nbas, 4);
    EXPECT_EQ(make_loc<COORDINATE_TYPE::SPH>(bas, nbas), (ivec{ 0, 1, 4, 9, 16 }));
    EXPECT_EQ(make_loc<COORDINATE_TYPE::CART>(bas, nbas), (ivec{ 0, 1, 4, 10, 20 }));
    EXPECT_EQ(Overlap2C_SPH::gen_loc(bas, nbas), (ivec{ 0, 1, 4, 9, 16 }));
    EXPECT_EQ(Overlap2C_CRT::gen_loc(bas, nbas), (ivec{ 0, 1, 4, 10, 20 }));
    EXPECT_EQ(Coulomb3C_CRT::gen_loc(bas, nbas), (ivec{ 0, 1, 4, 10, 20 }));
    EXPECT_EQ(Overlap3C_SPH::gen_loc(bas, nbas), (ivec{ 0, 1, 4, 9, 16 }));
    EXPECT_EQ(p.get_nao(), 16u);
}

// libcint's common factor for s and p is the Y00 / Y1m angular normalisation, 1 for higher l
TEST(LibcintDispatchTests, CommonFactorMatchesAngularNormalisation)
{
    EXPECT_NEAR(NOS_CINTcommon_fac_sp(0), std::sqrt(1.0 / (4.0 * constants::PI)), 1e-15);
    EXPECT_NEAR(NOS_CINTcommon_fac_sp(1), std::sqrt(3.0 / (4.0 * constants::PI)), 1e-15);
    EXPECT_EQ(NOS_CINTcommon_fac_sp(2), 1.0);
    EXPECT_EQ(NOS_CINTcommon_fac_sp(7), 1.0);
}

// s is the identity; p is the PYPZPX permutation (cart x, y, z -> sph y, z, x) of the libcint build,
// both scaled by the common factor when 'unnormalised' is requested
TEST(LibcintDispatchTests, Cart2SphLowLIsIdentityOrPermutation)
{
    dMatrix2 s = cart2sph(0, true), p = cart2sph(1, true);
    ASSERT_EQ(s.extent(0), 1u);
    ASSERT_EQ(p.extent(0), 3u);
    ASSERT_EQ(p.extent(1), 3u);
    EXPECT_EQ(s(0, 0), 1.0);
    const auto p_perm = [](int i, int j) { return (i == 1 && j == 0) || (i == 2 && j == 1) || (i == 0 && j == 2); };
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            EXPECT_EQ(p(i, j), p_perm(i, j) ? 1.0 : 0.0) << i << "," << j;

    dMatrix2 su = cart2sph(0, false), pu = cart2sph(1, false);
    EXPECT_NEAR(su(0, 0), 0.282094791773878143, 1e-15);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            EXPECT_NEAR(pu(i, j), p_perm(i, j) ? 0.488602511902919921 : 0.0, 1e-15) << i << "," << j;
}

// the d block carries libcint's real solid harmonic coefficients in cart order xx,xy,xz,yy,yz,zz and m = -2..2
TEST(LibcintDispatchTests, Cart2SphDBlockHoldsSolidHarmonicCoefficients)
{
    dMatrix2 d = cart2sph(2, true);
    ASSERT_EQ(d.extent(0), 6u);
    ASSERT_EQ(d.extent(1), 5u);
    const double c22 = std::sqrt(15.0 / (4.0 * constants::PI));
    const double c20 = std::sqrt(5.0 / (16.0 * constants::PI));
    const double c2p2 = std::sqrt(15.0 / (16.0 * constants::PI));
    // xy -> m=-2, yz -> m=-1, (2zz - xx - yy) -> m=0, xz -> m=1, (xx - yy) -> m=2; every other entry is 0
    const double expected_d[6][5] = {
        { 0.0, 0.0, -c20, 0.0, c2p2 },
        { c22, 0.0, 0.0, 0.0, 0.0 },
        { 0.0, 0.0, 0.0, c22, 0.0 },
        { 0.0, 0.0, -c20, 0.0, -c2p2 },
        { 0.0, c22, 0.0, 0.0, 0.0 },
        { 0.0, 0.0, 2.0 * c20, 0.0, 0.0 } };
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 5; j++)
            EXPECT_NEAR(d(i, j), expected_d[i][j], 1e-14) << "entry " << i << "," << j;
    // every real solid harmonic column is unit normalised over the sphere: sum of squares
    // weighted by the cartesian angular integrals equals 1 for the m = 0 column
    double n0 = 0.0;
    const double xx4 = 4.0 * constants::PI / 5.0, xxyy = 4.0 * constants::PI / 15.0;
    n0 = d(0, 2) * d(0, 2) * xx4 + d(3, 2) * d(3, 2) * xx4 + d(5, 2) * d(5, 2) * xx4
        + 2.0 * (d(0, 2) * d(3, 2) + d(0, 2) * d(5, 2) + d(3, 2) * d(5, 2)) * xxyy;
    EXPECT_NEAR(n0, 1.0, 1e-13);
    // f goes through the general libcint path; its m = 0 column is sqrt(7/16pi) (2zzz - 3xxz - 3yyz)
    // in cart order xxx,xxy,xxz,xyy,xyz,xzz,yyy,yyz,yzz,zzz
    dMatrix2 f = cart2sph(3, false);
    ASSERT_EQ(f.extent(0), 10u);
    ASSERT_EQ(f.extent(1), 7u);
    const double c30 = std::sqrt(7.0 / (16.0 * constants::PI));
    for (int i = 0; i < 10; i++)
    {
        const double expected = i == 9 ? 2.0 * c30 : (i == 2 || i == 7) ? -3.0 * c30 : 0.0;
        EXPECT_NEAR(f(i, 3), expected, 1e-14) << "cart " << i;
    }
}

// the wavefunction-wide matrix is block diagonal in shell order with 1-based types on the atoms
TEST(LibcintDispatchTests, WavefunctionCart2SphMatrixIsBlockDiagonal)
{
    // get_cart2sph_matrix reads type - 1 as l, so build the shells with 1-based types
    WFN w1 = make_wfn({ make_atom("H", 1, 0.0, 0.0, 0.0, { {0, 1.0}, {1, 0.8} }, 1), make_atom("He", 2, 0.0, 0.0, 1.0, { {2, 0.6} }, 1) });
    dMatrix2 C = get_cart2sph_matrix(w1, true);
    ASSERT_EQ(C.extent(0), 10u);
    ASSERT_EQ(C.extent(1), 9u);
    const dMatrix2 p = cart2sph(1, true), d = cart2sph(2, true);
    for (int i = 0; i < 10; i++)
        for (int j = 0; j < 9; j++)
        {
            double expected = 0.0;
            if (i == 0 && j == 0) expected = 1.0;
            else if (i >= 1 && i < 4 && j >= 1 && j < 4) expected = p(i - 1, j - 1);
            else if (i >= 4 && j >= 4) expected = d(i - 4, j - 4);
            EXPECT_NEAR(C(i, j), expected, 1e-14) << "entry " << i << "," << j;
        }
}

// a shell behind a d shell lands at the spherical column offset, not the cartesian one
TEST(LibcintDispatchTests, WavefunctionCart2SphMatrixOffsetsSphericalColumns)
{
    WFN w0 = make_wfn({ make_atom("H", 1, 0.0, 0.0, 0.0, { {2, 0.9}, {0, 1.1} }) });
    WFN w1 = make_wfn({ make_atom("H", 1, 0.0, 0.0, 0.0, { {2, 0.9}, {0, 1.1} }, 1) });
    const dMatrix2 C = get_cart2sph_matrix(w1, true), expected = block_cart2sph(w0);
    ASSERT_EQ(C.extent(0), 7u);
    ASSERT_EQ(C.extent(1), 6u);
    for (int i = 0; i < 7; i++)
        for (int j = 0; j < 6; j++)
            EXPECT_NEAR(C(i, j), expected(i, j), 1e-14) << "entry " << i << "," << j;
    const dMatrix2 p = cart2sph(1, true);
    EXPECT_EQ(p(1, 0), 1.0);
    EXPECT_EQ(p(2, 1), 1.0);
    EXPECT_EQ(p(0, 2), 1.0);
}

// the spherical overlap must be the cartesian overlap transformed with the same coefficients libcint uses
TEST(LibcintDispatchTests, SphericalOverlapIsTransformedCartesianOverlap)
{
    const shell_list h_shells = { {0, 1.1}, {1, 0.7}, {2, 0.9}, {3, 0.5} }, li_shells = { {0, 0.6}, {1, 0.4} };
    WFN w0 = make_wfn({ make_atom("H", 1, 0.0, 0.0, 0.0, h_shells), make_atom("Li", 3, 0.3, -0.4, 1.2, li_shells) });
    Int_Params p(w0);
    vec sph, crt;
    compute2C<Overlap2C_SPH>(p, sph);
    compute2C<Overlap2C_CRT>(p, crt);
    ASSERT_EQ(sph.size(), 20u * 20u);
    ASSERT_EQ(crt.size(), 24u * 24u);
    const vec transformed = transform(crt, block_cart2sph(w0));
    for (size_t i = 0; i < sph.size(); i++)
        EXPECT_NEAR(transformed[i], sph[i], 1e-12) << "element " << i;
    for (int i = 0; i < 20; i++) EXPECT_NEAR(sph[(size_t)i * 21], 1.0, 1e-12);
}

// the same identity for the Coulomb metric exercises the cartesian 2c2e kernel
TEST(LibcintDispatchTests, SphericalCoulombIsTransformedCartesianCoulomb)
{
    const shell_list h_shells = { {0, 1.1}, {1, 0.7}, {2, 0.9} }, li_shells = { {0, 0.6}, {1, 0.4} };
    WFN w0 = make_wfn({ make_atom("H", 1, 0.0, 0.0, 0.0, h_shells), make_atom("Li", 3, 0.3, -0.4, 1.2, li_shells) });
    Int_Params p(w0);
    vec sph, crt;
    compute2C<Coulomb2C_SPH>(p, sph);
    compute2C<Coulomb2C_CRT>(p, crt);
    ASSERT_EQ(sph.size(), 13u * 13u);
    ASSERT_EQ(crt.size(), 14u * 14u);
    const vec transformed = transform(crt, block_cart2sph(w0));
    for (size_t i = 0; i < sph.size(); i++)
        EXPECT_NEAR(transformed[i], sph[i], 1e-10) << "element " << i;
    for (int i = 0; i < 13; i++) EXPECT_GT(sph[(size_t)i * 14], 0.0);
}

// three-centre overlaps of same-centre s primitives are (pi / (a + b + c))^(3/2) times the norms; the
// result is left in Fortran order i + j * naoi + k * naoi * naoj
TEST(LibcintDispatchTests, Compute3cOverlapMatchesSameCentreFormula)
{
    const double e[2] = { 1.1, 0.4 }, f[2] = { 0.8, 2.0 };
    WFN qm = make_wfn({ make_atom("H", 1, 0.0, 0.0, 0.0, { {0, e[0]}, {0, e[1]} }) });
    WFN aux = make_wfn({ make_atom("He", 2, 0.0, 0.0, 0.0, { {0, f[0]}, {0, f[1]} }) });
    Int_Params pq(qm), pa(aux);
    vec res;
    compute3C<Overlap3C_SPH>(pq, pa, res);
    ASSERT_EQ(res.size(), 8u);
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
            {
                const double expected = s_norm(e[i]) * s_norm(e[j]) * s_norm(f[k]) * std::pow(constants::PI / (e[i] + e[j] + f[k]), 1.5);
                EXPECT_NEAR(res[i + j * 2 + k * 4], expected, 1e-12) << i << j << k;
            }
}

// the three-centre Coulomb integral of same-centre s primitives has the closed form 2 pi^(5/2) / (p q sqrt(p + q))
TEST(LibcintDispatchTests, Compute3cCoulombMatchesSameCentreFormula)
{
    const double e[2] = { 1.1, 0.4 }, f[2] = { 0.8, 2.0 };
    WFN qm = make_wfn({ make_atom("H", 1, 0.0, 0.0, 0.0, { {0, e[0]}, {0, e[1]} }) });
    WFN aux = make_wfn({ make_atom("He", 2, 0.0, 0.0, 0.0, { {0, f[0]}, {0, f[1]} }) });
    Int_Params pq(qm), pa(aux);
    vec res;
    compute3C<Coulomb3C_SPH>(pq, pa, res);
    ASSERT_EQ(res.size(), 8u);
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
            {
                const double p = e[i] + e[j], q = f[k];
                const double expected = s_norm(e[i]) * s_norm(e[j]) * s_norm(f[k]) * 2.0 * std::pow(constants::PI, 2.5) / (p * q * std::sqrt(p + q));
                EXPECT_NEAR(res[i + j * 2 + k * 4], expected, 1e-10) << i << j << k;
            }
}

// computeEri3c reorders the same libcint result into C order i * naoj * naok + j * naok + k
TEST(LibcintDispatchTests, Eri3cIsTheCOrderedCompute3cResult)
{
    WFN qm = make_wfn({ make_atom("H", 1, 0.0, 0.0, 0.0, { {0, 1.0}, {1, 0.6} }), make_atom("Li", 3, 0.0, 0.5, 1.4, { {0, 0.5} }) });
    WFN aux = make_wfn({ make_atom("He", 2, 0.0, 0.0, 0.0, { {0, 0.9}, {1, 0.7} }), make_atom("Be", 4, 0.0, 0.5, 1.4, { {0, 1.3}, {2, 0.8} }) });
    Int_Params pq(qm), pa(aux);
    const int nao = (int)pq.get_nao(), naok = (int)pa.get_nao();
    ASSERT_EQ(nao, 5);
    ASSERT_EQ(naok, 10);
    vec f_order, c_order;
    compute3C<Coulomb3C_SPH>(pq, pa, f_order);
    computeEri3c<Coulomb3C_SPH>(pq, pa, c_order);
    ASSERT_EQ(f_order.size(), (size_t)nao * nao * naok);
    ASSERT_EQ(c_order.size(), f_order.size());
    double max_abs = 0.0;
    for (int i = 0; i < nao; i++)
        for (int j = 0; j < nao; j++)
            for (int k = 0; k < naok; k++)
            {
                const double v = f_order[i + j * nao + k * nao * nao];
                EXPECT_EQ(c_order[((size_t)i * nao + j) * naok + k], v);
                EXPECT_NEAR(f_order[j + i * nao + k * nao * nao], v, 1e-12);
                max_abs = std::max(max_abs, std::abs(v));
            }
    EXPECT_GT(max_abs, 0.1);
}

// computeRho over atom-pair blocks with weight 2 for i != j equals the plain contraction sum_ij (ij|k) dm_ij
TEST(LibcintDispatchTests, ComputeRhoCoulombEqualsContractedEri3c)
{
    WFN qm = make_wfn({ make_atom("H", 1, 0.0, 0.0, 0.0, { {0, 1.0}, {1, 0.6} }), make_atom("Li", 3, 0.0, 0.5, 1.4, { {0, 0.5} }) });
    WFN aux = make_wfn({ make_atom("He", 2, 0.0, 0.0, 0.0, { {0, 0.9}, {1, 0.7} }), make_atom("Be", 4, 0.0, 0.5, 1.4, { {0, 1.3}, {2, 0.8} }) });
    Int_Params pq(qm), pa(aux);
    const int nao = (int)pq.get_nao(), naok = (int)pa.get_nao();
    dMatrix2 dm = symmetric_dm(nao);
    vec eri, rho;
    computeEri3c<Coulomb3C_SPH>(pq, pa, eri);
    computeRho<Coulomb3C_SPH>(pq, pa, dm, rho);
    ASSERT_EQ((int)rho.size(), naok);
    const vec expected = contract_eri(eri, dm, nao, naok);
    for (int k = 0; k < naok; k++) EXPECT_NEAR(rho[k], expected[k], 1e-10) << "k = " << k;
    EXPECT_GT(std::abs(rho[0]), 0.1);
}

// the same identity with the overlap metric, whose 3c result is only available in Fortran order
TEST(LibcintDispatchTests, ComputeRhoOverlapEqualsContractedCompute3c)
{
    WFN qm = make_wfn({ make_atom("H", 1, 0.0, 0.0, 0.0, { {0, 1.0}, {1, 0.6} }), make_atom("Li", 3, 0.0, 0.5, 1.4, { {0, 0.5} }) });
    WFN aux = make_wfn({ make_atom("He", 2, 0.0, 0.0, 0.0, { {0, 0.9}, {1, 0.7} }), make_atom("Be", 4, 0.0, 0.5, 1.4, { {0, 1.3} }) });
    Int_Params pq(qm), pa(aux);
    const int nao = (int)pq.get_nao(), naok = (int)pa.get_nao();
    dMatrix2 dm = symmetric_dm(nao);
    vec ovl, rho;
    compute3C<Overlap3C_SPH>(pq, pa, ovl);
    computeRho<Overlap3C_SPH>(pq, pa, dm, rho);
    ASSERT_EQ((int)rho.size(), naok);
    for (int k = 0; k < naok; k++)
    {
        double expected = 0.0;
        for (int i = 0; i < nao; i++)
            for (int j = 0; j < nao; j++)
                expected += ovl[i + j * nao + k * nao * nao] * dm(i, j);
        EXPECT_NEAR(rho[k], expected, 1e-10) << "k = " << k;
    }
    EXPECT_GT(std::abs(rho[0]), 0.1);
}

// with only s and p shells the cartesian and spherical functions coincide, so the cartesian kernel must reproduce the spherical rho
TEST(LibcintDispatchTests, ComputeRhoCartesianMatchesSphericalForSpBasis)
{
    WFN qm = make_wfn({ make_atom("H", 1, 0.0, 0.0, 0.0, { {0, 1.0}, {1, 0.6} }), make_atom("Li", 3, 0.0, 0.5, 1.4, { {0, 0.5} }) });
    WFN aux = make_wfn({ make_atom("He", 2, 0.0, 0.0, 0.0, { {0, 0.9}, {1, 0.7} }), make_atom("Be", 4, 0.0, 0.5, 1.4, { {0, 1.3}, {1, 0.4} }) });
    Int_Params pq(qm), pa(aux);
    const int nao = (int)pq.get_nao();
    dMatrix2 dm = symmetric_dm(nao);
    // for s and p the cart -> sph matrices are permutations (PYPZPX): dm_crt = Cq dm Cq^T, rho_sph = Ca^T rho_crt
    const dMatrix2 Cq = block_cart2sph(qm), Ca = block_cart2sph(aux);
    dMatrix2 dm_crt(nao, nao);
    for (int i = 0; i < nao; i++)
        for (int j = 0; j < nao; j++)
            for (int p = 0; p < nao; p++)
                for (int q = 0; q < nao; q++)
                    dm_crt(i, j) += Cq(i, p) * dm(p, q) * Cq(j, q);
    vec sph, crt;
    computeRho<Coulomb3C_SPH>(pq, pa, dm, sph);
    computeRho<Coulomb3C_CRT>(pq, pa, dm_crt, crt);
    ASSERT_EQ(sph.size(), 8u);
    ASSERT_EQ(crt.size(), sph.size());
    for (size_t k = 0; k < sph.size(); k++)
    {
        double from_crt = 0.0;
        for (size_t c = 0; c < crt.size(); c++) from_crt += Ca(c, k) * crt[c];
        EXPECT_NEAR(from_crt, sph[k], 1e-10) << "k = " << k;
    }
}

// the overlap screening flags a pair 60 bohr apart, keeps the diagonal, reports the largest block and clears on no atoms
TEST(LibcintDispatchTests, ScreeningFlagsFarPairsAndHandlesNoAtoms)
{
    WFN qm = make_wfn({ make_atom("H", 1, 0.0, 0.0, 0.0, { {0, 1.0}, {1, 0.6} }), make_atom("Li", 3, 0.0, 0.0, 60.0, { {0, 0.8} }) });
    Int_Params p(qm);
    ivec bas = p.get_bas();
    ivec aoloc = make_loc<COORDINATE_TYPE::SPH>(bas, (int)p.get_nbas());
    ivec per_atom = generate_bas_indices_per_atom(p);
    EXPECT_EQ(per_atom, (ivec{ 0, 2, 3 }));
    bvec2 screened;
    int max_ij = -1;
    calc_screend_functions_and_max_ij(p.get_atoms(), aoloc, per_atom, screened, max_ij);
    ASSERT_EQ(screened.size(), 2u);
    EXPECT_FALSE(screened[0][0]);
    EXPECT_TRUE(screened[0][1]);
    EXPECT_FALSE(screened[1][1]);
    EXPECT_EQ(max_ij, 16);

    calc_screend_functions_and_max_ij({}, aoloc, per_atom, screened, max_ij);
    EXPECT_TRUE(screened.empty());
    EXPECT_EQ(max_ij, 0);
}

// a screened pair contributes exp(-1600) to rho, so the blocked contraction still matches the full one;
// this only catches a pair screened wrongly (a lost contribution), a screen that never fires is
// invisible here and is pinned down by ScreeningFlagsFarPairsAndHandlesNoAtoms instead
TEST(LibcintDispatchTests, ComputeRhoWithScreenedPairMatchesFullContraction)
{
    WFN qm = make_wfn({ make_atom("H", 1, 0.0, 0.0, 0.0, { {0, 1.0}, {1, 0.6} }), make_atom("Li", 3, 0.0, 0.0, 60.0, { {0, 0.8} }) });
    WFN aux = make_wfn({ make_atom("He", 2, 0.0, 0.0, 0.0, { {0, 0.9} }), make_atom("Be", 4, 0.0, 0.0, 60.0, { {0, 1.3}, {1, 0.5} }) });
    Int_Params pq(qm), pa(aux);
    const int nao = (int)pq.get_nao(), naok = (int)pa.get_nao();
    dMatrix2 dm = symmetric_dm(nao);
    vec eri, rho;
    computeEri3c<Coulomb3C_SPH>(pq, pa, eri);
    computeRho<Coulomb3C_SPH>(pq, pa, dm, rho);
    const vec expected = contract_eri(eri, dm, nao, naok);
    ASSERT_EQ((int)rho.size(), naok);
    for (int k = 0; k < naok; k++) EXPECT_NEAR(rho[k], expected[k], 1e-10) << "k = " << k;
}

// AO values on three points against normalised real solid harmonics times the Gaussian radial part
TEST(LibcintDispatchTests, EvalGtoSphMatchesSolidHarmonicsOnThreePoints)
{
    const double as = 0.9, ap = 0.7, ad = 0.5, af = 0.4, ag = 0.3;
    WFN w = make_wfn({ make_atom("H", 1, 0.0, 0.0, 0.0, { {0, as}, {1, ap}, {2, ad}, {3, af}, {4, ag} }) });
    Int_Params p(w);
    const int nao = (int)p.get_nao();
    ASSERT_EQ(nao, 25);
    vec2 grid = { { 0.3, 0.0, 0.6 }, { -0.2, 0.0, 0.0 }, { 0.5, 0.8, 0.0 } };
    ivec slice;
    vec v = eval_GTO_sph(p, grid, slice);
    EXPECT_EQ(slice, (ivec{ 0, 5 }));
    ASSERT_EQ(v.size(), (size_t)nao * 3);
    const double pi4 = 4.0 * constants::PI;

    // point 0: s and the three p functions in libcint order x, y, z
    {
        const double x = 0.3, y = -0.2, z = 0.5, r2 = x * x + y * y + z * z;
        EXPECT_NEAR(v[0], s_norm(as) * std::exp(-as * r2), 1e-12);
        const double np = radial_norm(1, ap) * std::sqrt(3.0 / pi4) * std::exp(-ap * r2);
        EXPECT_NEAR(v[1], np * x, 1e-12);
        EXPECT_NEAR(v[2], np * y, 1e-12);
        EXPECT_NEAR(v[3], np * z, 1e-12);
    }
    // point 1 on the z axis: only m = 0 of d, f and g survive
    {
        const double z = 0.8, z2 = z * z;
        const double *row = v.data() + nao;
        EXPECT_NEAR(row[0], s_norm(as) * std::exp(-as * z2), 1e-12);
        EXPECT_NEAR(row[3], radial_norm(1, ap) * std::sqrt(3.0 / pi4) * z * std::exp(-ap * z2), 1e-12);
        EXPECT_NEAR(row[4 + 2], radial_norm(2, ad) * std::sqrt(5.0 / pi4) * z2 * std::exp(-ad * z2), 1e-12);
        EXPECT_NEAR(row[9 + 3], radial_norm(3, af) * std::sqrt(7.0 / pi4) * z2 * z * std::exp(-af * z2), 1e-12);
        EXPECT_NEAR(row[16 + 4], radial_norm(4, ag) * std::sqrt(9.0 / pi4) * z2 * z2 * std::exp(-ag * z2), 1e-12);
        for (int m : { 4, 5, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18, 19, 21, 22, 23, 24 })
            EXPECT_NEAR(row[m], 0.0, 1e-14) << "function " << m;
    }
    // point 2 on the x axis: d0 = -sqrt(5/16pi) x^2, d2 = sqrt(15/16pi) x^2
    {
        const double x = 0.6, x2 = x * x;
        const double *row = v.data() + 2 * nao;
        const double nd = radial_norm(2, ad) * std::exp(-ad * x2);
        EXPECT_NEAR(row[4 + 2], -nd * std::sqrt(5.0 / (16.0 * constants::PI)) * x2, 1e-12);
        EXPECT_NEAR(row[4 + 4], nd * std::sqrt(15.0 / (16.0 * constants::PI)) * x2, 1e-12);
        EXPECT_NEAR(row[4], 0.0, 1e-14);
        EXPECT_NEAR(row[5], 0.0, 1e-14);
        EXPECT_NEAR(row[7], 0.0, 1e-14);
        EXPECT_NEAR(row[2], 0.0, 1e-14);
        EXPECT_NEAR(row[3], 0.0, 1e-14);
    }
}

// the point count is the length of the coordinate rows, not the number of rows
TEST(LibcintDispatchTests, EvalGtoSphHandlesTwoGridPoints)
{
    const double as = 0.9;
    WFN w = make_wfn({ make_atom("H", 1, 0.0, 0.0, 0.0, { {0, as} }) });
    Int_Params p(w);
    vec2 grid = { { 0.3, 1.1 }, { -0.2, 0.4 }, { 0.5, -0.7 } };
    ivec slice;
    vec v = eval_GTO_sph(p, grid, slice);
    ASSERT_EQ(v.size(), 2u);
    for (int j = 0; j < 2; j++)
    {
        const double r2 = grid[0][j] * grid[0][j] + grid[1][j] * grid[1][j] + grid[2][j] * grid[2][j];
        EXPECT_NEAR(v[j], s_norm(as) * std::exp(-as * r2), 1e-12);
    }
}

// a tonto wavefunction (1-based types) gets sqrt(4 pi / (2l+1)!!) per shell plus sqrt(10 / 4 pi) on d
TEST(LibcintDispatchTests, TontoOriginRescalesCoefficients)
{
    WFN w = make_wfn({ make_atom("H", 1, 0.0, 0.0, 0.0, { {0, 1.3}, {2, 0.6} }, 1) }, e_origin::tonto);
    ASSERT_EQ(w.get_origin(), e_origin::tonto);
    Int_Params p(w);
    std::map<int, LibCintBasis> sets = p.get_basis_sets();
    ASSERT_EQ(sets.size(), 1u);
    const LibCintBasis &b = sets[1];
    ASSERT_EQ(b.coefficients.size(), 2u);
    EXPECT_EQ(b.shelltypes, (ivec{ 0, 2 }));
    EXPECT_EQ(b.exponents, (vec{ 1.3, 0.6 }));
    EXPECT_NEAR(b.coefficients[0], std::sqrt(4.0 * constants::PI), 1e-13);
    EXPECT_NEAR(b.coefficients[1], std::sqrt(4.0 * constants::PI / 15.0) * std::sqrt(10.0 / (4.0 * constants::PI)), 1e-13);
    EXPECT_EQ(p.get_nao(), 6u);
    // the s coefficient cancels libcint's 1 / (4 pi) common factor, so its self overlap is the raw Gaussian one
    vec S;
    compute2C<Overlap2C_SPH>(p, S);
    ASSERT_EQ(S.size(), 36u);
    EXPECT_NEAR(S[0], std::pow(constants::PI / 2.6, 1.5), 1e-12);
}

// a molden wavefunction stores l + 1 as the shell type: its shells are sorted by l without touching the
// primitives, so the F2 fixture (s, p on each F) comes back as shelltypes 0, 1 with the file's exponents
TEST(LibcintDispatchTests, MoldenOriginSortsShellsByAngularMomentum)
{
    const std::filesystem::path path = nos_test_repo_root() / "tests" / "molden_file" / "F2.molden";
    if (!std::filesystem::exists(path)) GTEST_SKIP() << "fixture missing: " << path;
    WFN w(e_origin::molden);
    std::ostringstream log;
    ASSERT_TRUE(w.read_molden(path, log, false));
    Int_Params p(w);
    EXPECT_EQ(p.get_nbas(), 4u);
    EXPECT_EQ(p.get_nao(), 8u);
    std::map<int, LibCintBasis> sets = p.get_basis_sets();
    ASSERT_EQ(sets.size(), 1u);
    const LibCintBasis &b = sets[9];
    EXPECT_EQ(b.shelltypes, (ivec{ 0, 1 }));
    EXPECT_EQ(b.shellcount, (ivec{ 4, 4 }));
    ASSERT_EQ(b.exponents.size(), 8u);
    ASSERT_EQ(b.coefficients.size(), 8u);
    EXPECT_NEAR(b.exponents[0], 67.8191594740697, 1e-10);
    EXPECT_NEAR(b.exponents[3], 0.357670001280168, 1e-12);
    EXPECT_NEAR(b.exponents[4], 9.58240403549692, 1e-11);
    EXPECT_NEAR(b.exponents[7], 0.348706630431394, 1e-12);

    // a p shell listed before an s shell is moved behind it
    WFN r = make_wfn({ make_atom("H", 1, 0.0, 0.0, 0.0, { {1, 0.8}, {0, 1.2} }, 1) }, e_origin::molden);
    Int_Params pr(r);
    sets = pr.get_basis_sets();
    EXPECT_EQ(sets[1].shelltypes, (ivec{ 0, 1 }));
    EXPECT_EQ(sets[1].exponents, (vec{ 1.2, 0.8 }));
    EXPECT_EQ(pr.get_nao(), 4u);
}

// an origin without a normalisation rule (xtb) leaves the coefficients untouched; the self overlap is
// then the raw Gaussian one times libcint's s common factor squared
TEST(LibcintDispatchTests, UnknownOriginLeavesCoefficientsUntouched)
{
    WFN w = make_wfn({ make_atom("H", 1, 0.0, 0.0, 0.0, { {0, 1.3} }, 1) }, e_origin::xtb);
    Int_Params p(w);
    std::map<int, LibCintBasis> sets = p.get_basis_sets();
    ASSERT_EQ(sets.size(), 1u);
    EXPECT_EQ(sets[1].coefficients, (vec{ 1.0 }));
    EXPECT_EQ(sets[1].shelltypes, (ivec{ 0 }));
    vec S;
    compute2C<Overlap2C_SPH>(p, S);
    ASSERT_EQ(S.size(), 1u);
    EXPECT_NEAR(S[0], std::pow(constants::PI / 2.6, 1.5) / (4.0 * constants::PI), 1e-12);
}

// print_data writes the three tables to <name>.txt, one atom / shell per line
TEST(LibcintDispatchTests, PrintDataWritesAtmBasEnvSections)
{
    WFN w = make_wfn({ make_atom("H", 1, 0.0, 0.0, 0.0, { {0, 1.2} }), make_atom("He", 2, 0.0, 0.0, 1.5, { {0, 0.7} }) });
    Int_Params p(w);
    const std::filesystem::path stem = std::filesystem::temp_directory_path() / "LibcintDispatch_print_data";
    const std::filesystem::path file = stem.string() + ".txt";
    std::filesystem::remove(file);
    p.print_data(stem.string());
    ASSERT_TRUE(std::filesystem::exists(file));
    const std::string text = read_file(file);
    std::filesystem::remove(file);
    EXPECT_EQ(text.rfind("ATM:\n", 0), 0u);
    EXPECT_NE(text.find("1 20 1 23 0 0 \n2 24 1 27 0 0 \n"), std::string::npos);
    EXPECT_NE(text.find("BAS:\n0 0 1 1 0 28 29 0 \n1 0 1 1 0 30 31 0 \n"), std::string::npos);
    const size_t env_pos = text.find("ENV:\n");
    ASSERT_NE(env_pos, std::string::npos);
    std::istringstream env(text.substr(env_pos + 5));
    vec values;
    double d;
    while (env >> d) values.push_back(d);
    ASSERT_EQ(values.size(), 32u);
    EXPECT_EQ(values[26], 1.5);
    EXPECT_EQ(values[28], 1.2);
}

// run_app with no job writes the banner and the help to the log named by -out and returns 0
TEST(LibcintDispatchIoTests, NoJobWritesHelpToTheRequestedLog)
{
    const std::filesystem::path log = std::filesystem::temp_directory_path() / "LibcintDispatch_nojob.log";
    std::filesystem::remove(log);
    std::vector<std::string> args = { "NoSpherA2", "-out", log.string(), "-no_date" };
    std::vector<char *> argv;
    for (std::string &a : args) argv.push_back(a.data());
    argv.push_back(nullptr);
    const int rc = run_app((int)args.size(), argv.data());
    EXPECT_EQ(rc, 0);
    ASSERT_TRUE(std::filesystem::exists(log));
    const std::string text = read_file(log);
    std::filesystem::remove(log);
    EXPECT_NE(text.find("Did not understand the task to perform!"), std::string::npos);
    EXPECT_NE(text.find("-h, --h, -help, --help"), std::string::npos);
}
