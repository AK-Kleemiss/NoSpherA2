#include "pch.h"

#include "core/nbo.h"

//NBO search on densities whose accepted orbital set is known by hand; corner arithmetic and the
//level-scheduled Gauss-Seidel sweep are equivalences, so results must be exact and thread-independent.

namespace {

    //exit(-1) is status 255 on POSIX and -1 on Windows; accept either, but only a clean exit, not a signal.
    struct ExitedWithErrCheckfCode
    {
        bool operator()(int status) const
        {
            return ::testing::ExitedWithCode(255)(status) || ::testing::ExitedWithCode(-1)(status);
        }
    };

    //Several s NAOs per atom: with one, every block is 1 or 2 wide and any eigensolve looks correct.
    NAOResult h_chain_shells(const int na, const int per_atom)
    {
        NAOResult nao;
        const size_t n = static_cast<size_t>(na) * per_atom;
        nao.C = dMatrix2(n, n);
        for (int a = 0; a < na; a++) {
            for (int s = 0; s < per_atom; s++) {
                NAO o;
                o.atom = a;
                o.l = 0;
                o.m = 0;
                o.shell = s;
                o.n = s + 1;
                o.type = NAOClass::Valence;
                o.occupation = 0.0;
                nao.orbitals.push_back(o);
            }
            NAOAtom at;
            at.index = a;
            at.Z = 1;
            at.label = "H" + std::to_string(a + 1);
            at.Z_eff = 1.0;
            nao.atoms.push_back(at);
        }
        return nao;
    }

    bvec2 chain_bondable(const int na)
    {
        bvec2 b(na, bvec(na, false));
        for (int a = 0; a + 1 < na; a++) {
            b[a][a + 1] = true;
            b[a + 1][a] = true;
        }
        return b;
    }

    //Gamma = occ * sum_k v_k v_k^T over orthonormal vectors.
    dMatrix2 density_of(const std::vector<vec>& v, const double occ)
    {
        const size_t n = v.front().size();
        dMatrix2 g(n, n);
        for (const vec& u : v)
            for (size_t i = 0; i < n; i++)
                for (size_t j = 0; j < n; j++)
                    g(i, j) += occ * u[i] * u[j];
        return g;
    }

    //Bonds 0-1, 1-2 and 3-4, 4-5: within a fragment they share an atom, across fragments nothing,
    //giving two dependency levels of two. Each bond is half on each centre; a fragment's two vectors
    //are orthogonal through the sign on the shared atom's second shell.
    std::vector<vec> two_fragment_bonds()
    {
        const int n = 12;
        std::vector<vec> out;
        for (const int base : { 0, 6 }) {
            vec v1(n, 0.0), v2(n, 0.0);
            v1[base + 0] = 0.5;  v1[base + 1] = 0.5;  v1[base + 2] = 0.5;  v1[base + 3] = 0.5;
            v2[base + 2] = 0.5;  v2[base + 3] = -0.5; v2[base + 4] = 0.5;  v2[base + 5] = 0.5;
            out.push_back(v1);
            out.push_back(v2);
        }
        return out;
    }

    NboLewis search_at(const int threads)
    {
        NboOptions opt;
        opt.threads = threads;
        return nbo_search(h_chain_shells(6, 2), density_of(two_fragment_bonds(), 2.0),
                          chain_bondable(6), 4, 2.0, opt);
    }

}  //namespace

//Checks the corner arithmetic: updates and eigensolve run on the 4 x 4 corner of a 12 x 12 matrix,
//so wrong index bookkeeping shows as a wrong occupancy.
TEST(NboSearchTests, SweepRecoversTheBondsThatAreInTheDensity)
{
    const NboLewis L = search_at(1);
    ASSERT_EQ(L.n_lewis, 4);
    for (int j = 0; j < L.n_lewis; j++) {
        EXPECT_EQ(L.orbitals[j].type, "BD");
        EXPECT_EQ(L.orbitals[j].centers.size(), 2u);
        EXPECT_NEAR(L.orbitals[j].occupancy, 2.0, 1e-9);
        EXPECT_NEAR(L.orbitals[j].center_weight[0], 0.5, 1e-9);
        EXPECT_NEAR(L.orbitals[j].center_weight[1], 0.5, 1e-9);
    }
    EXPECT_NEAR(L.rho_nl, 0.0, 1e-9);
    EXPECT_EQ(L.topo[0][1], 1);
    EXPECT_EQ(L.topo[1][2], 1);
    EXPECT_EQ(L.topo[3][4], 1);
    EXPECT_EQ(L.topo[4][5], 1);
}

//Gauss-Seidel threaded by dependency level is order-preserving, so results must be bitwise equal (EXPECT_EQ on doubles).
TEST(NboSearchTests, ThreadCountCannotMoveASingleBit)
{
    const NboLewis a = search_at(1);
    const NboLewis b = search_at(8);
    ASSERT_EQ(a.orbitals.size(), b.orbitals.size());
    ASSERT_EQ(a.n_lewis, b.n_lewis);
    EXPECT_EQ(a.rho_nl, b.rho_nl);
    EXPECT_EQ(a.threshold, b.threshold);
    for (size_t j = 0; j < a.orbitals.size(); j++) {
        const NboFunction& fa = a.orbitals[j];
        const NboFunction& fb = b.orbitals[j];
        EXPECT_EQ(fa.type, fb.type) << "orbital " << j;
        EXPECT_EQ(fa.centers, fb.centers) << "orbital " << j;
        EXPECT_EQ(fa.occupancy, fb.occupancy) << "orbital " << j;
        ASSERT_EQ(fa.coefficients.size(), fb.coefficients.size());
        for (size_t i = 0; i < fa.coefficients.size(); i++)
            EXPECT_EQ(fa.coefficients[i], fb.coefficients[i]) << "orbital " << j << " coefficient " << i;
        EXPECT_EQ(fa.center_weight, fb.center_weight) << "orbital " << j;
        EXPECT_EQ(fa.center_lchar, fb.center_lchar) << "orbital " << j;
    }
    EXPECT_EQ(a.topo, b.topo);
}

//0.3 electrons on a pair is below the lowest threshold, so no orbital is accepted and the OWSO would face a 0 x 0 eigensolve.
TEST(NboSearchTests, ADensityWithNoLewisStructureSaysSoInsteadOfCrashing)
{
    EXPECT_EXIT(
        {
            std::vector<vec> v(1, vec(4, 0.0));
            v[0][0] = v[0][1] = v[0][2] = v[0][3] = 0.5;
            NboOptions opt;
            nbo_search(h_chain_shells(2, 2), density_of(v, 0.3), chain_bondable(2), 1, 2.0, opt);
        },
        ExitedWithErrCheckfCode(), ".*");
}

//One determinant: every acceptor is filled by its donors alone (s = 0), and the acceptors' q sum to rho(NL).
TEST(NboSearchTests, DensityTransferOfAnIdempotentDensityHasNoRemainder)
{
    const NboLewis L = search_at(1);
    ASSERT_EQ(L.orbitals.size(), 12u);
    std::vector<NboTransfer> t;
    std::vector<NboAcceptorBudget> b;
    std::vector<NboLocalNoon> l;
    density_transfer(L, density_of(two_fragment_bonds(), 2.0), 2.0, 0.0, t, b, l);
    ASSERT_EQ(b.size(), 12u - L.n_lewis);
    double q = 0.0, occ = 0.0;
    for (const NboAcceptorBudget& x : b) {
        EXPECT_NEAR(x.s, 0.0, 1e-10) << "acceptor " << x.index;
        EXPECT_NEAR(x.q_lewis + x.q_nonlewis, x.occupancy, 1e-10) << "acceptor " << x.index;
        q += x.q_lewis + x.q_nonlewis;
        occ += x.occupancy;
    }
    EXPECT_NEAR(occ, L.rho_nl, 1e-10);
    EXPECT_NEAR(q, L.rho_nl, 1e-10);
    for (const NboLocalNoon& x : l) EXPECT_NEAR(x.b_eff, 1.0, 1e-8) << "bond " << x.bond_index;
}

//H2 in two NAOs, natural orbitals rotated by theta off the BD/BD* pair: the 2 x 2 block returns the natural
//occupations, occupancy = q + s, and at n = 1/1 (the dissociated singlet) q = 0, s = 1, b_eff = 0.
TEST(NboSearchTests, DensityTransferOfACorrelatedBondKeepsTheBudget)
{
    const double r = std::sqrt(0.5);
    NboLewis L;
    L.n_lewis = 1;
    NboFunction bd, anti;
    bd.type = "BD", anti.type = "BD*";
    bd.centers = anti.centers = { 0, 1 };
    bd.coefficients = { r, r };
    anti.coefficients = { r, -r };
    L.orbitals = { bd, anti };
    for (const auto& [ng, theta] : std::vector<std::pair<double, double>>{ { 1.9, 0.0 }, { 1.6, 0.3 }, { 1.0, 0.0 } }) {
        const double nu = 2.0 - ng, c = std::cos(theta), s = std::sin(theta);
        //natural orbitals c BD + s BD*, -s BD + c BD* in the NAO basis
        const vec g = { r * (c + s), r * (c - s) }, u = { r * (c - s), -r * (c + s) };
        dMatrix2 gamma(2, 2);
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++) gamma(i, j) = ng * g[i] * g[j] + nu * u[i] * u[j];
        std::vector<NboTransfer> t;
        std::vector<NboAcceptorBudget> b;
        std::vector<NboLocalNoon> l;
        density_transfer(L, gamma, 2.0, 0.0, t, b, l);
        ASSERT_EQ(b.size(), 1u);
        ASSERT_EQ(l.size(), 1u);
        EXPECT_NEAR(b[0].q_lewis + b[0].s, b[0].occupancy, 1e-12);
        EXPECT_NEAR(l[0].n_plus, ng, 1e-12);
        EXPECT_NEAR(l[0].n_minus, nu, 1e-12);
        EXPECT_NEAR(l[0].b_eff, 0.5 * (ng - nu), 1e-12);
        if (theta == 0.0) EXPECT_NEAR(b[0].s, nu, 1e-12) << "n = " << ng;
    }
}

//One determinant in two NAOs, a donor mixed into an acceptor by F_kj: E2(D) = q dE and E2 = q_k F^2/dE agree to
//second order in F/dE (here 2 %), and E2(D) is 0 without a Fock.
TEST(NboSearchTests, DensityTransferPricesOneDeterminantAtE2)
{
    NboLewis L;
    L.n_lewis = 1;
    NboFunction donor, acceptor;
    donor.type = "LP", acceptor.type = "LV";
    donor.centers = acceptor.centers = { 0 };
    donor.coefficients = { 1.0, 0.0 };
    acceptor.coefficients = { 0.0, 1.0 };
    L.orbitals = { donor, acceptor };
    const double ek = -0.5, ej = 0.3, f = 0.04;
    dMatrix2 fock(2, 2);
    fock(0, 0) = ek, fock(1, 1) = ej, fock(0, 1) = fock(1, 0) = f;
    //the occupied eigenvector of the 2 x 2 Fock, (cos, -sin) with tan 2 theta = 2 F/dE, doubly occupied
    const double th = 0.5 * std::atan(2.0 * f / (ej - ek)), c = std::cos(th), s = std::sin(th);
    dMatrix2 gamma(2, 2);
    gamma(0, 0) = 2.0 * c * c, gamma(1, 1) = 2.0 * s * s, gamma(0, 1) = gamma(1, 0) = -2.0 * c * s;
    L.orbitals[0].occupancy = gamma(0, 0);
    L.orbitals[1].occupancy = gamma(1, 1);
    const std::vector<NboE2Entry> e2 = nbo_e2(L, fock, 0.0);
    ASSERT_EQ(e2.size(), 1u);
    std::vector<NboTransfer> t, bare;
    std::vector<NboAcceptorBudget> b;
    std::vector<NboLocalNoon> l;
    density_transfer(L, gamma, 2.0, 0.0, t, b, l, fock);
    ASSERT_EQ(t.size(), 1u);
    EXPECT_NEAR(t[0].e2d_kcal, e2[0].energy_kcal, 0.02 * e2[0].energy_kcal);
    density_transfer(L, gamma, 2.0, 0.0, bare, b, l);
    ASSERT_EQ(bare.size(), 1u);
    EXPECT_EQ(bare[0].e2d_kcal, 0.0);
}

//|F| at or above the gap: second order exceeds the exact two-level gain q|F|, flagged; at ten times the gap left out.
TEST(NboSearchTests, E2FlagsAndDropsPairsBeyondSecondOrder)
{
    NboLewis L;
    L.n_lewis = 1;
    NboFunction donor, acceptor;
    donor.type = "LP", acceptor.type = "LV";
    donor.centers = acceptor.centers = { 0 };
    donor.coefficients = { 1.0, 0.0 };
    acceptor.coefficients = { 0.0, 1.0 };
    donor.occupancy = 2.0, acceptor.occupancy = 0.0;
    L.orbitals = { donor, acceptor };
    dMatrix2 fock(2, 2);
    fock(0, 0) = -0.5, fock(1, 1) = -0.495, fock(0, 1) = fock(1, 0) = 0.1;
    int beyond = 0;
    EXPECT_TRUE(nbo_e2(L, fock, 0.0, &beyond).empty());
    EXPECT_EQ(beyond, 1);
    fock(1, 1) = -0.45;
    beyond = 0;
    std::vector<NboE2Entry> e2 = nbo_e2(L, fock, 0.0, &beyond);
    ASSERT_EQ(e2.size(), 1u);
    EXPECT_TRUE(e2[0].beyond_second_order);
    EXPECT_EQ(beyond, 0);
    fock(1, 1) = 0.3;
    e2 = nbo_e2(L, fock, 0.0, &beyond);
    ASSERT_EQ(e2.size(), 1u);
    EXPECT_FALSE(e2[0].beyond_second_order);
    EXPECT_EQ(beyond, 0);
}
