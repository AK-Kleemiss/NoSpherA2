#include "pch.h"

#include "core/nbo.h"

//The NBO search, on densities small enough that the accepted orbital set is known by hand.  The
//search is where the run spends its time, so it is also where it is optimised - corner arithmetic
//in place of n x n bookkeeping, and the Gauss-Seidel sweep run one dependency level at a time.  Both
//of those are equivalences, not approximations, and an equivalence needs a check that fails when it
//stops holding: the sweep has to find the orbitals that are actually in the density, and the number
//of threads must not move a single bit of what it reports.
//
//The corpus comparison against the 22 stored NBO 7 references lives in tests/nbo_reference and is
//run by compare_nbo.py; this file is the part that fails when the algebra breaks.

namespace {

    //err_checkf ends the process with exit(-1).  POSIX reports that as wait status 255, Windows
    //hands gtest the raw -1, so pinning one number passes on one platform and fails on the other -
    //which is exactly what CI's Windows Release and Windows GPU Release jobs were failing on while
    //Linux and macOS were green.  Accept either code, but still only a CLEAN exit: this test exists
    //to show the search refuses instead of crashing, so an access violation or a signal must not
    //satisfy it.  ExitedWithCode does that platform check for us, so delegate to it twice rather
    //than reimplementing WIFEXITED here.
    struct ExitedWithErrCheckfCode
    {
        bool operator()(int status) const
        {
            return ::testing::ExitedWithCode(255)(status) || ::testing::ExitedWithCode(-1)(status);
        }
    };

    //Several s NAOs per atom, so an atom block is wider than the orbital taken out of it and a pair
    //block is wider still - with one NAO per atom every block is 1 or 2 wide and any way of
    //computing its leading eigenpair looks correct.
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

    //Two three-atom fragments, two s NAOs per atom, and four bonds: 0-1 and 1-2 on the first
    //fragment, 3-4 and 4-5 on the second.  The two bonds of a fragment share their middle atom, so
    //the sweep has to keep their order; bonds on different fragments share nothing, so it may run
    //them together.  That is exactly the dependency structure the level scheduling exploits, and
    //here it is two levels of two.  Every vector is half on each of its centres, which keeps it a
    //bond rather than a polarised near-lone-pair, and the two vectors of a fragment are orthogonal
    //through the sign on the shared atom's second shell.
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

//The density is four doubly occupied two-centre orbitals, so the search has to come back with those
//four and nothing outside them.  This is the check on the corner arithmetic: the rank-one updates,
//the difference the eigensolve sees and the eigensolve itself all happen on the 4 x 4 corner of a
//12 x 12 matrix, and getting the index bookkeeping of that corner wrong shows up here as a wrong
//occupancy rather than as a crash.
TEST(NboSearchTests, SweepRecoversTheBondsThatAreInTheDensity)
{
    const NboLewis L = search_at(1);
    ASSERT_EQ(L.n_lewis, 4);
    for (int j = 0; j < L.n_lewis; j++) {
        EXPECT_EQ(L.orbitals[j].type, "BD");
        EXPECT_EQ(L.orbitals[j].centers.size(), 2u);
        EXPECT_NEAR(L.orbitals[j].occupancy, 2.0, 1e-9);
        //half on each centre, which is what makes it a bond
        EXPECT_NEAR(L.orbitals[j].center_weight[0], 0.5, 1e-9);
        EXPECT_NEAR(L.orbitals[j].center_weight[1], 0.5, 1e-9);
    }
    EXPECT_NEAR(L.rho_nl, 0.0, 1e-9);
    EXPECT_EQ(L.topo[0][1], 1);
    EXPECT_EQ(L.topo[1][2], 1);
    EXPECT_EQ(L.topo[3][4], 1);
    EXPECT_EQ(L.topo[4][5], 1);
}

//The sweep is Gauss-Seidel: orbital j reads the density with every other accepted orbital already
//removed, including the ones updated earlier in the same sweep.  It is threaded by level, so
//orbitals that share a centre still run in order and only independent ones run together.  If that
//is right the thread count cannot change anything, down to the last bit - not "to 1e-12", bitwise,
//which is why this compares with EXPECT_EQ on doubles.  A number that moves here is a number the
//published tables lose.
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

//A density with nothing in it: 0.3 electrons spread over a pair, below the bottom rung of the
//threshold ladder, so not one orbital is accepted.  tests/molden_file/F2.molden and epoxide.molden
//arrive at the search exactly like this - their FILE47 comes out with a fraction of the electrons
//the molecule has - and they used to segfault inside a 0 x 0 eigensolve in the OWSO.  There is no
//Lewis structure in such an input, so the search has to say which input it was rather than crash.
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
