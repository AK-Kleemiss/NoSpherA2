#include "pch.h"

#include "core/nbo.h"

//Natural resonance theory, checked on inputs small enough that the answer is known by hand rather
//than by a golden file.  The internals of nrt.cpp are in an anonymous namespace on purpose, so the
//checks drive the public entry point: a synthetic NAO basis, a synthetic density and a parent Lewis
//structure go in, and what comes out has to satisfy the identities the printed NBO tables satisfy.
//
//The corpus comparison against the 22 stored NBO 7 references lives in tests/nbo_reference and is
//run by compare_nbo.py; this file is the part that fails when the algebra breaks, not the chemistry.

namespace {

    //One s-type valence NAO per atom, which is all these checks need: the topology bookkeeping, the
    //OWSO orbital construction and the QP never look at l.
    NAOResult h_chain(const int na)
    {
        NAOResult nao;
        nao.C = dMatrix2(na, na);
        for (int a = 0; a < na; a++) {
            NAO o;
            o.atom = a;
            o.l = 0;
            o.m = 0;
            o.shell = 0;
            o.n = 1;
            o.type = NAOClass::Valence;
            o.occupation = 0.0;
            nao.orbitals.push_back(o);
            NAOAtom at;
            at.index = a;
            at.Z = 1;
            at.label = "H" + std::to_string(a + 1);
            at.Z_eff = 1.0;
            nao.atoms.push_back(at);
        }
        return nao;
    }

    //Gamma = 2 v v^T: one doubly occupied orbital, so the parent structure that puts exactly this
    //orbital where v lives reproduces the density and leaves nothing outside it.
    dMatrix2 rank_one(const vec& v)
    {
        const size_t n = v.size();
        dMatrix2 g(n, n);
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                g(i, j) = 2.0 * v[i] * v[j];
        return g;
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

    //A single two-centre bond over atoms 0 and 1, no cores and no lone pairs.
    NboLewis one_bond(const dMatrix2& gamma, const int na)
    {
        NboLewis L;
        L.gamma = gamma;
        NboFunction f;
        f.centers = { 0, 1 };
        f.type = "BD";
        f.multiplicity = 1;
        f.occupancy = 2.0;
        L.orbitals.push_back(f);
        L.n_lewis = 1;
        L.topo.assign(na, ivec(na, 0));
        L.topo[0][1] = 1;
        L.topo[1][0] = 1;
        return L;
    }

    const NboValency& valency_of(const NboNrt& nrt, const int atom)
    {
        for (const NboValency& v : nrt.valencies)
            if (v.atom == atom) return v;
        throw std::runtime_error("no valency for atom " + std::to_string(atom));
    }

    double bond_total(const NboNrt& nrt, const int a, const int b)
    {
        for (const NboBondOrder& o : nrt.bond_orders)
            if (o.atom1 == std::min(a, b) && o.atom2 == std::max(a, b) && (o.diagonal == (a == b)))
                return o.total;
        return 0.0;
    }

}  //namespace

//The density is exactly the parent's own orbital set, so the residual has to vanish: one structure
//at 100 %, D(w) = 0, no electrons outside the Lewis set, bond order 1 and two electrons on each H.
TEST(NrtTests, ParentThatSpansTheDensityTakesAllTheWeight)
{
    const NAOResult nao = h_chain(2);
    const double r = 1.0 / std::sqrt(2.0);
    const dMatrix2 gamma = rank_one({ r, r });
    const NboLewis lewis = one_bond(gamma, 2);
    NboOptions opt;
    opt.nrt = true;
    NboNrt nrt;
    std::ostringstream log;
    native_nrt(nrt, nao, lewis, {}, chain_bondable(2), opt, "", 2.0, log);

    ASSERT_TRUE(nrt.present);
    EXPECT_EQ(nrt.structures_found, 1);
    EXPECT_EQ(nrt.structures_used, 1);
    ASSERT_EQ(nrt.weights.size(), 1u);
    EXPECT_NEAR(nrt.weights[0].weight_percent, 100.0, 1e-6);
    EXPECT_NEAR(nrt.d_w, 0.0, 1e-8);
    EXPECT_NEAR(nrt.candidates[0].rho_nl, 0.0, 1e-8);

    EXPECT_NEAR(bond_total(nrt, 1, 2), 1.0, 1e-8);
    EXPECT_NEAR(bond_total(nrt, 1, 1), 0.0, 1e-8);
    const NboValency& v = valency_of(nrt, 1);
    EXPECT_NEAR(v.valency, 1.0, 1e-8);
    EXPECT_NEAR(v.electron_count, 2.0, 1e-8);
    ASSERT_EQ(nrt.leading_topo.size(), 1u);
    EXPECT_EQ(nrt.leading_topo[0].matrix[0][1], 1);
}

//The ionic share of a bond order is |i| b with i = c_A^2 - c_B^2, linear in |i|, which is what
//acetylene's printed table says (i = 0.2334, ionic/total = 0.2262/0.9693 = 0.23336) and what a
//plausible-looking i^2 would get wrong.  Here the bond orbital is (sqrt(0.8), sqrt(0.2)) by
//construction, so i = 0.6 exactly: ionic 0.6, covalent 0.4, and i^2 would give 0.36.
TEST(NrtTests, IonicBondOrderIsLinearInThePolarity)
{
    const NAOResult nao = h_chain(2);
    const dMatrix2 gamma = rank_one({ std::sqrt(0.8), std::sqrt(0.2) });
    const NboLewis lewis = one_bond(gamma, 2);
    NboOptions opt;
    opt.nrt = true;
    NboNrt nrt;
    std::ostringstream log;
    native_nrt(nrt, nao, lewis, {}, chain_bondable(2), opt, "", 2.0, log);

    ASSERT_EQ(nrt.bond_orders.size(), 3u);  //one pair plus the two diagonals
    const NboBondOrder& o = nrt.bond_orders[0];
    EXPECT_FALSE(o.diagonal);
    EXPECT_NEAR(o.total, 1.0, 1e-8);
    EXPECT_NEAR(o.ionic, 0.6, 1e-6);
    EXPECT_NEAR(o.covalent, 0.4, 1e-6);
    EXPECT_NEAR(valency_of(nrt, 1).electrovalency, 0.6, 1e-6);
    EXPECT_NEAR(valency_of(nrt, 1).covalency, 0.4, 1e-6);
}

//With more than one candidate the answer is no longer known by hand, but the identities are: the
//weights are a probability vector, the minimiser cannot do worse than the parent alone, every
//atom's valency is its bond-order row sum and its electron count is 2 (lone pairs + valency).
//The exhaustive mode is used so the check does not depend on the arrow search or on an E2 table.
TEST(NrtTests, DerivedQuantitiesAreConsistentOverAMultiStructureFit)
{
    const NAOResult nao = h_chain(3);
    const double r = 1.0 / std::sqrt(3.0);
    const dMatrix2 gamma = rank_one({ r, r, r });
    NboLewis lewis = one_bond(gamma, 3);
    NboOptions opt;
    opt.nrt = true;
    opt.nrt_exhaustive = true;
    NboNrt nrt;
    std::ostringstream log;
    native_nrt(nrt, nao, lewis, {}, chain_bondable(3), opt, "", 2.0, log);

    ASSERT_TRUE(nrt.present);
    EXPECT_GT(nrt.structures_found, 1);
    double sum = 0.0;
    for (const NboResonanceWeight& w : nrt.weights) {
        EXPECT_GE(w.weight_percent, 0.0);
        sum += w.weight_percent;
    }
    EXPECT_NEAR(sum, 100.0, 1e-3);          //the dropped sub-floor weights are below 5e-3 %
    EXPECT_LE(nrt.d_w, nrt.d_0 + 1e-9);

    for (int a = 1; a <= 3; a++) {
        double row = 0.0;
        for (int b = 1; b <= 3; b++)
            if (b != a) row += bond_total(nrt, a, b);
        const NboValency& v = valency_of(nrt, a);
        EXPECT_NEAR(v.valency, row, 1e-8);
        EXPECT_NEAR(v.valency, v.covalency + v.electrovalency, 1e-8);
        EXPECT_NEAR(v.electron_count, 2.0 * (bond_total(nrt, a, a) + v.valency), 1e-8);
    }
    //one electron pair, distributed over the chain: the total bond order plus the lone pairs is 1
    double pairs = 0.0;
    for (const NboBondOrder& o : nrt.bond_orders) pairs += o.total;
    EXPECT_NEAR(pairs, 1.0, 1e-3);
}
