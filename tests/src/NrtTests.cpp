#include "pch.h"

#include "core/constants.h"
#include "core/nbo.h"

//Natural resonance theory on inputs small enough that the answer is known by hand. nrt.cpp's internals sit
//in an anonymous namespace, so the checks drive the public entry point with a synthetic NAO basis, density
//and parent Lewis structure. The NBO 7 corpus comparison is tests/nbo_reference (compare_nbo.py).

namespace {

    //One s-type valence NAO per atom: the topology bookkeeping, the OWSO construction and the QP never look at l.
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

    //Several s NAOs per atom, so a block is wider than the orbital taken out of it; with one per atom every
    //block is 1 or 2 wide and any way of computing the leading eigenpair looks correct.
    NAOResult h_chain_shells(const int na, const int per_atom)
    {
        NAOResult nao;
        nao.C = dMatrix2(static_cast<size_t>(na) * per_atom, static_cast<size_t>(na) * per_atom);
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

    //Gamma = 2 (v1 v1^T + v2 v2^T) for orthonormal v1, v2: two doubly occupied orbitals.
    dMatrix2 rank_two(const vec& v1, const vec& v2)
    {
        const size_t n = v1.size();
        dMatrix2 g(n, n);
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                g(i, j) = 2.0 * (v1[i] * v1[j] + v2[i] * v2[j]);
        return g;
    }

    vec normalised(vec v)
    {
        double n = 0.0;
        for (const double x : v) n += x * x;
        n = std::sqrt(n);
        for (double& x : v) x /= n;
        return v;
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

//The density is exactly the parent's orbital set, so the residual vanishes and the parent takes all the weight.
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

//The ionic share of a bond order is |i| b with i = c_A^2 - c_B^2, linear in |i| as NBO's acetylene table
//shows. Here i = 0.6 exactly: ionic 0.6, covalent 0.4, where i^2 would give 0.36.
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

//An open-shell NRT runs once per spin and its integer unit is a single electron (a pair closed shell).
//Bond orders count pairs either way, so one alpha electron shared by two H is half a bond and the channel
//holds one electron. ParentThatSpansTheDensityTakesAllTheWeight is the closed-shell twin, with exactly
//twice these numbers and `scale` the only difference.
TEST(NrtTests, AnOpenShellSpinChannelCountsItsBondOrdersInPairsNotElectrons)
{
    const NAOResult nao = h_chain(2);
    const double r = 1.0 / std::sqrt(2.0);
    //Gamma = 1 v v^T, not rank_one()'s 2 v v^T: one electron in this spin channel, as the Lewis bond's occupancy says.
    dMatrix2 gamma(2, 2);
    for (size_t i = 0; i < 2; i++)
        for (size_t j = 0; j < 2; j++)
            gamma(i, j) = r * r;
    NboLewis lewis = one_bond(gamma, 2);
    lewis.orbitals[0].occupancy = 1.0;
    NboOptions opt;
    opt.nrt = true;
    NboNrt nrt;
    std::ostringstream log;
    native_nrt(nrt, nao, lewis, {}, chain_bondable(2), opt, "alpha", 1.0, log);

    ASSERT_TRUE(nrt.present);
    EXPECT_NEAR(bond_total(nrt, 1, 2), 0.5, 1e-8)
        << "one electron shared between two centres is half a bond; 1.0 would be counting the spin's "
           "single electron as a pair";
    const NboValency& v = valency_of(nrt, 1);
    EXPECT_NEAR(v.valency, 0.5, 1e-8);
    EXPECT_NEAR(v.electron_count, 1.0, 1e-8)
        << "the alpha channel of this system holds exactly one electron, so no bookkeeping derived "
           "from it may report two";
    //every bond order and lone pair of a spin channel sums to half its electron count: the weights are a
    //probability vector over topologies that each place the same number of units
    double pairs = 0.0;
    for (const NboBondOrder& o : nrt.bond_orders) pairs += o.total;
    EXPECT_NEAR(pairs, 0.5, 1e-8)
        << "the bond orders and lone pairs of a one-electron spin channel must sum to 0.5 pairs";
}

//A closed-shell density split into two identical spin channels must reproduce the closed-shell answer: the
//same density runs as gamma (occupancy 2, scale 2) and twice as gamma/2 (occupancy 1, scale 1). The totals
//must add up (the unit conversion, scale/2 in native_nrt), and the ionic fraction of each channel must equal
//the closed-shell one, as an orbital's polarity does not depend on its occupation. Two bonds on purpose: for
//one orbital the OWSO step is the identity (M w (w S w)^-1/2 = M for k = 1), so the fraction checks could not
//fail; two bonds sharing atom 2 overlap, and the weighting decides how. The invariant holds to roundoff,
//hence the 1e-12 fraction tolerances.
TEST(NrtTests, TwoIdenticalSpinChannelsAddUpToTheClosedShellAnswer)
{
    const NAOResult nao = h_chain(3);
    //Two polar bond orbitals, 1-2 and 2-3, not orthogonal to each other: they share atom 2, which is
    //what gives the OWSO step something to do.
    const vec u1 = normalised({ std::sqrt(0.8), std::sqrt(0.2), 0.0 });
    const vec u2 = normalised({ 0.0, std::sqrt(0.3), std::sqrt(0.7) });
    const dMatrix2 gamma = rank_two(u1, u2);   //2 (u1 u1^T + u2 u2^T)

    //The parent: one bond on each bondable pair of the chain.
    auto two_bond_lewis = [](const dMatrix2& g, const double occ) {
        NboLewis L;
        L.gamma = g;
        for (const std::vector<int>& c : { std::vector<int>{ 0, 1 }, std::vector<int>{ 1, 2 } }) {
            NboFunction f;
            f.centers = c;
            f.type = "BD";
            f.multiplicity = 1;
            f.occupancy = occ;
            L.orbitals.push_back(f);
        }
        L.n_lewis = 2;
        L.topo.assign(3, ivec(3, 0));
        L.topo[0][1] = L.topo[1][0] = 1;
        L.topo[1][2] = L.topo[2][1] = 1;
        return L;
    };

    auto bond_of = [](const NboNrt& n, const int a, const int b) {
        for (const NboBondOrder& o : n.bond_orders)
            if (!o.diagonal && o.atom1 == a && o.atom2 == b) return o;
        throw std::runtime_error("no bond order for the pair");
    };

    NboOptions opt;
    opt.nrt = true;
    opt.nrt_exhaustive = true;   //so both routes see the same candidate set by construction
    std::ostringstream log;

    NboNrt closed;
    NboLewis lc = two_bond_lewis(gamma, 2.0);
    native_nrt(closed, nao, lc, {}, chain_bondable(3), opt, "", 2.0, log);
    ASSERT_TRUE(closed.present);

    dMatrix2 half(3, 3);
    for (size_t i = 0; i < 3; i++)
        for (size_t j = 0; j < 3; j++)
            half(i, j) = 0.5 * gamma(i, j);

    double spin_total[2] = { 0.0, 0.0 }, spin_ionic[2] = { 0.0, 0.0 }, spin_electrons = 0.0;
    const int pair_a[2] = { 1, 2 }, pair_b[2] = { 2, 3 };
    for (const char* spin : { "alpha", "beta" }) {
        NboNrt s;
        NboLewis ls = two_bond_lewis(half, 1.0);
        native_nrt(s, nao, ls, {}, chain_bondable(3), opt, spin, 1.0, log);
        ASSERT_TRUE(s.present) << spin;
        ASSERT_EQ(s.spin_d.size(), 1u);
        EXPECT_NEAR(s.spin_d[0].d_0, closed.spin_d[0].d_0, 1e-10) << spin << ": NBO's D is per spin, Gamma / 2 is each spin of the closed shell";
        EXPECT_NEAR(s.spin_d[0].d_w, closed.spin_d[0].d_w, 1e-10) << spin;
        for (int k = 0; k < 2; k++) {
            const NboBondOrder o = bond_of(s, pair_a[k], pair_b[k]);
            const NboBondOrder c = bond_of(closed, pair_a[k], pair_b[k]);
            ASSERT_GT(o.total, 0.0) << spin;
            EXPECT_NEAR(o.total, 0.5 * c.total, 1e-8)
                << spin << " bond " << pair_a[k] << "-" << pair_b[k]
                << ": one electron of this spin where the closed-shell run has a pair is half the "
                   "bond order, not the same bond order";
            EXPECT_NEAR(o.ionic / o.total, c.ionic / c.total, 1e-12)
                << spin << " bond " << pair_a[k] << "-" << pair_b[k]
                << ": the ionic FRACTION is a property of the orbital's polarity and cannot depend on "
                   "whether one electron or two occupy it";
            spin_total[k] += o.total;
            spin_ionic[k] += o.ionic;
        }
        spin_electrons += valency_of(s, 2).electron_count;
    }

    for (int k = 0; k < 2; k++) {
        const NboBondOrder c = bond_of(closed, pair_a[k], pair_b[k]);
        EXPECT_NEAR(spin_total[k], c.total, 1e-8)
            << "the two spin channels of an unpolarised density must add up to the closed-shell bond "
               "order on " << pair_a[k] << "-" << pair_b[k] << ": " << spin_total[k] << " against "
            << c.total;
        EXPECT_NEAR(spin_ionic[k], c.ionic, 1e-12)
            << "and to its ionic part on " << pair_a[k] << "-" << pair_b[k] << ": " << spin_ionic[k]
            << " against " << c.ionic;
    }
    EXPECT_NEAR(spin_electrons, valency_of(closed, 2).electron_count, 1e-8)
        << "each channel's electron count is the closed-shell one halved, so the two add back to it "
           "and not to twice it";
}

//With several candidates the identities still hold: the weights are a probability vector, the fit is no
//worse than the parent alone, each atom's valency is its bond-order row sum and its electron count is 2.
//Exhaustive mode, so the check needs no arrow search or E2 table.
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

//The self-consistency sweep must find the exact orbitals, which greedy filling cannot: gamma =
//2(v1 v1^T + v2 v2^T) with v1 on atoms 1-2 and v2 on 2-3, so the first 4-wide block also holds part of v2
//and its leading eigenvector mixes the two. Only removing the other orbital and re-solving recovers v1 and
//v2, and then D(w) = 0. Two NAOs per atom make it a real eigenproblem for leading_block's power iteration.
TEST(NrtTests, TheSweepRecoversTheExactOrbitalsOutOfOverlappingBlocks)
{
    const NAOResult nao = h_chain_shells(3, 2);
    //orthogonal by construction: they meet only on atom 2, where v2 is antisymmetric and v1 is not
    const vec v1 = normalised({ 1.0, 1.0, 0.5, 0.5, 0.0, 0.0 });
    const vec v2 = normalised({ 0.0, 0.0, 0.5, -0.5, 1.0, 1.0 });
    const dMatrix2 gamma = rank_two(v1, v2);

    NboLewis lewis;
    lewis.gamma = gamma;
    for (const std::pair<int, int> b : { std::make_pair(0, 1), std::make_pair(1, 2) }) {
        NboFunction f;
        f.centers = { b.first, b.second };
        f.type = "BD";
        f.multiplicity = 1;
        f.occupancy = 2.0;
        lewis.orbitals.push_back(f);
    }
    lewis.n_lewis = 2;
    lewis.topo.assign(3, ivec(3, 0));
    lewis.topo[0][1] = lewis.topo[1][0] = 1;
    lewis.topo[1][2] = lewis.topo[2][1] = 1;

    NboOptions opt;
    opt.nrt = true;
    NboNrt nrt;
    std::ostringstream log;
    native_nrt(nrt, nao, lewis, {}, chain_bondable(3), opt, "", 2.0, log);

    ASSERT_TRUE(nrt.present);
    ASSERT_FALSE(nrt.candidates.empty());
    //The floor is the sweep's own fixed point, not the eigenpair: 50 sweeps of a linearly convergent iteration
    //leave D(w) near 5e-8.
    EXPECT_NEAR(nrt.d_w, 0.0, 1e-7);                      //the parent reproduces gamma
    EXPECT_NEAR(nrt.candidates[0].rho_nl, 0.0, 1e-7);     //nothing left outside the Lewis set
    EXPECT_NEAR(bond_total(nrt, 1, 2), 1.0, 1e-6);
    EXPECT_NEAR(bond_total(nrt, 2, 3), 1.0, 1e-6);
    //four electrons, two pairs: bond orders plus lone pairs
    double pairs = 0.0;
    for (const NboBondOrder& o : nrt.bond_orders) pairs += o.total;
    EXPECT_NEAR(pairs, 2.0, 1e-3);
}

//A candidate limit returns a subset of the answer: the 2-3 bond forms under a 25 kcal/mol interaction and
//the 3-4 bond under a 2 kcal/mol one, so a budget of four goes to the first and the cheap structure is lost.
TEST(NrtTests, ASmallCandidateBudgetKeepsTheExpensiveArrowsAndMoreThanTheParent)
{
    const NAOResult nao = h_chain(4);
    //Gamma must contain the delocalised structure, as nrt.candidates holds only structures above the weight
    //floor: 85 % BD(1,2) + LP(3) + LP(4) and 15 % LP(1) + BD(2,3) + LP(4), the 25 kcal/mol arrow's structure.
    //Both are three pairs over four NAOs, so the trace stays at six electrons.
    const double r = 1.0 / std::sqrt(2.0);
    dMatrix2 gamma(4, 4);
    const std::vector<std::pair<double, std::vector<vec>>> mix = {
        { 0.85, { vec{ r, r, 0.0, 0.0 }, vec{ 0.0, 0.0, 1.0, 0.0 }, vec{ 0.0, 0.0, 0.0, 1.0 } } },
        { 0.15, { vec{ 1.0, 0.0, 0.0, 0.0 }, vec{ 0.0, r, r, 0.0 }, vec{ 0.0, 0.0, 0.0, 1.0 } } },
    };
    for (const auto& part : mix)
        for (const vec& v : part.second)
            for (size_t i = 0; i < 4; i++)
                for (size_t j = 0; j < 4; j++)
                    gamma(i, j) += part.first * 2.0 * v[i] * v[j];

    NboLewis lewis;
    lewis.gamma = gamma;
    NboFunction bd;
    bd.centers = { 0, 1 };
    bd.type = "BD";
    bd.multiplicity = 1;
    bd.occupancy = 2.0;
    lewis.orbitals.push_back(bd);
    for (const int a : { 2, 3 }) {
        NboFunction lp;
        lp.centers = { a };
        lp.type = "LP";
        lp.multiplicity = 1;
        lp.occupancy = 2.0;
        lewis.orbitals.push_back(lp);
    }
    lewis.n_lewis = 3;
    lewis.topo.assign(4, ivec(4, 0));
    lewis.topo[0][1] = lewis.topo[1][0] = 1;
    lewis.topo[2][2] = lewis.topo[3][3] = 1;

    //LP(3)->BD(1,2) covers atoms 1,2,3 and marks the 1-2 and 2-3 bonds at 25, LP(4)->LP(3) marks 3-4 at 2;
    //both are above the 1 kcal/mol resonance threshold.
    std::vector<NboE2Entry> e2(2);
    e2[0].donor_index = 2;      //LP on atom 3
    e2[0].acceptor_index = 1;   //BD 1-2
    e2[0].energy_kcal = 25.0;
    e2[1].donor_index = 3;      //LP on atom 4
    e2[1].acceptor_index = 2;   //LP on atom 3
    e2[1].energy_kcal = 2.0;

    const auto pair_seen = [](const NboNrt& n, const int a, const int b) {
        for (const NboNrtCandidate& c : n.candidates)
            if (c.topo[a][b] > 0) return true;
        return false;
    };

    NboOptions opt;
    opt.nrt = true;
    opt.nrt_e2_kcal = 1.0;
    opt.nrt_max_set = true;     //obey the number, do not let the size guard pick one
    std::ostringstream log;

    opt.nrt_max_candidates = 10000;
    NboNrt full;
    native_nrt(full, nao, lewis, e2, chain_bondable(4), opt, "", 2.0, log);
    ASSERT_TRUE(full.present);
    ASSERT_GT(full.structures_found, 2);
    ASSERT_TRUE(pair_seen(full, 1, 2));   //the 25 kcal/mol structure carries weight, so it is visible

    opt.nrt_max_candidates = 4;
    NboNrt tight;
    native_nrt(tight, nao, lewis, e2, chain_bondable(4), opt, "", 2.0, log);
    ASSERT_TRUE(tight.present);
    EXPECT_GT(tight.structures_found, 1);                      //not the parent alone
    EXPECT_LE(tight.structures_found, full.structures_found);
    EXPECT_TRUE(pair_seen(tight, 1, 2));                       //the 25 kcal/mol arrows survive
    EXPECT_LE(tight.d_w, tight.d_0 + 1e-9);

    //`arrows` in the nbo_json holds arrows only; the budget and the dropped half-arrow intermediates are notes.
    for (const std::string& a : tight.arrows)
        EXPECT_EQ(a.rfind("ARROWS", 0), 0u) << a;
    for (const std::string& a : full.arrows)
        EXPECT_EQ(a.rfind("ARROWS", 0), 0u) << a;
    EXPECT_FALSE(tight.notes.empty());

    //-nrt_target 0: walk again from every weighted structure until none is left; never a worse fit, and the
    //notes say why it stopped.
    opt.nrt_max_candidates = 10000;
    opt.nrt_target = 0.0;
    NboNrt grown;
    native_nrt(grown, nao, lewis, e2, chain_bondable(4), opt, "", 2.0, log);
    ASSERT_TRUE(grown.present);
    EXPECT_GE(grown.structures_found, full.structures_found);
    EXPECT_LE(grown.d_w, full.d_w + 1e-6);
    bool stop_note = false;
    for (const std::string& n : grown.notes)
        stop_note = stop_note || n.find("D(w) target") != std::string::npos;
    EXPECT_TRUE(stop_note);
}

TEST(NrtTests, EveryResonanceRowTheResultCarriesIsAlsoPrinted)
{
    const NAOResult nao = h_chain(3);
    const double r = 1.0 / std::sqrt(3.0);
    const dMatrix2 gamma = rank_one({ r, r, r });
    NboLewis lewis = one_bond(gamma, 3);
    NboOptions opt;
    opt.nrt = true;
    opt.nrt_exhaustive = true;
    NboResults res;
    std::ostringstream log;
    native_nrt(res.nrt, nao, lewis, {}, chain_bondable(3), opt, "", 2.0, log);
    ASSERT_TRUE(res.nrt.present);
    ASSERT_FALSE(res.nrt.weights.empty());
    ASSERT_FALSE(res.nrt.bond_orders.empty());
    ASSERT_FALSE(res.nrt.valencies.empty());

    std::ostringstream out;
    print_nrt(res, out);
    const std::string text = out.str();
    const auto fixed_str = [](const double v, const int prec) {
        std::ostringstream s;
        s << std::fixed << std::setprecision(prec) << v;
        return s.str();
    };
    const auto printed = [&text](const std::string& needle) {
        return text.find(needle) != std::string::npos;
    };

    EXPECT_TRUE(printed("Natural resonance theory")) << text;
    EXPECT_TRUE(printed(fixed_str(res.nrt.d_w, 5))) << text;
    for (const NboResonanceWeight& w : res.nrt.weights)
        EXPECT_TRUE(printed(fixed_str(w.weight_percent, 2)))
            << "weight of structure " << w.structure << " missing from\n" << text;
    for (const NboBondOrder& b : res.nrt.bond_orders)
        EXPECT_TRUE(printed(fixed_str(b.total, 4)))
            << "bond order " << b.atom1 << "-" << b.atom2 << " missing from\n" << text;
    for (const NboValency& v : res.nrt.valencies)
        EXPECT_TRUE(printed(fixed_str(v.electron_count, 4)))
            << "electron count of atom " << v.atom << " missing from\n" << text;
    for (const std::string& s : res.nrt.notes)
        EXPECT_TRUE(printed(s)) << "note missing from\n" << text;
    //the element labels come from the populations table, which only a full native run fills; the
    //valencies carry them too, so the printer must not fall back to "?" here
    EXPECT_TRUE(printed("H1")) << text;
    EXPECT_FALSE(printed("?")) << text;

    //and nothing when the analysis did not run
    NboResults empty;
    std::ostringstream none;
    print_nrt(empty, none);
    EXPECT_TRUE(none.str().empty()) << none.str();
}

//The parent goes out ahead of the NRT section as a matrix and as the JSON an NRT input reads back,
//and into the .native.nbo.json under "parent_lewis".
TEST(NrtTests, TheParentIsPrintedAsMatrixAndJsonBeforeTheResonanceSection)
{
    const NAOResult nao = h_chain(4);
    const double r = 1.0 / std::sqrt(2.0);
    dMatrix2 gamma(4, 4);
    for (const vec& v : { vec{ r, r, 0.0, 0.0 }, vec{ 0.0, 0.0, 1.0, 0.0 }, vec{ 0.0, 0.0, 0.0, 1.0 } })
        for (size_t i = 0; i < 4; i++)
            for (size_t j = 0; j < 4; j++)
                gamma(i, j) += 2.0 * v[i] * v[j];
    NboLewis lewis = one_bond(gamma, 4);
    for (const int a : { 2, 3 }) {
        NboFunction lp;
        lp.centers = { a };
        lp.type = "LP";
        lp.multiplicity = 1;
        lp.occupancy = 2.0;
        lewis.orbitals.push_back(lp);
        lewis.topo[a][a] = 1;
    }
    lewis.n_lewis = 3;
    NboOptions opt;
    opt.nrt = true;
    opt.nrt_exhaustive = true;
    NboResults res;
    res.source = "D:/somewhere/h4.wfn";
    std::ostringstream log;
    native_nrt(res.nrt, nao, lewis, {}, chain_bondable(4), opt, "", 2.0, log);
    ASSERT_TRUE(res.nrt.present);
    ASSERT_EQ(res.nrt.parent_topo.size(), 1u);
    EXPECT_EQ(res.nrt.parent_topo[0].matrix, (ivec2{ { 0, 1, 0, 0 }, { 1, 0, 0, 0 }, { 0, 0, 1, 0 }, { 0, 0, 0, 1 } }));

    std::ostringstream out;
    print_nrt(res, out);
    const std::string text = out.str();
    const std::string json = "{\n"
                             "  \"format\": \"NoSpherA2 Lewis structure\",\n"
                             "  \"version\": 1,\n"
                             "  \"wavefunction\": \"h4.wfn\",\n"
                             "  \"atoms\": 4,\n"
                             "  \"lone_pairs\": { \"H3\": 1, \"H4\": 1 },\n"
                             "  \"bonds\": [ [\"H1\", \"H2\", 1] ]\n"
                             "}";
    const size_t matrix = text.find("\nParent Lewis structure:\n");
    const size_t block = text.find("\nParent Lewis structure as NRT input (JSON):\n" + json + "\n");
    const size_t nrt = text.find("\nNatural resonance theory (in house):");
    ASSERT_NE(matrix, std::string::npos) << text;
    ASSERT_NE(block, std::string::npos) << text;
    EXPECT_LT(matrix, block);
    EXPECT_LT(block, nrt);
    EXPECT_NE(text.find("            H1    H2    H3    H4\n"
                        "  H1         .     1     .     .\n"
                        "  H2         1     .     .     .\n"
                        "  H3         .     .     1     .\n"
                        "  H4         .     .     .     1\n"), std::string::npos) << text;

    const std::filesystem::path file = std::filesystem::temp_directory_path() / "NrtTests_parent.native.nbo.json";
    write_nbo_json(res, file);
    std::ifstream in(file);
    const std::string stored((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    in.close();
    std::filesystem::remove(file);
    std::string nested = json;
    for (size_t p = nested.find('\n'); p != std::string::npos; p = nested.find('\n', p + 1)) nested.insert(p + 1, "  ");
    EXPECT_NE(stored.find("  },\n  \"parent_lewis\": " + nested + "\n}\n"), std::string::npos) << stored;
}

//An open shell carries one parent per spin, and past 40 atoms the matrix becomes one row per atom.
TEST(NrtTests, AnOpenShellParentSplitsBySpinAndALargeOneIsPrintedAsRows)
{
    NboResults res;
    res.source = "chain.gbw";
    res.nrt.present = true;
    const int na = 41;
    ivec2 m(na, ivec(na, 0));
    for (int a = 0; a + 1 < na; a++) m[a][a + 1] = m[a + 1][a] = 1;
    m[0][0] = 3;
    res.nrt.parent_topo = { NboTopo{ "alpha", m }, NboTopo{ "beta", m } };
    for (int a = 0; a < na; a++) {
        NboAtomPopulation p;
        p.element = a ? "C" : "Cl";
        p.index = a + 1;
        res.npa.push_back(p);
    }
    std::ostringstream out;
    print_nrt(res, out);
    const std::string text = out.str();
    EXPECT_NE(text.find("\nParent Lewis structure (alpha spin):\n"), std::string::npos) << text;
    EXPECT_NE(text.find("\nParent Lewis structure (beta spin):\n"), std::string::npos) << text;
    EXPECT_NE(text.find("  Cl1             3   C2:1\n"), std::string::npos) << text;
    EXPECT_NE(text.find("  C2              0   Cl1:1 C3:1\n"), std::string::npos) << text;

    const std::string json = parent_lewis_json(res);
    EXPECT_EQ(json.rfind("{\n  \"format\": \"NoSpherA2 Lewis structure\",\n  \"version\": 1,\n"
                         "  \"wavefunction\": \"chain.gbw\",\n  \"atoms\": 41,\n  \"alpha\": {\n"
                         "    \"lone_pairs\": { \"Cl1\": 3 },\n"
                         "    \"bonds\": [ [\"Cl1\", \"C2\", 1], [\"C2\", \"C3\", 1], ", 0), 0u) << json;
    //eight bonds to a line
    EXPECT_NE(json.find("[\"C8\", \"C9\", 1],\n      [\"C9\", \"C10\", 1]"), std::string::npos) << json;
    EXPECT_NE(json.find("[\"C40\", \"C41\", 1] ]\n  },\n  \"beta\": {\n"), std::string::npos) << json;
    EXPECT_EQ(json.substr(json.size() - 7), "]\n  }\n}");
}

//A hexagon of one-NAO atoms whose density lies halfway between its two Kekule structures, with no E2 to walk on:
//the second structure is only reachable as the parent's image under the hexagon's rotation, and the two then fit
//the density exactly with half a bond on every ring pair.
TEST(NrtTests, TheParentsSymmetryImageCompletesTheHexagon)
{
    const int na = 6;
    const NAOResult nao = h_chain(na);
    bvec2 ring(na, bvec(na, false));
    vec2 dist(na, vec(na, 0.0));
    dMatrix2 gamma(na, na);
    for (int a = 0; a < na; a++) {
        const int b = (a + 1) % na;
        ring[a][b] = ring[b][a] = true;
        gamma(a, a) += 0.5;
        gamma(b, b) += 0.5;
        gamma(a, b) += 0.5;
        gamma(b, a) += 0.5;
        for (int c = 0; c < na; c++)
            dist[a][c] = 2.0 * std::abs(std::sin(constants::PI * (a - c) / na));
    }
    NboLewis L;
    L.gamma = gamma;
    L.topo.assign(na, ivec(na, 0));
    for (int a = 0; a < na; a += 2) {
        NboFunction f;
        f.centers = { a, a + 1 };
        f.type = "BD";
        f.multiplicity = 1;
        f.occupancy = 1.5;
        L.orbitals.push_back(f);
        L.topo[a][a + 1] = L.topo[a + 1][a] = 1;
    }
    L.n_lewis = 3;
    NboOptions opt;
    opt.nrt = true;
    NboNrt nrt;
    std::ostringstream log;
    native_nrt(nrt, nao, L, {}, ring, opt, "", 2.0, log, &dist);
    ASSERT_TRUE(nrt.present);
    EXPECT_NEAR(nrt.d_w, 0.0, 1e-8) << "the second Kekule structure is missing";
    for (int a = 0; a < na; a++)
        EXPECT_NEAR(bond_total(nrt, a + 1, (a + 1) % na + 1), 0.5, 1e-6) << "ring pair " << a + 1;
}

//A double bond whose two orbitals lean to opposite atoms, i = c_A^2 - c_B^2 = +0.6 and -0.6, like Ni(CO)4's donor
//and back-donor Ni-C orbitals: the ionic share is the weighted |i| of each orbital, 0.6 of the bond order, not
//|0.6 - 0.6| = 0. Occupancies 2 and 1.8 keep the two orbitals apart.
TEST(NrtTests, OppositelyPolarisedOrbitalsOfABondDoNotCancelInItsIonicShare)
{
    const NAOResult nao = h_chain_shells(2, 2);
    const vec v1 = { std::sqrt(0.8), 0.0, std::sqrt(0.2), 0.0 };
    const vec v2 = { 0.0, std::sqrt(0.2), 0.0, std::sqrt(0.8) };
    dMatrix2 gamma(4, 4);
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            gamma(i, j) = 2.0 * v1[i] * v1[j] + 1.8 * v2[i] * v2[j];
    NboLewis L;
    L.gamma = gamma;
    NboFunction f;
    f.centers = { 0, 1 };
    f.type = "BD";
    f.occupancy = 2.0;
    f.multiplicity = 1;
    L.orbitals.push_back(f);
    f.multiplicity = 2;
    f.occupancy = 1.8;
    L.orbitals.push_back(f);
    L.n_lewis = 2;
    L.topo = { { 0, 2 }, { 2, 0 } };
    NboOptions opt;
    opt.nrt = true;
    NboNrt nrt;
    std::ostringstream log;
    native_nrt(nrt, nao, L, {}, chain_bondable(2), opt, "", 2.0, log);
    ASSERT_TRUE(nrt.present);
    for (const NboBondOrder& o : nrt.bond_orders)
        if (!o.diagonal) {
            EXPECT_NEAR(o.total, 2.0, 1e-12);
            EXPECT_NEAR(o.ionic, 1.2, 1e-6);
        }
}
