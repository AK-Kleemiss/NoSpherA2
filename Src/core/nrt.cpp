#include "pch.h"
#include "nbo.h"
#include "constants.h"

#include <Eigen/Dense>
#include <chrono>
#include <functional>
#include <map>

using Eigen::MatrixXd;
using Eigen::VectorXd;

//Natural resonance theory, after Glendening and Weinhold, J. Comput. Chem. 19 (1998) 593, as what it
//mathematically is: a convex quadratic program on the probability simplex,
//
//    minimise  D(w) = || Gamma - sum_a w_a Gamma_a ||_F   subject to  w >= 0, sum w = 1,
//
//over a set of candidate resonance structures.  Gamma is the density in the NAO basis and Gamma_a is
//the idealised density of candidate a, s V_a V_a^T over the orthonormal orbital set V_a that its
//integer topology prescribes (s = 2 closed shell, 1 per spin of an open one).
//
//Nothing here ever materialises a Gamma_a.  Expanding the norm gives
//
//    D(w)^2 = Tr(Gamma^2) - 2 g.w + w^T G w,   g_a = s Tr(V_a^T Gamma V_a),
//                                              G_ab = s^2 || V_a^T V_b ||_F^2,
//
//so one n x k matrix per candidate and one k x k product per pair is the whole cost, and the pair
//loop is embarrassingly parallel.  NBO 7 is serial - its binaries carry no OpenMP or pthread symbols
//- so this is where the speed comes from, not from a cleverer minimiser.
//
//Two honest deviations from NBO, both measured and both reported:
//
//  * the value of D.  NBO's printed D is scaled by something this code does not reproduce (acetylene
//    prints D(0) = 0.01884235 where the Frobenius norm of the same difference is 0.10108, a ratio of
//    5.37 that is not the electron count, the orbital count or any power of either).  D here is the
//    plain Frobenius norm defined above.  It is a residual, so its absolute value carries no physics;
//    what is compared against the reference is the weights by rank, the bond orders, the valencies
//    and the retained-structure count, all of which are D-normalisation independent.
//  * the argmin.  The candidate Gramians are strongly near-collinear - two resonance structures that
//    differ by one arrow overlap in almost every orbital - so the residual minimum is unique but the
//    weight vector attaining it need not be.  This is not an implementation artefact: gennbo itself
//    reproduces D(w), every bond order and every valency of TiCl4 to all eight printed decimals at
//    NRTE2 = 10 and at 20 while individual weights move by 7.1 percentage points.  Weights are
//    therefore compared by rank and never by label.

namespace
{
    MatrixXd to_eigen(const dMatrix2& m)
    {
        const int r = static_cast<int>(m.extent(0)), c = static_cast<int>(m.extent(1));
        MatrixXd out(r, c);
        for (int i = 0; i < r; i++)
            for (int j = 0; j < c; j++)
                out(i, j) = m(i, j);
        return out;
    }

    //Leading eigenpair of the block of R over idx, embedded back into the full space.  A copy of
    //nbo.cpp's helper; ten lines duplicated is cheaper than a header two files under parallel
    //development both have to include.
    double leading_block(const MatrixXd& R, const ivec& idx, VectorXd& v)
    {
        const int k = static_cast<int>(idx.size());
        MatrixXd B(k, k);
        for (int i = 0; i < k; i++)
            for (int j = 0; j < k; j++)
                B(i, j) = R(idx[i], idx[j]);
        Eigen::SelfAdjointEigenSolver<MatrixXd> es(B);
        v = VectorXd::Zero(R.rows());
        for (int i = 0; i < k; i++) v(idx[i]) = es.eigenvectors()(i, k - 1);
        return es.eigenvalues()(k - 1);
    }

    //--------------------------------------------------------------------------------------
    // resonance structures
    //--------------------------------------------------------------------------------------

    //A resonance structure is an integer symmetric matrix: bond multiplicities off the diagonal,
    //lone pairs on it.  Cores are not on the diagonal - they are in every structure and so cannot
    //distinguish two of them.  That convention is the reference's: acetylene's leading_topo prints 0
    //on both carbons while both carry a 1s core pair.
    struct Topology {
        int n = 0;
        ivec t;

        Topology() = default;
        explicit Topology(const int na) : n(na), t(static_cast<size_t>(na) * na, 0) {}
        int& at(const int a, const int b) { return t[static_cast<size_t>(a) * n + b]; }
        int at(const int a, const int b) const { return t[static_cast<size_t>(a) * n + b]; }
        void set(const int a, const int b, const int v) { at(a, b) = v; at(b, a) = v; }

        //orbitals the atom carries: lone pairs plus one per bond, a double bond counting twice
        int used(const int a) const
        {
            int u = at(a, a);
            for (int b = 0; b < n; b++) if (b != a) u += at(a, b);
            return u;
        }
        //electrons the Lewis structure assigns to the atom, both electrons of a shared bond counted
        int electrons(const int a) const
        {
            int e = 2 * at(a, a);
            for (int b = 0; b < n; b++) if (b != a) e += at(a, b);
            return e;
        }
        int pairs() const
        {
            int p = 0;
            for (int a = 0; a < n; a++)
                for (int b = a; b < n; b++) p += at(a, b);
            return p;
        }
        //upper triangle including the diagonal: the deduplication key
        ivec key() const
        {
            ivec k;
            k.reserve(static_cast<size_t>(n) * (n + 1) / 2);
            for (int a = 0; a < n; a++)
                for (int b = a; b < n; b++) k.push_back(at(a, b));
            return k;
        }
        ivec2 matrix() const
        {
            ivec2 m(n, ivec(n, 0));
            for (int a = 0; a < n; a++)
                for (int b = 0; b < n; b++) m[a][b] = at(a, b);
            return m;
        }
    };

    //Invariant under a relabelling of equivalent atoms: the sorted multiset of per-atom signatures
    //(element, lone pairs, sorted neighbour element/multiplicity list).  Two structures with the same
    //string are the same Lewis structure drawn on differently numbered atoms.
    std::string canonical_form(const Topology& s, const ivec& Z)
    {
        std::vector<std::string> sig(s.n);
        for (int a = 0; a < s.n; a++) {
            std::vector<std::string> nb;
            for (int b = 0; b < s.n; b++)
                if (b != a && s.at(a, b))
                    nb.push_back(std::string(constants::atnr2letter(Z[b])) +
                                 std::to_string(s.at(a, b)));
            std::sort(nb.begin(), nb.end());
            sig[a] = std::string(constants::atnr2letter(Z[a])) + "|" +
                     std::to_string(s.at(a, a)) + "|";
            for (const std::string& x : nb) sig[a] += x + ",";
        }
        std::sort(sig.begin(), sig.end());
        std::string out;
        for (const std::string& x : sig) out += x + ";";
        return out;
    }

    //What a candidate did to the parent, in the spirit of NBO's "Added(Removed)" column.
    std::string change_string(const Topology& s, const Topology& p, const ivec& Z)
    {
        const auto label = [&](const int a) {
            return std::string(constants::atnr2letter(Z[a])) + std::to_string(a + 1);
        };
        std::vector<std::string> add, rem;
        for (int a = 0; a < s.n; a++)
            for (int b = a; b < s.n; b++) {
                const int d = s.at(a, b) - p.at(a, b);
                if (!d) continue;
                const std::string what = (a == b) ? ("LP " + label(a))
                                                  : (label(a) + "-" + label(b));
                for (int k = 0; k < std::abs(d); k++) (d > 0 ? add : rem).push_back(what);
            }
        std::string out;
        for (size_t i = 0; i < add.size(); i++) out += (i ? ", " : "") + add[i];
        if (!rem.empty()) {
            out += "(";
            for (size_t i = 0; i < rem.size(); i++) out += (i ? ", " : "") + rem[i];
            out += ")";
        }
        return out.empty() ? std::string("parent") : out;
    }

    //--------------------------------------------------------------------------------------
    // a-priori screens
    //--------------------------------------------------------------------------------------

    //Which atom pairs the resonance search may move electrons between: exactly the pairs a
    //second-order interaction above nrt_e2_kcal connects.  This is screen (b) of the brief - E2 as the
    //pricing score - used as a hard gate rather than as a ranking, and it is the screen that decides
    //how many candidates exist at all.  Seeding it with the pairs the parent already bonds, which
    //looks harmless, is what made water generate nine ionic structures where NBO 7 generates one
    //structure and finds one: NBO's candidate count for a molecule with no delocalisation above the
    //threshold is one, and that is a consequence of this gate, not of a weight floor.
    bvec2 delocalisation_graph(const NboLewis& lewis, const std::vector<NboE2Entry>& e2,
                               const double kcal, const int na, const bvec2& bondable)
    {
        bvec2 g(na, bvec(na, false));
        for (const NboE2Entry& e : e2) {
            if (e.energy_kcal < kcal) continue;
            ivec c = lewis.orbitals[e.donor_index - 1].centers;
            const ivec& ac = lewis.orbitals[e.acceptor_index - 1].centers;
            c.insert(c.end(), ac.begin(), ac.end());
            for (size_t i = 0; i < c.size(); i++)
                for (size_t j = i + 1; j < c.size(); j++)
                    if (c[i] != c[j] && bondable[c[i]][c[j]])
                        g[c[i]][c[j]] = g[c[j]][c[i]] = true;
        }
        return g;
    }

    //Screen (a): connected components of that graph.  Two moves in different components describe
    //independent resonance, whose joint weight factorises, so the product of the two candidate sets
    //is combinatorial waste.
    ivec components(const bvec2& g)
    {
        const int n = static_cast<int>(g.size());
        ivec comp(n, -1);
        int c = 0;
        for (int a = 0; a < n; a++) {
            if (comp[a] >= 0) continue;
            ivec stack{ a };
            comp[a] = c;
            while (!stack.empty()) {
                const int x = stack.back();
                stack.pop_back();
                for (int y = 0; y < n; y++)
                    if (g[x][y] && comp[y] < 0) { comp[y] = c; stack.push_back(y); }
            }
            c++;
        }
        return comp;
    }

    //--------------------------------------------------------------------------------------
    // candidate generation
    //--------------------------------------------------------------------------------------

    struct Candidate {
        Topology topo;
        int depth = 0;
        int component = -1;       //which component of the delocalisation graph the moves touched
        //filled by the orbital construction
        bool feasible = false;
        double g = 0.0;           //s Tr(V^T Gamma V)
        double rho_nl = 0.0;      //electrons the structure leaves outside its own orbital set
        MatrixXd V;
        std::vector<std::array<double, 3>> polarity;  //{a, b, c_a^2 - c_b^2} per two-centre orbital
    };

    struct Limits {
        ivec cap;            //orbitals an atom may carry
        bvec free_atom;      //may a move touch this atom at all (subspace screen)
        bvec2 bondable;      //geometry screen
        bvec2 deloc;         //E2 screen
        ivec comp;           //component of each atom
        int max_charge = 2;  //how far a candidate may move an atom's Lewis electron count
        ivec parent_electrons;
        bool ion = true;
        bool use_components = true;
    };

    bool acceptable(const Topology& s, const Limits& L)
    {
        for (int a = 0; a < s.n; a++) {
            if (s.used(a) > L.cap[a]) return false;
            if (std::abs(s.electrons(a) - L.parent_electrons[a]) > L.max_charge) return false;
        }
        return true;
    }

    //One arrow: a single electron pair moves from one slot of the structure to another.  The
    //charge-neutral NRT arrow - a bond shifts and a lone pair takes its place - is the depth-two
    //composition of two of these, which is why nrt_max_arrows defaults to 2.
    void expand(const Candidate& c, const Limits& L, std::vector<Candidate>& out)
    {
        const int n = c.topo.n;
        const auto emit = [&](Topology s, const int x, const int y) {
            if (!L.free_atom[x] || !L.free_atom[y]) return;
            if (!acceptable(s, L)) return;
            Candidate d;
            d.topo = std::move(s);
            d.depth = c.depth + 1;
            const int comp = L.comp[x];
            if (L.use_components && c.component >= 0 && comp != c.component) return;
            d.component = comp;
            out.push_back(std::move(d));
        };
        for (int a = 0; a < n; a++) {
            //a lone pair becomes a bond, or moves to another atom
            if (c.topo.at(a, a) > 0) {
                for (int b = 0; b < n; b++) {
                    if (b == a || !L.deloc[a][b]) continue;
                    Topology s = c.topo;
                    s.at(a, a)--;
                    s.set(a, b, s.at(a, b) + 1);
                    emit(std::move(s), a, b);
                    if (L.ion) {
                        Topology i = c.topo;
                        i.at(a, a)--;
                        i.at(b, b)++;
                        emit(std::move(i), a, b);
                    }
                }
            }
            //a bond becomes a lone pair on either end, or shifts to a neighbouring pair
            for (int b = 0; b < n; b++) {
                if (b == a || c.topo.at(a, b) <= 0 || !L.deloc[a][b]) continue;
                if (L.ion) {
                    Topology s = c.topo;
                    s.set(a, b, s.at(a, b) - 1);
                    s.at(a, a)++;
                    emit(std::move(s), a, b);
                }
                for (int d = 0; d < n; d++) {
                    if (d == a || d == b || !L.deloc[a][d]) continue;
                    Topology s = c.topo;
                    s.set(a, b, s.at(a, b) - 1);
                    s.set(a, d, s.at(a, d) + 1);
                    emit(std::move(s), a, d);
                }
            }
        }
    }

    //Every topology the constraints admit, not only the ones an arrow walk reaches.  This is the
    //reference mode the screens are measured against; it is exponential in the number of bondable
    //pairs and is capped by nrt_max_candidates.
    void enumerate(const Topology& parent, const Limits& L, const int n_pairs, const size_t cap,
                   std::vector<Candidate>& out)
    {
        struct Slot { int a, b; };
        std::vector<Slot> slots;
        for (int a = 0; a < parent.n; a++) {
            if (L.free_atom[a]) slots.push_back({ a, a });
            for (int b = a + 1; b < parent.n; b++)
                if (L.bondable[a][b] && L.free_atom[a] && L.free_atom[b]) slots.push_back({ a, b });
        }
        //the slots a move may not touch keep their parent value and do not enter the recursion
        Topology fixed(parent.n);
        int fixed_pairs = 0;
        for (int a = 0; a < parent.n; a++)
            for (int b = a; b < parent.n; b++) {
                const bool movable = (a == b) ? L.free_atom[a]
                                              : (L.bondable[a][b] && L.free_atom[a] && L.free_atom[b]);
                if (!movable) { fixed.set(a, b, parent.at(a, b)); fixed_pairs += parent.at(a, b); }
            }
        Topology s = fixed;
        std::function<void(size_t, int)> rec = [&](const size_t i, const int left) {
            if (out.size() >= cap) return;
            if (i == slots.size()) {
                if (left == 0 && acceptable(s, L)) {
                    Candidate d;
                    d.topo = s;
                    d.depth = -1;
                    out.push_back(std::move(d));
                }
                return;
            }
            const int a = slots[i].a, b = slots[i].b;
            //an orbital count above 3 is not a Lewis structure of any element in the reference set
            for (int m = 0; m <= std::min(3, left); m++) {
                s.set(a, b, m);
                if (s.used(a) <= L.cap[a] && s.used(b) <= L.cap[b]) rec(i + 1, left - m);
                if (out.size() >= cap) break;
            }
            s.set(a, b, 0);
        };
        rec(0, n_pairs - fixed_pairs);
    }

    //--------------------------------------------------------------------------------------
    // the orbitals of one candidate
    //--------------------------------------------------------------------------------------

    MatrixXd sym_power(const MatrixXd& M, const double p, const double rel_floor = 1e-10)
    {
        Eigen::SelfAdjointEigenSolver<MatrixXd> es(M);
        VectorXd w = es.eigenvalues();
        const double cut = rel_floor * std::max(w.maxCoeff(), 1e-300);
        VectorXd f(w.size());
        for (int i = 0; i < w.size(); i++)
            f(i) = (w(i) > cut) ? std::pow(w(i), p) : 0.0;
        return es.eigenvectors() * f.asDiagonal() * es.eigenvectors().transpose();
    }

    struct Blocks {
        std::vector<ivec> atom;     //every NAO index of the atom, cores included
        ivec core;                  //core NAO indices, one fixed orbital each in every candidate
        std::map<std::pair<int, int>, ivec> pair;
        vec score_atom;             //leading eigenvalue of the undepleted block, the slot order
        std::map<std::pair<int, int>, double> score_pair;
    };

    /**
     * A candidate's orthonormal orbital set, built by the same three steps the NBO search uses -
     * greedy fill, self-consistency sweep, occupancy-weighted symmetric orthogonalisation - with the
     * topology prescribed instead of searched.  The sweep is not optional and it is not a refinement:
     * without it every candidate looks like a decent fit, because a greedy pass leaves each orbital
     * carrying the bias of the ones accepted after it.  On water the parent came out at rho_NL = 0.73
     * instead of the search's 0.02, and the residual was then so insensitive to the topology that
     * water's ionic structures took 23 % of the weight against the reference's 0.  With the sweep a
     * wrong topology cannot hide: it is the difference between measuring the structures and
     * measuring the construction.
     */
    void build_orbitals(Candidate& c, const MatrixXd& G0, const Blocks& B, const int max_sweeps)
    {
        struct Slot { int a, b, mult; double score; };
        std::vector<Slot> slots;
        for (int a = 0; a < c.topo.n; a++)
            if (c.topo.at(a, a) > 0) slots.push_back({ a, a, c.topo.at(a, a), B.score_atom[a] });
        for (int a = 0; a < c.topo.n; a++)
            for (int b = a + 1; b < c.topo.n; b++)
                if (c.topo.at(a, b) > 0) {
                    const auto it = B.score_pair.find({ a, b });
                    if (it == B.score_pair.end()) return;   //not a bondable pair: infeasible
                    slots.push_back({ a, b, c.topo.at(a, b), it->second });
                }
        std::sort(slots.begin(), slots.end(), [](const Slot& x, const Slot& y) {
            if (x.score != y.score) return x.score > y.score;
            return std::tie(x.a, x.b) < std::tie(y.a, y.b);
        });

        const int n = static_cast<int>(G0.rows());
        int k = static_cast<int>(B.core.size());
        for (const Slot& sl : slots) k += sl.mult;
        if (k > n) return;
        std::vector<VectorXd> v(k);
        std::vector<const ivec*> blk(k, nullptr);
        std::vector<std::pair<int, int>> owner(k, { -1, -1 });
        vec occ(k, 0.0);

        //1. the cores: one NAO each, in every candidate, never swept
        MatrixXd sum = MatrixXd::Zero(n, n);
        int col = 0;
        for (const int i : B.core) {
            v[col] = VectorXd::Zero(n);
            v[col](i) = 1.0;
            occ[col] = G0(i, i);
            sum(i, i) += occ[col];
            col++;
        }
        const int first_valence = col;

        //2. greedy fill out of the density the cores and the earlier orbitals have been taken from
        MatrixXd R = G0 - sum;
        for (const Slot& sl : slots) {
            const ivec& idx = (sl.a == sl.b) ? B.atom[sl.a] : B.pair.at({ sl.a, sl.b });
            if (static_cast<int>(idx.size()) < sl.mult) return;
            for (int m = 0; m < sl.mult; m++) {
                VectorXd x;
                const double lam = leading_block(R, idx, x);
                R -= lam * x * x.transpose();
                v[col] = x;
                blk[col] = &idx;
                owner[col] = { sl.a, sl.b };
                occ[col] = x.dot(G0 * x);
                sum += occ[col] * x * x.transpose();
                col++;
            }
        }

        //3. self consistency: every orbital against the density with all the others removed
        for (int s = 0; s < max_sweeps; s++) {
            double change = 0.0;
            for (int j = first_valence; j < k; j++) {
                sum -= occ[j] * v[j] * v[j].transpose();
                VectorXd x;
                leading_block(MatrixXd(G0 - sum), *blk[j], x);
                if (x.dot(v[j]) < 0.0) x = -x;
                change = std::max(change, (x - v[j]).norm());
                v[j] = x;
                occ[j] = x.dot(G0 * x);
                sum += occ[j] * x * x.transpose();
            }
            if (change < 1e-9) break;
        }

        //4. OWSO.  Two bonds at the same atom overlap by about a fifth; weighting by occupancy lets
        //the occupied orbitals keep their shape at the expense of the empty ones and keeps two
        //equivalent bonds equivalent.  A topology whose orbitals are linearly dependent - two lone
        //pairs asked of an atom whose block the bonds have already used up - drops out here.
        MatrixXd M(n, k);
        for (int j = 0; j < k; j++) M.col(j) = v[j];
        VectorXd wt(k);
        for (int j = 0; j < k; j++) wt(j) = std::max(occ[j], 1e-6);
        const MatrixXd S = M.transpose() * M;
        if (Eigen::SelfAdjointEigenSolver<MatrixXd>(S).eigenvalues()(0) < 1e-8) return;
        c.V = M * wt.asDiagonal() *
              sym_power(MatrixXd(wt.asDiagonal() * S * wt.asDiagonal()), -0.5);

        for (int j = first_valence; j < k; j++) {
            if (owner[j].first == owner[j].second) continue;
            const int a = owner[j].first, b = owner[j].second;
            double wa = 0.0, wb = 0.0;
            for (const int i : B.atom[a]) wa += c.V(i, j) * c.V(i, j);
            for (const int i : B.atom[b]) wb += c.V(i, j) * c.V(i, j);
            c.polarity.push_back({ static_cast<double>(a), static_cast<double>(b), wa - wb });
        }
        c.feasible = true;
    }

    //--------------------------------------------------------------------------------------
    // the quadratic program
    //--------------------------------------------------------------------------------------

    //Euclidean projection onto the probability simplex, Duchi et al., ICML 2008.
    void project_simplex(VectorXd& w)
    {
        const int n = static_cast<int>(w.size());
        if (!n) return;
        VectorXd u = w;
        std::sort(u.data(), u.data() + n, std::greater<double>());
        double css = 0.0, theta = 0.0;
        for (int i = 0; i < n; i++) {
            css += u(i);
            const double t = (css - 1.0) / (i + 1);
            if (u(i) - t > 0.0) theta = t;
        }
        for (int i = 0; i < n; i++) w(i) = std::max(0.0, w(i) - theta);
    }

    double objective(const MatrixXd& G, const VectorXd& g, const double trg2, const VectorXd& w)
    {
        return trg2 - 2.0 * g.dot(w) + w.dot(G * w);
    }

    //Accelerated projected gradient (FISTA with adaptive restart).  The Hessian 2G is badly
    //conditioned by construction - candidates that differ by one arrow are near-collinear - so this
    //runs to identify the support and the support problem is then solved again on its own.
    double solve_qp(const MatrixXd& G, const VectorXd& g, const double trg2, VectorXd& w,
                    const int maxit, const vec* rho, const std::string& spin,
                    std::vector<NboQpIteration>* trace)
    {
        const int n = static_cast<int>(G.rows());
        if (!n) return trg2;
        VectorXd z = VectorXd::Constant(n, 1.0 / std::sqrt(static_cast<double>(n)));
        double lam = G.diagonal().maxCoeff();
        for (int it = 0; it < 100; it++) {
            const VectorXd y = G * z;
            const double nn = y.norm();
            if (nn <= 0.0) break;
            lam = nn;
            z = y / nn;
        }
        const double L = std::max(2.0 * lam, 1e-12);
        VectorXd y = w, wp = w;
        double t = 1.0, f = objective(G, g, trg2, w);
        for (int it = 1; it <= maxit; it++) {
            VectorXd wn = y - (2.0 / L) * (G * y - g);
            project_simplex(wn);
            const double fn = objective(G, g, trg2, wn);
            if (fn > f) { y = wn; t = 1.0; }          //restart: the momentum overshot
            else {
                const double tn = 0.5 * (1.0 + std::sqrt(1.0 + 4.0 * t * t));
                y = wn + ((t - 1.0) / tn) * (wn - wp);
                t = tn;
            }
            const double step = (wn - wp).cwiseAbs().maxCoeff();
            const double df = std::abs(fn - f);
            wp = wn;
            f = fn;
            if (trace && (it <= 10 || it % 50 == 0)) {
                NboQpIteration q;
                q.iteration = it;
                q.structures = static_cast<int>((wn.array() > 1e-6).count());
                q.d_w = std::sqrt(std::max(0.0, fn));
                q.spin = spin;
                if (rho)
                    for (int i = 0; i < n; i++) q.rho_nl += wn(i) * (*rho)[i];
                trace->push_back(q);
            }
            if (df < 1e-14 * std::max(1.0, std::abs(f)) && step < 1e-12) break;
        }
        w = wp;
        return f;
    }
}

void native_nrt(NboNrt& nrt, const NAOResult& nao, const NboLewis& lewis,
                const std::vector<NboE2Entry>& e2, const bvec2& bondable,
                const NboOptions& options, const std::string& spin, const double scale,
                std::ostream& log)
{
    const auto clock = [] { return std::chrono::steady_clock::now(); };
    const auto secs = [](const std::chrono::steady_clock::time_point a,
                         const std::chrono::steady_clock::time_point b) {
        return std::chrono::duration<double>(b - a).count();
    };
    const auto t_start = clock();

    const int na = static_cast<int>(nao.atoms.size());
    const int nn = static_cast<int>(nao.orbitals.size());
    ivec Z(na, 0);
    for (int a = 0; a < na; a++) Z[a] = nao.atoms[a].Z;

    //--- the parent structure, from the NBO search -------------------------------------------------
    ivec ncore(na, 0);
    for (const NboFunction& f : lewis.orbitals)
        if (f.type == "CR") ncore[f.centers[0]]++;
    Topology parent(na);
    for (int a = 0; a < na; a++) {
        parent.at(a, a) = lewis.topo[a][a] - ncore[a];
        for (int b = 0; b < na; b++)
            if (b != a) parent.at(a, b) = lewis.topo[a][b];
    }
    const int n_pairs = parent.pairs();

    Limits L;
    L.bondable = bondable;  //the resonance screen, nrt_bond_scale, not the search's (see nbo.h)
    L.ion = options.nrt_ion;
    L.use_components = options.nrt_components;
    L.deloc = delocalisation_graph(lewis, e2, options.nrt_e2_kcal, na, L.bondable);
    L.comp = components(L.deloc);
    L.cap.assign(na, 0);
    L.parent_electrons.assign(na, 0);
    L.free_atom.assign(na, options.nrt_subspace.empty());
    for (const int a : options.nrt_subspace)
        if (a >= 1 && a <= na) L.free_atom[a - 1] = true;
    for (int a = 0; a < na; a++) L.parent_electrons[a] = parent.electrons(a);
    {
        ivec nval(na, 0);
        for (const NAO& o : nao.orbitals)
            if (o.type == NAOClass::Valence) nval[o.atom]++;
        //An atom carries at most as many orbitals as it has valence NAOs, but never fewer than the
        //parent already gives it: the NBO search relaxes that cap for a genuinely hypervalent
        //density, and a resonance search that could not reproduce its own parent is nonsense.
        for (int a = 0; a < na; a++) L.cap[a] = std::max(nval[a], parent.used(a));
    }

    //--- candidates --------------------------------------------------------------------------------
    std::vector<Candidate> cands;
    std::map<ivec, int> seen;
    {
        Candidate p;
        p.topo = parent;
        cands.push_back(p);
        seen[parent.key()] = 0;
    }
    const auto t_search0 = clock();
    if (options.nrt_exhaustive) {
        std::vector<Candidate> all;
        enumerate(parent, L, n_pairs, static_cast<size_t>(options.nrt_max_candidates), all);
        for (Candidate& c : all)
            if (seen.emplace(c.topo.key(), static_cast<int>(cands.size())).second)
                cands.push_back(std::move(c));
        nrt.arrows.push_back("exhaustive enumeration yields " + std::to_string(cands.size()) +
                             " feasible topologies");
    }
    else {
        size_t level_begin = 0;
        for (int depth = 1; depth <= options.nrt_max_arrows; depth++) {
            const size_t level_end = cands.size();
            std::vector<Candidate> made;
            for (size_t i = level_begin; i < level_end; i++) expand(cands[i], L, made);
            size_t added = 0;
            for (Candidate& c : made) {
                if (cands.size() >= static_cast<size_t>(options.nrt_max_candidates)) break;
                if (seen.emplace(c.topo.key(), static_cast<int>(cands.size())).second) {
                    cands.push_back(std::move(c));
                    added++;
                }
            }
            nrt.arrows.push_back("ARROWS depth " + std::to_string(depth) + " generates " +
                                 std::to_string(added) + " new structures from " +
                                 std::to_string(level_end - level_begin));
            level_begin = level_end;
            if (!added) break;
        }
        //A single pair-move is never a resonance structure of its own, only half of one.  Count the
        //entries in the "changes" column of all 588 structures the 22 references list: 78 carry none
        //(they are the reference structures themselves), 239 carry four, 52 six, 193 eight, 23 ten,
        //and not one carries two.  Four entries is one add and one removal per move for two moves, so
        //NBO's elementary step is a coupled pair - a bond shifts while a lone pair takes its place -
        //and the intermediate is not a candidate.  Keeping the depth-one half-moves is what made LiF
        //bond its ionic parent back together with 99.65 % weight and water invent ionic structures.
        const size_t before = cands.size();
        cands.erase(std::remove_if(cands.begin() + 1, cands.end(),
                                   [](const Candidate& c) { return c.depth == 1; }),
                    cands.end());
        nrt.arrows.push_back("half-arrow intermediates dropped: " +
                             std::to_string(before - cands.size()) + " of " +
                             std::to_string(before));
    }
    const double search_seconds = secs(t_search0, clock());

    //--- the orbital sets and the Gram matrix ------------------------------------------------------
    const auto t_gram0 = clock();
    const MatrixXd gamma = to_eigen(lewis.gamma);
    double electrons = 0.0;
    for (int i = 0; i < nn; i++) electrons += gamma(i, i);

    Blocks B;
    B.atom.assign(na, ivec());
    for (int i = 0; i < nn; i++) {
        B.atom[nao.orbitals[i].atom].push_back(i);
        if (nao.orbitals[i].type == NAOClass::Core) B.core.push_back(i);
    }
    //the slot order is candidate independent: the leading eigenvalue of the slot's own block of the
    //density with the cores taken out, so a core cannot be mistaken for a lone pair
    MatrixXd g_nocore = gamma;
    for (const int i : B.core) g_nocore(i, i) -= gamma(i, i);
    B.score_atom.assign(na, 0.0);
    for (int a = 0; a < na; a++) {
        if (B.atom[a].empty()) continue;
        VectorXd v;
        B.score_atom[a] = leading_block(g_nocore, B.atom[a], v);
    }
    for (int a = 0; a < na; a++)
        for (int b = a + 1; b < na; b++) {
            if (!bondable[a][b] || B.atom[a].empty() || B.atom[b].empty()) continue;
            ivec idx = B.atom[a];
            idx.insert(idx.end(), B.atom[b].begin(), B.atom[b].end());
            VectorXd v;
            B.score_pair[{ a, b }] = leading_block(g_nocore, idx, v);
            B.pair[{ a, b }] = std::move(idx);
        }

#ifdef _OPENMP
    const int nthreads = options.threads > 0 ? options.threads : omp_get_max_threads();
#else
    const int nthreads = 1;
#endif
    const int nc0 = static_cast<int>(cands.size());
#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int i = 0; i < nc0; i++) {
        build_orbitals(cands[i], gamma, B, 50);
        if (!cands[i].feasible) continue;
        const MatrixXd GV = gamma * cands[i].V;
        double tr = 0.0;
        for (int c = 0; c < cands[i].V.cols(); c++) tr += cands[i].V.col(c).dot(GV.col(c));
        cands[i].g = scale * tr;
        cands[i].rho_nl = electrons - tr;
    }

    //drop the infeasible ones, keeping the parent first
    std::vector<Candidate> keep;
    for (Candidate& c : cands)
        if (c.feasible) keep.push_back(std::move(c));
    err_checkf(!keep.empty() && keep[0].topo.key() == parent.key(),
               "NRT: the parent Lewis structure is not representable in its own NAO basis", log);
    cands = std::move(keep);
    const int nc = static_cast<int>(cands.size());

    MatrixXd G(nc, nc);
    VectorXd g(nc);
    vec rho(nc, 0.0);
    for (int i = 0; i < nc; i++) { g(i) = cands[i].g; rho[i] = cands[i].rho_nl; }
#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int i = 0; i < nc; i++)
        for (int j = i; j < nc; j++) {
            const double v = scale * scale *
                             (cands[i].V.transpose() * cands[j].V).squaredNorm();
            G(i, j) = v;
            G(j, i) = v;
        }
    const double trg2 = gamma.squaredNorm();
    const double gram_seconds = secs(t_gram0, clock());

    //--- minimise ---------------------------------------------------------------------------------
    const auto t_min0 = clock();
    VectorXd w = VectorXd::Zero(nc);
    w(0) = 1.0;                       //start at the parent vertex
    const double f0 = objective(G, g, trg2, w);
    double f = solve_qp(G, g, trg2, w, 5000, &rho, spin, &nrt.qp_iterations);
    {
        NboNrtCycle cy;
        cy.cycle = 1;
        cy.spin = spin;
        cy.structures_found = nc;
        cy.structures_used = static_cast<int>((w.array() > options.nrt_weight_floor).count());
        cy.d_w = std::sqrt(std::max(0.0, f));
        nrt.cycles.push_back(cy);
    }
    //Polish: drop the structures below the reporting floor and solve again on the support alone.  The
    //dropped weights are not noise-free - the argmin is not unique - but the residual is, and this is
    //what makes the reported weights a solution of the problem that is actually reported.
    ivec support;
    for (int i = 0; i < nc; i++)
        if (w(i) > options.nrt_weight_floor) support.push_back(i);
    if (!support.empty() && static_cast<int>(support.size()) < nc) {
        const int ns = static_cast<int>(support.size());
        MatrixXd Gs(ns, ns);
        VectorXd gs(ns), ws(ns);
        vec rs(ns, 0.0);
        for (int i = 0; i < ns; i++) {
            gs(i) = g(support[i]);
            rs[i] = rho[support[i]];
            ws(i) = w(support[i]);
            for (int j = 0; j < ns; j++) Gs(i, j) = G(support[i], support[j]);
        }
        project_simplex(ws);
        const double fs = solve_qp(Gs, gs, trg2, ws, 200000, &rs, spin, nullptr);
        if (fs <= f + 1e-12) {
            w.setZero();
            for (int i = 0; i < ns; i++) w(support[i]) = ws(i);
            f = fs;
        }
        NboNrtCycle cy;
        cy.cycle = 2;
        cy.spin = spin;
        cy.structures_found = ns;
        cy.structures_used = static_cast<int>((w.array() > options.nrt_weight_floor).count());
        cy.d_w = std::sqrt(std::max(0.0, f));
        nrt.cycles.push_back(cy);
    }

    //Screen (d), applied after the fact: average the weights over classes of structures that are the
    //same Lewis structure on differently numbered atoms.  The argmin set of a convex function is
    //convex, so averaging over members that are all optimal stays optimal - and if the class is not a
    //real symmetry of the density the average is worse, which is checked and rejected.
    int orbits = 0, in_orbits = 0;
    if (options.nrt_symmetry) {
        //Every candidate enters an orbit, not only the ones above the weight floor, and the orbit
        //weight is not averaged after the fact but solved for: substituting w_i = u_p / |p| into the
        //objective gives the same quadratic program over the orbits, with G and g block-averaged, and
        //its solution is the symmetric optimum - the other structures readjust, which a post-hoc
        //average does not let them do.  Averaging instead kept formate's two equivalent structures at
        //47.98 and 21.62 where the reference has 33.29 and 33.29, and the rise it cost (0.029 in D)
        //then looked like evidence against the symmetry rather than against the method.
        //The canonical form is a graph invariant, and a graph invariant is not a symmetry of the
        //molecule: ozone's ring structure, O 1- O 3 across 2.24 A, is isomorphic to the open one with
        //its 1.28 A bond, and forcing the two to equal weight is simply wrong.  Splitting each
        //canonical class by Tr(V^T Gamma V) - equal to 1e-5, single-link - is the numerical
        //"close enough" that tells a real symmetry from an accidental isomorphism, and it needs no
        //coordinates: two structures that fit this density equally well and have the same graph are
        //the same structure on renumbered atoms.
        std::map<std::string, ivec> classes;
        for (int i = 0; i < nc; i++)
            if (cands[i].feasible) classes[canonical_form(cands[i].topo, Z)].push_back(i);
        std::vector<ivec> member;
        for (auto& kv : classes) {
            ivec& m = kv.second;
            std::sort(m.begin(), m.end(),
                      [&](const int a, const int b) { return cands[a].g < cands[b].g; });
            for (size_t k = 0; k < m.size(); k++) {
                if (k && cands[m[k]].g - cands[m[k - 1]].g <= 1e-5)
                    member.back().push_back(m[k]);
                else
                    member.push_back(ivec{ m[k] });
            }
        }
        const int no = static_cast<int>(member.size());
        for (const ivec& m : member)
            if (m.size() > 1) { orbits++; in_orbits += static_cast<int>(m.size()); }
        if (orbits) {
            MatrixXd Gr(no, no);
            VectorXd gr(no);
            for (int p = 0; p < no; p++) {
                double s = 0.0;
                for (const int i : member[p]) s += g(i);
                gr(p) = s / static_cast<double>(member[p].size());
                for (int q = 0; q < no; q++) {
                    double t = 0.0;
                    for (const int i : member[p])
                        for (const int j : member[q]) t += G(i, j);
                    Gr(p, q) = t / static_cast<double>(member[p].size() * member[q].size());
                }
            }
            VectorXd u = VectorXd::Constant(no, 1.0 / static_cast<double>(no));
            const double fu = solve_qp(Gr, gr, trg2, u, 200000, nullptr, spin, nullptr);
            const double rise = std::sqrt(std::max(0.0, fu)) - std::sqrt(std::max(0.0, f));
            //A rise means the orbits are not a symmetry of this density after all; 1 % of D(w) is the
            //"close enough" the brief asks for, and the number itself is reported either way.
            if (rise <= 0.01 * std::sqrt(std::max(1e-12, f))) {
                w.setZero();
                for (int p = 0; p < no; p++)
                    for (const int i : member[p]) w(i) = u(p) / static_cast<double>(member[p].size());
                f = fu;
                nrt.symmetry = std::to_string(orbits) + " graph-invariant orbit(s) covering " +
                               std::to_string(in_orbits) + " structure(s), " + std::to_string(no) +
                               " orbit weights solved for, D(w) rise " + std::to_string(rise);
            }
            else {
                nrt.symmetry = "graph-invariant orbit constraint rejected, it raised D(w) by " +
                               std::to_string(rise);
                orbits = 0;
            }
        }
    }
    if (nrt.symmetry.empty())
        nrt.symmetry = options.nrt_symmetry
            ? (std::to_string(orbits) + " graph-invariant orbit(s) covering " +
               std::to_string(in_orbits) + " structure(s), weights averaged inside each")
            : "symmetry screen off";
    const double minimize_seconds = secs(t_min0, clock());

    //--- what the weights imply -------------------------------------------------------------------
    const auto t_other0 = clock();
    vec2 bo(na, vec(na, 0.0));        //bond orders, diagonal = lone pairs
    vec2 pol(na, vec(na, 0.0));       //weight-summed c_A^2 - c_B^2 of the bonds on the pair
    for (int i = 0; i < nc; i++) {
        if (w(i) <= 0.0) continue;
        for (int a = 0; a < na; a++)
            for (int b = a; b < na; b++)
                bo[a][b] += w(i) * cands[i].topo.at(a, b);
        for (const std::array<double, 3>& p : cands[i].polarity) {
            const int a = static_cast<int>(p[0]), b = static_cast<int>(p[1]);
            pol[std::min(a, b)][std::max(a, b)] += w(i) * ((a < b) ? p[2] : -p[2]);
        }
    }

    //The ionic share of a bond order is |i| b with i = c_A^2 - c_B^2, linear in |i| and not in i^2.
    //Acetylene's printed tables settle it: C-H carries i = 0.2334 and the table gives
    //ionic/total = 0.2262/0.9693 = 0.23336.
    vec valency(na, 0.0), covalency(na, 0.0), electrovalency(na, 0.0);
    for (int a = 0; a < na; a++)
        for (int b = a + 1; b < na; b++) {
            const double t = bo[a][b];
            if (t <= 0.0) continue;
            const double ion = std::min(1.0, std::abs(pol[a][b]) / std::max(t, 1e-30));
            NboBondOrder o;
            o.atom1 = a + 1;
            o.atom2 = b + 1;
            o.spin = spin;
            o.total = t;
            o.ionic = ion * t;
            o.covalent = t - o.ionic;
            nrt.bond_orders.push_back(o);
            valency[a] += t;            valency[b] += t;
            covalency[a] += o.covalent; covalency[b] += o.covalent;
            electrovalency[a] += o.ionic; electrovalency[b] += o.ionic;
        }
    for (int a = 0; a < na; a++) {
        NboBondOrder o;
        o.atom1 = o.atom2 = a + 1;
        o.diagonal = true;
        o.spin = spin;
        o.total = bo[a][a];
        nrt.bond_orders.push_back(o);
        NboValency v;
        v.atom = a + 1;
        v.element = constants::atnr2letter(Z[a]);
        v.spin = spin;
        v.valency = valency[a];
        v.covalency = covalency[a];
        v.electrovalency = electrovalency[a];
        //the octet count: both electrons of a shared bond counted for each partner, cores excluded.
        //Acetylene's carbon prints 7.9657 and 2 (0.0198 + 3.9631) = 7.9658.
        v.electron_count = 2.0 * (bo[a][a] + valency[a]);
        nrt.valencies.push_back(v);
    }

    ivec order;
    for (int i = 0; i < nc; i++)
        if (w(i) > options.nrt_weight_floor) order.push_back(i);
    std::stable_sort(order.begin(), order.end(), [&](const int a, const int b) {
        return w(a) > w(b);
    });
    for (size_t r = 0; r < order.size(); r++) {
        const int i = order[r];
        NboResonanceWeight x;
        x.structure = static_cast<int>(r) + 1;
        x.idxres = i + 1;
        x.spin = spin;
        x.weight_percent = 100.0 * w(i);
        x.weight_fraction = w(i);
        x.changes = change_string(cands[i].topo, parent, Z);
        nrt.weights.push_back(x);
        NboNrtCandidate c;
        c.structure = x.structure;
        c.idxres = x.idxres;
        c.rho_nl = cands[i].rho_nl;
        c.spin = spin;
        c.topo = cands[i].topo.matrix();
        nrt.candidates.push_back(c);
    }
    //Every structure that ties the leading weight, not just the first of them: ozone's two Lewis
    //structures both come out at 32.96 % and which one is printed first is arbitrary - in NBO too,
    //whose own leading pair is 25.03 % twice.  A single printed matrix invites a comparison that
    //mistakes the tie-break for a disagreement.
    for (size_t r = 0; r < order.size(); r++) {
        if (w(order[r]) < w(order[0]) - 5.0e-4) break;
        NboTopo t;
        t.spin = spin;
        t.matrix = cands[order[r]].topo.matrix();
        nrt.leading_topo.push_back(t);
    }

    nrt.present = true;
    nrt.structures_used += static_cast<int>(order.size());
    nrt.structures_found += nc;
    nrt.d_w = std::sqrt(std::max(0.0, f));
    nrt.d_0 = std::sqrt(std::max(0.0, f0));
    nrt.deloc_threshold_kcal = options.nrt_e2_kcal;
    nrt.max_search_cycles = static_cast<int>(nrt.cycles.size());
    nrt.initial_topo = nc;
    nrt.search_seconds += search_seconds;
    nrt.gram_seconds += gram_seconds;
    nrt.minimize_seconds += minimize_seconds;
    nrt.other_seconds += secs(t_other0, clock()) + secs(t_start, t_search0);
    if (options.debug)
        log << "NRT" << (spin.empty() ? "" : " " + spin) << ": " << nc << " candidates, "
            << order.size() << " retained, D(0) = " << nrt.d_0 << ", D(w) = " << nrt.d_w
            << " (search " << search_seconds << " s, gram " << gram_seconds << " s, minimise "
            << minimize_seconds << " s)\n";
}
