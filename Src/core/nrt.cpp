#include "pch.h"
#include "nbo.h"
#include "constants.h"

#include <Eigen/Dense>
#include <chrono>
#include <functional>
#include <limits>
#include <map>

using Eigen::MatrixXd;
using Eigen::VectorXd;

//NRT minimizes ||Gamma - sum_a w_a Gamma_a||_F over simplex weights w.

namespace
{
	//Per-slot candidate budget
	constexpr int NRT_PER_SLOT = 64;

	MatrixXd to_eigen(const dMatrix2& m)
	{
		const int r = static_cast<int>(m.extent(0)), c = static_cast<int>(m.extent(1));
		MatrixXd out(r, c);
		for (int i = 0; i < r; i++)
			for (int j = 0; j < c; j++)
				out(i, j) = m(i, j);
		return out;
	}

	//Leading eigenpair on block idx by warm-started power iteration; a negative dominant eigenvalue
	//falls back to the full symmetric solve
	double leading_block(const MatrixXd& R, const ivec& idx, VectorXd& v,
						 const MatrixXd* minus = nullptr, const VectorXd* warm = nullptr)
	{
		const int k = static_cast<int>(idx.size());
		MatrixXd B(k, k);
		for (int i = 0; i < k; i++)
			for (int j = 0; j < k; j++)
				B(i, j) = minus ? R(idx[i], idx[j]) - (*minus)(idx[i], idx[j])
								: R(idx[i], idx[j]);

		const auto embed = [&](const VectorXd& x) {
			v = VectorXd::Zero(R.rows());
			for (int i = 0; i < k; i++) v(idx[i]) = x(i);
		};

		if (warm) {
			VectorXd x(k);
			for (int i = 0; i < k; i++) x(i) = (*warm)(idx[i]);
			const double nx = x.norm();
			if (nx > 0.0) {
				x /= nx;
				for (int it = 0; it < 30; it++) {
					VectorXd y = B * x;
					const double ny = y.norm();
					if (ny == 0.0) break;
					y /= ny;
					//Up to sign: a dominant negative eigenvalue flips y each step and must reach the fallback
					const double d = std::min((y - x).norm(), (y + x).norm());
					x = y;
					if (d < 1e-12) {
						const double lam = x.dot(B * x);
						if (lam > 0.0) { embed(x); return lam; }
						break;
					}
				}
			}
		}

		Eigen::SelfAdjointEigenSolver<MatrixXd> es(B);
		embed(es.eigenvectors().col(k - 1));
		return es.eigenvalues()(k - 1);
	}

	//x^T A x on the nonzero block idx
	double block_quad(const MatrixXd& A, const ivec& idx, const VectorXd& x)
	{
		const int k = static_cast<int>(idx.size());
		double s = 0.0;
		for (int i = 0; i < k; i++) {
			double r = 0.0;
			for (int j = 0; j < k; j++) r += A(idx[i], idx[j]) * x(idx[j]);
			s += x(idx[i]) * r;
		}
		return s;
	}

	//A += s x x^T on the nonzero block idx
	void rank1_block(MatrixXd& A, const ivec& idx, const double s, const VectorXd& x)
	{
		const int k = static_cast<int>(idx.size());
		for (int i = 0; i < k; i++) {
			const double xi = s * x(idx[i]);
			for (int j = 0; j < k; j++) A(idx[i], idx[j]) += xi * x(idx[j]);
		}
	}

	//Resonance structure: bond multiplicities off the diagonal, lone pairs on it.  Cores are excluded,
	//as in NBO's own output: they are in every structure and distinguish none
	struct Topology {
		int n = 0;
		ivec t;

		Topology() = default;
		explicit Topology(const int na) : n(na), t(static_cast<size_t>(na) * na, 0) {}
		int& at(const int a, const int b) { return t[static_cast<size_t>(a) * n + b]; }
		int at(const int a, const int b) const { return t[static_cast<size_t>(a) * n + b]; }
		void set(const int a, const int b, const int v) { at(a, b) = v; at(b, a) = v; }

		//orbitals on the atom: lone pairs plus bond multiplicities
		int used(const int a) const
		{
			int u = at(a, a);
			for (int b = 0; b < n; b++) if (b != a) u += at(a, b);
			return u;
		}
		//Lewis electron count, both electrons of a shared bond counted
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
		//upper triangle with diagonal: the deduplication key
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

	//Relabelling-invariant key: sorted per-atom signatures (element, lone pairs, sorted neighbour/multiplicity)
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

	//NBO's "Added(Removed)" column
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

	//Pair moves allowed only where an E2 interaction exceeds nrt_e2_kcal
	vec2 delocalisation_graph(const NboLewis& lewis, const std::vector<NboE2Entry>& e2,
							  const double kcal, const int na, const bvec2& bondable)
	{
		vec2 g(na, vec(na, 0.0));
		for (const NboE2Entry& e : e2) {
			if (e.energy_kcal < kcal) continue;
			ivec c = lewis.orbitals[e.donor_index - 1].centers;
			const ivec& ac = lewis.orbitals[e.acceptor_index - 1].centers;
			c.insert(c.end(), ac.begin(), ac.end());
			for (size_t i = 0; i < c.size(); i++)
				for (size_t j = i + 1; j < c.size(); j++)
					if (c[i] != c[j] && bondable[c[i]][c[j]]) {
						const double p = std::max(g[c[i]][c[j]], e.energy_kcal);
						g[c[i]][c[j]] = g[c[j]][c[i]] = p;
					}
		}
		return g;
	}

	//Moves in different components are independent resonance whose weights factorise, so they are not combined
	ivec components(const vec2& g)
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
					if (g[x][y] > 0.0 && comp[y] < 0) { comp[y] = c; stack.push_back(y); }
			}
			c++;
		}
		return comp;
	}

	struct Candidate {
		Topology topo;
		int depth = 0;
		int component = -1;       //delocalisation-graph component the moves touched
		double score = 0.0;       //summed E2 price of its arrows, kcal/mol
		bool feasible = false;
		double g = 0.0;           //s Tr(V^T Gamma V)
		double rho_nl = 0.0;      //electrons outside the structure's own orbital set
		MatrixXd V;
		std::vector<std::array<double, 3>> polarity;  //{a, b, c_a^2 - c_b^2} per two-centre orbital
		int sweeps = 0;
		double change = 0.0;      //largest vector change of the last sweep
		bool converged = false;   //valence span stopped moving before the sweep cap
	};

	struct Limits {
		ivec cap;            //orbitals an atom may carry
		bvec free_atom;      //subspace screen
		bvec2 bondable;      //geometry screen
		vec2 deloc;          //E2 screen and kcal/mol price of each open pair
		ivec comp;
		int max_charge = 2;  //max change of an atom's Lewis electron count
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

	//One arrow moves one electron pair between slots; the charge-neutral NRT arrow is two of these,
	//hence nrt_max_arrows = 2
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
			d.score = c.score + L.deloc[x][y];
			out.push_back(std::move(d));
		};
		for (int a = 0; a < n; a++) {
			//a lone pair becomes a bond, or moves to another atom
			if (c.topo.at(a, a) > 0) {
				for (int b = 0; b < n; b++) {
					if (b == a || L.deloc[a][b] <= 0.0) continue;
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
				if (b == a || c.topo.at(a, b) <= 0 || L.deloc[a][b] <= 0.0) continue;
				if (L.ion) {
					Topology s = c.topo;
					s.set(a, b, s.at(a, b) - 1);
					s.at(a, a)++;
					emit(std::move(s), a, b);
				}
				for (int d = 0; d < n; d++) {
					if (d == a || d == b || L.deloc[a][d] <= 0.0) continue;
					Topology s = c.topo;
					s.set(a, b, s.at(a, b) - 1);
					s.set(a, d, s.at(a, d) + 1);
					emit(std::move(s), a, d);
				}
			}
		}
	}

	//Every admissible topology, not only those an arrow walk reaches; exponential in the bondable
	//pairs, capped by nrt_max_candidates
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
		//immovable slots keep their parent value outside the recursion
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
			//bond multiplicity above 3 is not a Lewis structure; lone pairs run to the atom's cap (F-, Ne: 4)
			const int mmax = (a == b) ? L.cap[a] : 3;
			for (int m = 0; m <= std::min(mmax, left); m++) {
				s.set(a, b, m);
				if (s.used(a) <= L.cap[a] && s.used(b) <= L.cap[b]) rec(i + 1, left - m);
				if (out.size() >= cap) break;
			}
			s.set(a, b, 0);
		};
		rec(0, n_pairs - fixed_pairs);
	}

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
		ivec core;                  //core NAOs, one fixed orbital each in every candidate
		std::map<std::pair<int, int>, ivec> pair;
		vec score_atom;             //leading eigenvalue of the undepleted block: the slot order
		std::map<std::pair<int, int>, double> score_pair;
	};

	//Orbitals of a prescribed topology by the NBO search's steps: greedy fill, self-consistency sweep,
	//occupancy-weighted orthogonalisation.  The sweep is required: a greedy pass leaves each orbital
	//biased by later ones and the residual stops distinguishing topologies
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
		if (k == 0) return;
		std::vector<VectorXd> v(k);
		std::vector<const ivec*> blk(k, nullptr);
		std::vector<std::pair<int, int>> owner(k, { -1, -1 });
		vec occ(k, 0.0);

		//cores: one NAO each, never swept
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

		//greedy fill from the density depleted by cores and earlier orbitals
		MatrixXd R = G0 - sum;
		for (const Slot& sl : slots) {
			const ivec& idx = (sl.a == sl.b) ? B.atom[sl.a] : B.pair.at({ sl.a, sl.b });
			if (static_cast<int>(idx.size()) < sl.mult) return;
			for (int m = 0; m < sl.mult; m++) {
				VectorXd x;
				const double lam = leading_block(R, idx, x);
				rank1_block(R, idx, -lam, x);
				v[col] = x;
				blk[col] = &idx;
				owner[col] = { sl.a, sl.b };
				occ[col] = block_quad(G0, idx, x);
				rank1_block(sum, idx, occ[col], x);
				col++;
			}
		}

		//Self consistency: each orbital against the density with all others removed, warm-started.
		//Converged on the valence span, which is all G, g and rho_nl depend on; the vectors keep wobbling
		//inside it, so a vector criterion never fires.  Anderson mixing starts only near convergence:
		//from the first sweep it can land on a different fixed point
		const int nv = k - first_valence;
		const auto orth = [&]() {
			MatrixXd Mv(n, nv);
			for (int j = first_valence; j < k; j++) Mv.col(j - first_valence) = v[j];
			Eigen::HouseholderQR<MatrixXd> qr(Mv);
			return MatrixXd(qr.householderQ() * MatrixXd::Identity(n, nv));
		};
		ivec off(nv + 1, 0);
		for (int j = 0; j < nv; j++) off[j + 1] = off[j] + static_cast<int>(blk[first_valence + j]->size());
		const auto pack = [&]() {
			VectorXd x(off[nv]);
			for (int j = 0; j < nv; j++) {
				const ivec& idx = *blk[first_valence + j];
				for (size_t t = 0; t < idx.size(); t++) x(off[j] + t) = v[first_valence + j](idx[t]);
			}
			return x;
		};
		constexpr int anderson_depth = 5;
		constexpr double anderson_start = 1e-3, span_start = 1e-6, span_tol = 1e-10;
		MatrixXd Q;  //span the sweep started from; a QR costs about a sweep, so only kept near the end
		std::deque<VectorXd> dX, dF;
		VectorXd x_prev, f_prev;
		c.converged = nv == 0;
		for (int s = 0; s < max_sweeps && !c.converged; s++) {
			const VectorXd x0 = pack();
			double change = 0.0;
			for (int j = first_valence; j < k; j++) {
				rank1_block(sum, *blk[j], -occ[j], v[j]);
				VectorXd x;
				leading_block(G0, *blk[j], x, &sum, &v[j]);
				if (x.dot(v[j]) < 0.0) x = -x;
				change = std::max(change, (x - v[j]).norm());
				v[j] = x;
				occ[j] = block_quad(G0, *blk[j], x);
				rank1_block(sum, *blk[j], occ[j], x);
			}
			c.sweeps = s + 1;
			c.change = change;
			if (change < span_start) {
				const MatrixXd Qn = orth();
				if (Q.size() && (Qn - Q * (Q.transpose() * Qn)).norm() < span_tol) {
					c.converged = true;
					break;
				}
				Q = Qn;
			}
			else Q.resize(0, 0);
			if (x_prev.size() == 0 && change > anderson_start) continue;
			//Anderson (type II) on the stacked block vectors
			const VectorXd F = pack(), f = F - x0;
			if (x_prev.size()) {
				dX.push_back(x0 - x_prev);
				dF.push_back(f - f_prev);
				if (static_cast<int>(dX.size()) > anderson_depth) { dX.pop_front(); dF.pop_front(); }
			}
			x_prev = x0;
			f_prev = f;
			if (dX.empty()) continue;
			const int h = static_cast<int>(dX.size());
			MatrixXd DX(F.size(), h), DF(F.size(), h);
			for (int i = 0; i < h; i++) { DX.col(i) = dX[i]; DF.col(i) = dF[i]; }
			const VectorXd xn = F - (DX + DF) * DF.colPivHouseholderQr().solve(f);
			sum.setZero();
			for (int j = 0; j < first_valence; j++) sum(B.core[j], B.core[j]) += occ[j];
			for (int j = 0; j < nv; j++) {
				const int o = first_valence + j;
				const ivec& idx = *blk[o];
				v[o].setZero();
				for (size_t t = 0; t < idx.size(); t++) v[o](idx[t]) = xn(off[j] + t);
				v[o].normalize();
				occ[o] = block_quad(G0, idx, v[o]);
				rank1_block(sum, idx, occ[o], v[o]);
			}
			if (Q.size()) Q = orth();  //the next sweep starts from the mixed iterate
		}

		//Occupancy-weighted symmetric orthogonalisation: V = M W (W S W)^(-1/2)
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

	//G v over the nonzero weights only: simplex iterates are sparse
	void gather_mv(const MatrixXd& G, const VectorXd& v, VectorXd& out)
	{
		const int n = static_cast<int>(v.size());
		out.setZero(n);
		for (int j = 0; j < n; j++) {
			const double vj = v(j);
			if (vj == 0.0) continue;
			out.noalias() += vj * G.col(j);
		}
	}

	//Equality-constrained KKT solve on a support, dropping the most negative weight until none is
	double solve_qp_support(const MatrixXd& G, const VectorXd& g, const double trg2, VectorXd& w)
	{
		const int n = static_cast<int>(G.rows());
		ivec act(n);
		for (int i = 0; i < n; i++) act[i] = i;
		for (int sweep = 0; sweep < n && !act.empty(); sweep++) {
			const int m = static_cast<int>(act.size());
			MatrixXd K = MatrixXd::Zero(m + 1, m + 1);
			VectorXd rhs(m + 1);
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < m; j++) K(i, j) = G(act[i], act[j]);
				K(i, m) = 1.0;
				K(m, i) = 1.0;
				rhs(i) = g(act[i]);
			}
			rhs(m) = 1.0;
			const VectorXd sol = K.completeOrthogonalDecomposition().solve(rhs);
			if (!sol.allFinite()) break;
			int worst = -1;
			double least = 0.0;
			for (int i = 0; i < m; i++)
				if (sol(i) < least) { least = sol(i); worst = i; }
			if (worst < 0) {
				w.setZero();
				for (int i = 0; i < m; i++) w(act[i]) = sol(i);
				return objective(G, g, trg2, w);
			}
			act.erase(act.begin() + worst);
		}
		return std::numeric_limits<double>::infinity();
	}

	//FISTA on the simplex with adaptive restart, step 1/L, L = 2 lambda_max(G)
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
		//Windowed progress test: a single restart step can stall without convergence
		double f_window = f;
		//G y follows by linearity from G wn and G wp, one product per iteration
		VectorXd Gy(n), Gwn(n), Gwp(n);
		gather_mv(G, y, Gy);
		Gwp = Gy;                                     //wp == y == w on entry
		for (int it = 1; it <= maxit; it++) {
			VectorXd wn = y - (2.0 / L) * (Gy - g);
			project_simplex(wn);
			gather_mv(G, wn, Gwn);
			const double fn = trg2 - 2.0 * g.dot(wn) + wn.dot(Gwn);
			if (fn > f) { y = wn; t = 1.0; Gy = Gwn; }  //restart: the momentum overshot
			else {
				const double tn = 0.5 * (1.0 + std::sqrt(1.0 + 4.0 * t * t));
				const double c = (t - 1.0) / tn;
				y = wn + c * (wn - wp);
				Gy = Gwn + c * (Gwn - Gwp);
				t = tn;
			}
			const double step = (wn - wp).cwiseAbs().maxCoeff();
			const double df = std::abs(fn - f);
			wp = wn;
			Gwp = Gwn;
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
			if (it % 500 == 0) {
				if (f_window - f < 1e-11 * std::max(1.0, std::abs(f))) break;
				f_window = f;
			}
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
	L.bondable = bondable;  //nrt_bond_scale screen, not the search's
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
		//Cap = valence NAOs, but never below the parent's count: a hypervalent parent must stay reachable
		for (int a = 0; a < na; a++) L.cap[a] = std::max(nval[a], parent.used(a));
	}

	//Budget from octet slots and a Gram-matrix work guard (k^2 n_NAO per pair); -nrt_max overrides both
	int budget = options.nrt_max_candidates;
	if (!options.nrt_max_set) {
		int slots = 0, active = 0;
		for (int a = 0; a < na; a++) {
			bool deloc = false;
			for (int b = 0; b < na && !deloc; b++) deloc = L.deloc[a][b] > 0.0;
			if (!deloc) continue;
			active++;
			slots += 1 + std::max(0, L.cap[a] - parent.used(a));
		}
		int ncore_total = 0;
		for (const int c : ncore) ncore_total += c;
		const double k = n_pairs + ncore_total;
		const int afford = static_cast<int>(std::sqrt(1.0e14 / std::max(k * k * nn, 1.0)));
		const int chem = std::max(64, NRT_PER_SLOT * slots);
		const int want = std::min(chem, std::max(64, afford));
		if (want < budget) {
			budget = want;
			//The budget changes the answer, so it is reported
			log << "NRT" << (spin.empty() ? "" : " " + spin) << ": candidate budget " << budget
				<< " of " << options.nrt_max_candidates << ": " << active << " delocalising atom(s), "
				<< slots << " octet slot(s) -> " << chem << ", machine guard " << afford << " at "
				<< nn << " NAOs (-nrt_max overrides both)\n";
			nrt.notes.push_back("candidate budget " + std::to_string(budget) + " of " +
								 std::to_string(options.nrt_max_candidates) + ": " +
								 std::to_string(active) + " delocalising atom(s), " +
								 std::to_string(slots) + " octet slot(s) -> " + std::to_string(chem) +
								 ", machine guard " + std::to_string(afford) + " at " +
								 std::to_string(nn) + " NAOs (-nrt_max overrides both)");
		}
	}

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
		enumerate(parent, L, n_pairs, static_cast<size_t>(budget), all);
		for (Candidate& c : all)
			if (seen.emplace(c.topo.key(), static_cast<int>(cands.size())).second)
				cands.push_back(std::move(c));
		nrt.notes.push_back("exhaustive enumeration yields " + std::to_string(cands.size()) +
							 " feasible topologies");
	}
	else {
		//Rank intermediates by summed E2 price before applying the budget.
		const auto richer = [](const Candidate& a, const Candidate& b) { return a.score > b.score; };
		const auto shortlist = [&richer](std::vector<Candidate>& v, const size_t keep) {
			std::stable_sort(v.begin(), v.end(), richer);
			if (v.size() > keep) v.resize(keep);
		};
		size_t level_begin = 0, paid = 0;
		for (int depth = 1; depth <= options.nrt_max_arrows; depth++) {
			const size_t level_end = cands.size();
			const size_t level_cap = (depth == 1) ? 4 * static_cast<size_t>(budget)
												  : static_cast<size_t>(budget) - paid;
			std::vector<Candidate> made;
			for (size_t i = level_begin; i < level_end; i++) {
				expand(cands[i], L, made);
				//O(n^2) children per parent: prune before the level grows an order of magnitude past its cap
				if (made.size() > 8 * level_cap + 1024) shortlist(made, level_cap);
			}
			shortlist(made, level_cap);
			size_t added = 0;
			for (Candidate& c : made) {
				if (seen.emplace(c.topo.key(), static_cast<int>(cands.size())).second) {
					if (c.depth > 1) paid++;
					cands.push_back(std::move(c));
					added++;
				}
			}
			nrt.arrows.push_back("ARROWS depth " + std::to_string(depth) + " generates " +
								 std::to_string(added) + " new structures from " +
								 std::to_string(level_end - level_begin));
			level_begin = level_end;
			if (!added || paid >= static_cast<size_t>(budget)) break;
		}
		//A resonance arrow couples two pair moves; depth-one states are intermediates.
		const size_t before = cands.size();
		cands.erase(std::remove_if(cands.begin() + 1, cands.end(),
								   [](const Candidate& c) { return c.depth == 1; }),
					cands.end());
		nrt.notes.push_back("half-arrow intermediates dropped: " +
							 std::to_string(before - cands.size()) + " of " +
							 std::to_string(before));
	}
	const double search_seconds = secs(t_search0, clock());

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
	//Candidate-independent slot order: leading eigenvalue of the slot block with cores removed, so a
	//core is never taken for a lone pair
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
		build_orbitals(cands[i], gamma, B, 200);
		if (!cands[i].feasible) continue;
		const MatrixXd GV = gamma * cands[i].V;
		double tr = 0.0;
		for (int c = 0; c < cands[i].V.cols(); c++) tr += cands[i].V.col(c).dot(GV.col(c));
		cands[i].g = scale * tr;
		cands[i].rho_nl = electrons - tr;
	}

	if (options.debug) {
		int capped = 0, feas = 0;
		long total = 0;
		double worst = 0.0;
		for (const Candidate& c : cands) {
			if (!c.feasible) continue;
			feas++;
			total += c.sweeps;
			if (!c.converged) { capped++; worst = std::max(worst, c.change); }
		}
		log << "NRT" << (spin.empty() ? "" : " " + spin) << ": orbital sweeps " << total << " over "
			<< feas << " candidates, " << capped << " stopped at the cap (worst change " << worst
			<< ", parent " << cands[0].change << " after " << cands[0].sweeps << "); last change by decade:";
		ivec dec(13, 0);
		for (const Candidate& c : cands)
			if (c.feasible)
				dec[std::min(12, std::max(0, static_cast<int>(-std::floor(std::log10(std::max(c.change, 1e-300))))))]++;
		for (int d = 0; d < 13; d++)
			if (dec[d]) log << " 1e-" << d << ":" << dec[d];
		log << "\n";
	}
	//parent stays at index 0
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
	const double orbital_seconds = secs(t_gram0, clock());
	const auto t_pairs0 = clock();

	//G_ij = s^2 ||V_i^T V_j||_F^2 = s^2 <P_i, P_j>: pair loop or projector product, whichever is cheaper
	int kmax = 0;
	for (const Candidate& c : cands) kmax = std::max(kmax, static_cast<int>(c.V.cols()));
	const double dnc = static_cast<double>(nc), dnn = static_cast<double>(nn);
	const double cost_pairs = 0.5 * dnc * dnc * static_cast<double>(kmax) * kmax * dnn;
	const double cost_proj = 0.25 * dnc * dnc * dnn * dnn + 0.5 * dnc * dnn * dnn * kmax;
	const bool projector = cost_proj < cost_pairs;
	if (projector) {
		//One (R, C) NAO block per pass, sized so Z (nc columns of wr*wc rows) stays near 128 MB
		int nbs = static_cast<int>(std::sqrt(134217728.0 / (8.0 * std::max(1, nc))));
		nbs = std::max(16, std::min(nbs, nn));
		const int nblk = (nn + nbs - 1) / nbs;
		std::vector<std::pair<int, int>> blocks;
		for (int R = 0; R < nblk; R++)
			for (int C = R; C < nblk; C++) blocks.push_back({ R, C });
		vec Zbuf(static_cast<size_t>(nbs) * nbs * nc);
		MatrixXd Goff = MatrixXd::Zero(nc, nc);      //blocks with R < C, counted twice below
		G.setZero();                                 //blocks with R == C, counted once
		//Upper-triangle candidate tiles; diagonal tiles also fill their lower half to stay a plain GEMM,
		//since Eigen does not parallelise the triangular kernel
		const int tile = 128;
		std::vector<std::pair<int, int>> tiles;
		for (int I = 0; I < nc; I += tile)
			for (int J = I; J < nc; J += tile) tiles.push_back({ I, J });
		const int ntiles = static_cast<int>(tiles.size());
		for (const std::pair<int, int>& blk : blocks) {
			const int r0 = blk.first * nbs, c0 = blk.second * nbs;
			const int wr = std::min(nbs, nn - r0), wc = std::min(nbs, nn - c0);
			const size_t rows = static_cast<size_t>(wr) * wc;
#pragma omp parallel for schedule(static) num_threads(nthreads)
			for (int a = 0; a < nc; a++) {
				Eigen::Map<MatrixXd> Y(Zbuf.data() + static_cast<size_t>(a) * rows, wr, wc);
				Y.noalias() = cands[a].V.middleRows(r0, wr) *
							  cands[a].V.middleRows(c0, wc).transpose();
			}
			const Eigen::Map<const MatrixXd> Z(Zbuf.data(), static_cast<Eigen::Index>(rows), nc);
			MatrixXd& acc = (blk.first == blk.second) ? G : Goff;
#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
			for (int t = 0; t < ntiles; t++) {
				const int I = tiles[t].first, J = tiles[t].second;
				const int wi = std::min(tile, nc - I), wj = std::min(tile, nc - J);
				acc.block(I, J, wi, wj).noalias() +=
					Z.middleCols(I, wi).transpose() * Z.middleCols(J, wj);
			}
		}
		const double s2 = scale * scale;
		for (int i = 0; i < nc; i++)
			for (int j = i; j < nc; j++) {
				G(i, j) = s2 * (G(i, j) + 2.0 * Goff(i, j));
				G(j, i) = G(i, j);
			}

		//Probe the tiled product against the direct formula
		const int probes = 96;
		//nc * nc passes INT_MAX from about 46000 candidates
		const size_t ncc = static_cast<size_t>(nc) * nc;
		const size_t stride = std::max<size_t>(1, ncc / probes);
		double worst = 0.0, scale_ref = 1e-300;
		int bad_i = -1, bad_j = -1;
		for (size_t p = 0; p < ncc; p += stride) {
			const int i = static_cast<int>(p / nc), j = static_cast<int>(p % nc);
			if (j < i) continue;
			const double direct = s2 * (cands[i].V.transpose() * cands[j].V).squaredNorm();
			const double d = std::abs(direct - G(i, j));
			scale_ref = std::max(scale_ref, std::abs(direct));
			if (d > worst) { worst = d; bad_i = i; bad_j = j; }
		}
		err_checkf(worst <= 1e-8 * scale_ref,
				   "NRT: the projector route disagrees with the pair product at G(" +
					   std::to_string(bad_i) + "," + std::to_string(bad_j) + ") by " +
					   std::to_string(worst),
				   log);
	}
	else {
#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
		for (int i = 0; i < nc; i++)
			for (int j = i; j < nc; j++) {
				const double v = scale * scale *
								 (cands[i].V.transpose() * cands[j].V).squaredNorm();
				G(i, j) = v;
				G(j, i) = v;
			}
	}
	const double pair_seconds = secs(t_pairs0, clock());
	const double trg2 = gamma.squaredNorm();
	const double gram_seconds = secs(t_gram0, clock());

	const auto t_min0 = clock();
	VectorXd w = VectorXd::Zero(nc);
	w(0) = 1.0;
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
	//Re-solve on the support above the reporting floor, so the reported weights solve the reported
	//problem; the argmin is not unique, only the residual is
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
		VectorXd wq = ws;
		double fs = solve_qp_support(Gs, gs, trg2, wq);
		if (fs <= f + 1e-12) ws = wq;
		else fs = solve_qp(Gs, gs, trg2, ws, 200000, &rs, spin, nullptr);   //singular: iterate instead
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

	int orbits = 0, in_orbits = 0;
	if (options.nrt_symmetry) {
		//Orbits: graph-isomorphic structures with equal g; one weight per orbit, shared equally
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
			//A rise means the orbits are no symmetry of this density; accepted within 1 % of D(w)
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

	const auto t_other0 = clock();
	vec2 bo(na, vec(na, 0.0));        //bond orders, diagonal = lone pairs
	vec2 pol(na, vec(na, 0.0));       //weight-summed c_A^2 - c_B^2 of the bonds on the pair
	//Topology units: electron pairs closed shell, single electrons per spin open shell
	const double unit = scale / 2.0;
	for (int i = 0; i < nc; i++) {
		if (w(i) <= 0.0) continue;
		for (int a = 0; a < na; a++)
			for (int b = a; b < na; b++)
				bo[a][b] += unit * w(i) * cands[i].topo.at(a, b);
		for (const std::array<double, 3>& p : cands[i].polarity) {
			const int a = static_cast<int>(p[0]), b = static_cast<int>(p[1]);
			pol[std::min(a, b)][std::max(a, b)] += unit * w(i) * ((a < b) ? p[2] : -p[2]);
		}
	}

	//Ionic share of a bond order is |i| b, i = c_A^2 - c_B^2: linear in |i|, not i^2 (NBO 7 convention)
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
		//octet count: both electrons of a shared bond counted for each partner, cores excluded
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
	//every structure tied for the leading weight
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
			<< " (search " << search_seconds << " s, gram " << gram_seconds << " s [orbitals "
			<< orbital_seconds << " s, pairs " << pair_seconds << " s, "
			<< (projector ? "projector" : "pair loop") << "], minimise "
			<< minimize_seconds << " s)\n";
}
