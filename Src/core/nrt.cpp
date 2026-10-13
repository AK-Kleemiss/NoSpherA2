#include "pch.h"
#include "nbo.h"
#include "constants.h"
#include "wfn_class.h"

#include <Eigen/Dense>
#include <chrono>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>

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

	//Leading eigenpair on block idx by warm-started power iteration; where that stalls or meets a negative dominant
	//eigenvalue, inverse iteration from the warm vector, and the full symmetric solve without one
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

		if (warm && k > 2) {
			//lambda_max from the tridiagonal form, then sigma - B is positive definite for sigma just above it
			Eigen::Tridiagonalization<MatrixXd> tri(B);
			Eigen::SelfAdjointEigenSolver<MatrixXd> ev;
			ev.computeFromTridiagonal(tri.diagonal(), tri.subDiagonal(), Eigen::EigenvaluesOnly);
			const double lam = ev.eigenvalues()(k - 1);
			MatrixXd A = -B;
			A.diagonal().array() += lam + 1e-9 * std::max(1.0, std::abs(lam));
			const Eigen::LLT<MatrixXd> llt(A);
			VectorXd x(k);
			for (int i = 0; i < k; i++) x(i) = (*warm)(idx[i]);
			if (llt.info() == Eigen::Success && x.norm() > 0.0) {
				x.normalize();
				for (int it = 0; it < 20; it++) {
					VectorXd y = llt.solve(x);
					y.normalize();
					if (y.dot(x) < 0.0) y = -y;
					const double d = (y - x).norm();
					x = y;
					if (d < 1e-13) { embed(x); return lam; }
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
		vec Mp;                   //orbitals V = M diag(wt) X, kept as M packed over each column's block, rebuilt by orbitals()
		VectorXd wt;
		MatrixXd X;
		std::vector<std::array<double, 3>> polarity;  //{a, b, c_a^2 - c_b^2} per two-centre orbital
		int sweeps = 0;
		double change = 0.0;      //largest vector change of the last sweep
		bool converged = false;   //valence span stopped moving before the sweep cap
		std::array<int, 3> moved{ -1, -1, -1 };  //atoms the last pair move touched
		int root = -1;                             //the weighted structure whose walk produced it
		MatrixXd M;                                //vectors given per slot: a start, after a build those before the orthogonalisation
		std::vector<std::pair<int, int>> owner;    //slot of each column of M, and after a build of V; {-1, -1} for a core
		bool warm = false;                         //started from a symmetry image's rotated vectors
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
		const auto emit = [&](Topology s, const int x, const int y, const int z = -1) {
			if (!L.free_atom[x] || !L.free_atom[y]) return;
			if (!acceptable(s, L)) return;
			//Two moves that share no atom are two separate charge separations, not one arrow
			const auto touched = [&](const int t) {
				return t >= 0 && std::find(c.moved.begin(), c.moved.end(), t) != c.moved.end();
			};
			if (c.depth > 0 && !touched(x) && !touched(y) && !touched(z)) return;
			Candidate d;
			d.topo = std::move(s);
			d.root = c.root;
			d.depth = c.depth + 1;
			d.moved = { x, y, z };
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
					emit(std::move(s), a, d, b);
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

	//Weisfeiler-Lehman colours of the bondable skeleton, starting from Z, three rounds
	ivec skeleton_colours(const ivec& Z, const bvec2& adj)
	{
		const int n = static_cast<int>(Z.size());
		ivec col = Z;
		for (int r = 0; r < 3; r++) {
			std::map<std::pair<int, ivec>, int> id;
			ivec next(n);
			for (int a = 0; a < n; a++) {
				ivec nb;
				for (int b = 0; b < n; b++)
					if (b != a && adj[a][b]) nb.push_back(col[b]);
				std::sort(nb.begin(), nb.end());
				next[a] = id.emplace(std::make_pair(col[a], nb), static_cast<int>(id.size())).first->second;
			}
			col = next;
		}
		return col;
	}

	//An automorphism pi of the coloured skeleton with pi[u] = v, atoms placed breadth first from u; with distances
	//(bohr) an isometry to 0.01 bohr
	bool skeleton_map(const int u, const int v, const ivec& col, const bvec2& adj, const vec2* dist, ivec& pi, long budget)
	{
		const int n = static_cast<int>(col.size());
		ivec order, from(n, -1);
		bvec in(n, false);
		for (int r = -1; r < n; r++) {
			const int s = r < 0 ? u : r;
			if (in[s]) continue;
			in[s] = true;
			order.push_back(s);
			for (size_t k = order.size() - 1; k < order.size(); k++)
				for (int b = 0; b < n; b++)
					if (!in[b] && adj[order[k]][b]) { in[b] = true; from[b] = order[k]; order.push_back(b); }
		}
		pi.assign(n, -1);
		bvec used(n, false);
		std::function<bool(int)> place = [&](const int k) {
			if (k == n) return true;
			if (--budget < 0) return false;
			const int a = order[k];
			for (int b = 0; b < n; b++) {
				if (used[b] || col[b] != col[a] || (k == 0 && b != v) || (from[a] >= 0 && !adj[pi[from[a]]][b])) continue;
				bool ok = true;
				for (int j = 0; j < k && ok; j++)
					ok = dist ? std::abs((*dist)[a][order[j]] - (*dist)[b][pi[order[j]]]) < 0.01 : adj[a][order[j]] == adj[b][pi[order[j]]];
				if (!ok) continue;
				pi[a] = b;
				used[b] = true;
				if (place(k + 1)) return true;
				pi[a] = -1;
				used[b] = false;
			}
			return false;
		};
		return place(0);
	}

	//Isometries of the coloured skeleton, one per pair of atoms the maps found so far do not already join into one orbit
	std::vector<ivec> symmetry_maps(const ivec& Z, const Limits& L, const vec2* dist)
	{
		const int n = static_cast<int>(Z.size());
		const ivec col = skeleton_colours(Z, L.bondable);
		std::vector<ivec> maps;
		ivec pi, orbit(n);
		for (int a = 0; a < n; a++) orbit[a] = a;
		const auto root = [&](int a) {
			while (orbit[a] != a) a = orbit[a] = orbit[orbit[a]];
			return a;
		};
		for (int u = 0; u < n; u++)
			for (int v = u + 1; v < n; v++)
				if (col[v] == col[u] && root(u) != root(v) && skeleton_map(u, v, col, L.bondable, dist, pi, 20000)) {
					maps.push_back(pi);
					for (int a = 0; a < n; a++) orbit[root(a)] = root(pi[a]);
				}
		return maps;
	}

	//The images of s under the maps, closed under composition up to cap images, without those that overfill an atom or
	//move a fixed atom. The charge screen is the parent's and not symmetric: an image passes it for the parent's image
	std::vector<Topology> images_of(const Topology& s, const std::vector<ivec>& maps, const Limits& L, const size_t cap, bool& capped)
	{
		const int n = s.n;
		std::set<ivec> seen{ s.key() };
		std::vector<Topology> all{ s };
		capped = false;
		for (size_t i = 0; i < all.size() && !capped; i++)
			for (const ivec& p : maps) {
				Topology t(n);
				for (int a = 0; a < n; a++)
					for (int b = 0; b < n; b++) t.at(p[a], p[b]) = all[i].at(a, b);
				if (!seen.insert(t.key()).second) continue;
				if (all.size() > cap) { capped = true; break; }
				all.push_back(std::move(t));
			}
		std::vector<Topology> out;
		for (size_t i = 1; i < all.size(); i++) {
			bool keep = true;
			for (int a = 0; a < n && keep; a++) keep = all[i].used(a) <= L.cap[a];
			for (int a = 0; a < n && keep; a++)
				for (int b = 0; b < n && keep && !L.free_atom[a]; b++) keep = all[i].at(a, b) == s.at(a, b);
			if (keep) out.push_back(all[i]);
		}
		return out;
	}

	//Real spherical harmonics of order l at the unit vector u, in FILE47's m = 0, +1, -1, +2, -2, ... order, without
	//the Condon-Shortley phase: P_l^m by the upward recursion in l
	VectorXd real_harmonics(const int l, const Eigen::Vector3d& u)
	{
		const double x = std::max(-1.0, std::min(1.0, u(2))), st = std::sqrt(std::max(0.0, 1.0 - x * x));
		const double phi = std::atan2(u(1), u(0));
		VectorXd y(2 * l + 1);
		double ratio = 1.0;  //(l - m)! / (l + m)!
		for (int m = 0; m <= l; m++) {
			if (m) ratio /= (l - m + 1.0) * (l + m);
			double pmm = 1.0;
			for (int i = 1; i <= m; i++) pmm *= (2 * i - 1) * st;
			double p0 = pmm, p1 = l > m ? x * (2 * m + 1) * pmm : pmm;
			for (int ll = m + 2; ll <= l; ll++) {
				const double p2 = ((2 * ll - 1) * x * p1 - (ll + m - 1) * p0) / (ll - m);
				p0 = p1;
				p1 = p2;
			}
			const double norm = std::sqrt((2 * l + 1) / (4.0 * constants::PI) * ratio) * p1;
			if (m == 0) y(0) = norm;
			else {
				y(2 * m - 1) = std::sqrt(2.0) * norm * std::cos(m * phi);
				y(2 * m) = std::sqrt(2.0) * norm * std::sin(m * phi);
			}
		}
		return y;
	}

	//c' = D c carries sum_m c_m Y_m(r) to sum_m c_m Y_m(R^T r); fitted on Fibonacci directions
	MatrixXd wigner_real(const int l, const Eigen::Matrix3d& R)
	{
		const int k = 4 * (2 * l + 1);
		MatrixXd A(k, 2 * l + 1), Y(k, 2 * l + 1);
		for (int i = 0; i < k; i++) {
			const double z = 1.0 - (2.0 * i + 1.0) / k, r = std::sqrt(1.0 - z * z), phi = 2.399963229728653 * i;
			const Eigen::Vector3d u(r * std::cos(phi), r * std::sin(phi), z);
			A.row(i) = real_harmonics(l, u);
			Y.row(i) = real_harmonics(l, R.transpose() * u);
		}
		return A.colPivHouseholderQr().solve(Y);
	}

	//NAOs of one (atom, l, shell), idx[c] the NAO of component c
	struct NaoGroup {
		int atom = 0, l = 0, shell = 0;
		ivec idx;
	};

	//An isometry in the NAO basis: group g goes to group target[g] through D[g], which carries the group's sign;
	//miss = ||U Gamma U^T - Gamma|| / ||Gamma||
	struct SymOp {
		ivec perm, target;
		std::vector<MatrixXd> D;
		double miss = 0.0;
	};

	MatrixXd apply_op(const SymOp& op, const std::vector<NaoGroup>& groups, const MatrixXd& V, const bool inverse)
	{
		MatrixXd out = MatrixXd::Zero(V.rows(), V.cols());
		for (size_t g = 0; g < groups.size(); g++) {
			const ivec& from = inverse ? groups[op.target[g]].idx : groups[g].idx;
			const ivec& to = inverse ? groups[g].idx : groups[op.target[g]].idx;
			const int w = static_cast<int>(from.size());
			MatrixXd x(w, V.cols());
			for (int c = 0; c < w; c++) x.row(c) = V.row(from[c]);
			const MatrixXd y = inverse ? MatrixXd(op.D[g].transpose() * x) : MatrixXd(op.D[g] * x);
			for (int c = 0; c < w; c++) out.row(to[c]) = y.row(c);
		}
		return out;
	}

	//The NAO-basis image of the atom map pi: R from the coordinates (improper allowed), D^l(R) per group, and per group
	//the sign the NAO construction left arbitrary, taken from the density along a maximum spanning tree of the group
	//couplings. Without and with ORCA's sign on |m| = 3, 4, whichever U Gamma U^T matches better; empty when that
	//misses Gamma by 1 %
	std::optional<SymOp> nao_op(const ivec& pi, const vec2& xyz, const std::vector<NaoGroup>& groups, const MatrixXd& gamma)
	{
		const int na = static_cast<int>(pi.size()), ng = static_cast<int>(groups.size());
		Eigen::Vector3d c0 = Eigen::Vector3d::Zero();
		for (int a = 0; a < na; a++) c0 += Eigen::Vector3d(xyz[a][0], xyz[a][1], xyz[a][2]) / na;
		Eigen::Matrix3d H = Eigen::Matrix3d::Zero();
		for (int a = 0; a < na; a++)
			H += (Eigen::Vector3d(xyz[a][0], xyz[a][1], xyz[a][2]) - c0) *
				 (Eigen::Vector3d(xyz[pi[a]][0], xyz[pi[a]][1], xyz[pi[a]][2]) - c0).transpose();
		const Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
		const Eigen::Matrix3d R = svd.matrixV() * svd.matrixU().transpose();
		for (int a = 0; a < na; a++)
			if ((R * (Eigen::Vector3d(xyz[a][0], xyz[a][1], xyz[a][2]) - c0) -
				 (Eigen::Vector3d(xyz[pi[a]][0], xyz[pi[a]][1], xyz[pi[a]][2]) - c0)).norm() > 0.01) return std::nullopt;
		SymOp op;
		op.perm = pi;
		op.target.assign(ng, -1);
		for (int g = 0; g < ng; g++)
			for (int h = 0; h < ng && op.target[g] < 0; h++)
				if (groups[h].atom == pi[groups[g].atom] && groups[h].l == groups[g].l && groups[h].shell == groups[g].shell)
					op.target[g] = h;
		for (int g = 0; g < ng; g++)
			if (op.target[g] < 0 || groups[op.target[g]].idx.size() != groups[g].idx.size()) return std::nullopt;
		const auto block = [&](const int g, const int h) {
			MatrixXd b(groups[g].idx.size(), groups[h].idx.size());
			for (int i = 0; i < b.rows(); i++)
				for (int j = 0; j < b.cols(); j++) b(i, j) = gamma(groups[g].idx[i], groups[h].idx[j]);
			return b;
		};
		std::map<int, MatrixXd> wig;
		std::optional<SymOp> best;
		for (int orca = 0; orca < 2; orca++) {
			op.D.assign(ng, MatrixXd());
			for (int g = 0; g < ng; g++) {
				const int l = groups[g].l;
				if (!wig.count(l)) wig[l] = wigner_real(l, R);
				op.D[g] = wig[l];
				if (orca)
					for (int c = 1; c < 2 * l + 1; c++)
						if (orca_pure_sign_flips((c + 1) / 2)) { op.D[g].row(c) *= -1.0; op.D[g].col(c) *= -1.0; }
			}
			vec sign(ng, 0.0), link(ng, 0.0);
			ivec from(ng, -1);
			for (int left = ng; left > 0; left--) {
				int g = -1;
				for (int h = 0; h < ng; h++)
					if (sign[h] == 0.0 && (g < 0 || link[h] > link[g])) g = h;
				sign[g] = 1.0;
				if (from[g] >= 0) {
					const int f = from[g];
					sign[g] = sign[f] * ((block(op.target[f], op.target[g]).cwiseProduct(op.D[f] * block(f, g) * op.D[g].transpose())).sum() < 0.0 ? -1.0 : 1.0);
				}
				for (int h = 0; h < ng; h++)
					if (sign[h] == 0.0) {
						const double c = block(g, h).norm();
						if (c > link[h]) { link[h] = c; from[h] = g; }
					}
			}
			for (int g = 0; g < ng; g++) op.D[g] *= sign[g];
			double num = 0.0, den = 0.0;
			for (int g = 0; g < ng; g++)
				for (int h = 0; h < ng; h++) {
					num += (block(op.target[g], op.target[h]) - op.D[g] * block(g, h) * op.D[h].transpose()).squaredNorm();
					den += block(g, h).squaredNorm();
				}
			op.miss = std::sqrt(num / std::max(den, 1e-300));
			if (!best || op.miss < best->miss) best = op;
		}
		if (best && best->miss < 1e-2) return best;
		return std::nullopt;
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

	//V of a built candidate, the product build_orbitals formed
	MatrixXd orbitals(const Candidate& c, const Blocks& B, const int n)
	{
		const int k = static_cast<int>(c.owner.size());
		MatrixXd M = MatrixXd::Zero(n, k);
		size_t p = 0;
		for (int j = 0; j < k; j++) {
			const std::pair<int, int>& o = c.owner[j];
			if (o.first < 0) M(B.core[j], j) = c.Mp[p++];
			else
				for (const int i : o.first == o.second ? B.atom[o.first] : B.pair.at(o)) M(i, j) = c.Mp[p++];
		}
		return M * c.wt.asDiagonal() * c.X;
	}

	//Orbitals of a prescribed topology by the NBO search's steps: greedy fill, self-consistency sweep,
	//occupancy-weighted orthogonalisation.  The sweep is required: a greedy pass leaves each orbital
	//biased by later ones and the residual stops distinguishing topologies.  Returns V, empty if infeasible
	MatrixXd build_orbitals(Candidate& c, const MatrixXd& G0, const Blocks& B, const int max_sweeps, const bool keep = false,
							const int stall = 0)
	{
		struct Slot { int a, b, mult; double score; };
		std::vector<Slot> slots;
		for (int a = 0; a < c.topo.n; a++)
			if (c.topo.at(a, a) > 0) slots.push_back({ a, a, c.topo.at(a, a), B.score_atom[a] });
		for (int a = 0; a < c.topo.n; a++)
			for (int b = a + 1; b < c.topo.n; b++)
				if (c.topo.at(a, b) > 0) {
					const auto it = B.score_pair.find({ a, b });
					if (it == B.score_pair.end()) return {};   //not a bondable pair: infeasible
					slots.push_back({ a, b, c.topo.at(a, b), it->second });
				}
		//Lone pairs first, as in the search: a bond slot filled earlier takes a metal d pair and keeps
		//it at every sweep, leaving the ligand pair it should hold unclaimed
		std::sort(slots.begin(), slots.end(), [](const Slot& x, const Slot& y) {
			if ((x.a == x.b) != (y.a == y.b)) return x.a == x.b;
			if (x.score != y.score) return x.score > y.score;
			return std::tie(x.a, x.b) < std::tie(y.a, y.b);
		});

		const int n = static_cast<int>(G0.rows());
		int k = static_cast<int>(B.core.size());
		for (const Slot& sl : slots) k += sl.mult;
		if (k > n) return {};
		if (k == 0) return {};
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

		//a slot keeps the vectors M gives it (a symmetry image's, rotated, or the walk root's); the others are filled
		//greedily in slot order from the density depleted by cores and every kept vector
		const int given = static_cast<int>(c.M.cols());
		bvec taken(given, false);
		ivec open;
		for (const Slot& sl : slots) {
			const ivec& idx = (sl.a == sl.b) ? B.atom[sl.a] : B.pair.at({ sl.a, sl.b });
			if (static_cast<int>(idx.size()) < sl.mult) return {};
			for (int m = 0; m < sl.mult; m++) {
				int t = 0;
				while (t < given && (taken[t] || c.owner[t] != std::make_pair(sl.a, sl.b))) t++;
				blk[col] = &idx;
				owner[col] = { sl.a, sl.b };
				if (t < given) {
					taken[t] = true;
					v[col] = c.M.col(t).normalized();
					occ[col] = block_quad(G0, idx, v[col]);
					rank1_block(sum, idx, occ[col], v[col]);
				}
				else open.push_back(col);
				col++;
			}
		}
		c.warm = given > 0 && open.empty();
		MatrixXd R = G0 - sum;
		for (const int j : open) {
			VectorXd x;
			const double lam = leading_block(R, *blk[j], x);
			rank1_block(R, *blk[j], -lam, x);
			v[j] = x;
			occ[j] = block_quad(G0, *blk[j], x);
			rank1_block(sum, *blk[j], occ[j], x);
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
		double least = 0.0;   //-nrt_stall: smallest change since mixing started, and its sweep
		int least_at = -1;
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
			if (stall > 0 && (x_prev.size() || change <= anderson_start) && (least_at < 0 || change < 0.999 * least)) {
				least = change;
				least_at = s;
			}
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
			if (stall > 0 && least_at >= 0 && s - least_at >= stall) break;
		}

		//Occupancy-weighted symmetric orthogonalisation: V = M W (W S W)^(-1/2)
		MatrixXd M(n, k);
		for (int j = 0; j < k; j++) M.col(j) = v[j];
		VectorXd wt(k);
		for (int j = 0; j < k; j++) wt(j) = std::max(occ[j], 1e-6);
		c.owner = owner;
		if (keep) c.M = M;
		else c.M.resize(0, 0);
		const MatrixXd S = M.transpose() * M;
		if (Eigen::SelfAdjointEigenSolver<MatrixXd>(S).eigenvalues()(0) < 1e-8) return {};
		c.wt = wt;
		c.X = sym_power(MatrixXd(wt.asDiagonal() * S * wt.asDiagonal()), -0.5);
		c.Mp.clear();
		for (int j = 0; j < k; j++)
			if (j < first_valence) c.Mp.push_back(M(B.core[j], j));
			else
				for (const int i : *blk[j]) c.Mp.push_back(M(i, j));
		MatrixXd V = orbitals(c, B, n);

		for (int j = first_valence; j < k; j++) {
			if (owner[j].first == owner[j].second) continue;
			const int a = owner[j].first, b = owner[j].second;
			double wa = 0.0, wb = 0.0;
			for (const int i : B.atom[a]) wa += V(i, j) * V(i, j);
			for (const int i : B.atom[b]) wb += V(i, j) * V(i, j);
			c.polarity.push_back({ static_cast<double>(a), static_cast<double>(b), wa - wb });
		}
		c.feasible = true;
		return V;
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

	//Frank-Wolfe gap grad.w - min_i grad_i, grad = 2 (G w - g): an upper bound on f(w) - min f over the
	//simplex, zero exactly at the optimum.  Gw = G w.
	double fw_gap(const VectorXd& Gw, const VectorXd& g, const VectorXd& w)
	{
		const VectorXd grad = 2.0 * (Gw - g);
		return std::max(0.0, grad.dot(w) - grad.minCoeff());
	}

	//FISTA on the simplex with adaptive restart, step 1/L, L = 2 lambda_max(G).  status: a progress line
	//every 500 iterations.  gap_out: the Frank-Wolfe gap of the returned weights.
	double solve_qp(const MatrixXd& G, const VectorXd& g, const double trg2, VectorXd& w,
					const int maxit, const vec* rho, const std::string& spin,
					std::vector<NboQpIteration>* trace, std::ostream* status = nullptr,
					double* gap_out = nullptr, int* it_out = nullptr, const std::function<double(double)>& nbo = {})
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
		int it = 1;
		for (; it <= maxit; it++) {
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
			if (it % 50 == 0 && fw_gap(Gwn, g, wn) <= 1e-10 * std::max(1e-6, f)) break;
			if (it % 500 == 0) {
				if (status)
					*status << "NRT" << (spin.empty() ? "" : " " + spin) << ": iteration " << it << ", D(w) "
							<< std::sqrt(std::max(0.0, f))
							<< (nbo ? " (NBO " + std::to_string(nbo(std::sqrt(std::max(0.0, f)))) + ")" : "")
							<< ", " << (wn.array() > 1e-6).count()
							<< " weighted, gap " << fw_gap(Gwn, g, wn) << std::endl;
				if (f_window - f < 1e-11 * std::max(1.0, std::abs(f))) break;
				f_window = f;
			}
		}
		w = wp;
		if (gap_out) gap_out[0] = fw_gap(Gwp, g, wp);
		if (it_out) it_out[0] = std::min(it, maxit);
		return f;
	}
}

void native_nrt(NboNrt& nrt, const NAOResult& nao, const NboLewis& lewis,
				const std::vector<NboE2Entry>& e2, const bvec2& bondable,
				const NboOptions& options, const std::string& spin, const double scale,
				std::ostream& log, const vec2* xyz)
{
	const auto clock = [] { return std::chrono::steady_clock::now(); };
	const auto secs = [](const std::chrono::steady_clock::time_point a,
						 const std::chrono::steady_clock::time_point b) {
		return std::chrono::duration<double>(b - a).count();
	};
	const auto t_start = clock();

	const int na = static_cast<int>(nao.atoms.size());
	const int nn = static_cast<int>(nao.orbitals.size());
	//NBO 7's D per spin is sqrt((||Gamma_s - sum_a w_a P_a||^2 + 0.001) / n_NAO), Gamma_s = Gamma / scale: one structure
	//leaves 2 rho_NL(s) in the norm, and n D^2 - 2 rho_NL(s) = 0.0010 in every spin of the 41 nbo_reference_v2 goldens
	const auto nbo_units = [&](const double d) { return std::sqrt((d * d / (scale * scale) + 0.001) / nn); };
	const auto nbo_d = [&](const double d) { return " (NBO " + std::to_string(nbo_units(d)) + ")"; };
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
	nrt.parent_topo.push_back(NboTopo{ spin, parent.matrix() });
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
	//-nrt_target growth stops at this many candidates: the Gram matrix is nc^2
	size_t pool_cap = static_cast<size_t>(options.nrt_max_candidates);
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
		pool_cap = std::min(pool_cap, static_cast<size_t>(std::max(64, afford)));
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
	std::set<ivec> seen;
	std::set<ivec> walked{ parent.key() };  //roots the arrow walk has started from
	{
		Candidate p;
		p.topo = parent;
		cands.push_back(p);
		seen.insert(parent.key());
	}
	//Arrow walk from roots, at most `want` structures past depth one.  A resonance arrow couples two
	//pair moves, so depth-one states are intermediates: walked through, counted, not returned
	size_t intermediates = 0;
	const auto grow = [&](std::vector<Candidate> level, const size_t want, const bool report) {
		//Rank intermediates by summed E2 price before applying the budget.
		const auto richer = [](const Candidate& a, const Candidate& b) { return a.score > b.score; };
		const auto shortlist = [&richer](std::vector<Candidate>& v, const size_t keep) {
			std::stable_sort(v.begin(), v.end(), richer);
			if (v.size() > keep) v.resize(keep);
		};
		std::vector<Candidate> found;
		size_t paid = 0;
		for (int depth = 1; depth <= options.nrt_max_arrows; depth++) {
			const size_t level_cap = (depth == 1) ? 4 * want : want - paid;
			std::vector<Candidate> made;
			for (const Candidate& c : level) {
				expand(c, L, made);
				//O(n^2) children per parent: prune before the level grows an order of magnitude past its cap
				if (made.size() > 8 * level_cap + 1024) shortlist(made, level_cap);
			}
			shortlist(made, level_cap);
			std::vector<Candidate> next;
			for (Candidate& c : made)
				if (seen.insert(c.topo.key()).second) {
					if (c.depth > 1) { paid++; found.push_back(c); }
					else intermediates++;
					next.push_back(std::move(c));
				}
			if (report)
				nrt.arrows.push_back("ARROWS depth " + std::to_string(depth) + " generates " +
									 std::to_string(next.size()) + " new structures from " +
									 std::to_string(level.size()));
			level = std::move(next);
			if (level.empty() || paid >= want) break;
		}
		return found;
	};
	const auto t_search0 = clock();
	vec2 dist;
	if (xyz) {
		dist.assign(na, vec(na, 0.0));
		for (int a = 0; a < na; a++)
			for (int b = 0; b < na; b++)
				dist[a][b] = std::hypot((*xyz)[a][0] - (*xyz)[b][0], (*xyz)[a][1] - (*xyz)[b][1], (*xyz)[a][2] - (*xyz)[b][2]);
	}
	const std::vector<ivec> maps = options.nrt_symmetry && !options.nrt_exhaustive ? symmetry_maps(Z, L, xyz ? &dist : nullptr) : std::vector<ivec>();
	//the walk's budget can cut through a symmetry orbit; the images of what it kept complete it
	const auto complete_orbits = [&](const size_t from) {
		bool capped = false;
		for (size_t i = from, end = cands.size(); i < end; i++) {
			const Candidate c0 = cands[i];
			for (Topology& t : images_of(c0.topo, maps, L, 64, capped))
				if (seen.insert(t.key()).second) {
					Candidate c;
					c.topo = std::move(t);
					c.depth = c0.depth;
					c.score = c0.score;
					c.root = c0.root;
					cands.push_back(std::move(c));
				}
		}
	};
	if (options.nrt_exhaustive) {
		std::vector<Candidate> all;
		enumerate(parent, L, n_pairs, static_cast<size_t>(budget), all);
		for (Candidate& c : all)
			if (seen.insert(c.topo.key()).second)
				cands.push_back(std::move(c));
		nrt.notes.push_back("exhaustive enumeration yields " + std::to_string(cands.size()) +
							 " feasible topologies");
	}
	else {
		for (Candidate& c : grow({ cands[0] }, static_cast<size_t>(budget), true))
			cands.push_back(std::move(c));
		nrt.notes.push_back("half-arrow intermediates dropped: " + std::to_string(intermediates) +
							 " of " + std::to_string(cands.size() + intermediates));
		//NBO's symmetry references: the structures the parent's own invariant matches are walked as well, so
		//a parent that breaks the molecule's symmetry (Ni(CO)4's one Ni-C bond) does not tilt the arrow walk
		//toward its own ligand
		if (options.nrt_symmetry) {
			const std::string pform = canonical_form(parent, Z);
			std::vector<Candidate> refs;
			for (size_t i = 1; i < cands.size(); i++)
				if (canonical_form(cands[i].topo, Z) == pform) {
					Candidate r;
					r.topo = cands[i].topo;
					refs.push_back(std::move(r));
				}
			//the images no arrow walk reaches, such as benzene's second Kekule structure three pair moves away
			bool capped = false;
			int images = 0;
			for (Topology& t : images_of(parent, maps, L, 64, capped))
				if (seen.insert(t.key()).second) {
					Candidate c;
					c.topo = t;
					cands.push_back(std::move(c));
					Candidate r;
					r.topo = std::move(t);
					refs.push_back(std::move(r));
					images++;
				}
			nrt.notes.push_back(std::to_string(images) + " symmetry image(s) of the parent added" +
								 (capped ? ", the closure stopped at 64 images" : ""));
			if (!refs.empty() && cands.size() < pool_cap) {
				const size_t before = cands.size(), n_refs = refs.size();
				for (const Candidate& r : refs) walked.insert(r.topo.key());
				const size_t room = std::min(n_refs * static_cast<size_t>(budget), pool_cap - cands.size());
				for (Candidate& c : grow(std::move(refs), room, false)) cands.push_back(std::move(c));
				nrt.notes.push_back(std::to_string(n_refs) + " symmetry reference(s) of the parent walked, " +
									 std::to_string(cands.size() - before) + " new structures");
			}
			const size_t before = cands.size();
			complete_orbits(1);
			nrt.notes.push_back(std::to_string(cands.size() - before) + " symmetry image(s) complete the walk's orbits");
		}
	}
	double search_seconds = secs(t_search0, clock());

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

	//the isometries in the NAO basis: an image of a built structure starts its sweeps from the source's vectors, rotated
	std::vector<NaoGroup> groups;
	std::vector<SymOp> ops;
	if (options.nrt_rotate && xyz && !maps.empty()) {
		std::map<std::array<int, 3>, std::vector<std::pair<int, int>>> members;
		for (int i = 0; i < nn; i++) members[{ nao.orbitals[i].atom, nao.orbitals[i].l, nao.orbitals[i].shell }].push_back({ nao.orbitals[i].m, i });
		for (auto& [k, v] : members) {
			std::sort(v.begin(), v.end());
			NaoGroup grp{ k[0], k[1], k[2], ivec() };
			for (const auto& mi : v) grp.idx.push_back(mi.second);
			groups.push_back(std::move(grp));
		}
		for (const ivec& pi : maps)
			if (std::optional<SymOp> op = nao_op(pi, *xyz, groups, gamma)) ops.push_back(std::move(*op));
		if (ops.size() < maps.size())
			nrt.notes.push_back(std::to_string(maps.size() - ops.size()) + " of " + std::to_string(maps.size()) +
								 " isometries miss the NAO density by 1 % or more; their images start from scratch");
	}
	int rotated = 0, stalled = 0;
#ifdef _OPENMP
	const int nthreads = options.threads > 0 ? options.threads : omp_get_max_threads();
#else
	const int nthreads = 1;
#endif
	const double trg2 = gamma.squaredNorm();
	double orbital_seconds = 0.0, pair_seconds = 0.0, gram_seconds = 0.0, minimize_seconds = 0.0;
	bool projector = false;
	int nc = 0;
	MatrixXd G;
	VectorXd g, w;
	vec rho;
	double f = 0.0, f0 = 0.0, fista_gap = 0.0;
	int fista_its = 0;
	//-nrt_target: walk arrows from every structure the fit weights and refit, until D(w) reaches the
	//target or no unwalked structure carries weight (NBO 7's multi-reference search)
	const size_t trace_mark = nrt.qp_iterations.size();
	size_t built = 0;
	auto t_cycle0 = t_gram0;
	//-nrt_screen: a walk's structures get orbitals from their root's untouched slots, a greedy fill and three sweeps;
	//from those, the reduced cost r_j = 2 (s Tr(V_j^T Gamma_w V_j) - g_j) - 2 (w.G.w - g.w), Gamma_w = sum_i w_i s V_i V_i^T,
	//admits the 16 most negative below the threshold with their symmetry orbits, built in full, refitted, until none is left
	std::vector<Candidate> deferred;
	int screen_rounds = 0;
	size_t screen_walked = 0, screen_built = 0;
	const auto approximate = [&](Candidate& c) {
		if (c.root >= 0) {
			const Candidate& r = cands[c.root];
			bvec touched(na, false);
			for (int a = 0; a < na; a++)
				for (int b = 0; b < na; b++)
					if (c.topo.at(a, b) != r.topo.at(a, b)) touched[a] = touched[b] = true;
			const MatrixXd rV = orbitals(r, B, nn);
			ivec keep;
			for (int q = 0; q < rV.cols(); q++)
				if (r.owner[q].first >= 0 && !touched[r.owner[q].first] && !touched[r.owner[q].second]) keep.push_back(q);
			c.M = MatrixXd::Zero(nn, keep.size());
			for (size_t q = 0; q < keep.size(); q++) {
				const std::pair<int, int> o = r.owner[keep[q]];
				for (const int i : o.first == o.second ? B.atom[o.first] : B.pair.at(o)) c.M(i, q) = rV(i, keep[q]);
				c.owner.push_back(o);
			}
		}
		const MatrixXd V = build_orbitals(c, gamma, B, 3);
		if (!c.feasible) return;
		const MatrixXd GV = gamma * V;
		double tr = 0.0;
		for (int y = 0; y < V.cols(); y++) tr += V.col(y).dot(GV.col(y));
		c.g = scale * tr;
	};
	const auto admit = [&]() {
		const int nd = static_cast<int>(deferred.size());
		if (!nd) return 0;
		MatrixXd gw = MatrixXd::Zero(nn, nn);
		for (int i = 0; i < nc; i++)
			if (w(i) > 0.0) {
				const MatrixXd V = orbitals(cands[i], B, nn);
				gw.noalias() += (w(i) * scale) * V * V.transpose();
			}
		VectorXd Gw;
		gather_mv(G, w, Gw);
		const double level = w.dot(Gw) - g.dot(w);
		vec r(nd);
#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
		for (int j = 0; j < nd; j++) {
			const MatrixXd V = orbitals(deferred[j], B, nn);
			const MatrixXd GV = gw * V;
			double tr = 0.0;
			for (int y = 0; y < V.cols(); y++) tr += V.col(y).dot(GV.col(y));
			r[j] = 2.0 * (scale * tr - deferred[j].g) - 2.0 * level;
		}
		ivec order(nd);
		for (int j = 0; j < nd; j++) order[j] = j;
		std::stable_sort(order.begin(), order.end(), [&](const int a, const int b) { return r[a] < r[b]; });
		std::set<ivec> pick;
		bool capped = false;
		for (int j = 0, seeds = 0; j < nd && seeds < 16 && r[order[j]] < options.nrt_screen; j++)
			if (pick.insert(deferred[order[j]].topo.key()).second) {
				seeds++;
				for (const Topology& t : images_of(deferred[order[j]].topo, maps, L, 64, capped)) pick.insert(t.key());
			}
		int joined = 0;
		std::vector<Candidate> rest;
		for (Candidate& c : deferred)
			if (pick.count(c.topo.key())) {
				Candidate d;
				d.topo = std::move(c.topo);
				d.depth = c.depth;
				d.component = c.component;
				d.score = c.score;
				d.root = c.root;
				cands.push_back(std::move(d));
				joined++;
			}
			else rest.push_back(std::move(c));
		deferred = std::move(rest);
		screen_built += joined;
		return joined;
	};
	for (int cycle = 1;;) {
		const int nc0 = static_cast<int>(cands.size());
		//progress output only where a step takes long enough to watch; -no_date silences it
		const bool show = !constants::hide_timings && nc0 >= 500;
		std::optional<ProgressBar> pb;
		//source[i] = {j, op, inverse}: breadth first through the ops from the structures that have orbitals, then from
		//each structure no op reaches, which is built
		std::vector<std::array<int, 3>> source(nc0, { -1, -1, 0 });
		ivec order;
		if (!ops.empty()) {
			std::map<ivec, int> index;
			for (int i = 0; i < nc0; i++) index.emplace(cands[i].topo.key(), i);
			bvec reached(nc0, false);
			std::deque<int> queue;
			const auto spread = [&]() {
				while (!queue.empty()) {
					const int j = queue.front();
					queue.pop_front();
					for (int k = 0; k < static_cast<int>(ops.size()); k++)
						for (int inv = 0; inv < 2; inv++) {
							const ivec& p = ops[k].perm;
							Topology t(na);
							for (int a = 0; a < na; a++)
								for (int b = 0; b < na; b++)
									if (inv) t.at(a, b) = cands[j].topo.at(p[a], p[b]);
									else t.at(p[a], p[b]) = cands[j].topo.at(a, b);
							const auto it = index.find(t.key());
							if (it == index.end() || reached[it->second]) continue;
							reached[it->second] = true;
							source[it->second] = { j, k, inv };
							order.push_back(it->second);
							queue.push_back(it->second);
						}
				}
			};
			for (int i = 0; i < static_cast<int>(built); i++) {
				reached[i] = true;
				queue.push_back(i);
			}
			spread();
			for (int i = static_cast<int>(built); i < nc0; i++)
				if (!reached[i]) {
					reached[i] = true;
					queue.push_back(i);
					spread();
				}
		}
		//levels of the breadth-first tree: a structure starts from its source's rotated vectors once the source is built,
		//and only from a source whose sweeps converged; one that stopped at the cap leaves its images their own start
		std::vector<ivec> levels(1);
		ivec depth(nc0, 0);
		bvec keep(nc0, false);
		for (int i = static_cast<int>(built); i < nc0; i++)
			if (source[i][0] < 0) levels[0].push_back(i);
		for (const int i : order) {
			keep[source[i][0]] = true;
			depth[i] = depth[source[i][0]] + 1;
			if (static_cast<int>(levels.size()) <= depth[i]) levels.emplace_back();
			levels[depth[i]].push_back(i);
		}
		for (const ivec& list : levels) {
			for (const int i : list) {
				if (source[i][0] < 0) continue;
				const Candidate& s = cands[source[i][0]];
				if (!s.feasible || !s.converged || !s.M.cols()) continue;
				const SymOp& op = ops[source[i][1]];
				const bool inv = source[i][2] != 0;
				ivec p = op.perm;
				if (inv)
					for (int a = 0; a < na; a++) p[op.perm[a]] = a;
				Candidate& c = cands[i];
				c.M = apply_op(op, groups, s.M, inv);
				c.owner = s.owner;
				for (std::pair<int, int>& o : c.owner)
					if (o.first >= 0) o = std::minmax(p[o.first], p[o.second]);
			}
			const int nb = static_cast<int>(list.size());
			if (show && nb >= 500)
				pb.emplace(static_cast<unsigned long long>(nb), 50, "=", " ", "NRT: orbitals of " + std::to_string(nb) + " structures", log);
#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
			for (int x = 0; x < nb; x++) {
				if (pb) pb->update();
				Candidate& c = cands[list[x]];
				const MatrixXd V = build_orbitals(c, gamma, B, 200, keep[list[x]], options.nrt_stall);
				if (!c.feasible) continue;
				const MatrixXd GV = gamma * V;
				double tr = 0.0;
				for (int y = 0; y < V.cols(); y++) tr += V.col(y).dot(GV.col(y));
				c.g = scale * tr;
				c.rho_nl = electrons - tr;
			}
			pb.reset();
			for (const int i : list)
				if (cands[i].feasible && !cands[i].converged && cands[i].sweeps < 200) stalled++;
		}
		for (int i = static_cast<int>(built); i < nc0; i++) {
			rotated += cands[i].feasible && cands[i].warm;
			cands[i].M.resize(0, 0);
		}

		if (options.debug) {
			int capped = 0, feas = 0, warm = 0;
			long total = 0, warm_sweeps = 0;
			double worst = 0.0;
			for (const Candidate& c : cands) {
				if (!c.feasible) continue;
				feas++;
				total += c.sweeps;
				if (c.warm) { warm++; warm_sweeps += c.sweeps; }
				if (!c.converged) { capped++; worst = std::max(worst, c.change); }
			}
			log << "NRT" << (spin.empty() ? "" : " " + spin) << ": orbital sweeps " << total << " over "
				<< feas << " candidates (" << warm_sweeps << " over " << warm << " rotated starts), " << capped << " stopped at the cap (worst change " << worst
				<< ", parent " << cands[0].change << " after " << cands[0].sweeps << "); last change by decade:";
			ivec dec(13, 0);
			for (const Candidate& c : cands)
				if (c.feasible)
					dec[std::min(12, std::max(0, static_cast<int>(-std::floor(std::log10(std::max(c.change, 1e-300))))))]++;
			for (int d = 0; d < 13; d++)
				if (dec[d]) log << " 1e-" << d << ":" << dec[d];
			log << "\n";
		}
		//parent stays at index 0, earlier cycles keep their indices
		size_t kept = built;
		for (size_t i = built; i < cands.size(); i++)
			if (cands[i].feasible) {
				if (i != kept) cands[kept] = std::move(cands[i]);
				kept++;
			}
		cands.resize(kept);
		err_checkf(!cands.empty() && cands[0].topo.key() == parent.key(),
				   "NRT: the parent Lewis structure is not representable in its own NAO basis", log);
		nc = static_cast<int>(cands.size());

		//structures of earlier cycles keep their indices, so their block of G stands
		const int n0 = static_cast<int>(built);
		G.conservativeResize(nc, nc);
		g.resize(nc);
		rho.assign(nc, 0.0);
		for (int i = 0; i < nc; i++) { g(i) = cands[i].g; rho[i] = cands[i].rho_nl; }
		orbital_seconds += secs(t_cycle0, clock());
		const auto t_pairs0 = clock();

		//G_ij = s^2 ||V_i^T V_j||_F^2 = s^2 <P_i, P_j>: pair loop or projector product, whichever is cheaper
		int kmax = 0;
		for (const Candidate& c : cands) kmax = std::max(kmax, static_cast<int>(c.X.cols()));
		const double dnc = static_cast<double>(nc), dnn = static_cast<double>(nn), dnew = dnc * dnc - static_cast<double>(n0) * n0;
		const double cost_pairs = 0.5 * dnew * static_cast<double>(kmax) * kmax * dnn;
		const double cost_proj = 0.25 * dnew * dnn * dnn + 0.5 * dnc * dnn * dnn * kmax;
		projector = cost_proj < cost_pairs;
		if (projector) {
			//One (R, C) NAO block per pass, sized so Z (nc columns of wr*wc rows) stays near 128 MB
			int nbs = static_cast<int>(std::sqrt(134217728.0 / (8.0 * std::max(1, nc))));
			nbs = std::max(16, std::min(nbs, nn));
			const int nblk = (nn + nbs - 1) / nbs;
			std::vector<std::pair<int, int>> blocks;
			for (int R = 0; R < nblk; R++)
				for (int C = R; C < nblk; C++) blocks.push_back({ R, C });
			vec Zbuf(static_cast<size_t>(nbs) * nbs * nc);
			MatrixXd Goff = MatrixXd::Zero(nblk > 1 ? nc : 0, nc - n0);  //blocks with R < C over the new columns, counted twice below
			G.rightCols(nc - n0).setZero();              //blocks with R == C, counted once
			G.bottomRows(nc - n0).setZero();
			//Tiles over the new structures' columns and every row up to them; a tile across the diagonal also fills
			//part of the lower half to stay a plain GEMM, since Eigen does not parallelise the triangular kernel
			const int tile = 128;
			std::vector<std::pair<int, int>> tiles;
			for (int J = n0; J < nc; J += tile)
				for (int I = 0; I < std::min(nc, J + tile); I += tile) tiles.push_back({ I, J });
			const int ntiles = static_cast<int>(tiles.size());
			if (show)
				pb.emplace(static_cast<unsigned long long>(blocks.size()), 50, "=", " ",
						   "NRT: overlaps of " + std::to_string(nc) + " structures", log);
			for (const std::pair<int, int>& blk : blocks) {
				if (pb) pb->update();
				const int r0 = blk.first * nbs, c0 = blk.second * nbs;
				const int wr = std::min(nbs, nn - r0), wc = std::min(nbs, nn - c0);
				const size_t rows = static_cast<size_t>(wr) * wc;
#pragma omp parallel for schedule(static) num_threads(nthreads)
				for (int a = 0; a < nc; a++) {
					const MatrixXd V = orbitals(cands[a], B, nn);
					Eigen::Map<MatrixXd> Y(Zbuf.data() + static_cast<size_t>(a) * rows, wr, wc);
					Y.noalias() = V.middleRows(r0, wr) * V.middleRows(c0, wc).transpose();
				}
				const Eigen::Map<const MatrixXd> Z(Zbuf.data(), static_cast<Eigen::Index>(rows), nc);
				const bool diag = blk.first == blk.second;
				MatrixXd& acc = diag ? G : Goff;
				const int shift = diag ? 0 : n0;
#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
				for (int t = 0; t < ntiles; t++) {
					const int I = tiles[t].first, J = tiles[t].second;
					const int wi = std::min(tile, nc - I), wj = std::min(tile, nc - J);
					acc.block(I, J - shift, wi, wj).noalias() +=
						Z.middleCols(I, wi).transpose() * Z.middleCols(J, wj);
				}
			}
			pb.reset();
			const double s2 = scale * scale;
			for (int i = 0; i < nc; i++)
				for (int j = std::max(i, n0); j < nc; j++) {
					G(i, j) = s2 * (nblk > 1 ? G(i, j) + 2.0 * Goff(i, j - n0) : G(i, j));
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
				const double direct = s2 * (orbitals(cands[i], B, nn).transpose() * orbitals(cands[j], B, nn)).squaredNorm();
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
			if (show)
				pb.emplace(static_cast<unsigned long long>(nc), 50, "=", " ",
						   "NRT: overlaps of " + std::to_string(nc) + " structures", log);
			std::vector<MatrixXd> Vn(nc - n0);
#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
			for (int j = n0; j < nc; j++) Vn[j - n0] = orbitals(cands[j], B, nn);
#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
			for (int i = 0; i < nc; i++) {
				if (pb) pb->update();
				const MatrixXd Vi = i >= n0 ? Vn[i - n0] : orbitals(cands[i], B, nn);
				for (int j = std::max(i, n0); j < nc; j++) {
					const double v = scale * scale *
									 (Vi.transpose() * Vn[j - n0]).squaredNorm();
					G(i, j) = v;
					G(j, i) = v;
				}
			}
			pb.reset();
		}
		pair_seconds += secs(t_pairs0, clock());
		gram_seconds += secs(t_cycle0, clock());

		const auto t_min0 = clock();
		//warm start: the previous cycle's weights, zero on the new structures
		const VectorXd w_prev = w;
		w = VectorXd::Zero(nc);
		if (built == 0) {
			w(0) = 1.0;
			f0 = objective(G, g, trg2, w);
		}
		else w.head(w_prev.size()) = w_prev;
		//Frank-Wolfe gap of the warm start over the grown set: f(w_prev) - min f, so this cycle cannot take D(w)
		//below sqrt(f(w_prev) - gap)
		double floor_dw = -1.0;
		if (built > 0) {
			VectorXd Gw;
			gather_mv(G, w, Gw);
			floor_dw = std::sqrt(std::max(0.0, objective(G, g, trg2, w) - fw_gap(Gw, g, w)));
		}
		built = static_cast<size_t>(nc);
		nrt.qp_iterations.resize(trace_mark);  //the trace describes the final fit
		f = solve_qp(G, g, trg2, w, 5000, &rho, spin, &nrt.qp_iterations, show ? &log : nullptr,
					 &fista_gap, &fista_its, nbo_units);
		minimize_seconds += secs(t_min0, clock());

		const double dw = std::sqrt(std::max(0.0, f));
		if (options.nrt_exhaustive || options.nrt_target < 0.0 || dw <= options.nrt_target) {
			if (cycle > 1)
				nrt.notes.push_back("D(w) target " + std::to_string(options.nrt_target) + " reached after " +
									 std::to_string(cycle) + " cycles, " + std::to_string(nc) + " structures");
			break;
		}
		t_cycle0 = clock();
		if (admit()) {
			screen_rounds++;
			continue;
		}
		deferred.clear();
		std::vector<Candidate> roots;
		for (int i = 0; i < nc; i++)
			if (w(i) > options.nrt_weight_floor && walked.insert(cands[i].topo.key()).second) {
				Candidate r;
				r.topo = cands[i].topo;
				r.root = i;
				roots.push_back(std::move(r));
			}
		const size_t n_roots = roots.size();
		std::string stop;
		if (roots.empty()) stop = "no unwalked structure carries weight";
		else if (cands.size() >= pool_cap)
			stop = "the candidate pool is full at " + std::to_string(pool_cap) + " (-nrt_max raises it)";
		else {
			const size_t room = std::min(static_cast<size_t>(budget), pool_cap - cands.size());
			const size_t before = cands.size();
			for (Candidate& c : grow(std::move(roots), room, false)) cands.push_back(std::move(c));
			complete_orbits(before);
			if (cands.size() == before) stop = "the weighted structures generate no new structure";
			else if (options.nrt_screen >= 0.0 && cands.size() - before >= static_cast<size_t>(options.nrt_screen_min)) {
				const int nw = static_cast<int>(cands.size() - before);
#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
				for (int x = 0; x < nw; x++) approximate(cands[before + x]);
				for (size_t i = before; i < cands.size(); i++)
					if (cands[i].feasible) deferred.push_back(std::move(cands[i]));
				cands.resize(before);
				screen_walked += nw;
				if (!admit()) stop = "no structure of the walk passes the -nrt_screen threshold";
			}
		}
		nrt.notes.push_back("cycle " + std::to_string(cycle) + ": D(w) = " + std::to_string(dw) + nbo_d(dw) +
							 (floor_dw >= 0.0 ? ", bound before the fit " + std::to_string(floor_dw) + nbo_d(floor_dw) : std::string()) + " over " +
							 std::to_string(nc) + " structures, " + std::to_string(n_roots) +
							 " weighted structure(s) walked, " + std::to_string(cands.size() - static_cast<size_t>(nc)) + " new");
		if (!constants::hide_timings) log << "NRT" << (spin.empty() ? "" : " " + spin) << ": " << nrt.notes.back() << std::endl;
		search_seconds += secs(t_cycle0, clock());
		if (!stop.empty()) {
			nrt.notes.push_back("D(w) target " + std::to_string(options.nrt_target) + " not reached: " + stop);
			break;
		}
		t_cycle0 = clock();
		cycle++;
	}
	if (screen_walked)
		nrt.notes.push_back("-nrt_screen " + std::to_string(options.nrt_screen) + ": " + std::to_string(screen_built) + " of " +
							 std::to_string(screen_walked) + " walked structures built, over " + std::to_string(screen_rounds) + " refit(s)");
	if (rotated)
		nrt.notes.push_back(std::to_string(rotated) + " structure(s) started from a symmetry image's orbitals, rotated");
	if (stalled)
		nrt.notes.push_back("-nrt_stall " + std::to_string(options.nrt_stall) + ": " + std::to_string(stalled) +
							 " structure(s) stopped their sweeps unconverged before the cap of 200");
	const auto t_min1 = clock();
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
	//Active-set finish: FISTA stops on a stall or its cap, not on optimality.  Add the unweighted structures
	//whose gradient says they lower D(w), re-solve exactly on the grown support, until the Frank-Wolfe gap
	//over all nc structures certifies the optimum.
	{
		const double tol = 1e-10 * std::max(1e-6, f);
		double gap = fw_gap(G * w, g, w);
		const double gap0 = gap;
		int added = 0, pass = 0;
		for (; pass < 100 && gap > tol; pass++) {
			const VectorXd grad = 2.0 * (G * w - g);
			const double level = grad.dot(w);
			std::vector<std::pair<double, int>> viol;
			for (int i = 0; i < nc; i++)
				if (w(i) <= 0.0 && grad(i) < level) viol.push_back({ grad(i), i });
			if (viol.empty()) break;
			std::sort(viol.begin(), viol.end());
			ivec sup;
			for (int i = 0; i < nc; i++)
				if (w(i) > 0.0) sup.push_back(i);
			const int nadd = static_cast<int>(std::min<size_t>(10, viol.size()));
			for (int k = 0; k < nadd; k++) sup.push_back(viol[k].second);
			const int ns = static_cast<int>(sup.size());
			MatrixXd Gs(ns, ns);
			VectorXd gs(ns), ws(ns);
			for (int i = 0; i < ns; i++) {
				gs(i) = g(sup[i]);
				ws(i) = w(sup[i]);
				for (int j = 0; j < ns; j++) Gs(i, j) = G(sup[i], sup[j]);
			}
			VectorXd wq = ws;
			double fs = solve_qp_support(Gs, gs, trg2, wq);
			if (fs < f) ws = wq;
			else fs = solve_qp(Gs, gs, trg2, ws, 200000, nullptr, spin, nullptr);   //singular: iterate
			if (!(fs < f)) break;
			w.setZero();
			for (int i = 0; i < ns; i++) w(sup[i]) = ws(i);
			f = fs;
			added += nadd;
			gap = fw_gap(G * w, g, w);
		}
		std::ostringstream os;
		os << "minimiser: FISTA " << fista_its << " iterations, Frank-Wolfe gap " << fista_gap
		   << "; active-set finish " << pass << " pass(es), " << added << " structure(s) added, gap "
		   << gap0 << " -> " << gap << (gap > tol ? " (not certified)" : "");
		nrt.notes.push_back(os.str());
		if (!constants::hide_timings) log << "NRT" << (spin.empty() ? "" : " " + spin) << ": " << os.str() << std::endl;
	}

	int orbits = 0, in_orbits = 0;
	if (options.nrt_symmetry) {
		//Orbits: structures one symmetry map takes into another, or without maps graph-isomorphic structures with
		//equal g; one weight per orbit, shared equally
		std::map<std::string, ivec> classes;
		std::vector<ivec> member;
		bool capped = false;
		for (int i = 0; i < nc; i++) {
			if (!cands[i].feasible) continue;
			if (maps.empty()) { classes[canonical_form(cands[i].topo, Z)].push_back(i); continue; }
			ivec least = cands[i].topo.key();
			for (const Topology& t : images_of(cands[i].topo, maps, L, 64, capped)) least = std::min(least, t.key());
			std::string k(least.begin(), least.end());
			classes[k].push_back(i);
		}
		for (auto& kv : classes) {
			ivec& m = kv.second;
			if (!maps.empty()) { member.push_back(m); continue; }
			std::sort(m.begin(), m.end(),
					  [&](const int a, const int b) { return cands[a].g < cands[b].g; });
			for (size_t k = 0; k < m.size(); k++) {
				if (k && cands[m[k]].g - cands[m[k - 1]].g <= 1e-3)
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
	minimize_seconds += secs(t_min1, clock());

	const auto t_other0 = clock();
	vec2 bo(na, vec(na, 0.0));        //bond orders, diagonal = lone pairs
	vec2 pol(na, vec(na, 0.0));       //weight-summed |c_A^2 - c_B^2| of the bonds on the pair
	//Topology units: electron pairs closed shell, single electrons per spin open shell
	const double unit = scale / 2.0;
	for (int i = 0; i < nc; i++) {
		if (w(i) <= 0.0) continue;
		for (int a = 0; a < na; a++)
			for (int b = a; b < na; b++)
				bo[a][b] += unit * w(i) * cands[i].topo.at(a, b);
		for (const std::array<double, 3>& p : cands[i].polarity) {
			const int a = static_cast<int>(p[0]), b = static_cast<int>(p[1]);
			pol[std::min(a, b)][std::max(a, b)] += unit * w(i) * std::abs(p[2]);
		}
	}

	//Bond orders the near-optimal fits allow: f(w') >= f + sum_j w'_j r_j with reduced costs r_j >= 0, so within
	//D(w')^2 <= 1.001^2 f only structures with r_j <= eps carry weight; on them min and max of t.w' by bisection on the
	//multiplier mu of f(w') + mu t.w'
	const double eps = 0.002001 * std::max(f, 0.0);
	ivec nearw;
	{
		VectorXd Gw;
		gather_mv(G, w, Gw);
		const VectorXd grad = 2.0 * (Gw - g);
		const double level = grad.dot(w);
		for (int i = 0; i < nc; i++)
			if (w(i) > 0.0 || grad(i) - level <= eps) nearw.push_back(i);
	}
	const int nr = static_cast<int>(nearw.size());
	MatrixXd Gn(nr, nr);
	VectorXd gn(nr), wn(nr);
	for (int i = 0; i < nr; i++) {
		gn(i) = g(nearw[i]);
		wn(i) = w(nearw[i]);
		for (int j = 0; j < nr; j++) Gn(i, j) = G(nearw[i], nearw[j]);
	}
	std::map<std::pair<int, int>, std::pair<double, double>> ranges;
	const auto bond_range = [&](const int a, const int b) {
		//one bisection per pair orbit of the symmetry maps
		std::set<std::pair<int, int>> orbit{ { a, b } };
		std::vector<std::pair<int, int>> todo{ { a, b } };
		while (!todo.empty()) {
			const std::pair<int, int> q = todo.back();
			todo.pop_back();
			for (const ivec& p : maps) {
				const std::pair<int, int> r{ std::min(p[q.first], p[q.second]), std::max(p[q.first], p[q.second]) };
				if (orbit.insert(r).second) todo.push_back(r);
			}
		}
		const auto known = ranges.find(*orbit.begin());
		if (known != ranges.end()) return known->second;
		VectorXd t(nr);
		for (int i = 0; i < nr; i++) t(i) = unit * cands[nearw[i]].topo.at(a, b);
		double ext[2];
		for (int s = 0; s < 2; s++) {
			VectorXd x = wn, best = wn;
			const auto feasible = [&](const double mu) {
				solve_qp(Gn, gn - (s ? -0.5 : 0.5) * mu * t, trg2, x, 3000, nullptr, spin, nullptr);
				return objective(Gn, gn, trg2, x) <= f + eps;
			};
			double lo = 0.0, hi = 1e-3;
			while (hi < 1e6 && feasible(hi)) { best = x; lo = hi; hi *= 4.0; }
			for (int it = 0; it < 30; it++) {
				const double mid = 0.5 * (lo + hi);
				if (feasible(mid)) { best = x; lo = mid; }
				else hi = mid;
			}
			ext[s] = t.dot(best);
		}
		return ranges[*orbit.begin()] = { ext[0], ext[1] };
	};

	//Ionic share of a bond order is |i| b per structure, i = c_A^2 - c_B^2: linear in |i|, not i^2 (NBO 7 convention);
	//Ni(CO)4's donor and back-donor Ni-C orbitals have opposite i and must not cancel
	vec valency(na, 0.0), covalency(na, 0.0), electrovalency(na, 0.0);
	for (int a = 0; a < na; a++)
		for (int b = a + 1; b < na; b++) {
			const double t = bo[a][b];
			if (t <= 0.0) continue;
			const double ion = std::min(1.0, pol[a][b] / std::max(t, 1e-30));
			NboBondOrder o;
			o.atom1 = a + 1;
			o.atom2 = b + 1;
			o.spin = spin;
			o.total = t;
			o.ionic = ion * t;
			o.covalent = t - o.ionic;
			const std::pair<double, double> range = bond_range(a, b);
			o.range_lo = range.first;
			o.range_hi = range.second;
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
	const double d_w = std::sqrt(std::max(0.0, f)), d_0 = std::sqrt(std::max(0.0, f0));
	nrt.d_w = std::sqrt(nrt.d_w * nrt.d_w + d_w * d_w);
	nrt.d_0 = std::sqrt(nrt.d_0 * nrt.d_0 + d_0 * d_0);
	nrt.spin_d.push_back({ spin, nbo_units(d_0), nbo_units(d_w) });
	nrt.deloc_threshold_kcal = options.nrt_e2_kcal;
	nrt.max_search_cycles = static_cast<int>(nrt.cycles.size());
	nrt.initial_topo = nc;
	nrt.search_seconds += search_seconds;
	nrt.gram_seconds += gram_seconds;
	nrt.minimize_seconds += minimize_seconds;
	nrt.other_seconds += secs(t_other0, clock()) + secs(t_start, t_search0);
	if (options.debug)
		log << "NRT" << (spin.empty() ? "" : " " + spin) << ": " << nc << " candidates, "
			<< order.size() << " retained, D(0) = " << d_0 << ", D(w) = " << d_w
			<< " (search " << search_seconds << " s, gram " << gram_seconds << " s [orbitals "
			<< orbital_seconds << " s, pairs " << pair_seconds << " s, "
			<< (projector ? "projector" : "pair loop") << "], minimise "
			<< minimize_seconds << " s)\n";
}
