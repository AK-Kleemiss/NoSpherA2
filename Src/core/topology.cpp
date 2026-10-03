#include "pch.h"
#include "topology.h"
#include "constants.h"
#include "citations.h"
#include <iomanip>
#include <iostream>
#include <functional>
#include <queue>
#include <sstream>

namespace topology
{
	const char* kind_name(const cp_kind k)
	{
		switch (k) {
		case cp_kind::attractor: return "(3,-3) attractor";
		case cp_kind::bond: return "(3,-1) bond";
		case cp_kind::ring: return "(3,+1) ring";
		case cp_kind::cage: return "(3,+3) cage";
		default: return "degenerate";
		}
	}

	const char* seed_class_name(const seed_class s)
	{
		switch (s) {
		case seed_class::nucleus: return "nucleus";
		case seed_class::bond: return "bond";
		case seed_class::ring: return "ring";
		case seed_class::cage: return "cage";
		default: return "grid";
		}
	}

	//Same adjugate/determinant inverse b2c.cpp uses for its Newton step, and for the same reason:
	//a fixed 3x3 system solved thousands of times, where a LAPACK call is all overhead.  Kept here
	//rather than exported from b2c.cpp so that file needs no edit.
	bool invert_3x3_topology(const double m[9], double inv[9])
	{
		const double det =
			m[0] * (m[4] * m[8] - m[5] * m[7]) -
			m[1] * (m[3] * m[8] - m[5] * m[6]) +
			m[2] * (m[3] * m[7] - m[4] * m[6]);
		const double scale = std::max({ std::abs(m[0]), std::abs(m[1]), std::abs(m[2]),
										std::abs(m[3]), std::abs(m[4]), std::abs(m[5]),
										std::abs(m[6]), std::abs(m[7]), std::abs(m[8]) });
		const double threshold = std::max(1E-300, scale * scale * scale * 1E-12);
		if (!std::isfinite(det) || std::abs(det) < threshold)
			return false;
		const double id = 1.0 / det;
		inv[0] = (m[4] * m[8] - m[5] * m[7]) * id;
		inv[1] = (m[2] * m[7] - m[1] * m[8]) * id;
		inv[2] = (m[1] * m[5] - m[2] * m[4]) * id;
		inv[3] = (m[5] * m[6] - m[3] * m[8]) * id;
		inv[4] = (m[0] * m[8] - m[2] * m[6]) * id;
		inv[5] = (m[2] * m[3] - m[0] * m[5]) * id;
		inv[6] = (m[3] * m[7] - m[4] * m[6]) * id;
		inv[7] = (m[1] * m[6] - m[0] * m[7]) * id;
		inv[8] = (m[0] * m[4] - m[1] * m[3]) * id;
		return true;
	}

	d3 centroid(const std::vector<d3>& p)
	{
		d3 c{ 0.0, 0.0, 0.0 };
		if (p.empty())
			return c;
		for (const d3& q : p) for (int k = 0; k < 3; k++) c[k] += q[k];
		for (int k = 0; k < 3; k++) c[k] /= (double)p.size();
		return c;
	}

	cp_kind classify_hessian(const double H[9], d3& eigenvalues, int& negative, int& positive, int& zero, const double eigen_tolerance)
	{
		negative = positive = 0;
		zero = 3;
		eigenvalues = { 0.0, 0.0, 0.0 };
		vec A(H, H + 9), W(3);
		if (!try_make_Eigenvalues(A, W))
			return cp_kind::degenerate; //no eigenvalues, so no signature: not silently binned
		std::sort(W.begin(), W.end());
		eigenvalues = { W[0], W[1], W[2] };
		const double max_abs = std::max({ std::abs(W[0]), std::abs(W[1]), std::abs(W[2]) });
		//Relative to the largest curvature, with an absolute floor so a vacuum point where all
		//three eigenvalues are ~0 is degenerate rather than accidentally signed
		const double tol = std::max(1E-12, max_abs * eigen_tolerance);
		zero = 0;
		for (int i = 0; i < 3; i++) {
			if (W[i] < -tol) negative++;
			else if (W[i] > tol) positive++;
			else zero++;
		}
		if (zero > 0)
			return cp_kind::degenerate;
		if (negative == 3) return cp_kind::attractor;
		if (negative == 2) return cp_kind::bond;
		if (negative == 1) return cp_kind::ring;
		return cp_kind::cage;
	}

	std::vector<nucleus> nuclei_of(const WFN& wavy)
	{
		std::vector<nucleus> n;
		n.reserve(wavy.get_ncen());
		for (int a = 0; a < wavy.get_ncen(); a++)
			n.push_back({ wavy.get_atom_pos(a), wavy.get_atom_charge(a) });
		return n;
	}

	//A bond critical point lies on the path between the two nuclei it connects, so the two nearest
	//nuclei are that pair.  Cheap and standard; it only mis-assigns where two candidate pairs are
	//degenerate in distance, which would also make the ring seed indistinguishable.
	std::vector<std::vector<int>> bond_graph_cycles(result& r, const std::vector<nucleus>& nuclei)
	{
		const int V = (int)nuclei.size();
		r.graph_vertices = V;
		r.graph_edges = 0;
		r.graph_components = 0;
		r.required_ring_minus_cage = 0;
		if (V == 0)
			return {};

		std::vector<std::vector<int>> adj(V);
		std::vector<std::array<int, 2>> edges;
		for (const cp& point : r.points) {
			if (point.kind != cp_kind::bond)
				continue;
			int first = -1, second = -1;
			double d1 = std::numeric_limits<double>::infinity(), d2 = d1;
			for (int a = 0; a < V; a++) {
				const double d = array_length(point.position, nuclei[a].pos);
				if (d < d1) { d2 = d1; second = first; d1 = d; first = a; }
				else if (d < d2) { d2 = d; second = a; }
			}
			if (first < 0 || second < 0 || first == second)
				continue;
			const int u = std::min(first, second), v = std::max(first, second);
			bool known = false;
			for (const std::array<int, 2>& e : edges) if (e[0] == u && e[1] == v) { known = true; break; }
			if (known)
				continue;
			edges.push_back({ u, v });
			adj[u].push_back(v);
			adj[v].push_back(u);
		}
		r.graph_edges = (int)edges.size();

		//Spanning forest by BFS.  Its non-tree edges number E - V + C, the graph's cycle rank, which
		//is exactly how many independent rings the molecule has - and so, by Poincare-Hopf with
		//n_NCP = V and n_BCP = E, how many ring points minus cage points there must be
		std::vector<int> parent(V, -2);
		std::vector<std::array<int, 2>> non_tree;
		std::vector<std::pair<int, int>> tree_edges;
		for (int root = 0; root < V; root++) {
			if (parent[root] != -2)
				continue;
			r.graph_components++;
			parent[root] = -1;
			std::queue<int> q;
			q.push(root);
			while (!q.empty()) {
				const int u = q.front();
				q.pop();
				for (const int v : adj[u]) {
					if (parent[v] == -2) { parent[v] = u; tree_edges.push_back({ std::min(u, v), std::max(u, v) }); q.push(v); }
				}
			}
		}
		for (const std::array<int, 2>& e : edges) {
			bool in_tree = false;
			for (const std::pair<int, int>& t : tree_edges)
				if (t.first == e[0] && t.second == e[1]) { in_tree = true; break; }
			if (!in_tree)
				non_tree.push_back(e);
		}
		r.required_ring_minus_cage = r.graph_edges - r.graph_vertices + r.graph_components;

		std::vector<std::vector<int>> cycles;

		//Every simple cycle up to six vertices, enumerated exhaustively.  The cycle rank alone is
		//not enough to seed from: an octahedron has rank 12 - 6 + 1 = 7 but eight triangular faces,
		//because the eight are linearly dependent in the cycle space, so a basis of seven can only
		//ever reach seven of the eight ring points.  Enumerating the small cycles instead costs
		//nothing - a seed that is not inside a ring is rejected by the search, and two seeds inside
		//the same ring merge - and covers the 3- to 6-membered rings that are almost all of them.
		//ponytail: length capped at 6 and the list capped at max_cycles; the fundamental cycles
		//below still give one seed per independent ring of any size, so a macrocycle is not lost.
		{
			const size_t max_cycles = 5000;
			std::vector<int> path;
			std::vector<char> on_path(V, 0);
			//start only from the smallest vertex of the cycle, and fix the direction by requiring the
			//second vertex to be smaller than the last, so each cycle is enumerated exactly once
			std::function<void(const int, const int)> walk = [&](const int start, const int u) {
				if (cycles.size() >= max_cycles)
					return;
				for (const int v : adj[u]) {
					if (v == start) {
						if (path.size() >= 3 && path[1] < path.back())
							cycles.push_back(path);
						continue;
					}
					if (v < start || on_path[v] || path.size() >= 6)
						continue;
					on_path[v] = 1;
					path.push_back(v);
					walk(start, v);
					path.pop_back();
					on_path[v] = 0;
				}
			};
			for (int start = 0; start < V && cycles.size() < max_cycles; start++) {
				path.assign(1, start);
				std::fill(on_path.begin(), on_path.end(), (char)0);
				on_path[start] = 1;
				walk(start, start);
			}
		}

		//For each non-tree edge, the smallest cycle through it: the shortest path between its ends
		//in the graph with that edge cut.  This is what still gives a seed for a ring larger than
		//the six-vertex cap above.
		for (const std::array<int, 2>& e : non_tree) {
			std::vector<int> prev(V, -1);
			std::vector<char> visited(V, 0);
			std::queue<int> q;
			q.push(e[0]);
			visited[e[0]] = 1;
			while (!q.empty()) {
				const int u = q.front();
				q.pop();
				if (u == e[1])
					break;
				for (const int v : adj[u]) {
					if (visited[v])
						continue;
					if ((u == e[0] && v == e[1]) || (u == e[1] && v == e[0]))
						continue; //the cut edge
					visited[v] = 1;
					prev[v] = u;
					q.push(v);
				}
			}
			if (!visited[e[1]])
				continue; //no cycle through this edge after all
			std::vector<int> cycle;
			for (int at = e[1]; at != -1; at = prev[at]) cycle.push_back(at);
			if (cycle.size() >= 3)
				cycles.push_back(cycle);
		}
		return cycles;
	}

	bool covalently_bonded(const nucleus& a, const nucleus& b, const options& opt)
	{
		const double ra = a.Z > 0 && a.Z < 114 ? constants::covalent_radii[a.Z] : 1.5;
		const double rb = b.Z > 0 && b.Z < 114 ? constants::covalent_radii[b.Z] : 1.5;
		return array_length(a.pos, b.pos) <= constants::ang2bohr(opt.bond_scale * (ra + rb));
	}

	int covalent_fragment_count(const std::vector<nucleus>& nuclei, const options& opt)
	{
		const int V = (int)nuclei.size();
		if (V <= 1)
			return 1;
		//Union-find over the pair criterion.  Quadratic in the nuclei, which is what the bond seeding
		//already is, and it runs once
		std::vector<int> parent(V);
		for (int a = 0; a < V; a++) parent[a] = a;
		auto root = [&parent](int x) { while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; };
		for (int a = 0; a < V; a++)
			for (int b = a + 1; b < V; b++)
				if (covalently_bonded(nuclei[a], nuclei[b], opt)) {
					const int ra = root(a), rb = root(b);
					if (ra != rb) parent[ra] = rb;
				}
		int components = 0;
		for (int a = 0; a < V; a++) if (root(a) == a) components++;
		return components;
	}

	void tally(result& r, const options& opt)
	{
		r.n_attractor = r.n_bond = r.n_ring = r.n_cage = r.n_degenerate = r.n_nna = 0;
		for (const cp& point : r.points) {
			switch (point.kind) {
			case cp_kind::attractor: r.n_attractor++; if (point.is_nna) r.n_nna++; break;
			case cp_kind::bond: r.n_bond++; break;
			case cp_kind::ring: r.n_ring++; break;
			case cp_kind::cage: r.n_cage++; break;
			default: r.n_degenerate++; break;
			}
		}
		r.sum = r.n_attractor - r.n_bond + r.n_ring - r.n_cage;
		r.target = opt.poincare_hopf_target;
		//The molecular form's 1 is ONE isolated molecule's index sum.  C of them, with rho -> 0
		//between them, sum to C: tests/TFVC/water.gbw is a water with a helium 13.2 bohr away, its 4
		//attractors and 2 bond points sum to 2, and it was refused for a deficit of -1 while its own
		//diagnosis read "the bond graph falls into 2 covalent fragments".
		//
		//C is taken from the bond paths that were FOUND, not from the covalent geometry, because a
		//closed-shell contact makes a bond path where the covalent criterion sees none: a
		//hydrogen-bonded dimer is covalently two fragments and its sum is 1, so a geometric count
		//would refuse it.  What that costs is stated and paid for separately: dropping a bridging bond
		//point splits a component and raises the sum, moving both sides of sum == C together, so the
		//sum alone can no longer see it.  covalent_fragments is what sees it - the found graph must
		//not be MORE disconnected than the geometry, checked with graph_consistent below.
		//
		//Only the molecular default is raised.  An explicit poincare_hopf_target - the periodic Morse
		//form's 0 among them - is the caller's choice and is left alone.
		if (opt.poincare_hopf_target == 1 && r.graph_components > 1)
			r.target = r.graph_components;
		r.balanced = r.sum == r.target;
		r.found_ring_minus_cage = r.n_ring - r.n_cage;
		//r.complete is set below, once the bond-graph consistency it now depends on has been worked out

		std::ostringstream d;
		if (r.n_degenerate > 0)
			d << r.n_degenerate << " degenerate critical point(s): an eigenvalue of the Hessian is below "
			  << "tolerance there, so the point has no (3,s) signature and the Poincare-Hopf sum is "
			  << "undefined however it adds up. Tighten the wavefunction or inspect those points. ";
		if (!r.balanced) {
			const int deficit = r.target - r.sum;
			d << "Poincare-Hopf sum is " << r.sum << ", not " << r.target << " (deficit " << deficit << "). ";
			//Which class.  Adding an attractor or a ring point raises the sum by one, a bond or cage
			//point lowers it, so the sign of the deficit already halves the search
			if (deficit > 0)
				d << "The search is short of " << deficit << " attractor or ring point(s), or has that many spurious bond or cage points. ";
			else
				d << "The search is short of " << -deficit << " bond or cage point(s), or has that many spurious attractor or ring points. ";
		}
		//The bond-graph accounting used to sit inside the block above, so it was only ever consulted
		//once the sum had already failed - and the sum is a necessary condition, not a sufficient one.
		//A spurious bond point and a spurious ring point cancel in the alternating sum: HgH2 (3 nuclei,
		//2 bonds, no ring) came out at 3 - 4 + 2 - 0 = 1, balanced, and was reported COMPLETE with two
		//bond points and two ring points that cannot exist.  The graph rank is what catches that, and
		//the code already computed it.
		{
			const bool rings_ok = r.found_ring_minus_cage == r.required_ring_minus_cage;
			const bool attractors_ok = r.n_attractor == r.graph_vertices + r.n_nna;
			//The bond paths cannot be MORE disconnected than the covalent geometry: a closed-shell
			//contact joins fragments the covalent criterion keeps apart, never the other way round.
			//This is what pays for taking the target from the found graph - a missing bridging bond
			//point moves the sum and the component count together and so is invisible to the sum, but
			//it leaves one covalent fragment split in two, which is visible here.
			const bool connected_ok = r.graph_components <= r.covalent_fragments;
			r.graph_consistent = r.graph_vertices == 0 || (rings_ok && attractors_ok && connected_ok);
			if (r.graph_vertices > 0 && !(r.balanced && r.graph_consistent)) {
				d << "The bond graph has " << r.graph_vertices << " nuclei, " << r.graph_edges
				  << " bonds and " << r.graph_components << " component(s), cycle rank "
				  << r.required_ring_minus_cage << "; found n_ring - n_cage = " << r.found_ring_minus_cage << ". ";
				if (r.found_ring_minus_cage < r.required_ring_minus_cage)
					d << "Ring seeding is the likely gap: " << (r.required_ring_minus_cage - r.found_ring_minus_cage)
					  << " independent cycle(s) of the bond graph produced no ring point. ";
				else if (r.found_ring_minus_cage > r.required_ring_minus_cage)
					d << "More ring than cage points than the bond graph allows: either a spurious ring point, or "
					  << (r.found_ring_minus_cage - r.required_ring_minus_cage) << " cage point(s) were not found - cage seeding is the likely gap. ";
				if (r.n_attractor != r.graph_vertices + r.n_nna)
					d << "Attractors (" << r.n_attractor << ") do not match nuclei plus non-nuclear attractors ("
					  << r.graph_vertices << " + " << r.n_nna << "): nuclear seeding is the likely gap. ";
				if (r.graph_components > r.covalent_fragments)
					d << "The bond paths fall into " << r.graph_components << " component(s) where the covalent "
					  << "geometry has only " << r.covalent_fragments << ": " << (r.graph_components - r.covalent_fragments)
					  << " covalently bonded pair(s) carry no bond critical point, so a bridging bond point is "
					  << "missing - and that is the one gap the alternating sum cannot see, because losing a "
					  << "bridge raises the sum and splits a component at the same time. ";
			}
		}
		//Every sentence above blames a seeding class, and on an ECP wavefunction all of them are
		//blaming the search for something it cannot find.  Measured over the five topology cells:
		//ECP_SF/Au2Br2.gbw, ECP_SF/malbac.gbw and ELI_heavy/hgh2_ecp.gbw are INCOMPLETE and all three
		//carry coreless nuclei; Fe_gbw/Fe.gbw and TFVC/water.gbw are COMPLETE and carry none.  On
		//Au2Br2 the two nuclei with no attractor are exactly its two Au.
		//d is non-empty only when something above has already failed, which is the condition for
		//saying any of this: on a complete set there is no verdict to explain
		if (!r.coreless_nuclei.empty() && d.tellp() > 0) {
			d << r.coreless_nuclei.size() << " nucleus/nuclei carry no core density (rho below "
			  << opt.core_rho_floor << " at a nucleus of Z >= 5), atom number";
			for (size_t i = 0; i < r.coreless_nuclei.size() && i < 8; i++)
				d << (i ? "," : " ") << r.coreless_nuclei[i] + 1;
			if (r.coreless_nuclei.size() > 8)
				d << " and " << r.coreless_nuclei.size() - 8 << " more";
			d << ": a pseudopotential replaced the core there, so the density has no cusp at those "
			  << "nuclei and no readable shell structure around them. A missing nuclear attractor on "
			  << "one of them, or a spurious ring or cage point within about a bohr of one, is a "
			  << "property of this wavefunction and not a gap in the seeding - QTAIM is not defined on "
			  << "a pseudo-density inside the core radius. ";
		}
		r.complete = r.balanced && r.n_degenerate == 0 && r.graph_consistent;
		r.diagnosis = d.str();
	}

	void report_topology(const result& r, const std::vector<nucleus>& nuclei, std::ostream& log, const options& opt)
	{
		log << "\n---------------- Topological analysis of rho ----------------\n";
		citations::cite(citations::Method::QTAIM, log);
		log << "Newton-Raphson on grad rho = 0 with the analytic Hessian. Positions in bohr.\n"
			//A rejected point is as much a result as an accepted one, so the thresholds that did the
			//rejecting are printed: a critical point below the density floor would not appear above
			<< "Accepted when |grad rho| <= " << opt.gradient_tolerance << " * max(1, rho) and <= "
			<< opt.relative_gradient_tolerance << " * rho, with rho > " << opt.density_floor
			<< "; points merged below " << opt.merge_distance << " bohr; an eigenvalue below "
			<< opt.eigen_tolerance << " * max|lambda| counts as zero.\n\n"
			<< std::setw(4) << "#" << std::setw(19) << "type" << std::setw(11) << "seed"
			<< std::setw(11) << "x" << std::setw(11) << "y" << std::setw(11) << "z"
			<< std::setw(13) << "rho" << std::setw(13) << "lap(rho)" << std::setw(11) << "|grad|"
			<< std::setw(7) << "it" << std::setw(8) << "near" << std::setw(9) << "d(near)" << "\n";
		log << std::string(138, '-') << "\n";
		int i = 0;
		for (const cp& p : r.points) {
			log << std::setw(4) << ++i << std::setw(19) << kind_name(p.kind) << std::setw(11) << seed_class_name(p.from)
				<< std::fixed << std::setprecision(5)
				<< std::setw(11) << p.position[0] << std::setw(11) << p.position[1] << std::setw(11) << p.position[2]
				<< std::scientific << std::setprecision(4)
				<< std::setw(13) << p.density << std::setw(13) << p.laplacian << std::setw(11) << p.gradient_norm
				<< std::setw(7) << p.iterations;
			if (p.nearest_nucleus >= 0 && p.nearest_nucleus < (int)nuclei.size())
				log << std::setw(8) << (std::string(constants::atnr2letter(nuclei[p.nearest_nucleus].Z)) + std::to_string(p.nearest_nucleus + 1));
			else
				log << std::setw(8) << "-";
			log << std::fixed << std::setprecision(4) << std::setw(9) << p.nearest_nucleus_distance << "\n";
		}
		log << std::defaultfloat << std::setprecision(6) << "\n";

		log << "Counts: " << r.n_attractor << " attractor (3,-3), " << r.n_bond << " bond (3,-1), "
			<< r.n_ring << " ring (3,+1), " << r.n_cage << " cage (3,+3)";
		if (r.n_degenerate)
			log << ", " << r.n_degenerate << " DEGENERATE";
		log << "\n";
		//Molecular form.  A periodic density in a unit cell obeys the Morse relation, whose sum is 0
		log << "Poincare-Hopf (molecular form, isolated molecule): n_NCP - n_BCP + n_RCP - n_CCP = "
			<< r.n_attractor << " - " << r.n_bond << " + " << r.n_ring << " - " << r.n_cage
			<< " = " << r.sum << " (expected " << r.target << ")\n";
		//Where a target above 1 comes from, said out loud: an unexplained "expected 2" is worse than
		//the wrong "expected 1" it replaced
		if (r.target > 1 && opt.poincare_hopf_target == 1)
			log << "The bond paths fall into " << r.graph_components << " disconnected component(s), so the "
				<< "expected sum is " << r.target << " and not 1: each separated fragment contributes 1. The "
				<< "covalent geometry has " << r.covalent_fragments << " fragment(s), and a component count "
				<< "above that would mean a bond path is missing rather than that the system is separated.\n";
		if (r.escalated) {
			int from_grid = 0;
			for (const cp& p : r.points) if (p.from == seed_class::grid) from_grid++;
			log << "Topology-driven seeding did not close the sum; a coarse grid of extra seeds was added and "
				<< "contributed " << from_grid << " point(s).\n";
		}
		if (r.complete)
			log << "The set of critical points is COMPLETE: the relation holds, no point is degenerate, and "
				"n_ring - n_cage matches the bond graph's cycle rank of " << r.required_ring_minus_cage << ".\n";
		else
			log << "The set of critical points is INCOMPLETE. " << r.diagnosis << "\n";

		if (r.n_nna) {
			log << "\nNon-nuclear attractors: " << r.n_nna << ". A (3,-3) point away from every nucleus owns a\n"
				<< "basin of its own, so any basin integration that ignores it mis-assigns that charge.\n";
			for (const cp& p : r.points) {
				if (!p.is_nna)
					continue;
				log << "  NNA at " << std::fixed << std::setprecision(5) << p.position[0] << " " << p.position[1] << " " << p.position[2]
					<< "  rho = " << std::scientific << std::setprecision(6) << p.density
					<< "  nearest nucleus " << std::defaultfloat;
				if (p.nearest_nucleus >= 0 && p.nearest_nucleus < (int)nuclei.size())
					log << constants::atnr2letter(nuclei[p.nearest_nucleus].Z) << p.nearest_nucleus + 1;
				log << " at " << std::fixed << std::setprecision(4) << p.nearest_nucleus_distance << " bohr\n";
			}
			log << std::defaultfloat << std::setprecision(6);
		}
		log << "------------------------------------------------------------\n";
	}

	bool report(const std::filesystem::path& wfn_path, std::ostream& log)
	{
		err_checkf(std::filesystem::exists(wfn_path), "Could not find " + wfn_path.string(), log);
		WFN wavy(wfn_path);
		const std::vector<nucleus> n = nuclei_of(wavy);
		err_checkf(!n.empty(), "No nuclei in " + wfn_path.string(), log);
		log << "Topological analysis of " << wfn_path.string() << ": " << n.size() << " nuclei\n";
		const options opt{};
		const result r = analyze_topology(wavy, n, opt);
		report_topology(r, n, log, opt);
		//the verdict is a result, not decoration: the caller exits non-zero when it is false
		return r.complete;
	}
}
