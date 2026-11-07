#pragma once
//All critical points of a density and the Poincare-Hopf completeness test.  Unlike b2c.cpp (cube
//voxels, no ring or cage seeding) the search is cube-free and seeded from the topology: nuclei,
//bonded pairs, cycles of the bond-path graph, and centroids of the ring points those give.
//Newton-Raphson on grad rho = 0 with the analytic Hessian, because Newton converges on a stationary
//point of any index while gradient following only reaches maxima.
//Signature from the eigenvalue signs: (3,-3) attractor (an NNA away from every nucleus), (3,-1) bond,
//(3,+1) ring, (3,+3) cage.  An eigenvalue below tolerance makes the point degenerate, never binned,
//since the Poincare-Hopf sum is then undefined.
//Poincare-Hopf, molecular (rho -> 0 at infinity): n_NCP - n_BCP + n_RCP - n_CCP = 1; the periodic
//Morse form over a cell sums to 0.
#include "convenience.h"
#include "density_source.h"
#include "nos_math.h"
#include "b2c.h" //critical_point
#include <iosfwd>
#include <string>
#include <vector>

namespace topology
{
	struct nucleus
	{
		d3 pos;
		int Z;
	};

	enum class cp_kind
	{
		attractor, //(3,-3)
		bond,      //(3,-1)
		ring,      //(3,+1)
		cage,      //(3,+3)
		degenerate //an eigenvalue below tolerance
	};
	const char* kind_name(const cp_kind k);
	//Seeding class that produced a point, so a gap in the sum can be blamed on the class that is short
	enum class seed_class
	{
		nucleus,
		bond,
		ring,
		cage,
		grid //escalation fallback when the sum does not close
	};
	const char* seed_class_name(const seed_class s);

	struct cp
	{
		d3 position{ 0.0, 0.0, 0.0 };
		d3 gradient{ 0.0, 0.0, 0.0 };
		d3 eigenvalues{ 0.0, 0.0, 0.0 }; //ascending
		double density = 0.0;
		double laplacian = 0.0;
		double gradient_norm = 0.0;
		double ellipticity = std::numeric_limits<double>::quiet_NaN(); //lambda1/lambda2 - 1, bond points only
		cp_kind kind = cp_kind::degenerate;
		seed_class from = seed_class::grid;
		int negative = 0, positive = 0, zero = 0;
		int iterations = 0;
		int nearest_nucleus = -1;
		double nearest_nucleus_distance = 0.0;
		bool is_nna = false; //attractor further than options::nna_distance from every nucleus
	};

	struct options
	{
		//|grad rho| at an accepted point, a.u., scaled by rho where rho > 1 (see newton_to_critical_point)
		double gradient_tolerance = 1E-7;
		//In the density tail |grad rho| is of the order of rho, so an absolute tolerance accepts any
		//tail seed without a step; |grad rho| / rho is large there and zero at a real critical point
		double relative_gradient_tolerance = 1E-3; //|grad rho| <= this * rho, bohr^-1
		double trust_radius = 0.3;                 //longest Newton step, bohr
		int max_iterations = 80;
		double merge_distance = 0.05;  //two accepted points closer than this are one point, bohr
		double nna_distance = 0.3;     //an attractor beyond this from every nucleus is an NNA, bohr
		//Below this rho is vacuum, far under a real cage point (~1E-3 a.u.)
		double density_floor = 1E-7;
		//rho at a Z >= 5 nucleus below this means an ECP removed the core: all-electron boron already
		//has rho ~ 40 there and the cusp grows as Z^3, a pseudo-density is orders of magnitude lower
		double core_rho_floor = 1.0;
		double bond_scale = 1.3;       //pair is bonded if d <= bond_scale * (r_cov,a + r_cov,b)
		double eigen_tolerance = 1E-6; //relative to max|lambda|: below it an eigenvalue is zero
		bool escalate_on_mismatch = true; //add a coarse grid of seeds when the sum does not close
		int escalation_grid = 7;          //points per axis over the nuclear bounding box
		int poincare_hopf_target = 1;     //1 molecular, 0 for the periodic Morse form
	};

	struct result
	{
		std::vector<cp> points;
		int n_attractor = 0, n_bond = 0, n_ring = 0, n_cage = 0, n_degenerate = 0, n_nna = 0;
		int sum = 0;                //n_NCP - n_BCP + n_RCP - n_CCP
		int target = 1;
		bool balanced = false;      //sum == target
		bool complete = false;      //balanced and no degenerate point
		int graph_vertices = 0;
		int graph_edges = 0;        //bond critical points that link two nuclei
		int graph_components = 0;
		int required_ring_minus_cage = 0; //E - V + C
		int found_ring_minus_cage = 0;
		//n_ring - n_cage equals the bond graph's cycle rank and attractors equal nuclei plus NNAs.
		//Poincare-Hopf is necessary, not sufficient: a spurious bond and ring point cancel in the sum
		bool graph_consistent = true;
		//Components of the covalent graph from geometry alone.  The found bond-path graph may not have
		//more: a missing bridging bond point raises both the sum and its target, invisible to the sum
		int covalent_fragments = 1;
		bool escalated = false;
		//Z >= 5 nuclei with rho below options::core_rho_floor: an ECP pseudo-density has no cusp and
		//no shell structure there, so missing or spurious points near them are not a seeding gap
		std::vector<int> coreless_nuclei;
		std::string diagnosis;  //empty when complete, otherwise what is missing and where
	};

	//d <= bond_scale * (r_cov,a + r_cov,b).  Shared by bond seeding and fragment count, which must agree
	bool covalently_bonded(const nucleus& a, const nucleus& b, const options& opt);
	//Components of that graph; an empty system gives 1, not a deficit
	int covalent_fragment_count(const std::vector<nucleus>& nuclei, const options& opt);

	//Signature from a row-major 3x3 Hessian; eigen_tolerance is relative to the largest |lambda|
	cp_kind classify_hessian(const double H[9], d3& eigenvalues, int& negative, int& positive, int& zero, const double eigen_tolerance = 1E-6);

	//Poincare-Hopf accounting and diagnosis from the counts alone
	void tally(result& r, const options& opt = {});

	std::vector<nucleus> nuclei_of(const WFN& wavy);

	//Helpers declared before the templates: non-dependent calls are looked up at definition, not instantiation.
	//Closed-form inverse; false (inv untouched) when m is singular
	bool invert_3x3_topology(const double m[9], double inv[9]);
	d3 centroid(const std::vector<d3>& p);
	//Smallest independent cycles of the bond-critical-point graph as nucleus index lists; fills r's graph counts
	std::vector<std::vector<int>> bond_graph_cycles(result& r, const std::vector<nucleus>& nuclei);

	template <class S>
	cp describe_cp(const S& source, const d3& p, const std::vector<nucleus>& nuclei, const options& opt)
	{
		cp c;
		c.position = p;
		double H[9]{};
		c.density = calculate_hessian(source, p, c.gradient, H);
		c.gradient_norm = array_length(c.gradient);
		c.laplacian = H[0] + H[4] + H[8];
		c.kind = classify_hessian(H, c.eigenvalues, c.negative, c.positive, c.zero, opt.eigen_tolerance);
		if (c.kind == cp_kind::bond && std::abs(c.eigenvalues[1]) > 0.0)
			c.ellipticity = c.eigenvalues[0] / c.eigenvalues[1] - 1.0;
		c.nearest_nucleus_distance = std::numeric_limits<double>::infinity();
		for (size_t a = 0; a < nuclei.size(); a++) {
			const double d = array_length(p, nuclei[a].pos);
			if (d < c.nearest_nucleus_distance) { c.nearest_nucleus_distance = d; c.nearest_nucleus = (int)a; }
		}
		c.is_nna = c.kind == cp_kind::attractor && c.nearest_nucleus_distance > opt.nna_distance;
		return c;
	}

	//Measured from rho at the nucleus, not a flag: the file does not always record an ECP
	template <class S>
	void note_coreless_nuclei(result& r, const S& source, const std::vector<nucleus>& nuclei, const options& opt)
	{
		r.coreless_nuclei.clear();
		for (size_t a = 0; a < nuclei.size(); a++) {
			if (nuclei[a].Z < 5)
				continue;
			d3 grad{ 0.0, 0.0, 0.0 };
			double H[9]{};
			if (calculate_hessian(source, nuclei[a].pos, grad, H) < opt.core_rho_floor)
				r.coreless_nuclei.push_back((int)a);
		}
	}

	//S is any density_source.h source implementing hessian(p, grad, H).
	//Step -H^-1 grad, capped at opt.trust_radius and halved until |grad rho| does not grow; accepted
	//only when the gradient vanishes both absolutely and relative to rho
	template <class S>
	bool newton_to_critical_point(const S& source, d3& p, const options& opt, int& iterations)
	{
		//The absolute bound scales with rho above 1: at a heavy nucleus (rho 1E3-1E4, Hessian 1E8-1E9)
		//cancellation among primitive terms leaves a gradient floor far above 1E-7 and nuclei go missing
		auto converged = [&opt](const double gnorm, const double rho) {
			return rho > opt.density_floor && gnorm <= opt.gradient_tolerance * (rho > 1.0 ? rho : 1.0)
				&& gnorm <= opt.relative_gradient_tolerance * rho;
		};
		d3 grad{ 0.0, 0.0, 0.0 };
		double H[9]{};
		double rho = calculate_hessian(source, p, grad, H);
		double gnorm = array_length(grad);
		iterations = 0;
		for (; iterations < opt.max_iterations; iterations++) {
			if (!std::isfinite(gnorm))
				return false;
			if (converged(gnorm, rho))
				return true;
			if (rho <= opt.density_floor)
				return false;
			double Hinv[9];
			if (!invert_3x3_topology(H, Hinv))
				return false;
			d3 step{
				-(Hinv[0] * grad[0] + Hinv[1] * grad[1] + Hinv[2] * grad[2]),
				-(Hinv[3] * grad[0] + Hinv[4] * grad[1] + Hinv[5] * grad[2]),
				-(Hinv[6] * grad[0] + Hinv[7] * grad[1] + Hinv[8] * grad[2])
			};
			if (!std::isfinite(step[0]) || !std::isfinite(step[1]) || !std::isfinite(step[2]))
				return false;
			const double snorm = array_length(step);
			if (snorm > opt.trust_radius)
				for (int k = 0; k < 3; k++) step[k] *= opt.trust_radius / snorm;

			bool accepted = false;
			for (int attempt = 0; attempt < 12 && !accepted; attempt++) {
				const double damping = std::pow(0.5, attempt);
				const d3 q{ p[0] + damping * step[0], p[1] + damping * step[1], p[2] + damping * step[2] };
				d3 tg{ 0.0, 0.0, 0.0 };
				double tH[9]{};
				const double trho = calculate_hessian(source, q, tg, tH);
				const double tn = array_length(tg);
				if (!std::isfinite(tn) || tn > gnorm)
					continue;
				p = q;
				grad = tg;
				gnorm = tn;
				rho = trho;
				std::copy(tH, tH + 9, H);
				accepted = true;
			}
			if (!accepted)
				return false; //stalled: not a critical point
		}
		return converged(gnorm, rho);
	}

	template <class S>
	result analyze_topology(const S& source, const std::vector<nucleus>& nuclei, const options& opt = {})
	{
		result r;
		std::vector<d3> seeds;
		std::vector<seed_class> seed_of;
		auto add_seed = [&](const d3& p, const seed_class c) { seeds.push_back(p); seed_of.push_back(c); };

		auto run = [&](const size_t first) {
			//Searches in parallel, merge serial in seed order so the point set is deterministic
			const int n = static_cast<int>(seeds.size() - first);
			std::vector<cp> found(n);
			std::vector<char> ok(n, 0);
#pragma omp parallel for schedule(dynamic, 1)
			for (int i = 0; i < n; i++) {
				d3 p = seeds[first + i];
				int it = 0;
				if (!newton_to_critical_point(source, p, opt, it))
					continue;
				found[i] = describe_cp(source, p, nuclei, opt);
				found[i].from = seed_of[first + i];
				found[i].iterations = it;
				ok[i] = 1;
			}
			for (int i = 0; i < n; i++) {
				if (!ok[i])
					continue;
				const cp &point = found[i];
				bool duplicate = false;
				//Same kind first: an attractor and a saddle within merge_distance are two points, and
				//merging by distance alone lets a bond point replace a nuclear attractor
				for (cp& existing : r.points)
					if (existing.kind == point.kind && array_length(existing.position, point.position) <= opt.merge_distance) {
						duplicate = true;
						if (point.gradient_norm < existing.gradient_norm) {
							const seed_class keep = existing.from; //the class that found it first
							existing = point;
							existing.from = keep;
						}
						break;
					}
				if (!duplicate)
					r.points.push_back(point);
			}
		};

		for (const nucleus& n : nuclei) add_seed(n.pos, seed_class::nucleus);
		run(0);

		//Several seeds per bond, not the midpoint: a polar bond's point sits off centre and the trust radius is short
		{
			const size_t first = seeds.size();
			for (size_t a = 0; a < nuclei.size(); a++)
				for (size_t b = a + 1; b < nuclei.size(); b++) {
					if (!covalently_bonded(nuclei[a], nuclei[b], opt))
						continue;
					for (int i = 1; i < 10; i++) {
						const double t = 0.1 * i;
						add_seed({ nuclei[a].pos[0] + t * (nuclei[b].pos[0] - nuclei[a].pos[0]),
								   nuclei[a].pos[1] + t * (nuclei[b].pos[1] - nuclei[a].pos[1]),
								   nuclei[a].pos[2] + t * (nuclei[b].pos[2] - nuclei[a].pos[2]) },
							seed_class::bond);
					}
				}
			run(first);
		}

		//Ring seeds: one per non-tree edge of a spanning forest, at the centroid of the smallest cycle through it
		{
			const size_t first = seeds.size();
			for (const std::vector<int>& cycle : bond_graph_cycles(r, nuclei)) {
				d3 c{ 0.0, 0.0, 0.0 };
				for (const int a : cycle) for (int k = 0; k < 3; k++) c[k] += nuclei[a].pos[k];
				for (int k = 0; k < 3; k++) c[k] /= (double)cycle.size();
				add_seed(c, seed_class::ring);
			}
			run(first);
		}

		//Cage seeds: centroid of all ring points (one cage) and of each with its three nearest (fused cages)
		{
			const size_t first = seeds.size();
			std::vector<d3> rings;
			for (const cp& point : r.points) if (point.kind == cp_kind::ring) rings.push_back(point.position);
			if (rings.size() >= 4) {
				add_seed(centroid(rings), seed_class::cage);
				for (size_t i = 0; i < rings.size(); i++) {
					std::vector<size_t> order(rings.size());
					for (size_t j = 0; j < order.size(); j++) order[j] = j;
					std::sort(order.begin(), order.end(), [&](const size_t x, const size_t y) {
						return array_length(rings[i], rings[x]) < array_length(rings[i], rings[y]);
					});
					std::vector<d3> quad;
					for (int k = 0; k < 4; k++) quad.push_back(rings[order[k]]);
					add_seed(centroid(quad), seed_class::cage);
				}
			}
			run(first);
		}

		r.covalent_fragments = covalent_fragment_count(nuclei, opt);
		note_coreless_nuclei(r, source, nuclei, opt);
		tally(r, opt);

		//A coarse grid over the nuclear bounding box assumes nothing about what is missing
		if (!r.balanced && opt.escalate_on_mismatch && opt.escalation_grid > 1 && !nuclei.empty()) {
			const size_t first = seeds.size();
			d3 lo = nuclei.front().pos, hi = nuclei.front().pos;
			for (const nucleus& n : nuclei)
				for (int k = 0; k < 3; k++) { lo[k] = std::min(lo[k], n.pos[k]); hi[k] = std::max(hi[k], n.pos[k]); }
			for (int k = 0; k < 3; k++) { lo[k] -= 1.0; hi[k] += 1.0; }
			const int n = opt.escalation_grid;
			for (int ix = 0; ix < n; ix++)
				for (int iy = 0; iy < n; iy++)
					for (int iz = 0; iz < n; iz++)
						add_seed({ lo[0] + (hi[0] - lo[0]) * ix / (n - 1.0),
								   lo[1] + (hi[1] - lo[1]) * iy / (n - 1.0),
								   lo[2] + (hi[2] - lo[2]) * iz / (n - 1.0) },
							seed_class::grid);
			run(first);
			r.escalated = true;
			//The grid can add bond points, so the graph counts are rebuilt before the second tally
			bond_graph_cycles(r, nuclei);
			tally(r, opt);
		}
		return r;
	}

	void report_topology(const result& r, const std::vector<nucleus>& nuclei, std::ostream& log, const options& opt = {});
	//-topology <wfn>: true when the critical point set is complete, so the caller can exit non-zero
	bool report(const std::filesystem::path& wfn_path, std::ostream& log);
}
