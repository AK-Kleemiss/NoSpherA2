#pragma once
//The complete set of critical points of a density and the one test that says whether it is
//complete: Poincare-Hopf.  b2c.cpp finds critical points from a cube's voxels and seeds along
//bonds, which covers nuclei and bond critical points; it has no ring or cage seeding at all, so
//a ring critical point only turns up when a voxel happens to bracket it.  Here the search is
//cube-free and seeded from the topology itself: nuclei, bonded pairs, the cycles of the bond
//graph the bond critical points define, and the centroids of the ring critical points those
//cycles produce.
//
//Newton-Raphson on grad rho = 0 with the analytic Hessian from density_source.h.  Newton is what
//finds a saddle: it converges on a stationary point of any index, where following the gradient
//can only reach a maximum.  Nothing here uses a finite difference.
//
//Signature by the signs of the three eigenvalues of grad grad rho at a point where grad rho = 0:
//  (3,-3) attractor - a maximum, at a nucleus normally; away from every nucleus a non-nuclear
//                     attractor (NNA), which owns a basin of its own in any integration
//  (3,-1) bond       - two negative, one positive
//  (3,+1) ring       - one negative, two positive
//  (3,+3) cage       - all positive, a minimum
//An eigenvalue whose modulus falls below the tolerance makes the point degenerate.  A degenerate
//point is reported as such and never binned into one of the four: the Poincare-Hopf sum is then
//undefined, and saying so is worth more than a sum that happens to come out right.
//
//Poincare-Hopf, molecular form (an isolated molecule, rho -> 0 at infinity):
//    n_NCP - n_BCP + n_RCP - n_CCP = 1
//The periodic (Morse) form over a unit cell sums to 0 instead.  NoSpherA2 reads molecular
//wavefunctions, so 1 is what analyze_topology assumes and what report_topology prints; the
//assumption is named in the output.
#include "convenience.h"
#include "density_source.h"
#include "nos_math.h"
#include "b2c.h" //critical_point
#include <iosfwd>
#include <string>
#include <vector>

namespace topology
{
	//A nucleus the search seeds from and measures NNA distances against
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
		degenerate //an eigenvalue below tolerance: not classifiable
	};
	const char* kind_name(const cp_kind k);
	//Which seeding class produced a point, so a gap in the sum can be blamed on the class that
	//is short rather than reported as a bare wrong number
	enum class seed_class
	{
		nucleus,
		bond,
		ring,
		cage,
		grid //the escalation fallback, only used when the sum does not close
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
		bool is_nna = false; //an attractor further than options::nna_distance from every nucleus
	};

	struct options
	{
		double gradient_tolerance = 1E-7; //|grad rho| at an accepted point, atomic units
		//An absolute gradient tolerance alone is not a convergence test.  Far out in the density
		//tail |grad rho| is itself of the order of rho, so |grad rho| <= 1E-7 holds at every point
		//out there and a seed placed in the tail reports a critical point without taking a step -
		//measured on an octahedral test density, an unguarded coarse-grid escalation invented 128
		//ring points and a Poincare-Hopf sum of 121.  The guard is the dimensionless logarithmic
		//derivative |grad rho| / rho (bohr^-1), which is large in the tail and zero at a real
		//critical point; with it the same six test densities all give the exact topology.
		double relative_gradient_tolerance = 1E-3; //|grad rho| <= this * rho, bohr^-1
		double trust_radius = 0.3;                 //longest Newton step, bohr
		int max_iterations = 80;
		double merge_distance = 0.05;  //two accepted points closer than this are one point, bohr
		double nna_distance = 0.3;     //an attractor beyond this from every nucleus is an NNA, bohr
		//Below this rho is vacuum.  Four orders of magnitude under the rho of a real cage point
		//(~1E-3 a.u.), so it cannot hide one; the relative gradient above is the actual guard
		double density_floor = 1E-7;
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
		int target = 1;             //what it is compared against
		bool balanced = false;      //sum == target
		bool complete = false;      //balanced and no degenerate point
		int graph_vertices = 0;     //nuclei the bond graph connects
		int graph_edges = 0;        //bond critical points that link two nuclei
		int graph_components = 0;
		int required_ring_minus_cage = 0; //E - V + C: what Poincare-Hopf needs of the two classes
		int found_ring_minus_cage = 0;
		bool escalated = false; //a coarse grid of seeds was added
		std::string diagnosis;  //empty when complete, otherwise what is missing and where
	};

	//Classification of a stationary point from its Hessian alone (row-major 3x3).  Split out so
	//the signature logic can be checked against a matrix whose eigenvalues are known exactly,
	//without a density anywhere near it.  eigen_tolerance is relative to the largest |lambda|.
	cp_kind classify_hessian(const double H[9], d3& eigenvalues, int& negative, int& positive, int& zero, const double eigen_tolerance = 1E-6);

	//Poincare-Hopf accounting and the diagnosis, from counts alone.  Public so the reporting can
	//be checked without running a search.
	void tally(result& r, const options& opt = {});

	std::vector<nucleus> nuclei_of(const WFN& wavy);

	//----- helpers.  Declared before the templates that call them: the calls take plain doubles and
	//vectors, so they are looked up where the template is defined, not where it is instantiated ---

	//Closed-form 3x3 inverse for the Newton step; false (inv untouched) when m is singular
	bool invert_3x3_topology(const double m[9], double inv[9]);
	d3 centroid(const std::vector<d3>& p);
	//Smallest independent cycles of the graph whose edges are the accepted bond critical points,
	//each as a list of nucleus indices.  Also fills r's graph counts.
	std::vector<std::vector<int>> bond_graph_cycles(result& r, const std::vector<nucleus>& nuclei);

	//rho, its gradient, the Hessian signature and the distance to the nearest nucleus at p
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

	//----- the templated search.  S is any density_source.h source: a WFN, a fitted density, a
	//Centred<> atom model, or a test source that implements hessian(p, grad, H) -----

	//Newton-Raphson from p onto a stationary point of rho.  The step is -H^-1 grad, capped at
	//opt.trust_radius and backtracked by halving until |grad rho| does not grow; the point is
	//accepted only when the gradient has vanished both absolutely and relative to rho, so neither a
	//Newton iterate that stalls against the trust radius nor a seed sitting in the density tail is
	//returned as a critical point.
	template <class S>
	bool newton_to_critical_point(const S& source, d3& p, const options& opt, int& iterations)
	{
		auto converged = [&opt](const double gnorm, const double rho) {
			return rho > opt.density_floor && gnorm <= opt.gradient_tolerance && gnorm <= opt.relative_gradient_tolerance * rho;
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
				return false; //stalled: not a critical point, whatever the iterate looks like
		}
		return converged(gnorm, rho);
	}

	//The whole analysis: seed, converge, classify, count, check Poincare-Hopf.  nuclei supplies
	//both the nuclear seeds and the reference the NNA test measures against.
	template <class S>
	result analyze_topology(const S& source, const std::vector<nucleus>& nuclei, const options& opt = {})
	{
		result r;
		std::vector<d3> seeds;
		std::vector<seed_class> seed_of;
		auto add_seed = [&](const d3& p, const seed_class c) { seeds.push_back(p); seed_of.push_back(c); };

		auto run = [&](const size_t first) {
			for (size_t s = first; s < seeds.size(); s++) {
				d3 p = seeds[s];
				int it = 0;
				if (!newton_to_critical_point(source, p, opt, it))
					continue;
				cp point = describe_cp(source, p, nuclei, opt);
				point.from = seed_of[s];
				point.iterations = it;
				bool duplicate = false;
				for (cp& existing : r.points)
					if (array_length(existing.position, point.position) <= opt.merge_distance) {
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

		//1. nuclei: every nucleus is an attractor seed
		for (const nucleus& n : nuclei) add_seed(n.pos, seed_class::nucleus);
		run(0);

		//2. bonded pairs.  Several points along the vector, not just the midpoint: a polar bond's
		//critical point sits well off centre and the trust radius is short
		{
			const size_t first = seeds.size();
			for (size_t a = 0; a < nuclei.size(); a++)
				for (size_t b = a + 1; b < nuclei.size(); b++) {
					const double ra = nuclei[a].Z > 0 && nuclei[a].Z < 114 ? constants::covalent_radii[nuclei[a].Z] : 1.5;
					const double rb = nuclei[b].Z > 0 && nuclei[b].Z < 114 ? constants::covalent_radii[nuclei[b].Z] : 1.5;
					if (array_length(nuclei[a].pos, nuclei[b].pos) > constants::ang2bohr(opt.bond_scale * (ra + rb)))
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

		//3. rings: the cycles of the bond graph the bond critical points just built.  One seed per
		//independent cycle (a spanning forest's non-tree edges), at the centroid of the smallest
		//cycle through that edge - which is the ring the critical point belongs to
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

		//4. cages: centroids of the ring critical points.  A cage is enclosed by its rings, so the
		//centroid of all of them and of each ring point with its three nearest neighbours cover
		//both a single cage and several fused ones
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

		tally(r, opt);

		//5. the sum did not close.  A coarse grid over the nuclear bounding box is the only seeding
		//that makes no assumption about what is missing; whether it closed the gap is reported
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
			tally(r, opt);
		}
		return r;
	}

	//Every critical point of a wavefunction's density with the Poincare-Hopf verdict, printed
	void report_topology(const result& r, const std::vector<nucleus>& nuclei, std::ostream& log, const options& opt = {});
	//-topology <wfn>: load, analyze, report
	void report(const std::filesystem::path& wfn_path, std::ostream& log);
}
