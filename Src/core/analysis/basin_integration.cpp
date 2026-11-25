#include "pch.h"
#include "tuning.h"
#include "wfn_class.h"
#include "convenience.h"
#include "constants.h"
#include "b2c.h"
#include "nos_math.h"
#include "GridManager.h"
#include "fiber_loop.h"
#include "eli_family.h"
#include <limits>

static int g_basin_gpu = -1;
void basin_gpu_set_mode(const int mode) { g_basin_gpu = mode < 0 ? -1 : mode > 0; }
//Auto keeps the host below this cost per field evaluation: a device-fed thread pays a fixed host cost per evaluation
//(park, batch, copies, driver contention) that the probe does not see (-tune NOS_BASIN_GPU_MIN_US overrides).
//ponytail: one break-even taken at 24 host threads; fewer threads break even lower, a per-thread-count model if it matters
static constexpr double basin_gpu_min_us = 15.0;

//One per host thread of a device-fed loop. Each fiber (fiber_loop) holds a trajectory and parks its evaluations here;
//once all have parked or finished the round goes to the device as one batch on this feed's own context and stream.
struct field_feed {
	struct req { d3 p; double v = 0.0, rho = 0.0; d3 g{}; bool ok = false; };
	const WFN &wavy;
	const bool eli;
	void *ctx;
	bool failed = false;
	std::vector<req *> reqs;
	vec pts, val, grad, rho;
	long long batches = 0, points = 0;
	double device_s = 0.0; //inside field_gpu_run; the rest of a round is the host's turn
	const int fibers; //per host thread, at most one batch
	field_feed(const WFN &w, const bool e, const int max_points) : wavy(w), eli(e), ctx(w.field_gpu_open(e, max_points)), fibers(max_points) {}
	~field_feed() { WFN::field_gpu_close(ctx); }
	field_feed(const field_feed &) = delete;
	field_feed &operator=(const field_feed &) = delete;
	void flush()
	{
		const size_t n = reqs.size();
		if (n == 0) return;
		pts.resize(3 * n); val.resize(n); grad.resize(3 * n); rho.resize(n);
		for (size_t i = 0; i < n; i++) std::copy_n(reqs[i]->p.data(), 3, &pts[3 * i]);
		const auto t0 = std::chrono::steady_clock::now();
		if (!failed && !WFN::field_gpu_run(ctx, static_cast<int>(n), pts.data(), eli ? val.data() : nullptr, grad.data(), rho.data()))
			failed = true;
		device_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
		if (!failed)
			for (size_t i = 0; i < n; i++) {
				req &r = *reqs[i];
				r.v = val[i]; r.rho = rho[i]; std::copy_n(&grad[3 * i], 3, r.g.data()); r.ok = true;
			}
		batches++;
		points += static_cast<long long>(n);
		reqs.clear();
	}
	//ELI-D value (eli only), gradient and density at p, as computeELIGrad / computeGrad give them
	void eval(const d3 &p, double &v, d3 &g, double &rho_p)
	{
		req r;
		r.p = p;
		if (!failed) {
			reqs.push_back(&r);
			fiber_yield();
		}
		if (r.ok) { v = r.v; g = r.g; rho_p = r.rho; return; }
		//the device declined: every fiber evaluates its own point from here on
		if (eli) wavy.computeELIGrad(p, v, g, &rho_p);
		else wavy.computeGrad(p, g, &rho_p);
	}
};
//The feed of the calling thread's fibers; null outside a device-fed loop
static thread_local field_feed *t_batch = nullptr;
//One feed per host thread for a device-fed loop, or none; why gets the reason
static std::vector<std::unique_ptr<field_feed>> open_field_feeds(const WFN &wavy, bool eli, const vec &sample,
	const std::function<void(const d3 &)> &host_eval, int threads, std::ostringstream &why);
//body(i) for i in [0, n) on the calling thread's fibers; every thread of a parallel region calls it
static void feed_loop(field_feed &f, std::atomic<int> &next, const int n, const std::function<void(int)> &body)
{
	t_batch = &f;
	fiber_loop(f.fibers, next, n, body, [&f]() { f.flush(); });
	t_batch = nullptr;
}

//Monotone ascent of a spin-resolved ELI-D on rho >= basin_density_cutoff. Where one orbital carries a spin's density the
//same-spin field rises out to the isosurface, so the attractor is its maximum there: a step leaving the domain is pulled
//back along grad rho and the walk slides along the surface. False if p cannot be put on the domain.
bool bounded_eli_ascent(const WFN &wavy, const eli_spin_field &eli, d3 &p, double &f)
{
	//just inside, so the integration's own isosurface test accepts the attractor it lands on
	const double target = 1.001 * basin_density_cutoff;
	//Newton on ln rho, which the exponential tail makes nearly linear; a plain Newton on rho overshoots
	//inward and leaves the walk off the surface it is meant to slide along
	auto onto = [&](d3 &q) {
		for (int k = 0; k < 12; k++) {
			d3 n;
			double rho = 0.0;
			wavy.computeGrad(q, n, &rho);
			if (rho >= basin_density_cutoff) return true;
			const double n2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
			if (!(n2 > 1e-30) || !(rho > 0.0)) return false;
			const double s = std::log(target / rho) * rho / n2;
			for (int i = 0; i < 3; i++) q[i] += s * n[i];
		}
		return false;
	};
	if (!onto(p)) return false;
	d3 g;
	eli(p, f, g, nullptr);
	if (!std::isfinite(f)) return false;
	double step = 0.05;
	for (int it = 0; it < 500 && step > 1e-5; it++) {
		//On the surface with the field pointing out, only its tangential part moves the walk; the full gradient projected back stops short
		d3 d = g, n;
		double rho = 0.0;
		wavy.computeGrad(p, n, &rho);
		const double gdn = g[0] * n[0] + g[1] * n[1] + g[2] * n[2], n2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
		if (rho < 1.05 * basin_density_cutoff && gdn > 0.0 && n2 > 1e-30)
			for (int i = 0; i < 3; i++) d[i] -= gdn / n2 * n[i];
		const double gn = array_length(d);
		if (!std::isfinite(gn) || gn < 1e-12 || gn < 1e-9 * array_length(g)) break;
		d3 t{ p[0] + step * d[0] / gn, p[1] + step * d[1] / gn, p[2] + step * d[2] / gn };
		if (!onto(t)) { step *= 0.5; continue; }
		double ft;
		d3 gt;
		eli(t, ft, gt, nullptr);
		if (std::isfinite(ft) && ft > f) { p = t; f = ft; g = gt; step = std::min(1.5 * step, 0.3); }
		else step *= 0.5;
	}
	return true;
}

//One linearly implicit Euler step of the gradient flow dx/dt = g: (I - tau H) dx = tau g. A TM outer-core shell is a
//curved ridge whose radial curvature is ~1e4 times the tangential one (Ni(CO)4: 1.4092-1.4145 over the whole shell), so an
//explicit step longer than g_t/H_rr overshoots the crest and the walk stops rising on a steep slope. The implicit step
//lands on the crest along the stiff direction and moves tau g along the soft ones; tau -> infinity is the Newton step onto
//a maximum. H from central differences of the analytic gradient. tau is the largest (doubling from len/|g|) whose step
//stays within len and keeps I - tau H positive definite, then halved until the field rises. vg(p, g) returns the value.
//False when no step rises: p is a maximum to the field's precision.
template <class VG>
static bool implicit_ascent_step(const VG &vg, d3 &p, double &f, d3 &g, const double len)
{
	const double gn = array_length(g);
	if (!(gn > 0.0) || !std::isfinite(gn)) return false;
	constexpr double h = 1e-4;
	double H[3][3];
	for (int k = 0; k < 3; k++) {
		d3 a = p, b = p, ga, gb;
		a[k] += h; b[k] -= h;
		vg(a, ga); vg(b, gb);
		for (int i = 0; i < 3; i++) H[i][k] = (ga[i] - gb[i]) / (2.0 * h);
	}
	for (int i = 0; i < 3; i++)
		for (int k = 0; k < i; k++) H[i][k] = H[k][i] = 0.5 * (H[i][k] + H[k][i]);
	//(I - tau H) d = tau g by Cramer; false unless the matrix is positive definite (Sylvester) and d within len
	auto solve = [&](const double tau, d3 &d) {
		double A[3][3];
		for (int i = 0; i < 3; i++)
			for (int k = 0; k < 3; k++) A[i][k] = (i == k) - tau * H[i][k];
		const double m2 = A[0][0] * A[1][1] - A[0][1] * A[1][0];
		const double det = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
		if (!(A[0][0] > 0.0) || !(m2 > 0.0) || !(det > 0.0)) return false;
		const d3 b{ tau * g[0], tau * g[1], tau * g[2] };
		for (int c = 0; c < 3; c++) {
			double M[3][3];
			for (int i = 0; i < 3; i++)
				for (int k = 0; k < 3; k++) M[i][k] = k == c ? b[i] : A[i][k];
			d[c] = (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0])) / det;
		}
		return std::isfinite(d[0] + d[1] + d[2]) && array_length(d) <= len;
	};
	double tau = len / gn, best = -1.0;
	d3 d, dbest{};
	for (int k = 0; k < 60 && solve(tau, d); k++) { best = tau; dbest = d; tau *= 2.0; }
	for (tau = 0.5 * len / gn; best < 0.0 && tau > 1e-12 * len / gn; tau *= 0.5)
		if (solve(tau, d)) { best = tau; dbest = d; }
	for (int k = 0; k < 30 && best > 0.0; k++) {
		const d3 t{ p[0] + dbest[0], p[1] + dbest[1], p[2] + dbest[2] };
		d3 gt;
		const double ft = vg(t, gt);
		if (std::isfinite(ft) && ft > f) { p = t; f = ft; g = gt; return true; }
		best *= 0.5;
		if (!solve(best, dbest)) break;
	}
	return false;
}

std::vector<d4> analytic_eli_maxima(const WFN &wavy, const bool debug, const eli_spin_field *eli)
{
	//rho comes out of computeELIGrad's own orbital pass; the spin field leaves it unset (nullptr)
	auto eli_grad = [&](const d3 &p, double &f, d3 &g, double *rho = nullptr) {
		if (eli) (*eli)(p, f, g, nullptr);
		else if (t_batch) { double r; t_batch->eval(p, f, g, r); if (rho) *rho = r; }
		else wavy.computeELIGrad(p, f, g, rho);
	};
	auto sphere = [](const int n) {
		std::vector<d3> d(n);
		for (int i = 0; i < n; i++) {
			const double z = 1.0 - 2.0 * (i + 0.5) / n;
			const double s = std::sqrt(std::max(0.0, 1.0 - z * z));
			const double phi = 2.39996322972865332 * i;
			d[i] = d3{ s * std::cos(phi), s * std::sin(phi), z };
		}
		return d;
	};
	const std::vector<d3> inner = sphere(26), outer = sphere(110);
	//bohr; the core shells of a fourth-row atom sit inside 0.1, a lone pair or bond maximum 0.6-2.5 out
	const double radii[] = { 0.03, 0.06, 0.1, 0.15, 0.22, 0.3, 0.4, 0.55, 0.75, 1.0, 1.3, 1.7, 2.2, 2.8, 3.5 };
	const int ncen = wavy.get_ncen();
	std::vector<d3> seeds;
	for (int a = 0; a < ncen; a++) {
		const d3 pa = wavy.get_atom_pos(a);
		seeds.push_back(pa);
		for (const double r : radii)
			for (const d3 &u : (r < 0.5 ? inner : outer)) {
				const d3 p{ pa[0] + r * u[0], pa[1] + r * u[1], pa[2] + r * u[2] };
				//Every seed belongs to the atom nearest to it, so no region is seeded twice
				bool own = true;
				for (int b = 0; b < ncen && own; b++)
					if (b != a && array_length(p, wavy.get_atom_pos(b)) < r) own = false;
				if (own) seeds.push_back(p);
			}
	}
	//A spin field's surface attractor can sit behind an interior maximum and a minimum where no seed above reaches, so a
	//spin field also starts one bounded ascent per outer direction on the isosurface itself (first point out from 3.5 bohr
	//in 0.5 bohr steps with rho below the cutoff)
	const size_t n_interior = seeds.size();
	if (eli)
		for (int a = 0; a < ncen; a++) {
			const d3 pa = wavy.get_atom_pos(a);
			for (const d3 &u : outer)
				for (double r = 3.5; r < 25.0; r += 0.5) {
					const d3 p{ pa[0] + r * u[0], pa[1] + r * u[1], pa[2] + r * u[2] };
					if (wavy.compute_dens(p) >= basin_density_cutoff) continue;
					bool own = true;
					for (int b = 0; b < ncen && own; b++)
						if (b != a && array_length(p, wavy.get_atom_pos(b)) < r) own = false;
					if (own) seeds.push_back(p);
					break;
				}
		}
	//Monotone ascent along the normalised gradient: a step is taken only if ELI-D rises, grows after a success and halves
	//after a failure. No Hessian test: a hydrogen's valence maximum sits on the cusp of the nucleus
	auto climb = [&](d3 p) {
		d3 g;
		double f;
		if (wavy.compute_dens(p) < basin_density_cutoff) return d4{ 0.0, 0.0, 0.0, -1.0 };
		eli_grad(p, f, g);
		if (!std::isfinite(f)) return d4{ 0.0, 0.0, 0.0, -1.0 };
		bool inside = true;
		double step = 0.05;
		for (int it = 0; it < 500 && step > 1e-5; it++) {
			const double gn = array_length(g);
			if (!std::isfinite(gn) || gn < 1e-12) break;
			const d3 t{ p[0] + step * g[0] / gn, p[1] + step * g[1] / gn, p[2] + step * g[2] / gn };
			double ft, rt = 0.0;
			d3 gt;
			eli_grad(t, ft, gt, eli ? nullptr : &rt);
			if (std::isfinite(ft) && ft > f) {
				p = t; f = ft; g = gt;
				step = std::min(1.5 * step, 0.3);
				//ELI-D is undefined in the tail and the basins end at the same isosurface; the closed-shell rho is the one
				//computeELIGrad just summed
				if ((eli ? wavy.compute_dens(p) : rt) < basin_density_cutoff) {
					//A spin field climbing out has its attractor on the isosurface; the spin-summed field drops such climbs
					inside = eli != nullptr && bounded_eli_ascent(wavy, *eli, p, f);
					break;
				}
			}
			else step *= 0.5;
		}
		return inside ? d4{ p[0], p[1], p[2], f } : d4{ 0.0, 0.0, 0.0, -1.0 };
	};
	const int ns = static_cast<int>(seeds.size());
	std::vector<d4> ends(ns, d4{ 0.0, 0.0, 0.0, -1.0 });
	//Closed-shell climbs go to the device like the point loops, one seed per fiber; auto samples the seeds inside the isosurface
	std::vector<std::unique_ptr<field_feed>> feeds;
	std::ostringstream why;
	if (eli) why << "the spin-resolved ELI-D";
	else {
		vec sample;
		for (int s = 0; s < ns && g_basin_gpu < 0; s += std::max(1, ns / 1024))
			if (wavy.compute_dens(seeds[s]) >= basin_density_cutoff) sample.insert(sample.end(), seeds[s].begin(), seeds[s].end());
		int threads = 1;
#ifdef _OPENMP
		threads = omp_get_max_threads();
#endif
		feeds = open_field_feeds(wavy, true, sample, [&](const d3 &p) { double f; d3 g; wavy.computeELIGrad(p, f, g); }, threads, why);
	}
	std::cout << "  ELI-D maxima search on the " << (feeds.empty() ? "host" : "device") << ": " << why.str() << std::endl;
	if (!feeds.empty()) {
		std::atomic<int> next{ 0 };
#pragma omp parallel
		{
			int tid = 0;
#ifdef _OPENMP
			tid = omp_get_thread_num();
#endif
			feed_loop(*feeds[tid], next, ns, [&](const int s) { ends[s] = climb(seeds[s]); });
		}
	}
	else {
#pragma omp parallel for schedule(dynamic, 16)
		for (int s = 0; s < ns; s++) {
			if (static_cast<size_t>(s) < n_interior) { ends[s] = climb(seeds[s]); continue; }
			d3 p = seeds[s];
			double f;
			if (bounded_eli_ascent(wavy, *eli, p, f)) ends[s] = d4{ p[0], p[1], p[2], f };
		}
	}
	feeds.clear();
	//Highest first; an end within 0.1 bohr of a kept one is that maximum reached again
	std::vector<int> order;
	for (int s = 0; s < ns; s++) if (ends[s][3] >= 0.0) order.push_back(s);
	std::stable_sort(order.begin(), order.end(), [&](const int x, const int y) { return ends[x][3] > ends[y][3]; });
	std::vector<d4> maxima;
	for (const int s : order) {
		bool seen = false;
		for (const d4 &m : maxima)
			if (std::pow(m[0] - ends[s][0], 2) + std::pow(m[1] - ends[s][1], 2) + std::pow(m[2] - ends[s][2], 2) < 0.01) { seen = true; break; }
		if (!seen) maxima.push_back(ends[s]);
	}
	//The explicit climb stops short on a shell's crest, 0.1 bohr from the true maximum with a gradient left (Ni(CO)4: 4 of 14
	//shell maxima, each stealing up to 0.6 e). Implicit steps finish each walk; then the dedup runs again. A spin field's surface
	//attractor stays where bounded_eli_ascent left it
	{
		auto vg = [&](const d3 &p, d3 &g) { double f; eli_grad(p, f, g); return f; };
#pragma omp parallel for schedule(dynamic, 1)
		for (int k = 0; k < static_cast<int>(maxima.size()); k++) {
			d3 p{ maxima[k][0], maxima[k][1], maxima[k][2] }, g;
			if (wavy.compute_dens(p) < 2.0 * basin_density_cutoff) continue;
			double f = vg(p, g);
			for (int it = 0; it < 200 && implicit_ascent_step(vg, p, f, g, 0.05); it++) {}
			if (wavy.compute_dens(p) >= basin_density_cutoff && f >= maxima[k][3]) maxima[k] = d4{ p[0], p[1], p[2], f };
		}
		std::stable_sort(maxima.begin(), maxima.end(), [](const d4 &x, const d4 &y) { return x[3] > y[3]; });
		std::vector<d4> polished;
		for (const d4 &e : maxima) {
			bool seen = false;
			for (const d4 &m : polished)
				if (std::pow(m[0] - e[0], 2) + std::pow(m[1] - e[1], 2) + std::pow(m[2] - e[2], 2) < 0.01) { seen = true; break; }
			if (!seen) polished.push_back(e);
		}
		if (debug) std::cout << "ELI-D maxima after implicit polish: " << polished.size() << " of " << maxima.size() << std::endl;
		maxima.swap(polished);
	}
	//A seed on a mirror plane climbs inside it and can stop on a saddle of the full field. As in the density's accept, six
	//axis-displaced starts must come back; a saddle loses one to a clearly higher end elsewhere. Not tested: a nuclear
	//maximum, a spin field's surface attractor and the near-equal maxima on a core shell's ridge, which the core/shell
	//unification merges anyway
	//ponytail: a saddle within 0.5 % of the maximum it falls to is kept
	std::vector<char> saddle(maxima.size(), 0);
#pragma omp parallel for schedule(dynamic, 1)
	for (int k = 0; k < static_cast<int>(maxima.size()); k++) {
		const d3 m{ maxima[k][0], maxima[k][1], maxima[k][2] };
		double rn = std::numeric_limits<double>::max();
		for (int a = 0; a < ncen; a++) rn = std::min(rn, array_length(m, wavy.get_atom_pos(a)));
		const double delta = std::min(0.05, 0.3 * rn);
		if (delta < 1e-3 || wavy.compute_dens(m) < 2.0 * basin_density_cutoff) continue;
		for (int t = 0; t < 6 && !saddle[k]; t++) {
			d3 p = m;
			p[t / 2] += t % 2 ? -delta : delta;
			const d4 e = climb(p);
			saddle[k] = e[3] > 1.005 * maxima[k][3] && std::pow(e[0] - m[0], 2) + std::pow(e[1] - m[1], 2) + std::pow(e[2] - m[2], 2) >= 0.01;
		}
	}
	size_t kept = 0;
	for (size_t k = 0; k < maxima.size(); k++) {
		if (!saddle[k]) { maxima[kept++] = maxima[k]; continue; }
		if (debug) std::cout << "Dropped an ELI-D saddle found as a maximum at " << maxima[k][0] << " " << maxima[k][1] << " " << maxima[k][2] << " (ELI-D " << maxima[k][3] << ")" << std::endl;
	}
	maxima.resize(kept);
	if (debug)
		std::cout << "ELI-D maxima from " << ns << " seeds by analytic ascent: " << order.size() << " climbs finished, " << maxima.size() << " distinct maxima" << std::endl;
	return maxima;
}

static bool g_beta_spheres = true;
void beta_spheres_set_enabled(const bool on) { g_beta_spheres = on; }
bool beta_spheres_enabled() { return g_beta_spheres; }
//The adaptive step uses field turning to grow on straight trajectories.
static constexpr double adp_cap_default = 8.0;    //at most this many times the validated floor step
static constexpr double adp_grow_default = 0.99999;   //midpoint cosine that earns a doubling
static constexpr double adp_keep_default = 0.999;   //below this the step is thrown away and retaken at the floor
static constexpr double adp_reach_default = 0.25;  //fraction of the distance to the nearest maximum
static double g_adp_cap = adp_cap_default;
static double g_adp_grow = adp_grow_default;
static double g_adp_keep = adp_keep_default;
static double g_adp_reach = adp_reach_default;
static double tuning_double(const char *name, const double fallback, const double upper = 0.0)
{
	const char *v = tuning(name);
	if (v == nullptr || *v == '\0') return fallback;
	try {
		const double d = std::stod(v);
		if (d > 0.0 && std::isfinite(d) && (upper <= 0.0 || d <= upper)) return d;
		std::cout << "Ignoring " << name << "=" << v << ": not a positive finite number";
		if (upper > 0.0) std::cout << " of at most " << upper;
		std::cout << std::endl;
	}
	catch (const std::exception &) {
		std::cout << "Ignoring " << name << "=" << v << ": not a number" << std::endl;
	}
	return fallback;
}
//Auto (-1) times host_eval on the sample (xyz) on one thread and the device on it in one batch after a warm-up; the device
//is taken when a host evaluation costs at least -tune NOS_BASIN_GPU_MIN_US and the device outruns all host threads together.
static std::vector<std::unique_ptr<field_feed>> open_field_feeds(const WFN &wavy, const bool eli, const vec &sample,
	const std::function<void(const d3 &)> &host_eval, const int threads, std::ostringstream &why)
{
	std::vector<std::unique_ptr<field_feed>> feeds;
	if (g_basin_gpu == 0) { why << "-no_basin_gpu"; return feeds; }
	const int fibers = std::max(1, static_cast<int>(tuning_double("NOS_BASIN_GPU_FIBERS", 256.0)));
	const int ns = static_cast<int>(sample.size() / 3);
	auto at = [&](const int i) { return d3{ sample[3 * i], sample[3 * i + 1], sample[3 * i + 2] }; };
	double host_us = 0.0;
	if (g_basin_gpu < 0 && ns > 1) {
		host_eval(at(0)); //sizes the scratch and is not counted
		const auto t0 = std::chrono::steady_clock::now();
		for (int i = 1; i < ns; i++) host_eval(at(i));
		host_us = 1e6 * std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() / (ns - 1);
	}
	const double min_us = tuning_double("NOS_BASIN_GPU_MIN_US", basin_gpu_min_us);
	std::ostringstream us;
	us << std::setprecision(3) << "the host field costs " << host_us << " us a point";
	if (g_basin_gpu < 0 && host_us < min_us) { why << "auto, " << us.str() << ", below " << min_us; return feeds; }
	for (int t = 0; t < threads; t++) {
		feeds.push_back(std::make_unique<field_feed>(wavy, eli, fibers));
		if (!feeds.back()->ctx) { feeds.clear(); why << "no usable device"; return feeds; }
	}
	if (g_basin_gpu > 0) { why << "-basin_gpu"; return feeds; }
	double dev_us = 1e30;
	void *ctx = wavy.field_gpu_open(eli, ns);
	vec v(ns), g(3 * static_cast<size_t>(ns)), r(ns);
	auto run = [&]() { return WFN::field_gpu_run(ctx, ns, sample.data(), eli ? v.data() : nullptr, g.data(), r.data()); };
	if (ctx && run()) {
		const auto t0 = std::chrono::steady_clock::now();
		if (run()) dev_us = 1e6 * std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() / ns;
	}
	WFN::field_gpu_close(ctx);
	us << " on " << threads << " threads, the device " << dev_us;
	if (host_us < threads * dev_us) { feeds.clear(); why << "auto, " << us.str() << " (slower than the host)"; }
	else why << "auto, " << us.str() << ", at least " << min_us;
	return feeds;
}
static void adp_knobs_from_tuning()
{
	//From the defaults every time, so removing a knob restores its default
	g_adp_cap = tuning_double("NOS_ADP_CAP", adp_cap_default);
	g_adp_grow = tuning_double("NOS_ADP_GROW", adp_grow_default);
	g_adp_keep = tuning_double("NOS_ADP_KEEP", adp_keep_default);
	g_adp_reach = tuning_double("NOS_ADP_REACH", adp_reach_default);
}
//Keep beta spheres within the smallest sampled safe radius.
static constexpr double beta_margin_default = 0.9;
double basin_beta_margin() { return tuning_double("NOS_BETA_MARGIN", beta_margin_default, 1.0); }
void basin_adaptive_step_knobs(double &cap, double &grow, double &keep, double &reach)
{
	cap = g_adp_cap; grow = g_adp_grow; keep = g_adp_keep; reach = g_adp_reach;
}
static bool g_adaptive_step = false;
static std::atomic<long long> g_adp_steps{ 0 }, g_adp_tries{ 0 }, g_adp_turn{ 0 }, g_adp_fall{ 0 }, g_adp_shrink{ 0 };
//Trajectories that ran their step budget out
static std::atomic<long long> g_adp_exhaust{ 0 };
//ELI-D stalls handed to the implicit ridge walk, and how many of those reached an attractor
static std::atomic<long long> g_ridge_walks{ 0 }, g_ridge_done{ 0 };
//Steps a single trajectory may take before it is given up on (NOS_BASIN_STEP_CAP)
static constexpr int step_cap_default = 2000;
static int g_step_cap = step_cap_default;
//Factor by which a shortened step grows back toward the base step
static constexpr double step_relax_default = 1.5;
static double g_step_relax = step_relax_default;
//Whether a floor step that fails to rise may halve and try again (NOS_BASIN_SHRINK=0 turns it off)
static bool g_step_shrink = true;
//Fraction of the base step below which a failing step stops halving (NOS_BASIN_SHRINK_FLOOR).
//ponytail: 1/64 clears more stalls but splits near-degenerate lone pairs apart, so the default stays 1/16
static constexpr double shrink_floor_default = 0.0625;
static double g_shrink_floor = shrink_floor_default;
//Account separately for density that never reaches a basin.
static std::atomic<long long> g_stall_vacuum{ 0 }, g_stall_field{ 0 }, g_stall_far{ 0 };
static std::atomic<double> g_stall_far_rho{ 0.0 };
//Bucket stalled gradients to distinguish tail and interior failures.
static constexpr int g_stall_basins = 8;
static std::atomic<long long> g_stall_gn[4] = {};
//The >=1e-2 bucket again, cleared when an integration starts so basin_stalls_on_a_slope can read it after
static std::atomic<long long> g_stall_slope{ 0 };
static std::atomic<long long> g_stall_to[g_stall_basins + 1] = {};
//A negative distance marks a stall below the density cutoff.
static void basin_stall_seen(const double dist, const double rho, const double gn = -1.0)
{
	if (dist < 0.0) { g_stall_vacuum.fetch_add(1, std::memory_order_relaxed); return; }
	g_stall_field.fetch_add(1, std::memory_order_relaxed);
	if (gn >= 0.0) {
		const int b = gn < 1e-6 ? 0 : gn < 1e-4 ? 1 : gn < 1e-2 ? 2 : 3;
		g_stall_gn[b].fetch_add(1, std::memory_order_relaxed);
		if (b == 3) g_stall_slope.fetch_add(1, std::memory_order_relaxed);
	}
	if (dist <= 1.0) return;
	g_stall_far.fetch_add(1, std::memory_order_relaxed);
	double cur = g_stall_far_rho.load(std::memory_order_relaxed);
	while (rho > cur && !g_stall_far_rho.compare_exchange_weak(cur, rho, std::memory_order_relaxed)) {}
}
long long basin_stalls_on_a_slope()
{
	return g_stall_slope.load();
}
static void basin_stall_gave_to(const int basin)
{
	g_stall_to[basin >= 0 && basin < g_stall_basins ? basin : g_stall_basins].fetch_add(1, std::memory_order_relaxed);
}
static void basin_stall_counters(long long &vacuum, long long &in_field, long long &beyond_a_bohr, double &worst_rho)
{
	vacuum = g_stall_vacuum.load();
	in_field = g_stall_field.load();
	beyond_a_bohr = g_stall_far.load();
	worst_rho = g_stall_far_rho.load();
}
static void basin_stall_counters_reset()
{
	g_stall_vacuum = 0; g_stall_field = 0; g_stall_far = 0; g_stall_far_rho = 0.0;
	for (auto &c : g_stall_gn) c = 0;
	for (auto &c : g_stall_to) c = 0;
}
static inline void adp_count(std::atomic<long long> &c) { c.fetch_add(1, std::memory_order_relaxed); }
void basin_adaptive_step_counters(long long &steps, long long &proposed, long long &turned_back, long long &fell_back)
{
	steps = g_adp_steps.load();
	proposed = g_adp_tries.load();
	turned_back = g_adp_turn.load();
	fell_back = g_adp_fall.load();
}
void basin_adaptive_step_counters_reset()
{
	g_adp_steps = 0; g_adp_tries = 0; g_adp_turn = 0; g_adp_fall = 0; g_adp_shrink = 0; g_adp_exhaust = 0;
}
void basin_adaptive_step_set_enabled(const bool on) { g_adaptive_step = on; if (on) adp_knobs_from_tuning(); }
bool basin_adaptive_step_enabled() { return g_adaptive_step; }
static double g_basin_step_scale = 1.0;
void basin_step_scale_set(const double f) { g_basin_step_scale = f > 0.0 ? f : 1.0; }
double basin_step_scale() { return g_basin_step_scale; }
static bool g_basin_timing = false;
void basin_timing_set_enabled(const bool on) { g_basin_timing = on; }
bool basin_timing_enabled() { return g_basin_timing; }
void basin_flags_reset()
{
	g_basin_gpu = -1; g_beta_spheres = true; g_adaptive_step = false; g_basin_step_scale = 1.0; g_basin_timing = false;
}
void basin_stage_timer::lap(const std::string &what) {
	const auto now = std::chrono::steady_clock::now();
	const double s = std::chrono::duration<double>(now - t).count();
	t = now;
	if (g_basin_timing) std::cout << "  [timing] " << what << ": " << std::fixed << std::setprecision(2) << s << " s" << std::endl;
}

//Sized for the atoms' box plus 4 bohr each side at half load, so a molecule smaller than the cap does not pay for it
basin_memo::basin_memo(const d3 &lo, const d3 &hi, const double voxel, const size_t megabytes) : inv(1.0 / voxel)
{
	double n = 2.0;
	for (int d = 0; d < 3; d++) {
		origin[d] = 0.5 * (lo[d] + hi[d]) - (1 << (bits - 1)) * voxel;
		n *= (hi[d] - lo[d] + 8.0) / voxel;
	}
	const size_t cap = std::max<size_t>(megabytes << 20, 1 << 19) / sizeof(uint64_t);
	size_t size = 1 << 16;
	while (size * 2 <= cap && static_cast<double>(size) < n) size <<= 1;
	table = std::vector<std::atomic<uint64_t>>(size);
	mask = size - 1;
	//Linear probing degrades past three quarters full; inserts stop there and reads go on
	//ponytail: a full table just stops learning; grow-and-rehash if a large system hits the cap
	limit = size / 4 * 3;
}

uint64_t basin_memo::key(const d3 &p) const
{
	uint64_t k = 0;
	for (int d = 0; d < 3; d++) {
		const double t = std::floor((p[d] - origin[d]) * inv);
		//one voxel of margin both sides, so every neighbour of a key is a key
		if (!(t >= 1.0 && t <= static_cast<double>((1 << bits) - 2))) return 0;
		k |= static_cast<uint64_t>(t) << (d * bits);
	}
	return k;
}

static inline uint64_t memo_hash(uint64_t k)
{
	k ^= k >> 33; k *= 0xff51afd7ed558ccdull; k ^= k >> 33; k *= 0xc4ceb9fe1a85ec53ull; k ^= k >> 33;
	return k;
}

uint64_t basin_memo::find(const uint64_t k) const
{
	for (size_t i = memo_hash(k) & mask, n = 0; n <= mask; n++, i = (i + 1) & mask) {
		const uint64_t e = table[i].load(std::memory_order_relaxed);
		if (e == 0) return 0;
		if ((e & key_mask) == k) return e;
	}
	return 0;
}

void basin_memo::add(const uint64_t k, const uint64_t label)
{
	for (size_t i = memo_hash(k) & mask, n = 0; n <= mask; n++, i = (i + 1) & mask) {
		uint64_t e = table[i].load(std::memory_order_relaxed);
		if (e == 0) {
			if (used.load(std::memory_order_relaxed) >= limit) return;
			if (table[i].compare_exchange_strong(e, k | label << 42 | one, std::memory_order_relaxed)) {
				used.fetch_add(1, std::memory_order_relaxed);
				return;
			}
			//e is now whoever won the slot, which may be this very voxel
		}
		if ((e & key_mask) != k) continue;
		for (;;) {
			if (e & conflict) return;
			uint64_t next = e;
			if (((e >> 42) & label_mask) != label) next |= conflict;
			else if ((e & ~conflict) >> 55 < 255) next += one;
			else return;
			if (table[i].compare_exchange_weak(e, next, std::memory_order_relaxed)) return;
		}
	}
}

void basin_memo::write(const std::vector<uint64_t> &path, const int label)
{
	if (label <= 0 || static_cast<uint64_t>(label) > label_mask) return;
	for (const uint64_t k : path) add(k, static_cast<uint64_t>(label));
}

int basin_memo::settled(const d3 &p) const
{
	const uint64_t k = key(p);
	if (k == 0) return 0;
	const uint64_t c = find(k);
	if (c == 0 || (c & conflict) || (c >> 55) < 2) return 0;
	const uint64_t label = (c >> 42) & label_mask;
	for (int dz = -1; dz <= 1; dz++)
		for (int dy = -1; dy <= 1; dy++)
			for (int dx = -1; dx <= 1; dx++) {
				if (dx == 0 && dy == 0 && dz == 0) continue;
				//every index stays inside its 14 bits, so the offsets add field by field
				const uint64_t e = find(k + static_cast<uint64_t>(dx + dy * (int64_t(1) << bits) + dz * (int64_t(1) << 2 * bits)));
				if (e == 0 || (e & conflict) || ((e >> 42) & label_mask) != label) return 0;
			}
	return static_cast<int>(label);
}

vec integrate_basins_on_atomic_grids(const cube *cub, const cubei *basin_cube, const std::vector<d4> &maxima, const WFN &wavy, const int accuracy, const bool eli_field, vec &volumes, double &outside, const std::function<double(const d3&)> *core_density, const std::function<void(const d3&, d3&)> *core_gradient, const int grid_boost, const density_field *field, basin_overlaps *ovl, const ivec *maximum_basin, const eli_spin_field *eli, vec2 *spin_pop, const std::vector<shell_split> *split, const std::function<void(const std::function<int(const d3&)>&)> *probe, const bool probe_only, std::vector<basin_point> *points)
{
	if (!eli_field) eli = nullptr;
	if (!eli) spin_pop = nullptr;
	//The filled core steers the trajectories only. An ECP atom's grid is built for its
	//valence basis and cannot integrate a 1s at Z = 80, so the core electrons are added to
	//the nucleus's basin by count once the valence density is integrated; a Thakkar core
	//lies whole inside its atom's basin
	auto valence = [&](const d3 &p) -> double {
		if (t_batch) { double v, r; d3 g; t_batch->eval(p, v, g, r); return r; }
		return field ? field->rho(p) : wavy.compute_dens(p);
	};
	//Streaming: no cube and no basin cube; the maxima are the whole topology and every point walks the field to its basin
	const bool streaming = cub == nullptr || basin_cube == nullptr;
	//Cleared here and not at the print, so the count belongs to this call and survives it
	g_stall_slope = 0;
	//From the defaults every time, so a knob an earlier run in this process set does not stick
	g_step_cap = step_cap_default; g_step_relax = step_relax_default; g_step_shrink = true; g_shrink_floor = shrink_floor_default;
	if (const char *e = tuning("NOS_BASIN_STEP_CAP")) {
		const int v = std::atoi(e);
		if (v >= 100) g_step_cap = v;
	}
	if (const char *e = tuning("NOS_BASIN_STEP_RELAX")) {
		const double v = std::atof(e);
		if (v >= 1.0 && v <= 4.0) g_step_relax = v;
	}
	if (const char *e = tuning("NOS_BASIN_SHRINK")) g_step_shrink = std::string(e) != "0";
	if (const char *e = tuning("NOS_BASIN_SHRINK_FLOOR")) {
		const double v = std::atof(e);
		if (v > 0.0 && v <= 1.0) g_shrink_floor = v;
	}
	//The basin a maximum belongs to, 1-based; without a map that is the maximum's own index
	auto final_of = [&](const size_t m) { return maximum_basin ? (*maximum_basin)[m + 1] : static_cast<int>(m) + 1; };
	int nb = 0;
	if (!streaming) nb = basin_cube->max_value();
	else if (!maximum_basin) nb = static_cast<int>(maxima.size());
	else for (size_t m = 0; m < maxima.size(); m++) nb = std::max(nb, final_of(m));
	//A split shoulder (shell_split) answers walks as basin nb + 1 + k; settle hands the point to the exit whose saddle is
	//nearest in angle about the atom. Exits are found by position: every later merge keeps a member's exact position
	struct split_at { d3 o; std::vector<std::pair<d3, int>> exits; int home; };
	std::vector<split_at> splits;
	ivec split_of(maxima.size(), -1);
	if (split && streaming && eli_field) {
		auto raw = [&](const d3 &x) {
			for (size_t m = 0; m < maxima.size(); m++)
				if (std::pow(x[0] - maxima[m][0], 2) + std::pow(x[1] - maxima[m][1], 2) + std::pow(x[2] - maxima[m][2], 2) < 1e-12) return static_cast<int>(m);
			return -1;
		};
		for (const shell_split &s : *split) {
			const int m = raw(s.at);
			if (m < 0 || split_of[m] >= 0) continue;
			split_at t{ s.o, {}, final_of(m) };
			bool several = false;
			for (const auto &e : s.exits) {
				const int c = raw(e.second);
				if (c < 0) continue;
				t.exits.emplace_back(e.first, final_of(c));
				several |= final_of(c) != t.exits.front().second;
			}
			if (!several) continue;
			split_of[m] = static_cast<int>(splits.size());
			splits.push_back(std::move(t));
		}
	}
	auto basin_of = [&](const size_t m) { return split_of[m] >= 0 ? nb + 1 + split_of[m] : final_of(m); };
	auto settle = [&](const int n, const d3 &p) {
		if (n <= nb) return n;
		const split_at &s = splits[n - nb - 1];
		const d3 v{ p[0] - s.o[0], p[1] - s.o[1], p[2] - s.o[2] };
		if (!(array_length(v) > 1e-6)) return s.home;
		int best = s.home;
		double cos_best = -2.0;
		for (const auto &e : s.exits) {
			const double c = v[0] * e.first[0] + v[1] * e.first[1] + v[2] * e.first[2];
			if (c > cos_best) { cos_best = c; best = e.second; }
		}
		return best;
	};
	//The overlap matrices use the populations' points and weights: a point's density is sum_i occ_i phi_i^2, so their
	//diagonal sums to the population exactly
	if (ovl && field) ovl = nullptr;
	if (ovl) {
		ovl->mo_index.clear();
		for (int m = 0; m < wavy.get_nmo(); m++)
			if (std::abs(wavy.get_MO_occ(m)) > 1e-8) ovl->mo_index.push_back(m);
		ovl->nmo = static_cast<int>(ovl->mo_index.size());
		ovl->S.assign(nb, vec(ovl->triangle(), 0.0));
		ovl->outside.assign(ovl->triangle(), 0.0);
	}
	vec pop(nb, 0.0);
	volumes.assign(nb, 0.0);
	if (spin_pop) spin_pop->assign(nb, vec(3, 0.0));
	if (points) points->clear();
	outside = 0.0;
	double cutoff_outside = 0.0, unresolved_outside = 0.0;
	const int nx = streaming ? 0 : cub->get_size(0), ny = streaming ? 0 : cub->get_size(1), nz = streaming ? 0 : cub->get_size(2);
	//Without a cube the trajectory keeps the step of a 0.1 A grid
	d3 h{ constants::ang2bohr(0.1), constants::ang2bohr(0.1), constants::ang2bohr(0.1) };
	if (cub)
		for (int i = 0; i < 3; i++)
			h[i] = std::sqrt(cub->get_vector(0, i) * cub->get_vector(0, i) + cub->get_vector(1, i) * cub->get_vector(1, i) + cub->get_vector(2, i) * cub->get_vector(2, i));
	//A third of a voxel within 1.5 bohr of a nucleus, a whole voxel beyond
	const double voxel = std::min({ h[0], h[1], h[2] });
	const std::vector<atom> atoms = wavy.get_atoms();
	const double sscale = basin_step_scale();
	auto step_at = [&](const d3 &p) {
		double d2 = std::numeric_limits<double>::max();
		for (const atom &at : atoms) {
			const d3 ap = at.get_pos();
			d2 = std::min(d2, std::pow(p[0] - ap[0], 2) + std::pow(p[1] - ap[1], 2) + std::pow(p[2] - ap[2], 2));
		}
		if (d2 < 2.25) return sscale * 0.3 * voxel;
		//The density cutoff bounds the tail outside the basins.
		if (streaming && d2 > 36.0) return sscale * std::min(0.25 * std::sqrt(d2), 4.0);
		return sscale * voxel;
	};
	const double step = 0.3 * voxel;
	//Streaming only: a gridded run's basin cube already answers what the memo would. Half a step voxel, so a 0.19 bohr
	//step crosses about two memo voxels and the midpoint fills the gap.
	//ponytail: 13-bit labels, so more than 8191 basins runs without the memo; widen the label field if that happens
	std::unique_ptr<basin_memo> memo;
	{
		const char *e = tuning("NOS_BASIN_MEMO");
		if (streaming && !atoms.empty() && nb < (1 << 13) && !(e != nullptr && std::string(e) == "0")) {
			d3 lo = atoms[0].get_pos(), hi = lo;
			for (const atom &at : atoms)
				for (int d = 0; d < 3; d++) { lo[d] = std::min(lo[d], at.get_pos()[d]); hi[d] = std::max(hi[d], at.get_pos()[d]); }
			memo = std::make_unique<basin_memo>(lo, hi, 0.5 * voxel, static_cast<size_t>(tuning_double("NOS_BASIN_MEMO_MB", 256.0)));
		}
	}
	//Trajectories by who asked for them (outer probe, inner probe, cell centre, bisection), their step iterations and how many a settled memo block ended
	enum { climb_outer, climb_inner, climb_centre, climb_bisect };
	std::atomic<long long> climbs[4]{}, climb_steps[4]{}, memo_hits{ 0 };
	//Cube cell of a position and the position within it; false outside the cube
	auto cell = [&](const d3 &p, int *c, d3 &f) {
		if (streaming) return false;
		const int sz[3] = { nx, ny, nz };
		for (int d = 0; d < 3; d++) {
			const double t = (p[d] - cub->get_origin(d)) / cub->get_vector(d, d);
			c[d] = static_cast<int>(std::floor(t));
			if (c[d] < 0 || c[d] >= sz[d] - 1) return false;
			f[d] = t - c[d];
		}
		return true;
	};
	//Basin of the nearest grid point, and whether every voxel within band of the cell agrees
	const int band = 1;
	auto lookup = [&](const d3 &p, bool &settled) {
		int c[3]; d3 f;
		settled = false;
		if (!cell(p, c, f)) return 0;
		const int b = basin_cube->get_value(c[0] + (f[0] > 0.5), c[1] + (f[1] > 0.5), c[2] + (f[2] > 0.5));
		settled = b != 0;
		for (int dx = -band; dx < band + 2 && settled; dx++)
			for (int dy = -band; dy < band + 2 && settled; dy++)
				for (int dz = -band; dz < band + 2; dz++) {
					const int x = c[0] + dx, y = c[1] + dy, z = c[2] + dz;
					if (x < 0 || y < 0 || z < 0 || x >= nx || y >= ny || z >= nz) continue;
					if (basin_cube->get_value(x, y, z) != b) { settled = false; break; }
				}
		return b;
	};
	//The basin whose maximum lies within reach of p, 0 when none does. Two voxels serve an
	//ELI-D maximum, which is broad; a nucleus gets a tenth of a bohr, since a hydroxyl
	//hydrogen's basin is 0.4 bohr thick and a wider net catches the oxygen's electrons
	const double catch2 = eli_field ? std::pow(2.0 * std::max({ h[0], h[1], h[2] }), 2) : 0.01;
	//Store each attractor centre and squared beta-sphere radius; -1 marks a core maximum that
	//answers through its nucleus's sphere and has no catch radius of its own.
	vec beta2(maxima.size(), 0.0);
	std::vector<d3> bcen(maxima.size());
	for (size_t m = 0; m < maxima.size(); m++) bcen[m] = d3{ maxima[m][0], maxima[m][1], maxima[m][2] };
	auto at_maximum = [&](const d3 &p) {
		for (size_t m = 0; m < maxima.size(); m++) {
			//a measured beta sphere replaces the catch radius, which can reach past a core shell minimum
			if (beta2[m] > 0.0) {
				if (std::pow(p[0] - bcen[m][0], 2) + std::pow(p[1] - bcen[m][1], 2) + std::pow(p[2] - bcen[m][2], 2) < beta2[m])
					return basin_of(m);
			}
			else if (beta2[m] == 0.0 && std::pow(p[0] - maxima[m][0], 2) + std::pow(p[1] - maxima[m][1], 2) + std::pow(p[2] - maxima[m][2], 2) < catch2)
				return basin_of(m);
		}
		return 0;
	};
	//The basin of the maximum nearest p, 0 when further than reach: a stalled trajectory has run out of slope and the point
	//still has to belong to somebody
	auto nearest_maximum = [&](const d3 &p, const double reach) {
		int best = 0;
		double d2 = reach * reach;
		for (size_t m = 0; m < maxima.size(); m++) {
			const double q = std::pow(p[0] - maxima[m][0], 2) + std::pow(p[1] - maxima[m][1], 2) + std::pow(p[2] - maxima[m][2], 2);
			if (q < d2) { d2 = q; best = basin_of(m); }
		}
		return best;
	};
	//Half the distance to the nearest maximum: a grown step must still land inside the catch radius (or beta sphere) of the
	//attractor it walks into, since at_maximum only looks at where the walk landed
	auto reach_limit = [&](const d3 &p) {
		double d2 = std::numeric_limits<double>::max();
		for (size_t m = 0; m < maxima.size(); m++)
			d2 = std::min(d2, std::pow(p[0] - maxima[m][0], 2) + std::pow(p[1] - maxima[m][1], 2) + std::pow(p[2] - maxima[m][2], 2));
		return g_adp_reach * std::sqrt(d2);
	};
	//Evaluate density and gradient in the same orbital pass.
	auto gradient = [&](const d3 &p, d3 &g, double *val = nullptr) {
		if (!eli_field) {
			if (field) { field->grad(p, g); if (val) *val = field->rho(p); }
			else if (t_batch) { double v, r; t_batch->eval(p, v, g, r); if (val) *val = r; }
			else wavy.computeGrad(p, g, val);
			if (core_gradient) { d3 c; (*core_gradient)(p, c); for (int k = 0; k < 3; k++) g[k] += c[k]; }
			if (val && core_density) *val += (*core_density)(p);
			return;
		}
		double e;
		if (eli) (*eli)(p, e, g, nullptr);
		else if (t_batch) { double r; t_batch->eval(p, e, g, r); }
		else wavy.computeELIGrad(p, e, g);
		if (val) *val = e;
	};
	//The field's value and gradient in one call, so the climb sees when it stops rising; computeELIGrad builds both from one orbital pass
	auto value_and_gradient = [&](const d3 &p, d3 &g) {
		double v = 0.0;
		gradient(p, g, &v);
		return v;
	};
	//A beta sphere terminates trajectories that cannot leave its attractor basin.
	basin_stage_timer T;
	const std::string fieldname = eli_field ? "ELI-D " : "QTAIM ";
	double margin = 0.0;
	if (streaming && beta_spheres_enabled() && !maxima.empty()) {
		constexpr int ndir = 302;
		static const std::vector<d3> dirs = [] {
			std::vector<d3> d(ndir);
			for (int i = 0; i < ndir; i++) {
				const double z = 1.0 - 2.0 * (i + 0.5) / ndir;
				const double s = std::sqrt(std::max(0.0, 1.0 - z * z));
				const double phi = 2.39996322972865332 * i;   //golden angle: no two samples line up
				d[i] = d3{ s * std::cos(phi), s * std::sin(phi), z };
			}
			return d;
		}();
		const double march = 0.05;   //bohr; the radius is only ever needed to within a step
		const int nm = static_cast<int>(maxima.size());
		//Once, outside the parallel loop: every sphere has to be drawn at the same margin
		margin = basin_beta_margin();
#pragma omp parallel for schedule(dynamic)
		for (int m = 0; m < nm; m++) {
			//Ascend onto the attractor first: a cube maximum is a voxel centre, and past the true top
			//the radial derivative already points back in
			d3 c{ maxima[m][0], maxima[m][1], maxima[m][2] };
			{
				d3 g;
				double f = value_and_gradient(c, g), sl = 0.5 * voxel;
				for (int it = 0; it < 60 && sl > 1e-4; it++) {
					const double gn = std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
					if (gn < 1e-12) break;
					const d3 t{ c[0] + sl * g[0] / gn, c[1] + sl * g[1] / gn, c[2] + sl * g[2] / gn };
					d3 gt;
					const double ft = value_and_gradient(t, gt);
					if (ft > f) { c = t; f = ft; g = gt; }
					else sl *= 0.5;
				}
			}
			//Walked further than a voxel and a half: that is a different attractor and a sphere would answer for the wrong basin,
			//so the maximum keeps the plain catch radius and no sphere
			if (std::pow(c[0] - maxima[m][0], 2) + std::pow(c[1] - maxima[m][1], 2) + std::pow(c[2] - maxima[m][2], 2) > std::pow(1.5 * voxel, 2))
				continue;
			double cap = 3.0;
			for (int n = 0; n < nm; n++) {
				if (n == m || basin_of(n) == basin_of(m)) continue;
				const double d = std::sqrt(std::pow(c[0] - maxima[n][0], 2) + std::pow(c[1] - maxima[n][1], 2) + std::pow(c[2] - maxima[n][2], 2));
				cap = std::min(cap, 0.45 * d);
			}
			//A maximum on a nucleus heads the atom's unified core: the shells out to its furthest core
			//maximum are its own basin, so a ray crosses their minima and stops at the first beyond
			double inner = 0.0;
			for (const atom &a : atoms) {
				const d3 ap = a.get_pos();
				if (std::pow(c[0] - ap[0], 2) + std::pow(c[1] - ap[1], 2) + std::pow(c[2] - ap[2], 2) > march * march) continue;
				for (int n = 0; n < nm; n++)
					if (n != m && basin_of(n) == basin_of(m))
						inner = std::max(inner, std::sqrt(std::pow(c[0] - maxima[n][0], 2) + std::pow(c[1] - maxima[n][1], 2) + std::pow(c[2] - maxima[n][2], 2)));
			}
			double r = cap;
			for (const d3 &u : dirs) {
				double rr = inner + march;
				bool fell = inner == 0.0;
				for (; rr <= cap + 1e-12; rr += march) {
					const d3 q{ c[0] + rr * u[0], c[1] + rr * u[1], c[2] + rr * u[2] };
					d3 g;
					gradient(q, g);
					if (g[0] * u[0] + g[1] * u[1] + g[2] * u[2] < 0.0) fell = true;
					else if (fell) break;
				}
				//never past the outer core shell's top within the cap: nothing measured
				if (!fell) { r = 0.0; break; }
				r = std::min(r, rr - march);
				if (r <= march) break;
			}
			if (r <= 2.0 * march) continue;
			bcen[m] = c;
			beta2[m] = std::pow(margin * r, 2);
		}
		//A core maximum inside its nucleus's sphere answers through that sphere; its own catch radius
		//reaches past the shell minimum the sphere stops at
		for (int m = 0; m < nm; m++) {
			if (beta2[m] > 0.0) continue;
			for (int n = 0; n < nm; n++)
				if (beta2[n] > 0.0 && basin_of(n) == basin_of(m) && std::pow(maxima[m][0] - bcen[n][0], 2) + std::pow(maxima[m][1] - bcen[n][1], 2) + std::pow(maxima[m][2] - bcen[n][2], 2) < beta2[n]) { beta2[m] = -1.0; break; }
		}
	}
	T.lap(fieldname + "beta spheres");
	if (basin_timing_enabled() && margin > 0.0) {
		//A maximum inside a same-basin sphere (-1) is covered by it; on a heavy atom that is most of them
		size_t with = 0, via = 0;
		for (const double b : beta2) { if (b > 0.0) with++; else if (b < 0.0) via++; }
		std::cout << "  [timing] " << fieldname << "beta spheres: margin " << margin << ", "
			<< with << " of " << beta2.size() << " maxima carry one, " << via << " more answer through a same-basin sphere, "
			<< beta2.size() - with - via << " have only the catch radius" << std::endl;
	}
	//At least level 3; strongly ionic centres can need level 4
	GridConfiguration config;
	config.accuracy = std::max(accuracy, 3);
	config.alpha_max_scale = static_cast<double>(grid_boost) * grid_boost;
	config.radial_step_scale = grid_boost;
	config.angular_boost = grid_boost - 1;
	config.partition_type = PartitionType::Becke;
	config.no_density_eval = true;
	const char *rot_env = tuning("NOS_BASIN_ROTATE");
	config.rotate_angular = rot_env == nullptr || std::string(rot_env) != "0";
	GridManager grids(config);
	ivec every_atom(wavy.get_ncen());
	std::iota(every_atom.begin(), every_atom.end(), 0);
	grids.setup3DGridsForMolecule(wavy, every_atom);
	const GridData &gd = grids.getGridData();
	T.lap(fieldname + "atomic quadrature grids");
	int host_threads = 1;
#ifdef _OPENMP
	host_threads = omp_get_max_threads();
#endif
	//On the device each host thread runs -tune NOS_BASIN_GPU_FIBERS trajectories as fibers and batches their field calls
	//on its own feed. Plain wavefunction fields only: a density_source field and the spin ELI-D stay on the host.
	//Auto decides on a sample of the grid points inside the density isosurface (open_field_feeds).
	std::vector<std::unique_ptr<field_feed>> feeds;
	std::ostringstream why;
	if (!streaming) why << "the basin cube decides";
	else if (field) why << "a fitted density source";
	else if (eli) why << "the spin-resolved ELI-D";
	else {
		vec sample;
		if (g_basin_gpu < 0) {
			size_t total = 0;
			for (const auto &a : gd.atomic_grids) total += a[GridData::X].size();
			const size_t stride = std::max<size_t>(1, total / 1024);
			for (const auto &a : gd.atomic_grids)
				for (size_t i = 0; i < a[GridData::X].size(); i += stride) {
					const d3 p{ a[GridData::X][i], a[GridData::Y][i], a[GridData::Z][i] };
					if (a[GridData::BECKE_WEIGHT][i] != 0.0 && valence(p) >= basin_density_cutoff) sample.insert(sample.end(), p.begin(), p.end());
				}
		}
		feeds = open_field_feeds(wavy, eli_field, sample, [&](const d3 &p) { d3 g; double v; gradient(p, g, &v); }, host_threads, why);
		T.lap(fieldname + "device probe and feeds");
	}
	std::cout << "  " << fieldname << "basin field on the " << (feeds.empty() ? "host" : "device") << ": " << why.str() << std::endl;
	//body(i) for the calling thread's share of [0, np): its fibers on its feed, or the host's dynamic schedule
	std::atomic<int> next_point{ 0 };
	auto point_loop = [&](const int np, const auto &body) {
		if (feeds.empty()) {
#pragma omp for schedule(dynamic, 16) nowait
			for (int i = 0; i < np; i++) body(i);
			return;
		}
		int tid = 0;
#ifdef _OPENMP
		tid = omp_get_thread_num();
#endif
		feed_loop(*feeds[tid], next_point, np, body);
	};
	//No reach limit: the density isosurface bounds gridded and streaming basins alike
	const double stall_reach = 1e30;
	//The trajectory and quadrature share the density isosurface.
	const double stall_floor = basin_density_cutoff;
	//An ELI-D walk stalls on a flat ridge such as a lone-pair torus, where the nearest maximum can be
	//the core behind the shell minimum. Take the maximum, among those up to twice the nearest
	//distance away, whose straight path from p stays highest. d2 is the nearest squared.
	auto highest_path = [&](const d3 &p, const double d2) {
		int best = 0;
		double best_low = -1.0, best_q = 0.0;
		for (size_t m = 0; m < maxima.size(); m++) {
			const d3 v{ maxima[m][0] - p[0], maxima[m][1] - p[1], maxima[m][2] - p[2] };
			const double q = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
			if (q > 4.0 * d2) continue;
			double low = std::numeric_limits<double>::max();
			for (int k = 1; k < 8; k++)
			{
				const d3 x{ p[0] + k / 8.0 * v[0], p[1] + k / 8.0 * v[1], p[2] + k / 8.0 * v[2] };
				double e;
				d3 g;
				if (eli) (*eli)(x, e, g, nullptr);
				else e = wavy.computeELI(x);
				low = std::min(low, e);
			}
			if (low > best_low || (low == best_low && q < best_q)) { best = basin_of(m); best_low = low; best_q = q; }
		}
		return best;
	};
	//Resolve a stalled density trajectory from its starting side of the separatrix.
	auto stalled =[&](const d3 &start, const d3 &r, const double gn = -1.0) {
		const double rho = valence(r);
		if (rho < stall_floor) { basin_stall_seen(-1.0, 0.0); return 0; }
		double d2 = std::numeric_limits<double>::max();
		for (size_t m = 0; m < maxima.size(); m++)
			d2 = std::min(d2, std::pow(r[0] - maxima[m][0], 2) + std::pow(r[1] - maxima[m][1], 2) + std::pow(r[2] - maxima[m][2], 2));
		basin_stall_seen(std::sqrt(d2), rho, gn);
		const int b = eli_field ? highest_path(r, d2) : nearest_maximum(start, stall_reach);
		basin_stall_gave_to(b);
		return b;
	};
	//An ELI-D walk that stops rising on a steep slope straddles a ridge crest, not a maximum (Ni(CO)4: 4623 of 4726 stalls at
	//|g| >= 1e-2, all on the Ni shell). Implicit steps follow the crest to its attractor; 0 when they stop rising first
	auto ridge_walk = [&](d3 r, const double len) {
		if (!eli_field) return 0;
		adp_count(g_ridge_walks);
		d3 g;
		double f = value_and_gradient(r, g);
		for (int s = 0; s < 400; s++) {
			if (const int m = at_maximum(r)) { adp_count(g_ridge_done); return m; }
			if (valence(r) < basin_density_cutoff || !implicit_ascent_step(value_and_gradient, r, f, g, len)) return 0;
		}
		return 0;
	};
	//kind says who asked (climb_outer...). Every clean ending is written to the memo, but only
	//the edge probes read it: a cell centre or a bisection probe is climbed because it sits within
	//a cell of the separatrix, which is exactly where a 3x3x3 block can be wrong
	auto climb_walk = [&](const d3 &p, long long &lb, long long &ll, const int kind) {
		if (valence(p) < basin_density_cutoff) return 0;
		bool settled;
		int b = lookup(p, settled);
		//A gridded ELI-D takes a settled cell straight from the cube and leaves the crop outside;
		//streaming has no cube to ask and walks from every point
		if (eli_field && !streaming && (settled || b == 0)) return b;
		int c[3]; d3 f;
		//Off the cube there is nothing to integrate; streaming has no cube to be off
		if (!streaming && !cell(p, c, f)) return 0;
		//Already inside a beta sphere (or on a maximum): the answer needs no trajectory, so it is
		//not counted as one either
		if (const int m0 = at_maximum(p)) return m0;
		lb++;
		climbs[kind].fetch_add(1, std::memory_order_relaxed);
		struct tally { std::atomic<long long> &to; long long n = 0; ~tally() { to.fetch_add(n, std::memory_order_relaxed); } } steps{ climb_steps[kind] };
		const bool read_memo = memo && (kind == climb_outer || kind == climb_inner);
		//The voxels this walk has held, and the midpoint of the step in flight, which joins them
		//only once the point it led to has been accepted
		std::vector<uint64_t> path;
		uint64_t pending = 0;
		auto note = [&](const uint64_t k) { if (k != 0 && (path.empty() || path.back() != k)) path.push_back(k); };
		//Field at r_prev, kept so a step that is taken back does not evaluate it a second time
		d3 g_prev{};
		double here_prev = 0.0;
		bool have_prev = false, reuse = false;
		d3 r = p, g;
		double last_value = -1.0;
		//Grow only steps whose field direction stays nearly straight.
		double mult = 1.0;
		//Shorten a base step if it would decrease the field.
		double shrink = 1.0;
		//Whether the step that reached r was grown: mult > 1 at the top of an iteration only says the next step may be,
		//and only a grown last step is grounds for throwing a point away
		bool grown_last = false;
		d3 r_prev = p;
		const bool grow = basin_adaptive_step_enabled();
		double value_prev = -1.0;
		for (int s = 0; s < g_step_cap; s++) {
			steps.n++;
			//Two steps of its own first, so a probe that starts on a boundary does not read the
			//block it starts in
			if (read_memo && s >= 2)
				if (const int mb = memo->settled(r)) { memo->write(path, mb); memo_hits.fetch_add(1, std::memory_order_relaxed); return mb; }
			//A streaming ELI-D step evaluates the field at r anyway; that orbital pass carries the
			//density too, so the isosurface test costs no separate density evaluation
			bool have_vg = false;
			double vg = 0.0;
			if (reuse) {
				//r is r_prev again, whose field (and density test) the accepted step already has
				g = g_prev; vg = here_prev; have_vg = true; reuse = false;
			}
			else if (eli_field) {
				double rho_r;
				if (streaming && eli) { double aux[4]; (*eli)(r, vg, g, aux); rho_r = aux[0]; have_vg = true; }
				else if (streaming && !field) {
					if (t_batch) t_batch->eval(r, vg, g, rho_r);
					else wavy.computeELIGrad(r, vg, g, &rho_r);
					have_vg = true;
				}
				else rho_r = valence(r);
				if (rho_r < basin_density_cutoff) {
					//A spin field's walk that leaves the domain goes on along the isosurface to the attractor analytic_eli_maxima found there
					if (!(streaming && eli)) return 0;
					d3 q = r;
					double fq;
					if (!bounded_eli_ascent(wavy, *eli, q, fq)) return 0;
					if (const int mb = at_maximum(q)) return mb;
					//A near-degenerate surface stops the ascent anywhere on it; the nearest maximum of the same value within
					//unify_boundary_basins' 3 bohr linkage is its basin, and that pass then merges the shell
					int best = 0;
					double d2 = 9.0;
					for (size_t m = 0; m < maxima.size(); m++) {
						if (std::abs(maxima[m][3] - fq) > 1e-3 * std::max(maxima[m][3], fq)) continue;
						const double dq = std::pow(q[0] - maxima[m][0], 2) + std::pow(q[1] - maxima[m][1], 2) + std::pow(q[2] - maxima[m][2], 2);
						if (dq < d2) { d2 = dq; best = basin_of(m); }
					}
					//ponytail: a few exits stall just short of an axial surface attractor and end here with 0; a constrained Newton step if that ever carries weight
					return best ? best : nearest_maximum(q, 0.1);
				}
			}
			if (grow) adp_count(g_adp_steps);
			const int m = at_maximum(r);
			if (m) {
				if (memo) { note(pending); note(memo->key(r)); memo->write(path, m); }
				return m;
			}
			const double floor_step = step_at(r);
			//The gridded ELI-D climb is steered by the cube below and never asked whether it is
			//still rising; every streaming walk is, since nothing else can stop it
			double here = last_value;
			if (!eli_field || streaming) {
				here = have_vg ? vg : value_and_gradient(r, g);
				if (here <= last_value) {
					//Retry a grown step at the base length after excessive turning.
					if (grown_last) { adp_count(g_adp_fall); r = r_prev; last_value = value_prev; mult = 1.0; grown_last = false; reuse = have_prev; pending = 0; continue; }
					mult = 1.0;
					if (g_step_shrink && shrink > g_shrink_floor) {
						adp_count(g_adp_shrink);
						shrink *= 0.5;
						r = r_prev;
						last_value = value_prev;
						grown_last = false;
						reuse = have_prev;
						pending = 0;
						continue;
					}
					int n = ridge_walk(r_prev, floor_step);
					if (!n) n = stalled(p, r, std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]));
					if (n) return n;
					break;
				}
			}
			else gradient(r, g);
			if (memo) { note(pending); pending = 0; note(memo->key(r)); }
			double gn = std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
			if (gn < 1e-12) {
				//A critical point: a streaming density walk has no cube behind it, so it goes to the
				//stall tie-break instead of outside every basin
				if (!eli_field) { const int n = stalled(p, r, gn); if (n) return n; }
				break;
			}
			//Preserve the midpoint expression order used for validated basin populations.
			const d3 g0{ g[0], g[1], g[2] };
			const double gn0 = gn;
			const d3 dir{ g0[0] / gn0, g0[1] / gn0, g0[2] / gn0 };
			//Retry a step at the base length when its field direction turns too far.
			const double base = floor_step * shrink;
			double sl = base;
			double cosine = 0.0;
			bool stepped = false;
			d3 mid;
			for (int attempt = 0; attempt < 2 && !stepped; attempt++) {
				sl = mult > 1.0 ? std::min(base * mult, std::max(base, reach_limit(r))) : base;
				if (mult > 1.0) adp_count(g_adp_tries);
				for (int k = 0; k < 3; k++) mid[k] = r[k] + 0.5 * sl * g0[k] / gn0;
				gradient(mid, g);
				gn = std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
				if (gn < 1e-12) break;
				cosine = (g[0] * dir[0] + g[1] * dir[1] + g[2] * dir[2]) / gn;
				if (mult > 1.0 && cosine < g_adp_keep) { adp_count(g_adp_turn); mult = 1.0; continue; }
				stepped = true;
			}
			//Same case one level in: the attempt loop gives up only when the midpoint gradient vanished too
			if (!stepped) {
				if (!eli_field) { const int n = stalled(p, r, gn); if (n) return n; }
				break;
			}
			r_prev = r;
			g_prev = g0;
			here_prev = here;
			have_prev = true;
			if (memo) pending = memo->key(mid);
			value_prev = last_value;
			last_value = here;
			grown_last = sl > base;
			for (int k = 0; k < 3; k++) r[k] += sl * g[k] / gn;
			if (shrink < 1.0) shrink = std::min(1.0, shrink * g_step_relax);
			else if (cosine > g_adp_grow && mult < g_adp_cap && grow) mult *= 2.0;
			if (!streaming) {
				const int b2 = lookup(r, settled);
				if (b2 == 0 && eli_field) { ll++; break; }
				if (b2) b = b2;
				if (settled && eli_field) break;
			}
		}
		//Assign a stalled streaming trajectory by its starting point.
		if (streaming) {
			adp_count(g_adp_exhaust);
			if (const int n = ridge_walk(r, step_at(r))) return n;
			return stalled(p, r);
		}
		return b;
	};
	auto climb = [&](const d3 &p, long long &lb, long long &ll, const int kind) { return settle(climb_walk(p, lb, ll, kind), p); };
	//Only the caller's points, no quadrature: the memo then fills from the caller's own climbs
	if (probe && probe_only) {
		(*probe)([&](const d3 &p) { long long lb = 0, ll = 0; return climb(p, lb, ll, climb_outer); });
		return pop;
	}
	//Bisect a radial cell only when its edge trajectories reach different basins.
	constexpr int bisections = 6;
	long long boundary_points = 0, lost = 0;
	//NOS_BASIN_SKIP: electrons per molecule the lightest points may move instead of being climbed; 0 climbs every point.
	//The default, half the last printed digit of a population, cannot move a basin by a printed unit
	const char *skip_env = tuning("NOS_BASIN_SKIP");
	const double skip_budget = (skip_env != nullptr && std::string(skip_env) == "0") ? 0.0 : tuning_double("NOS_BASIN_SKIP", 5e-5);
	long long skipped_points = 0;
	double skipped_mass = 0.0;
	for (size_t a = 0; a < gd.atomic_grids.size(); a++) {
		const vec &X = gd.atomic_grids[a][GridData::X], &Y = gd.atomic_grids[a][GridData::Y], &Z = gd.atomic_grids[a][GridData::Z], &W = gd.atomic_grids[a][GridData::BECKE_WEIGHT];
		const int np = static_cast<int>(X.size());
		const d3 centre = atoms[a].get_pos();
		vec radius(np), shells;
		for (int i = 0; i < np; i++) radius[i] = std::sqrt(std::pow(X[i] - centre[0], 2) + std::pow(Y[i] - centre[1], 2) + std::pow(Z[i] - centre[2], 2));
		shells = radius;
		std::sort(shells.begin(), shells.end());
		shells.erase(std::unique(shells.begin(), shells.end(), [](double x, double y) { return std::abs(x - y) < 1e-8; }), shells.end());
		//The cell of point i as a shell segment: the two radii halfway to the neighbouring
		//shells. False for a cell with no radial extent, which is given to one basin whole.
		auto cell_edges = [&](const int i, double &in, double &out) {
			const size_t k = std::lower_bound(shells.begin(), shells.end(), radius[i] - 1e-8) - shells.begin();
			in = k > 0 ? 0.5 * (shells[k - 1] + shells[k]) : 0.0;
			out = k + 1 < shells.size() ? 0.5 * (shells[k] + shells[k + 1]) : shells[k];
			return radius[i] > 1e-8 && out > in + 1e-8;
		};
		//A point at radius r along point i's direction
		auto along_i = [&](const int i, const double r) {
			const double f = r / radius[i];
			return d3{ centre[0] + (X[i] - centre[0]) * f, centre[1] + (Y[i] - centre[1]) * f, centre[2] + (Z[i] - centre[2]) * f };
		};
		//Reuse probes only for adjacent shells with matching directions.
		ivec partner(np, -1);
		{
			ivec start;
			for (int i = 0; i < np; i++)
				if (i == 0 || std::abs(radius[i] - radius[i - 1]) > 1e-8) start.push_back(i);
			start.push_back(np);
			for (size_t s = 1; s + 1 < start.size(); s++) {
				const int a0 = start[s - 1], a1 = start[s], a2 = start[s + 1];
				if (a1 - a0 != a2 - a1) continue;
				//Only a cell whose inner edge is the previous cell's outer edge, which needs the shells to be neighbours in shells[];
				//a pair closer than the 1e-8 of the unique pass shares one entry and has no edge between them
				double i0, o0, i1, o1;
				if (!cell_edges(a0, i0, o0) || !cell_edges(a1, i1, o1) || o0 != i1) continue;
				for (int j = 0; j < a1 - a0; j++) {
					const int lo = a0 + j, hi = a1 + j;
					const double dot = ((X[lo] - centre[0]) * (X[hi] - centre[0]) + (Y[lo] - centre[1]) * (Y[hi] - centre[1]) + (Z[lo] - centre[2]) * (Z[hi] - centre[2])) / (radius[lo] * radius[hi]);
					if (dot > 1.0 - 1e-12) partner[hi] = lo;
				}
			}
		}
		//The points carrying the least |w| * rho, summed from the smallest up to this grid's share of the budget, are not
		//climbed but go to the nearest maximum, moving at most the budget in electrons and budget / cutoff in volume; that is
		//the Becke weight's long tail. |w|, since the rule has negative weights that would sort first and pull the sum down.
		vec rho_at;
		std::vector<char> skip;
		if (skip_budget > 0.0) {
			rho_at.assign(np, 0.0);
			bool on_device = false;
			if (!feeds.empty()) {
				ivec at;
				vec q, gq, rq;
				for (int i = 0; i < np; i++)
					if (W[i] != 0.0) { at.push_back(i); q.insert(q.end(), { X[i], Y[i], Z[i] }); }
				gq.resize(3 * at.size());
				rq.resize(at.size());
				on_device = at.empty() || WFN::field_gpu_run(feeds[0]->ctx,static_cast<int>(at.size()), q.data(), nullptr, gq.data(), rq.data());
				if (on_device)
					for (size_t k = 0; k < at.size(); k++) rho_at[at[k]] = rq[k];
			}
			if (!on_device) {
#pragma omp parallel for schedule(dynamic, 64)
				for (int i = 0; i < np; i++)
					if (W[i] != 0.0) rho_at[i] = valence(d3{ X[i], Y[i], Z[i] });
			}
			ivec order;
			for (int i = 0; i < np; i++)
				if (W[i] != 0.0 && rho_at[i] >= basin_density_cutoff) order.push_back(i);
			std::sort(order.begin(), order.end(), [&](const int x, const int y) { return std::abs(W[x]) * rho_at[x] < std::abs(W[y]) * rho_at[y]; });
			skip.assign(np, 0);
			double sum = 0.0;
			const double share = skip_budget / static_cast<double>(gd.atomic_grids.size());
			for (const int i : order) {
				if (sum + std::abs(W[i]) * rho_at[i] > share) break;
				sum += std::abs(W[i]) * rho_at[i];
				skip[i] = 1;
				skipped_points++;
			}
			skipped_mass += sum;
		}
		//A skipped cell's outer probe is still wanted when a climbed cell reads it as its inner edge
		std::vector<char> probe_wanted(np, 1);
		if (!skip.empty()) {
			for (int i = 0; i < np; i++) probe_wanted[i] = !skip[i];
			for (int i = 0; i < np; i++)
				if (partner[i] >= 0 && !skip[i]) probe_wanted[partner[i]] = 1;
		}
		//Every cell's outer probe, climbed once. A gridded ELI-D takes most cells from the cube without probing, so it keeps
		//the on-demand path and this pass is skipped.
		ivec outer_probe(np, -1);
		if (!(eli_field && !streaming)) {
			next_point = 0;
#pragma omp parallel num_threads(host_threads)
			{
				long long lb = 0, ll = 0;
				point_loop(np, [&](const int i) {
					double in, out;
					if (W[i] == 0.0 || !probe_wanted[i] || !cell_edges(i, in, out)) return;
					outer_probe[i] = climb(along_i(i, out), lb, ll, climb_outer);
				});
#pragma omp critical
				{ boundary_points += lb; lost += ll; }
			}
		}
		next_point = 0;
#pragma omp parallel num_threads(host_threads)
		{
			vec lp(nb, 0.0), lv(nb, 0.0), lspin(spin_pop ? 3 * (size_t)nb : 0, 0.0);
			std::vector<basin_point> lpts;
			double lo = 0.0, lc = 0.0, lu = 0.0;
			long long lb = 0, ll = 0;
			//ponytail: one triangle per basin per thread, nb * nmo^2 / 2 doubles each; the caller
			//sizes the job, a molecule big enough to hurt here has other limits first
			vec2 ls;
			vec ls_out;
			vec phi;
			vec2 dbuf;
			int phi_at = -1; //the point phi belongs to: the thread's fibers share it across their climbs
			if (ovl) {
				ls.assign(nb, vec(ovl->triangle(), 0.0));
				ls_out.assign(ovl->triangle(), 0.0);
				phi.resize(wavy.get_nmo(), 0.0);
				dbuf.assign(wavy.get_ncen(), vec(16, 0.0));
			}
			//Rank-1 update of one basin's triangle: the point's share of the quadrature weight
			//times the outer product of the orbitals it sees
			auto accumulate = [&](const int b, const double wq) {
				if (!ovl || wq == 0.0) return;
				double *Sb = b > 0 ? ls[b - 1].data() : ls_out.data();
				const int *idx = ovl->mo_index.data();
				for (int a2 = 0; a2 < ovl->nmo; a2++) {
					const double pa = wq * phi[idx[a2]];
					if (pa == 0.0) continue;
					double *row = Sb + (size_t)a2 * (a2 + 1) / 2;
					for (int b2 = 0; b2 <= a2; b2++) row[b2] += pa * phi[idx[b2]];
				}
			};
			point_loop(np, [&](const int i) {
				const double w = W[i];
				if (w == 0.0) return;
				const d3 p{ X[i], Y[i], Z[i] };
				d3 rec_at = p;
				bool settled;
				int b = lookup(p, settled);
				//The same density, taken from the orbital pass that also hands out phi
				const double rho = ovl ? (phi_at = i, wavy.compute_dens(p, dbuf, phi)) : (rho_at.empty() ? valence(p) : rho_at[i]);
				//rho, rho_alpha, rho_beta, rho_s Y_q at the point; only basins carry them, so not below the cutoff
				double aux[4]{ 0.0, 0.0, 0.0, 0.0 };
				if (spin_pop && rho >= basin_density_cutoff) { double e; d3 g; (*eli)(p, e, g, aux); }
				//A basin's share of the cell's quadrature weight; the weight itself stays with
				//the rule, only who gets it is decided here
				auto give = [&](const int bb, const double fr) {
					if (fr <= 0.0) return;
					if (ovl && phi_at != i) { wavy.compute_dens(p, dbuf, phi); phi_at = i; }
					if (bb == 0) {
						lo += w * rho * fr;
						if (rho < basin_density_cutoff) lc += w * rho * fr; else lu += w * rho * fr;
						accumulate(0, w * fr);
						return;
					}
					lp[bb - 1] += w * rho * fr;
					lv[bb - 1] += w * fr;
					if (spin_pop) for (int k = 0; k < 3; k++) lspin[3 * (size_t)(bb - 1) + k] += w * fr * aux[1 + k];
					if (points) lpts.push_back({ rec_at, w * rho * fr, bb });
					accumulate(bb, w * fr);
				};
				if (rho < basin_density_cutoff) { give(0, 1.0); return; }
				if (!skip.empty() && skip[i]) { give(settle(nearest_maximum(p, stall_reach), p), 1.0); return; }
				//A gridded ELI-D takes a cell whose neighbourhood agrees from the grid; streaming has no grid and refines every cell
				if (eli_field && !streaming && (b == 0 || settled)) { give(b, 1.0); return; }
				//Climb the cell centre only when the edge probes do not decide its basin.
				int bc = -1;
				auto centre_basin = [&]() { if (bc < 0) bc = climb(p, lb, ll, climb_centre); return bc; };
				double inner, outer;
				if (!cell_edges(i, inner, outer)) { give(centre_basin(), 1.0); return; }
				auto along = [&](const double r) { return along_i(i, r); };
				//A failed edge probe defers to the cell centre.
				const int pin = partner[i];
				int bi = pin >= 0 && outer_probe[pin] >= 0 ? outer_probe[pin] : climb(along(inner), lb, ll, climb_inner);
				int bo = outer_probe[i] >= 0 ? outer_probe[i] : climb(along(outer), lb, ll, climb_outer);
				if (bi == 0) bi = centre_basin();
				if (bo == 0) bo = centre_basin();
				if (bi == bo) { give(bi, 1.0); return; }
				//ponytail: one crossing per cell. Three basins meeting inside a single quadrature
				//cell is a smaller thing than the rule's own error; bisect for more if it is not
				double lo_r = inner, hi_r = outer;
				for (int it = 0; it < bisections; it++) {
					const double mid = 0.5 * (lo_r + hi_r);
					int bm = climb(along(mid), lb, ll, climb_bisect);
					if (bm == 0) bm = centre_basin();
					if (bm == bi) lo_r = mid; else hi_r = mid;
				}
				const double rc = 0.5 * (lo_r + hi_r);
				//Each side gets the density it carries over its own part of the shell segment,
				//whose volume goes as r^3
				const double wi = valence(along(0.5 * (inner + rc))) * (rc * rc * rc - inner * inner * inner);
				const double wo = valence(along(0.5 * (rc + outer))) * (outer * outer * outer - rc * rc * rc);
				const double sum = wi + wo;
				if (sum > 0.0) {
					rec_at = along(0.5 * (inner + rc)); give(bi, wi / sum);
					rec_at = along(0.5 * (rc + outer)); give(bo, wo / sum);
				}
				else give(centre_basin(), 1.0);
			});
#pragma omp critical
			{
				for (int b = 0; b < nb; b++) { pop[b] += lp[b]; volumes[b] += lv[b]; }
				if (spin_pop) for (int b = 0; b < nb; b++) for (int k = 0; k < 3; k++) (*spin_pop)[b][k] += lspin[3 * (size_t)b + k];
				if (points) points->insert(points->end(), lpts.begin(), lpts.end());
				if (ovl) {
					for (int b = 0; b < nb; b++)
						for (size_t t = 0; t < ovl->S[b].size(); t++) ovl->S[b][t] += ls[b][t];
					for (size_t t = 0; t < ovl->outside.size(); t++) ovl->outside[t] += ls_out[t];
				}
				outside += lo;
				cutoff_outside += lc;
				unresolved_outside += lu;
				boundary_points += lb;
				lost += ll;
			}
		}
	}
	if (core_density)
		for (int a = 0; a < wavy.get_ncen(); a++) {
			const int ncore = wavy.get_atom_ECP_electrons(a);
			if (ncore <= 0) continue;
			bool settled;
			//Streaming has no basin cube to read the nucleus out of; it is one of the maxima by
			//construction, so the maximum it sits on is its basin
			const int b = streaming ? settle(at_maximum(wavy.get_atom_pos(a)), wavy.get_atom_pos(a)) : lookup(wavy.get_atom_pos(a), settled);
			if (b > 0) pop[b - 1] += ncore;
			else outside += ncore;
		}
	T.lap(fieldname + "point loop");
	if (g_basin_timing && !feeds.empty()) {
		long long points = 0, batches = 0;
		double device_s = 0.0;
		bool failed = false;
		for (const auto &f : feeds) { points += f->points; batches += f->batches; device_s += f->device_s; failed |= f->failed; }
		std::cout << "  [timing] " << fieldname << "GPU field: " << points << " points in " << batches << " batches from "
			<< host_threads << " host threads x " << feeds[0]->fibers << " fibers, " << std::fixed << std::setprecision(2) << device_s
			<< " thread-s waiting on the device" << (failed ? ", device failed part way, the host finished" : "") << std::endl;
	}
	{
		//Read and reset unconditionally so the counts never carry into the next field or molecule
		long long sv, sf, sfar; double srho;
		basin_stall_counters(sv, sf, sfar, srho);
		long long sgn[4], sto[g_stall_basins + 1];
		for (int i = 0; i < 4; i++) sgn[i] = g_stall_gn[i].load();
		for (int i = 0; i <= g_stall_basins; i++) sto[i] = g_stall_to[i].load();
		const long long sx = g_adp_exhaust.exchange(0);
		const long long rw = g_ridge_walks.exchange(0), rd = g_ridge_done.exchange(0);
		basin_stall_counters_reset();
		if (g_basin_timing && rw > 0) std::cout << "  " << fieldname << "ridge walks: " << rd << " of " << rw << " stopped walks reached an attractor by implicit steps" << std::endl;
		//Own stream, so the scientific format does not stay on std::cout
		if (g_basin_timing && sf + sv + sx > 0) {
			std::ostringstream line;
			line << std::scientific << std::setprecision(3);
			line << "  " << fieldname << "trajectories that stalled: " << sf << " in the field, "
				<< sv << " below " << stall_floor << " e/bohr^3; " << sfar
				<< " further than a bohr from every attractor, the densest at " << srho << " e/bohr^3";
			line << std::defaultfloat << "\n  " << fieldname << "stall gradients: ";
			const char *edge[4] = { "<1e-6", "<1e-4", "<1e-2", ">=1e-2" };
			for (int i = 0; i < 4; i++) line << (i ? ", " : "") << sgn[i] << " " << edge[i];
			line << "; " << sx << " ran the step budget out; handed to basin";
			for (int i = 0; i <= g_stall_basins; i++)
				if (sto[i]) line << " " << (i == g_stall_basins ? std::string("other") : std::to_string(i)) << ":" << sto[i];
			std::cout << line.str() << std::endl;
		}
	}
	if (g_basin_timing && skip_budget > 0.0) {
		const auto flags = std::cout.flags();
		std::cout << "  [timing] " << fieldname << "not climbed: " << skipped_points << " points carrying " << std::scientific << skipped_mass << " e, to their nearest maximum" << std::endl;
		std::cout.flags(flags);
	}
	if (g_basin_timing) std::cout << "  [timing] " << fieldname << "outside the density isosurface: " << cutoff_outside << " e, unresolved inside: " << unresolved_outside << " e" << std::endl;
	if (g_basin_timing) {
		const char *who[4] = { "outer", "inner", "centre", "bisection" };
		long long all = 0;
		std::cout << "  [timing] " << fieldname << "climbs (step iterations):";
		for (int k = 0; k < 4; k++) { std::cout << (k ? ", " : " ") << climbs[k] << " " << who[k] << " (" << climb_steps[k] << ")"; all += climb_steps[k]; }
		std::cout << "; " << all << " step iterations";
		if (memo) std::cout << "; memo: " << memo_hits << " ended in a settled block, " << memo->entries() << " voxels in "
			<< memo->capacity() << " slots";
		else std::cout << "; memo off";
		std::cout << std::endl;
	}
	if (g_basin_timing && g_adaptive_step) {
		const long long st = g_adp_steps.exchange(0), tr = g_adp_tries.exchange(0);
		const long long tu = g_adp_turn.exchange(0), fa = g_adp_fall.exchange(0);
		const long long sh = g_adp_shrink.exchange(0);
		std::cout << "  [timing] " << fieldname << "grown step: " << st << " steps, " << tr
			<< " proposals, " << tu << " turned back, " << fa << " fell back, " << sh << " shrunk below the floor, "
			<< std::fixed << std::setprecision(1)
			<< (tr ? 100.0 * static_cast<double>(tu + fa) / static_cast<double>(tr) : 0.0)
			<< " % of proposals wasted" << std::endl;
		std::cout << "  [timing] " << fieldname << "grown step knobs: cap " << std::setprecision(4) << g_adp_cap
			<< ", grow " << std::setprecision(8) << g_adp_grow << ", keep " << g_adp_keep
			<< ", reach " << std::setprecision(4) << g_adp_reach
			<< ", exp cutoff " << constants::exp_cutoff << std::endl;
	}
	std::cout << "Quadrature points sent along the field: " << boundary_points << ", left the grid: " << lost << std::endl;
	//The caller's own points walk the same trajectories to the same basins, after the memo is filled
	if (probe) (*probe)([&](const d3 &p) { long long lb = 0, ll = 0; return climb(p, lb, ll, climb_outer); });
	return pop;
}
