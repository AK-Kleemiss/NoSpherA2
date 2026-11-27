#include "pch.h"
#include "scattering_factors.h"
#include "AtomGrid.h"
#include "fourier_bessel.h"
#ifdef NOSPHERA2_USE_GPU
#include "sf_gpu.h"
#include "itensor_gpu.h"
#endif

//calc_SF phase recurrence: reflections per task, points per L1 tile, multiplies between sincos anchors, max row-jump steps
constexpr int SF_BLOCK = 64, SF_TILE = 256, SF_ANCHOR = 32, SF_NMAX = 4;
//calc_SF lanes: NEON on every ARM64 (Apple included), AVX or SSE2 on x86, plain doubles elsewhere (armv7). sf_fma and
//sf_fms are fused where the ISA has it: c + a b and c - a b.
#if defined(__aarch64__) || defined(_M_ARM64)
using sf_v = float64x2_t;
constexpr int SF_W = 2;
static inline sf_v sf_ld(const double* p) { return vld1q_f64(p); }
static inline void sf_st(double* p, const sf_v a) { vst1q_f64(p, a); }
static inline sf_v sf_set(const double a) { return vdupq_n_f64(a); }
static inline sf_v sf_add(const sf_v a, const sf_v b) { return vaddq_f64(a, b); }
static inline sf_v sf_mul(const sf_v a, const sf_v b) { return vmulq_f64(a, b); }
static inline sf_v sf_fma(const sf_v c, const sf_v a, const sf_v b) { return vfmaq_f64(c, a, b); }
static inline sf_v sf_fms(const sf_v c, const sf_v a, const sf_v b) { return vfmsq_f64(c, a, b); }
static inline double sf_hsum(const sf_v a) { return vgetq_lane_f64(a, 0) + vgetq_lane_f64(a, 1); }
#elif defined(__AVX__)
using sf_v = __m256d;
constexpr int SF_W = 4;
static inline sf_v sf_ld(const double* p) { return _mm256_loadu_pd(p); }
static inline void sf_st(double* p, const sf_v a) { _mm256_storeu_pd(p, a); }
static inline sf_v sf_set(const double a) { return _mm256_set1_pd(a); }
static inline sf_v sf_add(const sf_v a, const sf_v b) { return _mm256_add_pd(a, b); }
static inline sf_v sf_mul(const sf_v a, const sf_v b) { return _mm256_mul_pd(a, b); }
#if defined(__FMA__) || defined(__AVX2__)
static inline sf_v sf_fma(const sf_v c, const sf_v a, const sf_v b) { return _mm256_fmadd_pd(a, b, c); }
static inline sf_v sf_fms(const sf_v c, const sf_v a, const sf_v b) { return _mm256_fnmadd_pd(a, b, c); }
#else
static inline sf_v sf_fma(const sf_v c, const sf_v a, const sf_v b) { return _mm256_add_pd(c, _mm256_mul_pd(a, b)); }
static inline sf_v sf_fms(const sf_v c, const sf_v a, const sf_v b) { return _mm256_sub_pd(c, _mm256_mul_pd(a, b)); }
#endif
static inline double sf_hsum(const sf_v a) { alignas(32) double l[4]; _mm256_store_pd(l, a); return ((l[0] + l[1]) + l[2]) + l[3]; }
static inline sf_v sf_sub(const sf_v a, const sf_v b) { return _mm256_sub_pd(a, b); }
static inline sf_v sf_and(const sf_v a, const sf_v b) { return _mm256_and_pd(a, b); }
static inline sf_v sf_andnot(const sf_v a, const sf_v b) { return _mm256_andnot_pd(a, b); }
static inline sf_v sf_or(const sf_v a, const sf_v b) { return _mm256_or_pd(a, b); }
static inline sf_v sf_xor(const sf_v a, const sf_v b) { return _mm256_xor_pd(a, b); }
static inline bool sf_small(const sf_v a) { return _mm256_movemask_pd(_mm256_cmp_pd(_mm256_andnot_pd(_mm256_set1_pd(-0.0), a), _mm256_set1_pd(1e5), _CMP_LT_OQ)) == 15; }
static inline __m128i sf_round_i(const sf_v a) { return _mm256_cvtpd_epi32(a); }
static inline sf_v sf_from_i(const __m128i q) { return _mm256_cvtepi32_pd(q); }
static inline sf_v sf_wide(const __m128i m) { return _mm256_castsi256_pd(_mm256_insertf128_si256(_mm256_castsi128_si256(_mm_unpacklo_epi32(m, m)), _mm_unpackhi_epi32(m, m), 1)); }
#elif defined(__SSE2__) || defined(_M_X64)
#include <emmintrin.h>
using sf_v = __m128d;
constexpr int SF_W = 2;
static inline sf_v sf_ld(const double* p) { return _mm_loadu_pd(p); }
static inline void sf_st(double* p, const sf_v a) { _mm_storeu_pd(p, a); }
static inline sf_v sf_set(const double a) { return _mm_set1_pd(a); }
static inline sf_v sf_add(const sf_v a, const sf_v b) { return _mm_add_pd(a, b); }
static inline sf_v sf_mul(const sf_v a, const sf_v b) { return _mm_mul_pd(a, b); }
static inline sf_v sf_fma(const sf_v c, const sf_v a, const sf_v b) { return _mm_add_pd(c, _mm_mul_pd(a, b)); }
static inline sf_v sf_fms(const sf_v c, const sf_v a, const sf_v b) { return _mm_sub_pd(c, _mm_mul_pd(a, b)); }
static inline double sf_hsum(const sf_v a) { alignas(16) double l[2]; _mm_store_pd(l, a); return l[0] + l[1]; }
static inline sf_v sf_sub(const sf_v a, const sf_v b) { return _mm_sub_pd(a, b); }
static inline sf_v sf_and(const sf_v a, const sf_v b) { return _mm_and_pd(a, b); }
static inline sf_v sf_andnot(const sf_v a, const sf_v b) { return _mm_andnot_pd(a, b); }
static inline sf_v sf_or(const sf_v a, const sf_v b) { return _mm_or_pd(a, b); }
static inline sf_v sf_xor(const sf_v a, const sf_v b) { return _mm_xor_pd(a, b); }
static inline bool sf_small(const sf_v a) { return _mm_movemask_pd(_mm_cmplt_pd(_mm_andnot_pd(_mm_set1_pd(-0.0), a), _mm_set1_pd(1e5))) == 3; }
static inline __m128i sf_round_i(const sf_v a) { return _mm_cvtpd_epi32(a); }
static inline sf_v sf_from_i(const __m128i q) { return _mm_cvtepi32_pd(q); }
static inline sf_v sf_wide(const __m128i m) { return _mm_castsi128_pd(_mm_unpacklo_epi32(m, m)); }
#else
using sf_v = double;
constexpr int SF_W = 1;
static inline sf_v sf_ld(const double* p) { return *p; }
static inline void sf_st(double* p, const sf_v a) { *p = a; }
static inline sf_v sf_set(const double a) { return a; }
static inline sf_v sf_add(const sf_v a, const sf_v b) { return a + b; }
static inline sf_v sf_mul(const sf_v a, const sf_v b) { return a * b; }
static inline sf_v sf_fma(const sf_v c, const sf_v a, const sf_v b) { return c + a * b; }
static inline sf_v sf_fms(const sf_v c, const sf_v a, const sf_v b) { return c - a * b; }
static inline double sf_hsum(const sf_v a) { return a; }
#endif
#if !(defined(__aarch64__) || defined(_M_ARM64)) && (defined(__AVX__) || defined(__SSE2__) || defined(_M_X64))
//sincos_shared (convenience.h) on the SF_W x86 lanes: the same reduction and kernels, the quadrant from cvtpd_epi32
//(round to nearest, as nearbyint) as int32 masks widened to the 64-bit lanes. |x| < 1e5 only, sf_sincos_n tests it.
static inline void sf_sincos_v(const sf_v x, sf_v* s, sf_v* c)
{
	const __m128i qi = sf_round_i(sf_mul(x, sf_set(0.63661977236758134308)));
	const sf_v q = sf_from_i(qi);
	const sf_v r = sf_fms(sf_fms(sf_fms(x, q, sf_set(1.57079632673412561417e+00)), q, sf_set(6.07710050630396597660e-11)),
		q, sf_set(2.02226624871116645580e-21));
	const sf_v z = sf_mul(r, r);
	sf_v ps = sf_fma(sf_set(-2.50507602534068634195e-08), z, sf_set(1.58969099521155010221e-10));
	ps = sf_fma(sf_set(2.75573137070700676789e-06), z, ps);
	ps = sf_fma(sf_set(-1.98412698298579493134e-04), z, ps);
	ps = sf_fma(sf_set(8.33333333332248946124e-03), z, ps);
	ps = sf_fma(sf_set(-1.66666666666666324348e-01), z, ps);
	const sf_v sr = sf_fma(r, sf_mul(r, z), ps);
	sf_v pc = sf_fma(sf_set(2.08757232129817482790e-09), z, sf_set(-1.13596475577881948265e-11));
	pc = sf_fma(sf_set(-2.75573143513906633035e-07), z, pc);
	pc = sf_fma(sf_set(2.48015872894767294178e-05), z, pc);
	pc = sf_fma(sf_set(-1.38888888888741095749e-03), z, pc);
	pc = sf_fma(sf_set(4.16666666666666019037e-02), z, pc);
	const sf_v hz = sf_mul(sf_set(0.5), z), w = sf_sub(sf_set(1.0), hz);
	const sf_v cr = sf_add(w, sf_fma(sf_sub(sf_sub(sf_set(1.0), w), hz), sf_mul(z, z), pc));
	//quadrant q & 3: odd swaps sin and cos, q & 2 negates sin, (q + 1) & 2 negates cos
	const __m128i one = _mm_set1_epi32(1);
	const sf_v odd = sf_wide(_mm_cmpeq_epi32(_mm_and_si128(qi, one), one)), sign = sf_set(-0.0);
	const sf_v sgn_s = sf_and(sf_wide(_mm_slli_epi32(qi, 30)), sign);
	const sf_v sgn_c = sf_and(sf_wide(_mm_slli_epi32(_mm_add_epi32(qi, one), 30)), sign);
	*s = sf_xor(sf_or(sf_and(odd, cr), sf_andnot(odd, sr)), sgn_s);
	*c = sf_xor(sf_or(sf_and(odd, sr), sf_andnot(odd, cr)), sgn_c);
}
#endif
//n sines and cosines from one shared reduction: the NEON kernels on ARM64, sf_sincos_v on x86, else sincos_shared
//lane by lane (half the time of glibc's sincos)
static inline void sf_sincos_n(const int n, const double* x, double* s, double* c)
{
#if defined(__aarch64__) || defined(_M_ARM64)
	sincos_shared_n(n, x, s, c);
#else
	int p = 0;
#if defined(__AVX__) || defined(__SSE2__) || defined(_M_X64)
	for (; p + SF_W <= n; p += SF_W)
	{
		const sf_v v = sf_ld(x + p);
		if (sf_small(v)) { sf_v sv, cv; sf_sincos_v(v, &sv, &cv); sf_st(s + p, sv); sf_st(c + p, cv); }
		else for (int e = p; e < p + SF_W; e++) sincos_shared(x[e], s + e, c + e);
	}
#endif
	for (; p < n; p++) sincos_shared(x[p], s + p, c + p);
#endif
}
//M consecutive recurrence steps of calc_SF in one pass over a tile: Q stays in registers between steps, so the A72's
//~16 bytes/cycle of vector results go to multiplies instead of reloading Q. Two vectors run side by side to hide the
//chain latency. Every product and sum is the one a single step forms, in the same order: bit-identical to M = 1
//(MSVC /fp:fast may reassociate the intrinsic sums and differ by a few ulp).
template <int M>
static inline void sf_steps(const int nt, double* qr, double* qi, const double* tr, const double* ti, double* acc_r, double* acc_i)
{
	const int ntw = nt - nt % SF_W;
	sf_v sr[M], sm[M];
	for (int k = 0; k < M; k++) sr[k] = sm[k] = sf_set(0.0);
	int q = 0;
	for (; q + 2 * SF_W <= ntw; q += 2 * SF_W)
	{
		sf_v a0 = sf_ld(qr + q), b0 = sf_ld(qi + q), a1 = sf_ld(qr + q + SF_W), b1 = sf_ld(qi + q + SF_W);
		const sf_v c0 = sf_ld(tr + q), d0 = sf_ld(ti + q), c1 = sf_ld(tr + q + SF_W), d1 = sf_ld(ti + q + SF_W);
		for (int k = 0; k < M; k++)
		{
			const sf_v nr0 = sf_fms(sf_mul(a0, c0), b0, d0), ni0 = sf_fma(sf_mul(a0, d0), b0, c0);
			const sf_v nr1 = sf_fms(sf_mul(a1, c1), b1, d1), ni1 = sf_fma(sf_mul(a1, d1), b1, c1);
			sr[k] = sf_add(sf_add(sr[k], nr0), nr1); sm[k] = sf_add(sf_add(sm[k], ni0), ni1);
			a0 = nr0; b0 = ni0; a1 = nr1; b1 = ni1;
		}
		sf_st(qr + q, a0); sf_st(qi + q, b0); sf_st(qr + q + SF_W, a1); sf_st(qi + q + SF_W, b1);
	}
	if (q < ntw)
	{
		sf_v a = sf_ld(qr + q), b = sf_ld(qi + q);
		const sf_v c2 = sf_ld(tr + q), d = sf_ld(ti + q);
		for (int k = 0; k < M; k++)
		{
			const sf_v nr = sf_fms(sf_mul(a, c2), b, d), ni = sf_fma(sf_mul(a, d), b, c2);
			sr[k] = sf_add(sr[k], nr); sm[k] = sf_add(sm[k], ni);
			a = nr; b = ni;
		}
		sf_st(qr + q, a); sf_st(qi + q, b);
	}
	double re1[M] = {}, im1[M] = {};
	for (int e = ntw; e < nt; e++)
	{
		double a = qr[e], b = qi[e];
		for (int k = 0; k < M; k++)
		{
			const double nr = a * tr[e] - b * ti[e], ni = a * ti[e] + b * tr[e];
			a = nr; b = ni;
			re1[k] += a; im1[k] += b;
		}
		qr[e] = a; qi[e] = b;
	}
	for (int k = 0; k < M; k++)
	{
		acc_r[k] += sf_hsum(sr[k]) + re1[k];
		acc_i[k] += sf_hsum(sm[k]) + im1[k];
	}
}
//Row jump: Q = R e^{i db.r} (e^{i dc.r}, or its conjugate for n < 0)^|n| from the saved row start R, summed like a step.
static inline void sf_jump(const int nt, const int n, const double* rr, const double* ri, const double* br, const double* bi,
	const double* tr, const double* ti, double* qr, double* qi, double* acc_r, double* acc_i)
{
	const int ntw = nt - nt % SF_W, an = n < 0 ? -n : n;
	const double sg = n < 0 ? -1.0 : 1.0;
	const sf_v sgv = sf_set(sg);
	sf_v sr = sf_set(0.0), sm = sf_set(0.0);
	for (int q = 0; q < ntw; q += SF_W)
	{
		const sf_v a0 = sf_ld(rr + q), b0 = sf_ld(ri + q), c = sf_ld(br + q), d = sf_ld(bi + q);
		const sf_v c2 = sf_ld(tr + q), d2 = sf_mul(sgv, sf_ld(ti + q));
		sf_v a = sf_fms(sf_mul(a0, c), b0, d), b = sf_fma(sf_mul(a0, d), b0, c);
		for (int k = 0; k < an; k++)
		{
			const sf_v na = sf_fms(sf_mul(a, c2), b, d2);
			b = sf_fma(sf_mul(a, d2), b, c2);
			a = na;
		}
		sf_st(qr + q, a); sf_st(qi + q, b);
		sr = sf_add(sr, a); sm = sf_add(sm, b);
	}
	double re1 = 0.0, im1 = 0.0;
	for (int e = ntw; e < nt; e++)
	{
		double a = rr[e] * br[e] - ri[e] * bi[e], b = rr[e] * bi[e] + ri[e] * br[e];
		const double c2 = tr[e], d2 = sg * ti[e];
		for (int k = 0; k < an; k++)
		{
			const double na = a * c2 - b * d2;
			b = a * d2 + b * c2;
			a = na;
		}
		qr[e] = a; qi[e] = b;
		re1 += a; im1 += b;
	}
	*acc_r += sf_hsum(sr) + re1;
	*acc_i += sf_hsum(sm) + im1;
}
//Anchor: Q = rho e^{i k.r} from a sincos per point, summed like a step
static inline void sf_anchor(const int nt, const double k1, const double k2, const double k3, const double* x, const double* y,
	const double* z, const double* r, double* qr, double* qi, double* acc_r, double* acc_i)
{
	alignas(32) double ph[SF_TILE], sn[SF_TILE], cs[SF_TILE];
	const int ntw = nt - nt % SF_W;
	const sf_v k1v = sf_set(k1), k2v = sf_set(k2), k3v = sf_set(k3);
	int q = 0;
	for (; q < ntw; q += SF_W)
		sf_st(ph + q, sf_fma(sf_fma(sf_mul(k1v, sf_ld(x + q)), k2v, sf_ld(y + q)), k3v, sf_ld(z + q)));
	for (; q < nt; q++)
		ph[q] = k1 * x[q] + k2 * y[q] + k3 * z[q];
	sf_sincos_n(nt, ph, sn, cs);
	sf_v sr = sf_set(0.0), sm = sf_set(0.0);
	for (q = 0; q + 2 * SF_W <= ntw; q += 2 * SF_W)
	{
		const sf_v r0 = sf_ld(r + q), r1 = sf_ld(r + q + SF_W);
		const sf_v a0 = sf_mul(r0, sf_ld(cs + q)), b0 = sf_mul(r0, sf_ld(sn + q));
		const sf_v a1 = sf_mul(r1, sf_ld(cs + q + SF_W)), b1 = sf_mul(r1, sf_ld(sn + q + SF_W));
		sf_st(qr + q, a0); sf_st(qi + q, b0); sf_st(qr + q + SF_W, a1); sf_st(qi + q + SF_W, b1);
		sr = sf_add(sf_add(sr, a0), a1); sm = sf_add(sf_add(sm, b0), b1);
	}
	if (q < ntw)
	{
		const sf_v rv = sf_ld(r + q), a = sf_mul(rv, sf_ld(cs + q)), b = sf_mul(rv, sf_ld(sn + q));
		sf_st(qr + q, a); sf_st(qi + q, b);
		sr = sf_add(sr, a); sm = sf_add(sm, b);
	}
	double re1 = 0.0, im1 = 0.0;
	for (int e = ntw; e < nt; e++)
	{
		qr[e] = r[e] * cs[e]; qi[e] = r[e] * sn[e];
		re1 += qr[e]; im1 += qi[e];
	}
	*acc_r += sf_hsum(sr) + re1;
	*acc_i += sf_hsum(sm) + im1;
}


std::vector<char> calc_SF_fourier_bessel(const vec2& k_pt, const vec2& dens, const sf_shell_layout* layout,
	cvec2& sf, const int L_override, const bool exact, const bool force)
{
	const int na = static_cast<int>(dens.size());
	std::vector<char> done(na, 0);
	if (!layout || static_cast<int>(layout->proto.size()) != na)
		return done;
	const long long nk = static_cast<long long>(k_pt[0].size());
	const double* kx = k_pt[0].data(), * ky = k_pt[1].data(), * kz = k_pt[2].data();
	vec qk(nk);
	double qmax = 0.0;
	for (long long k = 0; k < nk; k++)
		qmax = std::max(qmax, qk[k] = std::sqrt(kx[k] * kx[k] + ky[k] * ky[k] + kz[k] * kz[k]));
	//reflections by |k|: a run of them on the same table nodes is one assembly GEMM
	std::vector<long long> ko(nk);
	std::iota(ko.begin(), ko.end(), 0LL);
	std::sort(ko.begin(), ko.end(), [&](const long long u, const long long v) { return qk[u] < qk[v]; });
	vec qs(nk);
	for (long long k = 0; k < nk; k++)
		qs[k] = qk[ko[k]];
	const std::vector<long long> run = fb_runs(qs);
	const long long nd = static_cast<long long>(run.size()) - 1;
	double ms_shells = 0.0, ms_proj = 0.0, ms_tab = 0.0, ms_asm = 0.0; //-gflops phase split
	constexpr size_t FB_MEM = size_t(1) << 24; //doubles per work array (128 MB): sets the atoms per column block
	constexpr int FB_BIN = 256, FB_ROWS = 64, FB_PB = 512; //reflections per assembly GEMM, table nodes per radial GEMM,
	//rule points per Y block of the projection
	//Direct-sum point x reflection steps per FB flop, timed on sucrose and rubredoxin (8 threads, AVX2, MKL)
	constexpr double FB_STEP_PER_FLOP = 0.05;
	_time_point t_ph = get_time();
	//A prototype's atoms are GEMM columns; it keeps its own shells, cut and Lebedev rules. All prototypes share the
	//|k| nodes (the finest any asks for) and Y(k^) up to the largest L, so the assembly evaluates Y once per reflection
	struct fb_rule { int s0, n, Lr; ivec sh; }; //first shell, points, largest Ls of its shells, its shells
	struct fb_part { int P = 0, ns = 0, L = 0; double as = 0.0, rad = 0.0; ivec mem, Ls, sl, rule_of; vec r; std::vector<fb_rule> rules; };
	std::vector<fb_part> parts;
	int L = 0;
	double h = 0.025;
	for (int P = 0; P < static_cast<int>(layout->protos.size()); P++)
	{
		const sf_shell_layout::prototype& pr = layout->protos[P];
		const ivec& ss = pr.shell_start;
		fb_part F;
		F.P = P;
		double direct = 0.0;
		for (int a = 0; a < na; a++)
		{
			const ivec& ix = layout->idx[a];
			if (layout->proto[a] != P || ix.empty() || ix.size() != dens[a].size())
				continue;
			F.mem.push_back(a);
			//shells up to the outermost one any atom kept a point on
			F.ns = std::max(F.ns, static_cast<int>(std::upper_bound(ss.begin(), ss.end(), *std::max_element(ix.begin(), ix.end())) - ss.begin()));
			direct += double(ix.size()) * double(nk);
		}
		if (F.mem.empty())
			continue;
		const int ns = F.ns;
		F.r.resize(ns);
		int most = 0;
		for (int s = 0; s < ns; s++)
		{
			const int p = ss[s];
			F.r[s] = std::sqrt(pr.x[p] * pr.x[p] + pr.y[p] * pr.y[p] + pr.z[p] * pr.z[p]);
			most = std::max(most, ss[s + 1] - ss[s]);
		}
		F.L = L_override >= 0 ? L_override : fb_lebedev_degree(most) / 2;
		//4-point Lagrange in |k|: error ~ (h r)^4 on a shell's term, so h shrinks with the outermost shell
		const double hP = std::min(0.025, 0.25 / std::max(F.r[ns - 1], 1e-3));
		const int nqP = static_cast<int>(qmax / hP) + 4;
		//j_l(x) < 1e-12 j_0 past x + 9.4 x^(1/3) + 6 (Debye's asymptotics, the 12-digit rule of fast multipole codes),
		//so a shell keeps only the l its |k| range can reach: exact to that level, and inner shells drop most of L
		F.Ls.resize(ns);
		F.as = 2.0 * double(nk) * double(F.L + 1) * (F.L + 1);
		double flops = 0.0;
		for (int s = 0; s < ns; s++)
		{
			const double x = qmax * F.r[s];
			F.Ls[s] = std::min(F.L, static_cast<int>(std::ceil(x + 9.4 * std::cbrt(x))) + 6);
			const double nls = double(F.Ls[s] + 1) * (F.Ls[s] + 1);
			F.rad += 2.0 * nls;
			flops += 2.0 * nls * double(ss[s + 1] - ss[s]);
		}
		//table: 4 nodes per reflection on nqP; exact: one per distinct |k|
		flops += std::min(exact ? 1e300 : 4.0 * F.as + F.rad * nqP, F.as + F.rad * double(nd));
		//the direct sum keeps a prototype it does faster; asking for L or exact forces Fourier-Bessel
		if (!force && L_override < 0 && !exact && FB_STEP_PER_FLOP * flops * double(F.mem.size()) > direct)
			continue;
		L = std::max(L, F.L);
		h = std::min(h, hP);
		parts.push_back(std::move(F));
	}
	if (parts.empty())
		return done;
	const int ne = fb_n_even(L), no = fb_n_odd(L), nl = ne + no;
	const int nq = static_cast<int>(qmax / h) + 4;
	double ct = 0.0, ce = 0.0;
	for (const fb_part& F : parts)
	{
		ct += double(F.mem.size()) * (4.0 * F.as + F.rad * nq);
		ce += double(F.mem.size()) * (F.as + F.rad * double(nd));
	}
	const bool table = !exact && ct < ce;
	const int S = table ? 4 : 1; //nodes per reflection; exact: the distinct |k| are the nodes
	const long long nn = table ? nq : nd;
	const int NC = table ? static_cast<int>(nn) : static_cast<int>(std::min<long long>(nd, 2048));
	const fb_ylm Y(L);
	int Mtot = 0;
	size_t widest = size_t(NC) * nl;
	for (fb_part& F : parts)
	{
		for (const int a : F.mem)
			done[a] = 1;
		const sf_shell_layout::prototype& pr = layout->protos[F.P];
		const ivec& ss = pr.shell_start;
		const int ns = F.ns;
		F.sl.assign(L + 1, ns); //first shell reaching l: J_l is zero on every shell before it
		for (int l = 0; l <= L; l++)
			for (int s = 0; s < ns && F.sl[l] == ns; s++)
				if (F.Ls[s] >= l)
					F.sl[l] = s;
		//Shells on the same Lebedev rule and rotation share Y_lm(u_p): the projection takes all of a rule's shells at once
		F.rule_of.resize(ns);
		for (int s = 0; s < ns; s++)
		{
			const int n = ss[s + 1] - ss[s], p1 = ss[s];
			const double i1 = F.r[s] > 0.0 ? 1.0 / F.r[s] : 0.0;
			int R = -1;
			for (int t = 0; t < static_cast<int>(F.rules.size()) && R < 0; t++)
			{
				if (F.rules[t].n != n)
					continue;
				const int p0 = ss[F.rules[t].s0];
				const double i0 = F.r[F.rules[t].s0] > 0.0 ? 1.0 / F.r[F.rules[t].s0] : 0.0;
				bool same = true;
				for (int i = 0; i < n && same; i++)
					same = std::abs(pr.x[p0 + i] * i0 - pr.x[p1 + i] * i1) < 1e-10 && std::abs(pr.y[p0 + i] * i0 - pr.y[p1 + i] * i1) < 1e-10
					&& std::abs(pr.z[p0 + i] * i0 - pr.z[p1 + i] * i1) < 1e-10;
				if (same)
					R = t;
			}
			if (R < 0)
			{
				R = static_cast<int>(F.rules.size());
				F.rules.push_back({ s, n, 0, {} });
			}
			F.rules[R].sh.push_back(s);
			F.rules[R].Lr = std::max(F.rules[R].Lr, F.Ls[s]);
			F.rule_of[s] = R;
		}
		Mtot += static_cast<int>(F.mem.size());
		widest = std::max({ widest, size_t(ns) * nl, size_t(ss[ns]) });
	}
	//J[l][n][s] = 4pi (+-1) j_l(q_n r_s), zero past Ls[s]: i^l = +-1 (even l) or +-i (odd), the i left to the parity
	vec J;
	const auto fill_J = [&](const fb_part& F, const long long n_lo, const int cnt) {
		const int ns = F.ns;
		J.assign(size_t(F.L + 1) * cnt * ns, 0.0);
#pragma omp parallel
		{
			vec jl(F.L + 1);
#pragma omp for schedule(static)
			for (int n = 0; n < cnt; n++)
			{
				const double q = table ? double(n_lo + n) * h : qs[run[n_lo + n]];
				for (int s = 0; s < ns; s++)
				{
					fb_sph_bessel(F.Ls[s], q * F.r[s], jl.data());
					for (int l = 0; l <= F.Ls[s]; l++)
						J[(size_t(l) * cnt + n) * ns + s] = constants::FOUR_PI * (((l >> 1) & 1) ? -1.0 : 1.0) * jl[l];
				}
			}
		}
	};
	//assembly runs: reflections (in |k| order) on the same nodes n0..n0+S-1, at most FB_BIN each
	struct fb_bin { long long k0; int nb, n0; };
	std::vector<fb_bin> bins;
	if (table)
		for (long long k = 0; k < nk;)
		{
			const int n0 = std::clamp(static_cast<int>(qk[ko[k]] / h) - 1, 0, nq - 4);
			long long e = k + 1;
			while (e < nk && e - k < FB_BIN && std::clamp(static_cast<int>(qk[ko[e]] / h) - 1, 0, nq - 4) == n0)
				e++;
			bins.push_back({ k, static_cast<int>(e - k), n0 });
			k = e;
		}
	const int Mc = std::max(1, static_cast<int>(std::min<size_t>(Mtot, FB_MEM / widest)));
	ms_shells += get_msec(t_ph, get_time());
	//a column block holds whole or partial prototypes: segment = part, its first member, members, first column
	struct fb_seg { int f, lo, mp, col; };
	vec W, Ce, Co, Te, To, Yb, Ct;
	for (int c0 = 0; c0 < Mtot; c0 += Mc)
	{
		const int mc = std::min(Mc, Mtot - c0);
		std::vector<fb_seg> seg;
		for (int f = 0, g = 0; f < static_cast<int>(parts.size()); g += static_cast<int>(parts[f].mem.size()), f++)
		{
			const int lo = std::max(c0, g), hi = std::min(c0 + mc, g + static_cast<int>(parts[f].mem.size()));
			if (lo < hi)
				seg.push_back({ f, lo - g, hi - lo, lo - c0 });
		}
		t_ph = get_time();
		//c_s,lm = sum_p Y_lm(u_p) W_p: per rule one GEMM over the rule's points against its shells x atoms side by side
		//(W zero where pruned), Y evaluated in blocks of points and the moments scattered to [shell][lm][atom]
		std::vector<size_t> ceo(seg.size() + 1, 0), coo(seg.size() + 1, 0);
		for (size_t j = 0; j < seg.size(); j++)
		{
			ceo[j + 1] = ceo[j] + size_t(parts[seg[j].f].ns) * ne * seg[j].mp;
			coo[j + 1] = coo[j] + size_t(parts[seg[j].f].ns) * no * seg[j].mp;
		}
		Ce.assign(ceo.back(), 0.0);
		Co.assign(coo.back(), 0.0);
		for (size_t j = 0; j < seg.size(); j++)
		{
			const fb_part& F = parts[seg[j].f];
			const sf_shell_layout::prototype& pr = layout->protos[F.P];
			const ivec& ss = pr.shell_start;
			const int mp = seg[j].mp, nr = static_cast<int>(F.rules.size());
			std::vector<size_t> wofs(nr + 1, 0);
			for (int t = 0; t < nr; t++)
				wofs[t + 1] = wofs[t] + size_t(F.rules[t].n) * F.rules[t].sh.size() * mp;
			ivec slot(F.ns), sop(ss[F.ns]); //a shell's place among its rule's shells; a prototype point's shell
			for (const fb_rule& R : F.rules)
				for (size_t q = 0; q < R.sh.size(); q++)
					slot[R.sh[q]] = static_cast<int>(q);
			for (int s = 0; s < F.ns; s++)
				std::fill(sop.begin() + ss[s], sop.begin() + ss[s + 1], s);
			W.assign(wofs.back(), 0.0); //W_t[i][slot][atom]
			for (int c = 0; c < mp; c++)
			{
				const int a = F.mem[seg[j].lo + c];
				const ivec& ix = layout->idx[a];
				for (size_t p = 0; p < ix.size(); p++)
				{
					const int s = sop[ix[p]], t = F.rule_of[s];
					W[wofs[t] + (size_t(ix[p] - ss[s]) * F.rules[t].sh.size() + slot[s]) * mp + c] = dens[a][p];
				}
			}
			for (int t = 0; t < nr; t++)
			{
				const fb_rule& R = F.rules[t];
				const int nsh = static_cast<int>(R.sh.size()), N = nsh * mp, re = fb_n_even(R.Lr), ro = fb_n_odd(R.Lr), rl = re + ro;
				const int p0 = ss[R.s0], nbk = (rl + FB_ROWS - 1) / FB_ROWS;
				//at r = 0 only l = 0 survives j_l(0)
				const double ir = F.r[R.s0] > 0.0 ? 1.0 / F.r[R.s0] : 0.0;
				Yb.resize(size_t(std::min(FB_PB, R.n)) * rl);
				Ct.resize(size_t(rl) * N);
				for (int b0 = 0; b0 < R.n; b0 += FB_PB)
				{
					const int nb = std::min(FB_PB, R.n - b0);
#pragma omp parallel
					{
#if !defined(__APPLE__) && !defined(NSA2_OPENBLAS)
						const int mkl_before = mkl_set_num_threads_local(1);
#endif
#pragma omp for schedule(static)
						for (int i = 0; i < nb; i++)
						{
							const int p = p0 + b0 + i;
							Y.eval(pr.x[p] * ir, pr.y[p] * ir, pr.z[p] * ir, Yb.data() + size_t(i) * rl, Yb.data() + size_t(i) * rl + re, R.Lr);
						}
#pragma omp for schedule(dynamic)
						for (int g = 0; g < nbk; g++)
						{
							const int r0 = g * FB_ROWS, rows = std::min(FB_ROWS, rl - r0);
							cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, rows, N, nb, 1.0, Yb.data() + r0, rl,
								W.data() + wofs[t] + size_t(b0) * N, N, b0 ? 1.0 : 0.0, Ct.data() + size_t(r0) * N, N);
						}
#if !defined(__APPLE__) && !defined(NSA2_OPENBLAS)
						mkl_set_num_threads_local(mkl_before);
#endif
					}
				}
#pragma omp parallel for schedule(static)
				for (int q = 0; q < nsh; q++)
				{
					const int s = R.sh[q];
					for (int i = 0; i < re; i++)
						std::copy_n(Ct.data() + size_t(i) * N + size_t(q) * mp, mp, Ce.data() + ceo[j] + (size_t(s) * ne + i) * mp);
					for (int i = 0; i < ro; i++)
						std::copy_n(Ct.data() + size_t(re + i) * N + size_t(q) * mp, mp, Co.data() + coo[j] + (size_t(s) * no + i) * mp);
				}
			}
		}
		ms_proj += get_msec(t_ph, get_time());
		for (long long n_lo = 0; n_lo < nn; n_lo += NC)
		{
			const int cnt = static_cast<int>(std::min<long long>(NC, nn - n_lo));
			t_ph = get_time();
			if (!table)
			{
				bins.clear();
				for (int j = 0; j < cnt; j++)
					for (long long k = run[n_lo + j]; k < run[n_lo + j + 1]; k += FB_BIN)
						bins.push_back({ k, static_cast<int>(std::min<long long>(FB_BIN, run[n_lo + j + 1] - k)), j });
			}
			//T[lm][n][atom] = sum_{s >= sl[l]} J_l[n][s] c_s,lm[atom]: one GEMM per l, block of nodes and segment, scattered
			//so that a reflection's S nodes are S mc consecutive columns of row lm. Zero past a prototype's own L.
			Te.assign(size_t(cnt) * ne * mc, 0.0);
			To.assign(size_t(cnt) * no * mc, 0.0);
			const int nbk = (cnt + FB_ROWS - 1) / FB_ROWS;
			for (size_t j = 0; j < seg.size(); j++)
			{
				const fb_part& F = parts[seg[j].f];
				const int ns = F.ns, mp = seg[j].mp;
				fill_J(F, n_lo, cnt);
#pragma omp parallel
				{
#if !defined(__APPLE__) && !defined(NSA2_OPENBLAS)
					const int mkl_before = mkl_set_num_threads_local(1);
#endif
					vec tb(size_t(FB_ROWS) * (2 * L + 1) * mp);
#pragma omp for schedule(dynamic)
					for (int t = 0; t < (L + 1) * nbk; t++)
					{
						const int l = L - t / nbk, n0 = (t % nbk) * FB_ROWS, rows = std::min(FB_ROWS, cnt - n0); //big l first
						if (F.sl[l] == ns)
							continue;
						const int width = (l & 1) ? no : ne, o = l * (l - 1) / 2, w = (2 * l + 1) * mp;
						const double* C = (l & 1) ? Co.data() + coo[j] : Ce.data() + ceo[j];
						cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rows, w, ns - F.sl[l], 1.0,
							J.data() + (size_t(l) * cnt + n0) * ns + F.sl[l], ns, C + (size_t(F.sl[l]) * width + o) * mp, width * mp,
							0.0, tb.data(), w);
						double* T = ((l & 1) ? To : Te).data();
						for (int m = 0; m < 2 * l + 1; m++)
							for (int n = 0; n < rows; n++)
								std::copy_n(tb.data() + size_t(n) * w + size_t(m) * mp, mp,
									T + ((size_t(o) + m) * cnt + n0 + n) * mc + seg[j].col);
					}
#if !defined(__APPLE__) && !defined(NSA2_OPENBLAS)
					mkl_set_num_threads_local(mkl_before);
#endif
				}
			}
			ms_tab += get_msec(t_ph, get_time());
			t_ph = get_time();
			//f(k) = sum_m w_m(|k|) Y(k^) . T[n0 + m]: Y rows times the S nodes x every column in one GEMM per run and
			//parity, the Lagrange weights applied after
			ivec colatom(mc);
			for (const fb_seg& G : seg)
				for (int c = 0; c < G.mp; c++)
					colatom[G.col + c] = parts[G.f].mem[G.lo + c];
#pragma omp parallel
			{
#if !defined(__APPLE__) && !defined(NSA2_OPENBLAS)
				const int mkl_before = mkl_set_num_threads_local(1);
#endif
				vec Ye(size_t(FB_BIN) * ne), Yo(size_t(FB_BIN) * no), Ze(size_t(FB_BIN) * S * mc), Zo(size_t(FB_BIN) * S * mc, 0.0),
					wt(size_t(FB_BIN) * 4);
#pragma omp for schedule(dynamic)
				for (int b = 0; b < static_cast<int>(bins.size()); b++)
				{
					const fb_bin& B = bins[b];
					for (int j = 0; j < B.nb; j++)
					{
						const long long k = ko[B.k0 + j];
						const double q = qk[k], iq = q > 0.0 ? 1.0 / q : 0.0;
						Y.eval(kx[k] * iq, ky[k] * iq, q > 0.0 ? kz[k] * iq : 1.0, Ye.data() + size_t(j) * ne, Yo.data() + size_t(j) * no);
						double* w = wt.data() + size_t(j) * 4;
						w[0] = 1.0;
						if (table)
						{
							const double u = q / h - B.n0;
							w[0] = -(u - 1) * (u - 2) * (u - 3) / 6.0; w[1] = u * (u - 2) * (u - 3) / 2.0;
							w[2] = -u * (u - 1) * (u - 3) / 2.0; w[3] = u * (u - 1) * (u - 2) / 6.0;
						}
					}
					const int N = S * mc;
					cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, B.nb, N, ne, 1.0, Ye.data(), ne,
						Te.data() + size_t(B.n0) * mc, size_t(cnt) * mc, 0.0, Ze.data(), N);
					if (no > 0)
						cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, B.nb, N, no, 1.0, Yo.data(), no,
							To.data() + size_t(B.n0) * mc, size_t(cnt) * mc, 0.0, Zo.data(), N);
					for (int j = 0; j < B.nb; j++)
					{
						const long long k = ko[B.k0 + j];
						const double* w = wt.data() + size_t(j) * 4, * ze = Ze.data() + size_t(j) * N, * zo = Zo.data() + size_t(j) * N;
						for (int c = 0; c < mc; c++)
						{
							double re = 0.0, im = 0.0;
							for (int m = 0; m < S; m++)
							{
								re += w[m] * ze[m * mc + c];
								im += w[m] * zo[m * mc + c];
							}
							sf[colatom[c]][k] = std::complex<double>(re, im);
						}
					}
				}
#if !defined(__APPLE__) && !defined(NSA2_OPENBLAS)
				mkl_set_num_threads_local(mkl_before);
#endif
			}
			ms_asm += get_msec(t_ph, get_time());
		}
	}
	throughput::record_time("FB shells + cost model", false, ms_shells);
	throughput::record_time("FB projection", false, ms_proj);
	throughput::record_time("FB radial table", false, ms_tab);
	throughput::record_time("FB assembly", false, ms_asm);
	return done;
}

/**
 * Calculates the scattering factors for a given set of parameters.
 *
 * @param points The number of points.
 * @param k_pt The vector of k points.
 * @param d1 The vector of distances in x.
 * @param d2 The vector of distances in y.
 * @param d3 The vector of distances in z.
 * @param dens The vector of density values.
 * @param sf The vector of scattering factors.
 * @param file The output stream to write the results.
 * @param start The start time of the calculation.
 * @param end1 The end time of the calculation.
 * @param debug Flag indicating whether to enable debug mode.
 * @param no_date Flag indicating whether to exclude the date in the output.
 */
void calc_SF(const int& points,
	const vec2& k_pt,
	const vec2& d1,
	const vec2& d2,
	const vec2& d3,
	const vec2& dens,
	cvec2& sf,
	std::ostream& file,
	_time_point& start,
	_time_point& end1,
	bool debug,
	bool no_date,
	bool do_XCW,
	bool use_gpu,
	bool gpu_fp64,
	bool gpu_fp32,
	const sf_shell_layout* layout)
{
	const long long int imax = static_cast<long long int>(dens.size());
	const long long int smax = static_cast<long long int>(k_pt[0].size());
	sf.reserve(imax * smax);
	sf.resize(imax);
#pragma omp parallel for
	for (int i = 0; i < imax; i++)
		sf[i].resize(smax, constants::cnull);
	using namespace std;
	if (debug)
		file << "Initialized FFs" << std::endl
		<< "asym atom list size: " << imax << " total grid size: " << points << endl;
#ifdef NOSPHERA2_USE_GPU
	if (use_gpu)
		sf_gpu_warmup_wait();
#endif
	end1 = get_time();

	if (!no_date)
	{
		long long int dur = get_sec(start, end1);
		if (dur < 1)
			file << "Time to prepare: " << fixed << setprecision(0) << get_msec(start, end1) << " ms" << endl << endl;
		else
			file << "Time to prepare: " << fixed << setprecision(0) << dur << " s" << endl << endl;
	}
	//Fourier-Bessel on the prototype shells unless -tune NOS_SF_DIRECT; atoms without a layout (cube grids) and those
	//the direct sum does faster go on to it, on the device when there is one
	std::vector<char> fb_done(imax, 0);
	if (layout && !tuning("NOS_SF_DIRECT")) {
		const _time_point fb_t0 = get_time();
		const char* fb_L = tuning("NOS_SF_FB_L");
		fb_done = calc_SF_fourier_bessel(k_pt, dens, layout, sf, fb_L ? std::atoi(fb_L) : -1,
			tuning("NOS_SF_FB_EXACT") != nullptr, tuning("NOS_SF_FB_FORCE") != nullptr);
		throughput::record_time("scattering-factor transform (Fourier-Bessel)", false, get_msec(fb_t0, get_time()));
		if (debug)
			file << "Fourier-Bessel transform for " << std::count(fb_done.begin(), fb_done.end(), 1) << " of " << imax
			<< " atoms, the direct sum for the rest" << endl;
		if (std::all_of(fb_done.begin(), fb_done.end(), [](const char c) { return c != 0; })) {
			if (!do_XCW) {
				ProgressBar fprogress(imax, 60, "=", " ", "Calculating Scattering Factors", file);
				for (int i = 0; i < imax; i++)
					fprogress.update();
			}
			return;
		}
	}
#ifdef NOSPHERA2_USE_GPU
	//-gpu_fp64 wins over -gpu_fp32 if both are given: the accurate one is the safer default.
	const sf_precision prec = gpu_fp64 ? sf_precision::FP64
		: gpu_fp32 ? sf_precision::FP32 : sf_precision::Auto;
	//An APU's own cores beat its device here unless a precision is asked for (Radeon 780M f32: 176 ms at d 0.5 sucrose,
	//its Ryzen 8700G 106 ms)
	if (use_gpu && sf_gpu_available() && (prec != sf_precision::Auto || !itensor_gpu_integrated())) {
		ivec left; //the atoms Fourier-Bessel left
		for (int i = 0; i < imax; i++)
			if (!fb_done[i])
				left.push_back(i);
		const int nleft = (int)left.size();
		ivec offs(nleft + 1, 0);
		for (int j = 0; j < nleft; j++)
			offs[j + 1] = offs[j] + (int)dens[left[j]].size();
		const long long tot = offs[nleft];
		vec fd1(tot), fd2(tot), fd3(tot), fde(tot);
		std::vector<double*> rows(nleft);
		for (int j = 0; j < nleft; j++) {
			const int i = left[j], lo = offs[j], n = (int)dens[i].size();
			for (int p = 0; p < n; p++) {
				fd1[lo + p] = d1[i][p]; fd2[lo + p] = d2[i][p];
				fd3[lo + p] = d3[i][p]; fde[lo + p] = dens[i][p];
			}
			rows[j] = reinterpret_cast<double*>(sf[i].data());
		}
		const _time_point sf_gpu_t0 = get_time();
		if (sf_gpu_run(nleft, smax, k_pt[0].data(), k_pt[1].data(), k_pt[2].data(),
			fd1.data(), fd2.data(), fd3.data(), fde.data(), offs.data(), tot,
			rows.data(), prec)) {
			//Rate includes the transfers: that is what the caller gets.
			throughput::record("scattering-factor transform", true,
				throughput::flops_ndft(static_cast<double>(tot), static_cast<double>(smax)),
				get_msec(sf_gpu_t0, get_time()));
			//The bar is what the reference logs expect, so draw it even though the work is done
			if (!do_XCW) {
				ProgressBar gprogress(imax, 60, "=", " ", "Calculating Scattering Factors", file);
				for (int i = 0; i < imax; i++)
					gprogress.update();
			}
			if (!no_date) {
				_time_point gend = get_time();
				const int ratio = sf_gpu_fp64_ratio();
				file << "GPU in use: scattering-factor Fourier transform on " << sf_gpu_backend() << ": " << get_msec(end1, gend) << " ms ("
					 << (sf_gpu_uses_fp32(prec) ? "centre-relative f32 phase and sincos" : "f64 sincos")
					 << ", fp32:fp64 ratio " << ratio << ")" << std::endl;
			}
			return;
		}
	}
#endif
	ProgressBar* progress = nullptr;
	if (!do_XCW) {
		progress = new ProgressBar(imax, 60, "=", " ", "Calculating Scattering Factors", file);
	}
	long long int pmax;
	const double* d1_local, * d2_local, * d3_local, * dens_local;

	// Pre-fetch k_pt data pointers for better cache locality
	const double* k1_data = k_pt[0].data();
	const double* k2_data = k_pt[1].data();
	const double* k3_data = k_pt[2].data();

	//Timed around the atom loop: a timer inside the omp parallel for would sum thread time beyond the wall clock.
	const _time_point sf_cpu_t0 = get_time();
	double sf_cpu_points = 0.0;

	//Phase recurrence: the reflections are a sorted hkl set, so most consecutive k differ by one step dc (c*), and
	//Q_p = rho_p e^{i k.r_p} continues as Q_p e^{i dc.r_p}, one complex multiply instead of a sincos. Q is re-anchored
	//with a sincos at block starts, row breaks no row jump (below) reaches and every SF_ANCHOR multiplies, so the
	//drift stays at rounding level (1e-14 relative on sucrose against the direct sum).
	//dc is taken where two successive steps agree, at the smallest |k| for the least rounding; cont[s] marks
	//k[s] = k[s-1] + dc. An unsorted k list leaves cont all false and every reflection anchors.
	constexpr signed char SF_STEP = 100, SF_ANC = 101; //sf_op values; anything else is a row jump by that many dc steps
	double dc[3] = { 0.0, 0.0, 0.0 }, db[3] = { 0.0, 0.0, 0.0 };
	bool jumps = false;
	std::vector<char> cont(smax, 0);
	std::vector<signed char> sf_op(smax, SF_ANC);
	{
		double kmax = 0.0, best = std::numeric_limits<double>::max();
		for (long long t = 0; t < smax; t++)
			kmax = std::max({ kmax, std::abs(k1_data[t]), std::abs(k2_data[t]), std::abs(k3_data[t]) });
		const double tol = 1e-13 * (1.0 + kmax);
		const auto steps_by = [&](const long long t, const double* d) {
			return std::abs(k1_data[t] - k1_data[t - 1] - d[0]) <= tol && std::abs(k2_data[t] - k2_data[t - 1] - d[1]) <= tol
				&& std::abs(k3_data[t] - k3_data[t - 1] - d[2]) <= tol;
		};
		for (long long t = 2; t < smax; t++) {
			const double d[3] = { k1_data[t - 1] - k1_data[t - 2], k2_data[t - 1] - k2_data[t - 2], k3_data[t - 1] - k3_data[t - 2] };
			const double n = k1_data[t - 1] * k1_data[t - 1] + k2_data[t - 1] * k2_data[t - 1] + k3_data[t - 1] * k3_data[t - 1];
			if (n < best && steps_by(t, d)) { best = n; std::copy(d, d + 3, dc); }
		}
		if (best < std::numeric_limits<double>::max())
			for (long long t = 1; t < smax; t++)
				cont[t] = steps_by(t, dc);
		//Row jumps: a row start that is the saved start of the previous row (or the block's first reflection) plus
		//db and n whole dc steps, |n| <= SF_NMAX, is reached as R e^{i db.r} e^{+-i dc.r}^|n| instead of anchored.
		//db is the commonest difference between consecutive row starts, reduced by whole dc steps. Every complex
		//multiply since the last sincos counts towards SF_ANCHOR, along the chain of saved row starts too.
		const double dcc = dc[0] * dc[0] + dc[1] * dc[1] + dc[2] * dc[2];
		const auto reduce = [&](const long long t, const long long u, const double* base, double* v) {
			v[0] = k1_data[t] - k1_data[u] - base[0]; v[1] = k2_data[t] - k2_data[u] - base[1]; v[2] = k3_data[t] - k3_data[u] - base[2];
			const long long n = std::llround((v[0] * dc[0] + v[1] * dc[1] + v[2] * dc[2]) / dcc);
			v[0] -= n * dc[0]; v[1] -= n * dc[1]; v[2] -= n * dc[2];
			return n;
		};
		const auto same = [&](const double* v, const double* w) {
			return std::abs(v[0] - w[0]) <= tol && std::abs(v[1] - w[1]) <= tol && std::abs(v[2] - w[2]) <= tol;
		};
		const double zero[3] = { 0.0, 0.0, 0.0 };
		if (dcc > 0.0)
		{
			std::vector<std::array<double, 3>> cand;
			for (long long t = 1, u = 0; t < smax; t++)
				if (!cont[t]) { double v[3]; reduce(t, u, zero, v); cand.push_back({ v[0], v[1], v[2] }); u = t; }
			size_t most = 1;
			for (size_t c = 0; c < std::min<size_t>(cand.size(), 32); c++)
			{
				const size_t hits = std::count_if(cand.begin(), cand.end(), [&](const std::array<double, 3>& w) { return same(cand[c].data(), w.data()); });
				if (hits > most) { most = hits; std::copy(cand[c].begin(), cand[c].end(), db); jumps = true; }
			}
		}
		for (long long b0 = 0; b0 < smax; b0 += SF_BLOCK)
		{
			const long long b1 = std::min<long long>(b0 + SF_BLOCK, smax);
			long long rs = b0;
			int since = 0, rs_since = 0;
			for (long long t = b0; t < b1; t++)
			{
				if (t > b0 && cont[t] && since < SF_ANCHOR) { sf_op[t] = SF_STEP; since++; continue; }
				since = 0;
				if (jumps && t > b0 && !cont[t])
				{
					double v[3];
					const long long n = reduce(t, rs, db, v);
					if (same(v, zero) && std::llabs(n) <= SF_NMAX && rs_since + 1 + std::llabs(n) <= SF_ANCHOR)
					{
						sf_op[t] = static_cast<signed char>(n);
						since = rs_since + 1 + static_cast<int>(std::llabs(n));
					}
				}
				if (t == b0 || !cont[t]) { rs = t; rs_since = since; }
			}
		}
	}
	vec e_ph, e_re, e_im, b_re, b_im;

	for (int i = 0; i < imax; i++)
	{
		if (fb_done[i]) {
			if (!do_XCW)
				progress->update();
			continue;
		}
		pmax = static_cast<long long int>(dens[i].size());
		sf_cpu_points += static_cast<double>(pmax);
		dens_local = dens[i].data();
		d1_local = d1[i].data();
		d2_local = d2[i].data();
		d3_local = d3[i].data();

		e_ph.resize(pmax); e_re.resize(pmax); e_im.resize(pmax);
		for (long long t = 0; t < pmax; t++)
			e_ph[t] = dc[0] * d1_local[t] + dc[1] * d2_local[t] + dc[2] * d3_local[t];
		sf_sincos_n(static_cast<int>(pmax), e_ph.data(), e_im.data(), e_re.data());
		if (jumps)
		{
			b_re.resize(pmax); b_im.resize(pmax);
			for (long long t = 0; t < pmax; t++)
				e_ph[t] = db[0] * d1_local[t] + db[1] * d2_local[t] + db[2] * d3_local[t];
			sf_sincos_n(static_cast<int>(pmax), e_ph.data(), b_im.data(), b_re.data());
		}
		const double* er = e_re.data(), * ei = e_im.data(), * bre = b_re.data(), * bim = b_im.data();
		//Blocks of SF_BLOCK reflections per task; per SF_TILE points (L1-resident Q) the block's reflections in order
#pragma omp parallel for schedule(dynamic, 1)
		for (long long b0 = 0; b0 < smax; b0 += SF_BLOCK)
		{
			const int nb = static_cast<int>(std::min<long long>(SF_BLOCK, smax - b0));
			double qr[SF_TILE], qi[SF_TILE], rr[SF_TILE], ri[SF_TILE], acc_r[SF_BLOCK] = {}, acc_i[SF_BLOCK] = {};
			for (long long t0 = 0; t0 < pmax; t0 += SF_TILE)
			{
				const int nt = static_cast<int>(std::min<long long>(SF_TILE, pmax - t0));
				const double* r = dens_local + t0, * x = d1_local + t0, * y = d2_local + t0, * z = d3_local + t0;
				const double* tr = er + t0, * ti = ei + t0;
				for (int j = 0; j < nb;)
				{
					const long long sj = b0 + j;
					if (sf_op[sj] == SF_STEP)
					{
						//steps j .. j+m-1 all continue the row, none reaching the anchor limit
						int m = 1;
						while (m < 4 && j + m < nb && sf_op[sj + m] == SF_STEP) m++;
						switch (m)
						{
						case 1: sf_steps<1>(nt, qr, qi, tr, ti, acc_r + j, acc_i + j); break;
						case 2: sf_steps<2>(nt, qr, qi, tr, ti, acc_r + j, acc_i + j); break;
						case 3: sf_steps<3>(nt, qr, qi, tr, ti, acc_r + j, acc_i + j); break;
						default: sf_steps<4>(nt, qr, qi, tr, ti, acc_r + j, acc_i + j); break;
						}
						j += m;
						continue;
					}
					if (sf_op[sj] != SF_ANC)
						sf_jump(nt, sf_op[sj], rr, ri, bre + t0, bim + t0, tr, ti, qr, qi, acc_r + j, acc_i + j);
					else
						sf_anchor(nt, k1_data[sj], k2_data[sj], k3_data[sj], x, y, z, r, qr, qi, acc_r + j, acc_i + j);
					//a row start is kept for the next row's jump
					if (j == 0 || !cont[sj])
					{
						std::copy(qr, qr + nt, rr);
						std::copy(qi, qi + nt, ri);
					}
					j++;
				}
			}
			complex<double>* out = sf[i].data() + b0;
			for (int j = 0; j < nb; j++)
				out[j] = complex<double>(acc_r[j], acc_i[j]);
		}
		if (!do_XCW) {
			progress->update();
		}
	}
	throughput::record("scattering-factor transform", false,
		throughput::flops_ndft(sf_cpu_points, static_cast<double>(smax)),
		get_msec(sf_cpu_t0, get_time()));
	if (!do_XCW) {
		delete (progress);
	}
}
