#pragma once
#include "pch.h"
#include <bit>

// Fast exp approximation for negative values
inline double fast_exp_neg(double x) {
	// For x in [-42, 0], use a fast approximation
	if (x < -42.0) return 0.0;
	if (x > -0.693147) { // ln(0.5) - use standard exp for values close to 0
		return exp(x);
	}
	// Power of 2 approximation: exp(x) ≈ (1 + x/1024)^1024
	x = 1.0 + x / 1024.0;
	x *= x; x *= x; x *= x; x *= x; x *= x; // x^32
	x *= x; x *= x; x *= x; x *= x; x *= x; // x^1024
	return x;
}

// sin and cos from one shared argument reduction, where the CRT reduces twice. Cody-Waite reduction by pi/2 and the
// FreeBSD k_sin/k_cos kernels: <= 2.2e-16 absolute error for |x| < 1e5; larger or non-finite x goes to the CRT (which
// MSVC /fp:fast swaps for a scalar sin good to ~4e-8 at |x| ~ 1e8). calc_SF runs the same kernels on vector lanes.
inline void sincos_shared(const double x, double* s, double* c) {
	if (!(std::abs(x) < 1e5)) { *s = std::sin(x); *c = std::cos(x); return; }
#ifdef __ANDROID__
	//armv7 has no inline nearbyint or double->int64: both are calls, 2.4x the kernel on the A53. The 1.5*2^52 shift rounds to nearest without -ffast-math
	const double q = (x * 0.63661977236758134308 + 6755399441055744.0) - 6755399441055744.0;
	const int quad = static_cast<int>(q);
#else
	const double q = std::nearbyint(x * 0.63661977236758134308);
	const long long quad = static_cast<long long>(q);
#endif
	const double r = ((x - q * 1.57079632673412561417e+00) - q * 6.07710050630396597660e-11) - q * 2.02226624871116645580e-21;
	const double z = r * r;
	const double sr = r + r * z * (-1.66666666666666324348e-01 + z * (8.33333333332248946124e-03 + z * (-1.98412698298579493134e-04 +
		z * (2.75573137070700676789e-06 + z * (-2.50507602534068634195e-08 + z * 1.58969099521155010221e-10)))));
	const double hz = 0.5 * z, w = 1.0 - hz;
	const double cr = w + (((1.0 - w) - hz) + z * z * (4.16666666666666019037e-02 + z * (-1.38888888888741095749e-03 +
		z * (2.48015872894767294178e-05 + z * (-2.75573143513906633035e-07 + z * (2.08757232129817482790e-09 +
		z * -1.13596475577881948265e-11))))));
	switch (quad & 3) {
	case 0: *s = sr; *c = cr; break;
	case 1: *s = cr; *c = -sr; break;
	case 2: *s = -sr; *c = -cr; break;
	default: *s = -cr; *c = sr; break;
	}
}

#if defined(__aarch64__) || defined(_M_ARM64)
// sincos_shared on two NEON lanes: the same reduction and kernels with fused multiply-adds, the quadrant switch as a
// lane select plus sign bits. Branch-free and for |x| < 1e5 only; sincos_shared2/4 test the range.
inline void sincos_kernel2(const float64x2_t x, float64x2_t* s, float64x2_t* c) {
	const float64x2_t q = vrndnq_f64(vmulq_f64(x, vdupq_n_f64(0.63661977236758134308)));
	float64x2_t r = vfmsq_f64(x, q, vdupq_n_f64(1.57079632673412561417e+00));
	r = vfmsq_f64(r, q, vdupq_n_f64(6.07710050630396597660e-11));
	r = vfmsq_f64(r, q, vdupq_n_f64(2.02226624871116645580e-21));
	const float64x2_t z = vmulq_f64(r, r);
	float64x2_t ps = vfmaq_f64(vdupq_n_f64(-2.50507602534068634195e-08), z, vdupq_n_f64(1.58969099521155010221e-10));
	ps = vfmaq_f64(vdupq_n_f64(2.75573137070700676789e-06), z, ps);
	ps = vfmaq_f64(vdupq_n_f64(-1.98412698298579493134e-04), z, ps);
	ps = vfmaq_f64(vdupq_n_f64(8.33333333332248946124e-03), z, ps);
	ps = vfmaq_f64(vdupq_n_f64(-1.66666666666666324348e-01), z, ps);
	const float64x2_t sr = vfmaq_f64(r, vmulq_f64(r, z), ps);
	float64x2_t pc = vfmaq_f64(vdupq_n_f64(2.08757232129817482790e-09), z, vdupq_n_f64(-1.13596475577881948265e-11));
	pc = vfmaq_f64(vdupq_n_f64(-2.75573143513906633035e-07), z, pc);
	pc = vfmaq_f64(vdupq_n_f64(2.48015872894767294178e-05), z, pc);
	pc = vfmaq_f64(vdupq_n_f64(-1.38888888888741095749e-03), z, pc);
	pc = vfmaq_f64(vdupq_n_f64(4.16666666666666019037e-02), z, pc);
	const float64x2_t hz = vmulq_f64(vdupq_n_f64(0.5), z);
	const float64x2_t w = vsubq_f64(vdupq_n_f64(1.0), hz);
	const float64x2_t cr = vaddq_f64(w, vfmaq_f64(vsubq_f64(vsubq_f64(vdupq_n_f64(1.0), w), hz), vmulq_f64(z, z), pc));
	//quadrant q & 3: odd swaps sin and cos, q & 2 negates sin, (q + 1) & 2 negates cos
	const int64x2_t qi = vcvtq_s64_f64(q);
	const uint64x2_t odd = vtstq_s64(qi, vdupq_n_s64(1));
	const uint64x2_t sgn_s = vshlq_n_u64(vreinterpretq_u64_s64(vandq_s64(qi, vdupq_n_s64(2))), 62);
	const uint64x2_t sgn_c = vshlq_n_u64(vreinterpretq_u64_s64(vandq_s64(vaddq_s64(qi, vdupq_n_s64(1)), vdupq_n_s64(2))), 62);
	*s = vreinterpretq_f64_u64(veorq_u64(vreinterpretq_u64_f64(vbslq_f64(odd, cr, sr)), sgn_s));
	*c = vreinterpretq_f64_u64(veorq_u64(vreinterpretq_u64_f64(vbslq_f64(odd, sr, cr)), sgn_c));
}

// A pair with a lane at |x| >= 1e5 or non-finite goes lane by lane through sincos_shared.
inline void sincos_shared2(const float64x2_t x, float64x2_t* s, float64x2_t* c) {
	const uint64x2_t in_range = vcaltq_f64(x, vdupq_n_f64(1e5));
	if ((vgetq_lane_u64(in_range, 0) & vgetq_lane_u64(in_range, 1)) == 0) {
		double xs[2], ss[2], cs[2];
		vst1q_f64(xs, x);
		sincos_shared(xs[0], ss, cs);
		sincos_shared(xs[1], ss + 1, cs + 1);
		*s = vld1q_f64(ss); *c = vld1q_f64(cs);
		return;
	}
	sincos_kernel2(x, s, c);
}

// Two pairs, lane for lane the results of sincos_shared2. One range test puts both kernels in one basic block: each
// kernel is a long dependency chain, and the A72 (Pi 4) only overlaps two of them when they sit side by side, 1.44x.
inline void sincos_shared4(const float64x2_t x0, const float64x2_t x1, float64x2_t* s0, float64x2_t* c0, float64x2_t* s1, float64x2_t* c1) {
	const uint64x2_t in_range = vandq_u64(vcaltq_f64(x0, vdupq_n_f64(1e5)), vcaltq_f64(x1, vdupq_n_f64(1e5)));
	if ((vgetq_lane_u64(in_range, 0) & vgetq_lane_u64(in_range, 1)) == 0) {
		sincos_shared2(x0, s0, c0);
		sincos_shared2(x1, s1, c1);
		return;
	}
	sincos_kernel2(x0, s0, c0);
	sincos_kernel2(x1, s1, c1);
}

// sincos_shared4/2 over n angles, an odd last one through sincos_shared
inline void sincos_shared_n(const int n, const double* x, double* s, double* c) {
	int p = 0;
	for (; p + 3 < n; p += 4) {
		float64x2_t s0, c0, s1, c1;
		sincos_shared4(vld1q_f64(x + p), vld1q_f64(x + p + 2), &s0, &c0, &s1, &c1);
		vst1q_f64(s + p, s0); vst1q_f64(s + p + 2, s1);
		vst1q_f64(c + p, c0); vst1q_f64(c + p + 2, c1);
	}
	for (; p + 1 < n; p += 2) {
		float64x2_t sv, cv;
		sincos_shared2(vld1q_f64(x + p), &sv, &cv);
		vst1q_f64(s + p, sv);
		vst1q_f64(c + p, cv);
	}
	if (p < n) sincos_shared(x[p], s + p, c + p);
}
#else
// No NEON double lanes (armv7, or an OpenBLAS build off ARM64, which XCW's sincos needs): angle by angle
inline void sincos_shared_n(const int n, const double* x, double* s, double* c) {
	for (int p = 0; p < n; p++) sincos_shared(x[p], s + p, c + p);
}
#endif

// exp by fdlibm's e_exp: k = round(x / ln2) by the 1.5*2^52 shift, Cody-Waite reduction, the rational kernel, 2^k into
// the exponent bits; within 1 ulp of the CRT. Selects are bit masks, so a #pragma omp simd loop vectorises it without
// -fno-trapping-math. x below -708.39 gives 0; for x <= 0 only, 2^k overflows above 709.
// MSVC takes the CRT: /fp:fast folds (x/ln2 + 1.5*2^52) - 1.5*2^52 back to x/ln2 (MBIS diverged), and with float_control
// precise the kernel was 15 % slower than its vectorised CRT exp, which these loops get anyway (sucrose MBIS 0.98 -> 0.82 s).
#ifdef _MSC_VER
inline double exp_kernel(const double x) { return std::exp(x); }
#else
inline double exp_kernel(const double x) {
	const uint64_t keep = x < -708.39 ? 0 : ~uint64_t(0);
	const double xc = std::bit_cast<double>((std::bit_cast<uint64_t>(x) & keep) | (std::bit_cast<uint64_t>(-708.39) & ~keep));
	const double kk = xc * 1.44269504088896338700e+00 + 6755399441055744.0;
	const double k = kk - 6755399441055744.0;
	const double hi = xc - k * 6.93147180369123816490e-01, lo = k * 1.90821492927058770002e-10;
	const double r = hi - lo, t = r * r;
	const double c = r - t * (1.66666666666666019037e-01 + t * (-2.77777777770155933842e-03 + t * (6.61375632143793436117e-05 +
		t * (-1.65339022054652515390e-06 + t * 4.13813679705723846039e-08))));
	const double y = 1.0 - ((lo - (r * c) / (2.0 - c)) - hi);
	return y * std::bit_cast<double>(((std::bit_cast<uint64_t>(kk) + 1023) << 52) & keep);
}
#endif
