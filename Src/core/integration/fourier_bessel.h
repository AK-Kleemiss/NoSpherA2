#pragma once
#include "convenience.h"

//Fourier-Bessel transform of an atom-centred shell grid. e^{ik.r} = 4pi sum_lm i^l j_l(kr) Y_lm(k^) Y_lm(r^) turns
//f(k) = sum_p W_p e^{ik.d_p} into 4pi sum_lm i^l Y_lm(k^) sum_s j_l(k r_s) c_s,lm with c_s,lm = sum_{p in s} W_p Y_lm(u_p).
//Cut at L = D/2 for the largest shell's Lebedev degree D, the projection is the exact quadrature of the band-limited
//part and the cut drops what the shells alias (grid lab 01: 3-40x closer than the direct sum to the converged f).
//Real Y_lm are kept by parity of l, entry l(l+1)/2 + m of the even or the odd block: i^l puts even l in Re, odd in Im.
inline int fb_block(const int top) { return top < 0 ? 0 : (top + 1) * (top + 2) / 2; } //entries for l <= top, one parity
inline int fb_n_even(const int L) { return fb_block(L - (L & 1)); }
inline int fb_n_odd(const int L) { return fb_block((L & 1) ? L : L - 1); }

//Real orthonormal Y_lm by the normalised associated-Legendre recurrence; (x+iy)^m carries the sin^m theta and phi.
struct fb_ylm
{
	int L;
	vec kmm, emm, a, b; //norm of P_m^m (sqrt 2 for m > 0), sqrt(2m+3), a_lm and 1/a_{l-1,m} at l(l+1)/2 + m
	explicit fb_ylm(const int L_) : L(L_), kmm(L_ + 1), emm(L_ + 1), a((L_ + 1) * (L_ + 2) / 2, 0.0), b(a)
	{
		kmm[0] = 1.0 / std::sqrt(constants::FOUR_PI);
		for (int m = 1; m <= L; m++)
			kmm[m] = kmm[m - 1] * std::sqrt((2.0 * m + 1.0) / (2.0 * m)) * (m == 1 ? constants::SQRT2 : 1.0);
		for (int m = 0; m <= L; m++)
			emm[m] = std::sqrt(2.0 * m + 3.0);
		for (int l = 2; l <= L; l++)
			for (int m = 0; m <= l - 2; m++)
			{
				a[l * (l + 1) / 2 + m] = std::sqrt((4.0 * l * l - 1.0) / (double(l) * l - double(m) * m));
				b[l * (l + 1) / 2 + m] = std::sqrt((double(l - 1) * (l - 1) - double(m) * m) / (4.0 * (l - 1) * (l - 1) - 1.0));
			}
	}
	//ye/yo: fb_n_even(top)/fb_n_odd(top) entries for l <= top <= L; the offsets do not depend on top
	void eval(const double x, const double y, const double z, double* ye, double* yo, int top = -1) const
	{
		if (top < 0)
			top = L;
		double c = 1.0, s = 0.0;
		for (int m = 0; m <= top; m++)
		{
			if (m > 0)
			{
				const double cn = c * x - s * y;
				s = s * x + c * y;
				c = cn;
			}
			const auto put = [&](const int l, const double p) {
				double* Y = (l & 1) ? yo : ye;
				const int o = l * (l + 1) / 2;
				if (m == 0)
					Y[o] = p;
				else
				{
					Y[o + m] = p * c;
					Y[o - m] = p * s;
				}
			};
			double p2 = kmm[m];
			put(m, p2);
			if (m == top)
				break;
			double p1 = emm[m] * z * p2;
			put(m + 1, p1);
			for (int l = m + 2; l <= top; l++)
			{
				const int i = l * (l + 1) / 2 + m;
				const double p = a[i] * (z * p1 - b[i] * p2);
				p2 = p1;
				p1 = p;
				put(l, p);
			}
		}
	}
};

//j_l(x), l = 0..L. std::sph_bessel is missing from libc++. Upward recurrence where it is stable (x > L), Miller's
//downward recurrence below, normalised on the larger of j0, j1.
inline void fb_sph_bessel(const int L, const double x, double* j)
{
	if (x < 1e-6)
	{
		double t = 1.0;
		j[0] = 1.0 - x * x / 6.0;
		for (int l = 1; l <= L; l++)
			j[l] = (t *= x / (2.0 * l + 1.0)); //next term is x^2/(4l+6) relative
		return;
	}
	const double sx = std::sin(x), cx = std::cos(x), j0 = sx / x, j1 = (j0 - cx) / x;
	if (x > L)
	{
		j[0] = j0;
		if (L > 0)
			j[1] = j1;
		for (int l = 1; l < L; l++)
			j[l + 1] = (2.0 * l + 1.0) / x * j[l] - j[l - 1];
		return;
	}
	const int ls = L + 20 + static_cast<int>(std::sqrt(10.0 * L));
	double fp = 0.0, f = 1e-30;
	for (int l = ls; l >= 1; l--)
	{
		const double fm = (2.0 * l + 1.0) / x * f - fp;
		fp = f;
		f = fm;
		if (l - 1 <= L)
			j[l - 1] = f;
		if (std::abs(f) > 1e200)
		{
			f *= 1e-200;
			fp *= 1e-200;
			for (int t = std::max(l - 1, 0); t <= L; t++)
				j[t] *= 1e-200;
		}
	}
	const double sc = std::abs(j0) > std::abs(j1) ? j0 / j[0] : j1 / j[1];
	for (int l = 0; l <= L; l++)
		j[l] *= sc;
}

//Lebedev degree of the smallest rule holding n points
inline int fb_lebedev_degree(const long long n)
{
	for (int i = 0; i < constants::max_LT - 1; i++)
		if (constants::lebedev_table[i] >= n)
			return i < 15 ? 2 * i + 3 : 35 + 6 * (i - 15);
	return 131;
}
