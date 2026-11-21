#include "pch.h"
#include "convenience.h"
#include "constants.h"

bool generate_cart2sph_mat(vec2 &d, vec2 &f, vec2 &g, vec2 &h)
{
	//
	// From 5D: D 0, D + 1, D - 1, D + 2, D - 2
	// To 6D : 1  2  3  4  5  6
	// XX, YY, ZZ, XY, XZ, YZ
	//
	d.resize(6);
	for (int i = 0; i < 6; i++)
	{
		d[i].resize(5, 0.0);
	}
	// D0 = -0.5 * XX - 0.5 * YY + ZZ
	d[0][0] = -0.5;
	d[1][0] = -0.5;
	d[2][0] = 1.0;
	// D + 1 = XZ
	d[4][1] = 1.0;
	// D - 1 = YZ
	d[5][2] = 1.0;
	// D + 2 = SQRT(3) / 2 * (XX - YY)
	d[0][3] = sqrt(3.0) / 2.0;
	d[1][3] = -sqrt(3.0) / 2.0;
	// D - 2 = XY
	d[3][4] = 1.0;

	// From 7F: F 0, F + 1, F - 1, F + 2, F - 2, F + 3, F - 3
	// To 10F : 1   2   3   4   5   6   7   8   9  10
	// XXX, YYY, ZZZ, XYY, XXY, XXZ, XZZ, YZZ, YYZ, XYZ(Gaussian sequence, not identical to Multiwfn)
	//
	f.resize(10);
	for (int i = 0; i < 10; i++)
	{
		f[i].resize(7, 0.0);
	}
	// F 0 = -3 / (2 * sqrt5) * (XXZ + YYZ) + ZZZ
	f[2][0] = 1.0;
	f[5][0] = -1.5 / sqrt(5.0);
	f[8][0] = -1.5 / sqrt(5.0);
	// F + 1 = -sqrt(3 / 8) * XXX - sqrt(3 / 40) * XYY + sqrt(6 / 5) * XZZ
	f[0][1] = -sqrt(3.0 / 8.0);
	f[3][1] = -sqrt(3.0 / 40.0);
	f[6][1] = sqrt(6.0 / 5.0);
	// F - 1 = -sqrt(3 / 40) * XXY - sqrt(3 / 8) * YYY + sqrt(6 / 5) * YZZ
	f[1][2] = -sqrt(3.0 / 8.0);
	f[4][2] = -sqrt(3.0 / 40.0);
	f[7][2] = sqrt(6.0 / 5.0);
	// F + 2 = sqrt3 / 2 * (XXZ - YYZ)
	f[5][3] = sqrt(3.0) / 2.0;
	f[8][3] = -sqrt(3.0) / 2.0;
	// F - 2 = XYZ
	f[9][4] = 1.0;
	// F + 3 = sqrt(5 / 8) * XXX - 3 / sqrt8 * XYY
	f[0][5] = sqrt(5.0 / 8.0);
	f[3][5] = -3.0 / sqrt(8.0);
	// F - 3 = 3 / sqrt8 * XXY - sqrt(5 / 8) * YYY
	f[1][6] = -sqrt(5.0 / 8.0);
	f[4][6] = 3.0 / sqrt(8.0);

	// From 9G: G 0, G + 1, G - 1, G + 2, G - 2, G + 3, G - 3, G + 4, G - 4
	// To 15G : 1    2    3    4    5    6    7    8
	// ZZZZ, YZZZ, YYZZ, YYYZ, YYYY, XZZZ, XYZZ, XYYZ
	// 9   10   11   12   13   14   15
	// XYYY, XXZZ, XXYZ, XXYY, XXXZ, XXXY, XXXX
	//
	g.resize(15);
	for (int i = 0; i < 15; i++)
	{
		g[i].resize(9, 0.0);
	}
	// G 0 = ZZZZ + 3 / 8 * (XXXX + YYYY) - 3 * sqrt(3 / 35) * (XXZZ + YYZZ - 1 / 4 * XXYY)
	g[0][0] = 1.0;
	g[2][0] = -3.0 * sqrt(3.0 / 35.0);
	g[4][0] = 3.0 / 8.0;
	g[9][0] = -3.0 * sqrt(3.0 / 35.0);
	g[11][0] = 3.0 / 4.0 * sqrt(3.0 / 35.0);
	g[14][0] = 3.0 / 8.0;
	// G + 1 = 2 * sqrt(5 / 14) * XZZZ - 3 / 2 * sqrt(5 / 14) * XXXZ - 3 / 2 / sqrt14 * XYYZ
	g[5][1] = 2.0 * sqrt(5.0 / 14.0);
	g[7][1] = -1.5 / sqrt(14.0);
	g[12][1] = -1.5 * sqrt(5.0 / 14.0);
	// G - 1 = 2 * sqrt(5 / 14) * YZZZ - 3 / 2 * sqrt(5 / 14) * YYYZ - 3 / 2 / sqrt14 * XXYZ
	g[1][2] = 2.0 * sqrt(5.0 / 14.0);
	g[3][2] = -1.5 * sqrt(5.0 / 14.0);
	g[10][2] = -1.5 / sqrt(14.0);
	// G + 2 = 3 * sqrt(3 / 28) * (XXZZ - YYZZ) - sqrt5 / 4 * (XXXX - YYYY)
	g[2][3] = -3.0 * sqrt(3.0 / 28.0);
	g[4][3] = sqrt(5.0) / 4.0;
	g[9][3] = 3.0 * sqrt(3.0 / 28.0);
	g[14][3] = -sqrt(5.0) / 4.0;
	// G - 2 = 3 / sqrt7 * XYZZ - sqrt(5 / 28) * (XXXY + XYYY)
	g[6][4] = 3.0 / sqrt(7.0);
	g[8][4] = -sqrt(5.0 / 28.0);
	g[13][4] = -sqrt(5.0 / 28.0);
	// G + 3 = sqrt(5 / 8) * XXXZ - 3 / sqrt8 * XYYZ
	g[7][5] = -3.0 / sqrt(8.0);
	g[12][5] = sqrt(5.0 / 8.0);
	// G - 3 = -sqrt(5 / 8) * YYYZ + 3 / sqrt8 * XXYZ
	g[3][6] = -sqrt(5.0 / 8.0);
	g[10][6] = 3.0 / sqrt(8.0);
	// G + 4 = sqrt35 / 8 * (XXXX + YYYY) - 3 / 4 * sqrt3 * XXYY
	g[4][7] = sqrt(35.0) / 8.0;
	g[11][7] = -3.0 / 4.0 * sqrt(3.0);
	g[14][7] = sqrt(35.0) / 8.0;
	// G - 4 = sqrt5 / 2 * (XXXY - XYYY)
	g[8][8] = -sqrt(5.0) / 2.0;
	g[13][8] = sqrt(5.0) / 2.0;

	// From 11H: H 0, H + 1, H - 1, H + 2, H - 2, H + 3, H - 3, H + 4, H - 4, H + 5, H - 5
	// To 21H : 1     2     3     4     5     6     7     8     9    10
	// ZZZZZ YZZZZ YYZZZ YYYZZ YYYYZ YYYYY XZZZZ XYZZZ XYYZZ XYYYZ
	// 11    12    13    14    15    16    17    18    19    20    21
	// XYYYY XXZZZ XXYZZ XXYYZ XXYYY XXXZZ XXXYZ XXXYY XXXXZ XXXXY XXXXX
	//
	h.resize(21);
	for (int i = 0; i < 21; i++)
	{
		h[i].resize(11);
		std::fill(h[i].begin(), h[i].end(), 0.0);
	}
	// H 0 = ZZZZZ - 5 / sqrt21 * (XXZZZ + YYZZZ) + 5 / 8 * (XXXXZ + YYYYZ) + sqrt(15 / 7) / 4 * XXYYZ
	h[0][0] = 1.0;
	h[11][0] = -5.0 / sqrt(21.0);
	h[2][0] = -5.0 / sqrt(21.0);
	h[18][0] = 5.0 / 8.0;
	h[4][0] = 5.0 / 8.0;
	h[13][0] = sqrt(15.0 / 7.0) / 4.0;
	// H + 1 = sqrt(5 / 3) * XZZZZ - 3 * sqrt(5 / 28) * XXXZZ - 3 / sqrt28 * XYYZZ + sqrt15 / 8 * XXXXX + sqrt(5 / 3) / 8 * XYYYY + sqrt(5 / 7) / 4 * XXXYY
	h[6][1] = sqrt(5.0 / 3.0);
	h[15][1] = -3.0 * sqrt(5.0 / 28.0);
	h[8][1] = -3.0 / sqrt(28.0);
	h[20][1] = sqrt(15.0) / 8.0;
	h[10][1] = sqrt(5.0 / 3.0) / 8.0;
	h[17][1] = sqrt(5.0 / 7.0) / 4.0;
	// H - 1 = sqrt(5 / 3) * YZZZZ - 3 * sqrt(5 / 28) * YYYZZ - 3 / sqrt28 * XXYZZ + sqrt15 / 8 * YYYYY + sqrt(5 / 3) / 8 * XXXXY + sqrt(5 / 7) / 4 * XXYYY
	h[1][2] = sqrt(5.0 / 3.0);
	h[3][2] = -3.0 * sqrt(5.0 / 28.0);
	h[12][2] = -3.0 / sqrt(28.0);
	h[5][2] = sqrt(15.0) / 8.0;
	h[19][2] = sqrt(5.0 / 3.0) / 8.0;
	h[14][2] = sqrt(5.0 / 7.0) / 4.0;
	// H + 2 = sqrt5 / 2 * (XXZZZ - YYZZZ) - sqrt(35 / 3) / 4 * (XXXXZ - YYYYZ)
	h[11][3] = sqrt(5.0) / 2.0;
	h[2][3] = -sqrt(5.0) / 2.0;
	h[18][3] = -sqrt(35.0 / 3.0) / 4.0;
	h[4][3] = sqrt(35.0 / 3.0) / 4.0;
	// H - 2 = sqrt(5 / 3) * XYZZZ - sqrt(5 / 12) * (XXXYZ + XYYYZ)
	h[7][4] = sqrt(5.0 / 3.0);
	h[16][4] = -sqrt(5.0 / 12.0);
	h[9][4] = -sqrt(5.0 / 12.0);
	// H + 3 = sqrt(5 / 6) * XXXZZ - sqrt(3 / 2) * XYYZZ - sqrt(35 / 2) / 8 * (XXXXX - XYYYY) + sqrt(5 / 6) / 4 * XXXYY
	h[15][5] = sqrt(5.0 / 6.0);
	h[8][5] = -sqrt(1.5);
	h[20][5] = -sqrt(17.5) / 8.0;
	h[10][5] = sqrt(17.5) / 8.0;
	h[17][5] = sqrt(5.0 / 6.0) / 4.0;
	// H - 3 = -sqrt(5 / 6) * YYYZZ + sqrt(3 / 2) * XXYZZ - sqrt(35 / 2) / 8 * (XXXXY - YYYYY) - sqrt(5 / 6) / 4 * XXYYY
	h[3][6] = -sqrt(5.0 / 6.0);
	h[12][6] = sqrt(1.5);
	h[19][6] = -sqrt(17.5) / 8.0;
	h[5][6] = sqrt(17.5) / 8.0;
	h[14][6] = -sqrt(5.0 / 6.0) / 4.0;
	// H + 4 = sqrt35 / 8 * (XXXXZ + YYYYZ) - 3 / 4 * sqrt3 * XXYYZ
	h[18][7] = sqrt(35.0) / 8.0;
	h[4][7] = sqrt(35.0) / 8.0;
	h[13][7] = -0.75 * sqrt(3.0);
	// H - 4 = sqrt5 / 2 * (XXXYZ - XYYYZ)
	h[16][8] = sqrt(5.0) / 2.0;
	h[9][8] = -sqrt(5.0) / 2.0;
	// H + 5 = 3 / 8 * sqrt(7 / 2) * XXXXX + 5 / 8 * sqrt(7 / 2) * XYYYY - 5 / 4 * sqrt(3 / 2) * XXXYY
	h[20][9] = 3.0 / 8.0 * sqrt(3.5);
	h[10][9] = 5.0 / 8.0 * sqrt(3.5);
	h[17][9] = -1.25 * sqrt(1.5);
	// H - 5 = 3 / 8 * sqrt(7 / 2) * YYYYY + 5 / 8 * sqrt(7 / 2) * XXXXY - 5 / 4 * sqrt(3 / 2) * XXYYY
	h[5][10] = 3.0 / 8.0 * sqrt(3.5);
	h[19][10] = 5.0 / 8.0 * sqrt(3.5);
	h[14][10] = -1.25 * sqrt(1.5);
	return true;
}


double get_bessel_ratio(const double nu, const double x)
{
	const double RECUR_BIG = DBL_MAX;
	const double RECUR_SMALL = DBL_MIN;
	const int maxiter = 10000;
	int n = 1;
	double Anm2 = 1.0;
	double Bnm2 = 0.0;
	double Anm1 = 0.0;
	double Bnm1 = 1.0;
	double a1 = x / (2.0 * (nu + 1.0));
	double An = Anm1 + a1 * Anm2;
	double Bn = Bnm1 + a1 * Bnm2;
	double an;
	double fn = An / Bn;
	double dn = a1;
	double s = 1.0;

	while (n < maxiter)
	{
		double old_fn;
		double del;
		n++;
		Anm2 = Anm1;
		Bnm2 = Bnm1;
		Anm1 = An;
		Bnm1 = Bn;
		an = -x * x / (4.0 * (nu + n - 1.0) * (nu + n));
		An = Anm1 + an * Anm2;
		Bn = Bnm1 + an * Bnm2;

		if (fabs(An) > RECUR_BIG || fabs(Bn) > RECUR_BIG)
		{
			An /= RECUR_BIG;
			Bn /= RECUR_BIG;
			Anm1 /= RECUR_BIG;
			Bnm1 /= RECUR_BIG;
			Anm2 /= RECUR_BIG;
		}
		else if (fabs(An) < RECUR_SMALL || fabs(Bn) < RECUR_SMALL)
		{
			An /= RECUR_SMALL;
			Bn /= RECUR_SMALL;
			Anm1 /= RECUR_SMALL;
			Bnm1 /= RECUR_SMALL;
			Anm2 /= RECUR_SMALL;
			Bnm2 /= RECUR_SMALL;
		}

		old_fn = fn;
		fn = An / Bn;
		del = old_fn / fn;

		dn = 1.0 / (2.0 * (nu + n) / x - dn);
		if (dn < 0.0)
			s = -s;

		if (fabs(del - 1.0) < 2.0 * DBL_EPSILON)
			break;
	}

	return fn;
}

double bessel_first_kind(const int l, const double x)
{
	if (l < 0 || x < 0.0)
	{
		err_not_impl_f("This is not implemented, pelase dont do this to me!", std::cout);
		return -1000;
	}
	else if (x == 0.0)
	{
		return (l > 0 ? 0.0 : 1.0);
	}
	else if (l == 0)
	{
		return sin(x) / x;
	}
	else if (l == 1)
	{
		return (sin(x) / x - cos(x)) / x;
	}
	else if (l == 2)
	{
		const double f = (3.0 / (x * x) - 1.0);
		return (f * sin(x) - 3.0 * cos(x) / x) / x;
	}
	else if (l == 3)
	{
		double x2 = x * x;
		const double f1 = (x2 - 15.0);
		const double f2 = (6. * x2 - 15.);
		return (-f2 * sin(x) + f1 * cos(x) * x) / pow(x, 4);
	}
	else if (l == 4)
	{
		double x2 = x * x;
		const double f1 = (10. * x2 - 105.0);
		const double f2 = (x2 * x2 - 45. * x2 + 105.);
		return (f2 * sin(x) + f1 * cos(x) * x) / pow(x, 5);
	}
	else if (l == 5)
	{
		double x2 = x * x;
		double x4 = x2 * x2;
		const double f1 = (-x4 + 105.0 * x2 - 945.);
		const double f2 = (15. * x4 - 420. * x2 + 945.);
		return (f2 * sin(x) + f1 * cos(x) * x) / (x4 * x2);
	}
	else if (l == 6)
	{
		double x2 = x * x;
		double x4 = x2 * x2;
		const double f1 = 21. * (x4 - 60.0 * x2 + 495.);
		const double f2 = (-x4 * x2 + 210. * x4 - 4725 * x2 + 10395.);
		return (f2 * sin(x) - f1 * cos(x) * x) / pow(x, 7);
	}
	else
	{
		double ratio = get_bessel_ratio(l + 0.5, x);
		const double smallest = DBL_MIN / DBL_EPSILON;
		double jellp1 = smallest * ratio;
		double jell = smallest;
		double jellm1;
		int ell;
		for (ell = l; ell > 0; ell--)
		{
			jellm1 = -jellp1 + (2 * ell + 1) / x * jell;
			jellp1 = jell;
			jell = jellm1;
		}

		if (fabs(jell) > fabs(jellp1))
		{
			double pre = smallest / jell;
			return bessel_first_kind(0, x) * pre;
		}
		else
		{
			double pre = smallest / jellp1;
			return bessel_first_kind(1, x) * pre;
		}
	}
}


cdouble hypergeometric(double a, double b, double c, cdouble x)
{
	const double TOLERANCE = 1.0e-10;
	cdouble term = a * b * x / c;
	cdouble value = 1.0 + term;
	int n = 1;

	while (std::abs(term) > TOLERANCE)
	{
		a++, b++, c++, n++;
		term *= a * b * x / c / static_cast<double>(n);
		value += term;
	}

	return value;
}

double hypergeometric(double a, double b, double c, double x)
{
	const double TOLERANCE = 1.0e-10;
	double term = a * b * x / c;
	double value = 1.0 + term;
	int n = 1;

	while (std::abs(term) > TOLERANCE)
	{
		a++, b++, c++, n++;
		term *= a * b * x / c / n;
		value += term;
	}

	return value;
}
