#pragma once
#include <limits>

//constants::sqrt and the real spherical harmonics of one l contracted with their 2l+1
//coefficients, x, y, z normalised; the l <= 8 switch of constants::spherical_harmonic(l, d, coefs)
//in constants.cpp, which forwards here. Templated on the number type so aux_density.h
//differentiates it with dual, and free of project includes so aux_density_gpu.cu compiles it
//under nvcc and hipcc: those see the host-device attribute, every other compiler plain inline
//functions. The prefactors are constants::sqrt, the same constexpr values as c_* in constants.h.
#if defined(__CUDACC__) || defined(__HIPCC__)
#define AUX_HD __host__ __device__
#else
#define AUX_HD
#endif

namespace constants
{
	constexpr double PI = 3.1415926535897932384626433832795028;
	AUX_HD double constexpr sqrtNewtonRaphson(double x, double curr, double prev)
	{
		return curr == prev
			? curr
			: sqrtNewtonRaphson(x, 0.5 * (curr + x / curr), curr);
	}

	/*
	 * Constexpr version of the square root
	 * Return value:
	 *   - For a finite and non-negative value of "x", returns an approximation for the square root of "x"
	 *   - Otherwise, returns NaN
	 * Taken from https://stackoverflow.com/questions/8622256/in-c11-is-sqrt-defined-as-constexpr
	 */
	AUX_HD double constexpr sqrt(double x)
	{
		return x >= 0 && x < std::numeric_limits<double>::infinity()
			? sqrtNewtonRaphson(x, x, 0)
			: std::numeric_limits<double>::quiet_NaN();
	}
	//sqrt(n / (d pi)) as c_n_dp, the prefactors of the real spherical harmonics; constexpr
	//so the Newton-Raphson sqrt runs at compile time, as a plain call inside the expressions
	//below MSVC evaluated it per point, 10 % of a sucrose ELI run (VTune, 19 Sep 2026)
	constexpr double c_1_4p = sqrt(1.0 / (4.0 * PI));
	constexpr double c_3_4p = sqrt(3.0 / (4.0 * PI));
	constexpr double c_5_4p = sqrt(5.0 / (4.0 * PI));
	constexpr double c_5_16p = sqrt(5.0 / (16.0 * PI));
	constexpr double c_7_4p = sqrt(7.0 / (4.0 * PI));
	constexpr double c_7_16p = sqrt(7.0 / (16.0 * PI));
	constexpr double c_9_4p = sqrt(9.0 / (4.0 * PI));
	constexpr double c_9_256p = sqrt(9.0 / (256.0 * PI));
	constexpr double c_11_4p = sqrt(11.0 / (4.0 * PI));
	constexpr double c_11_256p = sqrt(11.0 / (256.0 * PI));
	constexpr double c_13_p = sqrt(13.0 / PI);
	constexpr double c_13_4p = sqrt(13.0 / (4.0 * PI));
	constexpr double c_13_1024p = sqrt(13.0 / (1024.0 * PI));
	constexpr double c_15_p = sqrt(15.0 / PI);
	constexpr double c_15_4p = sqrt(15.0 / (4.0 * PI));
	constexpr double c_15_16p = sqrt(15.0 / (16.0 * PI));
	constexpr double c_17_p = sqrt(17.0 / PI);
	constexpr double c_17_4p = sqrt(17.0 / (4.0 * PI));
	constexpr double c_19_4p = sqrt(19.0 / (4.0 * PI));
	constexpr double c_21_32p = sqrt(21.0 / (32.0 * PI));
	constexpr double c_35_p = sqrt(35.0 / PI);
	constexpr double c_35_2p = sqrt(35.0 / (2 * PI));
	constexpr double c_35_32p = sqrt(35.0 / (32.0 * PI));
	constexpr double c_45_16p = sqrt(45.0 / (16.0 * PI));
	constexpr double c_45_32p = sqrt(45.0 / (32.0 * PI));
	constexpr double c_45_64p = sqrt(45.0 / (64.0 * PI));
	constexpr double c_91_p = sqrt(91.0 / (PI));
	constexpr double c_105_p = sqrt(105.0 / PI);
	constexpr double c_105_4p = sqrt(105.0 / (4.0 * PI));
	constexpr double c_105_16p = sqrt(105.0 / (16.0 * PI));
	constexpr double c_165_256p = sqrt(165.0 / (256.0 * PI));
	constexpr double c_273_p = sqrt(273.0 / (PI));
	constexpr double c_273_256p = sqrt(273.0 / (256.0 * PI));
	constexpr double c_315_16p = sqrt(315.0 / (16.0 * PI));
	constexpr double c_315_32p = sqrt(315.0 / (32.0 * PI));
	constexpr double c_315_256p = sqrt(315.0 / (256.0 * PI));
	constexpr double c_385_p = sqrt(385.0 / PI);
	constexpr double c_385_512p = sqrt(385.0 / (512.0 * PI));
	constexpr double c_595_2p = sqrt(595.0 / (2 * PI));
	constexpr double c_693_2048p = sqrt(693.0 / (2048.0 * PI));
	constexpr double c_715_p = sqrt(715.0 / PI);
	constexpr double c_1001_2p = sqrt(1001.0 / (2 * PI));
	constexpr double c_1155_64p = sqrt(1155.0 / (64.0 * PI));
	constexpr double c_1309_p = sqrt(1309.0 / PI);
	constexpr double c_1365_2p = sqrt(1365.0 / (2 * PI));
	constexpr double c_3003_2p = sqrt(3003.0 / (2 * PI));
	constexpr double c_3003_2048p = sqrt(3003.0 / (2048.0 * PI));
	constexpr double c_3465_256p = sqrt(3465.0 / (256.0 * PI));
	constexpr double c_5005_2p = sqrt(5005.0 / (2 * PI));
	constexpr double c_7293_2p = sqrt(7293.0 / (2 * PI));
	constexpr double c_12155_p = sqrt(12155.0 / PI);
	constexpr double c_17017_p = sqrt(17017.0 / PI);
	constexpr double c_19635_p = sqrt(19635.0 / PI);
	template <class T>
	AUX_HD inline T spherical_harmonic(const int l, const T x, const T y, const T z, const double* coefs)
	{
		switch (l)
		{
		case 0:
			return c_1_4p * coefs[0];
		case 1:
			return c_3_4p * (coefs[0] * y + coefs[1] * z + coefs[2] * x);
		case 2:
			return c_15_4p * (y * x * coefs[0] + y * z * coefs[1] + x * z * coefs[3]) + c_5_16p * (3 * z * z - 1.0) * coefs[2] + c_15_16p * (x * x - y * y) * coefs[4];
		case 3:
		{
			const T y2 = y * y, x2 = x * x, z2 = z * z;
			return c_35_32p * (y * (3 * x2 - y2) * coefs[0] + x * (x2 - 3 * y2) * coefs[6]) +
				c_105_4p * x * y * z * coefs[1] +
				c_21_32p * (y * (5 * z2 - 1.0) * coefs[2] + x * (5 * z2 - 1.0) * coefs[4]) +
				c_7_16p * (5 * z2 * z - 3 * z) * coefs[3] +
				c_105_16p * ((x2 - y2) * z) * coefs[5];
		}
		case 4:
		{
			const T x2 = x * x, y2 = y * y, z2 = z * z;
			return c_315_16p * x * y * (x2 - y2) * coefs[0] +
				c_315_32p * (y * (3 * x2 - y2) * z * coefs[1] + x * (x2 - 3 * y2) * z * coefs[7]) +
				c_45_16p * x * y * (7 * z2 - 1.0) * coefs[2] +
				c_45_32p * (y * (7 * z2 * z - 3 * z) * coefs[3] + x * (7 * z2 * z - 3 * z) * coefs[5]) +
				c_9_256p * (35 * z2 * z2 - 30 * z2 + 3.0) * coefs[4] +
				c_45_64p * (x2 - y2) * (7 * z2 - 1.0) * coefs[6] +
				c_315_256p * ((x2 * (x2 - 3 * y2)) - (y2 * (3 * x2 - y2))) * coefs[8];
		}
		case 5:
		{
			const T x2 = x * x, y2 = y * y, z2 = z * z;
			return c_693_2048p * (2 * y2 * y2 * y - 20 * x2 * y2 * y + 10 * y * x2 * x2) * coefs[0] +
				-c_3465_256p * z * ((4 * x * y2 * y - 4 * x2 * x * y) * coefs[1] - (x2 * x2 - 6 * x2 * y2 + y2 * y2) * coefs[9]) +
				c_385_512p * (9 * z2 - 1.0) * (y * (3 * x2 - y2) * coefs[2] + x * (x2 - 3 * y2) * coefs[8]) +
				c_1155_64p * (3 * z2 * z - z) * (2 * x * y * coefs[3] + (x2 - y2) * coefs[7]) +
				c_165_256p * (21 * z2 * z2 - 14 * z2 + 1.0) * (y * coefs[4] + x * coefs[6]) +
				c_11_256p * (63 * z2 * z2 * z - 70 * z2 * z + 15 * z) * coefs[5] +
				c_693_2048p * (2 * x2 * x2 * x - 20 * x2 * x * y2 + 10 * x * y2 * y2) * coefs[10];
		}
		case 6:
		{
			const T x2 = x * x, y2 = y * y, z2 = z * z;
			const T x4 = x2 * x2, y4 = y2 * y2, z4 = z2 * z2;
			return (1.0 / 32.0) * c_13_p * (z2 * (z2 * (231.0 * z2 - 315.0) + 105.0) - 5.0) * coefs[6] +
				(1.0 / 16.0) * c_273_p * (33 * z4 - 30 * z2 + 5) * (x * z * coefs[7] + y * z * coefs[5]) +
				(1.0 / 32.0) * c_1365_2p * (33 * z4 - 18 * z2 + 1) * ((x2 - y2) * coefs[8] + 2 * x * y * coefs[4]) +
				(1.0 / 16.0) * c_1365_2p * z * (11 * z2 - 3) * (x * (x2 - 3 * y2) * coefs[9] + y * (3 * x2 - y2) * coefs[3]) +
				(3.0 / 32.0) * c_91_p * (11 * z2 - 1) * ((x2 * (x2 - 6 * y2) + y4) * coefs[10] + x * y * 4 * (x2 - y2) * coefs[2]) +
				(3.0 / 16.0) * c_1001_2p * z * ((x4 * x - 10 * x2 * x * y2 + 5 * x * y4) * coefs[11] + (5 * x4 * y - 10 * x2 * y2 * y + y4 * y) * coefs[1]) +
				(1.0 / 32.0) * c_3003_2p * (-y4 * y2 + 15 * y4 * x2 - 15 * y2 * x4 + x4 * x2) * coefs[12] +
				(1.0 / 32.0) * c_3003_2p * (6 * x4 * x * y - 20 * x2 * x * y2 * y + 6 * x * y4 * y) * coefs[0];
		}
		case 7:
		{
			const T x2 = x * x, y2 = y * y, z2 = z * z;
			const T x4 = x2 * x2, y4 = y2 * y2;
			return (1.0 / 32.0) * c_15_p * z * (z2 * (z2 * (429.0 * z2 - 693.0) + 315.0) - 35.0) * coefs[7] +
				(1.0 / 64.0) * c_105_p * (z2 * (z2 * (429.0 * z2 - 495.0) + 135.0) - 5.0) * (x * coefs[8] + y * coefs[6]) +
				(3.0 / 32.0) * c_35_2p * z * (z2 * (143.0 * z2 - 110.0) + 15.0) * ((x2 - y2) * coefs[9] + 2 * x * y * coefs[5]) +
				(3.0 / 64.0) * c_35_p * (z2 * (143.0 * z2 - 66.0) + 3.0) * (x * (x2 - 3 * y2) * coefs[10] + y * (3 * x2 - y2) * coefs[4]) +
				(3.0 / 32.0) * c_385_p * z * (13.0 * z2 - 3.0) * ((x2 * (x2 - 6 * y2) + y4) * coefs[11] + x * y * 4 * (x2 - y2) * coefs[3]) +
				(3.0 / 64.0) * c_385_p * (13.0 * z2 - 1.0) * (x * (x4 + 5 * y2 * (-2 * x2 + y2)) * coefs[12] + y * (5 * x2 * (x2 - 2 * y2) + y4) * coefs[2]) +
				(3.0 / 32.0) * c_5005_2p * z * ((y2 * (y2 * (-y2 + 15 * x2) - 15 * x4) + x4 * x2) * coefs[13] + x * y * (x2 * (6 * x2 - 20 * y2) + 6 * y4) * coefs[1]) +
				(3.0 / 64.0) * c_715_p * x * (x4 * x2 + 7 * y2 * (-3 * x4 + y2 * (5 * x2 - y2))) * coefs[14] +
				(3.0 / 64.0) * c_715_p * y * (x2 * (x2 * (7 * x2 - 35 * y2) + 21 * y4) - y4 * y2) * coefs[0];
		}
		case 8:
		{
			const T x2 = x * x, y2 = y * y, z2 = z * z;
			const T x4 = x2 * x2, y4 = y2 * y2, z4 = z2 * z2;
			return (1.0 / 256.0) * c_17_p * (z2 * (z2 * (z2 * (6435.0 * z2 - 12012.0) + 6930.0) - 1260.0) + 35.0) * coefs[8] +
				(3.0 / 64.0) * c_17_p * z * (715.0 * z4 * z2 - 1001.0 * z4 + 385.0 * z2 - 35.0) * (x * coefs[9] + y * coefs[7]) +
				(3.0 / 64.0) * c_595_2p * (143 * z4 * z2 - 143 * z4 + 33 * z2 - 1) * ((x2 - y2) * coefs[10] + 2 * x * y * coefs[6]) +
				(1.0 / 64.0) * c_19635_p * z * (39 * z4 - 26 * z2 + 3) * ((x2 * x - 3 * x * y2) * coefs[11] + (3 * x2 * y - y2 * y) * coefs[5]) +
				(3.0 / 128.0) * c_1309_p * (65 * z4 - 26 * z2 + 1) * ((x2 * (x2 - 6 * y2) + y4) * coefs[12] + x * y * 4 * (x2 - y2) * coefs[4]) +
				(3.0 / 64.0) * c_17017_p * z * (5 * z2 - 1) * (x * (x4 + 5 * y2 * (-2 * x2 + y2)) * coefs[13] + y * (5 * x2 * (x2 - 2 * y2) + y4) * coefs[3]) +
				(1.0 / 64.0) * c_7293_2p * (15 * z2 - 1) * ((y2 * (y2 * (-y2 + 15 * x2) - 15 * x4) + x4 * x2) * coefs[14] + x * y * (x2 * (6 * x2 - 20 * y2) + 6 * y4) * coefs[2]) +
				(3.0 / 64.0) * c_12155_p * z * (x * (x4 * x2 + 7 * y2 * (-3 * x4 + y2 * (5 * x2 - y2))) * coefs[15] + y * (x2 * (x2 * (7 * x2 - 35 * y2) + 21 * y4) - y4 * y2) * coefs[1]) +
				(3.0 / 256.0) * c_12155_p * (y2 * (y2 * (y4 - 28 * x2 * y2 + 70 * x4) - 28 * x4 * x2) + x4 * x4) * coefs[16] +
				(3.0 / 256.0) * c_12155_p * y * x * (y2 * (y2 * (-8 * y2 + 56 * x2) - 56 * x4) + 8 * x4 * x2) * coefs[0];
		}
		}
		return T(0.0);
	}
}
