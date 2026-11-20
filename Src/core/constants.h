#pragma once

#include "pch.h"
#include "spherical_harmonic.h"

namespace constants
{
	extern double exp_cutoff;
	extern bool hide_gpu_notes;
	extern bool hide_timings; //-no_date: wall-clock lines would make the log non-reproducible
	//Bumped at the start of every run_app. Olex2 and the test suite call run_app repeatedly in one
	//process, where a plain static "said it once" latch stays set and the next run's log loses the line.
	extern unsigned run_id;
	//True on the first call per run for this latch (a static std::atomic<unsigned> latch{0} at the call site)
	inline bool first_this_run(std::atomic<unsigned> &latch) { return latch.exchange(run_id) != run_id; }
	//x^lx y^ly z^lz = sum_m sph2cart(l)[cart * n_spher(l) + m] R_lm, cartesians in WFN type order,
	//R_lm the real solid harmonics in ORCA/Gaussian order m = 0,+1,-1,+2,-2,...
	constexpr int n_cart(const int l) { return (l + 1) * (l + 2) / 2; };
	constexpr int n_spher(const int l) { return 2 * l + 1; };
	constexpr int first_type[11] = { 1, 2, 5, 11, 21, 36, 57, 85, 121, 166, 221 };
	const double* sph2cart(const int l);
	//ORCA stores g over sqrt(3) and h and up over sqrt((2l-1)!!) of the physical primitive norm; the tables follow ORCA, occ does not
	constexpr double sph2cart_norm2[11] = { 1, 1, 1, 1, 3, 945, 10395, 135135, 2027025, 34459425, 654729075 };
	constexpr double density_accuracy = 5.0e-5; // SQRT of the desired accuracy for density calculations
	//Electrons an RI (ab|P) triplet may move onto an atom before it is computed (see computeRho), sized for tsc within 1e-5
	//of max|f| and charges within 1e-4 e of the unscreened fit; -tune NOS_RI_SCREEN overrides it for one run.
	constexpr double ri_screen_threshold = 3e-7;
	//Overlap left between an orbital primitive pair and an aux shell, as erfc(sqrt(pq/(p+q)) R), below which computeRho
	//takes (ab|P) from the local expansion at the aux atom; -tune NOS_RI_FAR overrides it, 0 switches it off.
	constexpr double ri_far_erfc = 1e-14;
	constexpr int grid_max_no_flip = 50;
	int constexpr const_abs(int x)
	{
		return x < 0 ? -x : x;
	}

	template <typename K, typename V, std::size_t N>
	struct ConstexprMap {
		std::array<std::pair<K, V>, N> data;

		constexpr V at(const K &key) const {
			for (auto &[k, v] : data) {
				if (k == key) return v;
			}
			throw "Key not found";
		}
	};

	constexpr double constexpr_pow(double base, int exponent) {
		double result = 1.0;
		for (int i = 0; i < exponent; i++) result *= base;
		return result;
	}

	consteval double pow_const(const double &base, const int &exponent) {
		double result = 1.0;
		for (int i = 0; i < exponent; i++) result *= base;
		return result;
	}

	// Constants for later use
	constexpr double SQRT2 = sqrt(2.0);
	constexpr double SQRT3 = sqrt(3.0);
	constexpr double INV_SQRT2 = 1.0 / SQRT2;
	constexpr double INV_SQRT3 = 1.0 / SQRT3;
	constexpr int hardness = 3;
	constexpr double cutoff = 1.0e-20;
	constexpr double far_away = 10.0;
	constexpr double far_away_sq = far_away * far_away;
	constexpr double INV_PI = 1.0 / PI;
	constexpr double TWO_PI = 2 * PI;
	constexpr double FOUR_PI = 4 * PI;
	constexpr double INV_FOUR_PI = 1.0 / FOUR_PI;
	constexpr double EIGHT_PI = 8 * PI;
	constexpr double INV_EIGHT_PI = 1.0 / EIGHT_PI;
	constexpr double C0 = SQRT2 * FOUR_PI;
	constexpr double sqr_pi = sqrt(PI);
	constexpr double PI2 = pow_const(PI, 2);
	constexpr double PI3_2 = PI * sqr_pi;
	constexpr double PI3 = pow_const(PI, 3);
	constexpr double PI_180 = PI / 180.0;
	constexpr double INV_PI_180 = 1.0 / PI_180;
	constexpr double TG32 = 0.88622692545275801365;
	constexpr double ED_fact = 0.023934;
	constexpr int max_LT = 33;
	constexpr int MAG = 5810;
	//                       3,     5     7,    9,    11,   13,   15,   17
	//                      19,    21
	constexpr int lebedev_table[33] = { 6, 14, 26, 38, 50, 74, 86, 110,
									   146, 170, 194, 230, 266, 302, 350, 434,
									   590, 770, 974, 1202, 1454, 1730, 2030, 2354,
									   2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810 };
	constexpr long long int double_ft[25]{ 1, 1, 2, 3, 8, 15, 48, 105, 384, 945, 3840, 10395, 46080, 135135, 645120, 2027025, 10321920, 34459425, 185794560, 654729075, 3715891200, 13749310575, 81749606400, 316234143225, 1961990553600 };
	//Norm of the axial cartesian primitive x^l exp(-a r^2): (2a/PI)^(3/4) (4a)^(l/2) / sqrt((2l-1)!!), any l
	inline double axial_prim_norm(const int l, const double a)
	{
		return pow(pow(2, 3 + 4 * l) * pow(a, 2 * l + 3) / PI3 / pow(double_ft[std::max(2 * l - 1, 0)], 2), 0.25);
	}
	constexpr double alpha_coef = 0.1616204596739954813316614;
	constexpr double c_13 = 1.0 / 3.0;
	constexpr double c_43 = 4.0 / 3.0;
	constexpr double c_38 = 3.0 / 8.0;
	constexpr double c_16 = 1.0 / 6.0;
	constexpr double c_23 = 2.0 / 3.0;
	constexpr double c_53 = 5.0 / 3.0;
	constexpr double c_1_21 = 1.0 / 21.0;
	constexpr double c_1_35 = 1.0 / 35.0;
	constexpr double c_1_15 = 1.0 / 15.0;
	constexpr double c_3_40 = 3.0 / 40.0;
	constexpr double c_4_105 = 4.0 / 105.0;
	constexpr double c_1_105 = 1.0 / 105.0;
	constexpr double c_9_280 = 9.0 / 280.0;
	constexpr double c_m53 = -5.0 / 3.0;
	constexpr double keV_per_hartree = 0.027211386245988;
	constexpr double kcal_mol_per_hartree = 627.5094740631;
	constexpr double e_A2Debye = 0.2081943;                // in e*Angstrom
	constexpr double a0 = 0.529177210903E-10;              // in m
	constexpr std::complex<double> cnull = std::complex<double>(0.0, 0.0);

	constexpr double ctelf = 0.21940410461360746;

	constexpr int pow_2[] = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576 };

	constexpr size_t sod = sizeof(double);
	constexpr size_t soi = sizeof(int);
	constexpr size_t soli = sizeof(long int);

	constexpr double exp_approx(const double &x, int n = 25) {
		double result = 1.0;
		double term = 1.0;
		for (int i = 1; i <= n; i++) {
			term *= x / i;
			result += term;
		}
		return result;
	}


	//Conversion between ORCA and PySCF Orbital orderings
	// ORCA: 0,1,-1,2,-2,3,-3...
	// PySCF: -3,-2,-1,0,1,2,3...
	//MAP Layout: {l {m,new_index}}
	constexpr std::optional<std::size_t> orca_2_pySCF(unsigned int l, int m_idx) {
		switch (l) {
		case 0: { constexpr std::array<std::size_t, 1>  map = { 0 }; return map[m_idx]; }
			  //case 1: { constexpr std::array<std::size_t, 3>  map = { 2,0,1 }; return map[m_idx]; }  //For LibCint without -DPYPZPX="ON"
		case 1: { constexpr std::array<std::size_t, 3>  map = { 1,2,0 }; return map[m_idx]; }
		case 2: { constexpr std::array<std::size_t, 5>  map = { 2,3,1,4,0 }; return map[m_idx]; }
		case 3: { constexpr std::array<std::size_t, 7>  map = { 3,4,2,5,1,6,0 }; return map[m_idx]; }
		case 4: { constexpr std::array<std::size_t, 9>  map = { 4,5,3,6,2,7,1,8,0 }; return map[m_idx]; }
		case 5: { constexpr std::array<std::size_t, 11> map = { 5,6,4,7,3,8,2,9,1,10,0 }; return map[m_idx]; }
		case 6: { constexpr std::array<std::size_t, 13> map = { 6,7,5,8,4,9,3,10,2,11,1,12,0 }; return map[m_idx]; }
		case 7: { constexpr std::array<std::size_t, 15> map = { 7,8,6,9,5,10,4,11,3,12,2,13,1,14,0 }; return map[m_idx]; }
		case 8: { constexpr std::array<std::size_t, 17> map = { 8,9,7,10,6,11,5,12,4,13,3,14,2,15,1,16,0 }; return map[m_idx]; }
		case 9: { constexpr std::array<std::size_t, 19> map = { 9,10,8,11,7,12,6,13,5,14,4,15,3,16,2,17,1,18,0 }; return map[m_idx]; }
		case 10: { constexpr std::array<std::size_t, 21> map = { 10,11,9,12,8,13,7,14,6,15,5,16,4,17,3,18,2,19,1,20,0 }; return map[m_idx]; }
		default:
			return std::nullopt;
		}
	}

	constexpr unsigned int type_2_nbo(const int &type) {
		if (type == 1) return 1;
		//p types
		else if (type == 2) return 101;
		else if (type == 3) return 102;
		else if (type == 4) return 103;
		//d-types
		else if (type == 5) return 201;
		else if (type == 6) return 204;
		else if (type == 7) return 206;
		else if (type == 8) return 202;
		else if (type == 9) return 203;
		else if (type == 10) return 205;
		//f-types
		/*
		# 11 FXXX
		# 12 FYYY
		# 13 FZZZ
		# 14 FXXY
		# 15 FXXZ
		# 16 FYYZ
		# 17 FXYY
		# 18 FXZZ
		# 19 FYZZ
		# 20 FXYZ
		should become
		f (301-310) xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz
		*/
		else if (type == 11) return 301;
		else if (type == 12) return 307;
		else if (type == 13) return 310;
		else if (type == 14) return 302;
		else if (type == 15) return 303;
		else if (type == 16) return 308;
		else if (type == 17) return 304;
		else if (type == 18) return 306;
		else if (type == 19) return 309;
		else if (type == 20) return 305;
		//g-types
		/*
		# 21 GXXXX
		# 22 GYYYY
		# 23 GZZZZ
		# 24 GXXXY
		# 25 GXXXZ
		# 26 GXYYY
		# 27 GYYYZ
		# 28 GXZZZ
		# 29 GYZZZ
		# 30 GXXYY
		# 31 GXXZZ
		# 32 GYYZZ
		# 33 GXXYZ
		# 34 GXYYZ
		# 35 GXYZZ
		should become
		g (401-415) xxxx, xxxy, xxxz, xxyy, xxyz, xxzz, xyyy, xyyz, xyzz, xzzz, yyyy,  yyyz, yyzz, yzzz, zzzz
		*/
		else if (type == 21) return 401;
		else if (type == 22) return 411;
		else if (type == 23) return 415;
		else if (type == 24) return 402;
		else if (type == 25) return 403;
		else if (type == 26) return 407;
		else if (type == 27) return 412;
		else if (type == 28) return 414;
		else if (type == 29) return 413;
		else if (type == 30) return 404;
		else if (type == 31) return 406;
		else if (type == 32) return 410;
		else if (type == 33) return 405;
		else if (type == 34) return 408;
		else if (type == 35) return 409;
		else return 0; //error
	};

	//constexpr long long int ft[21]{ 1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200, 1307674368000, 20922789888000, 355687428096000, 6402373705728000, 121645100408832000, 2432902008176640000 };
	constexpr int MAX_FACTORIAL = 30;
	template <typename T>
	constexpr std::array<T, MAX_FACTORIAL> factorial_table() {
		std::array<T, MAX_FACTORIAL> f{};
		f[0] = 1;
		for (int i = 1; i < MAX_FACTORIAL; i++) f[i] = f[i - 1] * i;
		return f;
	}
	constexpr std::array<size_t, MAX_FACTORIAL> ft = factorial_table<size_t>();
	// ft[n] wraps silently from 21! on (> 2^64); use this one wherever n can exceed 20
	constexpr std::array<double, MAX_FACTORIAL> ftd = factorial_table<double>();

	constexpr double log_approx(const double &x, int n = 25) {
		if (x <= 0.0) return -1.0; // log is not defined for non-positive values
		double result = 0.0;
		double term = (x - 1) / (x + 1);
		double term_squared = term * term;
		double current_term = term;
		for (int i = 1; i <= n; i += 2) {
			result += current_term / i;
			current_term *= term_squared;
		}
		return 2 * result;
	}

	constexpr long long int ft_fun(const int &nr)
	{
		if (nr >= 0 && nr <= 20)
			return ft[nr];
		else if (nr < 0)
			return 0;
		else
			return ft_fun(nr - 1) * nr;
	}

	constexpr const double bohr2ang(const double &inp)
	{
		return inp * a0 * 1E10;
	}

	constexpr const double bohr2ang_p(const double &inp, const int &p)
	{
		if (p == 0)
			return 1.0;
		else if (p == 1)
			return bohr2ang(inp);
		else
			return bohr2ang_p(bohr2ang(inp), p - 1);
	}

	constexpr const double ang2bohr(const double &inp)
	{
		return inp / a0 * 1E-10;
	}
	constexpr const double ang2bohr_p(const double &inp, const int &p)
	{
		if (p == 0)
			return 1.0;
		else if (p == 1)
			return ang2bohr(inp);
		else
			return ang2bohr_p(ang2bohr(inp), p - 1);
	}

	constexpr const double cubic_ang2bohr(const double &inp)
	{
		return inp / ((a0 * a0 * a0) * 1E30);
	}

	constexpr const double cubic_bohr2ang(const double &inp)
	{
		return inp * (a0 * a0 * a0) * 1E30;
	}

	//------------------general functions for easy use of terminal input--------------------
	constexpr double bragg_angstrom[114]{
		0.00, // DUMMY LINE
		0.35, 0.35,
		1.45, 1.05, 0.85, 0.70, 0.65, 0.60, 0.50, 0.45,
		1.80, 1.50, 1.25, 1.10, 1.00, 1.00, 1.00, 1.00,
		2.20, 1.80, 1.60, 1.40, 1.35, 1.40, 1.40, 1.40, 1.35, 1.35, 1.35, 1.35, 1.30, 1.25, 1.15, 1.15, 1.15, 1.10,
		2.35, 2.00, 1.80, 1.55, 1.45, 1.45, 1.35, 1.30, 1.35, 1.40, 1.60, 1.55, 1.55, 1.45, 1.45, 1.40, 1.40, 1.40,
		2.60, 2.15, 1.95, 1.85, 1.85, 1.85, 1.85, 1.85, 1.85, 1.80, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.55, 1.45, 1.35, 1.30, 1.30, 1.35, 1.35, 1.35, 1.50, 1.90, 1.75, 1.60, 1.90, 1.50, 1.50,
		2.80, 2.35, 2.15, 2.05, 2.05, 2.05, 2.05, 2.05, 2.05, 2.00, 1.95, 1.95, 1.95, 1.95, 1.95, 1.95 };

	// Covalent Radii according to the CSD
	constexpr double covalent_radii[114]{
		0.0,
		0.23, 1.5,
		1.28, 0.96, 0.83, 0.68, 0.68, 0.68, 0.64, 1.5,
		1.66, 1.41, 1.21, 1.2, 1.05, 1.02, 0.99, 1.51,
		2.03, 1.76, 1.7, 1.6, 1.53, 1.39, 1.61, 1.52, 1.26, 1.24, 1.32, 1.22, 1.22, 1.17, 1.21, 1.22, 1.21, 1.5,
		2.2, 1.95, 1.9, 1.75, 1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.54, 1.42, 1.39, 1.39, 1.47, 1.4, 1.5,
		2.44, 2.15, 2.07, 2.04, 2.03, 2.01, 1.99, 1.98, 1.98, 1.96, 1.94, 1.92, 1.92, 1.89, 1.9, 1.87, 1.87, 1.75, 1.7, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32, 1.45, 1.46, 1.48, 1.4, 1.21, 1.5,
		2.6, 2.21, 2.15, 2.06, 2.00, 1.96, 1.9, 1.87, 1.8, 1.69, 1.54, 1.83, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5 };

	// van der Waals radii in Angstrom for d_norm: Bondi 1964 with the Mantina 2009 main-group additions, H = 1.09 as CrystalExplorer
	// ponytail: elements without a Bondi value (most of the d and f block) get 2.0, Alvarez 2013 has them at 2.0-2.3
	constexpr double vdW_radii[114]{
		0.0,
		1.09, 1.40,
		1.82, 1.53, 1.92, 1.70, 1.55, 1.52, 1.47, 1.54,
		2.27, 1.73, 1.84, 2.10, 1.80, 1.80, 1.75, 1.88,
		2.75, 2.31, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 1.63, 1.40, 1.39, 1.87, 2.11, 1.85, 1.90, 1.85, 2.02,
		3.03, 2.49, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 1.63, 1.72, 1.58, 1.93, 2.17, 2.06, 2.06, 1.98, 2.16,
		3.43, 2.68, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 1.75, 1.66, 1.55, 1.96, 2.02, 2.07, 1.97, 2.02, 2.20,
		3.48, 2.83, 2.00, 2.00, 2.00, 1.86, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00 };

	// Integer atom masses
	constexpr unsigned int integer_masses[]{
		1, 4,
		7, 9, 11, 12, 14, 16, 19, 20,
		23, 24, 27, 28, 31, 32, 35, 40,
		39, 40, 45, 48, 51, 52, 55, 56, 59, 58, 63, 64, 69, 74, 75, 80, 79, 84,
		85, 87, 88, 91, 92, 96, 98, 101, 103, 106, 108, 112, 115, 119, 122, 128, 127, 131,
		132, 137, 139, 140, 141, 144, 145, 150, 152, 157, 159, 163, 165, 167, 169, 173, 175, 178, 181, 184, 186, 190, 192, 195, 197, 201, 204, 207, 209, 209, 210, 222 };

	constexpr double real_masses[]{
		1.0079, 4.0026,
		6.941, 9.0122, 10.811, 12.011, 14.007, 15.999, 18.998, 20.18,
		22.99, 24.305, 26.986, 28.086, 30.974, 32.065, 35.453, 39.948,
		39.098, 40.078, 44.956, 47.867, 50.942, 51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.38, 69.723, 72.64, 74.922, 78.96, 79.904, 83.798,
		85.468, 87.62, 88.906, 91.224, 92.906, 95.96, 97.90, 101.07, 102.91, 106.42, 107.87, 112.41, 114.82, 118.71, 121.76, 127.6, 126.9, 131.29,
		132.91, 137.33, 139.91, 140.12, 140.91, 144.24, 144.9, 150.36, 151.96, 157.25, 158.93, 162.5, 164.93, 167.26, 168.93, 173.05, 174.97, 178.49, 180.95, 183.84, 186.21, 190.23, 192.22, 195.08, 196.97, 200.59, 204.38, 207.2, 208.98, 208.9, 209.9, 222.0 };

	// Avoid dynamic initialization from this widely included header.

	constexpr int ECP_electrons[] = { 0, 0, 0,
							 0, 0, 0, 0, 0, 0, 0, 0,
							 0, 0, 0, 0, 0, 0, 0, 0,
							 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28,
							 46, 46, 46, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60 };

	constexpr int MBIS_function[] = { 0, 1, 1,
							 2, 2, 2, 2, 2, 2, 2, 2,
							 3, 3, 3, 3, 3, 3, 3, 3,
							 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
							 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
							 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 };

	constexpr int ECP_electrons_xTB[] = { 0, 0, 0,
									 2, 2, 2, 2, 2, 2, 2, 2,
									 10, 10, 10, 10, 10, 10, 10, 10,
									 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 28, 28, 28, 28, 28, 28, 28,
									 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 46, 46, 46, 46, 46, 46, 46,
									 54, 54, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 68, 68, 68, 68, 68, 68, 68, 68, 78, 78, 78, 78, 78, 78, 78 };

	constexpr int ECP_electrons_pTB[] = { 0, 0, 0,
									 0, 0, 2, 2, 2, 2, 2, 2,
									 2, 2, 10, 10, 10, 10, 10, 10,
									 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 28, 28, 28, 28, 28, 28,
									 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 46, 46, 46, 46, 46, 46,
									 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 60, 60, 60, 60, 60, 60, 60, 60, 60, 78, 78, 78, 78, 78, 78 };

	//The tables stop at Rn (Z = 86), where the def2, xTB and pTB core definitions stop; direct indexing
	//reads past the end for heavier atoms. Beyond the table no core is defined: 0, the all-electron reading.
	template <std::size_t N>
	constexpr int ECP_core_electrons(const int(&table)[N], const int Z) noexcept
	{
		return (Z >= 0 && static_cast<std::size_t>(Z) < N) ? table[Z] : 0;
	}

	//Heaviest element the tables describe, for callers that warn instead of taking the zero
	constexpr int heaviest_ECP_element = static_cast<int>(sizeof(ECP_electrons) / sizeof(int)) - 1;
	static_assert(sizeof(ECP_electrons_xTB) == sizeof(ECP_electrons) &&
		sizeof(ECP_electrons_pTB) == sizeof(ECP_electrons),
		"the three ECP core tables are indexed by the same atomic number and reported by the same "
		"heaviest_ECP_element, so they have to end at the same element");

	constexpr const char *Labels[] = { "DM", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr" };
	const char *atnr2letter(const int &nr);

	int get_closest_num_angular(const int &n);

	int get_angular_order(const int &n);

	constexpr int get_Z_from_label(const char *tmp)
	{
		return
			(tmp[0] == 'H' && tmp[1] == '\0') ? 0 :
			(tmp[0] == 'D' && tmp[1] == '\0') ? 0 :
			(tmp[0] == 'T' && tmp[1] == '\0') ? 0 :
			(tmp[0] == 'H' && tmp[1] == 'e' && tmp[2] == '\0') ? 1 :
			(tmp[0] == 'L' && tmp[1] == 'i' && tmp[2] == '\0') ? 2 :
			(tmp[0] == 'B' && tmp[1] == 'e' && tmp[2] == '\0') ? 3 :
			(tmp[0] == 'B' && tmp[1] == '\0') ? 4 :
			(tmp[0] == 'C' && tmp[1] == '\0') ? 5 :
			(tmp[0] == 'N' && tmp[1] == '\0') ? 6 :
			(tmp[0] == 'O' && tmp[1] == '\0') ? 7 :
			(tmp[0] == 'F' && tmp[1] == '\0') ? 8 :
			(tmp[0] == 'N' && tmp[1] == 'e' && tmp[2] == '\0') ? 9 :
			(tmp[0] == 'N' && tmp[1] == 'a' && tmp[2] == '\0') ? 10 :
			(tmp[0] == 'M' && tmp[1] == 'g' && tmp[2] == '\0') ? 11 :
			(tmp[0] == 'A' && tmp[1] == 'l' && tmp[2] == '\0') ? 12 :
			(tmp[0] == 'S' && tmp[1] == 'i' && tmp[2] == '\0') ? 13 :
			(tmp[0] == 'P' && tmp[1] == '\0') ? 14 :
			(tmp[0] == 'S' && tmp[1] == '\0') ? 15 :
			(tmp[0] == 'C' && tmp[1] == 'l' && tmp[2] == '\0') ? 16 :
			(tmp[0] == 'A' && tmp[1] == 'r' && tmp[2] == '\0') ? 17 :
			(tmp[0] == 'K' && tmp[1] == '\0') ? 18 :
			(tmp[0] == 'C' && tmp[1] == 'a' && tmp[2] == '\0') ? 19 :
			(tmp[0] == 'S' && tmp[1] == 'c' && tmp[2] == '\0') ? 20 :
			(tmp[0] == 'T' && tmp[1] == 'i' && tmp[2] == '\0') ? 21 :
			(tmp[0] == 'V' && tmp[1] == '\0') ? 22 :
			(tmp[0] == 'C' && tmp[1] == 'r' && tmp[2] == '\0') ? 23 :
			(tmp[0] == 'M' && tmp[1] == 'n' && tmp[2] == '\0') ? 24 :
			(tmp[0] == 'F' && tmp[1] == 'e' && tmp[2] == '\0') ? 25 :
			(tmp[0] == 'C' && tmp[1] == 'o' && tmp[2] == '\0') ? 26 :
			(tmp[0] == 'N' && tmp[1] == 'i' && tmp[2] == '\0') ? 27 :
			(tmp[0] == 'C' && tmp[1] == 'u' && tmp[2] == '\0') ? 28 :
			(tmp[0] == 'Z' && tmp[1] == 'n' && tmp[2] == '\0') ? 29 :
			(tmp[0] == 'G' && tmp[1] == 'a' && tmp[2] == '\0') ? 30 :
			(tmp[0] == 'G' && tmp[1] == 'e' && tmp[2] == '\0') ? 31 :
			(tmp[0] == 'A' && tmp[1] == 's' && tmp[2] == '\0') ? 32 :
			(tmp[0] == 'S' && tmp[1] == 'e' && tmp[2] == '\0') ? 33 :
			(tmp[0] == 'B' && tmp[1] == 'r' && tmp[2] == '\0') ? 34 :
			(tmp[0] == 'K' && tmp[1] == 'r' && tmp[2] == '\0') ? 35 :
			(tmp[0] == 'R' && tmp[1] == 'b' && tmp[2] == '\0') ? 36 :
			(tmp[0] == 'S' && tmp[1] == 'r' && tmp[2] == '\0') ? 37 :
			(tmp[0] == 'Y' && tmp[1] == '\0') ? 38 :
			(tmp[0] == 'Z' && tmp[1] == 'r' && tmp[2] == '\0') ? 39 :
			(tmp[0] == 'N' && tmp[1] == 'b' && tmp[2] == '\0') ? 40 :
			(tmp[0] == 'M' && tmp[1] == 'o' && tmp[2] == '\0') ? 41 :
			(tmp[0] == 'T' && tmp[1] == 'c' && tmp[2] == '\0') ? 42 :
			(tmp[0] == 'R' && tmp[1] == 'u' && tmp[2] == '\0') ? 43 :
			(tmp[0] == 'R' && tmp[1] == 'h' && tmp[2] == '\0') ? 44 :
			(tmp[0] == 'P' && tmp[1] == 'd' && tmp[2] == '\0') ? 45 :
			(tmp[0] == 'A' && tmp[1] == 'g' && tmp[2] == '\0') ? 46 :
			(tmp[0] == 'C' && tmp[1] == 'd' && tmp[2] == '\0') ? 47 :
			(tmp[0] == 'I' && tmp[1] == 'n' && tmp[2] == '\0') ? 48 :
			(tmp[0] == 'S' && tmp[1] == 'n' && tmp[2] == '\0') ? 49 :
			(tmp[0] == 'S' && tmp[1] == 'b' && tmp[2] == '\0') ? 50 :
			(tmp[0] == 'T' && tmp[1] == 'e' && tmp[2] == '\0') ? 51 :
			(tmp[0] == 'I' && tmp[1] == '\0') ? 52 :
			(tmp[0] == 'X' && tmp[1] == 'e' && tmp[2] == '\0') ? 53 :
			(tmp[0] == 'C' && tmp[1] == 's' && tmp[2] == '\0') ? 54 :
			(tmp[0] == 'B' && tmp[1] == 'a' && tmp[2] == '\0') ? 55 :
			(tmp[0] == 'L' && tmp[1] == 'a' && tmp[2] == '\0') ? 56 :
			(tmp[0] == 'C' && tmp[1] == 'e' && tmp[2] == '\0') ? 57 :
			(tmp[0] == 'P' && tmp[1] == 'r' && tmp[2] == '\0') ? 58 :
			(tmp[0] == 'N' && tmp[1] == 'd' && tmp[2] == '\0') ? 59 :
			(tmp[0] == 'P' && tmp[1] == 'm' && tmp[2] == '\0') ? 60 :
			(tmp[0] == 'S' && tmp[1] == 'm' && tmp[2] == '\0') ? 61 :
			(tmp[0] == 'E' && tmp[1] == 'u' && tmp[2] == '\0') ? 62 :
			(tmp[0] == 'G' && tmp[1] == 'd' && tmp[2] == '\0') ? 63 :
			(tmp[0] == 'T' && tmp[1] == 'b' && tmp[2] == '\0') ? 64 :
			(tmp[0] == 'D' && tmp[1] == 'y' && tmp[2] == '\0') ? 65 :
			(tmp[0] == 'H' && tmp[1] == 'o' && tmp[2] == '\0') ? 66 :
			(tmp[0] == 'E' && tmp[1] == 'r' && tmp[2] == '\0') ? 67 :
			(tmp[0] == 'T' && tmp[1] == 'm' && tmp[2] == '\0') ? 68 :
			(tmp[0] == 'Y' && tmp[1] == 'b' && tmp[2] == '\0') ? 69 :
			(tmp[0] == 'L' && tmp[1] == 'u' && tmp[2] == '\0') ? 70 :
			(tmp[0] == 'H' && tmp[1] == 'f' && tmp[2] == '\0') ? 71 :
			(tmp[0] == 'T' && tmp[1] == 'a' && tmp[2] == '\0') ? 72 :
			(tmp[0] == 'W' && tmp[1] == '\0') ? 73 :
			(tmp[0] == 'R' && tmp[1] == 'e' && tmp[2] == '\0') ? 74 :
			(tmp[0] == 'O' && tmp[1] == 's' && tmp[2] == '\0') ? 75 :
			(tmp[0] == 'I' && tmp[1] == 'r' && tmp[2] == '\0') ? 76 :
			(tmp[0] == 'P' && tmp[1] == 't' && tmp[2] == '\0') ? 77 :
			(tmp[0] == 'A' && tmp[1] == 'u' && tmp[2] == '\0') ? 78 :
			(tmp[0] == 'H' && tmp[1] == 'g' && tmp[2] == '\0') ? 79 :
			(tmp[0] == 'T' && tmp[1] == 'l' && tmp[2] == '\0') ? 80 :
			(tmp[0] == 'P' && tmp[1] == 'b' && tmp[2] == '\0') ? 81 :
			(tmp[0] == 'B' && tmp[1] == 'i' && tmp[2] == '\0') ? 82 :
			(tmp[0] == 'P' && tmp[1] == 'o' && tmp[2] == '\0') ? 83 :
			(tmp[0] == 'A' && tmp[1] == 't' && tmp[2] == '\0') ? 84 :
			(tmp[0] == 'R' && tmp[1] == 'n' && tmp[2] == '\0') ? 85 :
			(tmp[0] == 'F' && tmp[1] == 'r' && tmp[2] == '\0') ? 86 :
			(tmp[0] == 'R' && tmp[1] == 'a' && tmp[2] == '\0') ? 87 :
			(tmp[0] == 'A' && tmp[1] == 'c' && tmp[2] == '\0') ? 88 :
			(tmp[0] == 'T' && tmp[1] == 'h' && tmp[2] == '\0') ? 89 :
			(tmp[0] == 'P' && tmp[1] == 'a' && tmp[2] == '\0') ? 90 :
			(tmp[0] == 'U' && tmp[1] == '\0') ? 91 :
			(tmp[0] == 'N' && tmp[1] == 'p' && tmp[2] == '\0') ? 92 :
			(tmp[0] == 'P' && tmp[1] == 'u' && tmp[2] == '\0') ? 93 :
			(tmp[0] == 'A' && tmp[1] == 'm' && tmp[2] == '\0') ? 94 :
			(tmp[0] == 'C' && tmp[1] == 'm' && tmp[2] == '\0') ? 95 :
			(tmp[0] == 'B' && tmp[1] == 'k' && tmp[2] == '\0') ? 96 :
			(tmp[0] == 'C' && tmp[1] == 'f' && tmp[2] == '\0') ? 97 :
			(tmp[0] == 'E' && tmp[1] == 's' && tmp[2] == '\0') ? 98 :
			(tmp[0] == 'F' && tmp[1] == 'm' && tmp[2] == '\0') ? 99 :
			(tmp[0] == 'M' && tmp[1] == 'd' && tmp[2] == '\0') ? 100 :
			(tmp[0] == 'N' && tmp[1] == 'o' && tmp[2] == '\0') ? 101 :
			(tmp[0] == 'L' && tmp[1] == 'r' && tmp[2] == '\0') ? 102 :
			(tmp[0] == 'R' && tmp[1] == 'f' && tmp[2] == '\0') ? 103 :
			(tmp[0] == 'D' && tmp[1] == 'b' && tmp[2] == '\0') ? 104 :
			(tmp[0] == 'S' && tmp[1] == 'g' && tmp[2] == '\0') ? 105 :
			(tmp[0] == 'B' && tmp[1] == 'h' && tmp[2] == '\0') ? 106 :
			(tmp[0] == 'H' && tmp[1] == 's' && tmp[2] == '\0') ? 107 :
			(tmp[0] == 'M' && tmp[1] == 't' && tmp[2] == '\0') ? 108 :
			(tmp[0] == 'D' && tmp[1] == 's' && tmp[2] == '\0') ? 109 :
			(tmp[0] == 'R' && tmp[1] == 'g' && tmp[2] == '\0') ? 110 :
			-1;
	};

	constexpr int type_vector[858]{
	0, 0, 0,
	1, 0, 0,
	0, 1, 0,
	0, 0, 1,
	2, 0, 0,
	0, 2, 0,
	0, 0, 2,
	1, 1, 0,
	1, 0, 1,
	0, 1, 1,
	3, 0, 0,
	0, 3, 0,
	0, 0, 3,
	2, 1, 0,
	2, 0, 1,
	0, 2, 1,
	1, 2, 0,
	1, 0, 2,
	0, 1, 2,
	1, 1, 1,
	0, 0, 4,
	0, 1, 3,
	0, 2, 2,
	0, 3, 1,
	0, 4, 0,
	1, 0, 3,
	1, 1, 2,
	1, 2, 1,
	1, 3, 0,
	2, 0, 2,
	2, 1, 1,
	2, 2, 0,
	3, 0, 1,
	3, 1, 0,
	4, 0, 0,
	0, 0, 5,
	0, 1, 4,
	0, 2, 3,
	0, 3, 2,
	0, 4, 1,
	0, 5, 0,
	1, 0, 4,
	1, 1, 3,
	1, 2, 2,
	1, 3, 1,
	1, 4, 0,
	2, 0, 3,
	2, 1, 2,
	2, 2, 1,
	2, 3, 0,
	3, 0, 2,
	3, 1, 1,
	3, 2, 0,
	4, 0, 1,
	4, 1, 0,
	5, 0, 0,
	0, 0, 6,
	0, 1, 5,
	0, 2, 4,
	0, 3, 3,
	0, 4, 2,
	0, 5, 1,
	0, 6, 0,
	1, 0, 5,
	1, 1, 4,
	1, 2, 3,
	1, 3, 2,
	1, 4, 1,
	1, 5, 0,
	2, 0, 4,
	2, 1, 3,
	2, 2, 2,
	2, 3, 1,
	2, 4, 0,
	3, 0, 3,
	3, 1, 2,
	3, 2, 1,
	3, 3, 0,
	4, 0, 2,
	4, 1, 1,
	4, 2, 0,
	5, 0, 1,
	5, 1, 0,
	6, 0, 0,
	0, 0, 7,
	0, 1, 6,
	0, 2, 5,
	0, 3, 4,
	0, 4, 3,
	0, 5, 2,
	0, 6, 1,
	0, 7, 0,
	1, 0, 6,
	1, 1, 5,
	1, 2, 4,
	1, 3, 3,
	1, 4, 2,
	1, 5, 1,
	1, 6, 0,
	2, 0, 5,
	2, 1, 4,
	2, 2, 3,
	2, 3, 2,
	2, 4, 1,
	2, 5, 0,
	3, 0, 4,
	3, 1, 3,
	3, 2, 2,
	3, 3, 1,
	3, 4, 0,
	4, 0, 3,
	4, 1, 2,
	4, 2, 1,
	4, 3, 0,
	5, 0, 2,
	5, 1, 1,
	5, 2, 0,
	6, 0, 1,
	6, 1, 0,
	7, 0, 0,
	0, 0, 8,
	0, 1, 7,
	0, 2, 6,
	0, 3, 5,
	0, 4, 4,
	0, 5, 3,
	0, 6, 2,
	0, 7, 1,
	0, 8, 0,
	1, 0, 7,
	1, 1, 6,
	1, 2, 5,
	1, 3, 4,
	1, 4, 3,
	1, 5, 2,
	1, 6, 1,
	1, 7, 0,
	2, 0, 6,
	2, 1, 5,
	2, 2, 4,
	2, 3, 3,
	2, 4, 2,
	2, 5, 1,
	2, 6, 0,
	3, 0, 5,
	3, 1, 4,
	3, 2, 3,
	3, 3, 2,
	3, 4, 1,
	3, 5, 0,
	4, 0, 4,
	4, 1, 3,
	4, 2, 2,
	4, 3, 1,
	4, 4, 0,
	5, 0, 3,
	5, 1, 2,
	5, 2, 1,
	5, 3, 0,
	6, 0, 2,
	6, 1, 1,
	6, 2, 0,
	7, 0, 1,
	7, 1, 0,
	8, 0, 0,
	0, 0, 9,
	0, 1, 8,
	0, 2, 7,
	0, 3, 6,
	0, 4, 5,
	0, 5, 4,
	0, 6, 3,
	0, 7, 2,
	0, 8, 1,
	0, 9, 0,
	1, 0, 8,
	1, 1, 7,
	1, 2, 6,
	1, 3, 5,
	1, 4, 4,
	1, 5, 3,
	1, 6, 2,
	1, 7, 1,
	1, 8, 0,
	2, 0, 7,
	2, 1, 6,
	2, 2, 5,
	2, 3, 4,
	2, 4, 3,
	2, 5, 2,
	2, 6, 1,
	2, 7, 0,
	3, 0, 6,
	3, 1, 5,
	3, 2, 4,
	3, 3, 3,
	3, 4, 2,
	3, 5, 1,
	3, 6, 0,
	4, 0, 5,
	4, 1, 4,
	4, 2, 3,
	4, 3, 2,
	4, 4, 1,
	4, 5, 0,
	5, 0, 4,
	5, 1, 3,
	5, 2, 2,
	5, 3, 1,
	5, 4, 0,
	6, 0, 3,
	6, 1, 2,
	6, 2, 1,
	6, 3, 0,
	7, 0, 2,
	7, 1, 1,
	7, 2, 0,
	8, 0, 1,
	8, 1, 0,
	9, 0, 0,
	0, 0, 10,
	0, 1, 9,
	0, 2, 8,
	0, 3, 7,
	0, 4, 6,
	0, 5, 5,
	0, 6, 4,
	0, 7, 3,
	0, 8, 2,
	0, 9, 1,
	0, 10, 0,
	1, 0, 9,
	1, 1, 8,
	1, 2, 7,
	1, 3, 6,
	1, 4, 5,
	1, 5, 4,
	1, 6, 3,
	1, 7, 2,
	1, 8, 1,
	1, 9, 0,
	2, 0, 8,
	2, 1, 7,
	2, 2, 6,
	2, 3, 5,
	2, 4, 4,
	2, 5, 3,
	2, 6, 2,
	2, 7, 1,
	2, 8, 0,
	3, 0, 7,
	3, 1, 6,
	3, 2, 5,
	3, 3, 4,
	3, 4, 3,
	3, 5, 2,
	3, 6, 1,
	3, 7, 0,
	4, 0, 6,
	4, 1, 5,
	4, 2, 4,
	4, 3, 3,
	4, 4, 2,
	4, 5, 1,
	4, 6, 0,
	5, 0, 5,
	5, 1, 4,
	5, 2, 3,
	5, 3, 2,
	5, 4, 1,
	5, 5, 0,
	6, 0, 4,
	6, 1, 3,
	6, 2, 2,
	6, 3, 1,
	6, 4, 0,
	7, 0, 3,
	7, 1, 2,
	7, 2, 1,
	7, 3, 0,
	8, 0, 2,
	8, 1, 1,
	8, 2, 0,
	9, 0, 1,
	9, 1, 0,
	10, 0, 0 };

	// https://en.wikipedia.org/wiki/Electronegativities_of_the_elements_(data_page)
	constexpr double allen_electronegativities[] =
	{ 2.300, 4.160,
		0.912, 1.576, 2.051, 2.544, 3.066, 3.610, 4.193, 4.787,
		0.869, 1.293, 1.613, 1.916, 2.253, 2.589, 2.869, 3.242,
		0.734, 1.034, 1.19, 1.38, 1.53, 1.65, 1.75, 1.80, 1.84, 1.88, 1.85, 1.588, 1.756, 1.994, 2.211, 2.424, 2.685, 2.966,
		0.706, 0.963, 1.12, 1.32, 1.41, 1.47, 1.51, 1.54, 1.56, 1.58, 1.87, 1.521, 1.656, 1.824, 1.984, 2.158, 2.359, 2.582,
		0.659, 0.881, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.09, 1.16, 1.34, 1.47, 1.60, 1.65, 1.68, 1.72, 1.92, 1.765, 1.789, 1.854, 2.01, 2.19, 2.39, 2.60,
		0.67, 0.89
	};

	constexpr void type2vector(int index, int *vector)
	{
		if (index < 1 || index > 286)
		{
			vector[0] = -1;
			vector[1] = -1;
			vector[2] = -1;
			return;
		}
		const int temp = index - 1;
		vector[0] = constants::type_vector[temp * 3];
		vector[1] = constants::type_vector[temp * 3 + 1];
		vector[2] = constants::type_vector[temp * 3 + 2];
	};

	const double normgauss(const int &type, const double &exp);

	const double spherical_harmonic(const int &l, const int &m, const double *d);
	const double spherical_harmonic_max(const int& l, const int& m);
	const double spherical_harmonic(const int &l, const double *d, const double *coefs);

	constexpr int ASSOCIATED_LEGENDRE_MAX_L = 12;
	using SphericalNormsArray = std::array<std::array<double, 2 * ASSOCIATED_LEGENDRE_MAX_L + 1>, ASSOCIATED_LEGENDRE_MAX_L + 1>;
	// Define constexpr storage generator
	constexpr SphericalNormsArray generate_spherical_norms() {
		SphericalNormsArray norms = {};

		// l = 0 case
		norms[0][0] = constants::c_1_4p;
		for (int l = 1; l <= ASSOCIATED_LEGENDRE_MAX_L; l++) {
			for (int m = -l; m <= l; m++) {
				norms[l][m + l] = sqrt(((2 * l + 1) * constants::ftd[l - m]) / (constants::TWO_PI * constants::ftd[l + m]));
			}
			norms[l][l] = sqrt(((2 * l + 1) / constants::FOUR_PI));
		}
		return norms;
	}
	constexpr SphericalNormsArray spherical_norms = generate_spherical_norms();



	double associated_legendre_polynomial(const int &l, const int &m, const double &x);
	std::vector<double> cartesian_to_spherical(const double &x, const double &y, const double &z);
	std::pair<double, double> norm_cartesian_to_spherical(const double &x, const double &y, const double &z);
	//Original implementation after P. Coppens DOI: 10.1107/97809553602060000759 Eq. 1.2.7.2b
	//I omitted the abs(m) in the factorial as most other sources do not include it
	// N_l|m| P_l^|m|(cos theta) {cos, sin}(|m| phi); P through the tables and their recurrence
	// on every platform (libc++ has no std::assoc_legendre)
	inline double real_spherical(const int &l, const int &m, const double &theta, const double &phi) {
		const int am = std::abs(m);
		return constants::spherical_norms[l][l + am] * associated_legendre_polynomial(l, am, cos(theta)) * ((m >= 0) ? cos(am * phi) : sin(am * phi));
	}

}
