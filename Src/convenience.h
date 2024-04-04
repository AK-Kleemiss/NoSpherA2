#pragma once
#define WIN32_LEAN_AND_MEAN
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#ifdef __cplusplus__
#include <cstdlib>
#else
#include <stdlib.h>
#endif
#include <functional>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <regex>
#include <set>
#include <string>
#include <stdexcept>
#include <sstream>
#include <typeinfo>
#include <vector>

// Here are the system specific libaries
#ifdef _WIN32
#include <direct.h>
#define GetCurrentDir _getcwd(NULL, 0)
#include <io.h>
#else
#define GetCurrentDir getcwd
#include <unistd.h>
#include <cfloat>
#include <sys/wait.h>
#include <termios.h>
#include <cstring>
#endif

// Pre-definition of classes included later
class WFN;
class cell;
class atom;

typedef std::complex<double> cdouble;
typedef std::vector<double> vec;
typedef std::vector<int> ivec;
typedef std::vector<cdouble> cvec;
typedef std::chrono::high_resolution_clock::time_point time_point;

inline double vec_sum(vec &in)
{
	double res = 0.0;
	for (int i = 0; i < in.size(); i++)
		res += in[i];
	return res;
}

inline cdouble vec_sum(cvec &in)
{
	cdouble res = 0.0;
	for (int i = 0; i < in.size(); i++)
		res += in[i];
	return res;
}

inline const std::complex<double> c_one(0, 1.0);

std::string help_message();
std::string NoSpherA2_message();
std::string build_date();

namespace constants
{
#include <limits>

	namespace Detail
	{
		double constexpr sqrtNewtonRaphson(double x, double curr, double prev)
		{
			return curr == prev
					   ? curr
					   : sqrtNewtonRaphson(x, 0.5 * (curr + x / curr), curr);
		}
	}

	/*
	 * Constexpr version of the square root
	 * Return value:
	 *   - For a finite and non-negative value of "x", returns an approximation for the square root of "x"
	 *   - Otherwise, returns NaN
	 * Taken from https://stackoverflow.com/questions/8622256/in-c11-is-sqrt-defined-as-constexpr
	 */
	double constexpr sqrt(double x)
	{
		return x >= 0 && x < std::numeric_limits<double>::infinity()
				   ? Detail::sqrtNewtonRaphson(x, x, 0)
				   : std::numeric_limits<double>::quiet_NaN();
	}
	// Constants for later use
	constexpr double SQRT2 = sqrt(2.0);
	constexpr int hardness = 3;
	constexpr double cutoff = 1.0e-20;
	constexpr double PI = 3.14159265358979323846;
	constexpr double TWO_PI = 2 * PI;
	constexpr double FOUR_PI = 4 * PI;
	constexpr double C0 = SQRT2 * FOUR_PI;
	const double sqr_pi = sqrt(PI);
	constexpr double PI2 = PI * PI;
	constexpr double PI3 = PI * PI * PI;
	constexpr double PI_180 = PI / 180.0;
	const double TG32 = tgamma(3.0 / 2.0);
	constexpr double ED_fact = 0.023934;
	constexpr int max_LT = 33;
	constexpr int MAG = 5810;
	//                       3,     5     7,    9,    11,   13,   15,   17
	//                      19,    21
	constexpr int lebedev_table[33] = {6, 14, 26, 38, 50, 74, 86, 110,
									   146, 170, 194, 230, 266, 302, 350, 434,
									   590, 770, 974, 1202, 1454, 1730, 2030, 2354,
									   2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810};
	constexpr long long int ft[21]{1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200, 1307674368000, 20922789888000, 355687428096000, 6402373705728000, 121645100408832000, 2432902008176640000};
	constexpr double alpha_coef = 0.1616204596739954813316614;
	constexpr double c_43 = 4.0 / 3.0;
	constexpr double c_38 = 3.0 / 8.0;
	constexpr double c_m53 = -5.0 / 3.0;
	constexpr double barnsbohr = 2.80028520539078E+7;
	constexpr double fine_struct = 7.2973525693E-3;
	constexpr double inv_fine_struct = 1 / fine_struct;
	constexpr double fine_pi = inv_fine_struct / TWO_PI / PI;
	constexpr double inv_fine_mod = inv_fine_struct / FOUR_PI;
	constexpr double keV_per_hartree = 0.027211386245988;
	constexpr double angstrom2eV = 1.23984193 * 10000;
	constexpr double angstrom2keV = 12.3984193;
	constexpr double f_to_mu = 4208.031548;
	constexpr double barns_to_electrons = 1.43110541E-8;
	constexpr double a0 = 0.529177210903E-10;			   // in m
	constexpr double h = 6.62607015E-34 / 1.602176634E-19; // in eV*s
	constexpr double Ryd_ener = 13.6056923;				   // in eV
	constexpr double alpha = 0.0072973525693;			   // Sommerfeld fine structure constant
	constexpr double el_mass = 9.1093837015E-31;		   // in kg
	constexpr double el_charge = 1.602176634E-19;		   // in C
	constexpr double speed_of_light = 2.99792458E8;		   // m/s

	const double ctelf = 10 * pow(2, -2.0 / 3.0) * pow(3, c_m53) * pow(PI, -c_43);
	constexpr double c_1_4p = sqrt(1.0 / (FOUR_PI));
	constexpr double c_3_4p = sqrt(3.0 / (FOUR_PI));
	constexpr double c_5_16p = sqrt(5.0 / (16.0 * PI));
	constexpr double c_7_16p = sqrt(7.0 / (16.0 * PI));
	constexpr double c_9_256p = sqrt(9.0 / (256.0 * PI));
	constexpr double c_11_256p = sqrt(11.0 / (256.0 * PI));
	constexpr double c_13_1024p = sqrt(13.0 / (1024.0 * PI));
	constexpr double c_15_4p = sqrt(15.0 / (FOUR_PI));
	constexpr double c_15_16p = sqrt(15.0 / (16.0 * PI));
	constexpr double c_21_32p = sqrt(21.0 / (32.0 * PI));
	constexpr double c_35_32p = sqrt(35.0 / (32.0 * PI));
	constexpr double c_45_16p = sqrt(45.0 / (16.0 * PI));
	constexpr double c_45_32p = sqrt(45.0 / (32.0 * PI));
	constexpr double c_45_64p = sqrt(45.0 / (64.0 * PI));
	constexpr double c_105_4p = sqrt(105.0 / (FOUR_PI));
	constexpr double c_105_16p = sqrt(105.0 / (16.0 * PI));
	constexpr double c_165_256p = sqrt(165.0 / (256.0 * PI));
	constexpr double c_273_256p = sqrt(273.0 / (256.0 * PI));
	constexpr double c_315_16p = sqrt(315.0 / (16.0 * PI));
	constexpr double c_315_32p = sqrt(315.0 / (32.0 * PI));
	constexpr double c_315_256p = sqrt(315.0 / (256.0 * PI));
	constexpr double c_385_512p = sqrt(385.0 / (512.0 * PI));
	constexpr double c_693_2048p = sqrt(693.0 / (2048.0 * PI));
	constexpr double c_1155_64p = sqrt(1155.0 / (64.0 * PI));
	constexpr double c_3465_256p = sqrt(3465.0 / (256.0 * PI));

	inline const long long int ft_fun(const int &nr)
	{
		if (nr >= 0 && nr <= 20)
			return ft[nr];
		else if (nr < 0)
			return 0;
		else
			return ft_fun(nr - 1) * nr;
	}

	constexpr double bohr2ang(const double &inp)
	{
		return inp * 0.529177249;
	}

	constexpr double bohr2ang_p(const double &inp, const int &p)
	{
		if (p == 0)
			return 1.0;
		else if (p == 1)
			return bohr2ang(inp);
		else
			return bohr2ang_p(bohr2ang(inp), p - 1);
	}

	constexpr double ang2bohr(const double &inp)
	{
		return inp / 0.529177249;
	}
	inline const double ang2bohr_p(const double &inp, const int &p)
	{
		if (p == 0)
			return 1.0;
		else if (p == 1)
			return ang2bohr(inp);
		else
			return ang2bohr_p(ang2bohr(inp), p - 1);
	}

	constexpr double cubic_ang2bohr(const double &inp)
	{
		return inp / (0.529177249 * 0.529177249 * 0.529177249);
	}

	constexpr double cubic_bohr2ang(const double &inp)
	{
		return inp * (0.529177249 * 0.529177249 * 0.529177249);
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
		2.80, 2.35, 2.15, 2.05, 2.05, 2.05, 2.05, 2.05, 2.05, 2.00, 1.95, 1.95, 1.95, 1.95, 1.95, 1.95};

	// Covalent Radii according to the CSD
	constexpr double covalent_radii[114]{
		0.0,
		0.23, 1.5,
		1.28, 0.96, 0.83, 0.68, 0.68, 0.68, 0.64, 1.5,
		1.66, 1.41, 1.21, 1.2, 1.05, 1.02, 0.99, 1.51,
		2.03, 1.76, 1.7, 1.6, 1.53, 1.39, 1.61, 1.52, 1.26, 1.24, 1.32, 1.22, 1.22, 1.17, 1.21, 1.22, 1.21, 1.5,
		2.2, 1.95, 1.9, 1.75, 1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.54, 1.42, 1.39, 1.39, 1.47, 1.4, 1.5,
		2.44, 2.15, 2.07, 2.04, 2.03, 2.01, 1.99, 1.98, 1.98, 1.96, 1.94, 1.92, 1.92, 1.89, 1.9, 1.87, 1.87, 1.75, 1.7, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32, 1.45, 1.46, 1.48, 1.4, 1.21, 1.5,
		2.6, 2.21, 2.15, 2.06, 2.00, 1.96, 1.9, 1.87, 1.8, 1.69, 1.54, 1.83, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5};

	// Integer atom masses
	constexpr unsigned int integer_masses[]{
		1, 4,
		7, 9, 11, 12, 14, 16, 19, 20,
		23, 24, 27, 28, 31, 32, 35, 40,
		39, 40, 45, 48, 51, 52, 55, 56, 59, 58, 63, 64, 69, 74, 75, 80, 79, 84,
		85, 87, 88, 91, 92, 96, 98, 101, 103, 106, 108, 112, 115, 119, 122, 128, 127, 131,
		132, 137, 139, 140, 141, 144, 145, 150, 152, 157, 159, 163, 165, 167, 169, 173, 175, 178, 181, 184, 186, 190, 192, 195, 197, 201, 204, 207, 209, 209, 210, 222};

	constexpr double real_masses[]{
		1.0079, 4.0026,
		6.941, 9.0122, 10.811, 12.011, 14.007, 15.999, 18.998, 20.18,
		22.99, 24.305, 26.986, 28.086, 30.974, 32.065, 35.453, 39.948,
		39.098, 40.078, 44.956, 47.867, 50.942, 51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.38, 69.723, 72.64, 74.922, 78.96, 79.904, 83.798,
		85.468, 87.62, 88.906, 91.224, 92.906, 95.96, 97.90, 101.07, 102.91, 106.42, 107.87, 112.41, 114.82, 118.71, 121.76, 127.6, 126.9, 131.29,
		132.91, 137.33, 139.91, 140.12, 140.91, 144.24, 144.9, 150.36, 151.96, 157.25, 158.93, 162.5, 164.93, 167.26, 168.93, 173.05, 174.97, 178.49, 180.95, 183.84, 186.21, 190.23, 192.22, 195.08, 196.97, 200.59, 204.38, 207.2, 208.98, 208.9, 209.9, 222.0};
}
// bool yesno();
bool is_similar_rel(const double &first, const double &second, const double &tolerance);
bool is_similar(const double &first, const double &second, const double &tolerance);
bool is_similar_abs(const double &first, const double &second, const double &tolerance);
void cls();
std::string get_home_path(void);
void join_path(std::string &s1, std::string &s2);
inline char asciitolower(char in)
{
	if (in <= 'Z' && in >= 'A')
		return in - ('Z' - 'z');
	return in;
}
inline void error_check(const bool condition, const std::string &file, const int &line, const std::string &function, const std::string &error_mesasge, std::ostream &log_file = std::cout)
{
	if (!condition)
	{
		log_file << "Error in " << function << " at: " << file << " : " << line << " " << error_mesasge << std::endl;
		log_file.flush();
		exit(-1);
	}
};
inline void not_implemented(const std::string &file, const int &line, const std::string &function, const std::string &error_mesasge, std::ostream &log_file)
{
	log_file << function << " at: " << file << ":" << line << " " << error_mesasge << " not yet implemented!" << std::endl;
	log_file.flush();
	exit(-1);
};
#define err_checkf(condition, error_message, file) error_check(condition, __FILE__, __LINE__, __func__, error_message, file)
#define err_chkf(condition, error_message, file) error_check(condition, __FILE__, __LINE__, __func__, error_message, file)
#define err_chekf(condition, error_message, file) error_check(condition, __FILE__, __LINE__, __func__, error_message, file)
#define err_not_impl_f(error_message, file) not_implemented(__FILE__, __LINE__, __func__, error_message, file)

bool generate_sph2cart_mat(std::vector<vec> &p, std::vector<vec> &d, std::vector<vec> &f, std::vector<vec> &g);
bool generate_cart2sph_mat(std::vector<vec> &d, std::vector<vec> &f, std::vector<vec> &g, std::vector<vec> &h);
std::string go_get_string(std::ifstream &file, std::string search, bool rewind = true);

inline const int sht2nbas(const int &type)
{
	const int st2bas[6]{1, 3, 6, 10, 15, 21};
	const int nst2bas[6]{11, 9, 7, 5, 4, 1};
	if (type >= 0)
		return st2bas[type];
	else
		return nst2bas[5 + type];
};

inline const int shell2function(const int &type, const int &prim)
{
	switch (type)
	{
	case (-5):
		return -32 + prim;
	case (-4):
		return -21 + prim;
	case (-3):
		return -12 + prim;
	case (-2):
		return -5 + prim;
	case (-1):
		return 1 + prim;
	case (0):
		return 1;
	case (1):
		return 2 + prim;
	case (2):
		return 5 + prim;
	case (3):
		if (prim == 0)
			return 11;
		if (prim == 1)
			return 12;
		if (prim == 2)
			return 13;
		if (prim == 3)
			return 17;
		if (prim == 4)
			return 14;
		if (prim == 5)
			return 15;
		if (prim == 6)
			return 18;
		if (prim == 7)
			return 19;
		if (prim == 8)
			return 16;
		if (prim == 9)
			return 20;
		break;
	case (4):
		return 21 + prim;
	case (5):
		return 36 + prim;
	default:
		return 0;
	}
	return 0;
}

const double normgauss(const int &type, const double &exp);
const double spherical_harmonic(const int &l, const int &m, const double *d);

template <class T>
std::string toString(const T &t)
{
	std::ostringstream stream;
	stream << t;
	return stream.str();
}

template <class T>
T fromString(const std::string &s)
{
	std::istringstream stream(s);
	T t;
	stream >> t;
	return t;
}

template <typename T>
void shrink_vector(std::vector<T> &g)
{
	g.clear();
	std::vector<T>(g).swap(g);
}

template <class T>
std::vector<T> split_string(const std::string &input, const std::string delimiter)
{
	std::string input_copy = input + delimiter; // Need to add one delimiter in the end to return all elements
	std::vector<T> result;
	size_t pos = 0;
	while ((pos = input_copy.find(delimiter)) != std::string::npos)
	{
		result.push_back(fromString<T>(input_copy.substr(0, pos)));
		input_copy.erase(0, pos + delimiter.length());
	}
	return result;
};

inline void remove_empty_elements(std::vector<std::string> &input, const std::string &empty = " ")
{
	for (int i = (int)input.size() - 1; i >= 0; i--)
		if (input[i] == empty || input[i] == "")
			input.erase(input.begin() + i);
}

inline std::chrono::high_resolution_clock::time_point get_time()
{
	// gets the current time using std chrono library
	std::chrono::high_resolution_clock::time_point time = std::chrono::high_resolution_clock::now();
	return time;
}

inline int get_µsec(std::chrono::high_resolution_clock::time_point start, std::chrono::high_resolution_clock::time_point end)
{
  // gets the time difference in microseconds
  std::chrono::microseconds µsec = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  return µsec.count();
}

inline int get_msec(std::chrono::high_resolution_clock::time_point start, std::chrono::high_resolution_clock::time_point end)
{
	// gets the time difference in milliseconds
	std::chrono::milliseconds msec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	return msec.count();
}

inline int get_sec(std::chrono::high_resolution_clock::time_point start, std::chrono::high_resolution_clock::time_point end)
{
	// gets the time difference in seconds
	std::chrono::seconds sec = std::chrono::duration_cast<std::chrono::seconds>(end - start);
	return sec.count();
}

inline void write_timing_to_file(std::ostream &file,
						  time_point start,
						  time_point end,
						  time_point end_prototypes,
						  time_point end_becke,
						  time_point end_spherical,
						  time_point end_prune,
						  time_point end_aspherical,
						  time_point before_kpts,
						  time_point after_kpts,
						  time_point end1)
{
	// writes the timing of different things to a file
	using namespace std;
	int dur = get_sec(start, end);

	if (dur < 1)
		file << "Total Time: " << fixed << setprecision(0) << get_msec(start, end) << " ms\n";
	else if (dur < 60)
		file << "Total Time: " << fixed << setprecision(0) << dur << " s\n";
	else if (dur < 3600)
		file << "Total Time: " << fixed << setprecision(0) << floor(dur / 60) << " m " << dur % 60 << " s\n";
	else
		file << "Total Time: " << fixed << setprecision(0) << floor(dur / 3600) << " h " << (dur % 3600) / 60 << " m\n";
	file << endl;
	file << "Time Breakdown:" << endl;
	if (get_sec(start, end_prototypes) > 1)
		if (get_sec(start, end_prototypes) < 100)
			file << " ... for Prototype Grid setup:" << setw(6) << get_sec(start, end_prototypes) << " s " << get_msec(start, end_prototypes) % 1000 << " ms" << endl;
		else
			file << " ... for Prototype Grid setup:" << setw(6) << get_sec(start, end_prototypes) << " s" << endl;
	else
		file << " ... for Prototype Grid setup:" << setw(6) << get_msec(start, end_prototypes) << " ms" << endl;
	if (get_sec(end_prototypes, end_becke) > 1)
		if (get_sec(end_prototypes, end_becke) < 100)
      file << " ... for Becke Grid setup:    " << setw(6) << get_sec(end_prototypes, end_becke) << " s " << get_msec(end_prototypes, end_becke) % 1000 << " ms" << endl;
    else
      file << " ... for Becke Grid setup:    " << setw(6) << get_sec(end_prototypes, end_becke) << " s" << endl;
	else
		file << " ... for Becke Grid setup:    " << setw(6) << get_msec(end_prototypes, end_becke) << " ms" << endl;
	if (get_sec(end_becke, end_spherical) > 1)
		if (get_sec(end_becke, end_spherical) < 100)
      file << " ... for spherical density:   " << setw(6) << get_sec(end_becke, end_spherical) << " s " << get_msec(end_becke, end_spherical) % 1000 << " ms" << endl;
    else
      file << " ... for spherical density:   " << setw(6) << get_sec(end_becke, end_spherical) << " s" << endl;
	else
		file << " ... for spherical density:   " << setw(6) << get_msec(end_becke, end_spherical) << " ms" << endl;
	if (get_sec(end_spherical, end_prune) > 1)
		if (get_sec(end_spherical, end_prune) < 100)
      file << " ... for Grid Pruning:        " << setw(6) << get_sec(end_spherical, end_prune) << " s " << get_msec(end_spherical, end_prune) % 1000 << " ms" << endl;
    else
      file << " ... for Grid Pruning:        " << setw(6) << get_sec(end_spherical, end_prune) << " s" << endl;
	else
		file << " ... for Grid Pruning:        " << setw(6) << get_msec(end_spherical, end_prune) << " ms" << endl;
	if (get_sec(end_prune, end_aspherical) > 1)
		if (get_sec(end_prune, end_aspherical) < 100)
      file << " ... for aspherical density:  " << setw(6) << get_sec(end_prune, end_aspherical) << " s " << get_msec(end_prune, end_aspherical) % 1000 << " ms" << endl;
    else
      file << " ... for aspherical density:  " << setw(6) << get_sec(end_prune, end_aspherical) << " s" << endl;
	else
		file << " ... for aspherical density:  " << setw(6) << get_msec(end_prune, end_aspherical) << " ms" << endl;
	if (get_sec(end_aspherical, before_kpts) > 1)
		if (get_sec(end_aspherical, before_kpts) < 100)
      file << " ... for density vectors:     " << setw(6) << get_sec(end_aspherical, before_kpts) << " s " << get_msec(end_aspherical, before_kpts) % 1000 << " ms" << endl;
    else
      file << " ... for density vectors:     " << setw(6) << get_sec(end_aspherical, before_kpts) << " s" << endl;
	else
		file << " ... for density vectors:     " << setw(6) << get_msec(end_aspherical, before_kpts) << " ms" << endl;
	if (get_sec(before_kpts, after_kpts) > 1)
		if (get_sec(before_kpts, after_kpts) < 100)
      file << " ... for k-points preparation:" << setw(6) << get_sec(before_kpts, after_kpts) << " s " << get_msec(before_kpts, after_kpts) % 1000 << " ms" << endl;
    else
      file << " ... for k-points preparation:" << setw(6) << get_sec(before_kpts, after_kpts) << " s" << endl;
	else
		file << " ... for k-points preparation:" << setw(6) << get_msec(before_kpts, after_kpts) << " ms" << endl;
	if (get_sec(after_kpts, end1) > 1)
		if (get_sec(after_kpts, end1) < 100)
      file << " ... for final preparation:   " << setw(6) << get_sec(after_kpts, end1) << " s " << get_msec(after_kpts, end1) % 1000 << " ms" << endl;
    else
      file << " ... for final preparation:   " << setw(6) << get_sec(after_kpts, end1) << " s" << endl;
	else
		file << " ... for final preparation:   " << setw(6) << get_msec(after_kpts, end1) << " ms" << endl;
	if (get_sec(end1, end) > 1)
		if (get_sec(end1, end) < 100)
      file << " ... for tsc calculation:     " << setw(6) << get_sec(end1, end) << " s " << get_msec(end1, end) % 1000 << " ms" << endl;
    else
      file << " ... for tsc calculation:     " << setw(6) << get_sec(end1, end) << " s" << endl;
	else
		file << " ... for tsc calculation:     " << setw(6) << get_msec(end1, end) << " ms" << endl;
}

inline int CountWords(const char *str)
{
	if (str == NULL)
		return -1;

	bool inSpaces = true;
	int numWords = 0;

	while (*str != '\0')
	{
		if (std::isspace(*str))
		{
			inSpaces = true;
		}
		else if (inSpaces)
		{
			numWords++;
			inSpaces = false;
		}

		++str;
	}

	return numWords;
};

inline bool exists(const std::string &name)
{
	if (FILE *file = fopen(name.c_str(), "r"))
	{
		fclose(file);
		return true;
	}
	else
	{
		return false;
	}
};

std::string atnr2letter(const int &nr);
void copy_file(std::string &from, std::string &to);
std::string shrink_string(std::string &input);
std::string shrink_string_to_atom(std::string &input, const int &atom_number);
std::string get_filename_from_path(const std::string &input);
std::string get_foldername_from_path(const std::string &input);
std::string get_basename_without_ending(const std::string &input);
//------------------Functions to read from .fchk files----------------------------------
bool read_fchk_integer_block(std::ifstream &in, std::string heading, ivec &result, bool rewind = true);
bool read_fchk_double_block(std::ifstream &in, std::string heading, vec &result, bool rewind = true);
int read_fchk_integer(std::string in);
int read_fchk_integer(std::ifstream &in, std::string search, bool rewind = true);
double read_fchk_double(std::string in);
double read_fchk_double(std::ifstream &in, std::string search, bool rewind = true);
//------------------Functions to work with configuration files--------------------------
void write_template_confi();
int program_confi(std::string &gaussian_path, std::string &turbomole_path,
				  std::string &basis, int &ncpus, double &mem, bool debug = false, bool expert = false, unsigned int counter = 0);
bool check_bohr(WFN &wave, bool debug);
int filetype_identifier(std::string &file, bool debug = false);

/*bool open_file_dialog(std::string &path, bool debug, std::vector <std::string> filter);
bool save_file_dialog(std::string &path, bool debug, const std::vector<std::string> &endings, const std::string &filename_given);
bool save_file_dialog(std::string &path, bool debug, const std::vector<std::string> &endings);*/
void select_cubes(std::vector<std::vector<unsigned int>> &selection, std::vector<WFN> &wavy, unsigned int nr_of_cubes = 1, bool wfnonly = false, bool debug = false);
bool unsaved_files(std::vector<WFN> &wavy);
int get_Z_from_label(const char *tmp);

inline int sum_of_bools(const std::vector<bool> in)
{
	int result = 0;
	for (int i = 0; i < in.size(); i++)
		if (in[i])
			result++;
	return result;
}

inline std::string trim(const std::string &s)
{
	if (s == "")
		return "";
	auto start = s.begin();
	while (start != s.end() && std::isspace(*start))
	{
		start++;
	}

	auto end = s.end();
	do
	{
		end--;
	} while (std::distance(start, end) > 0 && std::isspace(*end));

	return std::string(start, end + 1);
}

//-------------------------Progress_bar--------------------------------------------------

class progress_bar
{
	static const auto overhead = sizeof " [100%]";

	std::ostream &os;
	const std::size_t bar_width;
	std::string message;
	const std::string full_bar;
	const double precision;

public:
	progress_bar(std::ostream &os, std::size_t line_width,
				 std::string message_, const char symbol = '=', const double p = 0.05)
		: os{os},
		  bar_width{line_width - overhead},
		  message{std::move(message_)},
		  full_bar{std::string(bar_width, symbol) + std::string(bar_width, ' ')},
		  precision{p}
	{
		if (message.size() + 1 >= bar_width || message.find('\n') != message.npos)
		{
			os << message << '\n';
			message.clear();
		}
		else
		{
			message += ' ';
		}
		write(0.0);
	}

	// not copyable
	progress_bar(const progress_bar &) = delete;
	progress_bar &operator=(const progress_bar &) = delete;

	~progress_bar()
	{
		write(1.0);
		os << '\n';
	}

	void write(double fraction);
};
/*

class cosinus_annaeherung
{
public:
	cosinus_annaeherung();
	inline double get(double x) const
	{
		double xa = abs(x);
		size_t pos = static_cast<size_t>((xa * mSize) / MPI2); // Stueststelle bestimmen (Wird fuer grosse X ungenau, aber passt fuer x
		double dx = xa - pos * mStepwidth;
		pos = pos % mSize; // Modulo, da sinus periodisch ist.
		double y1 = mBase_values[pos];
		double y2 = mBase_values[pos + 1];
		return y1 + dx * (y2 - y1) / mStepwidth;
	}

	void   resize(size_t size);
	double calculate_error_at(double x) const;
private:
	size_t mSize;
	double* mBase_values;
	double mStepwidth;
};
struct sinus
{
	sinus(cosinus_annaeherung& helper) : helper(helper) {};
	double get(double x) { return helper.get(x - 1.57079632679489661922l); }
	cosinus_annaeherung& helper;
};

struct cosinus
{
	cosinus(cosinus_annaeherung& helper) : helper(helper) {};
	double get(double x) { return helper.get(x); }
	cosinus_annaeherung& helper;
};
*/
void readxyzMinMax_fromWFN(
	WFN &wavy,
	double *CoordMinMax,
	int *NbSteps,
	double Radius,
	double Increments,
	bool no_bohr = false);

void readxyzMinMax_fromCIF(
	std::string cif,
	double *CoordMinMax,
	int *NbSteps,
	std::vector<std::vector<double>> &cm,
	double Resolution,
	std::ofstream &file,
	bool debug = false);

void type2vector(
	const int &index,
	int *vector);

bool read_fracs_ADPs_from_CIF(std::string cif, WFN &wavy, cell &unit_cell, std::ofstream &log3, bool debug);

inline double double_from_string_with_esd(std::string in)
{
	if (in.find('(') == std::string::npos)
		return stod(in);
	else
		return stod(in.substr(0, in.find('(')));
}

void swap_sort(ivec order, cvec &v);

void swap_sort_multi(ivec order, std::vector<ivec> &v);

// Given a 3x3 matrix in a single array of double will find and sort eigenvalues and return biggest eigenvalue
double get_lambda_1(double *a);

double get_decimal_precision_from_CIF_number(std::string &given_string);

template <typename numtype = int>
struct hashFunction
{
	size_t operator()(const std::vector<numtype> &myVector) const
	{
		std::hash<numtype> hasher;
		size_t answer = 0;
		for (numtype i : myVector)
		{
			answer ^= hasher(i) + 0x9e3779b9 + (answer << 6) + (answer >> 2);
		}
		return answer;
	}
};

template <typename numtype = int>
struct hkl_equal
{
	bool operator()(const std::vector<numtype> &vec1, const std::vector<numtype> &vec2) const
	{
		const int size = vec1.size();
		if (size != vec2.size())
			return false;
		int similar = 0;
		for (int i = 0; i < size; i++)
		{
			if (vec1[i] == vec2[i])
				similar++;
			else if (vec1[i] == -vec2[i])
				similar--;
		}
		if (abs(similar) == size)
			return true;
		else
			return false;
	}
};

template <typename numtype = int>
struct hkl_less
{
	bool operator()(const std::vector<numtype> &vec1, const std::vector<numtype> &vec2) const
	{
		if (vec1[0] < vec2[0])
		{
			return true;
		}
		else if (vec1[0] == vec2[0])
		{
			if (vec1[1] < vec2[1])
			{
				return true;
			}
			else if (vec1[1] == vec2[1])
			{
				if (vec1[2] < vec2[2])
				{
					return true;
				}
				else
					return false;
			}
			else
				return false;
		}
		else
			return false;
	}
};

inline unsigned int doublefactorial(int n)
{
	if (n <= 1)
		return 1;
	return n * doublefactorial(n - 2);
}

struct primitive
{
	int center, type;
	double exp, coefficient;
	double norm_const = -10;
	void normalize()
	{
		coefficient *= normalization_constant();
	};
	void unnormalize()
	{
		coefficient /= normalization_constant();
	};
	double normalization_constant()
	{
		// assuming type is equal to angular momentum
		return norm_const;
	}
	primitive() : center(0), type(0), exp(0.0), coefficient(0.0) {}
	primitive(int c, int t, double e, double coef) : center(c), type(t), exp(e), coefficient(coef)
	{
		norm_const = pow(
			pow(2, 7 + 4 * type) * pow(exp, 3 + 2 * type) / constants::PI / pow(doublefactorial(2 * type + 1), 2),
			0.25);
	}
};

struct ECP_primitive : primitive
{
	int n;
	ECP_primitive() : primitive(), n(0) {}
	ECP_primitive(int c, int t, double e, double coef, int n) : primitive(c, t, e, coef), n(n) {}
};

struct tonto_primitive
{
	int center, type;
	double exp, coefficient;
	double norm_const = -10;
	void normalize()
	{
		coefficient *= normalization_constant();
	};
	void unnormalize()
	{
		coefficient /= normalization_constant();
	};
	double normalization_constant()
	{
		// assuming type is equal to angular momentum
		return norm_const;
	}
	tonto_primitive() : center(0), type(0), exp(0.0), coefficient(0.0) {}
	tonto_primitive(int c, int t, double e, double coef) : center(c), type(t), exp(e), coefficient(coef)
	{
		norm_const = pow(constants::PI, -0.75) * pow(2.0, type + 0.75) * pow(exp, type * 0.5 + 0.75) / sqrt(doublefactorial(type));
	}
};

typedef std::set<ivec> hkl_list;
typedef std::set<ivec>::const_iterator hkl_list_it;

typedef std::set<vec> hkl_list_d;
typedef std::set<vec>::const_iterator hkl_list_it_d;

//---------------- Object for handling all input options -------------------------------
struct options
{
	options() : log_file(std::cout)
	{
		groups.resize(1);
	};
	options(int &argc, char **argv, std::ostream &log) : log_file(log)
	{
		groups.resize(1);
		look_for_debug(argc, argv);
	};
	std::ostream &log_file;
	void look_for_debug(int &argc, char **argv);
	void digest_options();
	int accuracy = 2;
	int threads = -1;
	int pbc = 0;
	double resolution = 0.1;
	double radius = 2.0;
	bool becke = false;
	bool electron_diffraction = false;
	bool ECP = false;
	bool set_ECPs = false;
	int ECP_mode = 0;
	bool calc = false;
	bool eli = false;
	bool esp = false;
	bool elf = false;
	bool lap = false;
	bool rdg = false;
	bool hdef = false;
	bool def = false;
	bool fract = false;
	bool hirsh = false;
	bool s_rho = false;
	bool SALTED = false, SALTED_BECKE = false;
	bool Olex2_1_3_switch = false;
	bool iam_switch = false;
	bool read_k_pts = false;
	bool save_k_pts = false;
	bool combined_tsc_calc = false;
	bool binary_tsc = true;
	bool cif_based_combined_tsc_calc = false;
	bool no_date = false;
	bool gbw2wfn = false;
	bool old_tsc = false;
	bool thakkar_d_plot = false;
	bool write_CIF = false;
	bool ML_test2 = false;
	double d_sfac_scan = 0.0;
	double sfac_diffuse = 0.0;
	double dmin = 99.0;
	int hirsh_number = 0;
	double MinMax[6]{0, 0, 0, 0, 0, 0};
	int NbSteps[3]{0, 0, 0};
	ivec MOs;
	std::vector<ivec> groups;
	std::vector<vec> twin_law;
	std::vector<ivec> combined_tsc_groups;
	bool all_mos = false;
	bool test = false;
	std::string wfn;
	std::string wfn2;
	std::string fchk;
	std::string basis_set;
	std::string hkl;
	std::string cif;
	std::string method;
	std::string xyz_file;
	std::string coef_file;
	std::string fract_name;
	std::vector<std::string> combined_tsc_calc_files;
	std::vector<std::string> combined_tsc_calc_cifs;
	std::vector<unsigned int> combined_tsc_calc_mult;
	std::vector<int> combined_tsc_calc_charge;
	std::string wavename;
	std::string gaussian_path;
	std::string turbomole_path;
	std::string basis_set_path;
	std::vector<std::string> arguments;
	std::vector<std::string> combine_mo;
	std::vector<std::string> Cations;
	std::vector<std::string> Anions;
	ivec cmo1;
	ivec cmo2;
	ivec ECP_nrs;
	ivec ECP_elcounts;
	int ncpus = 0;
	double mem = 0.0;
	unsigned int mult = 0;
	int charge = 0;
	bool debug = false;
	hkl_list m_hkl_list;

	options(int accuracy, 
		int threads, 
		int pbc, 
		double resolution, 
		double radius, 
		bool becke, 
		bool electron_diffraction, 
		bool ECP, 
		bool set_ECPs, 
		int ECP_mode, 
		bool calc, 
		bool eli, 
		bool esp, 
		bool elf, 
		bool lap, 
		bool rdg, 
		bool hdef, 
		bool def, 
		bool fract, 
		bool hirsh, 
		bool s_rho, 
		bool SALTED, 
		bool SALTED_BECKE, bool Olex2_1_3_switch, bool iam_switch, bool read_k_pts, bool save_k_pts, bool combined_tsc_calc, bool binary_tsc, bool cif_based_combined_tsc_calc, bool density_test_cube, bool no_date, bool gbw2wfn, bool old_tsc, bool thakkar_d_plot, bool spherical_harmonic, bool ML_test, double sfac_scan, double sfac_diffuse, double dmin, int hirsh_number, const ivec &MOs, const std::vector<ivec> &groups, const std::vector<vec> &twin_law, const std::vector<ivec> &combined_tsc_groups, bool all_mos, bool test, const std::string &wfn, const std::string &fchk, const std::string &basis_set, const std::string &hkl, const std::string &cif, const std::string &method, const std::string &xyz_file, const std::string &coef_file, const std::string &fract_name, const std::vector<std::string> &combined_tsc_calc_files, const std::vector<std::string> &combined_tsc_calc_cifs, const std::string &wavename, const std::string &gaussian_path, const std::string &turbomole_path, const std::string &basis_set_path, const std::vector<std::string> &arguments, const std::vector<std::string> &combine_mo, const std::vector<std::string> &Cations, const std::vector<std::string> &Anions, const ivec &cmo1, const ivec &cmo2, const ivec &ECP_nrs, const ivec &ECP_elcounts, int ncpus, double mem, unsigned int mult, bool debug, bool ML_test2, const hkl_list &m_hkl_list, std::ostream log)
		: accuracy(accuracy), threads(threads), pbc(pbc), resolution(resolution), radius(radius), becke(becke), electron_diffraction(electron_diffraction), ECP(ECP), set_ECPs(set_ECPs), ECP_mode(ECP_mode), calc(calc), eli(eli), esp(esp), elf(elf), lap(lap), rdg(rdg), hdef(hdef), def(def), fract(fract), hirsh(hirsh), s_rho(s_rho), SALTED(SALTED), SALTED_BECKE(SALTED_BECKE), Olex2_1_3_switch(Olex2_1_3_switch), iam_switch(iam_switch), read_k_pts(read_k_pts), save_k_pts(save_k_pts), combined_tsc_calc(combined_tsc_calc), binary_tsc(binary_tsc), cif_based_combined_tsc_calc(cif_based_combined_tsc_calc), no_date(no_date), gbw2wfn(gbw2wfn), old_tsc(old_tsc), thakkar_d_plot(thakkar_d_plot), d_sfac_scan(sfac_scan), sfac_diffuse(sfac_diffuse), dmin(dmin), hirsh_number(hirsh_number), MOs(MOs), groups(groups), twin_law(twin_law), combined_tsc_groups(combined_tsc_groups), all_mos(all_mos), test(test), wfn(wfn), fchk(fchk), basis_set(basis_set), hkl(hkl), cif(cif), method(method), xyz_file(xyz_file), coef_file(coef_file), fract_name(fract_name), combined_tsc_calc_files(combined_tsc_calc_files), combined_tsc_calc_cifs(combined_tsc_calc_cifs), wavename(wavename), gaussian_path(gaussian_path), turbomole_path(turbomole_path), basis_set_path(basis_set_path), arguments(arguments), combine_mo(combine_mo), Cations(Cations), Anions(Anions), cmo1(cmo1), cmo2(cmo2), ECP_nrs(ECP_nrs), ECP_elcounts(ECP_elcounts), ncpus(ncpus), mem(mem), mult(mult), debug(debug), m_hkl_list(m_hkl_list), log_file(log), ML_test2(ML_test2)
	{
	}
};

const double gaussian_radial(primitive &p, double &r);

const double calc_density_ML(double &x,
							 double &y,
							 double &z,
							 vec &coefficients,
							 std::vector<atom> &atoms,
							 const int &exp_coefs);
const double calc_density_ML(double &x,
							 double &y,
							 double &z,
							 vec &coefficients,
							 std::vector<atom> &atoms,
							 const int &exp_coefs,
							 const int &atom_nr);

int load_basis_into_WFN(WFN &wavy, const std::vector<std::vector<primitive>> &b);

inline const std::vector<std::vector<primitive>> TZVP_JKfit(
	{
		{
			// center  type  exp         coef
			{0, 0, 9.5302493327, 1.0},
			{0, 0, 1.9174506246, 1.0},
			{0, 0, 0.68424049142, 1.0},
			{0, 0, 0.28413255710, 1.0},
			{0, 1, 2.9133232035, 1.0},
			{0, 1, 1.2621205398, 1.0},
			{0, 1, 0.50199775874, 1.0},
			{0, 2, 2.3135329149, 1.0},
			{0, 2, 0.71290724024, 1.0},
			{0, 3, 1.6565726132, 1.0},
		},	// H
		{}, // He
		{}, // Li
		{}, // Be
		{}, // B
		{
			{0, 0, 1113.9867719, 1.0},
			{0, 0, 369.16234180, 1.0},
			{0, 0, 121.79275232, 1.0},
			{0, 0, 48.127114540, 1.0},
			{0, 0, 20.365074004, 1.0},
			{0, 0, 8.0883596856, 1.0},
			{0, 0, 2.5068656570, 1.0},
			{0, 0, 1.2438537380, 1.0},
			{0, 0, 0.48449899601, 1.0},
			{0, 0, 0.19185160296, 1.0},
			{0, 1, 102.99176249, 1.0},
			{0, 1, 28.132594009, 1.0},
			{0, 1, 9.8364318173, 1.0},
			{0, 1, 3.3490544980, 1.0},
			{0, 1, 1.4947618613, 1.0},
			{0, 1, 0.57690108899, 1.0},
			{0, 1, 0.20320063291, 1.0},
			{0, 2, 10.594068356, 1.0},
			{0, 2, 3.5997195366, 1.0},
			{0, 2, 1.3355691094, 1.0},
			{0, 2, 0.51949764954, 1.0},
			{0, 2, 0.19954125200, 1.0},
			{0, 3, 1.194866338369, 1.0},
			{0, 3, .415866338369, 1.0},
			{0, 4, .858866338369, 1.0}}, // C
		{
			{0, 0, 1102.8622453, 1.0},
			{0, 0, 370.98041153, 1.0},
			{0, 0, 136.73555938, 1.0},
			{0, 0, 50.755871924, 1.0},
			{0, 0, 20.535656095, 1.0},
			{0, 0, 7.8318737184, 1.0},
			{0, 0, 3.4784063855, 1.0},
			{0, 0, 1.4552856603, 1.0},
			{0, 0, 0.63068989071, 1.0},
			{0, 0, 0.27276596483, 1.0},
			{0, 1, 93.540954073, 1.0},
			{0, 1, 29.524019527, 1.0},
			{0, 1, 10.917502987, 1.0},
			{0, 1, 4.3449288991, 1.0},
			{0, 1, 1.8216912640, 1.0},
			{0, 1, 0.75792424494, 1.0},
			{0, 1, 0.28241469033, 1.0},
			{0, 2, 16.419378926, 1.0},
			{0, 2, 5.0104049385, 1.0},
			{0, 2, 1.9793971884, 1.0},
			{0, 2, 0.78495771518, 1.0},
			{0, 2, 0.28954065963, 1.0},
			{0, 3, 1.79354239843, 1.0},
			{0, 3, .60854239843, 1.0},
			{0, 4, 1.23254239843, 1.0}}, // N
		{
			{0, 0, 1517.8667506, 1.0},
			{0, 0, 489.67952008, 1.0},
			{0, 0, 176.72118665, 1.0},
			{0, 0, 63.792233137, 1.0},
			{0, 0, 25.366499130, 1.0},
			{0, 0, 9.9135491200, 1.0},
			{0, 0, 4.4645306584, 1.0},
			{0, 0, 1.8017743661, 1.0},
			{0, 0, 0.80789710965, 1.0},
			{0, 0, 0.33864326862, 1.0},
			{0, 1, 120.16030921, 1.0},
			{0, 1, 34.409622474, 1.0},
			{0, 1, 12.581148610, 1.0},
			{0, 1, 5.0663824249, 1.0},
			{0, 1, 2.0346927092, 1.0},
			{0, 1, 0.86092967212, 1.0},
			{0, 1, 0.36681356726, 1.0},
			{0, 2, 19.043062805, 1.0},
			{0, 2, 5.8060381104, 1.0},
			{0, 2, 2.1891841580, 1.0},
			{0, 2, 0.87794613558, 1.0},
			{0, 2, 0.35623646700, 1.0},
			{0, 3, 2.493914788135, 1.0},
			{0, 3, .824914788135, 1.0},
			{0, 4, 1.607914788135, 1.0}}, // O
									
		{},								  // F
		{}								  // Ne
	});

inline const std::vector<std::vector<primitive>> QZVP_JKfit(
	{
		{// center  type  exp         coef
		 {0, 0, 9.5302493327, 1.0},
		 {0, 0, 1.9174506246, 1.0},
		 {0, 0, 0.68424049142, 1.0},
		 {0, 0, 0.28413255710, 1.0},
		 {0, 1, 2.9133232035, 1.0},
		 {0, 1, 1.2621205398, 1.0},
		 {0, 1, 0.50199775874, 1.0},
		 {0, 2, 2.8832083931, 1.0},
		 {0, 2, 1.2801701725, 1.0},
		 {0, 2, 0.52511317770, 1.0},
		 {0, 3, 2.7489448439, 1.0},
		 {0, 3, 1.1900885456, 1.0},
		 {0, 4, 1.4752662714, 1.0}}, // H
		{},							 // He
		{},							 // Li
		{},							 // Be
		{},							 // B
		{
			{0, 0, 1113.9867719, 1.0},
			{0, 0, 369.16234180, 1.0},
			{0, 0, 121.79275232, 1.0},
			{0, 0, 48.127114540, 1.0},
			{0, 0, 20.365074004, 1.0},
			{0, 0, 8.0883596856, 1.0},
			{0, 0, 2.5068656570, 1.0},
			{0, 0, 1.2438537380, 1.0},
			{0, 0, 0.48449899601, 1.0},
			{0, 0, 0.19185160296, 1.0},
			{0, 1, 102.99176249, 1.0},
			{0, 1, 28.132594009, 1.0},
			{0, 1, 9.8364318173, 1.0},
			{0, 1, 3.3490544980, 1.0},
			{0, 1, 1.4947618613, 1.0},
			{0, 1, 0.57690108899, 1.0},
			{0, 1, 0.20320063291, 1.0},
			{0, 2, 10.594068356, 1.0},
			{0, 2, 3.5997195366, 1.0},
			{0, 2, 1.3355691094, 1.0},
			{0, 2, 0.51949764954, 1.0},
			{0, 2, 0.19954125200, 1.0},
			{0, 3, 1.95390000000, 1.0},
			{0, 3, .75490000000, 1.0},
			{0, 3, .33390000000, 1.0},
			{0, 4, 1.52490000000, 1.0},
			{0, 4, .59090000000, 1.0},
			{0, 5, 1.11690000000, 1.0}}, // C
		{
			{0, 0, 1102.8622453, 1.0},
			{0, 0, 370.98041153, 1.0},
			{0, 0, 136.73555938, 1.0},
			{0, 0, 50.755871924, 1.0},
			{0, 0, 20.535656095, 1.0},
			{0, 0, 7.8318737184, 1.0},
			{0, 0, 3.4784063855, 1.0},
			{0, 0, 1.4552856603, 1.0},
			{0, 0, 0.63068989071, 1.0},
			{0, 0, 0.27276596483, 1.0},
			{0, 1, 93.540954073, 1.0},
			{0, 1, 29.524019527, 1.0},
			{0, 1, 10.917502987, 1.0},
			{0, 1, 4.3449288991, 1.0},
			{0, 1, 1.8216912640, 1.0},
			{0, 1, 0.75792424494, 1.0},
			{0, 1, 0.28241469033, 1.0},
			{0, 2, 16.419378926, 1.0},
			{0, 2, 5.0104049385, 1.0},
			{0, 2, 1.9793971884, 1.0},
			{0, 2, 0.78495771518, 1.0},
			{0, 2, 0.28954065963, 1.0},
			{0, 3, 2.98600000000, 1.0},
			{0, 3, 1.11700000000, 1.0},
			{0, 3, .48400000000, 1.0},
			{0, 4, 2.17600000000, 1.0},
			{0, 4, .83400000000, 1.0},
			{0, 5, 1.57600000000, 1.0}}, // N
		{
			{0, 0, 1517.8667506, 1.0},
			{0, 0, 489.67952008, 1.0},
			{0, 0, 176.72118665, 1.0},
			{0, 0, 63.792233137, 1.0},
			{0, 0, 25.366499130, 1.0},
			{0, 0, 9.9135491200, 1.0},
			{0, 0, 4.4645306584, 1.0},
			{0, 0, 1.8017743661, 1.0},
			{0, 0, 0.80789710965, 1.0},
			{0, 0, 0.33864326862, 1.0},
			{0, 1, 120.16030921, 1.0},
			{0, 1, 34.409622474, 1.0},
			{0, 1, 12.581148610, 1.0},
			{0, 1, 5.0663824249, 1.0},
			{0, 1, 2.0346927092, 1.0},
			{0, 1, 0.86092967212, 1.0},
			{0, 1, 0.36681356726, 1.0},
			{0, 2, 19.043062805, 1.0},
			{0, 2, 5.8060381104, 1.0},
			{0, 2, 2.1891841580, 1.0},
			{0, 2, 0.87794613558, 1.0},
			{0, 2, 0.35623646700, 1.0},
			{0, 3, 3.96585000000, 1.0},
			{0, 3, 1.49085000000, 1.0},
			{0, 3, .63485000000, 1.0},
			{0, 4, 2.85685000000, 1.0},
			{0, 4, 1.04985000000, 1.0},
			{0, 5, 2.03685000000, 1.0}}, // O
		{},								 // F
		{},								 // Ne
		{},							     // Na
		{},								// Mg
		{},								// Al
		{},                             // Si
		{},								// P
		{
			{0, 0, 6402.4580816, 1.0},
			{0, 0, 1704.9249408, 1.0},
			{0, 0, 598.09038904, 1.0},
			{0, 0, 225.90920364, 1.0},
			{0, 0, 94.426261909, 1.0},
			{0, 0, 37.222420840, 1.0},
			{0, 0, 18.259106532, 1.0},
			{0, 0, 8.7937857291, 1.0},
			{0, 0, 4.4619424717, 1.0},
			{0, 0, 2.2735987703, 1.0},
			{0, 0, 0.78075394724, 1.0},
			{0, 0, 0.48191153763, 1.0},
			{0, 0, 0.21573555055, 1.0},
			{0, 1, 796.16481618, 1.0},
			{0, 1, 226.17928266, 1.0},
			{0, 1, 85.598645757, 1.0},
			{0, 1, 36.268978384, 1.0},
			{0, 1, 16.325207000, 1.0},
			{0, 1, 7.6145421568, 1.0},
			{0, 1, 3.3388129956, 1.0},
			{0, 1, 1.8665221526, 1.0},
			{0, 1, 0.94025159837, 1.0},
			{0, 1, 0.51521255747, 1.0},
			{0, 1, 0.26318380469, 1.0},
			{0, 2, 190.55658672, 1.0},
			{0, 2, 60.346801430, 1.0},
			{0, 2, 25.141967063, 1.0},
			{0, 2, 12.146768196, 1.0},
			{0, 2, 6.2993488132, 1.0},
			{0, 2, 3.0706106656, 1.0},
			{0, 2, 1.2911977798, 1.0},
			{0, 2, 0.56678602865, 1.0},
			{0, 2, 0.24192419172, 1.0},
			{0, 3, 3.30961330000, 1.0},
			{0, 3, 1.35367000000, 1.0},
			{0, 3, 0.60767000000, 1.0},
			{0, 3, 0.30667000000, 1.0},
			{0, 4, 0.97267000000, 1.0},
			{0, 4, 0.43867000000, 1.0},
			{0, 5, 0.78667000000, 1.0}}, //S
		{}								// Cl
	});

inline const std::vector<std::vector<primitive>> def2_SVP_JKFIT(
	{
		{
		{0, 0, 22.068343      , 0.1},
		{0, 0, 0.2717874      , 1.0},
		{0, 1, 1.8529979      , 1.0},
		{0, 1, 0.3881034      , 1.0},
		{0, 2, 2.5579933      , 1.0},
		{0, 2, 0.32926490000  , 1.0}},		//H
		{
		{0, 0, 66.205029      , 0.1},
		{0, 0, 13.1717136     , 0.4},
		{0, 0, 3.1622361      , 0.9},
		{0, 0, 0.8153622      , 1.0},
		{0, 1, 9.612881982    , 1.0},
		{0, 1, 3.20429397     , 1.0},
		{0, 1, 1.1643102      , 1.0},
		{0, 2, 7.6739799      , 1.0},
		{0, 2, 0.98779470000  , 1.0}},		//He
		{
		{0, 0, 188.5699865    , 0.1},
		{0, 0, 2.803303       , 0.4},
		{0, 0, 1.213845       , 1.0},
		{0, 0, 0.5608472      , 1.0},
		{0, 0, 0.2721716      , 1.0},
		{0, 0, 0.1361922      , 1.0},
		{0, 0, 0.0688638      , 1.0},
		{0, 1, 9.0816823      , 0.3},
		{0, 1, 0.6337328      , 0.4},
		{0, 1, 0.3015717      , -0.1},
		{0, 1, 0.1491141      , 1.0},
		{0, 1, 0.0747033      , 1.0},
		{0, 2, 1.3650595      , 0.6},
		{0, 2, 0.3388178      , -0.2},
		{0, 2, 0.1702257      , 0.3},
		{0, 3, 0.6059012      , 1.0},
		{0, 3, 0.17311460000  , 1.0}},		//Li
		{
		{0, 0, 145.425355     , 0.2},
		{0, 0, 2.3772683      , 1.0},
		{0, 0, 1.109202       , 1.0},
		{0, 0, 0.562724       , 1.0},
		{0, 0, 0.3044705      , 1.0},
		{0, 0, 0.1715858      , 1.0},
		{0, 0, 0.0980565      , 1.0},
		{0, 1, 13.9643253     , 0.4},
		{0, 1, 0.7198232      , 1.0},
		{0, 1, 0.329289       , 1.0},
		{0, 1, 0.1595195      , 1.0},
		{0, 1, 0.0788265      , 1.0},
		{0, 2, 2.4566484      , 0.1},
		{0, 2, 0.2962418      , 1.0},
		{0, 2, 0.1028739      , 1.0},
		{0, 3, 0.709008       , 1.0},
		{0, 3, 0.20257371400  , 1.0}},		//Be
		{
		{0, 0, 305.0512594    , 0.2},
		{0, 0, 8.7025338      , 0.3},
		{0, 0, 4.7737916      , 0.1},
		{0, 0, 2.6974657      , 1.0},
		{0, 0, 1.5640429      , 1.0},
		{0, 0, 0.9266701      , 1.0},
		{0, 0, 0.5585253      , 1.0},
		{0, 0, 0.3408434      , 1.0},
		{0, 0, 0.2095733      , 1.0},
		{0, 0, 0.129185       , 1.0},
		{0, 1, 94.469984      , 0.1},
		{0, 1, 6.7457081      , 0.7},
		{0, 1, 3.1084655      , 0.5},
		{0, 1, 1.4682066      , 1.0},
		{0, 1, 0.7130179      , 1.0},
		{0, 1, 0.353384       , 1.0},
		{0, 1, 0.1773294      , 1.0},
		{0, 1, 0.0902499      , 1.0},
		{0, 2, 5.366408       , 0.4},
		{0, 2, 0.8743697      , 1.0},
		{0, 2, 0.4079433      , 1.0},
		{0, 2, 0.201544       , 1.0},
		{0, 2, 0.1015676      , 1.0},
		{0, 3, 1.5043374      , 0.3},
		{0, 4, 0.87665570000  , 0.7}},		//B
		{
		{0, 0, 384.4438241    , 0.2},
		{0, 0, 9.5470634      , 0.2},
		{0, 0, 5.1584143      , 1.0},
		{0, 0, 2.8816701      , 1.0},
		{0, 0, 1.6573522      , 1.0},
		{0, 0, 0.9768102      , 1.0},
		{0, 0, 0.58702        , 1.0},
		{0, 0, 0.3577927      , 1.0},
		{0, 0, 0.219955       , 1.0},
		{0, 0, 0.1356077      , 1.0},
		{0, 1, 62.0679561     , 0.3},
		{0, 1, 5.6367045      , 0.6},
		{0, 1, 2.8744918      , 1.0},
		{0, 1, 1.4513791      , 1.0},
		{0, 1, 0.7328327      , 1.0},
		{0, 1, 0.3700256      , 1.0},
		{0, 1, 0.186836       , 1.0},
		{0, 1, 0.0952821      , 1.0},
		{0, 2, 11.2820057     , 0.1},
		{0, 2, 1.7517699      , 0.9},
		{0, 2, 0.7164852      , 1.0},
		{0, 2, 0.2970402      , 1.0},
		{0, 2, 0.1237087      , 1.0},
		{0, 3, 1.8786227      , 0.6},
		{0, 4, 1.14306020000  , 1.0}},		//C
		{
		{0, 0, 502.8608492    , 0.2},
		{0, 0, 10.5912894     , 0.1},
		{0, 0, 5.6186543      , 1.0},
		{0, 0, 3.0952689      , 1.0},
		{0, 0, 1.7624885      , 1.0},
		{0, 0, 1.0319737      , 1.0},
		{0, 0, 0.6178352      , 1.0},
		{0, 0, 0.3759343      , 1.0},
		{0, 0, 0.2310083      , 1.0},
		{0, 0, 0.1424224      , 1.0},
		{0, 1, 44.6074972     , 0.5},
		{0, 1, 5.7421344      , 0.2},
		{0, 1, 2.8993384      , 1.0},
		{0, 1, 1.4639486      , 1.0},
		{0, 1, 0.7391861      , 1.0},
		{0, 1, 0.3732351      , 1.0},
		{0, 1, 0.1884567      , 1.0},
		{0, 1, 0.0951571      , 1.0},
		{0, 2, 46.8440406     , 0.0},
		{0, 2, 2.1154199      , 1.0},
		{0, 2, 0.8923218      , 1.0},
		{0, 2, 0.3848043      , 1.0},
		{0, 3, 2.4229788      , 0.7},
		{0, 3, 0.8335441      , 0.7},
		{0, 4, 1.45867910000  , 1.0}},		//N
		{
		{0, 0, 625.2829811    , 0.2},
		{0, 0, 11.8077591     , 1.0},
		{0, 0, 6.1827814      , 1.0},
		{0, 0, 3.3709061      , 1.0},
		{0, 0, 1.9042805      , 1.0},
		{0, 0, 1.1085447      , 1.0},
		{0, 0, 0.6609886      , 1.0},
		{0, 0, 0.4010814      , 1.0},
		{0, 0, 0.2459769      , 1.0},
		{0, 0, 0.1513939      , 1.0},
		{0, 1, 77.6874838     , 0.4},
		{0, 1, 5.4848863      , 1.0},
		{0, 1, 2.9732983      , 1.0},
		{0, 1, 1.473526       , 1.0},
		{0, 1, 0.7360341      , 1.0},
		{0, 1, 0.3697414      , 1.0},
		{0, 1, 0.1863721      , 1.0},
		{0, 1, 0.0949906      , 1.0},
		{0, 2, 37.7071074     , 0.1},
		{0, 2, 2.3304365      , 1.0},
		{0, 2, 0.9328267      , 1.0},
		{0, 2, 0.3739285      , 1.0},
		{0, 3, 3.0293422      , 0.8},
		{0, 3, 0.924849       , 0.6},
		{0, 4, 1.69348090000  , 1.0}},		//O
		{
		{0, 0, 858.4098655    , 0.2},
		{0, 0, 12.9733047     , 1.0},
		{0, 0, 6.6262885      , 1.0},
		{0, 0, 3.5457034      , 1.0},
		{0, 0, 1.9769882      , 1.0},
		{0, 0, 1.1415471      , 1.0},
		{0, 0, 0.6779147      , 1.0},
		{0, 0, 0.4109404      , 1.0},
		{0, 0, 0.2522467      , 1.0},
		{0, 0, 0.1554898      , 1.0},
		{0, 1, 146.4330804    , 0.2},
		{0, 1, 10.5795792     , 0.5},
		{0, 1, 5.4280075      , 1.0},
		{0, 1, 2.9079399      , 1.0},
		{0, 1, 1.4550833      , 1.0},
		{0, 1, 0.7315239      , 1.0},
		{0, 1, 0.368805       , 1.0},
		{0, 1, 0.1861123      , 1.0},
		{0, 2, 42.931821      , 0.1},
		{0, 2, 2.5066049      , 1.0},
		{0, 2, 0.9724219      , 1.0},
		{0, 2, 0.3772451      , 1.0},
		{0, 3, 3.4749889      , 0.7},
		{0, 3, 1.2194812      , 0.7},
		{0, 4, 2.04590910000  , 1.0}},		//F
		{
		{0, 0, 1073.0123319   , 0.2},
		{0, 0, 16.216630875   , 1.0},
		{0, 0, 8.282860625    , 1.0},
		{0, 0, 4.43212925     , 1.0},
		{0, 0, 2.47123525     , 1.0},
		{0, 0, 1.426933875    , 1.0},
		{0, 0, 0.847393375    , 1.0},
		{0, 0, 0.5136755      , 1.0},
		{0, 0, 0.315308375    , 1.0},
		{0, 0, 0.19436225     , 1.0},
		{0, 1, 183.0413505    , 0.2},
		{0, 1, 13.224474      , 0.5},
		{0, 1, 6.785009375    , 1.0},
		{0, 1, 3.634924875    , 1.0},
		{0, 1, 1.818854125    , 1.0},
		{0, 1, 0.914404875    , 1.0},
		{0, 1, 0.46100625     , 1.0},
		{0, 1, 0.232640375    , 1.0},
		{0, 2, 53.66477625    , 0.1},
		{0, 2, 3.133256125    , 1.0},
		{0, 2, 1.215527375    , 1.0},
		{0, 2, 0.471556375    , 1.0},
		{0, 3, 4.343736125    , 0.7},
		{0, 3, 1.5243515      , 0.7},
		{0, 4, 2.55738637500  , 1.0}},		//Ne
		{
		{0, 0, 1671.8975852   , 0.1},
		{0, 0, 36.2193848     , 1.0},
		{0, 0, 16.0730137     , 1.0},
		{0, 0, 7.5289891      , 1.0},
		{0, 0, 3.7112547      , 1.0},
		{0, 0, 1.9179111      , 1.0},
		{0, 0, 1.0345305      , 1.0},
		{0, 0, 0.5794957      , 1.0},
		{0, 0, 0.3351543      , 1.0},
		{0, 0, 0.1988643      , 1.0},
		{0, 0, 0.120223       , 1.0},
		{0, 0, 0.0735129      , 1.0},
		{0, 0, 0.0451232      , 1.0},
		{0, 1, 349.7521385    , 0.1},
		{0, 1, 10.4662639     , 1.0},
		{0, 1, 4.7731711      , 1.0},
		{0, 1, 2.240985       , 1.0},
		{0, 1, 1.0781356      , 1.0},
		{0, 1, 0.5288554      , 1.0},
		{0, 1, 0.2631013      , 1.0},
		{0, 1, 0.1320154      , 1.0},
		{0, 1, 0.0664311      , 1.0},
		{0, 2, 83.7159358     , 0.1},
		{0, 2, 6.0067198      , 0.8},
		{0, 2, 2.6999775      , 1.0},
		{0, 2, 1.248474       , 1.0},
		{0, 2, 0.5893866      , 1.0},
		{0, 2, 0.2817818      , 1.0},
		{0, 2, 0.1352918      , 1.0},
		{0, 3, 15.819398      , 1.0},
		{0, 3, 5.2731318      , 1.0},
		{0, 3, 1.7577169      , 1.0},
		{0, 3, 0.58590333     , 1.0},
		{0, 4, 3.60581590000  , 0.7}},		//Na
		{
		{0, 0, 2128.2495853   , 0.1},
		{0, 0, 60.8157036     , 0.5},
		{0, 0, 28.4040429     , 1.0},
		{0, 0, 13.8997073     , 1.0},
		{0, 0, 7.1055414      , 1.0},
		{0, 0, 3.7812917      , 1.0},
		{0, 0, 2.0863957      , 1.0},
		{0, 0, 1.1882511      , 1.0},
		{0, 0, 0.6950359      , 1.0},
		{0, 0, 0.4152765      , 1.0},
		{0, 0, 0.2519901      , 1.0},
		{0, 0, 0.1543511      , 1.0},
		{0, 0, 0.0948428      , 1.0},
		{0, 1, 396.4815262    , 0.1},
		{0, 1, 11.2223478     , 1.0},
		{0, 1, 5.086089       , 1.0},
		{0, 1, 2.3792166      , 1.0},
		{0, 1, 1.1430544      , 1.0},
		{0, 1, 0.560954       , 1.0},
		{0, 1, 0.2795766      , 1.0},
		{0, 1, 0.1406542      , 1.0},
		{0, 1, 0.070986       , 1.0},
		{0, 2, 101.5710428    , 0.1},
		{0, 2, 7.8352029      , 1.0},
		{0, 2, 3.4489951      , 1.0},
		{0, 2, 1.5360109      , 1.0},
		{0, 2, 0.6898367      , 1.0},
		{0, 2, 0.311389       , 1.0},
		{0, 2, 0.1407984      , 1.0},
		{0, 3, 13.3708335     , 1.0},
		{0, 3, 4.45694454     , 1.0},
		{0, 3, 1.4856481      , 1.0},
		{0, 3, 0.49521606     , 1.0},
		{0, 4, 2.68877660000  , 0.8}},		//Mg
		{
		{0, 0, 3311.4623116   , 0.1},
		{0, 0, 67.3569565     , 0.5},
		{0, 0, 30.7115064     , 1.0},
		{0, 0, 14.8515319     , 1.0},
		{0, 0, 7.5180147      , 1.0},
		{0, 0, 3.9696099      , 1.0},
		{0, 0, 2.1772745      , 1.0},
		{0, 0, 1.2347125      , 1.0},
		{0, 0, 0.7201817      , 1.0},
		{0, 0, 0.429605       , 1.0},
		{0, 0, 0.2604941      , 1.0},
		{0, 0, 0.1611259      , 1.0},
		{0, 0, 0.1004982      , 1.0},
		{0, 1, 440.453628     , 0.1},
		{0, 1, 66.0333518     , 0.9},
		{0, 1, 27.4313958     , 1.0},
		{0, 1, 11.8880832     , 1.0},
		{0, 1, 5.3525289      , 1.0},
		{0, 1, 2.4920135      , 1.0},
		{0, 1, 1.1934791      , 1.0},
		{0, 1, 0.584616       , 1.0},
		{0, 1, 0.2911107      , 1.0},
		{0, 1, 0.1464136      , 1.0},
		{0, 1, 0.0738856      , 1.0},
		{0, 2, 100.7280583    , 0.1},
		{0, 2, 15.275436      , 0.9},
		{0, 2, 6.4821533      , 1.0},
		{0, 2, 2.8906045      , 1.0},
		{0, 2, 1.3449804      , 1.0},
		{0, 2, 0.6476825      , 1.0},
		{0, 2, 0.3198804      , 1.0},
		{0, 2, 0.1604487      , 1.0},
		{0, 2, 0.0809001      , 1.0},
		{0, 3, 2.9655902      , 0.9},
		{0, 3, 0.7221518      , 0.4},
		{0, 3, 0.3685716      , 1.0},
		{0, 3, 0.095697       , 1.0},
		{0, 4, 0.57653060000  , 1.0}},		//Al
		{
		{0, 0, 2858.338883    , 0.1},
		{0, 0, 73.7956599     , 1.0},
		{0, 0, 33.8600909     , 1.0},
		{0, 0, 16.3287613     , 1.0},
		{0, 0, 8.2510141      , 1.0},
		{0, 0, 4.3529745      , 1.0},
		{0, 0, 2.3876618      , 1.0},
		{0, 0, 1.3551863      , 1.0},
		{0, 0, 0.7916905      , 1.0},
		{0, 0, 0.4732753      , 1.0},
		{0, 0, 0.2877135      , 1.0},
		{0, 0, 0.176701       , 1.0},
		{0, 0, 0.1088943      , 1.0},
		{0, 1, 472.1002712    , 0.1},
		{0, 1, 70.9586352     , 0.9},
		{0, 1, 29.506044      , 1.0},
		{0, 1, 12.7970916     , 1.0},
		{0, 1, 5.7652107      , 1.0},
		{0, 1, 2.685288       , 1.0},
		{0, 1, 1.2864019      , 1.0},
		{0, 1, 0.6302347      , 1.0},
		{0, 1, 0.3138494      , 1.0},
		{0, 1, 0.1578529      , 1.0},
		{0, 1, 0.0796582      , 1.0},
		{0, 2, 110.8657122    , 0.2},
		{0, 2, 16.8153379     , 1.0},
		{0, 2, 7.1390557      , 1.0},
		{0, 2, 3.1858809      , 1.0},
		{0, 2, 1.4837982      , 1.0},
		{0, 2, 0.7153555      , 1.0},
		{0, 2, 0.3537618      , 1.0},
		{0, 2, 0.1776898      , 1.0},
		{0, 2, 0.0897205      , 1.0},
		{0, 3, 4.3695835      , 0.9},
		{0, 3, 1.1183819      , 0.4},
		{0, 3, 0.5386777      , 1.0},
		{0, 3, 0.165924       , 1.0},
		{0, 4, 0.68448660000  , 1.0}},		//Si
		{
		{0, 0, 3184.8068027   , 0.1},
		{0, 0, 127.4482163    , 0.6},
		{0, 0, 63.3114097     , 1.0},
		{0, 0, 32.6594173     , 1.0},
		{0, 0, 17.4474402     , 1.0},
		{0, 0, 9.6232184      , 1.0},
		{0, 0, 5.4612796      , 1.0},
		{0, 0, 3.1771054      , 1.0},
		{0, 0, 1.8870259      , 1.0},
		{0, 0, 1.1393644      , 1.0},
		{0, 0, 0.6961763      , 1.0},
		{0, 0, 0.4284584      , 1.0},
		{0, 0, 0.2643306      , 1.0},
		{0, 1, 516.9904005    , 0.3},
		{0, 1, 76.117926      , 1.0},
		{0, 1, 31.4077596     , 1.0},
		{0, 1, 13.539519      , 1.0},
		{0, 1, 6.0722814      , 1.0},
		{0, 1, 2.8196165      , 1.0},
		{0, 1, 1.3482686      , 1.0},
		{0, 1, 0.6599992      , 1.0},
		{0, 1, 0.3286482      , 1.0},
		{0, 1, 0.1653611      , 1.0},
		{0, 1, 0.0834931      , 1.0},
		{0, 2, 121.1307862    , 0.2},
		{0, 2, 18.0623226     , 1.0},
		{0, 2, 7.6218265      , 1.0},
		{0, 2, 3.3857032      , 1.0},
		{0, 2, 1.571738       , 1.0},
		{0, 2, 0.7561454      , 1.0},
		{0, 2, 0.3734617      , 1.0},
		{0, 2, 0.187449       , 1.0},
		{0, 2, 0.094597       , 1.0},
		{0, 3, 6.0884261      , 0.9},
		{0, 3, 1.7078374      , 0.5},
		{0, 3, 0.7297299      , 1.0},
		{0, 3, 0.2582646      , 1.0},
		{0, 4, 0.82094560000  , 1.0}},		//P
		{
		{0, 0, 3681.2546546   , 0.1},
		{0, 0, 139.8936849    , 0.6},
		{0, 0, 68.8574614     , 1.0},
		{0, 0, 35.2463716     , 1.0},
		{0, 0, 18.7104368     , 1.0},
		{0, 0, 10.2682007     , 1.0},
		{0, 0, 5.8052483      , 1.0},
		{0, 0, 3.3680961      , 1.0},
		{0, 0, 1.9969333      , 1.0},
		{0, 0, 1.2045075      , 1.0},
		{0, 0, 0.7356456      , 1.0},
		{0, 0, 0.4526978      , 1.0},
		{0, 0, 0.2792847      , 1.0},
		{0, 1, 559.0682164    , 0.3},
		{0, 1, 80.8537686     , 1.0},
		{0, 1, 33.136945      , 1.0},
		{0, 1, 14.2085391     , 1.0},
		{0, 1, 6.3466236      , 1.0},
		{0, 1, 2.9386601      , 1.0},
		{0, 1, 1.4026898      , 1.0},
		{0, 1, 0.6860093      , 1.0},
		{0, 1, 0.3415063      , 1.0},
		{0, 1, 0.1718518      , 1.0},
		{0, 1, 0.0867926      , 1.0},
		{0, 2, 141.1811358    , 0.2},
		{0, 2, 19.3114271     , 1.0},
		{0, 2, 8.5935318      , 1.0},
		{0, 2, 3.6964401      , 1.0},
		{0, 2, 1.7159205      , 1.0},
		{0, 2, 0.8256256      , 1.0},
		{0, 2, 0.4078913      , 1.0},
		{0, 2, 0.2048044      , 1.0},
		{0, 2, 0.1033961      , 1.0},
		{0, 3, 6.440697       , 1.0},
		{0, 3, 1.4844098      , 1.0},
		{0, 3, 0.7619737      , 1.0},
		{0, 3, 0.3060145      , 1.0},
		{0, 4, 0.85880130000  , 1.0}},		//S
		{
		{0, 0, 4314.492304    , 0.1},
		{0, 0, 154.0662704    , 0.6},
		{0, 0, 74.9669272     , 1.0},
		{0, 0, 37.9959262     , 1.0},
		{0, 0, 20.0023883     , 1.0},
		{0, 0, 10.9019372     , 1.0},
		{0, 0, 6.1295625      , 1.0},
		{0, 0, 3.5409446      , 1.0},
		{0, 0, 2.0925534      , 1.0},
		{0, 0, 1.2591195      , 1.0},
		{0, 0, 0.7676092      , 1.0},
		{0, 0, 0.4716926      , 1.0},
		{0, 0, 0.2906248      , 1.0},
		{0, 1, 608.0316901    , 0.2},
		{0, 1, 85.7504275     , 0.9},
		{0, 1, 34.8051637     , 1.0},
		{0, 1, 14.8077119     , 1.0},
		{0, 1, 6.5744146      , 1.0},
		{0, 1, 3.0306998      , 1.0},
		{0, 1, 1.4422826      , 1.0},
		{0, 1, 0.7040776      , 1.0},
		{0, 1, 0.3501619      , 1.0},
		{0, 1, 0.1761315      , 1.0},
		{0, 1, 0.088932       , 1.0},
		{0, 2, 195.8762344    , 0.2},
		{0, 2, 27.1137012     , 1.0},
		{0, 2, 12.1220909     , 1.0},
		{0, 2, 5.5641243      , 1.0},
		{0, 2, 2.5442859      , 1.0},
		{0, 2, 1.1902538      , 1.0},
		{0, 2, 0.5661465      , 1.0},
		{0, 2, 0.2720231      , 1.0},
		{0, 2, 0.131146       , 1.0},
		{0, 3, 8.531085       , 1.0},
		{0, 3, 2.4197055      , 1.0},
		{0, 3, 0.8922925      , 1.0},
		{0, 3, 0.3552024      , 1.0},
		{0, 4, 1.03198430000  , 1.0}},		//Cl
		{
		{0, 0, 5263.6806109   , 0.1},
		{0, 0, 187.96084989   , 0.6},
		{0, 0, 91.459651184   , 1.0},
		{0, 0, 46.355029964   , 1.0},
		{0, 0, 24.402913726   , 1.0},
		{0, 0, 13.300363384   , 1.0},
		{0, 0, 7.47806625     , 1.0},
		{0, 0, 4.319952412    , 1.0},
		{0, 0, 2.552915148    , 1.0},
		{0, 0, 1.53612579     , 1.0},
		{0, 0, 0.936483224    , 1.0},
		{0, 0, 0.575464972    , 1.0},
		{0, 0, 0.354562256    , 1.0},
		{0, 1, 741.79866192   , 0.2},
		{0, 1, 104.61552155   , 0.9},
		{0, 1, 42.462299714   , 1.0},
		{0, 1, 18.065408518   , 1.0},
		{0, 1, 8.020785812    , 1.0},
		{0, 1, 3.697453756    , 1.0},
		{0, 1, 1.759584772    , 1.0},
		{0, 1, 0.858974672    , 1.0},
		{0, 1, 0.427197518    , 1.0},
		{0, 1, 0.21488043     , 1.0},
		{0, 1, 0.10849704     , 1.0},
		{0, 2, 238.96900597   , 0.2},
		{0, 2, 33.078715464   , 1.0},
		{0, 2, 14.788950898   , 1.0},
		{0, 2, 6.788231646    , 1.0},
		{0, 2, 3.104028798    , 1.0},
		{0, 2, 1.452109636    , 1.0},
		{0, 2, 0.69069873     , 1.0},
		{0, 2, 0.331868182    , 1.0},
		{0, 2, 0.15999812     , 1.0},
		{0, 3, 10.4079237     , 1.0},
		{0, 3, 2.95204071     , 1.0},
		{0, 3, 1.08859685     , 1.0},
		{0, 3, 0.433346928    , 1.0},
		{0, 4, 1.25902084600  , 1.0}},		//Ar
		{
		{0, 0, 8556.7117611   , 0.1},
		{0, 0, 146.5057165    , 1.0},
		{0, 0, 60.8415375     , 1.0},
		{0, 0, 26.6304714     , 1.0},
		{0, 0, 12.2605057     , 1.0},
		{0, 0, 5.922457       , 1.0},
		{0, 0, 2.992628       , 1.0},
		{0, 0, 1.5762429      , 1.0},
		{0, 0, 0.8618648      , 1.0},
		{0, 0, 0.4869592      , 1.0},
		{0, 0, 0.2828475      , 1.0},
		{0, 0, 0.1679521      , 1.0},
		{0, 0, 0.1013417      , 1.0},
		{0, 0, 0.0617485      , 1.0},
		{0, 0, 0.0377474      , 1.0},
		{0, 1, 1052.6413867   , 0.1},
		{0, 1, 36.1442743     , 1.0},
		{0, 1, 16.5969607     , 1.0},
		{0, 1, 7.7864208      , 1.0},
		{0, 1, 3.7229718      , 1.0},
		{0, 1, 1.8094243      , 1.0},
		{0, 1, 0.891429       , 1.0},
		{0, 1, 0.4438892      , 1.0},
		{0, 1, 0.2227452      , 1.0},
		{0, 1, 0.1122955      , 1.0},
		{0, 1, 0.0567011      , 1.0},
		{0, 2, 245.9248263    , 0.1},
		{0, 2, 19.2095772     , 1.0},
		{0, 2, 8.7467662      , 1.0},
		{0, 2, 4.089444       , 1.0},
		{0, 2, 1.9563381      , 1.0},
		{0, 2, 0.9540003      , 1.0},
		{0, 2, 0.4723243      , 1.0},
		{0, 2, 0.2364281      , 1.0},
		{0, 2, 0.1191361      , 1.0},
		{0, 2, 0.0601664      , 1.0},
		{0, 3, 7.2556601      , 1.0},
		{0, 3, 1.8118642      , 1.0},
		{0, 3, 0.5621445      , 1.0},
		{0, 3, 0.2476053      , 1.0},
		{0, 3, 0.069923       , 1.0},
		{0, 4, 5.3806511      , 1.0},
		{0, 4, 1.7935531      , 1.0},
		{0, 4, 0.59785110000  , 1.0}},		//K
		{
		{0, 0, 8607.040797    , 0.1},
		{0, 0, 186.6761164    , 1.0},
		{0, 0, 85.4284093     , 1.0},
		{0, 0, 40.9597011     , 1.0},
		{0, 0, 20.2067221     , 1.0},
		{0, 0, 9.4283304      , 1.0},
		{0, 0, 5.5105726      , 1.0},
		{0, 0, 3.0283137      , 1.0},
		{0, 0, 1.7139817      , 1.0},
		{0, 0, 0.9952817      , 1.0},
		{0, 0, 0.5904908      , 1.0},
		{0, 0, 0.3563493      , 1.0},
		{0, 0, 0.2177206      , 1.0},
		{0, 0, 0.1353608      , 1.0},
		{0, 0, 0.0843695      , 1.0},
		{0, 1, 1138.8086874   , 0.1},
		{0, 1, 40.25882       , 1.0},
		{0, 1, 18.5692989     , 1.0},
		{0, 1, 8.7423026      , 1.0},
		{0, 1, 4.190985       , 1.0},
		{0, 1, 2.0406572      , 1.0},
		{0, 1, 1.0065546      , 1.0},
		{0, 1, 0.5015582      , 1.0},
		{0, 1, 0.2517606      , 1.0},
		{0, 1, 0.1269331      , 1.0},
		{0, 1, 0.0640921      , 1.0},
		{0, 2, 272.6361647    , 0.1},
		{0, 2, 21.8506384     , 1.0},
		{0, 2, 10.0043167     , 1.0},
		{0, 2, 4.6969297      , 1.0},
		{0, 2, 2.2536441      , 1.0},
		{0, 2, 1.101133       , 1.0},
		{0, 2, 0.5457906      , 1.0},
		{0, 2, 0.27335        , 1.0},
		{0, 2, 0.137765       , 1.0},
		{0, 2, 0.0695779      , 1.0},
		{0, 3, 20.7187048     , 0.8},
		{0, 3, 6.5897516      , 0.6},
		{0, 3, 2.1054414      , 1.0},
		{0, 3, 0.7167493      , 1.0},
		{0, 3, 0.2632969      , 1.0},
		{0, 4, 7.174917       , 1.0},
		{0, 4, 2.391639       , 1.0},
		{0, 4, 0.79721300000  , 1.0}},		//Ca
	});

inline double hypergeometric(double a, double b, double c, double x)
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

inline cdouble hypergeometric(double a, double b, double c, cdouble x)
{
	const double TOLERANCE = 1.0e-10;
	cdouble term = a * b * x / c;
	cdouble value = 1.0 + term;
	int n = 1;

	while (std::abs(term) > TOLERANCE)
	{
		a++, b++, c++, n++;
		term *= a * b * x / c / (double)n;
		value += term;
	}

	return value;
}

inline bool ends_with(const std::string &str, const std::string &suffix)
{
	if (str.length() >= suffix.length())
	{
		return str.compare(str.length() - suffix.length(), suffix.length(), suffix) == 0;
	}
	return false;
}

template <typename T>
inline bool is_nan(T in)
{
	if (typeid(in) == typeid(double))
	{
		return in != in;
	}
	else if (typeid(in) == typeid(float))
	{
		return in != in;
	}
	else if (typeid(in) == typeid(long double))
	{
		return in != in;
	}
	else
	{
		return false;
	}
}

#include "wfn_class.h"
#include "atoms.h"