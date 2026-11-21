#pragma once
#include "../pch.h"
#include <bit>
#include "throughput.h"
#include "extinction.h"
#include "tuning.h"
#include "nos_types.h"
#include "error_check.h"
#include "properties_options.h"
#include "hkl.h"

int vec_sum(const bvec& in);
int vec_sum(const ivec& in);
double vec_sum(const vec& in);
cdouble vec_sum(const cvec& in);
double vec_length(const vec& in);
template <typename array>
const double array_length(const array& in)
{
	return std::hypot(in[0], in[1], in[2]);
}
template <typename array>
const double array_length(const array& in, const array& in2)
{
	if (std::size(in) == 3 && std::size(in2) == 3)
	{
		return std::hypot(in[0] - in2[0], in[1] - in2[1], in[2] - in2[2]);
	}

	double sum = 0.0;
	auto it1 = std::begin(in);
	auto it2 = std::begin(in2);
	for (; it1 != std::end(in); it1++, it2++) {
		sum += (*it1 - *it2) * (*it1 - *it2);
	}
	return sqrt(sum);
}

d3 vec_diff(const d3& a, const d3& b);

d3 vec_cross(const d3& a, const d3& b);

double vec_dot(const d3& a, const d3& b);

constexpr const std::complex<double> c_one(0, 1.0);

extern std::string help_message;
std::string NoSpherA2_message(bool no_date = false);
std::string WAEmpiRe_message();
extern std::string build_date;

#include "simd_kernels.h"

#include "sha256.h"

bool is_similar_rel(const double& first, const double& second, const double& tolerance);
bool is_similar(const double& first, const double& second, const double& tolerance);
bool is_similar_abs(const double& first, const double& second, const double& tolerance);
std::filesystem::path get_home_path(void);
bool ensure_occ_data_path(const char* argv0);
char asciitolower(char in);

bool generate_cart2sph_mat(vec2& d, vec2& f, vec2& g, vec2& h);
std::string go_get_string(std::ifstream& file, std::string search, bool rewind = true);

const int sht2nbas(const int& type);

const int shell2function(const int& type, const int& prim);

template <class T>
std::string toString(const T& t)
{
	std::ostringstream stream;
	stream << t;
	return stream.str();
}

template <class T>
T fromString(const std::string& s)
{
	std::istringstream stream(s);
	T t;
	stream >> t;
	return t;
}

template <typename T>
void shrink_vector(std::vector<T>& g)
{
	g.clear();
	std::vector<T>(g).swap(g);
}

template <class T>
std::vector<T> split_string(const std::string& input, const std::string delimiter)
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

void remove_empty_elements(svec& input, const std::string& empty = " ");
std::chrono::high_resolution_clock::time_point get_time();
long long int get_musec(std::chrono::high_resolution_clock::time_point start, std::chrono::high_resolution_clock::time_point end);
long long int get_msec(std::chrono::high_resolution_clock::time_point start, std::chrono::high_resolution_clock::time_point end);
long long int get_sec(std::chrono::high_resolution_clock::time_point start, std::chrono::high_resolution_clock::time_point end);

void write_timing_to_file(std::ostream& file, std::vector<_time_point> time_points, std::vector<std::string> descriptions);

int CountWords(const char* str);

void copy_file(std::filesystem::path& from, std::filesystem::path& to);
std::string shrink_string(std::string& input);
std::string shrink_string_to_atom(std::string& input, const int& atom_number);
//------------------Functions to work with configuration files--------------------------
bool check_bohr(const WFN& wave, bool debug);

void select_cubes(std::vector<std::vector<unsigned int>>& selection, std::vector<WFN>& wavy, unsigned int nr_of_cubes = 1, bool wfnonly = false, bool debug = false);
bool unsaved_files(std::vector<WFN>& wavy);

std::string trim(const std::string& s);

/**
 * @brief Physical memory this process can get in bytes, capped by a scheduler or container limit
 * (job object, cgroup); 0 when unknown, and the caller keeps its default.
 */
size_t available_memory_bytes();
/** @brief Resident memory this process holds now in bytes; 0 when unknown. */
size_t resident_memory_bytes();

#include "string_io.h"

//How many equal-sized items to keep resident inside a memory budget; 0 means "hold all of them", as does a budget of 0
//The fits-test is a division rather than n_items * item_bytes because that product is exactly what overflows on the structures this serves
inline size_t items_within_budget(const size_t n_items, const size_t item_bytes, const size_t budget_bytes)
{
	if (budget_bytes == 0 || n_items == 0 || item_bytes == 0)
		return 0;
	const size_t n = budget_bytes / item_bytes;
	if (n_items <= n)
		return 0;
	//a single item larger than the whole budget still has to be processed, one at a time
	return n ? n : 1;
}

#include "progress_bar.h"

//even_steps rounds the point count up to even so the box centre lies on a grid plane. Pass false only when
//the step is the resolution, not (Max-Min)/NbSteps: there an extra point only enlarges the box.
void readxyzMinMax_fromWFN(
	const WFN& wavy,
	properties_options& opts,
	const bool even_steps = true);

void readxyzMinMax_fromCIF(
	std::filesystem::path cif,
	properties_options& opts,
	vec2& cm);

double read_CIF(std::istream& cif_input, const cell& unit_cell, int& ncen,
	std::vector<asym_atom>& asym_atoms, std::vector<vec2>& ADPs, const bool debug);

double double_from_string_with_esd(std::string in);

void swap_sort(ivec order, cvec& v);

void swap_sort_multi(ivec order, std::vector<ivec>& v);

// Given a 3x3 symmetric matrix in a single row-major array of double, returns the median (middle) eigenvalue
double get_lambda_1(double* a);

double get_decimal_precision_from_CIF_number(std::string& given_string);

double bessel_first_kind(int l, double x);

constexpr unsigned int doublefactorial(int n)
{
	if (n <= 1)
		return 1;
	return n * doublefactorial(n - 2);
}

template <typename T>
void removeElement(std::vector<T>& vec, const T& x)
{
	auto new_end = std::remove(vec.begin(), vec.end(), x);
	vec.erase(new_end, vec.end());
}

inline void Enter() {
	std::cout << "press ENTER to continue... " << std::flush;
	std::cin.ignore();
	std::cin.get();
};

inline void cls() {
#ifdef _WIN32
	// On modern Windows 10+ terminals, ANSI codes are often supported.
	std::cout << "\033[2J\033[H";
#else
	std::cout << "\033[2J\033[H";
#endif
	std::cout.flush();
}

inline bool yesno() {
	bool end = false;
	while (!end) {
		char dum = '?';
		std::cout << "(Y/N)?";
		std::cin >> dum;
		if (dum == 'y' || dum == 'Y') {
			std::cout << "Okay..." << std::endl;
			return true;
		}
		else if (dum == 'N' || dum == 'n') return false;
		else std::cout << "Sorry, i did not understand that!" << std::endl;
	}
	return false;
};

#include "primitive.h"

#include "options.h"

double hypergeometric(double a, double b, double c, double x);

cdouble hypergeometric(double a, double b, double c, cdouble x);

bool ends_with(const std::string& str, const std::string& suffix);

bool read_block_from_fortran_binary(std::ifstream& file, void* Target, const size_t capacity);
template <typename T>
bool read_block_from_fortran_binary(std::ifstream& file, std::vector<T>& Target);
