#pragma once
#define WIN32_LEAN_AND_MEAN
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <complex>
#include <iomanip>
#include <omp.h>
#include <regex>
#include <algorithm>
#include <set>
#include <numeric>
#include <complex>
#include <functional>
#include <stdexcept>
#include <typeinfo>
#ifdef _WIN32
#include <direct.h>
#define GetCurrentDir _getcwd(NULL, 0)
#include <io.h>
#else
#define GetCurrentDir getcwd
#include <unistd.h>
#include <cfloat>
#include <sys/wait.h>
#include <sys/time.h>
#include <termios.h>
#endif
#ifdef __cplusplus__
#include <cstdlib>
#else
#include <stdlib.h>
#endif

class WFN;
class cell;
class atom;

typedef std::complex<double> cdouble;
typedef std::vector<double> vec;
typedef std::vector<int> ivec;
typedef std::vector<cdouble> cvec;

inline double vec_sum(vec& in) {
  double res = 0.0;
  for (int i = 0; i < in.size(); i++)
    res += in[i];
  return res;
}

inline cdouble vec_sum(cvec& in) {
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
    //Constants for later use
    constexpr double SQRT2 = sqrt(2.0);
    constexpr int hardness = 3;
    constexpr double cutoff = 1.0e-20;
    constexpr double PI = 3.14159265358979323846;
    constexpr double TWO_PI = 2 * PI;
    constexpr double FOUR_PI = 4 * PI;
    constexpr double C0 = SQRT2 * FOUR_PI;
    constexpr double PI2 = PI * PI;
    constexpr double PI3 = PI * PI * PI;
    constexpr double PI_180 = PI / 180.0;
    const double TG32 = tgamma(3.0 / 2.0);
    constexpr double ED_fact = 0.023934;
    constexpr int max_LT = 33;
    constexpr int MAG = 5810;
    //                       3,     5     7,    9,    11,   13,   15,   17
    //                      19,    21
    constexpr int lebedev_table[33] = { 6,    14,   26,   38,   50,   74,   86,   110,
                 146,  170,  194,  230,  266,  302,  350,  434,
                 590,  770,  974,  1202, 1454, 1730, 2030, 2354,
                 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810 };
    constexpr long long int ft[21]{ 1,1,2,6,24,120,720,5040,40320,362880,3628800,39916800,479001600,6227020800,87178291200,1307674368000,20922789888000,355687428096000,6402373705728000,121645100408832000,2432902008176640000 };
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
    constexpr double a0 = 0.529177210903E-10; //in m
    constexpr double h = 6.62607015E-34 / 1.602176634E-19; //in eV*s
    constexpr double Ryd_ener = 13.6056923; //in eV
    constexpr double alpha = 0.0072973525693; //Sommerfeld fine structure constant
    constexpr double el_mass = 9.1093837015E-31; //in kg
    constexpr double el_charge = 1.602176634E-19; // in C
    constexpr double speed_of_light = 2.99792458E8; //m/s

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

    inline const long long int ft_fun(const int& nr) {
        if (nr >= 0 && nr <= 20)
            return ft[nr];
        else if (nr < 0)
            return 0;
        else
            return ft_fun(nr - 1) * nr;
    }

    constexpr double bohr2ang(const double& inp)
    {
        return inp * 0.529177249;
    }

    constexpr double bohr2ang_p(const double& inp, const int& p)
    {
        if (p == 0)
            return 1.0;
        else if (p == 1) 
            return bohr2ang(inp);
        else
            return bohr2ang_p(bohr2ang(inp), p - 1);
    }

    constexpr double ang2bohr(const double& inp)
    {
        return inp / 0.529177249;
    }
    inline const double ang2bohr_p(const double& inp, const int& p)
    {
        if (p == 0)
            return 1.0;
        else if (p == 1)
            return ang2bohr(inp);
        else
            return ang2bohr_p(ang2bohr(inp), p - 1);
    }


    constexpr double cubic_ang2bohr(const double& inp)
    {
        return inp / (0.529177249 * 0.529177249 * 0.529177249);
    }

    constexpr double cubic_bohr2ang(const double& inp)
    {
        return inp * (0.529177249 * 0.529177249 * 0.529177249);
    }

    //------------------general functions for easy use of terminal input--------------------
    constexpr double bragg_angstrom[114]{
      0.00, //DUMMY LINE
      0.35,																																														                                                                                          0.35,
      1.45, 1.05,																																																										                              0.85, 0.70, 0.65, 0.60, 0.50, 0.45,
      1.80, 1.50,																																																										                              1.25, 1.10, 1.00, 1.00, 1.00, 1.00,
      2.20, 1.80,																																											1.60, 1.40, 1.35, 1.40, 1.40, 1.40, 1.35, 1.35, 1.35, 1.35, 1.30, 1.25, 1.15, 1.15, 1.15, 1.10,
      2.35, 2.00,																																											1.80, 1.55, 1.45, 1.45, 1.35, 1.30, 1.35, 1.40, 1.60, 1.55, 1.55, 1.45, 1.45, 1.40, 1.40, 1.40,
      2.60, 2.15, 1.95, 1.85, 1.85, 1.85, 1.85, 1.85, 1.85, 1.80, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.55, 1.45, 1.35, 1.30, 1.30, 1.35, 1.35, 1.35, 1.50, 1.90, 1.75, 1.60, 1.90, 1.50, 1.50,
      2.80, 2.35, 2.15, 2.05, 2.05, 2.05, 2.05, 2.05, 2.05, 2.00, 1.95, 1.95, 1.95, 1.95, 1.95, 1.95 };

    //Covalent Radii according to the CSD
    constexpr double covalent_radii[114]{
       0.0,
      0.23,	                                                                                                                                                                                    1.5,
      1.28, 0.96,	                                                                                                                                                0.83,	0.68,	0.68,	0.68,	0.64, 1.5,
      1.66,	1.41,																																																																									1.21,	1.2,	1.05,	1.02,	0.99,	1.51,
      2.03,	1.76,																																											1.7,	1.6,	1.53,	1.39,	1.61,	1.52,	1.26,	1.24,	1.32,	1.22,	1.22,	1.17,	1.21,	1.22,	1.21,	1.5,
      2.2,	1.95,																																											1.9,	1.75,	1.64,	1.54,	1.47,	1.46,	1.42,	1.39,	1.45,	1.54,	1.42,	1.39,	1.39,	1.47,	1.4,	1.5,
      2.44,	2.15,	2.07,	2.04,	2.03,	2.01,	1.99,	1.98,	1.98,	1.96,	1.94,	1.92,	1.92,	1.89,	1.9,	1.87,	1.87,	1.75,	1.7,	1.62,	1.51,	1.44,	1.41,	1.36,	1.36,	1.32,	1.45,	1.46,	1.48,	1.4,	1.21,	1.5,
      2.6,	2.21,	2.15,	2.06,	2.00,	1.96,	1.9,	1.87,	1.8,	1.69,	1.54,	1.83,	1.5,	1.5,	1.5,	1.5,	1.5,	1.5,	1.5,	1.5,	1.5,	1.5,	1.5,	1.5
    };

    //Integer atom masses
    constexpr unsigned int integer_masses[]{
     1,																																																	4,
     7,  9,																																												11,12,14,16,19,20,
     23,24,																																												27,28,31,32,35,40,
     39,40,																																		45,48,51,52,55,56,59,58,63,64,			69,74,75,80,79,84,
     85, 87,																																	88, 91, 92, 96, 98, 101, 103, 106, 108, 112, 115, 119, 122, 128, 127, 131,
     132,137,139, 140, 141, 144, 145, 150, 152, 157, 159, 163, 165, 167, 169, 173, 175, 178, 181, 184, 186, 190, 192, 195, 197, 201, 204, 207, 209, 209, 210, 222 };

    constexpr double real_masses[]{
     1.0079,																																																																	4.0026,
     6.941,		9.0122,																																											                            			10.811,	12.011,	14.007,	15.999,	18.998,	20.18,
     22.99,		24.305,        																																											                  				26.986,	28.086,	30.974,	32.065,	35.453,	39.948,
     39.098,	40.078,																																44.956,	47.867,	50.942,	51.996,	54.938,	55.845,	58.933,	58.693,	63.546,	65.38,		69.723,	72.64,	74.922,	78.96,	79.904,	83.798,
     85.468,	87.62,																																88.906, 91.224,	92.906, 95.96,	97.90,	101.07,	102.91,	106.42,	107.87,	112.41,		114.82, 118.71,	121.76,	127.6,	126.9,	131.29,
     132.91,	137.33,		139.91, 140.12, 140.91, 144.24, 144.9, 150.36, 151.96, 157.25, 158.93, 162.5, 164.93, 167.26, 168.93, 173.05, 174.97,			178.49, 180.95, 183.84, 186.21, 190.23, 192.22, 195.08, 196.97, 200.59,		204.38, 207.2,	208.98, 208.9,	209.9,	222.0 };
}
//bool yesno();
bool is_similar_rel(const double& first, const double& second, const double& tolerance);
bool is_similar(const double& first, const double& second, const double& tolerance);
bool is_similar_abs(const double& first, const double& second, const double& tolerance);
void cls();
std::string get_home_path(void);
void join_path(std::string& s1, std::string& s2);
inline char asciitolower(char in) {
  if (in <= 'Z' && in >= 'A')
    return in - ('Z' - 'z');
  return in;
}
inline void error_check(const bool condition, const std::string& file, const int& line, const std::string& function, const std::string& error_mesasge, std::ostream& log_file = std::cout)
{
  if (!condition) {
    log_file << "Error in " << function << " at: " << file << " : " << line << " " << error_mesasge << std::endl;
    log_file.flush();
    exit(-1);
  }
};
inline void not_implemented(const std::string& file, const int& line, const std::string& function, const std::string& error_mesasge, std::ostream& log_file) {
  log_file << function << " at: " << file << ":" << line << " " << error_mesasge << " not yet implemented!" << std::endl;
  log_file.flush();
  exit(-1);
};
#define err_checkf(condition, error_message, file) error_check(condition, __FILE__, __LINE__, __func__, error_message, file)
#define err_chkf(condition, error_message, file) error_check(condition, __FILE__, __LINE__, __func__, error_message, file)
#define err_chekf(condition, error_message, file) error_check(condition, __FILE__, __LINE__, __func__, error_message, file)
#define err_not_impl_f(error_message, file) not_implemented(__FILE__, __LINE__, __func__, error_message, file)

bool generate_sph2cart_mat(std::vector<vec>& p, std::vector<vec>& d, std::vector<vec>& f, std::vector<vec>& g);
bool generate_cart2sph_mat(std::vector<vec>& d, std::vector<vec>& f, std::vector<vec>& g, std::vector<vec>& h);
std::string go_get_string(std::ifstream& file, std::string search, bool rewind = true);

inline const int sht2nbas(const int& type)
{
  const int st2bas[6]{ 1,3,6,10,15,21 };
  const int nst2bas[6]{ 11,9,7,5,4,1 };
  if (type >= 0)
    return st2bas[type];
  else
    return nst2bas[5 + type];
};

inline const int shell2function(const int& type, const int& prim)
{
  switch (type) {
  case (-5):
    return -32 + prim;
  case (-4):
    return -21 + prim;
  case(-3):
    return -12 + prim;
  case(-2):
    return -5 + prim;
  case(-1):
    return 1 + prim;
  case(0):
    return 1;
  case(1):
    return 2 + prim;
  case(2):
    return 5 + prim;
  case(3):
    if (prim == 0) return 11;
    if (prim == 1) return 12;
    if (prim == 2) return 13;
    if (prim == 3) return 17;
    if (prim == 4) return 14;
    if (prim == 5) return 15;
    if (prim == 6) return 18;
    if (prim == 7) return 19;
    if (prim == 8) return 16;
    if (prim == 9) return 20;
    break;
  case(4):
    return 21 + prim;
  case(5):
    return 36 + prim;
  default:
    return 0;
  }
  return 0;
}

const double normgauss(const int& type, const double& exp);
const double spherical_harmonic(const int& l, const int& m, const double* d);

template<class T> std::string toString(const T& t)
{
  std::ostringstream stream;
  stream << t;
  return stream.str();
}

template<class T> T fromString(const std::string& s)
{
  std::istringstream stream(s);
  T t;
  stream >> t;
  return t;
}

template<typename T> void shrink_vector(std::vector<T>& g)
{
  g.clear();
  std::vector<T>(g).swap(g);
}

template <class T> std::vector<T> split_string(const std::string& input, const std::string delimiter)
{
  std::string input_copy = input + delimiter; // Need to add one delimiter in the end to return all elements
  std::vector<T> result;
  size_t pos = 0;
  while ((pos = input_copy.find(delimiter)) != std::string::npos) {
    result.push_back(fromString<T>(input_copy.substr(0, pos)));
    input_copy.erase(0, pos + delimiter.length());
  }
  return result;
};

inline void remove_empty_elements(std::vector <std::string>& input, const std::string& empty = " ")
{
  for (int i = (int) input.size() - 1; i >= 0; i--)
    if (input[i] == empty || input[i] == "")
      input.erase(input.begin() + i);
}

inline int CountWords(const char* str)
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

inline bool exists(const std::string& name)
{
    if (FILE* file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    }
    else {
        return false;
    }
};

std::string atnr2letter(const int& nr);
void copy_file(std::string& from, std::string& to);
std::string shrink_string(std::string& input);
std::string shrink_string_to_atom(std::string& input, const int& atom_number);
std::string get_filename_from_path(const std::string& input);
std::string get_foldername_from_path(const std::string& input);
std::string get_basename_without_ending(const std::string& input);
//------------------Functions to read from .fchk files----------------------------------
bool read_fchk_integer_block(std::ifstream& in, std::string heading, ivec& result, bool rewind = true);
bool read_fchk_double_block(std::ifstream& in, std::string heading, vec& result, bool rewind = true);
int read_fchk_integer(std::string in);
int read_fchk_integer(std::ifstream& in, std::string search, bool rewind = true);
double read_fchk_double(std::string in);
double read_fchk_double(std::ifstream& in, std::string search, bool rewind = true);
//------------------Functions to work with configuration files--------------------------
void write_template_confi();
int program_confi(std::string& gaussian_path, std::string& turbomole_path,
  std::string& basis, int& ncpus, double& mem, bool debug = false, bool expert = false, unsigned int counter = 0);
bool check_bohr(WFN& wave, bool debug);
int filetype_identifier(std::string& file, bool debug = false);

/*bool open_file_dialog(std::string &path, bool debug, std::vector <std::string> filter);
bool save_file_dialog(std::string &path, bool debug, const std::vector<std::string> &endings, const std::string &filename_given);
bool save_file_dialog(std::string &path, bool debug, const std::vector<std::string> &endings);*/
void select_cubes(std::vector <std::vector <unsigned int> >& selection, std::vector<WFN>& wavy, unsigned int nr_of_cubes = 1, bool wfnonly = false, bool debug = false);
bool unsaved_files(std::vector<WFN>& wavy);
int get_Z_from_label(const char* tmp);

inline int sum_of_bools(const std::vector<bool> in)
{
  int result = 0;
  for (int i = 0; i < in.size(); i++)
    if (in[i]) result++;
  return result;
}

inline std::string trim(const std::string& s)
{
  if (s == "") return "";
  auto start = s.begin();
  while (start != s.end() && std::isspace(*start)) {
    start++;
  }

  auto end = s.end();
  do {
    end--;
  } while (std::distance(start, end) > 0 && std::isspace(*end));

  return std::string(start, end + 1);
}

//-------------------------Progress_bar--------------------------------------------------

class progress_bar
{
  static const auto overhead = sizeof " [100%]";

  std::ostream& os;
  const std::size_t bar_width;
  std::string message;
  const std::string full_bar;
  const double precision;

public:
  progress_bar(std::ostream& os, std::size_t line_width,
    std::string message_, const char symbol = '=', const double p = 0.05)
    : os{ os },
    bar_width{ line_width - overhead },
    message{ std::move(message_) },
    full_bar{ std::string(bar_width, symbol) + std::string(bar_width, ' ') },
    precision{ p }
  {
    if (message.size() + 1 >= bar_width || message.find('\n') != message.npos) {
      os << message << '\n';
      message.clear();
    }
    else {
      message += ' ';
    }
    write(0.0);
  }

  // not copyable
  progress_bar(const progress_bar&) = delete;
  progress_bar& operator=(const progress_bar&) = delete;

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
  WFN& wavy,
  double* CoordMinMax,
  int* NbSteps,
  double Radius,
  double Increments,
  bool no_bohr = false);

void readxyzMinMax_fromCIF(
  std::string cif,
  double* CoordMinMax,
  int* NbSteps,
  std::vector < std::vector < double > >& cm,
  double Resolution,
  std::ofstream& file,
  bool debug = false);

void type2vector(
  const int& index,
  int* vector);

bool read_fracs_ADPs_from_CIF(std::string cif, WFN& wavy, cell& unit_cell, std::ofstream& log3, bool debug);

inline double double_from_string_with_esd(std::string in)
{
  if (in.find('(') == std::string::npos)
    return stod(in);
  else
    return stod(in.substr(0, in.find('(')));
}

void swap_sort(ivec order, cvec& v);

void swap_sort_multi(ivec order, std::vector<ivec>& v);

//Given a 3x3 matrix in a single array of double will find and sort eigenvalues and return biggest eigenvalue
double get_lambda_1(double* a);

double get_decimal_precision_from_CIF_number(std::string& given_string);

template <typename numtype = int>
struct hashFunction {
  size_t operator()(const std::vector<numtype>& myVector) const {
    std::hash<numtype> hasher;
    size_t answer = 0;
    for (numtype i : myVector) {
      answer ^= hasher(i) + 0x9e3779b9 + (answer << 6) + (answer >> 2);
    }
    return answer;
  }
};

template <typename numtype = int>
struct hkl_equal
{
  bool operator()(const std::vector<numtype>& vec1, const std::vector<numtype>& vec2) const {
    const int size = vec1.size();
    if (size != vec2.size()) return false;
    int similar = 0;
    for (int i = 0; i < size; i++) {
      if (vec1[i] == vec2[i])       similar++;
      else if (vec1[i] == -vec2[i]) similar--;
    }
    if (abs(similar) == size) return true;
    else                      return false;
  }
};

template <typename numtype = int>
struct hkl_less
{
  bool operator()(const std::vector<numtype>& vec1, const std::vector<numtype>& vec2) const {
    if (vec1[0] < vec2[0]) {
      return true;
    }
    else if (vec1[0] == vec2[0]) {
      if (vec1[1] < vec2[1]) {
        return true;
      }
      else if (vec1[1] == vec2[1]) {
        if (vec1[2] < vec2[2]) {
          return true;
        }
        else return false;
      }
      else return false;
    }
    else return false;
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
  void normalize() {
    coefficient *= normalization_constant();
  };
  void unnormalize() {
    coefficient /= normalization_constant();
  };
  double normalization_constant() {
    //assuming type is equal to angular momentum
    return norm_const;
  }
  primitive() : center(0), type(0), exp(0.0), coefficient(0.0) {}
  primitive(int c, int t, double e, double coef) : center(c), type(t), exp(e), coefficient(coef) {
    norm_const = pow(
      pow(2, 7 + 4 * type) * pow(exp, 3 + 2 * type)
      / constants::PI / pow(doublefactorial(2 * type + 1), 2),
      0.25
    );
  }
};

typedef std::set<ivec> hkl_list;
typedef std::set<ivec>::const_iterator hkl_list_it;

typedef std::set<vec> hkl_list_d;
typedef std::set<vec>::const_iterator hkl_list_it_d;

//---------------- Object for handling all input options -------------------------------
struct options {
  options() {
    groups.resize(1);
  };
  options(int& argc, char** argv) {
    groups.resize(1);
    look_for_debug(argc, argv);
  };
  void look_for_debug(int& argc, char** argv);
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
  bool density_test_cube = false;
  bool no_date = false;
  bool gbw2wfn = false;
  bool old_tsc = false;
  bool thakkar_d_plot = false;
  bool spherical_harmonic = false;
  bool ML_test = false;
  double sfac_scan = 0.0;
  double sfac_diffuse = 0.0;
  double dmin = 99.0;
  int hirsh_number = 0;
  double MinMax[6]{ 0,0,0,0,0,0 };
  int NbSteps[3]{ 0,0,0 };
  ivec MOs;
  std::vector < ivec > groups;
  std::vector < vec > twin_law;
  std::vector < ivec > combined_tsc_groups;
  bool all_mos = false;
  bool test = false;
  std::string wfn;
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
  bool debug = false;
  hkl_list m_hkl_list;
};


const double gaussian_radial(primitive& p, double& r);

const double calc_density_ML(double& x,
  double& y,
  double& z,
  vec& coefficients,
  std::vector<atom>& atoms,
  const int& exp_coefs);
const double calc_density_ML(double& x,
  double& y,
  double& z,
  vec& coefficients,
  std::vector<atom>& atoms,
  const int& exp_coefs,
  const int& atom_nr);

int load_basis_into_WFN(WFN& wavy, const std::vector<std::vector<primitive>>& b);

inline const std::vector<std::vector<primitive>> TZVP_JKfit(
  {
    { //center  type  exp         coef
      {0,       0, 9.5302493327,  1.0},
      {0,       0, 1.9174506246,  1.0},
      {0,       0, 0.68424049142, 1.0},
      {0,       0, 0.28413255710, 1.0},
      {0,       1, 2.9133232035,  1.0},
      {0,       1, 1.2621205398,  1.0},
      {0,       1, 0.50199775874, 1.0},
      {0,       2, 2.3135329149,  1.0},
      {0,       2, 0.71290724024, 1.0},
      {0,       3, 1.6565726132 , 1.0},
    },//H
    {},//He
    {},//Li
    {},//Be
    {},//B
    {
      {0,        0, 1113.9867719, 1.0},
      {0,        0, 369.16234180, 1.0},
      {0,        0, 121.79275232, 1.0},
      {0,        0, 48.127114540, 1.0},
      {0,        0, 20.365074004, 1.0},
      {0,        0, 8.0883596856, 1.0},
      {0,        0, 2.5068656570, 1.0},
      {0,        0, 1.2438537380, 1.0},
      {0,        0,0.48449899601, 1.0},
      {0,        0,0.19185160296, 1.0},
      {0,        1, 102.99176249, 1.0},
      {0,        1, 28.132594009, 1.0},
      {0,        1, 9.8364318173, 1.0},
      {0,        1, 3.3490544980, 1.0},
      {0,        1, 1.4947618613, 1.0},
      {0,        1,0.57690108899, 1.0},
      {0,        1,0.20320063291, 1.0},
      {0,        2, 10.594068356, 1.0},
      {0,        2, 3.5997195366, 1.0},
      {0,        2, 1.3355691094, 1.0},
      {0,        2,0.51949764954, 1.0},
      {0,        2,0.19954125200, 1.0},
      {0,        3,1.194866338369, 1.0},
      {0,        3, .415866338369, 1.0},
      {0,        4, .858866338369, 1.0}
    },//C
    {
      {0,        0, 1102.8622453, 1.0},
      {0,        0, 370.98041153, 1.0},
      {0,        0, 136.73555938, 1.0},
      {0,        0, 50.755871924, 1.0},
      {0,        0, 20.535656095, 1.0},
      {0,        0, 7.8318737184, 1.0},
      {0,        0, 3.4784063855, 1.0},
      {0,        0, 1.4552856603, 1.0},
      {0,        0,0.63068989071, 1.0},
      {0,        0,0.27276596483, 1.0},
      {0,        1, 93.540954073, 1.0},
      {0,        1, 29.524019527, 1.0},
      {0,        1, 10.917502987, 1.0},
      {0,        1, 4.3449288991, 1.0},
      {0,        1, 1.8216912640, 1.0},
      {0,        1,0.75792424494, 1.0},
      {0,        1,0.28241469033, 1.0},
      {0,        2, 16.419378926, 1.0},
      {0,        2, 5.0104049385, 1.0},
      {0,        2, 1.9793971884, 1.0},
      {0,        2,0.78495771518, 1.0},
      {0,        2,0.28954065963, 1.0},
      {0,        3,1.79354239843, 1.0},
      {0,        3, .60854239843, 1.0},
      {0,        4,1.23254239843, 1.0}
    },//N
    {
      {0,        0, 1517.8667506, 1.0},
      {0,        0, 489.67952008, 1.0},
      {0,        0, 176.72118665, 1.0},
      {0,        0, 63.792233137, 1.0},
      {0,        0, 25.366499130, 1.0},
      {0,        0, 9.9135491200, 1.0},
      {0,        0, 4.4645306584, 1.0},
      {0,        0, 1.8017743661, 1.0},
      {0,        0,0.80789710965, 1.0},
      {0,        0,0.33864326862, 1.0},
      {0,        1, 120.16030921, 1.0},
      {0,        1, 34.409622474, 1.0},
      {0,        1, 12.581148610, 1.0},
      {0,        1, 5.0663824249, 1.0},
      {0,        1, 2.0346927092, 1.0},
      {0,        1,0.86092967212, 1.0},
      {0,        1,0.36681356726, 1.0},
      {0,        2, 19.043062805, 1.0},
      {0,        2, 5.8060381104, 1.0},
      {0,        2, 2.1891841580, 1.0},
      {0,        2,0.87794613558, 1.0},
      {0,        2,0.35623646700, 1.0},
      {0,        3,2.493914788135, 1.0},
      {0,        3, .824914788135, 1.0},
      {0,        4,1.607914788135, 1.0}
    },//O
    {},//F
    {} //Ne
  }
);

inline const std::vector<std::vector<primitive>> QZVP_JKfit(
  {
    { //center  type  exp         coef
      {0,       0, 9.5302493327,  1.0},
      {0,       0, 1.9174506246,  1.0},
      {0,       0, 0.68424049142, 1.0},
      {0,       0, 0.28413255710, 1.0},
      {0,       1, 2.9133232035,  1.0},
      {0,       1, 1.2621205398,  1.0},
      {0,       1, 0.50199775874, 1.0},
      {0,       2, 2.8832083931,  1.0},
      {0,       2, 1.2801701725 , 1.0},
      {0,       2, 0.52511317770, 1.0},
      {0,       3, 2.7489448439 , 1.0},
      {0,       3, 1.1900885456 , 1.0},
      {0,       4, 1.4752662714 , 1.0}
    },//H
    {},//He
    {},//Li
    {},//Be
    {},//B
    {
      {0,        0, 1113.9867719, 1.0},
      {0,        0, 369.16234180, 1.0},
      {0,        0, 121.79275232, 1.0},
      {0,        0, 48.127114540, 1.0},
      {0,        0, 20.365074004, 1.0},
      {0,        0, 8.0883596856, 1.0},
      {0,        0, 2.5068656570, 1.0},
      {0,        0, 1.2438537380, 1.0},
      {0,        0,0.48449899601, 1.0},
      {0,        0,0.19185160296, 1.0},
      {0,        1, 102.99176249, 1.0},
      {0,        1, 28.132594009, 1.0},
      {0,        1, 9.8364318173, 1.0},
      {0,        1, 3.3490544980, 1.0},
      {0,        1, 1.4947618613, 1.0},
      {0,        1,0.57690108899, 1.0},
      {0,        1,0.20320063291, 1.0},
      {0,        2, 10.594068356, 1.0},
      {0,        2, 3.5997195366, 1.0},
      {0,        2, 1.3355691094, 1.0},
      {0,        2,0.51949764954, 1.0},
      {0,        2,0.19954125200, 1.0},
      {0,        3,1.95390000000, 1.0},
      {0,        3, .75490000000, 1.0},
      {0,        3, .33390000000, 1.0},
      {0,        4,1.52490000000, 1.0},
      {0,        4, .59090000000, 1.0},
      {0,        5,1.11690000000, 1.0}
    },//C
    {
      {0,        0, 1102.8622453, 1.0},
      {0,        0, 370.98041153, 1.0},
      {0,        0, 136.73555938, 1.0},
      {0,        0, 50.755871924, 1.0},
      {0,        0, 20.535656095, 1.0},
      {0,        0, 7.8318737184, 1.0},
      {0,        0, 3.4784063855, 1.0},
      {0,        0, 1.4552856603, 1.0},
      {0,        0,0.63068989071, 1.0},
      {0,        0,0.27276596483, 1.0},
      {0,        1, 93.540954073, 1.0},
      {0,        1, 29.524019527, 1.0},
      {0,        1, 10.917502987, 1.0},
      {0,        1, 4.3449288991, 1.0},
      {0,        1, 1.8216912640, 1.0},
      {0,        1,0.75792424494, 1.0},
      {0,        1,0.28241469033, 1.0},
      {0,        2, 16.419378926, 1.0},
      {0,        2, 5.0104049385, 1.0},
      {0,        2, 1.9793971884, 1.0},
      {0,        2,0.78495771518, 1.0},
      {0,        2,0.28954065963, 1.0},
      {0,        3,2.98600000000, 1.0},
      {0,        3,1.11700000000, 1.0},
      {0,        3, .48400000000, 1.0},
      {0,        4,2.17600000000, 1.0},
      {0,        4, .83400000000, 1.0},
      {0,        5,1.57600000000, 1.0}
    },//N
    {
      {0,        0, 1517.8667506, 1.0},
      {0,        0, 489.67952008, 1.0},
      {0,        0, 176.72118665, 1.0},
      {0,        0, 63.792233137, 1.0},
      {0,        0, 25.366499130, 1.0},
      {0,        0, 9.9135491200, 1.0},
      {0,        0, 4.4645306584, 1.0},
      {0,        0, 1.8017743661, 1.0},
      {0,        0,0.80789710965, 1.0},
      {0,        0,0.33864326862, 1.0},
      {0,        1, 120.16030921, 1.0},
      {0,        1, 34.409622474, 1.0},
      {0,        1, 12.581148610, 1.0},
      {0,        1, 5.0663824249, 1.0},
      {0,        1, 2.0346927092, 1.0},
      {0,        1,0.86092967212, 1.0},
      {0,        1,0.36681356726, 1.0},
      {0,        2, 19.043062805, 1.0},
      {0,        2, 5.8060381104, 1.0},
      {0,        2, 2.1891841580, 1.0},
      {0,        2,0.87794613558, 1.0},
      {0,        2,0.35623646700, 1.0},
      {0,        3,3.96585000000, 1.0},
      {0,        3,1.49085000000, 1.0},
      {0,        3, .63485000000, 1.0},
      {0,        4,2.85685000000, 1.0},
      {0,        4,1.04985000000, 1.0},
      {0,        5,2.03685000000, 1.0}
    },//O
    {},//F
    {} //Ne
  }
);

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

inline bool ends_with(const std::string& str, const std::string& suffix) {
    if (str.length() >= suffix.length()) {
        return str.compare(str.length() - suffix.length(), suffix.length(), suffix) == 0;
    }
    return false;
}

template <typename T>
inline bool is_nan(T in) {
  if (typeid(in) == typeid(double)) {
    return in != in;
  }
  else if (typeid(in) == typeid(float)) {
    return in != in;
  }
  else if (typeid(in) == typeid(long double)) {
    return in != in;
  }
  else {
    return false;
  }
}

#include "wfn_class.h"
#include "atoms.h"