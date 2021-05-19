#ifndef __CONVENIENCE_H__
#define __CONVENIENCE_H__

#include <iostream>
#include <string>
#include <sstream>
#include <sys/stat.h>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cctype>

class WFN;
class cell;

//------------------general functions for easy use of terminal input--------------------
const double bragg_angstrom[114]{
0.00, 
	0.35, 0.35, 
	1.45, 1.05,																																					0.85, 0.70, 0.65, 0.60, 0.50, 0.45,
	1.80, 1.50,																																					1.25, 1.10, 1.00, 1.00, 1.00, 1.00,
	2.20, 1.80,																						1.60, 1.40, 1.35, 1.40, 1.40, 1.40, 1.35, 1.35, 1.35, 1.35, 1.30, 1.25, 1.15, 1.15, 1.15, 1.10,
	2.35, 2.00,																						1.80, 1.55, 1.45, 1.45, 1.35, 1.30, 1.35, 1.40, 1.60, 1.55, 1.55, 1.45, 1.45, 1.40, 1.40, 1.40,
	2.60, 2.15, 1.95, 1.85, 1.85, 1.85, 1.85, 1.85, 1.85, 1.80, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.55, 1.45, 1.35, 1.30, 1.30, 1.35, 1.35, 1.35, 1.50, 1.90, 1.75, 1.60, 1.90, 1.50, 1.50,
	2.80, 2.35, 2.15, 2.05, 2.05, 2.05, 2.05, 2.05, 2.05, 2.00, 1.95, 1.95, 1.95, 1.95, 1.95, 1.95};
bool yesno();
bool is_similar_rel(double first, double second, double tolerance);
bool is_similar(double first, double second, double tolerance);
void Enter();
void cls();
std::string get_home_path(void);
void join_path(std::string &s1, std::string &s2);

inline bool generate_sph2cart_mat(std::vector<std::vector<double>> &d, std::vector<std::vector<double>> &f, std::vector<std::vector<double>> &g, std::vector<std::vector<double>> &h) {
	//                                                   
	//From 5D: D 0, D + 1, D - 1, D + 2, D - 2           
	//To 6D : 1  2  3  4  5  6                           
	//XX, YY, ZZ, XY, XZ, YZ      
	// 
	d.resize(6);
#pragma omp parallel for
	for (int i = 0; i < 6; i++) {
		d[i].resize(5);
		std::fill(d[i].begin(), d[i].end(), 0.0);
	}
	//D0 = -0.5 * XX - 0.5 * YY + ZZ
	d[0][0] = -0.5;
	d[1][0] = -0.5;
	d[2][0] = 1.0;
	//D + 1 = XZ
	d[4][1] = 1.0;
	//D - 1 = YZ
	d[5][2] = 1.0;
	//D + 2 = SQRT(3) / 2 * (XX - YY)
	d[0][3] = sqrt(3.0) / 2.0;
	d[1][3] = -sqrt(3.0) / 2.0;
	//D - 2 = XY
	d[3][4] = 1.0;

	//From 7F: F 0, F + 1, F - 1, F + 2, F - 2, F + 3, F - 3
	//To 10F : 1   2   3   4   5   6   7   8   9  10
	//XXX, YYY, ZZZ, XYY, XXY, XXZ, XZZ, YZZ, YYZ, XYZ(Gaussian sequence, not identical to Multiwfn)
	//
	f.resize(10);
#pragma omp parallel for
	for (int i = 0; i < 10; i++) {
		f[i].resize(7);
		std::fill(f[i].begin(), f[i].end(), 0.0);
	}
	//F 0 = -3 / (2 * sqrt5) * (XXZ + YYZ) + ZZZ
	f[2][0] = 1.0;
	f[5][0] = -1.5 / sqrt(5.0);
	f[8][0] = -1.5 / sqrt(5.0);
	//F + 1 = -sqrt(3 / 8) * XXX - sqrt(3 / 40) * XYY + sqrt(6 / 5) * XZZ
	f[0][1] = -sqrt(3.0 / 8.0);
	f[3][1] = -sqrt(3.0 / 40.0);
	f[6][1] = sqrt(6.0 / 5.0);
	//F - 1 = -sqrt(3 / 40) * XXY - sqrt(3 / 8) * YYY + sqrt(6 / 5) * YZZ
	f[1][2] = -sqrt(3.0 / 8.0);
	f[4][2] = -sqrt(3.0 / 40.0);
	f[7][2] = sqrt(6.0 / 5.0);
	//F + 2 = sqrt3 / 2 * (XXZ - YYZ)
	f[5][3] = sqrt(3.0) / 2.0;
	f[8][3] = -sqrt(3.0) / 2.0;
	//F - 2 = XYZ
	f[9][4] = 1.0;
	//F + 3 = sqrt(5 / 8) * XXX - 3 / sqrt8 * XYY
	f[0][5] = sqrt(5.0 / 8.0);
	f[3][5] = -3.0 / sqrt(8.0);
	//F - 3 = 3 / sqrt8 * XXY - sqrt(5 / 8) * YYY
	f[1][6] = -sqrt(5.0 / 8.0);
	f[4][6] = 3.0 / sqrt(8.0);

	//From 9G: G 0, G + 1, G - 1, G + 2, G - 2, G + 3, G - 3, G + 4, G - 4
	//To 15G : 1    2    3    4    5    6    7    8
	//ZZZZ, YZZZ, YYZZ, YYYZ, YYYY, XZZZ, XYZZ, XYYZ
	//9   10   11   12   13   14   15
	//XYYY, XXZZ, XXYZ, XXYY, XXXZ, XXXY, XXXX
	//
	g.resize(15);
#pragma omp parallel for
	for (int i = 0; i < 15; i++) {
		g[i].resize(9);
		std::fill(g[i].begin(), g[i].end(), 0.0);
	}
	//G 0 = ZZZZ + 3 / 8 * (XXXX + YYYY) - 3 * sqrt(3 / 35) * (XXZZ + YYZZ - 1 / 4 * XXYY)
	g[0][0] = 1.0;
	g[2][0] = -3.0 * sqrt(3.0 / 35.0);
	g[4][0] = 3.0 / 8.0;
	g[9][0] = -3.0 * sqrt(3.0 / 35.0);
	g[11][0] = 3.0 / 4.0 * sqrt(3.0 / 35.0);
	g[14][0] = 3.0 / 8.0;
	//G + 1 = 2 * sqrt(5 / 14) * XZZZ - 3 / 2 * sqrt(5 / 14) * XXXZ - 3 / 2 / sqrt14 * XYYZ
	g[5][1] = 2.0 * sqrt(5.0 / 14.0);
	g[7][1] = -1.5 / sqrt(14.0);
	g[12][1] = -1.5 * sqrt(5.0 / 14.0);
	//G - 1 = 2 * sqrt(5 / 14) * YZZZ - 3 / 2 * sqrt(5 / 14) * YYYZ - 3 / 2 / sqrt14 * XXYZ
	g[1][2] = 2.0 * sqrt(5.0 / 14.0);
	g[3][2] = -1.5 * sqrt(5.0 / 14.0);
	g[10][2] = -1.5 / sqrt(14.0);
	//G + 2 = 3 * sqrt(3 / 28) * (XXZZ - YYZZ) - sqrt5 / 4 * (XXXX - YYYY)
	g[2][3] = -3.0 * sqrt(3.0 / 28.0);
	g[4][3] = sqrt(5.0) / 4.0;
	g[9][3] = 3.0 * sqrt(3.0 / 28.0);
	g[14][3] = -sqrt(5.0) / 4.0;
	//G - 2 = 3 / sqrt7 * XYZZ - sqrt(5 / 28) * (XXXY + XYYY)
	g[6][4] = 3.0 / sqrt(7.0);
	g[8][4] = -sqrt(5.0 / 28.0);
	g[13][4] = -sqrt(5.0 / 28.0);
	//G + 3 = sqrt(5 / 8) * XXXZ - 3 / sqrt8 * XYYZ
	g[7][5] = -3.0 / sqrt(8.0);
	g[12][5] = sqrt(5.0 / 8.0);
	//G - 3 = -sqrt(5 / 8) * YYYZ + 3 / sqrt8 * XXYZ
	g[3][6] = -sqrt(5.0 / 8.0);
	g[10][6] = 3.0 / sqrt(8.0);
	//G + 4 = sqrt35 / 8 * (XXXX + YYYY) - 3 / 4 * sqrt3 * XXYY
	g[4][7] = sqrt(35.0) / 8.0;
	g[11][7] = -3.0 / 4.0 * sqrt(3.0);
	g[14][7] = sqrt(35.0) / 8.0;
	//G - 4 = sqrt5 / 2 * (XXXY - XYYY)
	g[8][8] = -sqrt(5.0) / 2.0;
	g[13][8] = sqrt(5.0) / 2.0;

	//From 11H: H 0, H + 1, H - 1, H + 2, H - 2, H + 3, H - 3, H + 4, H - 4, H + 5, H - 5
	//To 21H : 1     2     3     4     5     6     7     8     9    10
	//ZZZZZ YZZZZ YYZZZ YYYZZ YYYYZ YYYYY XZZZZ XYZZZ XYYZZ XYYYZ
	//11    12    13    14    15    16    17    18    19    20    21
	//XYYYY XXZZZ XXYZZ XXYYZ XXYYY XXXZZ XXXYZ XXXYY XXXXZ XXXXY XXXXX
	//
	h.resize(21);
#pragma omp parallel for
	for (int i = 0; i < 21; i++) {
		h[i].resize(11);
		std::fill(h[i].begin(), h[i].end(), 0.0);
	}
	//H 0 = ZZZZZ - 5 / sqrt21 * (XXZZZ + YYZZZ) + 5 / 8 * (XXXXZ + YYYYZ) + sqrt(15 / 7) / 4 * XXYYZ
	h[0][0] = 1.0;
	h[11][0] = -5.0 / sqrt(21.0);
	h[2][0] = -5.0 / sqrt(21.0);
	h[18][0] = 5.0 / 8.0;
	h[4][0] = 5.0 / 8.0;
	h[13][0] = sqrt(15.0 / 7.0) / 4.0;
	//H + 1 = sqrt(5 / 3) * XZZZZ - 3 * sqrt(5 / 28) * XXXZZ - 3 / sqrt28 * XYYZZ + sqrt15 / 8 * XXXXX + sqrt(5 / 3) / 8 * XYYYY + sqrt(5 / 7) / 4 * XXXYY
	h[6][1] = sqrt(5.0 / 3.0);
	h[15][1] = -3.0 * sqrt(5.0 / 28.0);
	h[8][1] = -3.0 / sqrt(28.0);
	h[20][1] = sqrt(15.0) / 8.0;
	h[10][1] = sqrt(5.0 / 3.0) / 8.0;
	h[17][1] = sqrt(5.0 / 7.0) / 4.0;
	//H - 1 = sqrt(5 / 3) * YZZZZ - 3 * sqrt(5 / 28) * YYYZZ - 3 / sqrt28 * XXYZZ + sqrt15 / 8 * YYYYY + sqrt(5 / 3) / 8 * XXXXY + sqrt(5 / 7) / 4 * XXYYY
	h[1][2] = sqrt(5.0 / 3.0);
	h[3][2] = -3.0 * sqrt(5.0 / 28.0);
	h[12][2] = -3.0 / sqrt(28.0);
	h[5][2] = sqrt(15.0) / 8.0;
	h[19][2] = sqrt(5.0 / 3.0) / 8.0;
	h[14][2] = sqrt(5.0 / 7.0) / 4.0;
	//H + 2 = sqrt5 / 2 * (XXZZZ - YYZZZ) - sqrt(35 / 3) / 4 * (XXXXZ - YYYYZ)
	h[11][3] = sqrt(5.0) / 2.0;
	h[2][3] = -sqrt(5.0) / 2.0;
	h[18][3] = -sqrt(35.0 / 3.0) / 4.0;
	h[4][3] = sqrt(35.0 / 3.0) / 4.0;
	//H - 2 = sqrt(5 / 3) * XYZZZ - sqrt(5 / 12) * (XXXYZ + XYYYZ)
	h[7][4] = sqrt(5.0 / 3.0);
	h[16][4] = -sqrt(5.0 / 12.0);
	h[9][4] = -sqrt(5.0 / 12.0);
	//H + 3 = sqrt(5 / 6) * XXXZZ - sqrt(3 / 2) * XYYZZ - sqrt(35 / 2) / 8 * (XXXXX - XYYYY) + sqrt(5 / 6) / 4 * XXXYY
	h[15][5] = sqrt(5.0 / 6.0);
	h[8][5] = -sqrt(1.5);
	h[20][5] = -sqrt(17.5) / 8.0;
	h[10][5] = sqrt(17.5) / 8.0;
	h[17][5] = sqrt(5.0 / 6.0) / 4.0;
	//H - 3 = -sqrt(5 / 6) * YYYZZ + sqrt(3 / 2) * XXYZZ - sqrt(35 / 2) / 8 * (XXXXY - YYYYY) - sqrt(5 / 6) / 4 * XXYYY
	h[3][6] = -sqrt(5.0 / 6.0);
	h[12][6] = sqrt(1.5);
	h[19][6] = -sqrt(17.5) / 8.0;
	h[5][6] = sqrt(17.5) / 8.0;
	h[14][6] = -sqrt(5.0 / 6.0) / 4.0;
	//H + 4 = sqrt35 / 8 * (XXXXZ + YYYYZ) - 3 / 4 * sqrt3 * XXYYZ
	h[18][7] = sqrt(35.0) / 8.0;
	h[4][7] = sqrt(35.0) / 8.0;
	h[13][7] = -0.75 * sqrt(3.0);
	//H - 4 = sqrt5 / 2 * (XXXYZ - XYYYZ)
	h[16][8] = sqrt(5.0) / 2.0;
	h[9][8] = -sqrt(5.0) / 2.0;
	//H + 5 = 3 / 8 * sqrt(7 / 2) * XXXXX + 5 / 8 * sqrt(7 / 2) * XYYYY - 5 / 4 * sqrt(3 / 2) * XXXYY
	h[20][9] = 3.0 / 8.0 * sqrt(3.5);
	h[10][9] = 5.0 / 8.0 * sqrt(3.5);
	h[17][9] = -1.25 * sqrt(1.5);
	//H - 5 = 3 / 8 * sqrt(7 / 2) * YYYYY + 5 / 8 * sqrt(7 / 2) * XXXXY - 5 / 4 * sqrt(3 / 2) * XXYYY
	h[5][10] = 3.0 / 8.0 * sqrt(3.5);
	h[19][10] = 5.0 / 8.0 * sqrt(3.5);
	h[14][10] = -1.25 * sqrt(1.5);
	return true;
}

std::string go_get_string(std::ifstream& file, std::string search, bool rewind = true);

inline const int sht2nbas(const int type) {
	const int st2bas[6]{ 1,3,6,10,15,21 };
	const int nst2bas[6]{ 11,9,7,5,4,1 };
	if (type >= 0)
		return st2bas[type];
	else
		return nst2bas[5+type];
};

inline const int shell2function(int type, int prim) {
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

const double normgauss(const int type, const double exp);

template<class T> std::string toString(const T& t)
{
     std::ostringstream stream;
     stream << t;
     return stream.str();
}

template<class T> T fromString(const std::string& s)
{
     std::istringstream stream (s);
     T t;
     stream >> t;
     return t;
}

inline void copyright();

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

inline bool copyright(std::ofstream& file);

inline bool exists(const std::string &name){
	std::ifstream f(name.c_str());
	return f.good();
};

std::string atnr2letter(const int nr);
void copy_file(std::string from, std::string to);
std::string shrink_string(std::string &input);
std::string shrink_string_to_atom(std::string &input, const int atom_number);
std::vector<std::string> split_string(const std::string& input, const std::string delimiter);
std::string get_filename_from_path(const std::string &input);
std::string get_foldername_from_path(const std::string &input);
std::string get_basename_without_ending(const std::string &input);
//------------------Functions to handle .wfn files--------------------------------------
bool writewfn(WFN &wavefunction, const std::string &path, bool &debug, bool occ);
bool readwfn(WFN &wavefunction, const std::string &path, bool &debug);
bool readwfn(WFN& wavefunction, const std::string& path, bool& debug, std::ofstream &file);
//------------------Functions to read from .fchk files----------------------------------
bool read_fchk_integer_block(std::ifstream& in, std::string heading, std::vector<int>& result, bool rewind = true);
bool read_fchk_double_block(std::ifstream& in, std::string heading, std::vector<double>& result, bool rewind = true);
int read_fchk_integer(std::string in);
int read_fchk_integer(std::ifstream& in, std::string search, bool rewind = true);
double read_fchk_double(std::string in);
double read_fchk_double(std::ifstream& in, std::string search, bool rewind = true);
//------------------Functions to work with configuration files--------------------------
void write_template_confi();
int program_confi(std::string &gaussian_path, std::string &turbomole_path, 
                   std::string &basis, int &ncpus, float &mem, bool debug = false, bool expert = false, unsigned int counter = 0);
bool check_bohr(WFN &wave, bool interactive, bool debug);
int filetype_identifier(std::string &file, bool debug = false);

/*bool open_file_dialog(std::string &path, bool debug, std::vector <std::string> filter);
bool save_file_dialog(std::string &path, bool debug, const std::vector<std::string> &endings, const std::string &filename_given);
bool save_file_dialog(std::string &path, bool debug, const std::vector<std::string> &endings);*/
void select_cubes(std::vector <std::vector <unsigned int> > &selection, std::vector<WFN> &wavy, unsigned int nr_of_cubes=1, bool wfnonly=false, bool debug = false);
bool unsaved_files(std::vector<WFN> &wavy);
int get_Z_from_label(const char * tmp);

inline int sum_of_bools(const std::vector<bool> in) {
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

inline void err(const char message, std::ofstream log) {
	log << message << std::endl;
	exit(-1);
}

inline void err(std::string& message, std::ofstream log) {
	log << message << std::endl;
	exit(-1);
}

inline void err(std::string& message) {
	std::cout << message << std::endl;
	exit(-1);
}


//-------------------------Progress_bar--------------------------------------------------

class progress_bar
{
	static const auto overhead = sizeof " [100%]";

	std::ofstream& os;
	const std::size_t bar_width;
	std::string message;
	const std::string full_bar;

public:
	progress_bar(std::ofstream& os, std::size_t line_width,
		std::string message_, const char symbol = '=')
		: os{ os },
		bar_width{ line_width - overhead },
		message{ std::move(message_) },
		full_bar{ std::string(bar_width, symbol) + std::string(bar_width, ' ') }
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

const double MPI2 = 2 * 3.14159265358979323844;

class cosinus_annaeherung
{
public:
	cosinus_annaeherung();
	inline double get(double x) const
	{
////////////////double xa = abs(x);
////////////////size_t pos = static_cast<size_t>((xa * mSize) / MPI2); // Stueststelle bestimmen (Wird fuer grosse X ungenau, aber passt fuer x
////////////////double dx = xa - pos * mStepwidth;
////////////////pos = pos % mSize; // Modulo, da sinus periodisch ist.
////////////////double y1 = mBase_values[pos];
////////////////double y2 = mBase_values[pos + 1];
////////////////return y1 + dx * (y2 - y1) / mStepwidth;
		return 0.0;
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

void readxyzMinMax_fromWFN(
	WFN& wavy,
	double* CoordMinMax,
	double* NbSteps,
	double Radius,
	double Increments);
void readxyzMinMax_fromCIF(
	std::string cif,
	double* CoordMinMax,
	double* NbSteps,
	std::vector < std::vector < double > > &cm,
	double Resolution,
	std::ofstream& file,
	bool debug = false);
void type2vector(
	const int index,
	int* vector);

bool read_fracs_ADPs_from_CIF(std::string cif, WFN& wavy, cell& unit_cell, std::ofstream& log3, bool debug);

inline double double_from_string_with_esd(std::string in) {
	if (in.find('(') == std::string::npos)
		return stod(in);
	else
		return stod(in.substr(0, in.find('(')));
}

#include "wfn_class.h"

#endif
