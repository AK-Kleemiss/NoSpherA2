#pragma once
#include "pch.h"
// Pre-definition of classes included later
class WFN;
class cell;
class atom;
class BasisSet;

inline std::streambuf *coutbuf = std::cout.rdbuf(); // save old buf
void error_check(const bool condition, const std::source_location loc, const std::string &error_mesasge, std::ostream &log_file = std::cout);
void not_implemented(const std::source_location loc, const std::string &error_mesasge, std::ostream &log_file);
#define err_checkf(condition, error_message, file) error_check(condition, std::source_location::current(), error_message, file)
#define err(error_message, file) error_check(false, std::source_location::current(), error_message, file)
#define err_chkf(condition, error_message, file) error_check(condition, std::source_location::current(), error_message, file)
#define err_chekf(condition, error_message, file) error_check(condition, std::source_location::current(), error_message, file)
#define err_not_impl_f(error_message, file) not_implemented(std::source_location::current(), error_message, file)
#ifdef _WIN32
#define ivdep loop(ivdep)
#else
#define ivdep GCC ivdep
#endif

typedef std::complex<double> cdouble;
typedef std::vector<double> vec;
typedef std::vector<vec> vec2;
typedef std::vector<vec2> vec3;
typedef std::vector<int> ivec;
typedef std::vector<ivec> ivec2;
typedef std::vector<cdouble> cvec;
typedef std::vector<cvec> cvec2;
typedef std::vector<cvec2> cvec3;
typedef std::vector<std::vector<cvec2>> cvec4;
typedef std::vector<bool> bvec;
typedef std::vector<std::string> svec;
typedef std::vector<std::filesystem::path> pathvec;
typedef std::chrono::high_resolution_clock::time_point _time_point;
typedef Kokkos::Experimental::mdarray<double, Kokkos::extents<unsigned long long, std::dynamic_extent>> dMatrix1;
typedef Kokkos::Experimental::mdarray<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent>> cMatrix1;
typedef Kokkos::Experimental::mdarray<int, Kokkos::extents<unsigned long long, std::dynamic_extent>> iMatrix1;
typedef Kokkos::Experimental::mdarray<bool, Kokkos::extents<unsigned long long, std::dynamic_extent>> bMatrix1;
typedef Kokkos::Experimental::mdarray<double, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> dMatrix2;
typedef Kokkos::Experimental::mdarray<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> cMatrix2;
typedef Kokkos::Experimental::mdarray<int, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> iMatrix2;
typedef Kokkos::Experimental::mdarray<bool, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> bMatrix2;
typedef Kokkos::Experimental::mdarray<double, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> dMatrix3;
typedef Kokkos::Experimental::mdarray<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> cMatrix3;
typedef Kokkos::Experimental::mdarray<int, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> iMatrix3;
typedef Kokkos::Experimental::mdarray<bool, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> bMatrix3;
typedef Kokkos::Experimental::mdarray<double, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> dMatrix4;
typedef Kokkos::Experimental::mdarray<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> cMatrix4;
typedef Kokkos::Experimental::mdarray<int, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> iMatrix4;
typedef Kokkos::Experimental::mdarray<bool, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> bMatrix4;

typedef Kokkos::mdspan<double, Kokkos::extents<unsigned long long, std::dynamic_extent>> dMatrixRef1;
typedef Kokkos::mdspan<double, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> dMatrixRef2;

typedef Kokkos::mdspan<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent>> cMatrixRef1;
typedef Kokkos::mdspan<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> cMatrixRef2;
typedef Kokkos::mdspan<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> cMatrixRef3;
typedef Kokkos::mdspan<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> cMatrixRef4;
typedef Kokkos::mdspan<const cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> ccMatrixRef4;

int vec_sum(const bvec &in);
int vec_sum(const ivec &in);
double vec_sum(const vec &in);
cdouble vec_sum(const cvec &in);
double vec_length(const vec &in);
template <typename array>
const double array_length(const array &in)
{
    double sum = 0.0;
    for (double val : in)
    {
        sum += val * val;
    }
    return sqrt(sum);
}

// Function to compute cross product
std::array<double, 3> cross(const std::array<double, 3> &a, const std::array<double, 3> &b);

// Function to compute dot product
double a_dot(const std::array<double, 3> &a, const std::array<double, 3> &b);

constexpr const std::complex<double> c_one(0, 1.0);

extern std::string help_message;
std::string NoSpherA2_message(bool no_date = false);
extern std::string build_date;

namespace sha
{
// Rotate right operation
#define ROTR(x, n) ((x >> n) | (x << (32 - n)))

// Logical functions for SHA-256
#define CH(x, y, z) ((x & y) ^ (~x & z))
#define MAJ(x, y, z) ((x & y) ^ (x & z) ^ (y & z))
#define EP0(x) (ROTR(x, 2) ^ ROTR(x, 13) ^ ROTR(x, 22))
#define EP1(x) (ROTR(x, 6) ^ ROTR(x, 11) ^ ROTR(x, 25))
#define SIG0(x) (ROTR(x, 7) ^ ROTR(x, 18) ^ (x >> 3))
#define SIG1(x) (ROTR(x, 17) ^ ROTR(x, 19) ^ (x >> 10))

    constexpr uint64_t custom_bswap_64(uint64_t x)
    {
        return ((x & 0xFF00000000000000ull) >> 56) |
               ((x & 0x00FF000000000000ull) >> 40) |
               ((x & 0x0000FF0000000000ull) >> 24) |
               ((x & 0x000000FF00000000ull) >> 8) |
               ((x & 0x00000000FF000000ull) << 8) |
               ((x & 0x0000000000FF0000ull) << 24) |
               ((x & 0x000000000000FF00ull) << 40) |
               ((x & 0x00000000000000FFull) << 56);
    }

    constexpr uint32_t custom_bswap_32(uint32_t value)
    {
        return ((value & 0x000000FF) << 24) |
               ((value & 0x0000FF00) << 8) |
               ((value & 0x00FF0000) >> 8) |
               ((value & 0xFF000000) >> 24);
    }

    // Initial hash values
    constexpr uint32_t k[64] = {
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
        0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
        0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
        0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
        0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
        0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
        0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
        0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
        0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
        0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
        0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
        0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
        0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};

    // SHA-256 processing function
    void sha256_transform(uint32_t state[8], const uint8_t block[64]);

    // SHA-256 update function
    void sha256_update(uint32_t state[8], uint8_t buffer[64], const uint8_t *data, size_t len, uint64_t &bitlen);

    // SHA-256 padding and final hash computation
    void sha256_final(uint32_t state[8], uint8_t buffer[64], uint64_t bitlen, uint8_t hash[32]);

    // Function to calculate SHA-256 hash
    std::string sha256(const std::string &input);
}

bool is_similar_rel(const double &first, const double &second, const double &tolerance);
bool is_similar(const double &first, const double &second, const double &tolerance);
bool is_similar_abs(const double &first, const double &second, const double &tolerance);
std::filesystem::path get_home_path(void);
char asciitolower(char in);

bool generate_sph2cart_mat(vec2 &p, vec2 &d, vec2 &f, vec2 &g);
bool generate_cart2sph_mat(vec2 &d, vec2 &f, vec2 &g, vec2 &h);
std::string go_get_string(std::ifstream &file, std::string search, bool rewind = true);

const int sht2nbas(const int &type);

const int shell2function(const int &type, const int &prim);

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

void remove_empty_elements(svec &input, const std::string &empty = " ");
std::chrono::high_resolution_clock::time_point get_time();
long long int get_musec(std::chrono::high_resolution_clock::time_point start, std::chrono::high_resolution_clock::time_point end);
long long int get_msec(std::chrono::high_resolution_clock::time_point start, std::chrono::high_resolution_clock::time_point end);
long long int get_sec(std::chrono::high_resolution_clock::time_point start, std::chrono::high_resolution_clock::time_point end);

void write_timing_to_file(std::ostream &file, std::vector<_time_point> time_points, std::vector<std::string> descriptions);

int CountWords(const char *str);

void copy_file(std::filesystem::path &from, std::filesystem::path &to);
std::string shrink_string(std::string &input);
std::string shrink_string_to_atom(std::string &input, const int &atom_number);
//------------------Functions to work with configuration files--------------------------
bool check_bohr(WFN &wave, bool debug);

bool open_file_dialog(std::filesystem::path &path, bool debug, std::vector <std::string> filter, const std::string& current_path);
bool save_file_dialog(std::filesystem::path &path, bool debug, const svec &endings, const std::string &filename_given = "", const std::string& current_path = "");
void select_cubes(std::vector<std::vector<unsigned int>> &selection, std::vector<WFN> &wavy, unsigned int nr_of_cubes = 1, bool wfnonly = false, bool debug = false);
bool unsaved_files(std::vector<WFN> &wavy);

std::string trim(const std::string &s);

inline void print_centered_text(const std::string &text, int &bar_width)
{
    const int text_length = static_cast<int>(text.length());
    const int total_padding = bar_width - text_length;
    const int padding_left = total_padding / 2;
    const int padding_right = total_padding - padding_left;

    std::cout << "["
              << std::setw(padding_left) << std::setfill(' ') << ""
              << text
              << std::setw(padding_right) << std::setfill(' ') << ""
              << "]" << std::endl;
}

//-------------------------Progress_bar--------------------------------------------------
// LMS: My implementation of a progress bar, I would like it to stay within one line that is compatible with parallel loops
class ProgressBar
{
public:
    ~ProgressBar();

    ProgressBar(const unsigned long long &worksize, const int &bar_width = 60, const std::string &fill = "#", const std::string &remainder = " ", const std::string &status_text = "")
        : worksize_(worksize), bar_width_(bar_width), fill_(fill), remainder_(remainder), status_text_(status_text), workdone(0), progress_(0.0f), workpart_(100.0f / worksize), percent_((worksize / 100 > 1) ? worksize / 100 : 1)
    {
        int bw = bar_width_ + 2;
        // Write status text
        print_centered_text(status_text_, bw);
        linestart = std::cout.tellp();
#ifdef _WIN32
        initialize_taskbar_progress();
#endif
    }

    void set_progress()
    {
        progress_ = (float)workdone * workpart_;
    }

    void update(std::ostream &os = std::cout)
    {
#pragma omp critical
        {
            workdone += 1;
            if (workdone % percent_ == 0)
            {
                set_progress();
                write_progress(os);
            }
        }
    }

    void write_progress(std::ostream &os = std::cout);

private:
    const unsigned long long worksize_;
    const float workpart_;
    const unsigned long long percent_;
    int bar_width_;
    std::string fill_;
    std::string remainder_;
    std::string status_text_;
    std::atomic<unsigned long long> workdone;
    float progress_;
    std::streampos linestart;
#ifdef _WIN32
    ITaskbarList3 *taskbarList_;

    void initialize_taskbar_progress();
#endif
};

void readxyzMinMax_fromWFN(
    WFN &wavy,
    double *CoordMinMax,
    int *NbSteps,
    double Radius,
    double Increments,
    bool no_bohr = false);

void readxyzMinMax_fromCIF(
    std::filesystem::path cif,
    double *CoordMinMax,
    int *NbSteps,
    vec2 &cm,
    double Resolution);

bool read_fracs_ADPs_from_CIF(std::filesystem::path &cif, WFN &wavy, cell &unit_cell, std::ofstream &log3, bool debug);

double double_from_string_with_esd(std::string in);

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

constexpr unsigned int doublefactorial(int n)
{
    if (n <= 1)
        return 1;
    return n * doublefactorial(n - 2);
}

template <typename T>
void removeElement(std::vector<T> &vec, const T &x)
{
    // Use std::remove to shift elements and get the new end iterator
    auto new_end = std::remove(vec.begin(), vec.end(), x);
    // Erase the elements from the new end to the actual end
    vec.erase(new_end, vec.end());
}

inline void Enter(){
	std::cout << "press ENTER to continue... " << std::flush;
	std::cin.ignore();
	std::cin.get();
};

inline void cls(){
//   std::cout << string( 100, '\n' );
#ifdef _WIN32
	if(system("CLS")) std::cout << "this should not happen...!" << std::endl;
#else
	if(system("clear")) std::cout << "this should not happen...!" << std::endl;
#endif
};

inline bool yesno(){
	bool end=false;
	while (!end) {
		char dum ='?';
		std::cout << "(Y/N)?";
		std::cin >> dum;
		if(dum == 'y'||dum == 'Y'){
			std::cout << "Okay..." << std::endl;
			 return true;
		}
		else if(dum == 'N'||dum == 'n') return false;
		else std::cout << "Sorry, i did not understand that!" << std::endl;
	}
	return false;
};

class primitive
{
private:
    int center, type;
    double exp, coefficient;
    double norm_const = -10;

public:
    void normalize()
    {
        coefficient *= normalization_constant();
    };
    void unnormalize()
    {
        coefficient /= normalization_constant();
    };
    double normalization_constant() const
    {
        // assuming type is equal to angular momentum
        return norm_const;
    }
    primitive() : center(0), type(0), exp(0.0), coefficient(0.0) {};
    primitive(int c, int t, double e, double coef);
    bool operator==(const primitive &other) const
    {
        return center == other.center &&
               type == other.type &&
               exp == other.exp &&
               coefficient == other.coefficient;
    };
    int get_center() const
    {
        return center;
    };
    int get_type() const
    {
        return type;
    };
    double get_exp() const
    {
        return exp;
    };
    double get_coef() const
    {
        return coefficient;
    };
    void set_center(const int &c)
    {
        center = c;
    };
    void set_type(const int &t)
    {
        type = t;
    };
    void set_exp(const double &e)
    {
        exp = e;
    };
    void set_coef(const double &c)
    {
        coefficient = c;
    };
};

struct ECP_primitive : primitive
{
    int n;
    ECP_primitive() : primitive(), n(0) {}
    ECP_primitive(int c, int t, double e, double coef, int n) : primitive(c, t, e, coef), n(n) {}
};

typedef std::array<int, 3> hkl_t;
typedef std::set<hkl_t> hkl_list;
typedef std::set<hkl_t>::const_iterator hkl_list_it;

typedef std::array<double, 3> hkl_d;
typedef std::set<hkl_d> hkl_list_d;
typedef std::set<hkl_d>::const_iterator hkl_list_it_d;

//---------------- Object for handling all input options -------------------------------
struct options
/**
 * @brief The `options` class represents a collection of options and settings for a program.
 *
 * It contains various member variables that store different configuration parameters.
 * These parameters control the behavior and functionality of the program.
 *
 * The `options` class also provides constructors and member functions to initialize and manipulate these parameters.
 *
 * @note This class is used to configure the behavior of a specific program and may have different member variables and functions depending on the program's requirements.
 */
{
    std::ostream &log_file;
    double resolution = 0.1;
    double radius = 2.0;
    double d_sfac_scan = 0.0;
    double sfac_diffuse = 0.0;
    double dmin = 99.0;
    double mem = 1000.0; // In MB
    double efield = 0.005;
    double MinMax[6]{0, 0, 0, 0, 0, 0};
    ivec MOs;
    ivec2 groups;
    ivec2 hkl_min_max{{-100, 100}, {-100, 100}, {-100, 100}};
    vec2 twin_law;
    ivec2 combined_tsc_groups;
    pathvec combined_tsc_calc_files;
    pathvec combined_tsc_calc_cifs;
    std::vector<unsigned int> combined_tsc_calc_mult;
    ivec combined_tsc_calc_charge;
    ivec combined_tsc_calc_ECP;
    svec arguments;
    pathvec combine_mo;
    svec Cations;
    svec Anions;
    pathvec pol_wfns;
    ivec cmo1;
    ivec cmo2;
    ivec ignore;
    std::filesystem::path SALTED_DIR;
    std::string SALTED_DFBASIS;
    std::filesystem::path wfn;
    std::filesystem::path wfn2;
    std::filesystem::path fchk;
    std::string basis_set;
    std::filesystem::path hkl;
    std::filesystem::path cif;
    std::string method;
    std::filesystem::path xyz_file;
    std::filesystem::path coef_file;
    std::filesystem::path hirshfeld_surface;
    std::filesystem::path hirshfeld_surface2;
    std::string fract_name;
    std::filesystem::path wavename;
    std::filesystem::path gaussian_path;
    std::filesystem::path turbomole_path;
    std::filesystem::path basis_set_path;
    std::string cwd;
    bool debug = false;
    bool rho = false;
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
    bool SALTED = false;
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
    bool all_mos = false;
    bool test = false;
    bool becke = false;
    bool electron_diffraction = false;
    bool ECP = false;
    bool RI_FIT = false;
    bool needs_Thakkar_fill = false;
    bool qct = false;
    int hirsh_number = 0;
    int NbSteps[3]{0, 0, 0};
    int accuracy = 2;
    int threads = -1;
    int pbc = 0;
    int charge = 0;
    int ECP_mode = 0;
    unsigned int mult = 0;
    hkl_list m_hkl_list;

    /**
     * @brief Looks for debug mode in the command line arguments.
     *
     * This function searches for a specific debug flag in the command line arguments
     * and modifies the `argc` and `argv` parameters accordingly and stores them internally
     *
     * @param argc The number of command line arguments.
     * @param argv An array of C-style strings representing the command line arguments.
     *
     */
    void look_for_debug(int &argc, char **argv);
    /**
     * @brief Digests the options.
     *
     * This function is responsible for digesting the options and performing the necessary actions based on the selected options.
     *
     * @note Make sure to call this function after looking for debug.
     */
    void digest_options();

    options() : log_file(std::cout)
    {
        groups.resize(1);
    };
    options(int &argc, char **argv, std::ostream &log) : log_file(log)
    {
        groups.resize(1);
        look_for_debug(argc, argv);
    };
};

const double gaussian_radial(primitive &p, double &r);

int load_basis_into_WFN(WFN &wavy, const std::array<std::vector<primitive>, 118> &b);
int load_basis_into_WFN(WFN &wavy, BasisSet &b);

double hypergeometric(double a, double b, double c, double x);

cdouble hypergeometric(double a, double b, double c, cdouble x);

bool ends_with(const std::string &str, const std::string &suffix);

bool is_nan(double &in);
bool is_nan(float &in);
bool is_nan(long double &in);
bool is_nan(cdouble &in);

bool read_block_from_fortran_binary(std::ifstream &file, void *Target);
/*
#ifdef _WIN32
bool ExtractDLL(const std::filesystem::path &dllName);
bool check_OpenBLAS_DLL(const bool &debug = false);
#endif
*/