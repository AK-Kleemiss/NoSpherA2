#pragma once
#include "pch.h"
#include "throughput.h"
#include "tuning.h"



// Pre-definition of classes included later
class WFN;
class cell;
class atom;
class GridManager;
class BasisSet;
struct asym_atom;
enum PartitionType { Becke, TFVC, Hirshfeld, RI, MBIS, EMBIS };
enum class MultipoleScheme { TFVC, HIRSHFELD, MBIS, EMBIS, NUCLEAR, MULLIKEN, SANDERSON };
enum class RGBIOrbitalBasis { NAO, ANO };

void error_check(const bool condition, const std::source_location loc, const std::string& error_mesasge, std::ostream& log_file = std::cout);
void not_implemented(const std::source_location loc, const std::string& error_mesasge, std::ostream& log_file);
//The message is built only on failure: a literal past the SSO length cost a heap allocation per call, also in passing
//checks inside hot loops (get_atom_pos per grid point and atom)
#define err_checkf(condition, error_message, file) ((condition) ? (void)0 : error_check(false, std::source_location::current(), error_message, file))
#define err(error_message, file) error_check(false, std::source_location::current(), error_message, file)
#define err_not_impl_f(error_message, file) not_implemented(std::source_location::current(), error_message, file)

typedef std::complex<double> cdouble;
typedef std::vector<double> vec;
typedef std::vector<vec> vec2;
typedef std::vector<vec2> vec3;
typedef std::vector<int> ivec;
typedef std::vector<ivec> ivec2;
typedef std::vector<ivec2> ivec3;
typedef std::vector<cdouble> cvec;
typedef std::vector<cvec> cvec2;
typedef std::vector<cvec2> cvec3;
typedef std::vector<std::vector<cvec2>> cvec4;
typedef std::vector<bool> bvec;
typedef std::vector<bvec> bvec2;
typedef std::vector<bvec2> bvec3;
//operator[] fails naming the index instead of reading past the end, so a short line in a truncated file is an error, not UB
struct svec : std::vector<std::string>
{
	using std::vector<std::string>::vector;
	svec(const std::vector<std::string> &v) : std::vector<std::string>(v) {}
	svec(std::vector<std::string> &&v) : std::vector<std::string>(std::move(v)) {}
	const std::string &operator[](size_t i) const
	{
		if (i >= size())
			throw std::out_of_range("field " + std::to_string(i + 1) + " missing, the line has only " + std::to_string(size()));
		return std::vector<std::string>::operator[](i);
	}
	std::string &operator[](size_t i) { return const_cast<std::string &>(std::as_const(*this)[i]); }
};
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
typedef Kokkos::mdspan<double, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> dMatrixRef3;

typedef Kokkos::mdspan<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent>> cMatrixRef1;
typedef Kokkos::mdspan<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> cMatrixRef2;
typedef Kokkos::mdspan<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> cMatrixRef3;
typedef Kokkos::mdspan<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> cMatrixRef4;
typedef Kokkos::mdspan<const cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> ccMatrixRef4;

struct properties_options
{
	bool rho = false;
	bool eli = false;
	bool esp = false;
	bool elf = false;
	bool lap = false;
	bool rdg = false;
	bool hdef = false;
	bool def = false;
	bool hirsh = false;
	bool s_rho = false;
	bool all_mos = false;
	//Fukui functions f+/f-/f0 and the dual descriptor, frozen-orbital approximation
	bool fukui = false;
	//rho isosurface value (au) coloured by the ESP and written as obj; 0 = off
	double esp_isosurface = 0.0;
	double resolution = 0.1;
	double radius = 2.0;
	double integral_accuracy = -1;
	double promol_nci_rcut1 = 0.95;
	double promol_nci_rcut2 = 0.75;
	double promol_nci_rho_abs_max = 0.5;
	double promol_nci_rdg_max = 1.0;
	double promol_nci_colour_max = 0.015; // VMD/Olex2 colour range on sign(l2)rho in a.u., symmetric
	double promol_nci_iso = 0.5; // RDG of the _nci.obj surface
	//The _values.dat writer schedules grid points dynamically, so row order (not the values) varies between runs; forces single-threaded for golden files
	bool promol_nci_single_threaded = false;
	std::array<int, 3> NbSteps = { 0, 0, 0 };
	std::array<double, 6> MinMax = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	ivec MO_numbers;
	//-ibo_cube <list>: IBOs to write as cubes (ibo_selection syntax)
	std::string ibo_cube;
	int hirsh_number = 0;
	bool calc() const {
		return rho || eli || esp || elf || lap || rdg || hdef || def || hirsh || s_rho || all_mos || fukui || esp_isosurface > 0 || MO_numbers.size() > 0 || !ibo_cube.empty();
	}
	size_t n_grid_points() const {
		size_t result = static_cast<size_t>(NbSteps[0]) * NbSteps[1] * NbSteps[2];
		return result;
	}
};

typedef std::array<int, 3> i3;
typedef std::set<i3> hkl_list;
typedef std::set<i3>::const_iterator hkl_list_it;

typedef std::array<double, 3> d3;
typedef std::array<double, 4> d4;
typedef std::set<d3> hkl_list_d;
typedef std::set<d3>::const_iterator hkl_list_it_d;

struct I3Less {
	bool operator()(const i3& a, const i3& b) const noexcept {
		if (a[0] != b[0]) return a[0] < b[0];
		if (a[1] != b[1]) return a[1] < b[1];
		return a[2] < b[2];
	}
};

//indexed by the major grid point the additional points belong to, this map of RefinePoints has field value and the non-integer index in the grid as second
typedef std::multimap<i3, std::pair<double, d3>, I3Less> Refinepointmap;

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
extern std::string build_date;

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
	const double q = std::nearbyint(x * 0.63661977236758134308);
	const double r = ((x - q * 1.57079632673412561417e+00) - q * 6.07710050630396597660e-11) - q * 2.02226624871116645580e-21;
	const double z = r * r;
	const double sr = r + r * z * (-1.66666666666666324348e-01 + z * (8.33333333332248946124e-03 + z * (-1.98412698298579493134e-04 +
		z * (2.75573137070700676789e-06 + z * (-2.50507602534068634195e-08 + z * 1.58969099521155010221e-10)))));
	const double hz = 0.5 * z, w = 1.0 - hz;
	const double cr = w + (((1.0 - w) - hz) + z * z * (4.16666666666666019037e-02 + z * (-1.38888888888741095749e-03 +
		z * (2.48015872894767294178e-05 + z * (-2.75573143513906633035e-07 + z * (2.08757232129817482790e-09 +
		z * -1.13596475577881948265e-11))))));
	switch (static_cast<long long>(q) & 3) {
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
		0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2 };

	// SHA-256 processing function
	void sha256_transform(uint32_t state[8], const uint8_t block[64]);

	// SHA-256 update function
	void sha256_update(uint32_t state[8], uint8_t buffer[64], const uint8_t* data, size_t len, uint64_t& bitlen);

	// SHA-256 padding and final hash computation
	void sha256_final(uint32_t state[8], uint8_t buffer[64], uint64_t bitlen, uint8_t hash[32]);

	// Function to calculate SHA-256 hash
	std::string sha256(const std::string& input);
}

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

bool open_file_dialog(std::filesystem::path& path, bool debug, std::vector <std::string> filter, const std::string& current_path);
bool save_file_dialog(std::filesystem::path& path, bool debug, const svec& endings, const std::string& filename_given = "", const std::string& current_path = "");
void select_cubes(std::vector<std::vector<unsigned int>>& selection, std::vector<WFN>& wavy, unsigned int nr_of_cubes = 1, bool wfnonly = false, bool debug = false);
bool unsaved_files(std::vector<WFN>& wavy);

std::string trim(const std::string& s);

/**
 * @brief Physical memory this process can get in bytes, capped by a scheduler or container limit
 * (job object, cgroup); 0 when unknown, and the caller keeps its default.
 */
size_t available_memory_bytes();

/** std::getline that drops a CRLF file's trailing \r, which would otherwise spoil comparisons and the last field; every reader uses it. */
inline std::istream& getline_universal(std::istream& is, std::string& line)
{
	std::getline(is, line);
	if (!line.empty() && line.back() == '\r')
		line.pop_back();
	return is;
}

//Line readers for the text formats. At EOF getline leaves the line untouched, so a seek loop
//written as while (line.find(tag) == npos) getline(...) never ends on a truncated file
inline void read_line_or_fail(std::istream& is, std::string& line, const std::string& what, std::ostream& log)
{
	err_checkf(static_cast<bool>(getline_universal(is, line)), "File ends while reading " + what, log);
}
//Advances to the first line containing tag, which is left in line
inline void seek_line(std::istream& is, std::string& line, const std::string& tag, std::ostream& log)
{
	while (line.find(tag) == std::string::npos)
		read_line_or_fail(is, line, tag, log);
}
//Appends every number on line to out; a token that is not a number fails naming the block
template <typename T>
void append_numbers(const std::string& line, std::vector<T>& out, const std::string& what, std::ostream& log)
{
	std::istringstream is(line);
	T v;
	while (is >> v)
		out.push_back(v);
	err_checkf(is.eof(), "Not a number in " + what + ": '" + line + "'", log);
}

//Restores a stream's flags, precision and width on scope exit, including unwinding. std::fixed and
//setprecision are sticky, so every analysis entry point that formats output holds one of these.
struct ostream_format_guard
{
	std::ostream& stream;
	const std::ios_base::fmtflags flags;
	const std::streamsize precision;
	const std::streamsize width;
	explicit ostream_format_guard(std::ostream& s)
		: stream(s), flags(s.flags()), precision(s.precision()), width(s.width()) {}
	ostream_format_guard(const ostream_format_guard&) = delete;
	ostream_format_guard& operator=(const ostream_format_guard&) = delete;
	~ostream_format_guard()
	{
		stream.flags(flags);
		stream.precision(precision);
		stream.width(width);
	}
};

inline void print_centered_text(const std::string& text, int& bar_width, std::ostream& file = std::cout)
{
	const int text_length = static_cast<int>(text.length());
	const int total_padding = bar_width - text_length;
	const int padding_left = total_padding / 2;
	const int padding_right = (total_padding - padding_left) - 1;

	file << "["
		<< std::setw(padding_left) << std::setfill(' ') << ""
		<< text
		<< std::setw(padding_right) << std::setfill(' ') << ""
		<< "]" << std::endl;
}

inline void print_centered_message(const std::string& text, int bar_width, std::ostream& os = std::cout)
{
	const int text_length = static_cast<int>(text.length());
	const int total_padding = bar_width - text_length;
	const int padding_left = total_padding / 2;
	const int padding_right = (total_padding - padding_left) - 1;

	os
		<< std::setw(padding_left) << std::setfill(' ') << ""
		<< text
		<< std::setw(padding_right) << std::setfill(' ') << ""
		<< std::endl;
}

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

//-------------------------Progress_bar--------------------------------------------------
// LMS: My implementation of a progress bar, I would like it to stay within one line that is compatible with parallel loops
class ProgressBar
{
public:
	~ProgressBar();

	ProgressBar(const unsigned long long& worksize, const int& bar_width = 60, const std::string& fill = "#", const std::string& remainder = " ", const std::string& status_text = "", std::ostream& stream_ = std::cout)
		: worksize_(worksize), bar_width_(bar_width), fill_(fill), remainder_(remainder), status_text_(status_text), workdone(0), progress_(0.0f), workpart_(100.0f / worksize), percent_((worksize / 100 > 1) ? worksize / 100 : 1), stream_(stream_)
	{
		int bw = bar_width_ + 2;
		print_centered_text(status_text_, bw, stream_);
		linestart = stream_.tellp();
		barend_ = linestart;
#ifdef _WIN32
			initialize_taskbar_progress();
#endif
	}

	void set_progress()
	{
		progress_ = (float)workdone * workpart_;
	}

	//Called once per reflection; workdone is atomic, so only the call that crosses a reporting boundary takes the lock instead of serialising every worker on the redraw
	void update(const unsigned long long n = 1)
	{
		update_calls_.fetch_add(1, std::memory_order_relaxed);
		const unsigned long long before = workdone.fetch_add(n);
		const unsigned long long after = before + n;
		if (before / percent_ != after / percent_)
		{
#pragma omp critical
			{
				bar_writes_.fetch_add(1, std::memory_order_relaxed);
				set_progress();
				write_progress();
			}
		}
	}

	// How often callers asked, and how often that actually needed the lock.
	unsigned long long update_calls() const { return update_calls_.load(); }
	unsigned long long bar_writes() const { return bar_writes_.load(); }
	static bool report_counts;

	void write_progress();

private:
	std::ostream& stream_;
	const unsigned long long worksize_;
	const float workpart_;
	const unsigned long long percent_;
	int bar_width_;
	std::string fill_;
	std::string remainder_;
	std::string status_text_;
	std::atomic<unsigned long long> workdone;
	std::atomic<unsigned long long> update_calls_{0};
	std::atomic<unsigned long long> bar_writes_{0};
	float progress_;
	std::streampos linestart;
	//end of the bar's last write; a put position elsewhere means the loop printed in between
	std::streampos barend_{};
	bool finished_ = false;
#ifdef _WIN32
	//Assigned only inside initialize_taskbar_progress()'s SUCCEEDED checks; without the initialiser the destructor calls through stack garbage when COM refuses
	ITaskbarList3* taskbarList_ = nullptr;

	void initialize_taskbar_progress();
#endif
};

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

bool read_fracs_ADPs_from_CIF(const std::filesystem::path& cif, WFN& wavy, cell& unit_cell, std::ofstream& log3, const bool& debug);

bool read_fracs_ADPs_from_CIF(const std::filesystem::path& cif, WFN& wavy, std::ofstream& log3, const bool& debug, const bool& grown, const ivec3& symmetry_linking_list);

vec read_U_iso_from_CIF(const std::filesystem::path& cif, WFN& wavy, cell& unit_cell, std::ofstream& log3, const bool& debug);

double double_from_string_with_esd(std::string in);

void swap_sort(ivec order, cvec& v);

void swap_sort_multi(ivec order, std::vector<ivec>& v);

// Given a 3x3 symmetric matrix in a single row-major array of double, returns the median (middle) eigenvalue
double get_lambda_1(double* a);

double get_decimal_precision_from_CIF_number(std::string& given_string);

double bessel_first_kind(int l, double x);

template <typename numtype = int>
struct hashFunction
{
	size_t operator()(const std::vector<numtype>& myVector) const
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
	bool operator()(const std::vector<numtype>& vec1, const std::vector<numtype>& vec2) const
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
	bool operator()(const std::vector<numtype>& vec1, const std::vector<numtype>& vec2) const
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

struct SimplePrimitive {
	int center;     // Center index
	int type;       // Primitive type (0=s, 1=p, 2=d, etc.)
	double exp;     // Exponent
	double coefficient; // Coefficient
	int shell;
};

class primitive
{
private:
	int center, type;
	double exp, coefficient;
	double norm_const = -10;
	double exp_l_plus_3_2 = -10;
	double normalized_coefficient = -10;

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
	primitive(const SimplePrimitive& p);
	bool operator==(const primitive& other) const
	{
		return center == other.center &&
			type == other.type &&
			exp == other.exp &&
			coefficient == other.coefficient &&
			exp_l_plus_3_2 == other.exp_l_plus_3_2;
	};
	int get_center() const
	{
		return center;
	};
	int get_type() const
	{
		return type;
	};
	double get_exp_l_plus_3_2() const
	{
		return exp_l_plus_3_2;
	};
	double get_normalized_coefficient() const
	{
		return normalized_coefficient;
	};
	double get_exp() const
	{
		return exp;
	};
	double get_coef() const
	{
		return coefficient;
	};
	void set_center(const int& c)
	{
		center = c;
	};
	void set_type(const int& t)
	{
		type = t;
	};
	void set_exp(const double& e)
	{
		exp = e;
	};
	void set_coef(const double& c)
	{
		coefficient = c;
	};
	void set_norm_const(const double& nc)
	{
		norm_const = nc;
		normalized_coefficient = nc * coefficient;
	};
	double eval_gaussian(const double& r) const
	{
		return pow(r, type) * std::exp(-exp * r * r) * normalized_coefficient;
	};
	double eval_gaussian_unnormalized(const double& r) const
	{
		return pow(r, type) * std::exp(-exp * r * r) * coefficient;
	};
	inline double eval_gaussian_unnormalized(const double& rl, const double& r2) const
	{
		return rl * std::exp(-exp * r2) * coefficient;
	};
};

struct ECP_primitive : primitive
{
	int n;
	ECP_primitive() : primitive(), n(0) {}
	ECP_primitive(int c, int t, double e, double coef, int n) : primitive(c, t, e, coef), n(n) {}
};

//---------------- Object for handling all input options -------------------------------
struct options
	/** @brief All command line options and settings controlling a run. */
{
	std::ostream &log_file;
	double d_sfac_scan = 0.0;
	d3 sfac_diffuse = { 0.0, 0.0, 0.0 };
	double dmin = 99.0;
	double mem = 1000.0; // In MB
	//Set only when -mem was passed; only then is mem a budget the tsc block size and XCW I tensor window are fitted to
	bool mem_given = false;
	double efield = 0.005;
	ivec2 groups;
	ivec2 hkl_min_max{ {-100, 100}, {-100, 100}, {-100, 100} };
	vec2 twin_law;
	ivec2 combined_tsc_groups;
	pathvec combined_tsc_calc_files;
	pathvec combined_tsc_calc_cifs;
	std::vector<unsigned int> combined_tsc_calc_mult;
	ivec combined_tsc_calc_charge;
	ivec combined_tsc_calc_ECP;
	//Bounds-checked: the digesters read arguments[i + n] freely, so a flag that is last on the
	//line throws missing_argument (caught in digest_options with the flag's name) instead of
	//reading past the end
	struct missing_argument : std::out_of_range
	{
		using std::out_of_range::out_of_range;
	};
	struct checked_svec : svec
	{
		using svec::svec;
		std::string &operator[](size_t i) { return const_cast<std::string &>(std::as_const(*this)[i]); }
		const std::string &operator[](size_t i) const
		{
			if (i >= size())
				throw missing_argument("argument " + std::to_string(i));
			return svec::operator[](i);
		}
	};
	checked_svec arguments;
	pathvec combine_mo;
	svec Cations;
	svec Anions;
	pathvec pol_wfns;
	ivec cmo1;
	ivec cmo2;
	ivec ignore;
	std::filesystem::path salted_model_dir;
	//Every model given to -SALTED, in that order; salted_model_dir is the first of them
	pathvec salted_model_dirs;
	std::vector<std::shared_ptr<BasisSet>> aux_basis;
	std::filesystem::path wfn;
	std::filesystem::path wfn2;
	std::filesystem::path cube_density;
	std::filesystem::path fchk;
	std::string basis_set;
	std::filesystem::path hkl;
	std::filesystem::path cif;
	std::string method;
	std::filesystem::path xyz_file;
	std::filesystem::path coef_file;
	std::filesystem::path hirshfeld_surface;
	std::filesystem::path hirshfeld_surface2;
	std::filesystem::path fract_name;
	std::filesystem::path wavename;
	std::filesystem::path gaussian_path;
	std::filesystem::path turbomole_path;
	std::filesystem::path basis_set_path;
	std::filesystem::path anom_disp_path;
	std::string occ;
	std::filesystem::path occ_toml_path;
	std::filesystem::path cwd;
	std::filesystem::path profiling_tests_root = "tests";
	pathvec promol_nci_xyz; //Two or more fragments; a grid point is intermolecular when no single fragment dominates
	//Geometry-aid jobs (-calc_featomic_descriptor(s), -classify_atoms(_list)): the flags queue, run_app_impl runs geometry_aid::run and quits
	bool calc_featomic_descriptor = false;
	std::filesystem::path classify_atoms_out;
	std::filesystem::path geometry_aid_model;
	pathvec featomic_structures;
	pathvec classify_structures;
	double geometry_aid_cutoff = 3.5;
	bool geometry_aid_metals = false;
	double geometry_aid_center_weight = 1.0;
	std::filesystem::path interaction_energies_job;
	std::filesystem::path xcw_settings_path;
	properties_options properties;
	bool debug = false;
	//Set by the one-shot options (-merge, -dipole_moments, -convert_to_47, ...) that used to
	//exit(0) inside the parser: that killed the host when run_app is a library call (Olex2, the
	//in-process tests), and no option after them was read. run_app_impl returns 0 instead.
	bool finished = false;
	bool all_charges = false;
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
	bool label_tsc_output = false;
	bool write_CIF = false;
	bool test = false;
	bool electron_diffraction = false;
	bool ECP = false;
	bool RI_FIT = false;
	bool needs_Thakkar_fill = false;
	//Set around a spherical fill, where a disorder part already covered by an earlier one legitimately yields no atoms and must not read as a broken CIF
	bool allow_empty_asym = false;
	//Set while a spherical fill runs for somebody else: it must RETURN a block, not stream experimental.tscb out from under the table being built
	bool spherical_fill = false;
	// Per-atom EEQ charges for atoms handed to the spherical fill, as
	// {x, y, z, q} in the wavefunction's own coordinate units. Keyed by
	// POSITION rather than index because the fill rebuilds its wavefunction
	// from the original file, and that is how CIF and WFN atoms are matched
	// everywhere else here. Empty means "no charges known" -> neutral fill,
	// which is the previous behaviour.
	std::vector<std::array<double, 4>> spherical_fill_charges{};
	//Reflections per block when streaming the tsc; 0 restores the single allocation of the whole scatterers x reflections x 16 byte table
	//The exact block size is not performance-critical
	size_t tsc_block_size = 1000;
	//Set only when -tsc_block was passed, so an explicit block size wins over one derived from -mem
	bool tsc_block_given = false;
	//Reflections to hold at once for n_scat scatterers, 0 for the whole table; -tsc_block wins, then -mem, then the default
	//A block costs about three copies of n_scat * block * 16 bytes - producer, queue and writer - which is what the budget is spent against
	size_t tsc_block_for(const size_t n_refl, const size_t n_scat) const
	{
		if (tsc_block_given || !mem_given || mem <= 0.0)
			return tsc_block_size;
		const size_t item = 3 * (n_scat ? n_scat : 1) * sizeof(std::complex<double>);
		return items_within_budget(n_refl, item, static_cast<size_t>(mem * 1024.0 * 1024.0));
	}

	//set once a streamed run wrote the file itself, so the caller does not overwrite it with an empty one-shot block
	bool tsc_written_by_stream = false;
	bool qct = false;
	bool do_XCW = false;
	bool xcw_gaussian_halt = false;
	double xcw_strong_cutoff = 3.0;
	bool calc_F_calc = false;
	bool rgbi = false;
	//-npa: NAO/NPA, run in-process
	bool npa = false;
	//per-NAO occupancy table beside the NPA; -npa_summary turns it off
	bool npa_orbitals = true;
	//-ibo: IAO charges and intrinsic bond orbitals, run in-process
	bool ibo = false;
	bool rgbi_no_sym = false;
	bool rgbi_EVs = false;
	bool rgbi_theta = false;
	//-rgbi_legacy_cutoff: atomic subspace by occupation threshold (1/6 NAO, 1/14 ANO) instead of the
	//free-atom orbital count; the projector rank, and so every bond index, jumps when an occupation
	//crosses it. Only for reproducing older numbers.
	bool rgbi_legacy_cutoff = false;
	RGBIOrbitalBasis rgbi_orbital_basis = RGBIOrbitalBasis::ANO;
	ivec3 rgbi_group_sets;
	bool fract = false;
	//GPU scattering factors when a device is present; -no_gpu forces the CPU loop
	bool use_gpu = true;
	//-gpu_fp64 keeps the double sincos on a card that would otherwise pick the fp32 one
	bool gpu_fp64 = false;
	//Single-precision tiles for the CPU I tensor, as the device path runs; sgemm is twice
	//dgemm's rate. -no_cpu_itensor_fp32 keeps double.
	bool cpu_itensor_fp32 = true;
	//-itensor_hybrid lets the CPU threads take reflections alongside the device. Off by
	//default: the two sides differ in the last bits and a shared counter decides which rows
	//each takes, so the result changes from run to run. For a card slow in double only.
	bool itensor_hybrid = false;
	//Seed each lambda step from the density extrapolated through the two previous steps
	//rather than the last one alone; the step is small and the trajectory smooth.
	bool xcw_extrapolate = true;
	//Build the two-electron Fock matrix from the change of the density between iterations,
	//rather than from the whole density every time: the stored integrals skip the segments
	//the difference cannot reach (stored_eri::JK), the direct kernel skips shell quartets.
	//The direct build only with xcw_int_precision 1e-12: at 1e-10 the increments' screening
	//error accumulates to a gradient floor of 3e-5 and the SCF never meets its 1e-5, and the
	//full build at 1e-10 is the faster of the two anyway.
	bool xcw_incremental = false;
	//Integral screening threshold of the XCW Fock build; OCC's own default is 1e-12.
	double xcw_int_precision = 1e-10;
	//-gflops reports achieved GFLOP/s per stage for the CPU and GPU paths at the end of a
	//run. The thresholds deciding what goes to the device were calibrated on one machine;
	//this is how they get re-derived on another.
	bool track_gflops = false;
	//-gpu_fp32 forces the reduced-argument single-precision sincos on a card that would
	//otherwise keep the double one. It exists for the test suite: without it the precision
	//a run uses depends on the card, so neither path can be pinned. -gpu_fp64 wins if both
	//are given, the accurate path being the safer thing to fall back to.
	bool gpu_fp32 = false;
	//-no_gpu_itensor keeps the XCW I tensor GEMMs on the CPU. On by default: it is much the
	//largest of the device paths, and it moves the total energy only in the tenth
	//significant figure. Read together with use_gpu, so -no_gpu turns it off as well.
	bool gpu_itensor = true;
	//FP16 Tensor Core operands with FP32 accumulation for the I tensor. Off by default: the
	//half-precision operands move the XCW energies in the fourth decimal, plain FP32 GEMM
	//sits within 1e-8 Eh of double.
	bool gpu_itensor_tensor = false;
	//SALTED descriptor combination uses the device when one is available; -no_gpu_salted keeps it on the CPU.
	bool gpu_salted = true;
	//-no_gpu_grid keeps the Becke/TFVC integration weights on the CPU
	bool gpu_grid = true;
	//Owned by the caller; the scattering-factor grid is built in it instead of a local, so
	//a second table for the same geometry reuses the points and weights
	GridManager* grid_cache = nullptr;
	//-no_gpu_density keeps the fitted density of the Gordon-Kim repulsion grid and the spherical-atom grids on the CPU
	bool gpu_density = true;
	//-gpu_blas offers large dense GEMMs in nos_math to the device
	bool gpu_blas = false;
	//The I tensor GEMM goes through cuBLAS when the machine has it, and through the
	//built-in CUTLASS path otherwise. cuBLAS is 1.65x faster on a V100 and level with
	//CUTLASS within measurement noise on consumer cards, so preferring it costs nothing
	//where it does not help. -no_gpu_cublas pins CUTLASS, which is what the reference tests
	//do: the two differ in the last digits, so a test left to pick would pass or fail on
	//whether a CUDA toolkit happened to be installed.
	bool gpu_cublas = true;
	//Standalone conceptual-DFT reactivity analysis (-fukui_analysis), run from run_app_impl rather than at parse time so its output survives
	bool fukui_analysis_run = false;
	//EQC decomposition (-eqc), run from run_app_impl like -fukui_analysis; fragments are 0-based atoms, charge, mult
	bool eqc = false;
	struct eqc_fragment
	{
		ivec atoms;
		int charge = 0;
		int mult = 1;
	};
	std::vector<eqc_fragment> eqc_frags;
	std::vector<std::filesystem::path> eqc_wfns;
	std::string eqc_method = "hf";  //occ's method name: hf or a functional
	std::string eqc_basis;          //a basis-library name; empty takes the basis out of the gbw
	bool eqc_cold = false;          //also converge everything from occ's own guess, to count what seeding saves
	//Basin analysis (-eli_analysis), run from run_app_impl for the same reason
	bool eli_analysis_run = false;
	//-fba: RGBI, native NBO/NPA with NRT, bondwise Laplacian, QTAIM and ELI-D on one wavefunction
	bool fba = false;
	bool profiling = false;
	bool promol_nci = false;
	bool get_g = false;
	int accuracy = 2;
	//-basin_grid <n>: the quadrature of the basin analysis pulled into the core, tightest
	//exponent sharpened n^2-fold, radial step divided by n, Lebedev order up n - 1 entries
	int basin_grid = 1;
	//-basin_cube: QTAIM and ELI-D basins on the cube instead of from the analytic critical points and field
	bool basin_cube = false;
	//-no_spin_eli: skip the spin-resolved ELI-D basins of a spin-polarised wavefunction
	bool spin_eli = true;
	int threads = -1;
	int pbc = 0;
	int charge = 0;
	int ECP_mode = 0;
	PartitionType partition_type = PartitionType::Hirshfeld;
	//-multipole_moments: the RI fit is restrained to this scheme's atomic moments up to this order, -1 = unrestrained
	int multipole_lmax = -1;
	MultipoleScheme multipole_scheme = MultipoleScheme::HIRSHFELD;
	double multipole_strength = 1.0;
	//-repulsion_overlap: exchange-repulsion of -interaction_energy as K * Int rhoA rhoB, 0 = not included
	double repulsion_overlap = 0.0;
	//-repulsion_exchange: exchange functional of the Gordon-Kim repulsion, 0 Dirac, 1 PBE, 2 B88, 3 r2SCAN-L
	int repulsion_exchange = 0;
	unsigned int mult = 0;
	hkl_list m_hkl_list;

	/** @brief Finds the debug flag, removes it from argc/argv and stores the rest internally. */
	void look_for_debug(int& argc, char** argv);
	/** @brief Digests the options; call only after look_for_debug(). */
	//file, format and conversion options
	bool digest_io_options(const std::string &temp, int &i);
	//resources, accuracy and general run control
	bool digest_run_options(const std::string &temp, int &i);
	//partitioning schemes, scattering factors and tsc tables
	bool digest_partition_options(const std::string &temp, int &i);
	//cube and property evaluation
	bool digest_property_options(const std::string &temp, int &i);
	//RI fitting, featomic descriptors and atom classification
	bool digest_ri_options(const std::string &temp, int &i);
	//X-ray constrained wavefunction fitting
	bool digest_xcw_options(const std::string &temp, int &i);
	//development and test-only switches
	bool digest_dev_options(const std::string &temp, int &i);
	void digest_options();
	/** @brief The error for a requested analysis with no wavefunction to run on, "" when runnable;
	 *  without -wfn/-occ run_app_impl skips its wavefunction branch silently. */
	std::string unrunnable_analysis() const;
	/** @brief Refuses -rgbi/-npa with an analysis that ends the run before them, and an -nbo_/-nrt_
	 *  option on a line with no NBO analysis. Called last in digest_options(), so order does not matter. */
	void refuse_unread_bonding_options();

	options() : log_file(std::cout)
	{
		groups.resize(1);
	};
	options(int& argc, char** argv, std::ostream& log) : log_file(log)
	{
		groups.resize(1);
		look_for_debug(argc, argv);
	};
};

/** @brief The analysis whose options begin with this flag's prefix, nullptr if none;
 *  digest_options refuses an unclaimed flag in such a family. */
const char *owning_analysis(const std::string &flag);

/** @brief The -nbo_/-nrt_ options the NBO handlers read from their own tokens, not through a digester.
 *  An option added to a handler's scan loop needs an entry here (CliRefusalTests asserts it). */
const std::set<std::string> &nbo_family_suboptions();

void convert_tonto_XCW_lambda_steps(const std::string& str, const std::string& lambda_step, bool debug, options& opt);

double hypergeometric(double a, double b, double c, double x);

cdouble hypergeometric(double a, double b, double c, cdouble x);

bool ends_with(const std::string& str, const std::string& suffix);

bool read_block_from_fortran_binary(std::ifstream& file, void* Target, const size_t capacity);
template <typename T>
bool read_block_from_fortran_binary(std::ifstream& file, std::vector<T>& Target);
