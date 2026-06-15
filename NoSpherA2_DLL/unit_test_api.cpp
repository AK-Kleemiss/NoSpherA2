#include "../Src/pch.h"
#include "unit_test_api.h"
#include "../Src/convenience.h"
#include "../Src/fchk.h"
#include "../Src/constants.h"
#include <cstring>
#include <cmath>
#include <thread>
#include <chrono>

extern "C" {

    UT_API double ut_array_length3(double x, double y, double z) {
        return std::hypot(x, y, z);
    }

    UT_API double ut_array_distance3(double x1, double y1, double z1,
                                     double x2, double y2, double z2) {
        d3 a{ x1, y1, z1 }, b{ x2, y2, z2 };
        return array_length(a, b);
    }

    UT_API void ut_vec_diff(double ax, double ay, double az,
                            double bx, double by, double bz, double out[3]) {
        d3 r = vec_diff({ ax, ay, az }, { bx, by, bz });
        out[0] = r[0]; out[1] = r[1]; out[2] = r[2];
    }

    UT_API void ut_vec_cross(double ax, double ay, double az,
                             double bx, double by, double bz, double out[3]) {
        d3 r = vec_cross({ ax, ay, az }, { bx, by, bz });
        out[0] = r[0]; out[1] = r[1]; out[2] = r[2];
    }

    UT_API double ut_vec_dot(double ax, double ay, double az,
                             double bx, double by, double bz) {
        return vec_dot({ ax, ay, az }, { bx, by, bz });
    }

    UT_API int ut_is_nan(double x) {
        // is_nan() uses "x != x" which /fp:fast optimises away; use std::isnan.
        return std::isnan(x) ? 1 : 0;
    }

    UT_API int ut_is_similar_rel(double a, double b, double tol) {
        return is_similar_rel(a, b, tol) ? 1 : 0;
    }

    UT_API int ut_is_similar_abs(double a, double b, double tol) {
        return is_similar_abs(a, b, tol) ? 1 : 0;
    }

    UT_API double ut_fast_exp_neg(double x) {
        return fast_exp_neg(x);
    }

    UT_API int ut_read_fchk_int(const char* line) {
        return read_fchk_integer(std::string(line));
    }

    UT_API double ut_read_fchk_dbl(const char* line) {
        return read_fchk_double(std::string(line));
    }

    UT_API int ut_ends_with(const char* str, const char* suffix) {
        return ends_with(std::string(str), std::string(suffix)) ? 1 : 0;
    }

    UT_API int ut_shrink_string(const char* input, char* output, int bufsize) {
        std::string s(input);
        std::string r = shrink_string(s);
        if (static_cast<int>(r.size()) >= bufsize)
            return -1;
        std::memcpy(output, r.c_str(), r.size() + 1);
        return static_cast<int>(r.size());
    }

    UT_API int ut_sha256(const char* input, char* output, int bufsize) {
        std::string r = sha::sha256(std::string(input));
        if (static_cast<int>(r.size()) >= bufsize)
            return -1;
        std::memcpy(output, r.c_str(), r.size() + 1);
        return static_cast<int>(r.size());
    }

    // --- Vector aggregates ---

    UT_API int ut_vec_sum_bool(const int* arr, int n) {
        if (n <= 0) return 0;
        bvec v(n);
        for (int i = 0; i < n; ++i) v[i] = (arr[i] != 0);
        return vec_sum(v);
    }

    UT_API int ut_vec_sum_int(const int* arr, int n) {
        if (n <= 0) return 0;
        ivec v(arr, arr + n);
        return vec_sum(v);
    }

    UT_API double ut_vec_sum_double(const double* arr, int n) {
        if (n <= 0) return 0.0;
        vec v(arr, arr + n);
        return vec_sum(v);
    }

    UT_API double ut_vec_length(const double* arr, int n) {
        if (n <= 0) return 0.0;
        vec v(arr, arr + n);
        return vec_length(v);
    }

    // --- Additional string utilities ---

    UT_API int ut_trim(const char* input, char* output, int bufsize) {
        std::string r = trim(std::string(input));
        if (static_cast<int>(r.size()) >= bufsize) return -1;
        std::memcpy(output, r.c_str(), r.size() + 1);
        return static_cast<int>(r.size());
    }

    UT_API char ut_asciitolower(char c) {
        return asciitolower(c);
    }

    UT_API double ut_double_from_esd(const char* s) {
        return double_from_string_with_esd(std::string(s));
    }

    UT_API double ut_decimal_precision_cif(const char* s) {
        std::string str(s);
        return get_decimal_precision_from_CIF_number(str);
    }

    // --- Basis-type helpers ---

    UT_API int ut_sht2nbas(int type) {
        return sht2nbas(type);
    }

    UT_API unsigned int ut_doublefactorial(int n) {
        return doublefactorial(n);
    }

    // --- Unit conversions ---

    UT_API double ut_bohr2ang(double inp)       { return constants::bohr2ang(inp); }
    UT_API double ut_ang2bohr(double inp)        { return constants::ang2bohr(inp); }
    UT_API double ut_cubic_bohr2ang(double inp)  { return constants::cubic_bohr2ang(inp); }
    UT_API double ut_cubic_ang2bohr(double inp)  { return constants::cubic_ang2bohr(inp); }

    // --- Constexpr math ---

    UT_API int    ut_const_abs(int x)              { return constants::const_abs(x); }
    UT_API double ut_constexpr_pow(double b, int e){ return constants::constexpr_pow(b, e); }
    UT_API double ut_constants_sqrt(double x)      { return constants::sqrt(x); }
    UT_API double ut_exp_approx(double x, int n)   { return constants::exp_approx(x, n); }
    UT_API double ut_log_approx(double x, int n)   { return constants::log_approx(x, n); }
    UT_API long long ut_factorial(int n)           { return static_cast<long long>(constants::ft_fun(n)); }

    // --- Orbital index mappings ---

    UT_API int ut_orca_2_pyscf(unsigned int l, int m_idx) {
        auto r = constants::orca_2_pySCF(l, m_idx);
        return r.has_value() ? static_cast<int>(r.value()) : -1;
    }

    UT_API unsigned int ut_type_2_nbo(int type) {
        return constants::type_2_nbo(type);
    }

    // --- Spherical Bessel functions ---

    UT_API double ut_bessel_j(int l, double x) {
        return bessel_first_kind(l, x);
    }

    // --- is_similar (pow-10 tolerance) ---

    UT_API int ut_is_similar_pow10(double a, double b, double tolerance) {
        return is_similar(a, b, tolerance) ? 1 : 0;
    }

    // --- shell2function ---

    UT_API int ut_shell2function(int type, int prim) {
        return shell2function(type, prim);
    }

    // --- CountWords ---

    UT_API int ut_count_words(const char* str) {
        return CountWords(str);
    }

    // --- shrink_string_to_atom ---

    UT_API int ut_shrink_string_to_atom(const char* input, int atom_number,
                                         char* output, int bufsize) {
        std::string s(input);
        std::string r = shrink_string_to_atom(s, atom_number);
        if (static_cast<int>(r.size()) >= bufsize)
            return -1;
        std::memcpy(output, r.c_str(), r.size() + 1);
        return static_cast<int>(r.size());
    }

    // --- split_string ---
    // Splits input by delim, writes up to max_out tokens into out_tokens[].
    // Returns total token count (may exceed max_out).

    UT_API int ut_split_string(const char* input, const char* delim,
                                char** out_tokens, int max_out, int tok_bufsize) {
        svec tokens = split_string<std::string>(std::string(input), std::string(delim));
        int written = std::min(static_cast<int>(tokens.size()), max_out);
        for (int i = 0; i < written; ++i) {
            strncpy_s(out_tokens[i], tok_bufsize, tokens[i].c_str(), tok_bufsize - 1);
        }
        return static_cast<int>(tokens.size());
    }

    // --- Timing helper ---

    UT_API long long ut_sleep_and_measure_us(int sleep_ms) {
        auto t0 = get_time();
        std::this_thread::sleep_for(std::chrono::milliseconds(sleep_ms));
        auto t1 = get_time();
        return get_musec(t0, t1);
    }

    // --- hypergeometric ---

    UT_API double ut_hypergeometric(double a, double b, double c, double x) {
        return hypergeometric(a, b, c, x);
    }

    // --- Atomic number → element symbol ---

    UT_API int ut_atnr2letter(int nr, char* output, int bufsize) {
        const char* sym = constants::atnr2letter(nr);
        int len = static_cast<int>(std::strlen(sym));
        if (len >= bufsize) return -1;
        std::memcpy(output, sym, len + 1);
        return len;
    }

    // --- Cartesian exponent vector from basis type ---

    UT_API void ut_type2vector(int index, int nx_out[1], int ny_out[1], int nz_out[1]) {
        int v[3];
        constants::type2vector(index, v);
        nx_out[0] = v[0];
        ny_out[0] = v[1];
        nz_out[0] = v[2];
    }

    // --- Gaussian normalization constant ---

    UT_API double ut_normgauss(int type, double exponent) {
        return constants::normgauss(type, exponent);
    }

    // --- Associated Legendre polynomial P_l^m(x) ---

    UT_API double ut_assoc_legendre(int l, int m, double x) {
        return constants::associated_legendre_polynomial(l, m, x);
    }

    // --- Cartesian → spherical coordinates ---

    UT_API void ut_cartesian_to_spherical(double x, double y, double z, double out3[3]) {
        std::vector<double> r = constants::cartesian_to_spherical(x, y, z);
        out3[0] = r[0];
        out3[1] = r[1];
        out3[2] = r[2];
    }

    // --- Median eigenvalue of 3x3 symmetric matrix ---

    UT_API double ut_get_lambda_1(const double a[9]) {
        // get_lambda_1 takes double* and reads a[0..8]
        return get_lambda_1(const_cast<double*>(a));
    }

} // extern "C"
