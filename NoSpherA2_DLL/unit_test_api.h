#pragma once
// Lightweight export header for VS unit tests.
// Uses only primitive types so UnitTests.cpp needs no Src/ include paths.

#ifdef BUILDING_DLL
#define UT_API __declspec(dllexport)
#else
#define UT_API __declspec(dllimport)
#endif

extern "C" {
    // --- Geometry ---
    UT_API double ut_array_length3(double x, double y, double z);
    UT_API double ut_array_distance3(double x1, double y1, double z1,
                                     double x2, double y2, double z2);
    UT_API void   ut_vec_diff (double ax, double ay, double az,
                               double bx, double by, double bz, double out[3]);
    UT_API void   ut_vec_cross(double ax, double ay, double az,
                               double bx, double by, double bz, double out[3]);
    UT_API double ut_vec_dot  (double ax, double ay, double az,
                               double bx, double by, double bz);

    // --- Numeric checks ---
    UT_API int    ut_is_nan        (double x);
    UT_API int    ut_is_similar_rel(double a, double b, double tol);
    UT_API int    ut_is_similar_abs(double a, double b, double tol);
    UT_API double ut_fast_exp_neg  (double x);

    // --- FCHK line parsing ---
    // Input: a full FCHK header line; value begins at character 49.
    UT_API int    ut_read_fchk_int(const char* line);
    UT_API double ut_read_fchk_dbl(const char* line);

    // --- String utilities ---
    UT_API int    ut_ends_with    (const char* str, const char* suffix);
    // Strips digits and spaces from input; writes result to output (bufsize incl. NUL).
    // Returns result length, or -1 if buffer too small.
    UT_API int    ut_shrink_string(const char* input, char* output, int bufsize);

    // --- SHA-256 ---
    // Writes the 64-char lowercase hex digest + NUL into output (need bufsize >= 65).
    // Returns 64 on success, -1 if buffer too small.
    UT_API int    ut_sha256(const char* input, char* output, int bufsize);

    // --- Vector aggregates ---
    // Pass elements as an array + count. bools encoded as 0/1 ints.
    UT_API int    ut_vec_sum_bool  (const int* arr, int n);
    UT_API int    ut_vec_sum_int   (const int* arr, int n);
    UT_API double ut_vec_sum_double(const double* arr, int n);
    UT_API double ut_vec_length    (const double* arr, int n);

    // --- Additional string utilities ---
    UT_API int    ut_trim                   (const char* input, char* output, int bufsize);
    UT_API char   ut_asciitolower           (char c);
    UT_API double ut_double_from_esd        (const char* s);   // strips "(esd)"
    UT_API double ut_decimal_precision_cif  (const char* s);   // uncertainty from "1.234(5)"

    // --- Basis-type helpers ---
    UT_API int          ut_sht2nbas        (int type);     // shell type → Cartesian/spherical count
    UT_API unsigned int ut_doublefactorial (int n);        // n!! (double factorial)

    // --- Unit conversions (constants namespace) ---
    UT_API double ut_bohr2ang       (double inp);
    UT_API double ut_ang2bohr       (double inp);
    UT_API double ut_cubic_bohr2ang (double inp);
    UT_API double ut_cubic_ang2bohr (double inp);

    // --- Constexpr math (constants namespace) ---
    UT_API int    ut_const_abs     (int x);
    UT_API double ut_constexpr_pow (double base, int exponent);
    UT_API double ut_constants_sqrt(double x);
    UT_API double ut_exp_approx    (double x, int n);
    UT_API double ut_log_approx    (double x, int n);
    UT_API long long ut_factorial  (int n);              // n! from lookup table

    // --- Orbital index mappings (constants namespace) ---
    // orca_2_pySCF: returns index, or -1 when l is out of range.
    UT_API int          ut_orca_2_pyscf (unsigned int l, int m_idx);
    UT_API unsigned int ut_type_2_nbo   (int type);

    // --- Spherical Bessel functions ---
    UT_API double ut_bessel_j(int l, double x);
}
