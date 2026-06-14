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
}
