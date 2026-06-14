#include "../Src/pch.h"
#include "unit_test_api.h"
#include "../Src/convenience.h"
#include "../Src/fchk.h"
#include <cstring>
#include <cmath>

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

} // extern "C"
