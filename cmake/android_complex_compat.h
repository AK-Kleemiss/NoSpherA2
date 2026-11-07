/* Force-included into C code on Android below API 26 (legacy APK line, API 21): bionic's libm gains
 * the C99 complex functions only at API 23 (creal ... csqrt, cexp, ccos, csin) and 26 (clog, cpow).
 * Used for libcint and for the API-21 OpenBLAS (its f2c LAPACK helpers call conj, cabs, cexp, ...).
 * Each name maps to a static inline built from real libm calls that exist at API 21. */
#pragma once
#include <complex.h>
#include <math.h>
#if defined(__ANDROID__) && __ANDROID_API__ < 26
static inline double _Complex nsa2_c(double r, double i) { return __builtin_complex(r, i); }
static inline double nsa2_cabs(double _Complex z) { return hypot(__builtin_creal(z), __builtin_cimag(z)); }
static inline double nsa2_carg(double _Complex z) { return atan2(__builtin_cimag(z), __builtin_creal(z)); }
static inline double _Complex nsa2_cexp(double _Complex z) {
  double e = exp(__builtin_creal(z)), y = __builtin_cimag(z);
  return nsa2_c(e * cos(y), e * sin(y));
}
static inline double _Complex nsa2_clog(double _Complex z) { return nsa2_c(log(nsa2_cabs(z)), nsa2_carg(z)); }
static inline double _Complex nsa2_csqrt(double _Complex z) { /* principal root, branch cut on the negative real axis */
  double x = __builtin_creal(z), y = __builtin_cimag(z);
  if (x == 0.0 && y == 0.0) return nsa2_c(0.0, y);
  double t = sqrt((fabs(x) + hypot(x, y)) * 0.5);
  return x >= 0.0 ? nsa2_c(t, y / (2.0 * t)) : nsa2_c(fabs(y) / (2.0 * t), copysign(t, y));
}
static inline double _Complex nsa2_cpow(double _Complex a, double _Complex b) {
  if (__builtin_creal(a) == 0.0 && __builtin_cimag(a) == 0.0)
    return (__builtin_creal(b) == 0.0 && __builtin_cimag(b) == 0.0) ? nsa2_c(1.0, 0.0) : nsa2_c(0.0, 0.0);
  return nsa2_cexp(b * nsa2_clog(a));
}
static inline double _Complex nsa2_ccos(double _Complex z) {
  double x = __builtin_creal(z), y = __builtin_cimag(z);
  return nsa2_c(cos(x) * cosh(y), -sin(x) * sinh(y));
}
static inline double _Complex nsa2_csin(double _Complex z) {
  double x = __builtin_creal(z), y = __builtin_cimag(z);
  return nsa2_c(sin(x) * cosh(y), cos(x) * sinh(y));
}
#if __ANDROID_API__ < 23
#define creal(z) __builtin_creal(z)
#define cimag(z) __builtin_cimag(z)
#define conj(z) __builtin_conj(z)
#define crealf(z) __builtin_crealf(z)
#define cimagf(z) __builtin_cimagf(z)
#define conjf(z) __builtin_conjf(z)
#define cabs(z) nsa2_cabs(z)
#define cabsf(z) ((float)nsa2_cabs(z))
#define carg(z) nsa2_carg(z)
#define cargf(z) ((float)nsa2_carg(z))
#define cexp(z) nsa2_cexp(z)
#define cexpf(z) ((float _Complex)nsa2_cexp(z))
#define csqrt(z) nsa2_csqrt(z)
#define csqrtf(z) ((float _Complex)nsa2_csqrt(z))
#define ccos(z) nsa2_ccos(z)
#define ccosf(z) ((float _Complex)nsa2_ccos(z))
#define csin(z) nsa2_csin(z)
#define csinf(z) ((float _Complex)nsa2_csin(z))
#endif
#define clog(z) nsa2_clog(z)
#define clogf(z) ((float _Complex)nsa2_clog(z))
#define cpow(a, b) nsa2_cpow(a, b)
#define cpowf(a, b) ((float _Complex)nsa2_cpow(a, b))
#endif
