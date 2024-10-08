/* Library libcerf:
 *   compute complex error functions,
 *   along with Dawson, Faddeeva and Voigt functions
 *
 * File defs.h:
 *   Define CMPLX, NaN, for internal use, for when sources are compiled as C++ code
 *
 * Copyright:
 *   (C) 2012 Massachusetts Institute of Technology
 *   (C) 2013 Forschungszentrum Jülich GmbH
 *
 * Licence:
 *   MIT Licence.
 *   See ../COPYING
 *
 * Authors:
 *   Steven G. Johnson, Massachusetts Institute of Technology, 2012, core author
 *   Joachim Wuttke, Forschungszentrum Jülich, 2013, package maintainer
 *
 * Website:
 *   http://apps.jcns.fz-juelich.de/libcerf
 */
#pragma once
#include <cfloat>
#include <cmath>
#include <limits>
#include <complex>

// use std::numeric_limits, since 1./0. and 0./0. fail with some compilers (MS)
#define Inf std::numeric_limits<double>::infinity()
#define NaN std::numeric_limits<double>::quiet_NaN()

// Use C-like complex syntax, since the C syntax is more restrictive
#define cexp(z) std::exp(z)
#define creal(z) z.real()
#define cimag(z) z.imag()
#define cpolar(r, t) std::polar(r, t)

#define _cdouble std::complex<double>

double voigt_hwhm(double sigma, double gamma);
_cdouble w_of_z(_cdouble z);
double im_w_of_x(double x);
double re_w_of_z(double x, double y);
double erfcx(double x);
_cdouble cerf(_cdouble z);
double voigt(double x, double sigma, double gamma);

/******************************************************************************/
/*  auxiliary functions                                                       */
/******************************************************************************/

constexpr double sinc(double x, double sinx)
{
    // return sinc(x) = sin(x)/x, given both x and sin(x)
    // [since we only use this in cases where sin(x) has already been computed]
    return fabs(x) < 1e-4 ? 1 - (0.1666666666666666666667) * x * x : sinx / x;
}

constexpr double sinh_taylor(double x)
{
    // sinh(x) via Taylor series, accurate to machine precision for |x| < 1e-2
    return x * (1 + (x * x) * (0.1666666666666666666667 + 0.00833333333333333333333 * (x * x)));
}

constexpr double sqr(double x) { return x * x; }

// If we are using the gnulib <cmath> (e.g. in the GNU Octave sources),
// gnulib generates a link warning if we use ::floor instead of gnulib::floor.
// This warning is completely innocuous because the only difference between
// gnulib::floor and the system ::floor (and only on ancient OSF systems)
// has to do with floor(-0), which doesn't occur in the usage below, but
// the Octave developers prefer that we silence the warning.
#ifdef GNULIB_NAMESPACE
#define floor GNULIB_NAMESPACE::floor
#endif
