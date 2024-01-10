/* Copyright (c) 2012 Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

 /* Available at: http://ab-initio.mit.edu/Faddeeva

	Header file for Faddeeva.c; see Faddeeva.cc for more information. */

	//THIS VERSION HAS BEEN MODIFIED BY F. KLEEMISS FOR USE IN NOSPHERA2

#ifndef FADDEEVA_H
#define FADDEEVA_H 1

#include <complex>

	// compute w(z) = exp(-z^2) erfc(-iz) [ Faddeeva / scaled complex error func ]
extern std::complex<double> Faddeeva_w(std::complex<double> z, double relerr);
extern double Faddeeva_w_im(double x); // special-case code for Im[w(x)] of real x

// Various functions that we can compute with the help of w(z)

// compute erfcx(z) = exp(z^2) erfc(z)
extern std::complex<double> Faddeeva_erfcx(std::complex<double> z, double relerr);
extern double Faddeeva_erfcx_re(double x); // special case for real x

// compute erf(z), the error function of complex arguments
extern std::complex<double> Faddeeva_erf(std::complex<double> z, double relerr);
extern double Faddeeva_erf_re(double x); // special case for real x

// compute erfi(z) = -i erf(iz), the imaginary error function
extern std::complex<double> Faddeeva_erfi(std::complex<double> z, double relerr);
extern double Faddeeva_erfi_re(double x); // special case for real x

// compute erfc(z) = 1 - erf(z), the complementary error function
extern std::complex<double> Faddeeva_erfc(std::complex<double> z, double relerr);
extern double Faddeeva_erfc_re(double x); // special case for real x

// compute Dawson(z) = sqrt(pi)/2  *  exp(-z^2) * erfi(z)
extern std::complex<double> Faddeeva_Dawson(std::complex<double> z, double relerr);
extern double Faddeeva_Dawson_re(double x); // special case for real x


#endif // FADDEEVA_H
