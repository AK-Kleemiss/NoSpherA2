/*
 * File width.c:
 *   Compute voit_hwhm, half width at half maximum of the Voigt profile,
 *   using iterative regula falsi (Illinois method).
 *   Modified from code originally written for gnuplot.
 *
 * Copyright:
 *   Ethan A Merritt 2020, 2021; Joachim Wuttke 2021
 *
 * License:
 *   Permission is hereby granted, free of charge, to any person obtaining
 *   a copy of this software and associated documentation files (the
 *   "Software"), to deal in the Software without restriction, including
 *   without limitation the rights to use, copy, modify, merge, publish,
 *   distribute, sublicense, and/or sell copies of the Software, and to
 *   permit persons to whom the Software is furnished to do so, subject to
 *   the following conditions:
 *
 *   The above copyright notice and this permission notice shall be
 *   included in all copies or substantial portions of the Software.
 *
 *   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 *   LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *   OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 *   WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */

#include <cassert>
#include <cmath>
using std::isnan;
#include "cerf.h"

#ifndef DBL_EPSILON
#define DBL_EPSILON 2.2204460492503131E-16
#endif

/* This approximation claims accuracy of 0.02%
 * Olivero & Longbothum [1977]
 * Journal of Quantitative Spectroscopy and Radiative Transfer. 17:233
 */
static double hwhm0(double sigma, double gamma)
{
    return .5 * (1.06868 * gamma + sqrt(0.86743 * gamma * gamma + 4 * 2 * log(2) * sigma * sigma));
}

double voigt_hwhm(double sigma, double gamma)
{
    double HM;
    double a, b, c; /* 3 points used by regula falsi */
    double del_a, del_b, del_c;
    int k;
    int side = 0;

    if (sigma == 0 && gamma == 0)
        return 0;
    if (isnan(sigma) || isnan(gamma))
        return NAN;

    /* Reduce order of magnitude to prevent overflow */
    double prefac = 1.;
    double s = fabs(sigma);
    double g = fabs(gamma);
    while (s > 0x1p+320 || g > 0x1p+320)
    {
        prefac *= 0x1p+80;
        s *= 0x1p-80;
        g *= 0x1p-80;
    }

    /* Increase order of magnitude to prevent underflow */
    while (s < 0x1p-320 && g < 0x1p-320)
    {
        prefac *= 0x1p-80;
        s *= 0x1p+80;
        g *= 0x1p+80;
    }

    HM = voigt(0.0, s, g) / 2;

    /* Choose initial points a,b that bracket the expected root */
    c = hwhm0(s, g);
    a = c * 0.995;
    b = c * 1.005;
    del_a = voigt(a, s, g) - HM;
    del_b = voigt(b, s, g) - HM;

    /* Iteration using regula falsi (Illinois variant).
     * Empirically, this takes <5 iterations to converge to FLT_EPSILON
     * and <10 iterations to converge to DBL_EPSILON.
     * We have never seen convergence worse than k = 15.
     */
    for (k = 0; k < 30; k++)
    {
        if (fabs(del_a - del_b) < 2 * DBL_EPSILON * HM)
            return prefac * (a + b) / 2;
        c = (b * del_a - a * del_b) / (del_a - del_b);
        if (fabs(b - a) < 2 * DBL_EPSILON * fabs(b + a))
            return prefac * c;
        del_c = voigt(c, s, g) - HM;

        if (del_b * del_c > 0)
        {
            b = c;
            del_b = del_c;
            if (side < 0)
                del_a /= 2;
            side = -1;
        }
        else if (del_a * del_c > 0)
        {
            a = c;
            del_a = del_c;
            if (side > 0)
                del_b /= 2;
            side = 1;
        }
        else
        {
            return prefac * c;
        }
    }
    assert(0); /* One should never arrive here */
}
