/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <cmath>
#include <iostream>

#include "grid_radial.h"
# define M_PI           3.14159265358979323846  /* pi */
# define C0             4.0 * std::sqrt(2.0) * M_PI  /* C0 */
# define M_PI2          pow(M_PI,2)  /* pi^2 */
# define TG32           tgamma(3.0 / 2.0)
#include "parameters.h"

// TCA 106, 178 (2001), eq. 25
// we evaluate r_inner for s functions
double get_r_inner(const double max_error, const double alpha_inner)
{
    int m = 0;
    double d = 1.9;

    double r = d - std::log(1.0 / max_error);
    r = r * 2.0 / (m + 3.0);
    r = std::exp(r) / (alpha_inner);
    r = std::sqrt(r);

    return r;
}

// TCA 106, 178 (2001), eq. 19
double get_r_outer(const double max_error,
                   const double alpha_outer,
                   const int l,
                   const double guess)
{
    const int m = 2 * l;
    double r = guess;
    double r_old = 1.0e50;
    double c, a, e;
    double step = 0.5;
    double sign, sign_old;
    double f = 1.0e50;

    (f > max_error) ? (sign = 1.0) : (sign = -1.0);

    while (std::abs(r_old - r) > SMALL)
    {
        c = tgamma((m + 3.0) / 2.0);
        a = std::pow(alpha_outer * r * r, (m + 1.0) / 2.0);
        e = std::exp(-alpha_outer * r * r);
        f = c * a * e;

        sign_old = sign;
        (f > max_error) ? (sign = 1.0) : (sign = -1.0);
        if (r < 0.0)
            sign = 1.0;
        if (sign != sign_old)
            step *= 0.1;

        r_old = r;
        r += sign * step;
    }

    return r;
}

// TCA 106, 178 (2001), eqs. 17 and 18
double get_h(const double max_error, const int l, const double guess)
{
    const int m = 2 * l;
    double h = guess;
    double h_old = h*2;
    double step = 0.1 * guess;
    double sign, sign_old;
    double f = 1.0e50;
    double cm, pm, rd0, e0;

    (f > max_error) ? (sign = -1.0) : (sign = 1.0);

    while (std::abs(h_old - h) > SMALL)
    {
        cm = TG32 / tgamma((m + 3.0) / 2.0);
        e0 = std::exp(-M_PI2 / (2.0 * h));
        pm = std::pow(M_PI / h, m / 2.0);
        rd0 = C0 / h * e0;
        f = cm * pm * rd0;

        sign_old = sign;
        (f > max_error) ? (sign = -1.0) : (sign = 1.0);
        if (h < 0.0)
            sign = 1.0;
        if (sign != sign_old)
            step *= 0.1;

        h_old = h;
        h += sign * step;
        if (h < 0.007) h = 0.007;
    }

    return h;
}
