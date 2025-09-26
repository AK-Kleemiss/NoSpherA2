#include "pch.h"
#include "spherical_density.h"
#include "convenience.h"
#include "Thakkar_coefs.h"
#include "def2-ECPs_GA.h"
#include "ECPs_corrections.h"
#include "constants.h"

Thakkar::Thakkar(const int g_atom_number, const int ECP_m) : Spherical_Atom(g_atom_number, ECP_m)
{
    nex = &(Thakkar_nex[0]);
    _first_ex = first_ex();
    ns = &(Thakkar_ns[0]);
    np = &(Thakkar_np[0]);
    nd = &(Thakkar_nd[0]);
    nf = &(Thakkar_nf[0]);
    occ = &(Thakkar_occ[0]);
    n = &(Thakkar_n[0]);
    z = &(Thakkar_z[0]);
    c = &(Thakkar_c[0]);
    if (atomic_number == 1)
        _prev_coef = 0;
    else
        _prev_coef = previous_element_coef();
};
Thakkar::Thakkar() : Spherical_Atom()
{
    nex = &(Thakkar_nex[0]);
    _first_ex = first_ex();
    ns = &(Thakkar_ns[0]);
    np = &(Thakkar_np[0]);
    nd = &(Thakkar_nd[0]);
    nf = &(Thakkar_nf[0]);
    occ = &(Thakkar_occ[0]);
    n = &(Thakkar_n[0]);
    z = &(Thakkar_z[0]);
    c = &(Thakkar_c[0]);
    if (atomic_number == 1)
        _prev_coef = 0;
    else
        _prev_coef = previous_element_coef();
};

const int Spherical_Atom::first_ex()
{
    if (atomic_number == 1)
        return 0;
    else if (atomic_number > 113)
        return 200000000;
    int ex = 0;
    for (int i = 0; i < atomic_number - 1; i++)
        ex += nex[i];
    return ex;
};

const int Spherical_Atom::previous_element_coef()
{
    if (atomic_number <= 2)
        return 0;
    int counter = 0;
    for (int temp = atomic_number - 2; temp >= 0; temp--)
    {
        for (int m = 0; m < 7; m++)
            if (occ[temp * 19 + 0 + m] != 0)
                counter += ns[temp];
        for (int m = 0; m < 6; m++)
            if (occ[temp * 19 + 7 + m] != 0)
                counter += np[temp];
        for (int m = 0; m < 4; m++)
            if (occ[temp * 19 + 13 + m] != 0)
                counter += nd[temp];
        for (int m = 0; m < 2; m++)
            if (occ[temp * 19 + 17 + m] != 0)
                counter += nf[temp];
    }
    return counter;
};

void Thakkar::calc_orbs(
    int& nr_ex,
    int& nr_coef,
    const double& dist,
    const int& offset,
    const int* n_vector,
    const int lower_m,
    const int upper_m,
    double* Orb)
{
    double exponent;
    for (int ex = 0; ex < n_vector[atomic_number - 1]; ex++)
    {
        for (int m = lower_m; m < upper_m; m++)
        {
            if (occ[offset + m] == 0)
                continue;
            exponent = -z[nr_ex] * dist;
            if (exponent > -46.5)
            { // Corresponds to at least 1E-20
                if (n[nr_ex] == 1)
                    Orb[m] += c[nr_coef] * exp(exponent);
                else
                    Orb[m] += c[nr_coef] * pow(dist, n[nr_ex] - 1) * exp(exponent);
            }
            nr_coef++;
        }
        nr_ex++;
    }
}

void Thakkar::calc_custom_orbs(
    int& nr_ex,
    int& nr_coef,
    const double& dist,
    const int& offset,
    const int* n_vector,
    const int lower_m,
    const int upper_m,
    const int& max,
    const int& min,
    double* Orb)
{
    double exponent;
    for (int ex = 0; ex < n_vector[atomic_number - 1]; ex++)
    {
        for (int m = lower_m + min; m < upper_m; m++)
        {
            if (occ[offset + m] == 0)
                continue;
            if (m < lower_m + max)
            {
                exponent = -z[nr_ex] * dist;
                if (exponent > -46.5)
                { // Corresponds to at least 1E-20
                    if (n[nr_ex] == 1)
                        Orb[m] += c[nr_coef] * exp(exponent);
                    else
                        Orb[m] += c[nr_coef] * pow(dist, n[nr_ex] - 1) * exp(exponent);
                }
            }
            nr_coef++;
        }
        nr_ex++;
    }
}

const double Thakkar::get_radial_density(const double& dist)
{
    // Speedup things for H
    if (atomic_number == 1)
        return 6.0835 * exp(-2.3 * dist) / constants::FOUR_PI;

    double Rho = 0.0;
    if (_first_ex == 200000000)
        return -20;
    int nr_coef = _prev_coef;
    int nr_ex = _first_ex;

    double Orb[19] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    calc_orbs(nr_ex, nr_coef, dist, _offset, ns, 0, 7, Orb);
    calc_orbs(nr_ex, nr_coef, dist, _offset, np, 7, 13, Orb);
    calc_orbs(nr_ex, nr_coef, dist, _offset, nd, 13, 17, Orb);
    calc_orbs(nr_ex, nr_coef, dist, _offset, nf, 17, 19, Orb);

    for (int m = 0; m < 19; m++)
    {
        if (Orb[m] == 0 || occ[_offset + m] == 0)
            continue;
        Rho += occ[_offset + m] * pow(Orb[m], 2);
    }
    return Rho / (constants::FOUR_PI);
};

const double Thakkar::get_radial_custom_density(const double& dist,
    const int& max_s,
    const int& max_p,
    const int& max_d,
    const int& max_f,
    const int& min_s,
    const int& min_p,
    const int& min_d,
    const int& min_f)
{
    // Speedup things for H
    if (atomic_number == 1)
        return 6.0835 * exp(-2.3 * dist) / constants::FOUR_PI;

    double Rho = 0.0;
    if (_first_ex == 200000000)
        return -20;

    double Orb[19] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    int nr_coef = _prev_coef;
    int nr_ex = _first_ex;

    calc_custom_orbs(nr_ex, nr_coef, dist, _offset, ns, 0, 7, max_s, min_s, Orb);
    calc_custom_orbs(nr_ex, nr_coef, dist, _offset, np, 7, 13, max_p, min_p, Orb);
    calc_custom_orbs(nr_ex, nr_coef, dist, _offset, nd, 13, 17, max_d, min_d, Orb);
    calc_custom_orbs(nr_ex, nr_coef, dist, _offset, nf, 17, 19, max_f, min_f, Orb);

    for (int m = 0; m < 19; m++)
    {
        if (Orb[m] == 0 || occ[_offset + m] == 0)
            continue;
        Rho += occ[_offset + m] * pow(Orb[m], 2);
    }
    return Rho / (constants::FOUR_PI);
};

constexpr double cosinus_integral(const int N, const double z, const double k);

static constexpr double sinus_integral(const int N, const double z, const double k)
{
    // Calculates the integral 0 - inf r ^ N e ^ -zr sin(kr) dr through recursion using the general integral int e^ax sin(bx) dx = -e^-ax/(a^2+b^2) * (a sin(bx) + b cos(bx)) and partial integration
    if (N == 0)
        return k / (z * z + k * k);
    else
        return N / (z * z + k * k) * (z * sinus_integral(N - 1, z, k) + k * cosinus_integral(N - 1, z, k));
};

constexpr double cosinus_integral(const int N, const double z, const double k)
{
    // Calculates the integral 0 - inf r ^ N e ^ -zr cos(kr) dr through recursion using the general integral int e^ax cos(bx) dx = -e^-ax/(a^2+b^2) * (a cos(bx) - b sin(bx)) and partial integration
    if (N == 0)
        return z / (z * z + k * k);
    else
        return N / (z * z + k * k) * (z * cosinus_integral(N - 1, z, k) - k * sinus_integral(N - 1, z, k));
};

const double Thakkar::get_form_factor(const double& k_vector)
{
    return get_custom_form_factor(k_vector, 7, 6, 4, 2, 0, 0, 0, 0);
};

void set_core_counts(int* max_s, int* max_p, int* max_d, int* max_f, const int& core_els, const int& ECP_mode)
{
    if (ECP_mode == 1 || ECP_mode == 3)
    {
        // If there are no core electrons, there is nothing to set for core orbital counts, so return early.
        if (core_els == 0)
            return;
        if (core_els == 2)
        {
            *max_s = 1;
            *max_p = 0;
            *max_d = 0;
            *max_f = 0;
        }
        else if (core_els == 10)
        {
            *max_s = 2;
            *max_p = 1;
            *max_d = 0;
            *max_f = 0;
        }
        else if (core_els == 18)
        {
            *max_s = 3;
            *max_p = 2;
            *max_d = 0;
            *max_f = 0;
        }
        else if (core_els == 28)
        {
            *max_s = 3;
            *max_p = 2;
            *max_d = 1;
            *max_f = 0;
        }
        else if (core_els == 46)
        {
            *max_s = 4;
            *max_p = 3;
            *max_d = 2;
            *max_f = 0;
        }
        else if (core_els == 60)
        {
            *max_s = 4;
            *max_p = 3;
            *max_d = 2;
            *max_f = 1;
        }
        else if (core_els == 78)
        {
            *max_s = 5;
            *max_p = 4;
            *max_d = 3;
            *max_f = 1;
        }
        else
        {
            err_not_impl_f("PLEASE DO NOT MAKE MORE ELECTRONS CORE ELECTRONS!", std::cout);
            return;
        }
    }
    else if (ECP_mode == 2)
    {
        if (core_els == 0)
            return;
        if (core_els == 2)
        {
            *max_s = 1;
            *max_p = 0;
            *max_d = 0;
            *max_f = 0;
        }
        else if (core_els == 10)
        {
            *max_s = 2;
            *max_p = 1;
            *max_d = 0;
            *max_f = 0;
        }
        else if (core_els == 18)
        {
            *max_s = 3;
            *max_p = 2;
            *max_d = 0;
            *max_f = 0;
        }
        else if (core_els == 28)
        {
            *max_s = 3;
            *max_p = 2;
            *max_d = 1;
            *max_f = 0;
        }
        else if (core_els == 36)
        {
            *max_s = 4;
            *max_p = 3;
            *max_d = 1;
            *max_f = 0;
        }
        else if (core_els == 46)
        {
            *max_s = 4;
            *max_p = 3;
            *max_d = 2;
            *max_f = 0;
        }
        else if (core_els == 54)
        {
            *max_s = 5;
            *max_p = 4;
            *max_d = 2;
            *max_f = 0;
        }
        else if (core_els > 54 && core_els < 69)
        {
            *max_s = 5;
            *max_p = 4;
            *max_d = 2;
            *max_f = 1;
        }
        else if (core_els == 78)
        {
            *max_s = 5;
            *max_p = 4;
            *max_d = 3;
            *max_f = 1;
        }
        else
        {
            err_not_impl_f("PLEASE DO NOT MAKE MORE ELECTRONS CORE ELECTRONS!", std::cout);
            return;
        }
    }
    else
    {
        err_not_impl_f("PLEASE DO NOT MAKE MORE ELECTRONS CORE ELECTRONS!", std::cout);
        return;
    }
};

const double Thakkar::get_core_form_factor(const double& k_vector, const int& core_els)
{
    int max_s = 0, max_p = 0, max_d = 0, max_f = 0;
    set_core_counts(&max_s, &max_p, &max_d, &max_f, core_els, ECP_mode);
    return get_custom_form_factor(k_vector, max_s, max_p, max_d, max_f, 0, 0, 0, 0);
};

const double Thakkar::get_core_density(const double& dist, const int& core_els)
{
    int max_s = 0, max_p = 0, max_d = 0, max_f = 0;
    set_core_counts(&max_s, &max_p, &max_d, &max_f, core_els, ECP_mode);
    return get_radial_custom_density(dist, max_s, max_p, max_d, max_f, 0, 0, 0, 0);
};

static double calc_int(const int& occ, const double& coef, const double& exp, const int& radial_exp, const double& k_vector)
{
    return occ * coef * sinus_integral(radial_exp, exp, k_vector);
}

static double calc_int_at_k0(const int& occ, const double& coef, const double& exp, const int& radial_exp, const double&)
{
    return occ * coef * constants::ft[radial_exp] / pow(exp, radial_exp + 1);
}

double Thakkar::calc_type(
    int& nr_ex,
    int& nr_coef,
    const double& k_vector,
    const int& offset,
    const int* n_vector,
    const int lower_m,
    const int upper_m,
    const int& max,
    const int& min)
{

    std::function<double(const int&, const double&, const double&, const int&, const double&)> func;
    if (k_vector == 0)
        func = calc_int_at_k0;
    else
        func = calc_int;
    const int l_n = n_vector[atomic_number - 1];
    double temp, result = 0;
    int i_j_distance = 0;
    for (int m = lower_m; m < upper_m; m++)
        if (occ[offset + m] != 0)
            i_j_distance++;
    for (int m = lower_m + min; m < lower_m + max; m++)
    {
        if (occ[offset + m] == 0)
            continue;
        for (int i = 0; i < l_n; i++)
        {
            for (int j = 0; j < l_n - i; j++)
            {
                temp = func(occ[offset + m],
                    c[nr_coef + m - lower_m + i * i_j_distance] * c[nr_coef + m - lower_m + (i + j) * i_j_distance],
                    z[nr_ex + i] + z[nr_ex + i + j],
                    n[nr_ex + i] + n[nr_ex + i + j] - 1,
                    k_vector);
                if (j != 0)
                    result += 2 * temp;
                else
                    result += temp;
            }
        }
    }
    nr_coef += i_j_distance * l_n;
    nr_ex += l_n;
    return result;
}

const double Thakkar::get_custom_form_factor(
    const double& k_vector,
    const int& max_s,
    const int& max_p,
    const int& max_d,
    const int& max_f,
    const int& min_s,
    const int& min_p,
    const int& min_d,
    const int& min_f)
{

    double result(0.0);

    if (_first_ex == 200000000)
        return -20;
    int nr_coef = _prev_coef;
    int nr_ex = _first_ex;

    result += calc_type(nr_ex, nr_coef, k_vector, _offset, ns, 0, 7, max_s, min_s);
    result += calc_type(nr_ex, nr_coef, k_vector, _offset, np, 7, 13, max_p, min_p);
    result += calc_type(nr_ex, nr_coef, k_vector, _offset, nd, 13, 17, max_d, min_d);
    result += calc_type(nr_ex, nr_coef, k_vector, _offset, nf, 17, 19, max_f, min_f);

    if (k_vector != 0)
        return result / k_vector;
    else
        return result;
};

void Thakkar::make_interpolator(const double& incr, const double& min_dist) {
    lincr = log(incr);
    start = min_dist;
    double current = 1;
    double _dist = min_dist;
    while (current > 1E-12)
    {
        radial_dist.push_back(_dist);
        current = get_radial_density(_dist);
        radial_density.push_back(current);
        _dist *= incr;
    }
};

double Thakkar::get_interpolated_density(const double& dist) {
    double result = 0;
    if (dist > radial_dist.back())
        return 0;
    else if (dist < radial_dist[0])
        return radial_density[0];
    int nr = int(floor(log(dist / start) / lincr));
    result = radial_density[nr] + (radial_density[nr + 1] - radial_density[nr]) / (radial_dist[nr] - radial_dist[nr - 1]) * (dist - radial_dist[nr - 1]);
    if (result < 1E-16)
        return 0;
    return result;
};

Thakkar_Anion::Thakkar_Anion(int g_atom_number) : Thakkar(g_atom_number)
{
    if (g_atom_number != 1 && g_atom_number != 6 && g_atom_number != 8 && g_atom_number != 15 && g_atom_number != 17)
        err_not_impl_f("Only selected anions are currently defined!", std::cout);
    nex = &(Anion_nex[0]);
    _first_ex = first_ex();
    ns = &(Anion_ns[0]);
    np = &(Anion_np[0]);
    nd = &(Anion_nd[0]);
    nf = &(Thakkar_nf[0]);
    occ = &(Anion_occ[0]);
    n = &(Anion_n[0]);
    z = &(Anion_z[0]);
    c = &(Anion_c[0]);
    charge = -1;
    if (atomic_number == 1)
        _prev_coef = 0;
    else
        _prev_coef = previous_element_coef();
};

Thakkar_Cation::Thakkar_Cation(int g_atom_number) : Thakkar(g_atom_number)
{
    if (g_atom_number < 3 || g_atom_number > 29)
        err_not_impl_f("Atoms with Z < 3 or bigger than 29 are not yet done!", std::cout);
    nex = &(Cation_nex[0]);
    _first_ex = first_ex();
    ns = &(Cation_ns[0]);
    np = &(Cation_np[0]);
    nd = &(Cation_nd[0]);
    nf = &(Thakkar_nf[0]);
    occ = &(Cation_occ[0]);
    n = &(Cation_n[0]);
    z = &(Cation_z[0]);
    c = &(Cation_c[0]);
    charge = +1;
    if (atomic_number == 1)
        _prev_coef = 0;
    else
        _prev_coef = previous_element_coef();
};

const double gauss_cos_integral(const int& N, const double& exp, const double& k_vector);

// This function calcualtes the integral of Int_0^Inf r^N exp(-zr^2) sin(kr) dr using a recursion of sinus and cosinus integrals with lower exponents of r
static const double gauss_sin_integral(const int& N, const double& exp, const double& k_vector)
{
    if (N == 1)
    {
        return k_vector * std::exp(-k_vector * k_vector / (4. * exp)) * pow(constants::PI / exp, 3. / 2.);
    }
    else
        return k_vector / (2. * exp) * gauss_cos_integral(N - 1, exp, k_vector) + (N - 1) / (2. * exp) * gauss_sin_integral(N - 2, exp, k_vector);
};
const double gauss_cos_integral(const int& N, const double& exp, const double& k_vector)
{
    if (N == 0)
    {
        return constants::sqr_pi * std::exp(-k_vector * k_vector / 4. / exp) / 2. / pow(exp, 1. / 2.);
    }
    else
        return (N - 1) / (2. * exp) * gauss_cos_integral(N - 2, exp, k_vector) + k_vector / (2. * exp) * gauss_sin_integral(N - 1, exp, k_vector);
};

// the integral in case of a gaussian function should be 1/k int r^(n) e^(-exp * r^2) sin(kr) dr
static double calc_Gaussian_int(const int& occ, const double& coef, const double& exp, const int& radial_exp, const double& k_vector)
{
    return occ * coef * gauss_sin_integral(radial_exp, exp, k_vector) / k_vector;
}

static double calc_Gaussian_int_at_k0(const int& occ, const double& coef, const double& exp, const int& radial_exp, const double& k_vector)
{
    const int N = radial_exp;
    return -pow(2.0, -N - 2.5) * pow(exp, -N - 1.5) * tgamma(N + 1.5) * coef * occ;
    (void)k_vector;
}

Gaussian_Atom::Gaussian_Atom(int g_atom_number, std::string& basis) : Spherical_Atom(g_atom_number, 2)
{
    _offset = (atomic_number - 37) * 21;
    if (basis == "def2-ECP")
    {
        first_atomic_number = 37;
        nex = &(def2_nex[0]);
        _first_ex = first_ex();
        ns = &(def2_ns[0]);
        np = &(def2_np[0]);
        nd = &(def2_nd[0]);
        nf = &(def2_nf[0]);
        ng = &(def2_ng[0]);
        nh = &(def2_nh[0]);
#pragma omp single
        {
            if (def2_n.size() == 0)
            {
                for (int i = 36; i < 86; i++)
                {
                    for (int s = 0; s < def2_ns[i]; s++)
                        def2_n.push_back(0);
                    for (int p = 0; p < def2_np[i]; p++)
                        def2_n.push_back(1);
                    for (int d = 0; d < def2_nd[i]; d++)
                        def2_n.push_back(2);
                    for (int f = 0; f < def2_nf[i]; f++)
                        def2_n.push_back(3);
                    for (int g = 0; g < def2_ng[i]; g++)
                        def2_n.push_back(4);
                    for (int h = 0; h < def2_nh[i]; h++)
                        def2_n.push_back(5);
                }
            }
        }
        occ = &(def2_occ[0]);
        n = def2_n.data();
        z = &(def2_z[0]);
        c = &(def2_c[0]);
    }
    else
    {
        err_checkf(false, "Basis set not implemented", std::cout);
        // Just to silence intellisense
        first_atomic_number = 37;
        nex = &(def2_nex[0]);
        _first_ex = first_ex();
        ns = &(def2_ns[0]);
        np = &(def2_np[0]);
        nd = &(def2_nd[0]);
        nf = &(def2_nf[0]);
        ng = &(def2_ng[0]);
        nh = &(def2_nh[0]);
#pragma omp single
        {
            if (def2_n.size() == 0)
            {
                for (int i = 36; i < 86; i++)
                {
                    for (int s = 0; s < def2_ns[i]; s++)
                        def2_n.push_back(0);
                    for (int p = 0; p < def2_np[i]; p++)
                        def2_n.push_back(1);
                    for (int d = 0; d < def2_nd[i]; d++)
                        def2_n.push_back(2);
                    for (int f = 0; f < def2_nf[i]; f++)
                        def2_n.push_back(3);
                    for (int g = 0; g < def2_ng[i]; g++)
                        def2_n.push_back(4);
                    for (int h = 0; h < def2_nh[i]; h++)
                        def2_n.push_back(5);
                }
            }
        }
        occ = &(def2_occ[0]);
        n = def2_n.data();
        z = &(def2_z[0]);
        c = &(def2_c[0]);
    }
    if (atomic_number == first_atomic_number)
        _prev_coef = 0;
    else
        _prev_coef = previous_element_coef();
};

double Gaussian_Atom::calc_type(
    int& nr_ex,
    int& nr_coef,
    const double& k_vector,
    const int& offset,
    const int* n_vector,
    const int lower_m,
    const int upper_m,
    const int& max,
    const int& min)
{

    std::function<double(const int&, const double&, const double&, const int&, const double&)> func;
    if (k_vector == 0)
        func = calc_Gaussian_int_at_k0;
    else
        func = calc_Gaussian_int;
    const int l_n = n_vector[atomic_number - 1];
    double temp, result = 0;
    int i_j_distance = 0;
    for (int m = lower_m; m < upper_m; m++)
        if (occ[offset + m] != 0)
            i_j_distance++;
    for (int m = lower_m + min; m < lower_m + max; m++)
    {
        if (occ[offset + m] == 0)
            continue;
        for (int i = 0; i < l_n; i++)
        {
            for (int j = 0; j < l_n - i; j++)
            {
                temp = func(occ[offset + m],
                    c[nr_coef + m - lower_m + i * i_j_distance] * c[nr_coef + m - lower_m + (i + j) * i_j_distance],
                    z[nr_ex + i] + z[nr_ex + i + j],
                    n[nr_ex + i] + n[nr_ex + i + j] + 1,
                    k_vector);
                if (j != 0)
                    result += 2 * temp;
                else
                    result += temp;
            }
        }
    }
    nr_coef += i_j_distance * l_n;
    nr_ex += l_n;
    return result;
}

const double Gaussian_Atom::get_custom_form_factor(
    const double& k_vector,
    const int& max_s,
    const int& max_p,
    const int& max_d,
    const int& max_f,
    const int& min_s,
    const int& min_p,
    const int& min_d,
    const int& min_f)
{
    err_not_impl_SA();
    (void)k_vector;
    (void)max_s;
    (void)max_p;
    (void)max_d;
    (void)max_f;
    (void)min_s;
    (void)min_p;
    (void)min_d;
    (void)min_f;
    return -1;
};
const double Gaussian_Atom::get_custom_form_factor(
    const double& k_vector,
    const int& max_s,
    const int& max_p,
    const int& max_d,
    const int& max_f,
    const int& max_g,
    const int& max_h,
    const int& min_s,
    const int& min_p,
    const int& min_d,
    const int& min_f,
    const int& min_g,
    const int& min_h)
{

    if (k_vector != 0)
    {
        double result(0.0);

        if (_first_ex == 200000000)
            return -20;
        int nr_coef = _prev_coef;
        int nr_ex = _first_ex;

        result += calc_type(nr_ex, nr_coef, k_vector, _offset, ns, 0, 1, max_s, min_s);
        result += calc_type(nr_ex, nr_coef, k_vector, _offset, np, 1, 2, max_p, min_p);
        result += calc_type(nr_ex, nr_coef, k_vector, _offset, nd, 2, 3, max_d, min_d);
        result += calc_type(nr_ex, nr_coef, k_vector, _offset, nf, 3, 4, max_f, min_f);
        result += calc_type(nr_ex, nr_coef, k_vector, _offset, ng, 4, 5, max_g, min_g);
        result += calc_type(nr_ex, nr_coef, k_vector, _offset, nh, 5, 6, max_h, min_h);

        return result / k_vector;
    }
    else
        return this->get_custom_form_factor(1E-12, max_s, max_p, max_d, max_f, max_g, max_h, min_s, min_p, min_d, min_f, min_g, min_h);
};

const double Gaussian_Atom::get_form_factor(const double& k_vector)
{
    return get_custom_form_factor(k_vector, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0);
};

const double Gaussian_Atom::get_core_form_factor(const double& k_vector, const int& core_els)
{
    int max_s = 0, max_p = 1, max_d = 1, max_f = 1, max_h = 1, max_g = 1;
    return get_custom_form_factor(k_vector, max_s, max_p, max_d, max_f, max_g, max_h, 0, 0, 0, 0, 0, 0);
    (void)core_els;
};

void Gaussian_Atom::calc_orbs(
    int& nr_ex,
    int& nr_coef,
    const double& dist,
    const int& offset,
    const int* n_vector,
    const int lower_m,
    const int upper_m,
    double* Orb)
{
    double exponent;
    for (int ex = 0; ex < n_vector[atomic_number - 1]; ex++)
    {
        for (int m = lower_m; m < upper_m; m++)
        {
            if (occ[offset + m] == 0)
                continue;
            exponent = -z[nr_ex] * dist * dist;
            if (exponent > -46.5)
            { // Corresponds to at least 1E-20
                if (n[nr_ex] == 1)
                    Orb[m] += c[nr_coef] * exp(exponent);
                else
                    Orb[m] += c[nr_coef] * pow(dist, n[nr_ex]) * exp(exponent);
            }
            nr_coef++;
        }
        nr_ex++;
    }
}

const int Gaussian_Atom::previous_element_coef()
{
    if (atomic_number <= first_atomic_number)
        return 0;
    int counter = 0;
    for (int temp = atomic_number - first_atomic_number - 1; temp >= 0; temp--)
    {
        if (occ[temp * 6 + 0] != 0)
            counter += ns[temp];
        if (occ[temp * 6 + 1] != 0)
            counter += np[temp];
        if (occ[temp * 6 + 2] != 0)
            counter += nd[temp];
        if (occ[temp * 6 + 3] != 0)
            counter += nf[temp];
        if (occ[temp * 6 + 4] != 0)
            counter += ng[temp];
        if (occ[temp * 6 + 5] != 0)
            counter += nh[temp];
    }
    return counter;
};

const double Gaussian_Atom::get_radial_density(const double& dist)
{
    double Rho = 0.0;
    if (_first_ex == 200000000)
        return -20;

    double Orb[6] = { 0, 0, 0, 0, 0, 0 };
    int nr_coef = _prev_coef;
    int nr_ex = _first_ex;

    calc_orbs(nr_ex, nr_coef, dist, _offset, ns, 0, 1, Orb);
    calc_orbs(nr_ex, nr_coef, dist, _offset, np, 1, 2, Orb);
    calc_orbs(nr_ex, nr_coef, dist, _offset, nd, 2, 3, Orb);
    calc_orbs(nr_ex, nr_coef, dist, _offset, nf, 3, 4, Orb);
    calc_orbs(nr_ex, nr_coef, dist, _offset, ng, 4, 5, Orb);
    calc_orbs(nr_ex, nr_coef, dist, _offset, nh, 5, 6, Orb);

    for (int m = 0; m < 6; m++)
    {
        if (Orb[m] == 0 || occ[_offset + m] == 0)
            continue;
        Rho += occ[_offset + m] * pow(Orb[m], 2);
    }
    return Rho / (constants::FOUR_PI); // 4pi is the angular function
};

const double Spherical_Gaussian_Density::get_radial_density(const double& dist)
{
    double res = 0;
    double d2 = dist * dist;
    for (int i = 0; i < nex; i++)
    {
        res += c[i] * exp(-z[i] * d2);
    }
    return res;
}

const double Spherical_Gaussian_Density::get_form_factor(const double& k_vector)
{
    std::function<double(const int&, const double&, const double&, const int&, const double&)> func;
    if (k_vector == 0)
        func = calc_Gaussian_int_at_k0;
    else
        func = calc_Gaussian_int;
    double res = 0;
    for (int i = 0; i < nex; i++)
    {
        res += func(1, c[i], z[i], 1, k_vector);
    }
    return res;
}
