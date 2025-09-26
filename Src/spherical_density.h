#pragma once

#include <vector>
#include <iostream>

#include "ECPs_corrections.h"

inline void not_implemented_SA(const std::string &file, const int &line, const std::string &function, const std::string &error_mesasge, std::ostream &log_file)
{
    log_file << function << " at: " << file << ":" << line << " " << error_mesasge << std::endl;
    log_file.flush();
    exit(-1);
};
#define err_not_impl_SA() not_implemented_SA(__FILE__, __LINE__, __func__, "Virtual_function", std::cout);

inline double linear_interpolate_spherical_density(
    const vec &radial_dens,
    const vec &spherical_dist,
    const double dist,
    const double lincr,
    const double start)
{
    double result = 0;
    if (dist > spherical_dist[spherical_dist.size() - 1])
        return 0;
    else if (dist < spherical_dist[0])
        return radial_dens[0];
    int nr = int(floor(log(dist / start) / lincr));
    result = radial_dens[nr] + (radial_dens[nr + 1] - radial_dens[nr]) / (spherical_dist[nr] - spherical_dist[nr - 1]) * (dist - spherical_dist[nr - 1]);
    if (result < 1E-10)
        result = 0;
    return result;
}

class Spherical_Atom
{
protected:
    int atomic_number;
    int ECP_mode;
    const int first_ex();
    virtual const int previous_element_coef();
    const int *nex, *ns, *np, *nd, *nf, *occ, *n;
    const double *z, *c;
    int charge;
    int _prev_coef, _offset, _first_ex;
    vec radial_density, radial_dist;
    double lincr, start;
    virtual void calc_orbs(int &nr_ex,
                           int &nr_coef,
                           const double &dist,
                           const int &offset,
                           const int *n_vector,
                           const int lower_m,
                           const int upper_m,
                           double *Orb)
    {
        err_not_impl_SA();
        (void)nr_ex;
        (void)nr_coef;
        (void)dist;
        (void)offset;
        (void)n_vector;
        (void)lower_m;
        (void)upper_m;
        (void)Orb;
    };
    virtual double calc_type(
        int &nr_ex,
        int &nr_coef,
        const double &k_vector,
        const int &offset,
        const int *n_vector,
        const int lower_m,
        const int upper_m,
        const int &max,
        const int &min)
    {
        err_not_impl_SA();
        (void)nr_ex;
        (void)nr_coef;
        (void)k_vector;
        (void)offset;
        (void)n_vector;
        (void)lower_m;
        (void)upper_m;
        (void)max;
        (void)min;
        return -1;
    };

public:
    Spherical_Atom(const int g_atom_number, const int ECP_m = 1) : atomic_number(g_atom_number),
                                                                   _offset((atomic_number - 1) * 19), _first_ex(0), _prev_coef(0), c(NULL), n(NULL), nd(NULL), ns(NULL), np(NULL), nf(NULL), nex(NULL), occ(NULL), z(NULL)
    {
        ECP_mode = ECP_m;
        charge = 0;
    };
    Spherical_Atom() : _first_ex(0), _offset(0), _prev_coef(0), c(NULL), n(NULL), nd(NULL), ns(NULL), np(NULL), nf(NULL), nex(NULL), occ(NULL), z(NULL)
    {
        ECP_mode = 1;
        atomic_number = 1;
        charge = 0;
    };
    virtual const double get_radial_density(const double &dist)
    {
        err_not_impl_SA();
        (void)dist;
        return -1;
    };
    virtual const double get_form_factor(const double &k_vector)
    {
        err_not_impl_SA();
        (void)k_vector;
        return -1;
    };
    virtual const double get_core_form_factor(const double &k_vector, const int &core_els)
    {
        err_not_impl_SA();
        (void)k_vector;
        (void)core_els;
        return -1;
    };
    virtual const double get_custom_form_factor(
        const double &k_vector,
        const int &max_s,
        const int &max_p,
        const int &max_d,
        const int &max_f,
        const int &min_s,
        const int &min_p,
        const int &min_d,
        const int &min_f)
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
    virtual void make_interpolator(const double& incr, const double& min_dist) {
        err_not_impl_SA();
        (void)incr;
        (void)min_dist;
    };
    virtual double get_interpolated_density(const double& dist) {
        err_not_impl_SA();
        (void)dist;
        return -1.0;
    };
    const int get_atomic_number() const { return atomic_number; };
    const int get_charge() const { return charge; };
    
};

class Thakkar : public Spherical_Atom
{
protected:
    void calc_orbs(int &nr_ex,
                   int &nr_coef,
                   const double &dist,
                   const int &offset,
                   const int *n_vector,
                   const int lower_m,
                   const int upper_m,
                   double *Orb) override;
    void calc_custom_orbs(int &nr_ex,
                          int &nr_coef,
                          const double &dist,
                          const int &offset,
                          const int *n_vector,
                          const int lower_m,
                          const int upper_m,
                          const int &max,
                          const int &min,
                          double *Orb);
    double calc_type(
        int &nr_ex,
        int &nr_coef,
        const double &k_vector,
        const int &offset,
        const int *n_vector,
        const int lower_m,
        const int upper_m,
        const int &max,
        const int &min) override;

public:
    Thakkar(const int g_atom_number, const int ECP_mode = 1);
    Thakkar();
    const double get_radial_density(const double &dist) override;
    const double get_radial_custom_density(
        const double &dist,
        const int &max_s,
        const int &max_p,
        const int &max_d,
        const int &max_f,
        const int &min_s,
        const int &min_p,
        const int &min_d,
        const int &min_f);
    const double get_form_factor(const double &k_vector) override;
    const double get_core_form_factor(const double &k_vector, const int &core_els) override;
    const double get_core_density(const double &dist, const int &core_els);
    const double get_custom_form_factor(
        const double &k_vector,
        const int &max_s,
        const int &max_p,
        const int &max_d,
        const int &max_f,
        const int &min_s,
        const int &min_p,
        const int &min_d,
        const int &min_f) override;
    void make_interpolator(const double& incr, const double& min_dist) override;
    double get_interpolated_density(const double& dist) override;
};

class Thakkar_Anion : public Thakkar
{
public:
    Thakkar_Anion(const int g_atom_number);
};

class Thakkar_Cation : public Thakkar
{
public:
    Thakkar_Cation(const int g_atom_number);
};

class Gaussian_Atom : public Spherical_Atom
{
protected:
    const int *ng, *nh;
    int first_atomic_number;
    const int previous_element_coef() override;
    void calc_orbs(int &nr_ex,
                   int &nr_coef,
                   const double &dist,
                   const int &offset,
                   const int *n_vector,
                   const int lower_m,
                   const int upper_m,
                   double *Orb) override;
    double calc_type(
        int &nr_ex,
        int &nr_coef,
        const double &k_vector,
        const int &offset,
        const int *n_vector,
        const int lower_m,
        const int upper_m,
        const int &max,
        const int &min) override;

public:
    Gaussian_Atom(const int g_atom_number, std::string &basis);
    Gaussian_Atom() = default;
    const double get_radial_density(const double &dist) override;
    const double get_form_factor(const double &k_vector) override;
    const double get_core_form_factor(const double &k_vector, const int &core_els) override;
    const double get_custom_form_factor(
        const double &k_vector,
        const int &max_s,
        const int &max_p,
        const int &max_d,
        const int &max_f,
        const int &min_s,
        const int &min_p,
        const int &min_d,
        const int &min_f) override;
    const double get_custom_form_factor(
        const double &k_vector,
        const int &max_s,
        const int &max_p,
        const int &max_d,
        const int &max_f,
        const int &max_g,
        const int &max_h,
        const int &min_s,
        const int &min_p,
        const int &min_d,
        const int &min_f,
        const int &min_g,
        const int &min_h);
};

class Spherical_Gaussian_Density
{
protected:
    const int atomic_number;
    const int ECP_mode;
    int nex;
    const double *z, *c;
    int charge;

public:
    Spherical_Gaussian_Density(const int g_atom_number, const int ECP_m = 1) : atomic_number(g_atom_number),
                                                                               ECP_mode(ECP_m),
                                                                               charge(0)
    {
        int temp_Z;
        switch (ECP_m)
        {
        case 2:
            temp_Z = atomic_number - 2;
            if (temp_Z < 0) {
                nex = 0;
                z = NULL;
                c = NULL;
            }
            else {
                nex = static_cast<int>(xtb_corrections::c[temp_Z].size());
                z = xtb_corrections::z[temp_Z].data();
                c = xtb_corrections::c[temp_Z].data();
            }
            break;
        case 3:
            temp_Z = atomic_number - 2;
            if (temp_Z < 0) {
                nex = 0;
                z = NULL;
                c = NULL;
            }
            else {
                nex = static_cast<int>(ptb_corrections::c[temp_Z].size());
                z = ptb_corrections::z[temp_Z].data();
                c = ptb_corrections::c[temp_Z].data();
            }
            break;
        case 1:
            temp_Z = atomic_number - 37;
            if (temp_Z < 0) {
                nex = 0;
                z = NULL;
                c = NULL;
            }
            else {
                nex = static_cast<int>(def_corrections::c[temp_Z].size());
                z = def_corrections::z[temp_Z].data();
                c = def_corrections::c[temp_Z].data();
            }
            break;
        default:
            z = NULL;
            c = NULL;
            nex = 0;
        }
    };
    Spherical_Gaussian_Density() : c(NULL), nex(0), z(NULL), charge(0), atomic_number(1), ECP_mode(1){};
    virtual const double get_radial_density(const double &dist);
    virtual const double get_form_factor(const double &k_vector);
    const int get_atomic_number() const { return atomic_number; };
    const int get_charge() const { return charge; };
};