#pragma once
//One set of point functions for every density NoSpherA2 knows: a WFN with orbitals, a fitted Gaussian density
//(Gaussian_Atom, Gaussian_Molecule) and an atom model placed at a nucleus (Centred<Thakkar>, Centred<MBIS_Atom>,
//Centred<EMBIS_Atom>, ...). A source is anything with rho(p); the derived properties (Laplacian, reduced gradient,
//ELI-D) are written once, here, on top of calculate_density.
#include "convenience.h"
#include "wfn_class.h"
#include "aux_density.h"
#include "constants.h"
#include <concepts>
#include <functional>

template<class S> concept RadialModel = requires(const S& s, const double r) { { s.get_radial_density(r) } -> std::convertible_to<double>; };
template<class S> concept RadialModelWithDerivatives = requires(const S& s, const double r, double& d1, double& d2) { { s.get_radial_density(r, d1, d2) } -> std::convertible_to<double>; };
template<class S> concept RelativeModel = requires(const S& s, const d3& p) { { s.get_density(p) } -> std::convertible_to<double>; };
template<class S> concept RelativeModelWithDerivatives = requires(const S& s, const d3& p, d3& g, double& lap) { { s.get_density(p, g, lap) } -> std::convertible_to<double>; };
template<class S> concept PointDensity = requires(const S& s, const d3& p) { { s.rho(p) } -> std::convertible_to<double>; };
template<class S> concept PointDensityWithDerivatives = requires(const S& s, const d3& p, d3& g, double& lap) { { s.values(p, g, lap) } -> std::convertible_to<double>; };
template<class S> concept RelativeModelWithHessian = requires(const S& s, const d3& p, d3& g, double* H) { { s.get_density(p, g, H) } -> std::convertible_to<double>; };
template<class S> concept PointDensityWithHessian = requires(const S& s, const d3& p, d3& g, double* H) { { s.hessian(p, g, H) } -> std::convertible_to<double>; };

//A spherical (get_radial_density(r)) or EMBIS (get_density of the nucleus-relative position) atom model at a nucleus;
//holds a reference, so the model outlives it
template<class A> struct Centred
{
    const A& model;
    d3 centre;
    double rho(const d3& p) const requires RadialModel<A> { return model.get_radial_density(array_length(p, centre)); }
    double rho(const d3& p) const requires RelativeModel<A> { return model.get_density({ p[0] - centre[0], p[1] - centre[1], p[2] - centre[2] }); }
    double values(const d3& p, d3& grad, double& lap) const requires RelativeModelWithDerivatives<A> { return model.get_density({ p[0] - centre[0], p[1] - centre[1], p[2] - centre[2] }, grad, lap); }
    double hessian(const d3& p, d3& grad, double* H) const requires RelativeModelWithHessian<A> { return model.get_density({ p[0] - centre[0], p[1] - centre[1], p[2] - centre[2] }, grad, H); }
};

//rho at p
inline double calculate_density(const WFN& w, const d3& p) { return w.compute_dens(p); }
template<PointDensity S> double calculate_density(const S& s, const d3& p) { return s.rho(p); }

//rho at p with its gradient and Laplacian
inline double calculate_density(const WFN& w, const d3& p, d3& grad, double& lap)
{
    double rho, tau, hess[9];
    w.computeValues(p, rho, grad, hess, tau);
    lap = hess[0] + hess[4] + hess[8];
    return rho;
}
template<PointDensityWithDerivatives S> double calculate_density(const S& s, const d3& p, d3& grad, double& lap) { return s.values(p, grad, lap); }
//ponytail: central differences for a source that only has rho; analytic derivatives when a case needs them
template<PointDensity S> requires (!PointDensityWithDerivatives<S>)
double calculate_density(const S& s, const d3& p, d3& grad, double& lap)
{
    const double h = 1E-4, rho = s.rho(p);
    lap = 0.0;
    for (int k = 0; k < 3; k++) {
        d3 pp = p, pm = p;
        pp[k] += h;
        pm[k] -= h;
        const double rp = s.rho(pp), rm = s.rho(pm);
        grad[k] = (rp - rm) / (2 * h);
        lap += (rp + rm - 2 * rho) / (h * h);
    }
    return rho;
}
//rho(r), rho', rho'' of a radial model: analytic where the model has them, otherwise one central difference
//(the model is even in r)
template<RadialModelWithDerivatives A> double radial_derivatives(const A& m, const double r, double& d1, double& d2) { return m.get_radial_density(r, d1, d2); }
template<RadialModel A> requires (!RadialModelWithDerivatives<A>)
double radial_derivatives(const A& m, const double r, double& d1, double& d2)
{
    const double h = 1E-4, rho = m.get_radial_density(r);
    const double rp = m.get_radial_density(r + h), rm = m.get_radial_density(std::abs(r - h));
    d1 = (rp - rm) / (2 * h), d2 = (rp + rm - 2 * rho) / (h * h);
    return rho;
}
//A radial model: lap = rho'' + 2 rho' / r
template<RadialModel A> double calculate_density(const Centred<A>& s, const d3& p, d3& grad, double& lap)
{
    const double h = 1E-4;
    const d3 d{ p[0] - s.centre[0], p[1] - s.centre[1], p[2] - s.centre[2] };
    double d1, d2;
    const double r = array_length(d), rho = radial_derivatives(s.model, r, d1, d2);
    if (r < h) {
        grad = { 0.0, 0.0, 0.0 };
        lap = 3 * d2;
    }
    else {
        for (int k = 0; k < 3; k++) grad[k] = d1 * d[k] / r;
        lap = d2 + 2 * d1 / r;
    }
    return rho;
}

inline double calculate_laplacian(const WFN& w, const d3& p) { return w.computeLap(p); }
template<class S> double calculate_laplacian(const S& s, const d3& p)
{
    d3 g;
    double lap;
    calculate_density(s, p, g, lap);
    return lap;
}
//Reduced density gradient |grad rho| / (2 (3 pi^2)^(1/3) rho^(4/3)), the NCI descriptor
template<class S> double calculate_rdg(const S& s, const d3& p)
{
    d3 g;
    double lap;
    const double rho = calculate_density(s, p, g, lap);
    return rho > 0 ? constants::alpha_coef * std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]) / std::pow(rho, constants::c_43) : 0.0;
}
//ELI-D: exact from the orbitals, the PC07 orbital-free estimate for a density alone
inline double calculate_eli(const WFN& w, const d3& p) { return w.computeELI(p); }
template<PointDensity S> double calculate_eli(const S& s, const d3& p)
{
    d3 g;
    double lap;
    const double rho = calculate_density(s, p, g, lap);
    return aux_density::eli_from_density(rho, g[0] * g[0] + g[1] * g[1] + g[2] * g[2], lap);
}

//rho at n points, in parallel over the points; the WFN keeps its scratch per thread
template<class S> void calculate_density(const S& s, const int n, const double* x, const double* y, const double* z, double* rho)
{
#pragma omp parallel for schedule(dynamic, 64)
    for (int p = 0; p < n; p++) rho[p] = calculate_density(s, d3{ x[p], y[p], z[p] });
}
//Where a source sits, for the "within radius of a centre" cut of the cube functions, and its name for messages
inline std::vector<d3> source_positions(const WFN& w)
{
    std::vector<d3> c;
    for (const atom& a : w.get_atoms()) c.push_back(a.get_pos());
    return c;
}
template<class A> std::vector<d3> source_positions(const Centred<A>& s) { return { s.centre }; }
inline const char* source_name(const WFN&) { return "wavefunction"; }
template<class A> const char* source_name(const Centred<A>&) { return "atom model"; }
//A source that also carries a potential (the fitted density with its nuclei)
template<class S> concept PointPotential = requires(const S& s, const d3& p) { { s.esp(p) } -> std::convertible_to<double>; };

//Hessian of rho at p, row-major 3x3
inline void calculate_hessian(const WFN& w, const d3& p, double* H)
{
    double rho, tau;
    d3 g;
    w.computeValues(p, rho, g, H, tau);
}
template<PointDensityWithHessian S> void calculate_hessian(const S& s, const d3& p, double* H)
{
    d3 g;
    s.hessian(p, g, H);
}
//A radial model: H = rho'' u u^T + rho'/r (1 - u u^T); at the nucleus every direction is radial, H = rho'' 1
template<RadialModel A> void calculate_hessian(const Centred<A>& s, const d3& p, double* H)
{
    const double h = 1E-4;
    const d3 d{ p[0] - s.centre[0], p[1] - s.centre[1], p[2] - s.centre[2] };
    double d1, d2;
    const double r = array_length(d);
    radial_derivatives(s.model, r, d1, d2);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            const double delta = i == j ? 1.0 : 0.0;
            if (r < h) H[3 * i + j] = d2 * delta;
            else {
                const double uu = d[i] * d[j] / (r * r);
                H[3 * i + j] = d2 * uu + d1 / r * (delta - uu);
            }
        }
}
//ponytail: central differences of the gradient, six evaluations, for a source without an analytic Hessian (none today)
template<class S> requires (!PointDensityWithHessian<S>)
void calculate_hessian(const S& s, const d3& p, double* H)
{
    const double h = 1E-4;
    for (int k = 0; k < 3; k++) {
        d3 pp = p, pm = p, gp, gm;
        double lap;
        pp[k] += h;
        pm[k] -= h;
        calculate_density(s, pp, gp, lap);
        calculate_density(s, pm, gm, lap);
        for (int j = 0; j < 3; j++) H[3 * k + j] = (gp[j] - gm[j]) / (2 * h);
    }
}

//What a GridManager evaluates on its atomic grids; captures s by reference
using DensityBatch = std::function<void(int, const double*, const double*, const double*, double*)>;
template<class S> DensityBatch density_batch(const S& s)
{
    return [&s](const int n, const double* x, const double* y, const double* z, double* rho) { calculate_density(s, n, x, y, z, rho); };
}
