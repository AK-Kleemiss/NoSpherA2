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
//Every source has analytic derivatives (Thakkar, MBIS, EMBIS, the fitted density, the WFN); a model without
//them (HE_Spherical_Atom, Spherical_Gaussian_Density) does not compile here rather than getting a stencil
//A radial model: lap = rho'' + 2 rho' / r
template<RadialModelWithDerivatives A> double calculate_density(const Centred<A>& s, const d3& p, d3& grad, double& lap)
{
	const double h = 1E-4;
	const d3 d{ p[0] - s.centre[0], p[1] - s.centre[1], p[2] - s.centre[2] };
	double d1, d2;
	const double r = array_length(d), rho = s.model.get_radial_density(r, d1, d2);
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

//rho at p with its gradient and Hessian (row-major 3x3); the Hessian-only form drops the rest
inline double calculate_hessian(const WFN& w, const d3& p, d3& grad, double* H)
{
	double rho, tau;
	w.computeValues(p, rho, grad, H, tau);
	return rho;
}
template<PointDensityWithHessian S> double calculate_hessian(const S& s, const d3& p, d3& grad, double* H) { return s.hessian(p, grad, H); }
template<class S> void calculate_hessian(const S& s, const d3& p, double* H)
{
	d3 g;
	calculate_hessian(s, p, g, H);
}
//A radial model, one evaluation of rho', rho'': grad = rho' u, H = rho'' u u^T + rho'/r (1 - u u^T); at the nucleus
//every direction is radial, grad = 0 and H = rho'' 1
template<RadialModelWithDerivatives A> double calculate_hessian(const Centred<A>& s, const d3& p, d3& grad, double* H)
{
	const double h = 1E-4;
	const d3 d{ p[0] - s.centre[0], p[1] - s.centre[1], p[2] - s.centre[2] };
	double d1, d2;
	const double r = array_length(d), rho = s.model.get_radial_density(r, d1, d2);
	for (int i = 0; i < 3; i++) {
		grad[i] = r < h ? 0.0 : d1 * d[i] / r;
		for (int j = 0; j < 3; j++) {
			const double delta = i == j ? 1.0 : 0.0;
			if (r < h) H[3 * i + j] = d2 * delta;
			else {
				const double uu = d[i] * d[j] / (r * r);
				H[3 * i + j] = d2 * uu + d1 / r * (delta - uu);
			}
		}
	}
	return rho;
}

//rho with its gradient and Laplacian, and with its gradient and Hessian (9 per point, row-major), at n points in
//parallel over the points; the fitted density replaces both with its batch kernel (gaussian_atom.h), on the GPU
//when one is enabled
template<class S> void calculate_density(const S& s, const int n, const double* x, const double* y, const double* z, double* rho, double* gx, double* gy, double* gz, double* lap)
{
#pragma omp parallel for schedule(dynamic, 64)
	for (int p = 0; p < n; p++) {
		d3 g;
		rho[p] = calculate_density(s, d3{ x[p], y[p], z[p] }, g, lap[p]);
		gx[p] = g[0], gy[p] = g[1], gz[p] = g[2];
	}
}
template<class S> void calculate_hessian(const S& s, const int n, const double* x, const double* y, const double* z, double* rho, double* gx, double* gy, double* gz, double* H)
{
#pragma omp parallel for schedule(dynamic, 64)
	for (int p = 0; p < n; p++) {
		d3 g;
		rho[p] = calculate_hessian(s, d3{ x[p], y[p], z[p] }, g, H + 9 * (size_t)p);
		gx[p] = g[0], gy[p] = g[1], gz[p] = g[2];
	}
}

//What a GridManager evaluates on its atomic grids; captures s by reference
using DensityBatch = std::function<void(int, const double*, const double*, const double*, double*)>;
template<class S> DensityBatch density_batch(const S& s)
{
	return [&s](const int n, const double* x, const double* y, const double* z, double* rho) { calculate_density(s, n, x, y, z, rho); };
}
