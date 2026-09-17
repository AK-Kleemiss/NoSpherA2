#pragma once
#include <cmath>
#include "spherical_harmonic.h"

//The fitted density at one point from a flattened auxiliary basis. One text for the CPU loop
//in SALTED_utilities.cpp and the kernel in aux_density_gpu.cu, so the two cannot drift; only
//spherical_harmonic.h is included so the .cu stays parseable by nvcc and hipcc.
namespace aux_density
{
    //Value with its three derivatives, so constants::spherical_harmonic differentiates itself:
    //instantiated on double for the density, on dual for its gradient
    struct dual
    {
        double v, x, y, z;
        AUX_HD dual(const double v_ = 0.0, const double x_ = 0.0, const double y_ = 0.0, const double z_ = 0.0) : v(v_), x(x_), y(y_), z(z_) {}
    };
    AUX_HD inline dual operator+(const dual& a, const dual& b) { return dual(a.v + b.v, a.x + b.x, a.y + b.y, a.z + b.z); }
    AUX_HD inline dual operator-(const dual& a, const dual& b) { return dual(a.v - b.v, a.x - b.x, a.y - b.y, a.z - b.z); }
    AUX_HD inline dual operator-(const dual& a) { return dual(-a.v, -a.x, -a.y, -a.z); }
    AUX_HD inline dual operator*(const dual& a, const dual& b) { return dual(a.v * b.v, a.v * b.x + a.x * b.v, a.v * b.y + a.y * b.v, a.v * b.z + a.z * b.v); }
    //rho at (x, y, z). Atoms carry their centre, the squared distance beyond which their most
    //diffuse primitive is below 1e-20 and their shell range; shells carry l, the primitive range
    //and the offset of their 2l+1 coefficients; primitives the exponent and the normalised
    //contraction coefficient. The 1e-10 radial cutoff is the one calc_density_ML always had.
    AUX_HD inline double at(const double x, const double y, const double z, const int n_at,
        const double* cx, const double* cy, const double* cz, const double* r2_max,
        const int* sh_start, const int* sh_l, const int* pr_start, const int* coef_off,
        const double* pr_exp, const double* pr_norm, const double* coefs)
    {
        double dens = 0.0;
        for (int a = 0; a < n_at; a++) {
            const double dx = x - cx[a], dy = y - cy[a], dz = z - cz[a], r2 = dx * dx + dy * dy + dz * dz;
            if (r2 > r2_max[a]) continue;
            const double r = sqrt(r2), ux = dx / r, uy = dy / r, uz = dz / r;
            for (int s = sh_start[a]; s < sh_start[a + 1]; s++) {
                double radial = 0.0, rl = 1.0;
                for (int p = pr_start[s]; p < pr_start[s + 1]; p++) radial += exp(-pr_exp[p] * r2) * pr_norm[p];
                for (int i = 0; i < sh_l[s]; i++) rl *= r;
                radial *= rl;
                if (std::abs(radial) < 1E-10) continue;
                dens += radial * constants::spherical_harmonic(sh_l[s], ux, uy, uz, coefs + coef_off[s]);
            }
        }
        return dens;
    }
    //rho, its gradient and, with LAP, its Laplacian at (x, y, z). Each shell is R(r^2) r^l Y(u) with R the
    //contraction of Gaussians, so d/dr of the radial part is analytic and dY/du comes from the dual harmonic,
    //projected onto the sphere: grad = d (-2 sum(a n e) r^l + l R r^(l-2)) Y + R r^(l-1) (1 - u u^T) dY/du.
    //r^l Y is a solid harmonic, hence harmonic itself, so lap = (4 r^2 sum(a^2 n e) - 2 (2l+3) sum(a n e)) r^l Y
    template <bool LAP>
    AUX_HD inline double at_deriv(const double x, const double y, const double z, const int n_at,
        const double* cx, const double* cy, const double* cz, const double* r2_max,
        const int* sh_start, const int* sh_l, const int* pr_start, const int* coef_off,
        const double* pr_exp, const double* pr_norm, const double* coefs, double& gx, double& gy, double& gz, double& lap)
    {
        double dens = 0.0;
        gx = gy = gz = lap = 0.0;
        for (int a = 0; a < n_at; a++) {
            const double dx = x - cx[a], dy = y - cy[a], dz = z - cz[a], r2 = dx * dx + dy * dy + dz * dz;
            if (r2 > r2_max[a]) continue;
            const double r = sqrt(r2), ux = dx / r, uy = dy / r, uz = dz / r;
            for (int s = sh_start[a]; s < sh_start[a + 1]; s++) {
                const int l = sh_l[s];
                double radial = 0.0, dradial = 0.0, ddradial = 0.0, rl = 1.0;
                for (int p = pr_start[s]; p < pr_start[s + 1]; p++) {
                    const double e = exp(-pr_exp[p] * r2) * pr_norm[p];
                    radial += e, dradial += pr_exp[p] * e;
                    if (LAP) ddradial += pr_exp[p] * pr_exp[p] * e;
                }
                for (int i = 0; i < l; i++) rl *= r;
                radial *= rl;
                if (std::abs(radial) < 1E-10) continue;
                const dual Y = constants::spherical_harmonic(l, dual(ux, 1.0, 0.0, 0.0), dual(uy, 0.0, 1.0, 0.0), dual(uz, 0.0, 0.0, 1.0), coefs + coef_off[s]);
                const double fr = (l * radial / r2 - 2.0 * dradial * rl) * Y.v, fa = radial / r, ud = ux * Y.x + uy * Y.y + uz * Y.z;
                dens += radial * Y.v;
                gx += fr * dx + fa * (Y.x - ux * ud);
                gy += fr * dy + fa * (Y.y - uy * ud);
                gz += fr * dz + fa * (Y.z - uz * ud);
                if (LAP) lap += (4.0 * r2 * ddradial - 2.0 * (2 * l + 3) * dradial) * rl * Y.v;
            }
        }
        return dens;
    }
    AUX_HD inline double at_grad(const double x, const double y, const double z, const int n_at,
        const double* cx, const double* cy, const double* cz, const double* r2_max,
        const int* sh_start, const int* sh_l, const int* pr_start, const int* coef_off,
        const double* pr_exp, const double* pr_norm, const double* coefs, double& gx, double& gy, double& gz)
    {
        double lap;
        return at_deriv<false>(x, y, z, n_at, cx, cy, cz, r2_max, sh_start, sh_l, pr_start, coef_off, pr_exp, pr_norm, coefs, gx, gy, gz, lap);
    }
    AUX_HD inline double at_lap(const double x, const double y, const double z, const int n_at,
        const double* cx, const double* cy, const double* cz, const double* r2_max,
        const int* sh_start, const int* sh_l, const int* pr_start, const int* coef_off,
        const double* pr_exp, const double* pr_norm, const double* coefs, double& gx, double& gy, double& gz, double& lap)
    {
        return at_deriv<true>(x, y, z, n_at, cx, cy, cz, r2_max, sh_start, sh_l, pr_start, coef_off, pr_exp, pr_norm, coefs, gx, gy, gz, lap);
    }
}
