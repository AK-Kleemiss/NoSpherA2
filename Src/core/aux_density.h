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
    //Value, gradient and the six upper-triangle second derivatives (xx xy xz yy yz zz), so the harmonic differentiates
    //itself twice through u(d) = d / r
    struct hyperdual
    {
        double v, g[3], h[6];
        AUX_HD hyperdual(const double v_ = 0.0) : v(v_), g{ 0.0, 0.0, 0.0 }, h{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } {}
    };
    AUX_HD inline hyperdual operator+(const hyperdual& a, const hyperdual& b)
    {
        hyperdual r(a.v + b.v);
        for (int i = 0; i < 3; i++) r.g[i] = a.g[i] + b.g[i];
        for (int i = 0; i < 6; i++) r.h[i] = a.h[i] + b.h[i];
        return r;
    }
    AUX_HD inline hyperdual operator-(const hyperdual& a)
    {
        hyperdual r(-a.v);
        for (int i = 0; i < 3; i++) r.g[i] = -a.g[i];
        for (int i = 0; i < 6; i++) r.h[i] = -a.h[i];
        return r;
    }
    AUX_HD inline hyperdual operator-(const hyperdual& a, const hyperdual& b) { return a + (-b); }
    AUX_HD inline hyperdual operator*(const hyperdual& a, const hyperdual& b)
    {
        hyperdual r(a.v * b.v);
        for (int i = 0; i < 3; i++) r.g[i] = a.v * b.g[i] + a.g[i] * b.v;
        for (int i = 0, n = 0; i < 3; i++)
            for (int j = i; j < 3; j++, n++) r.h[n] = a.v * b.h[n] + a.g[i] * b.g[j] + a.g[j] * b.g[i] + a.h[n] * b.v;
        return r;
    }
    //rho, its gradient and its Hessian (row-major 3x3) at (x, y, z). Each shell is f = g(r) Y(u) with g = r^l R(r^2):
    //grad f = g' u Y + g grad Y, H f = (g'' u u^T + g'/r (1 - u u^T)) Y + g' (u grad Y^T + grad Y u^T) + g H Y,
    //where grad Y and H Y with respect to the point come from the hyperdual harmonic seeded with u(d)
    AUX_HD inline double at_hess(const double x, const double y, const double z, const int n_at,
        const double* cx, const double* cy, const double* cz, const double* r2_max,
        const int* sh_start, const int* sh_l, const int* pr_start, const int* coef_off,
        const double* pr_exp, const double* pr_norm, const double* coefs, double& gx, double& gy, double& gz, double* H)
    {
        double dens = 0.0, grad[3] = { 0.0, 0.0, 0.0 };
        for (int i = 0; i < 9; i++) H[i] = 0.0;
        for (int a = 0; a < n_at; a++) {
            const double d[3] = { x - cx[a], y - cy[a], z - cz[a] }, r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
            if (r2 > r2_max[a]) continue;
            const double r = sqrt(r2), u[3] = { d[0] / r, d[1] / r, d[2] / r };
            hyperdual hu[3];
            for (int i = 0; i < 3; i++) {
                hu[i].v = u[i];
                for (int j = 0; j < 3; j++) hu[i].g[j] = ((i == j ? 1.0 : 0.0) - u[i] * u[j]) / r;
                for (int j = 0, n = 0; j < 3; j++)
                    for (int k = j; k < 3; k++, n++)
                        hu[i].h[n] = (3 * u[i] * u[j] * u[k] - (i == j ? u[k] : 0.0) - (i == k ? u[j] : 0.0) - (j == k ? u[i] : 0.0)) / r2;
            }
            for (int s = sh_start[a]; s < sh_start[a + 1]; s++) {
                const int l = sh_l[s];
                double R = 0.0, Ra = 0.0, Raa = 0.0, rl = 1.0;
                for (int p = pr_start[s]; p < pr_start[s + 1]; p++) {
                    const double e = exp(-pr_exp[p] * r2) * pr_norm[p];
                    R += e, Ra += pr_exp[p] * e, Raa += pr_exp[p] * pr_exp[p] * e;
                }
                for (int i = 0; i < l; i++) rl *= r;
                const double g = rl * R;
                if (std::abs(g) < 1E-10) continue;
                const double g1 = l * rl / r * R - 2.0 * rl * r * Ra;
                const double g2 = l * (l - 1) * rl / r2 * R - 2.0 * (2 * l + 1) * rl * Ra + 4.0 * rl * r2 * Raa;
                const hyperdual Y = constants::spherical_harmonic(l, hu[0], hu[1], hu[2], coefs + coef_off[s]);
                dens += g * Y.v;
                for (int i = 0; i < 3; i++) grad[i] += g1 * u[i] * Y.v + g * Y.g[i];
                for (int i = 0, n = 0; i < 3; i++)
                    for (int j = i; j < 3; j++, n++) {
                        const double hg = g2 * u[i] * u[j] + g1 / r * ((i == j ? 1.0 : 0.0) - u[i] * u[j]);
                        const double hf = hg * Y.v + g1 * (u[i] * Y.g[j] + Y.g[i] * u[j]) + g * Y.h[n];
                        H[3 * i + j] += hf;
                        if (i != j) H[3 * j + i] += hf;
                    }
            }
        }
        gx = grad[0], gy = grad[1], gz = grad[2];
        return dens;
    }
    //PC07 iso-orbital indicator alpha = tau_P / tau_TF (Perdew and Constantin, Phys. Rev. B 75, 155109 (2007)) with the
    //a, b of r2SCAN-L, from the reduced gradient p = s^2 and the reduced Laplacian q: z = GE4M - F_W is the part of the
    //fourth-order gradient expansion beyond von Weizsaecker, f_ab(z) switches it off where it would go negative
    AUX_HD inline double pc07_alpha(const double p, const double q)
    {
        const double a = 1.784720, b = 0.258304;
        const double D = 8 * q * q / 81 - p * q / 9 + 8 * p * p / 243, fW = 5 * p / 3, GE4 = 1 + 5 * p / 27 + 20 * q / 9 + D;
        const double z = GE4 / sqrt(1 + D * D / ((1 + fW) * (1 + fW))) - fW;
        if (z >= 0.975 * a) return z;
        if (z <= 0.025 * a) return 0.0;
        return z * exp(-a * b / z) * pow(1 + exp(-a / (a - z)), b) / pow(exp(-a / z) + exp(-a / (a - z)), b);
    }
    //ELI-D from the density alone: WFN::computeELI's 0.5 rho (48 / (rho tau - |grad rho|^2 / 4))^(3/8) is, with its
    //tau = 2 tau_conv, 0.5 rho (24 / (rho tau_P))^(3/8) in the Pauli kinetic energy density tau_P = tau - tau_W, and
    //tau_P = alpha_PC07 tau_TF is the deorbitalised estimate. Zero where the density is negligible.
    //ponytail: alpha floored at 0.05 - PC07 switches to 0 in single-orbital-like tails where the exact ELI is merely high,
    //and the floor caps the estimate near the exact maximum (epoxide: 3.4 vs 4.7, valence corr 0.66); tune if a case needs it
    AUX_HD inline double eli_from_density(const double rho, const double g2, const double lap)
    {
        if (rho < 1E-10) return 0.0;
        const double kf2 = pow(3 * 3.1415926535897932384626433832795028 * 3.1415926535897932384626433832795028 * rho, 2.0 / 3.0);
        double alpha = pc07_alpha(g2 / (4 * kf2 * rho * rho), lap / (4 * kf2 * rho));
        if (alpha < 0.05) alpha = 0.05;
        return 0.5 * rho * pow(24.0 / (rho * alpha * 0.3 * kf2 * rho), 0.375);
    }
    //gamma(l+3/2, x) / x^(l+3/2), the lower incomplete gamma without its leading power so the
    //potential below has no r^-(l+1) to cancel: the series below a+1, the erf recurrence above
    AUX_HD inline double lower_gamma_scaled(const int l, const double x)
    {
        const double a = l + 1.5;
        if (x < a + 1.0) {
            double term = 1.0 / a, sum = term;
            for (int k = 1; k < 500 && term > sum * 1E-17; k++) term *= x / (a + k), sum += term;
            return exp(-x) * sum;
        }
        const double ex = exp(-x);
        double g = 1.7724538509055160273 * erf(sqrt(x)), xa = sqrt(x), aa = 0.5;
        for (int i = 0; i <= l; i++) g = aa * g - xa * ex, xa *= x, aa += 1.0;
        return g / xa;
    }
    //Electrostatic potential of the fitted density plus the nuclei at (x, y, z), in Hartree/e like
    //WFN::computeESP. Each shell n exp(-a r^2) r^l Y(u) has the closed-form potential
    //4 pi / (2l+1) Y(u) sum_p n_p [ r^(l+2) G(l+3/2, a_p r^2) / 2 + r^l exp(-a_p r^2) / (2 a_p) ]
    //with G the scaled incomplete gamma above, so no r2_max cutoff: the far field is the multipole tail.
    AUX_HD inline double esp_at(const double x, const double y, const double z, const int n_at,
        const double* cx, const double* cy, const double* cz, const int* Z,
        const int* sh_start, const int* sh_l, const int* pr_start, const int* coef_off,
        const double* pr_exp, const double* pr_norm, const double* coefs)
    {
        double esp = 0.0;
        for (int a = 0; a < n_at; a++) {
            const double dx = x - cx[a], dy = y - cy[a], dz = z - cz[a], r2 = dx * dx + dy * dy + dz * dz;
            const double r = sqrt(r2), ux = dx / r, uy = dy / r, uz = dz / r;
            esp += Z[a] / r;
            for (int s = sh_start[a]; s < sh_start[a + 1]; s++) {
                const int l = sh_l[s];
                double rl = 1.0;
                for (int i = 0; i < l; i++) rl *= r;
                double radial = 0.0;
                for (int p = pr_start[s]; p < pr_start[s + 1]; p++)
                    radial += pr_norm[p] * (0.5 * rl * r2 * lower_gamma_scaled(l, pr_exp[p] * r2) + 0.5 * rl * exp(-pr_exp[p] * r2) / pr_exp[p]);
                esp -= 4.0 * 3.1415926535897932384626433832795028 / (2 * l + 1) * radial * constants::spherical_harmonic(l, ux, uy, uz, coefs + coef_off[s]);
            }
        }
        return esp;
    }
}
