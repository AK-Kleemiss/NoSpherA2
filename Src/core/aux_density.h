#pragma once
#include <cmath>

//The fitted density at one point from a flattened auxiliary basis. One text for the CPU loop
//in SALTED_utilities.cpp and the kernel in aux_density_gpu.cu, so the two cannot drift: nvcc
//and hipcc see the host-device attribute, every other compiler sees plain inline functions.
//Nothing from the project is included here on purpose; the .cu must stay parseable by both.
#if defined(__CUDACC__) || defined(__HIPCC__)
#define AUX_HD __host__ __device__
#else
#define AUX_HD
#endif

namespace aux_density
{
    constexpr double PI = 3.1415926535897932384626433832795028;
    //Value with its three derivatives, so the harmonic text below differentiates itself:
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
    //The real spherical harmonics of one l contracted with their coefficients, the switch of
    //constants::spherical_harmonic(l, d, coefs) up to l = 8 with x, y, z already normalised
    template <class T>
    AUX_HD inline T harmonic(const int l, const T x, const T y, const T z, const double* coefs)
    {
        switch (l)
        {
        case 0:
            return sqrt(1.0 / (4 * PI)) * coefs[0];
        case 1:
            return sqrt(3.0 / (4 * PI)) * (coefs[0] * y + coefs[1] * z + coefs[2] * x);
        case 2:
            return sqrt(15.0 / (4 * PI)) * (y * x * coefs[0] + y * z * coefs[1] + x * z * coefs[3]) + sqrt(5.0 / (16.0 * PI)) * (3 * z * z - 1.0) * coefs[2] + sqrt(15.0 / (16.0 * PI)) * (x * x - y * y) * coefs[4];
        case 3:
        {
            const T y2 = y * y, x2 = x * x, z2 = z * z;
            return sqrt(35.0 / (32.0 * PI)) * (y * (3 * x2 - y2) * coefs[0] + x * (x2 - 3 * y2) * coefs[6]) +
                sqrt(105.0 / (4 * PI)) * x * y * z * coefs[1] +
                sqrt(21.0 / (32.0 * PI)) * (y * (5 * z2 - 1.0) * coefs[2] + x * (5 * z2 - 1.0) * coefs[4]) +
                sqrt(7.0 / (16.0 * PI)) * (5 * z2 * z - 3 * z) * coefs[3] +
                sqrt(105.0 / (16.0 * PI)) * ((x2 - y2) * z) * coefs[5];
        }
        case 4:
        {
            const T x2 = x * x, y2 = y * y, z2 = z * z;
            return sqrt(315.0 / (16.0 * PI)) * x * y * (x2 - y2) * coefs[0] +
                sqrt(315.0 / (32.0 * PI)) * (y * (3 * x2 - y2) * z * coefs[1] + x * (x2 - 3 * y2) * z * coefs[7]) +
                sqrt(45.0 / (16.0 * PI)) * x * y * (7 * z2 - 1.0) * coefs[2] +
                sqrt(45.0 / (32.0 * PI)) * (y * (7 * z2 * z - 3 * z) * coefs[3] + x * (7 * z2 * z - 3 * z) * coefs[5]) +
                sqrt(9.0 / (256.0 * PI)) * (35 * z2 * z2 - 30 * z2 + 3.0) * coefs[4] +
                sqrt(45.0 / (64.0 * PI)) * (x2 - y2) * (7 * z2 - 1.0) * coefs[6] +
                sqrt(315.0 / (256.0 * PI)) * ((x2 * (x2 - 3 * y2)) - (y2 * (3 * x2 - y2))) * coefs[8];
        }
        case 5:
        {
            const T x2 = x * x, y2 = y * y, z2 = z * z;
            return sqrt(693.0 / (2048.0 * PI)) * (2 * y2 * y2 * y - 20 * x2 * y2 * y + 10 * y * x2 * x2) * coefs[0] +
                -sqrt(3465.0 / (256.0 * PI)) * z * ((4 * x * y2 * y - 4 * x2 * x * y) * coefs[1] - (x2 * x2 - 6 * x2 * y2 + y2 * y2) * coefs[9]) +
                sqrt(385.0 / (512.0 * PI)) * (9 * z2 - 1.0) * (y * (3 * x2 - y2) * coefs[2] + x * (x2 - 3 * y2) * coefs[8]) +
                sqrt(1155.0 / (64.0 * PI)) * (3 * z2 * z - z) * (2 * x * y * coefs[3] + (x2 - y2) * coefs[7]) +
                sqrt(165.0 / (256.0 * PI)) * (21 * z2 * z2 - 14 * z2 + 1.0) * (y * coefs[4] + x * coefs[6]) +
                sqrt(11.0 / (256.0 * PI)) * (63 * z2 * z2 * z - 70 * z2 * z + 15 * z) * coefs[5] +
                sqrt(693.0 / (2048.0 * PI)) * (2 * x2 * x2 * x - 20 * x2 * x * y2 + 10 * x * y2 * y2) * coefs[10];
        }
        case 6:
        {
            const T x2 = x * x, y2 = y * y, z2 = z * z;
            const T x4 = x2 * x2, y4 = y2 * y2, z4 = z2 * z2;
            return (1.0 / 32.0) * sqrt(13.0 / PI) * (z2 * (z2 * (231.0 * z2 - 315.0) + 105.0) - 5.0) * coefs[6] +
                (1.0 / 16.0) * sqrt(273.0 / (PI)) * (33 * z4 - 30 * z2 + 5) * (x * z * coefs[7] + y * z * coefs[5]) +
                (1.0 / 32.0) * sqrt(1365.0 / (2 * PI)) * (33 * z4 - 18 * z2 + 1) * ((x2 - y2) * coefs[8] + 2 * x * y * coefs[4]) +
                (1.0 / 16.0) * sqrt(1365.0 / (2 * PI)) * z * (11 * z2 - 3) * (x * (x2 - 3 * y2) * coefs[9] + y * (3 * x2 - y2) * coefs[3]) +
                (3.0 / 32.0) * sqrt(91.0 / (PI)) * (11 * z2 - 1) * ((x2 * (x2 - 6 * y2) + y4) * coefs[10] + x * y * 4 * (x2 - y2) * coefs[2]) +
                (3.0 / 16.0) * sqrt(1001.0 / (2 * PI)) * z * ((x4 * x - 10 * x2 * x * y2 + 5 * x * y4) * coefs[11] + (5 * x4 * y - 10 * x2 * y2 * y + y4 * y) * coefs[1]) +
                (1.0 / 32.0) * sqrt(3003.0 / (2 * PI)) * (-y4 * y2 + 15 * y4 * x2 - 15 * y2 * x4 + x4 * x2) * coefs[12] +
                (1.0 / 32.0) * sqrt(3003.0 / (2 * PI)) * (6 * x4 * x * y - 20 * x2 * x * y2 * y + 6 * x * y4 * y) * coefs[0];
        }
        case 7:
        {
            const T x2 = x * x, y2 = y * y, z2 = z * z;
            const T x4 = x2 * x2, y4 = y2 * y2;
            return (1.0 / 32.0) * sqrt(15.0 / PI) * z * (z2 * (z2 * (429.0 * z2 - 693.0) + 315.0) - 35.0) * coefs[7] +
                (1.0 / 64.0) * sqrt(105.0 / PI) * (z2 * (z2 * (429.0 * z2 - 495.0) + 135.0) - 5.0) * (x * coefs[8] + y * coefs[6]) +
                (3.0 / 32.0) * sqrt(35.0 / (2 * PI)) * z * (z2 * (143.0 * z2 - 110.0) + 15.0) * ((x2 - y2) * coefs[9] + 2 * x * y * coefs[5]) +
                (3.0 / 64.0) * sqrt(35.0 / PI) * (z2 * (143.0 * z2 - 66.0) + 3.0) * (x * (x2 - 3 * y2) * coefs[10] + y * (3 * x2 - y2) * coefs[4]) +
                (3.0 / 32.0) * sqrt(385.0 / PI) * z * (13.0 * z2 - 3.0) * ((x2 * (x2 - 6 * y2) + y4) * coefs[11] + x * y * 4 * (x2 - y2) * coefs[3]) +
                (3.0 / 64.0) * sqrt(385.0 / PI) * (13.0 * z2 - 1.0) * (x * (x4 + 5 * y2 * (-2 * x2 + y2)) * coefs[12] + y * (5 * x2 * (x2 - 2 * y2) + y4) * coefs[2]) +
                (3.0 / 32.0) * sqrt(5005.0 / (2 * PI)) * z * ((y2 * (y2 * (-y2 + 15 * x2) - 15 * x4) + x4 * x2) * coefs[13] + x * y * (x2 * (6 * x2 - 20 * y2) + 6 * y4) * coefs[1]) +
                (3.0 / 64.0) * sqrt(715.0 / PI) * x * (x4 * x2 + 7 * y2 * (-3 * x4 + y2 * (5 * x2 - y2))) * coefs[14] +
                (3.0 / 64.0) * sqrt(715.0 / PI) * y * (x2 * (x2 * (7 * x2 - 35 * y2) + 21 * y4) - y4 * y2) * coefs[0];
        }
        case 8:
        {
            const T x2 = x * x, y2 = y * y, z2 = z * z;
            const T x4 = x2 * x2, y4 = y2 * y2, z4 = z2 * z2;
            return (1.0 / 256.0) * sqrt(17.0 / PI) * (z2 * (z2 * (z2 * (6435.0 * z2 - 12012.0) + 6930.0) - 1260.0) + 35.0) * coefs[8] +
                (3.0 / 64.0) * sqrt(17.0 / PI) * z * (715.0 * z4 * z2 - 1001.0 * z4 + 385.0 * z2 - 35.0) * (x * coefs[9] + y * coefs[7]) +
                (3.0 / 64.0) * sqrt(595.0 / (2 * PI)) * (143 * z4 * z2 - 143 * z4 + 33 * z2 - 1) * ((x2 - y2) * coefs[10] + 2 * x * y * coefs[6]) +
                (1.0 / 64.0) * sqrt(19635.0 / PI) * z * (39 * z4 - 26 * z2 + 3) * ((x2 * x - 3 * x * y2) * coefs[11] + (3 * x2 * y - y2 * y) * coefs[5]) +
                (3.0 / 128.0) * sqrt(1309.0 / PI) * (65 * z4 - 26 * z2 + 1) * ((x2 * (x2 - 6 * y2) + y4) * coefs[12] + x * y * 4 * (x2 - y2) * coefs[4]) +
                (3.0 / 64.0) * sqrt(17017.0 / PI) * z * (5 * z2 - 1) * (x * (x4 + 5 * y2 * (-2 * x2 + y2)) * coefs[13] + y * (5 * x2 * (x2 - 2 * y2) + y4) * coefs[3]) +
                (1.0 / 64.0) * sqrt(7293.0 / (2 * PI)) * (15 * z2 - 1) * ((y2 * (y2 * (-y2 + 15 * x2) - 15 * x4) + x4 * x2) * coefs[14] + x * y * (x2 * (6 * x2 - 20 * y2) + 6 * y4) * coefs[2]) +
                (3.0 / 64.0) * sqrt(12155.0 / PI) * z * (x * (x4 * x2 + 7 * y2 * (-3 * x4 + y2 * (5 * x2 - y2))) * coefs[15] + y * (x2 * (x2 * (7 * x2 - 35 * y2) + 21 * y4) - y4 * y2) * coefs[1]) +
                (3.0 / 256.0) * sqrt(12155.0 / PI) * (y2 * (y2 * (y4 - 28 * x2 * y2 + 70 * x4) - 28 * x4 * x2) + x4 * x4) * coefs[16] +
                (3.0 / 256.0) * sqrt(12155.0 / PI) * y * x * (y2 * (y2 * (-8 * y2 + 56 * x2) - 56 * x4) + 8 * x4 * x2) * coefs[0];
        }
        }
        return T(0.0);
    }
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
                if (radial < 1E-10) continue;
                dens += radial * harmonic(sh_l[s], ux, uy, uz, coefs + coef_off[s]);
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
                if (radial < 1E-10) continue;
                const dual Y = harmonic(l, dual(ux, 1.0, 0.0, 0.0), dual(uy, 0.0, 1.0, 0.0), dual(uz, 0.0, 0.0, 1.0), coefs + coef_off[s]);
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
