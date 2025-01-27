#include "int_g2e.h"
#include "rys_roots.h"

void CINTcart_comp(int* nx, int* ny, int* nz, const int lmax)
{
    int inc = 0;
    int lx, ly, lz;

    for (lx = lmax; lx >= 0; lx--) {
        for (ly = lmax - lx; ly >= 0; ly--) {
            lz = lmax - lx - ly;
            nx[inc] = lx;
            ny[inc] = ly;
            nz[inc] = lz;
            inc++;
        }
    }
}

void CINTg2e_index_xyz(int* idx, const CINTEnvVars* envs)
{
    const int i_l = envs->i_l;
    const int j_l = envs->j_l;
    const int k_l = envs->k_l;
    const int l_l = envs->l_l;
    const int nfi = envs->nfi;
    const int nfj = envs->nfj;
    const int nfk = envs->nfk;
    const int nfl = envs->nfl;
    const int di = envs->g_stride_i;
    const int dk = envs->g_stride_k;
    const int dl = envs->g_stride_l;
    const int dj = envs->g_stride_j;
    int i, j, k, l, n;
    int ofx, ofkx, oflx;
    int ofy, ofky, ofly;
    int ofz, ofkz, oflz;
    int i_nx[CART_MAX], i_ny[CART_MAX], i_nz[CART_MAX];
    int j_nx[CART_MAX], j_ny[CART_MAX], j_nz[CART_MAX];
    int k_nx[CART_MAX], k_ny[CART_MAX], k_nz[CART_MAX];
    int l_nx[CART_MAX], l_ny[CART_MAX], l_nz[CART_MAX];

    CINTcart_comp(i_nx, i_ny, i_nz, i_l);
    CINTcart_comp(j_nx, j_ny, j_nz, j_l);
    CINTcart_comp(k_nx, k_ny, k_nz, k_l);
    CINTcart_comp(l_nx, l_ny, l_nz, l_l);

    ofx = 0;
    ofy = envs->g_size;
    ofz = envs->g_size * 2;
    n = 0;
    for (j = 0; j < nfj; j++) {
        for (l = 0; l < nfl; l++) {
            oflx = ofx + dj * j_nx[j] + dl * l_nx[l];
            ofly = ofy + dj * j_ny[j] + dl * l_ny[l];
            oflz = ofz + dj * j_nz[j] + dl * l_nz[l];
            for (k = 0; k < nfk; k++) {
                ofkx = oflx + dk * k_nx[k];
                ofky = ofly + dk * k_ny[k];
                ofkz = oflz + dk * k_nz[k];
                switch (i_l) {
                case 0:
                    idx[n + 0] = ofkx;
                    idx[n + 1] = ofky;
                    idx[n + 2] = ofkz;
                    n += 3;
                    break;
                case 1:
                    idx[n + 0] = ofkx + di;
                    idx[n + 1] = ofky;
                    idx[n + 2] = ofkz;
                    idx[n + 3] = ofkx;
                    idx[n + 4] = ofky + di;
                    idx[n + 5] = ofkz;
                    idx[n + 6] = ofkx;
                    idx[n + 7] = ofky;
                    idx[n + 8] = ofkz + di;
                    n += 9;
                    break;
                case 2:
                    idx[n + 0] = ofkx + di * 2;
                    idx[n + 1] = ofky;
                    idx[n + 2] = ofkz;
                    idx[n + 3] = ofkx + di;
                    idx[n + 4] = ofky + di;
                    idx[n + 5] = ofkz;
                    idx[n + 6] = ofkx + di;
                    idx[n + 7] = ofky;
                    idx[n + 8] = ofkz + di;
                    idx[n + 9] = ofkx;
                    idx[n + 10] = ofky + di * 2;
                    idx[n + 11] = ofkz;
                    idx[n + 12] = ofkx;
                    idx[n + 13] = ofky + di;
                    idx[n + 14] = ofkz + di;
                    idx[n + 15] = ofkx;
                    idx[n + 16] = ofky;
                    idx[n + 17] = ofkz + di * 2;
                    n += 18;
                    break;
                default:
                    for (i = 0; i < nfi; i++) {
                        idx[n + 0] = ofkx + di * i_nx[i]; //(:,ix,kx,lx,jx,1)
                        idx[n + 1] = ofky + di * i_ny[i]; //(:,iy,ky,ly,jy,2)
                        idx[n + 2] = ofkz + di * i_nz[i]; //(:,iz,kz,lz,jz,3)
                        n += 3;
                    } // i
                }
            } // k
        } // l
    } // j
}


void CINTg1e_index_xyz(int* idx, const CINTEnvVars* envs)
{
    const int i_l = envs->i_l;
    const int j_l = envs->j_l;
    const int nfi = envs->nfi;
    const int nfj = envs->nfj;
    const int di = envs->g_stride_i;
    const int dj = envs->g_stride_j;
    int i, j, n;
    int ofx, ofjx;
    int ofy, ofjy;
    int ofz, ofjz;
    int i_nx[CART_MAX], i_ny[CART_MAX], i_nz[CART_MAX];
    int j_nx[CART_MAX], j_ny[CART_MAX], j_nz[CART_MAX];

    CINTcart_comp(i_nx, i_ny, i_nz, i_l);
    CINTcart_comp(j_nx, j_ny, j_nz, j_l);

    ofx = 0;
    ofy = envs->g_size;
    ofz = envs->g_size * 2;
    n = 0;
    for (j = 0; j < nfj; j++) {
        ofjx = ofx + dj * j_nx[j];
        ofjy = ofy + dj * j_ny[j];
        ofjz = ofz + dj * j_nz[j];
        for (i = 0; i < nfi; i++) {
            idx[n + 0] = ofjx + di * i_nx[i];
            idx[n + 1] = ofjy + di * i_ny[i];
            idx[n + 2] = ofjz + di * i_nz[i];
            n += 3;
        }
    }
}


/*
 * g[i,k,l,j] = < ik | lj > = ( i j | k l )
 */
int CINTg0_2e(double* g, double* rij, double* rkl, double cutoff, CINTEnvVars* envs)
{
    int irys;
    int nroots = envs->nrys_roots;
    double aij = envs->ai[0] + envs->aj[0];
    double akl = envs->ak[0] + envs->al[0];
    double a0, a1, fac1, x;
    double u[MXRYSROOTS];
    double* w = g + envs->g_size * 2; // ~ gz
    double xij_kl = rij[0] - rkl[0];
    double yij_kl = rij[1] - rkl[1];
    double zij_kl = rij[2] - rkl[2];
    double rr = xij_kl * xij_kl + yij_kl * yij_kl + zij_kl * zij_kl;

    a1 = aij * akl;
    a0 = a1 / (aij + akl);
    fac1 = sqrt(a0 / (a1 * a1 * a1)) * envs->fac[0];
    x = a0 * rr;
    const double omega = envs->env[PTR_RANGE_OMEGA];
    double theta = 0;
    if (omega == 0.) {
        CINTrys_roots(nroots, x, u, w);
    }
    else if (omega < 0.) {
        // short-range part of range-separated Coulomb
        theta = omega * omega / (omega * omega + a0);
        // very small erfc() leads to ~0 weights. They can cause
        // numerical issue in sr_rys_roots.
        if (theta * x > cutoff || theta * x > EXPCUTOFF_SR) {
            return 0;
        }
        int rorder = envs->rys_order;
        if (rorder == nroots) {
            CINTsr_rys_roots(nroots, x, sqrt(theta), u, w);
        }
        else {
            double sqrt_theta = -sqrt(theta);
            CINTrys_roots(rorder, x, u, w);
            CINTrys_roots(rorder, theta * x, u + rorder, w + rorder);
            if (envs->g_size == 2) {
                g[0] = g[1] = g[2] = g[3] = 1;
                g[4] *= fac1;
                g[5] *= fac1 * sqrt_theta;
                return 1;
            }
            for (irys = rorder; irys < nroots; irys++) {
                double ut = u[irys] * theta;
                u[irys] = ut / (u[irys] + 1. - ut);
                w[irys] *= sqrt_theta;
            }
        }
    }
    else { // omega > 0.
        // long-range part of range-separated Coulomb
        theta = omega * omega / (omega * omega + a0);
        x *= theta;
        fac1 *= sqrt(theta);
        CINTrys_roots(nroots, x, u, w);
        /* u[:] = tau^2 / (1 - tau^2)
         * omega^2u^2 = a0 * tau^2 / (theta^-1 - tau^2)
         * transform u[:] to theta^-1 tau^2 / (theta^-1 - tau^2)
         * so the rest code can be reused.
         */
        for (irys = 0; irys < nroots; irys++) {
            double ut = u[irys] * theta;
            u[irys] = ut / (u[irys] + 1. - ut);
        }
    }
    if (envs->g_size == 1) {
        g[0] = g[1] = 1;
        g[2] *= fac1;
        return 1;
    }

    double u2, tmp1, tmp2, tmp3, tmp4, tmp5;
    double rijrx = rij[0] - envs->rx_in_rijrx[0];
    double rijry = rij[1] - envs->rx_in_rijrx[1];
    double rijrz = rij[2] - envs->rx_in_rijrx[2];
    double rklrx = rkl[0] - envs->rx_in_rklrx[0];
    double rklry = rkl[1] - envs->rx_in_rklrx[1];
    double rklrz = rkl[2] - envs->rx_in_rklrx[2];
    Rys2eT bc;
    double* b00 = bc.b00;
    double* b10 = bc.b10;
    double* b01 = bc.b01;
    double* c00x = bc.c00x;
    double* c00y = bc.c00y;
    double* c00z = bc.c00z;
    double* c0px = bc.c0px;
    double* c0py = bc.c0py;
    double* c0pz = bc.c0pz;

    for (irys = 0; irys < nroots; irys++) {
        /*
         *u(irys) = t2/(1-t2)
         *t2 = u(irys)/(1+u(irys))
         *u2 = aij*akl/(aij+akl)*t2/(1-t2)
         */
        u2 = a0 * u[irys];
        tmp4 = .5 / (u2 * (aij + akl) + a1);
        tmp5 = u2 * tmp4;
        tmp1 = 2. * tmp5;
        tmp2 = tmp1 * akl;
        tmp3 = tmp1 * aij;
        b00[irys] = tmp5;
        b10[irys] = tmp5 + tmp4 * akl;
        b01[irys] = tmp5 + tmp4 * aij;
        c00x[irys] = rijrx - tmp2 * xij_kl;
        c00y[irys] = rijry - tmp2 * yij_kl;
        c00z[irys] = rijrz - tmp2 * zij_kl;
        c0px[irys] = rklrx + tmp3 * xij_kl;
        c0py[irys] = rklry + tmp3 * yij_kl;
        c0pz[irys] = rklrz + tmp3 * zij_kl;
        w[irys] *= fac1;
    }

    (*envs->f_g0_2d4d)(g, &bc, envs);

    return 1;
}


/*
 * g(nroots,0:nmax,0:mmax)
 */
void CINTg0_2e_2d(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    const int nroots = envs->nrys_roots;
    const int nmax = envs->li_ceil + envs->lj_ceil;
    const int mmax = envs->lk_ceil + envs->ll_ceil;
    const int dm = envs->g2d_klmax;
    const int dn = envs->g2d_ijmax;
    int i, j, m, n, off;
    DEF_GXYZ(double, g, gx, gy, gz);
    double* p0x, * p0y, * p0z;
    double* p1x, * p1y, * p1z;
    double nb1, mb0;

    for (i = 0; i < nroots; i++) {
        gx[i] = 1;
        gy[i] = 1;
        //gz[i] = w[i];
    }

    double s0x, s1x, s2x, t0x, t1x;
    double s0y, s1y, s2y, t0y, t1y;
    double s0z, s1z, s2z, t0z, t1z;
    double c00x, c00y, c00z, c0px, c0py, c0pz, b10, b01, b00;
    for (i = 0; i < nroots; i++) {
        c00x = bc->c00x[i];
        c00y = bc->c00y[i];
        c00z = bc->c00z[i];
        c0px = bc->c0px[i];
        c0py = bc->c0py[i];
        c0pz = bc->c0pz[i];
        b10 = bc->b10[i];
        b01 = bc->b01[i];
        b00 = bc->b00[i];
        if (nmax > 0) {
            // gx(irys,0,1) = c00(irys) * gx(irys,0,0)
            // gx(irys,0,n+1) = c00(irys)*gx(irys,0,n)
            // + n*b10(irys)*gx(irys,0,n-1)
            s0x = gx[i];
            s0y = gy[i];
            s0z = gz[i];
            s1x = c00x * s0x;
            s1y = c00y * s0y;
            s1z = c00z * s0z;
            gx[i + dn] = s1x;
            gy[i + dn] = s1y;
            gz[i + dn] = s1z;
            for (n = 1; n < nmax; ++n) {
                s2x = c00x * s1x + n * b10 * s0x;
                s2y = c00y * s1y + n * b10 * s0y;
                s2z = c00z * s1z + n * b10 * s0z;
                gx[i + (n + 1) * dn] = s2x;
                gy[i + (n + 1) * dn] = s2y;
                gz[i + (n + 1) * dn] = s2z;
                s0x = s1x;
                s0y = s1y;
                s0z = s1z;
                s1x = s2x;
                s1y = s2y;
                s1z = s2z;
            }
        }

        if (mmax > 0) {
            // gx(irys,1,0) = c0p(irys) * gx(irys,0,0)
            // gx(irys,m+1,0) = c0p(irys)*gx(irys,m,0)
            // + m*b01(irys)*gx(irys,m-1,0)
            s0x = gx[i];
            s0y = gy[i];
            s0z = gz[i];
            s1x = c0px * s0x;
            s1y = c0py * s0y;
            s1z = c0pz * s0z;
            gx[i + dm] = s1x;
            gy[i + dm] = s1y;
            gz[i + dm] = s1z;
            for (m = 1; m < mmax; ++m) {
                s2x = c0px * s1x + m * b01 * s0x;
                s2y = c0py * s1y + m * b01 * s0y;
                s2z = c0pz * s1z + m * b01 * s0z;
                gx[i + (m + 1) * dm] = s2x;
                gy[i + (m + 1) * dm] = s2y;
                gz[i + (m + 1) * dm] = s2z;
                s0x = s1x;
                s0y = s1y;
                s0z = s1z;
                s1x = s2x;
                s1y = s2y;
                s1z = s2z;
            }

            if (nmax > 0) {
                // gx(irys,1,1) = c0p(irys)*gx(irys,0,1) + b00(irys)*gx(irys,0,0)
                // gx(irys,m+1,1) = c0p(irys)*gx(irys,m,1)
                // + m*b01(irys)*gx(irys,m-1,1)
                // + b00(irys)*gx(irys,m,0)
                s0x = gx[i + dn];
                s0y = gy[i + dn];
                s0z = gz[i + dn];
                s1x = c0px * s0x + b00 * gx[i];
                s1y = c0py * s0y + b00 * gy[i];
                s1z = c0pz * s0z + b00 * gz[i];
                gx[i + dn + dm] = s1x;
                gy[i + dn + dm] = s1y;
                gz[i + dn + dm] = s1z;
                for (m = 1; m < mmax; ++m) {
                    s2x = c0px * s1x + m * b01 * s0x + b00 * gx[i + m * dm];
                    s2y = c0py * s1y + m * b01 * s0y + b00 * gy[i + m * dm];
                    s2z = c0pz * s1z + m * b01 * s0z + b00 * gz[i + m * dm];
                    gx[i + dn + (m + 1) * dm] = s2x;
                    gy[i + dn + (m + 1) * dm] = s2y;
                    gz[i + dn + (m + 1) * dm] = s2z;
                    s0x = s1x;
                    s0y = s1y;
                    s0z = s1z;
                    s1x = s2x;
                    s1y = s2y;
                    s1z = s2z;
                }
            }
        }

        // gx(irys,m,n+1) = c00(irys)*gx(irys,m,n)
        // + n*b10(irys)*gx(irys,m,n-1)
        // + m*b00(irys)*gx(irys,m-1,n)
        for (m = 1; m <= mmax; ++m) {
            off = m * dm;
            j = off + i;
            s0x = gx[j];
            s0y = gy[j];
            s0z = gz[j];
            s1x = gx[j + dn];
            s1y = gy[j + dn];
            s1z = gz[j + dn];
            for (n = 1; n < nmax; ++n) {
                s2x = c00x * s1x + n * b10 * s0x + m * b00 * gx[j + n * dn - dm];
                s2y = c00y * s1y + n * b10 * s0y + m * b00 * gy[j + n * dn - dm];
                s2z = c00z * s1z + n * b10 * s0z + m * b00 * gz[j + n * dn - dm];
                gx[j + (n + 1) * dn] = s2x;
                gy[j + (n + 1) * dn] = s2y;
                gz[j + (n + 1) * dn] = s2z;
                s0x = s1x;
                s0y = s1y;
                s0z = s1z;
                s1x = s2x;
                s1y = s2y;
                s1z = s2z;
            }
        }
    }
}
/*
 * g0[i,k,l,j] = < ik | lj > = ( i j | k l )
 */
 /* 2d is based on l,j */
void CINTg0_lj2d_4d(double* g, CINTEnvVars* envs)
{
    int li = envs->li_ceil;
    int lk = envs->lk_ceil;
    if (li == 0 && lk == 0) {
        return;
    }
    int nmax = envs->li_ceil + envs->lj_ceil;
    int mmax = envs->lk_ceil + envs->ll_ceil;
    //int ll = envs->ll_ceil;
    int lj = envs->lj_ceil;
    int nroots = envs->nrys_roots;
    int i, j, k, l, ptr, n;
    int di = envs->g_stride_i;
    int dk = envs->g_stride_k;
    int dl = envs->g_stride_l;
    int dj = envs->g_stride_j;
    double* rirj = envs->rirj;
    double* rkrl = envs->rkrl;
    DEF_GXYZ(double, g, gx, gy, gz);
    double* p1x, * p1y, * p1z, * p2x, * p2y, * p2z;
    double rx, ry, rz;

    // g(i,...,j) = rirj * g(i-1,...,j) +  g(i-1,...,j+1)
    rx = rirj[0];
    ry = rirj[1];
    rz = rirj[2];
    p1x = gx - di;
    p1y = gy - di;
    p1z = gz - di;
    p2x = gx - di + dj;
    p2y = gy - di + dj;
    p2z = gz - di + dj;
    for (i = 1; i <= li; i++) {
        for (j = 0; j <= nmax - i; j++) {
            for (l = 0; l <= mmax; l++) {
                ptr = j * dj + l * dl + i * di;
                for (n = ptr; n < ptr + nroots; n++) {
                    gx[n] = rx * p1x[n] + p2x[n];
                    gy[n] = ry * p1y[n] + p2y[n];
                    gz[n] = rz * p1z[n] + p2z[n];
                }
            }
        }
    }

    // g(...,k,l,..) = rkrl * g(...,k-1,l,..) + g(...,k-1,l+1,..)
    rx = rkrl[0];
    ry = rkrl[1];
    rz = rkrl[2];
    p1x = gx - dk;
    p1y = gy - dk;
    p1z = gz - dk;
    p2x = gx - dk + dl;
    p2y = gy - dk + dl;
    p2z = gz - dk + dl;
    for (j = 0; j <= lj; j++) {
        for (k = 1; k <= lk; k++) {
            for (l = 0; l <= mmax - k; l++) {
                ptr = j * dj + l * dl + k * dk;
                for (n = ptr; n < ptr + dk; n++) {
                    gx[n] = rx * p1x[n] + p2x[n];
                    gy[n] = ry * p1y[n] + p2y[n];
                    gz[n] = rz * p1z[n] + p2z[n];
                }
            }
        }
    }
}


void CINTg0_2e_lj2d4d(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    CINTg0_2e_2d(g, bc, envs);
    CINTg0_lj2d_4d(g, envs);
}

/* 2d is based on i,l */
void CINTg0_il2d_4d(double* g, CINTEnvVars* envs)
{
    int lk = envs->lk_ceil;
    int lj = envs->lj_ceil;
    if (lj == 0 && lk == 0) {
        return;
    }
    int nmax = envs->li_ceil + envs->lj_ceil;
    int mmax = envs->lk_ceil + envs->ll_ceil;
    //int li = envs->li_ceil;
    int ll = envs->ll_ceil;
    int nroots = envs->nrys_roots;
    int i, j, k, l, ptr, n;
    int di = envs->g_stride_i;
    int dk = envs->g_stride_k;
    int dl = envs->g_stride_l;
    int dj = envs->g_stride_j;
    double* rirj = envs->rirj;
    double* rkrl = envs->rkrl;
    DEF_GXYZ(double, g, gx, gy, gz);
    double* p1x, * p1y, * p1z, * p2x, * p2y, * p2z;
    double rx, ry, rz;

    // g(...,k,l,..) = rkrl * g(...,k-1,l,..) + g(...,k-1,l+1,..)
    rx = rkrl[0];
    ry = rkrl[1];
    rz = rkrl[2];
    p1x = gx - dk;
    p1y = gy - dk;
    p1z = gz - dk;
    p2x = gx - dk + dl;
    p2y = gy - dk + dl;
    p2z = gz - dk + dl;
    for (k = 1; k <= lk; k++) {
        for (l = 0; l <= mmax - k; l++) {
            for (i = 0; i <= nmax; i++) {
                ptr = l * dl + k * dk + i * di;
                for (n = ptr; n < ptr + nroots; n++) {
                    gx[n] = rx * p1x[n] + p2x[n];
                    gy[n] = ry * p1y[n] + p2y[n];
                    gz[n] = rz * p1z[n] + p2z[n];
                }
            }
        }
    }

    // g(i,...,j) = rirj * g(i,...,j-1) +  g(i+1,...,j-1)
    rx = rirj[0];
    ry = rirj[1];
    rz = rirj[2];
    p1x = gx - dj;
    p1y = gy - dj;
    p1z = gz - dj;
    p2x = gx - dj + di;
    p2y = gy - dj + di;
    p2z = gz - dj + di;
    for (j = 1; j <= lj; j++) {
        for (l = 0; l <= ll; l++) {
            for (k = 0; k <= lk; k++) {
                ptr = j * dj + l * dl + k * dk;
                for (n = ptr; n < ptr + dk - di * j; n++) {
                    gx[n] = rx * p1x[n] + p2x[n];
                    gy[n] = ry * p1y[n] + p2y[n];
                    gz[n] = rz * p1z[n] + p2z[n];
                }
            }
        }
    }
}


void CINTg0_2e_il2d4d(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    CINTg0_2e_2d(g, bc, envs);
    CINTg0_il2d_4d(g, envs);
}

static inline void _g0_2d4d_0000(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    g[0] = 1;
    g[1] = 1;
    //g[2] = w[0];
}
static inline void _g0_2d4d_0001(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    g[0] = 1;
    g[1] = cpx[0];
    g[2] = 1;
    g[3] = cpy[0];
    //g[4] = w[0];
    g[5] = cpz[0] * g[4];
}
static inline void _g0_2d4d_0002(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = cpx[0];
    g[3] = cpx[1];
    g[4] = cpx[0] * cpx[0] + b01[0];
    g[5] = cpx[1] * cpx[1] + b01[1];
    g[6] = 1;
    g[7] = 1;
    g[8] = cpy[0];
    g[9] = cpy[1];
    g[10] = cpy[0] * cpy[0] + b01[0];
    g[11] = cpy[1] * cpy[1] + b01[1];
    //g[12] = w[0];
    //g[13] = w[1];
    g[14] = cpz[0] * g[12];
    g[15] = cpz[1] * g[13];
    g[16] = cpz[0] * g[14] + b01[0] * g[12];
    g[17] = cpz[1] * g[15] + b01[1] * g[13];
}
static inline void _g0_2d4d_0003(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = cpx[0];
    g[3] = cpx[1];
    g[4] = cpx[0] * cpx[0] + b01[0];
    g[5] = cpx[1] * cpx[1] + b01[1];
    g[6] = cpx[0] * (g[4] + 2 * b01[0]);
    g[7] = cpx[1] * (g[5] + 2 * b01[1]);
    g[8] = 1;
    g[9] = 1;
    g[10] = cpy[0];
    g[11] = cpy[1];
    g[12] = cpy[0] * cpy[0] + b01[0];
    g[13] = cpy[1] * cpy[1] + b01[1];
    g[14] = cpy[0] * (g[12] + 2 * b01[0]);
    g[15] = cpy[1] * (g[13] + 2 * b01[1]);
    //g[16] = w[0];
    //g[17] = w[1];
    g[18] = cpz[0] * g[16];
    g[19] = cpz[1] * g[17];
    g[20] = cpz[0] * g[18] + b01[0] * g[16];
    g[21] = cpz[1] * g[19] + b01[1] * g[17];
    g[22] = cpz[0] * g[20] + 2 * b01[0] * g[18];
    g[23] = cpz[1] * g[21] + 2 * b01[1] * g[19];
}
static inline void _g0_2d4d_0010(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    g[0] = 1;
    g[1] = cpx[0];
    g[2] = 1;
    g[3] = cpy[0];
    //g[4] = w[0];
    g[5] = cpz[0] * g[4];
}
static inline void _g0_2d4d_0011(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b01 = bc->b01;
    double xkxl = envs->rkrl[0];
    double ykyl = envs->rkrl[1];
    double zkzl = envs->rkrl[2];
    g[0] = 1;
    g[1] = 1;
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[0] * (xkxl + cpx[0]) + b01[0];
    g[7] = cpx[1] * (xkxl + cpx[1]) + b01[1];
    g[2] = xkxl + cpx[0];
    g[3] = xkxl + cpx[1];
    g[12] = 1;
    g[13] = 1;
    g[16] = cpy[0];
    g[17] = cpy[1];
    g[18] = cpy[0] * (ykyl + cpy[0]) + b01[0];
    g[19] = cpy[1] * (ykyl + cpy[1]) + b01[1];
    g[14] = ykyl + cpy[0];
    g[15] = ykyl + cpy[1];
    //g[24] = w[0];
    //g[25] = w[1];
    g[28] = cpz[0] * g[24];
    g[29] = cpz[1] * g[25];
    g[30] = g[28] * (zkzl + cpz[0]) + b01[0] * g[24];
    g[31] = g[29] * (zkzl + cpz[1]) + b01[1] * g[25];
    g[26] = g[24] * (zkzl + cpz[0]);
    g[27] = g[25] * (zkzl + cpz[1]);
}
static inline void _g0_2d4d_0012(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b01 = bc->b01;
    double xkxl = envs->rkrl[0];
    double ykyl = envs->rkrl[1];
    double zkzl = envs->rkrl[2];
    g[0] = 1;
    g[1] = 1;
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[8] = cpx[0] * cpx[0] + b01[0];
    g[9] = cpx[1] * cpx[1] + b01[1];
    g[10] = g[8] * (xkxl + cpx[0]) + cpx[0] * 2 * b01[0];
    g[11] = g[9] * (xkxl + cpx[1]) + cpx[1] * 2 * b01[1];
    g[6] = cpx[0] * (xkxl + cpx[0]) + b01[0];
    g[7] = cpx[1] * (xkxl + cpx[1]) + b01[1];
    g[2] = xkxl + cpx[0];
    g[3] = xkxl + cpx[1];
    g[16] = 1;
    g[17] = 1;
    g[20] = cpy[0];
    g[21] = cpy[1];
    g[24] = cpy[0] * cpy[0] + b01[0];
    g[25] = cpy[1] * cpy[1] + b01[1];
    g[26] = g[24] * (ykyl + cpy[0]) + cpy[0] * 2 * b01[0];
    g[27] = g[25] * (ykyl + cpy[1]) + cpy[1] * 2 * b01[1];
    g[22] = cpy[0] * (ykyl + cpy[0]) + b01[0];
    g[23] = cpy[1] * (ykyl + cpy[1]) + b01[1];
    g[18] = ykyl + cpy[0];
    g[19] = ykyl + cpy[1];
    //g[32] = w[0];
    //g[33] = w[1];
    g[36] = cpz[0] * g[32];
    g[37] = cpz[1] * g[33];
    g[40] = cpz[0] * g[36] + b01[0] * g[32];
    g[41] = cpz[1] * g[37] + b01[1] * g[33];
    g[42] = g[40] * (zkzl + cpz[0]) + 2 * b01[0] * g[36];
    g[43] = g[41] * (zkzl + cpz[1]) + 2 * b01[1] * g[37];
    g[38] = g[36] * (zkzl + cpz[0]) + b01[0] * g[32];
    g[39] = g[37] * (zkzl + cpz[1]) + b01[1] * g[33];
    g[34] = g[32] * (zkzl + cpz[0]);
    g[35] = g[33] * (zkzl + cpz[1]);
}
static inline void _g0_2d4d_0020(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = cpx[0];
    g[3] = cpx[1];
    g[4] = cpx[0] * cpx[0] + b01[0];
    g[5] = cpx[1] * cpx[1] + b01[1];
    g[6] = 1;
    g[7] = 1;
    g[8] = cpy[0];
    g[9] = cpy[1];
    g[10] = cpy[0] * cpy[0] + b01[0];
    g[11] = cpy[1] * cpy[1] + b01[1];
    //g[12] = w[0];
    //g[13] = w[1];
    g[14] = cpz[0] * g[12];
    g[15] = cpz[1] * g[13];
    g[16] = cpz[0] * g[14] + b01[0] * g[12];
    g[17] = cpz[1] * g[15] + b01[1] * g[13];
}
static inline void _g0_2d4d_0021(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b01 = bc->b01;
    double xkxl = envs->rkrl[0];
    double ykyl = envs->rkrl[1];
    double zkzl = envs->rkrl[2];
    g[0] = 1;
    g[1] = 1;
    g[2] = cpx[0];
    g[3] = cpx[1];
    g[4] = cpx[0] * cpx[0] + b01[0];
    g[5] = cpx[1] * cpx[1] + b01[1];
    g[8] = xkxl + cpx[0];
    g[9] = xkxl + cpx[1];
    g[10] = cpx[0] * (xkxl + cpx[0]) + b01[0];
    g[11] = cpx[1] * (xkxl + cpx[1]) + b01[1];
    g[12] = g[4] * (xkxl + cpx[0]) + cpx[0] * 2 * b01[0];
    g[13] = g[5] * (xkxl + cpx[1]) + cpx[1] * 2 * b01[1];
    g[16] = 1;
    g[17] = 1;
    g[18] = cpy[0];
    g[19] = cpy[1];
    g[20] = cpy[0] * cpy[0] + b01[0];
    g[21] = cpy[1] * cpy[1] + b01[1];
    g[24] = ykyl + cpy[0];
    g[25] = ykyl + cpy[1];
    g[26] = cpy[0] * (ykyl + cpy[0]) + b01[0];
    g[27] = cpy[1] * (ykyl + cpy[1]) + b01[1];
    g[28] = g[20] * (ykyl + cpy[0]) + cpy[0] * 2 * b01[0];
    g[29] = g[21] * (ykyl + cpy[1]) + cpy[1] * 2 * b01[1];
    //g[32] = w[0];
    //g[33] = w[1];
    g[34] = cpz[0] * g[32];
    g[35] = cpz[1] * g[33];
    g[36] = cpz[0] * g[34] + b01[0] * g[32];
    g[37] = cpz[1] * g[35] + b01[1] * g[33];
    g[40] = g[32] * (zkzl + cpz[0]);
    g[41] = g[33] * (zkzl + cpz[1]);
    g[42] = g[34] * (zkzl + cpz[0]) + b01[0] * g[32];
    g[43] = g[35] * (zkzl + cpz[1]) + b01[1] * g[33];
    g[44] = g[36] * (zkzl + cpz[0]) + 2 * b01[0] * g[34];
    g[45] = g[37] * (zkzl + cpz[1]) + 2 * b01[1] * g[35];
}
static inline void _g0_2d4d_0030(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = cpx[0];
    g[3] = cpx[1];
    g[4] = cpx[0] * cpx[0] + b01[0];
    g[5] = cpx[1] * cpx[1] + b01[1];
    g[6] = cpx[0] * (g[4] + 2 * b01[0]);
    g[7] = cpx[1] * (g[5] + 2 * b01[1]);
    g[8] = 1;
    g[9] = 1;
    g[10] = cpy[0];
    g[11] = cpy[1];
    g[12] = cpy[0] * cpy[0] + b01[0];
    g[13] = cpy[1] * cpy[1] + b01[1];
    g[14] = cpy[0] * (g[12] + 2 * b01[0]);
    g[15] = cpy[1] * (g[13] + 2 * b01[1]);
    //g[16] = w[0];
    //g[17] = w[1];
    g[18] = cpz[0] * g[16];
    g[19] = cpz[1] * g[17];
    g[20] = cpz[0] * g[18] + b01[0] * g[16];
    g[21] = cpz[1] * g[19] + b01[1] * g[17];
    g[22] = cpz[0] * g[20] + 2 * b01[0] * g[18];
    g[23] = cpz[1] * g[21] + 2 * b01[1] * g[19];
}
static inline void _g0_2d4d_0100(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    g[0] = 1;
    g[1] = c0x[0];
    g[2] = 1;
    g[3] = c0y[0];
    //g[4] = w[0];
    g[5] = c0z[0] * g[4];
}
static inline void _g0_2d4d_0101(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    g[0] = 1;
    g[1] = 1;
    g[2] = cpx[0];
    g[3] = cpx[1];
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = cpx[0] * c0x[0] + b00[0];
    g[7] = cpx[1] * c0x[1] + b00[1];
    g[8] = 1;
    g[9] = 1;
    g[10] = cpy[0];
    g[11] = cpy[1];
    g[12] = c0y[0];
    g[13] = c0y[1];
    g[14] = cpy[0] * c0y[0] + b00[0];
    g[15] = cpy[1] * c0y[1] + b00[1];
    //g[16] = w[0];
    //g[17] = w[1];
    g[18] = cpz[0] * g[16];
    g[19] = cpz[1] * g[17];
    g[20] = c0z[0] * g[16];
    g[21] = c0z[1] * g[17];
    g[22] = cpz[0] * g[20] + b00[0] * g[16];
    g[23] = cpz[1] * g[21] + b00[1] * g[17];
}
static inline void _g0_2d4d_0102(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    double* b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = cpx[0];
    g[3] = cpx[1];
    g[6] = c0x[0];
    g[7] = c0x[1];
    g[4] = cpx[0] * cpx[0] + b01[0];
    g[5] = cpx[1] * cpx[1] + b01[1];
    g[8] = cpx[0] * c0x[0] + b00[0];
    g[9] = cpx[1] * c0x[1] + b00[1];
    g[10] = cpx[0] * (g[8] + b00[0]) + b01[0] * c0x[0];
    g[11] = cpx[1] * (g[9] + b00[1]) + b01[1] * c0x[1];
    g[12] = 1;
    g[13] = 1;
    g[14] = cpy[0];
    g[15] = cpy[1];
    g[18] = c0y[0];
    g[19] = c0y[1];
    g[16] = cpy[0] * cpy[0] + b01[0];
    g[17] = cpy[1] * cpy[1] + b01[1];
    g[20] = cpy[0] * c0y[0] + b00[0];
    g[21] = cpy[1] * c0y[1] + b00[1];
    g[22] = cpy[0] * (g[20] + b00[0]) + b01[0] * c0y[0];
    g[23] = cpy[1] * (g[21] + b00[1]) + b01[1] * c0y[1];
    //g[24] = w[0];
    //g[25] = w[1];
    g[26] = cpz[0] * g[24];
    g[27] = cpz[1] * g[25];
    g[30] = c0z[0] * g[24];
    g[31] = c0z[1] * g[25];
    g[28] = cpz[0] * g[26] + b01[0] * g[24];
    g[29] = cpz[1] * g[27] + b01[1] * g[25];
    g[32] = cpz[0] * g[30] + b00[0] * g[24];
    g[33] = cpz[1] * g[31] + b00[1] * g[25];
    g[34] = cpz[0] * g[32] + b01[0] * g[30] + b00[0] * g[26];
    g[35] = cpz[1] * g[33] + b01[1] * g[31] + b00[1] * g[27];
}
static inline void _g0_2d4d_0110(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    g[0] = 1;
    g[1] = 1;
    g[2] = cpx[0];
    g[3] = cpx[1];
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = cpx[0] * c0x[0] + b00[0];
    g[7] = cpx[1] * c0x[1] + b00[1];
    g[8] = 1;
    g[9] = 1;
    g[10] = cpy[0];
    g[11] = cpy[1];
    g[12] = c0y[0];
    g[13] = c0y[1];
    g[14] = cpy[0] * c0y[0] + b00[0];
    g[15] = cpy[1] * c0y[1] + b00[1];
    //g[16] = w[0];
    //g[17] = w[1];
    g[18] = cpz[0] * g[16];
    g[19] = cpz[1] * g[17];
    g[20] = c0z[0] * g[16];
    g[21] = c0z[1] * g[17];
    g[22] = cpz[0] * g[20] + b00[0] * g[16];
    g[23] = cpz[1] * g[21] + b00[1] * g[17];
}
static inline void _g0_2d4d_0111(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    double* b01 = bc->b01;
    double xkxl = envs->rkrl[0];
    double ykyl = envs->rkrl[1];
    double zkzl = envs->rkrl[2];
    g[0] = 1;
    g[1] = 1;
    g[12] = c0x[0];
    g[13] = c0x[1];
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[16] = cpx[0] * c0x[0] + b00[0];
    g[17] = cpx[1] * c0x[1] + b00[1];
    g[6] = cpx[0] * (xkxl + cpx[0]) + b01[0];
    g[7] = cpx[1] * (xkxl + cpx[1]) + b01[1];
    g[18] = g[16] * (xkxl + cpx[0]) + cpx[0] * b00[0] + b01[0] * c0x[0];
    g[19] = g[17] * (xkxl + cpx[1]) + cpx[1] * b00[1] + b01[1] * c0x[1];
    g[2] = xkxl + cpx[0];
    g[3] = xkxl + cpx[1];
    g[14] = c0x[0] * (xkxl + cpx[0]) + b00[0];
    g[15] = c0x[1] * (xkxl + cpx[1]) + b00[1];
    g[24] = 1;
    g[25] = 1;
    g[36] = c0y[0];
    g[37] = c0y[1];
    g[28] = cpy[0];
    g[29] = cpy[1];
    g[40] = cpy[0] * c0y[0] + b00[0];
    g[41] = cpy[1] * c0y[1] + b00[1];
    g[30] = cpy[0] * (ykyl + cpy[0]) + b01[0];
    g[31] = cpy[1] * (ykyl + cpy[1]) + b01[1];
    g[42] = g[40] * (ykyl + cpy[0]) + cpy[0] * b00[0] + b01[0] * c0y[0];
    g[43] = g[41] * (ykyl + cpy[1]) + cpy[1] * b00[1] + b01[1] * c0y[1];
    g[26] = ykyl + cpy[0];
    g[27] = ykyl + cpy[1];
    g[38] = c0y[0] * (ykyl + cpy[0]) + b00[0];
    g[39] = c0y[1] * (ykyl + cpy[1]) + b00[1];
    //g[48] = w[0];
    //g[49] = w[1];
    g[60] = c0z[0] * g[48];
    g[61] = c0z[1] * g[49];
    g[52] = cpz[0] * g[48];
    g[53] = cpz[1] * g[49];
    g[64] = cpz[0] * g[60] + b00[0] * g[48];
    g[65] = cpz[1] * g[61] + b00[1] * g[49];
    g[54] = g[52] * (zkzl + cpz[0]) + b01[0] * g[48];
    g[55] = g[53] * (zkzl + cpz[1]) + b01[1] * g[49];
    g[66] = g[64] * (zkzl + cpz[0]) + b01[0] * g[60] + b00[0] * g[52];
    g[67] = g[65] * (zkzl + cpz[1]) + b01[1] * g[61] + b00[1] * g[53];
    g[50] = g[48] * (zkzl + cpz[0]);
    g[51] = g[49] * (zkzl + cpz[1]);
    g[62] = g[60] * (zkzl + cpz[0]) + b00[0] * g[48];
    g[63] = g[61] * (zkzl + cpz[1]) + b00[1] * g[49];
}
static inline void _g0_2d4d_0120(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    double* b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = cpx[0];
    g[3] = cpx[1];
    g[6] = c0x[0];
    g[7] = c0x[1];
    g[4] = cpx[0] * cpx[0] + b01[0];
    g[5] = cpx[1] * cpx[1] + b01[1];
    g[8] = cpx[0] * c0x[0] + b00[0];
    g[9] = cpx[1] * c0x[1] + b00[1];
    g[10] = cpx[0] * (g[8] + b00[0]) + b01[0] * c0x[0];
    g[11] = cpx[1] * (g[9] + b00[1]) + b01[1] * c0x[1];
    g[12] = 1;
    g[13] = 1;
    g[14] = cpy[0];
    g[15] = cpy[1];
    g[18] = c0y[0];
    g[19] = c0y[1];
    g[16] = cpy[0] * cpy[0] + b01[0];
    g[17] = cpy[1] * cpy[1] + b01[1];
    g[20] = cpy[0] * c0y[0] + b00[0];
    g[21] = cpy[1] * c0y[1] + b00[1];
    g[22] = cpy[0] * (g[20] + b00[0]) + b01[0] * c0y[0];
    g[23] = cpy[1] * (g[21] + b00[1]) + b01[1] * c0y[1];
    //g[24] = w[0];
    //g[25] = w[1];
    g[26] = cpz[0] * g[24];
    g[27] = cpz[1] * g[25];
    g[30] = c0z[0] * g[24];
    g[31] = c0z[1] * g[25];
    g[28] = cpz[0] * g[26] + b01[0] * g[24];
    g[29] = cpz[1] * g[27] + b01[1] * g[25];
    g[32] = cpz[0] * g[30] + b00[0] * g[24];
    g[33] = cpz[1] * g[31] + b00[1] * g[25];
    g[34] = cpz[0] * g[32] + b01[0] * g[30] + b00[0] * g[26];
    g[35] = cpz[1] * g[33] + b01[1] * g[31] + b00[1] * g[27];
}
static inline void _g0_2d4d_0200(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[2] = c0x[0];
    g[3] = c0x[1];
    g[4] = c0x[0] * c0x[0] + b10[0];
    g[5] = c0x[1] * c0x[1] + b10[1];
    g[6] = 1;
    g[7] = 1;
    g[8] = c0y[0];
    g[9] = c0y[1];
    g[10] = c0y[0] * c0y[0] + b10[0];
    g[11] = c0y[1] * c0y[1] + b10[1];
    //g[12] = w[0];
    //g[13] = w[1];
    g[14] = c0z[0] * g[12];
    g[15] = c0z[1] * g[13];
    g[16] = c0z[0] * g[14] + b10[0] * g[12];
    g[17] = c0z[1] * g[15] + b10[1] * g[13];
}
static inline void _g0_2d4d_0201(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    double* b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[8] = c0x[0] * c0x[0] + b10[0];
    g[9] = c0x[1] * c0x[1] + b10[1];
    g[2] = cpx[0];
    g[3] = cpx[1];
    g[6] = cpx[0] * c0x[0] + b00[0];
    g[7] = cpx[1] * c0x[1] + b00[1];
    g[10] = c0x[0] * (g[6] + b00[0]) + b10[0] * cpx[0];
    g[11] = c0x[1] * (g[7] + b00[1]) + b10[1] * cpx[1];
    g[12] = 1;
    g[13] = 1;
    g[16] = c0y[0];
    g[17] = c0y[1];
    g[20] = c0y[0] * c0y[0] + b10[0];
    g[21] = c0y[1] * c0y[1] + b10[1];
    g[14] = cpy[0];
    g[15] = cpy[1];
    g[18] = cpy[0] * c0y[0] + b00[0];
    g[19] = cpy[1] * c0y[1] + b00[1];
    g[22] = c0y[0] * (g[18] + b00[0]) + b10[0] * cpy[0];
    g[23] = c0y[1] * (g[19] + b00[1]) + b10[1] * cpy[1];
    //g[24] = w[0];
    //g[25] = w[1];
    g[28] = c0z[0] * g[24];
    g[29] = c0z[1] * g[25];
    g[32] = c0z[0] * g[28] + b10[0] * g[24];
    g[33] = c0z[1] * g[29] + b10[1] * g[25];
    g[26] = cpz[0] * g[24];
    g[27] = cpz[1] * g[25];
    g[30] = cpz[0] * g[28] + b00[0] * g[24];
    g[31] = cpz[1] * g[29] + b00[1] * g[25];
    g[34] = c0z[0] * g[30] + b10[0] * g[26] + b00[0] * g[28];
    g[35] = c0z[1] * g[31] + b10[1] * g[27] + b00[1] * g[29];
}
static inline void _g0_2d4d_0210(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    double* b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[2] = cpx[0];
    g[3] = cpx[1];
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = cpx[0] * c0x[0] + b00[0];
    g[7] = cpx[1] * c0x[1] + b00[1];
    g[8] = c0x[0] * c0x[0] + b10[0];
    g[9] = c0x[1] * c0x[1] + b10[1];
    g[10] = c0x[0] * (g[6] + b00[0]) + b10[0] * cpx[0];
    g[11] = c0x[1] * (g[7] + b00[1]) + b10[1] * cpx[1];
    g[12] = 1;
    g[13] = 1;
    g[14] = cpy[0];
    g[15] = cpy[1];
    g[16] = c0y[0];
    g[17] = c0y[1];
    g[18] = cpy[0] * c0y[0] + b00[0];
    g[19] = cpy[1] * c0y[1] + b00[1];
    g[20] = c0y[0] * c0y[0] + b10[0];
    g[21] = c0y[1] * c0y[1] + b10[1];
    g[22] = c0y[0] * (g[18] + b00[0]) + b10[0] * cpy[0];
    g[23] = c0y[1] * (g[19] + b00[1]) + b10[1] * cpy[1];
    //g[24] = w[0];
    //g[25] = w[1];
    g[26] = cpz[0] * g[24];
    g[27] = cpz[1] * g[25];
    g[28] = c0z[0] * g[24];
    g[29] = c0z[1] * g[25];
    g[30] = cpz[0] * g[28] + b00[0] * g[24];
    g[31] = cpz[1] * g[29] + b00[1] * g[25];
    g[32] = c0z[0] * g[28] + b10[0] * g[24];
    g[33] = c0z[1] * g[29] + b10[1] * g[25];
    g[34] = c0z[0] * g[30] + b10[0] * g[26] + b00[0] * g[28];
    g[35] = c0z[1] * g[31] + b10[1] * g[27] + b00[1] * g[29];
}
static inline void _g0_2d4d_0300(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[2] = c0x[0];
    g[3] = c0x[1];
    g[4] = c0x[0] * c0x[0] + b10[0];
    g[5] = c0x[1] * c0x[1] + b10[1];
    g[6] = c0x[0] * (g[4] + 2 * b10[0]);
    g[7] = c0x[1] * (g[5] + 2 * b10[1]);
    g[8] = 1;
    g[9] = 1;
    g[10] = c0y[0];
    g[11] = c0y[1];
    g[12] = c0y[0] * c0y[0] + b10[0];
    g[13] = c0y[1] * c0y[1] + b10[1];
    g[14] = c0y[0] * (g[12] + 2 * b10[0]);
    g[15] = c0y[1] * (g[13] + 2 * b10[1]);
    //g[16] = w[0];
    //g[17] = w[1];
    g[18] = c0z[0] * g[16];
    g[19] = c0z[1] * g[17];
    g[20] = c0z[0] * g[18] + b10[0] * g[16];
    g[21] = c0z[1] * g[19] + b10[1] * g[17];
    g[22] = c0z[0] * g[20] + 2 * b10[0] * g[18];
    g[23] = c0z[1] * g[21] + 2 * b10[1] * g[19];
}
static inline void _g0_2d4d_1000(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    g[0] = 1;
    g[1] = c0x[0];
    g[2] = 1;
    g[3] = c0y[0];
    //g[4] = w[0];
    g[5] = c0z[0] * g[4];
}
static inline void _g0_2d4d_1001(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    g[0] = 1;
    g[1] = 1;
    g[2] = c0x[0];
    g[3] = c0x[1];
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[0] * c0x[0] + b00[0];
    g[7] = cpx[1] * c0x[1] + b00[1];
    g[8] = 1;
    g[9] = 1;
    g[10] = c0y[0];
    g[11] = c0y[1];
    g[12] = cpy[0];
    g[13] = cpy[1];
    g[14] = cpy[0] * c0y[0] + b00[0];
    g[15] = cpy[1] * c0y[1] + b00[1];
    //g[16] = w[0];
    //g[17] = w[1];
    g[18] = c0z[0] * g[16];
    g[19] = c0z[1] * g[17];
    g[20] = cpz[0] * g[16];
    g[21] = cpz[1] * g[17];
    g[22] = cpz[0] * g[18] + b00[0] * g[16];
    g[23] = cpz[1] * g[19] + b00[1] * g[17];
}
static inline void _g0_2d4d_1002(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    double* b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = c0x[0];
    g[3] = c0x[1];
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[0] * c0x[0] + b00[0];
    g[7] = cpx[1] * c0x[1] + b00[1];
    g[8] = cpx[0] * cpx[0] + b01[0];
    g[9] = cpx[1] * cpx[1] + b01[1];
    g[10] = cpx[0] * (g[6] + b00[0]) + b01[0] * c0x[0];
    g[11] = cpx[1] * (g[7] + b00[1]) + b01[1] * c0x[1];
    g[12] = 1;
    g[13] = 1;
    g[14] = c0y[0];
    g[15] = c0y[1];
    g[16] = cpy[0];
    g[17] = cpy[1];
    g[18] = cpy[0] * c0y[0] + b00[0];
    g[19] = cpy[1] * c0y[1] + b00[1];
    g[20] = cpy[0] * cpy[0] + b01[0];
    g[21] = cpy[1] * cpy[1] + b01[1];
    g[22] = cpy[0] * (g[18] + b00[0]) + b01[0] * c0y[0];
    g[23] = cpy[1] * (g[19] + b00[1]) + b01[1] * c0y[1];
    //g[24] = w[0];
    //g[25] = w[1];
    g[26] = c0z[0] * g[24];
    g[27] = c0z[1] * g[25];
    g[28] = cpz[0] * g[24];
    g[29] = cpz[1] * g[25];
    g[30] = cpz[0] * g[26] + b00[0] * g[24];
    g[31] = cpz[1] * g[27] + b00[1] * g[25];
    g[32] = cpz[0] * g[28] + b01[0] * g[24];
    g[33] = cpz[1] * g[29] + b01[1] * g[25];
    g[34] = cpz[0] * g[30] + b01[0] * g[26] + b00[0] * g[28];
    g[35] = cpz[1] * g[31] + b01[1] * g[27] + b00[1] * g[29];
}
static inline void _g0_2d4d_1010(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    g[0] = 1;
    g[1] = 1;
    g[2] = c0x[0];
    g[3] = c0x[1];
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[0] * c0x[0] + b00[0];
    g[7] = cpx[1] * c0x[1] + b00[1];
    g[8] = 1;
    g[9] = 1;
    g[10] = c0y[0];
    g[11] = c0y[1];
    g[12] = cpy[0];
    g[13] = cpy[1];
    g[14] = cpy[0] * c0y[0] + b00[0];
    g[15] = cpy[1] * c0y[1] + b00[1];
    //g[16] = w[0];
    //g[17] = w[1];
    g[18] = c0z[0] * g[16];
    g[19] = c0z[1] * g[17];
    g[20] = cpz[0] * g[16];
    g[21] = cpz[1] * g[17];
    g[22] = cpz[0] * g[18] + b00[0] * g[16];
    g[23] = cpz[1] * g[19] + b00[1] * g[17];
}
static inline void _g0_2d4d_1011(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    double* b01 = bc->b01;
    double xkxl = envs->rkrl[0];
    double ykyl = envs->rkrl[1];
    double zkzl = envs->rkrl[2];
    g[0] = 1;
    g[1] = 1;
    g[2] = c0x[0];
    g[3] = c0x[1];
    g[8] = cpx[0];
    g[9] = cpx[1];
    g[10] = cpx[0] * c0x[0] + b00[0];
    g[11] = cpx[1] * c0x[1] + b00[1];
    g[12] = cpx[0] * (xkxl + cpx[0]) + b01[0];
    g[13] = cpx[1] * (xkxl + cpx[1]) + b01[1];
    g[14] = g[10] * (xkxl + cpx[0]) + cpx[0] * b00[0] + b01[0] * c0x[0];
    g[15] = g[11] * (xkxl + cpx[1]) + cpx[1] * b00[1] + b01[1] * c0x[1];
    g[4] = xkxl + cpx[0];
    g[5] = xkxl + cpx[1];
    g[6] = c0x[0] * (xkxl + cpx[0]) + b00[0];
    g[7] = c0x[1] * (xkxl + cpx[1]) + b00[1];
    g[24] = 1;
    g[25] = 1;
    g[26] = c0y[0];
    g[27] = c0y[1];
    g[32] = cpy[0];
    g[33] = cpy[1];
    g[34] = cpy[0] * c0y[0] + b00[0];
    g[35] = cpy[1] * c0y[1] + b00[1];
    g[36] = cpy[0] * (ykyl + cpy[0]) + b01[0];
    g[37] = cpy[1] * (ykyl + cpy[1]) + b01[1];
    g[38] = g[34] * (ykyl + cpy[0]) + cpy[0] * b00[0] + b01[0] * c0y[0];
    g[39] = g[35] * (ykyl + cpy[1]) + cpy[1] * b00[1] + b01[1] * c0y[1];
    g[28] = ykyl + cpy[0];
    g[29] = ykyl + cpy[1];
    g[30] = c0y[0] * (ykyl + cpy[0]) + b00[0];
    g[31] = c0y[1] * (ykyl + cpy[1]) + b00[1];
    //g[48] = w[0];
    //g[49] = w[1];
    g[50] = c0z[0] * g[48];
    g[51] = c0z[1] * g[49];
    g[56] = cpz[0] * g[48];
    g[57] = cpz[1] * g[49];
    g[58] = cpz[0] * g[50] + b00[0] * g[48];
    g[59] = cpz[1] * g[51] + b00[1] * g[49];
    g[60] = g[56] * (zkzl + cpz[0]) + b01[0] * g[48];
    g[61] = g[57] * (zkzl + cpz[1]) + b01[1] * g[49];
    g[62] = g[58] * (zkzl + cpz[0]) + b01[0] * g[50] + b00[0] * g[56];
    g[63] = g[59] * (zkzl + cpz[1]) + b01[1] * g[51] + b00[1] * g[57];
    g[52] = g[48] * (zkzl + cpz[0]);
    g[53] = g[49] * (zkzl + cpz[1]);
    g[54] = g[50] * (zkzl + cpz[0]) + b00[0] * g[48];
    g[55] = g[51] * (zkzl + cpz[1]) + b00[1] * g[49];
}
static inline void _g0_2d4d_1020(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    double* b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = c0x[0];
    g[3] = c0x[1];
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[0] * c0x[0] + b00[0];
    g[7] = cpx[1] * c0x[1] + b00[1];
    g[8] = cpx[0] * cpx[0] + b01[0];
    g[9] = cpx[1] * cpx[1] + b01[1];
    g[10] = cpx[0] * (g[6] + b00[0]) + b01[0] * c0x[0];
    g[11] = cpx[1] * (g[7] + b00[1]) + b01[1] * c0x[1];
    g[12] = 1;
    g[13] = 1;
    g[14] = c0y[0];
    g[15] = c0y[1];
    g[16] = cpy[0];
    g[17] = cpy[1];
    g[18] = cpy[0] * c0y[0] + b00[0];
    g[19] = cpy[1] * c0y[1] + b00[1];
    g[20] = cpy[0] * cpy[0] + b01[0];
    g[21] = cpy[1] * cpy[1] + b01[1];
    g[22] = cpy[0] * (g[18] + b00[0]) + b01[0] * c0y[0];
    g[23] = cpy[1] * (g[19] + b00[1]) + b01[1] * c0y[1];
    //g[24] = w[0];
    //g[25] = w[1];
    g[26] = c0z[0] * g[24];
    g[27] = c0z[1] * g[25];
    g[28] = cpz[0] * g[24];
    g[29] = cpz[1] * g[25];
    g[30] = cpz[0] * g[26] + b00[0] * g[24];
    g[31] = cpz[1] * g[27] + b00[1] * g[25];
    g[32] = cpz[0] * g[28] + b01[0] * g[24];
    g[33] = cpz[1] * g[29] + b01[1] * g[25];
    g[34] = cpz[0] * g[30] + b01[0] * g[26] + b00[0] * g[28];
    g[35] = cpz[1] * g[31] + b01[1] * g[27] + b00[1] * g[29];
}
static inline void _g0_2d4d_1100(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* b10 = bc->b10;
    double xixj = envs->rirj[0];
    double yiyj = envs->rirj[1];
    double zizj = envs->rirj[2];
    g[0] = 1;
    g[1] = 1;
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = c0x[0] * (xixj + c0x[0]) + b10[0];
    g[7] = c0x[1] * (xixj + c0x[1]) + b10[1];
    g[2] = xixj + c0x[0];
    g[3] = xixj + c0x[1];
    g[12] = 1;
    g[13] = 1;
    g[16] = c0y[0];
    g[17] = c0y[1];
    g[18] = c0y[0] * (yiyj + c0y[0]) + b10[0];
    g[19] = c0y[1] * (yiyj + c0y[1]) + b10[1];
    g[14] = yiyj + c0y[0];
    g[15] = yiyj + c0y[1];
    //g[24] = w[0];
    //g[25] = w[1];
    g[28] = c0z[0] * g[24];
    g[29] = c0z[1] * g[25];
    g[30] = g[28] * (zizj + c0z[0]) + b10[0] * g[24];
    g[31] = g[29] * (zizj + c0z[1]) + b10[1] * g[25];
    g[26] = g[24] * (zizj + c0z[0]);
    g[27] = g[25] * (zizj + c0z[1]);
}
static inline void _g0_2d4d_1101(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    double* b10 = bc->b10;
    double xixj = envs->rirj[0];
    double yiyj = envs->rirj[1];
    double zizj = envs->rirj[2];
    g[0] = 1;
    g[1] = 1;
    g[8] = c0x[0];
    g[9] = c0x[1];
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[12] = cpx[0] * c0x[0] + b00[0];
    g[13] = cpx[1] * c0x[1] + b00[1];
    g[10] = c0x[0] * (xixj + c0x[0]) + b10[0];
    g[11] = c0x[1] * (xixj + c0x[1]) + b10[1];
    g[2] = xixj + c0x[0];
    g[3] = xixj + c0x[1];
    g[14] = g[12] * (xixj + c0x[0]) + c0x[0] * b00[0] + b10[0] * cpx[0];
    g[15] = g[13] * (xixj + c0x[1]) + c0x[1] * b00[1] + b10[1] * cpx[1];
    g[6] = cpx[0] * (xixj + c0x[0]) + b00[0];
    g[7] = cpx[1] * (xixj + c0x[1]) + b00[1];
    g[24] = 1;
    g[25] = 1;
    g[32] = c0y[0];
    g[33] = c0y[1];
    g[28] = cpy[0];
    g[29] = cpy[1];
    g[36] = cpy[0] * c0y[0] + b00[0];
    g[37] = cpy[1] * c0y[1] + b00[1];
    g[34] = c0y[0] * (yiyj + c0y[0]) + b10[0];
    g[35] = c0y[1] * (yiyj + c0y[1]) + b10[1];
    g[26] = yiyj + c0y[0];
    g[27] = yiyj + c0y[1];
    g[38] = g[36] * (yiyj + c0y[0]) + c0y[0] * b00[0] + b10[0] * cpy[0];
    g[39] = g[37] * (yiyj + c0y[1]) + c0y[1] * b00[1] + b10[1] * cpy[1];
    g[30] = cpy[0] * (yiyj + c0y[0]) + b00[0];
    g[31] = cpy[1] * (yiyj + c0y[1]) + b00[1];
    //g[48] = w[0];
    //g[49] = w[1];
    g[56] = c0z[0] * g[48];
    g[57] = c0z[1] * g[49];
    g[52] = cpz[0] * g[48];
    g[53] = cpz[1] * g[49];
    g[60] = cpz[0] * g[56] + b00[0] * g[48];
    g[61] = cpz[1] * g[57] + b00[1] * g[49];
    g[58] = g[56] * (zizj + c0z[0]) + b10[0] * g[48];
    g[59] = g[57] * (zizj + c0z[1]) + b10[1] * g[49];
    g[50] = g[48] * (zizj + c0z[0]);
    g[51] = g[49] * (zizj + c0z[1]);
    g[62] = g[60] * (zizj + c0z[0]) + b10[0] * g[52] + b00[0] * g[56];
    g[63] = g[61] * (zizj + c0z[1]) + b10[1] * g[53] + b00[1] * g[57];
    g[54] = zizj * g[52] + cpz[0] * g[56] + b00[0] * g[48];
    g[55] = zizj * g[53] + cpz[1] * g[57] + b00[1] * g[49];
}
static inline void _g0_2d4d_1110(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    double* b10 = bc->b10;
    double xixj = envs->rirj[0];
    double yiyj = envs->rirj[1];
    double zizj = envs->rirj[2];
    g[0] = 1;
    g[1] = 1;
    g[8] = c0x[0];
    g[9] = c0x[1];
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[12] = cpx[0] * c0x[0] + b00[0];
    g[13] = cpx[1] * c0x[1] + b00[1];
    g[10] = c0x[0] * (xixj + c0x[0]) + b10[0];
    g[11] = c0x[1] * (xixj + c0x[1]) + b10[1];
    g[2] = xixj + c0x[0];
    g[3] = xixj + c0x[1];
    g[14] = g[12] * (xixj + c0x[0]) + c0x[0] * b00[0] + b10[0] * cpx[0];
    g[15] = g[13] * (xixj + c0x[1]) + c0x[1] * b00[1] + b10[1] * cpx[1];
    g[6] = cpx[0] * (xixj + c0x[0]) + b00[0];
    g[7] = cpx[1] * (xixj + c0x[1]) + b00[1];
    g[24] = 1;
    g[25] = 1;
    g[32] = c0y[0];
    g[33] = c0y[1];
    g[28] = cpy[0];
    g[29] = cpy[1];
    g[36] = cpy[0] * c0y[0] + b00[0];
    g[37] = cpy[1] * c0y[1] + b00[1];
    g[34] = c0y[0] * (yiyj + c0y[0]) + b10[0];
    g[35] = c0y[1] * (yiyj + c0y[1]) + b10[1];
    g[26] = yiyj + c0y[0];
    g[27] = yiyj + c0y[1];
    g[38] = g[36] * (yiyj + c0y[0]) + c0y[0] * b00[0] + b10[0] * cpy[0];
    g[39] = g[37] * (yiyj + c0y[1]) + c0y[1] * b00[1] + b10[1] * cpy[1];
    g[30] = cpy[0] * (yiyj + c0y[0]) + b00[0];
    g[31] = cpy[1] * (yiyj + c0y[1]) + b00[1];
    //g[48] = w[0];
    //g[49] = w[1];
    g[56] = c0z[0] * g[48];
    g[57] = c0z[1] * g[49];
    g[52] = cpz[0] * g[48];
    g[53] = cpz[1] * g[49];
    g[60] = cpz[0] * g[56] + b00[0] * g[48];
    g[61] = cpz[1] * g[57] + b00[1] * g[49];
    g[58] = g[56] * (zizj + c0z[0]) + b10[0] * g[48];
    g[59] = g[57] * (zizj + c0z[1]) + b10[1] * g[49];
    g[50] = g[48] * (zizj + c0z[0]);
    g[51] = g[49] * (zizj + c0z[1]);
    g[62] = g[60] * (zizj + c0z[0]) + b10[0] * g[52] + b00[0] * g[56];
    g[63] = g[61] * (zizj + c0z[1]) + b10[1] * g[53] + b00[1] * g[57];
    g[54] = zizj * g[52] + cpz[0] * g[56] + b00[0] * g[48];
    g[55] = zizj * g[53] + cpz[1] * g[57] + b00[1] * g[49];
}
static inline void _g0_2d4d_1200(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* b10 = bc->b10;
    double xixj = envs->rirj[0];
    double yiyj = envs->rirj[1];
    double zizj = envs->rirj[2];
    g[0] = 1;
    g[1] = 1;
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[8] = c0x[0] * c0x[0] + b10[0];
    g[9] = c0x[1] * c0x[1] + b10[1];
    g[10] = g[8] * (xixj + c0x[0]) + c0x[0] * 2 * b10[0];
    g[11] = g[9] * (xixj + c0x[1]) + c0x[1] * 2 * b10[1];
    g[6] = c0x[0] * (xixj + c0x[0]) + b10[0];
    g[7] = c0x[1] * (xixj + c0x[1]) + b10[1];
    g[2] = xixj + c0x[0];
    g[3] = xixj + c0x[1];
    g[16] = 1;
    g[17] = 1;
    g[20] = c0y[0];
    g[21] = c0y[1];
    g[24] = c0y[0] * c0y[0] + b10[0];
    g[25] = c0y[1] * c0y[1] + b10[1];
    g[26] = g[24] * (yiyj + c0y[0]) + c0y[0] * 2 * b10[0];
    g[27] = g[25] * (yiyj + c0y[1]) + c0y[1] * 2 * b10[1];
    g[22] = c0y[0] * (yiyj + c0y[0]) + b10[0];
    g[23] = c0y[1] * (yiyj + c0y[1]) + b10[1];
    g[18] = yiyj + c0y[0];
    g[19] = yiyj + c0y[1];
    //g[32] = w[0];
    //g[33] = w[1];
    g[36] = c0z[0] * g[32];
    g[37] = c0z[1] * g[33];
    g[40] = c0z[0] * g[36] + b10[0] * g[32];
    g[41] = c0z[1] * g[37] + b10[1] * g[33];
    g[42] = g[40] * (zizj + c0z[0]) + 2 * b10[0] * g[36];
    g[43] = g[41] * (zizj + c0z[1]) + 2 * b10[1] * g[37];
    g[38] = g[36] * (zizj + c0z[0]) + b10[0] * g[32];
    g[39] = g[37] * (zizj + c0z[1]) + b10[1] * g[33];
    g[34] = g[32] * (zizj + c0z[0]);
    g[35] = g[33] * (zizj + c0z[1]);
}
static inline void _g0_2d4d_2000(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[2] = c0x[0];
    g[3] = c0x[1];
    g[4] = c0x[0] * c0x[0] + b10[0];
    g[5] = c0x[1] * c0x[1] + b10[1];
    g[6] = 1;
    g[7] = 1;
    g[8] = c0y[0];
    g[9] = c0y[1];
    g[10] = c0y[0] * c0y[0] + b10[0];
    g[11] = c0y[1] * c0y[1] + b10[1];
    //g[12] = w[0];
    //g[13] = w[1];
    g[14] = c0z[0] * g[12];
    g[15] = c0z[1] * g[13];
    g[16] = c0z[0] * g[14] + b10[0] * g[12];
    g[17] = c0z[1] * g[15] + b10[1] * g[13];
}
static inline void _g0_2d4d_2001(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    double* b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[2] = c0x[0];
    g[3] = c0x[1];
    g[4] = c0x[0] * c0x[0] + b10[0];
    g[5] = c0x[1] * c0x[1] + b10[1];
    g[6] = cpx[0];
    g[7] = cpx[1];
    g[8] = cpx[0] * c0x[0] + b00[0];
    g[9] = cpx[1] * c0x[1] + b00[1];
    g[10] = c0x[0] * (g[8] + b00[0]) + b10[0] * cpx[0];
    g[11] = c0x[1] * (g[9] + b00[1]) + b10[1] * cpx[1];
    g[12] = 1;
    g[13] = 1;
    g[14] = c0y[0];
    g[15] = c0y[1];
    g[16] = c0y[0] * c0y[0] + b10[0];
    g[17] = c0y[1] * c0y[1] + b10[1];
    g[18] = cpy[0];
    g[19] = cpy[1];
    g[20] = cpy[0] * c0y[0] + b00[0];
    g[21] = cpy[1] * c0y[1] + b00[1];
    g[22] = c0y[0] * (g[20] + b00[0]) + b10[0] * cpy[0];
    g[23] = c0y[1] * (g[21] + b00[1]) + b10[1] * cpy[1];
    //g[24] = w[0];
    //g[25] = w[1];
    g[26] = c0z[0] * g[24];
    g[27] = c0z[1] * g[25];
    g[28] = c0z[0] * g[26] + b10[0] * g[24];
    g[29] = c0z[1] * g[27] + b10[1] * g[25];
    g[30] = cpz[0] * g[24];
    g[31] = cpz[1] * g[25];
    g[32] = cpz[0] * g[26] + b00[0] * g[24];
    g[33] = cpz[1] * g[27] + b00[1] * g[25];
    g[34] = c0z[0] * g[32] + b10[0] * g[30] + b00[0] * g[26];
    g[35] = c0z[1] * g[33] + b10[1] * g[31] + b00[1] * g[27];
}
static inline void _g0_2d4d_2010(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    double* b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[2] = c0x[0];
    g[3] = c0x[1];
    g[4] = c0x[0] * c0x[0] + b10[0];
    g[5] = c0x[1] * c0x[1] + b10[1];
    g[6] = cpx[0];
    g[7] = cpx[1];
    g[8] = cpx[0] * c0x[0] + b00[0];
    g[9] = cpx[1] * c0x[1] + b00[1];
    g[10] = c0x[0] * (g[8] + b00[0]) + b10[0] * cpx[0];
    g[11] = c0x[1] * (g[9] + b00[1]) + b10[1] * cpx[1];
    g[12] = 1;
    g[13] = 1;
    g[14] = c0y[0];
    g[15] = c0y[1];
    g[16] = c0y[0] * c0y[0] + b10[0];
    g[17] = c0y[1] * c0y[1] + b10[1];
    g[18] = cpy[0];
    g[19] = cpy[1];
    g[20] = cpy[0] * c0y[0] + b00[0];
    g[21] = cpy[1] * c0y[1] + b00[1];
    g[22] = c0y[0] * (g[20] + b00[0]) + b10[0] * cpy[0];
    g[23] = c0y[1] * (g[21] + b00[1]) + b10[1] * cpy[1];
    //g[24] = w[0];
    //g[25] = w[1];
    g[26] = c0z[0] * g[24];
    g[27] = c0z[1] * g[25];
    g[28] = c0z[0] * g[26] + b10[0] * g[24];
    g[29] = c0z[1] * g[27] + b10[1] * g[25];
    g[30] = cpz[0] * g[24];
    g[31] = cpz[1] * g[25];
    g[32] = cpz[0] * g[26] + b00[0] * g[24];
    g[33] = cpz[1] * g[27] + b00[1] * g[25];
    g[34] = c0z[0] * g[32] + b10[0] * g[30] + b00[0] * g[26];
    g[35] = c0z[1] * g[33] + b10[1] * g[31] + b00[1] * g[27];
}
static inline void _g0_2d4d_2100(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* b10 = bc->b10;
    double xixj = envs->rirj[0];
    double yiyj = envs->rirj[1];
    double zizj = envs->rirj[2];
    g[0] = 1;
    g[1] = 1;
    g[2] = c0x[0];
    g[3] = c0x[1];
    g[4] = c0x[0] * c0x[0] + b10[0];
    g[5] = c0x[1] * c0x[1] + b10[1];
    g[12] = g[4] * (xixj + c0x[0]) + c0x[0] * 2 * b10[0];
    g[13] = g[5] * (xixj + c0x[1]) + c0x[1] * 2 * b10[1];
    g[10] = c0x[0] * (xixj + c0x[0]) + b10[0];
    g[11] = c0x[1] * (xixj + c0x[1]) + b10[1];
    g[8] = xixj + c0x[0];
    g[9] = xixj + c0x[1];
    g[16] = 1;
    g[17] = 1;
    g[18] = c0y[0];
    g[19] = c0y[1];
    g[20] = c0y[0] * c0y[0] + b10[0];
    g[21] = c0y[1] * c0y[1] + b10[1];
    g[28] = g[20] * (yiyj + c0y[0]) + c0y[0] * 2 * b10[0];
    g[29] = g[21] * (yiyj + c0y[1]) + c0y[1] * 2 * b10[1];
    g[26] = c0y[0] * (yiyj + c0y[0]) + b10[0];
    g[27] = c0y[1] * (yiyj + c0y[1]) + b10[1];
    g[24] = yiyj + c0y[0];
    g[25] = yiyj + c0y[1];
    //g[32] = w[0];
    //g[33] = w[1];
    g[34] = c0z[0] * g[32];
    g[35] = c0z[1] * g[33];
    g[36] = c0z[0] * g[34] + b10[0] * g[32];
    g[37] = c0z[1] * g[35] + b10[1] * g[33];
    g[44] = g[36] * (zizj + c0z[0]) + 2 * b10[0] * g[34];
    g[45] = g[37] * (zizj + c0z[1]) + 2 * b10[1] * g[35];
    g[42] = g[34] * (zizj + c0z[0]) + b10[0] * g[32];
    g[43] = g[35] * (zizj + c0z[1]) + b10[1] * g[33];
    g[40] = g[32] * (zizj + c0z[0]);
    g[41] = g[33] * (zizj + c0z[1]);
}
static inline void _g0_2d4d_3000(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[2] = c0x[0];
    g[3] = c0x[1];
    g[4] = c0x[0] * c0x[0] + b10[0];
    g[5] = c0x[1] * c0x[1] + b10[1];
    g[6] = c0x[0] * (g[4] + 2 * b10[0]);
    g[7] = c0x[1] * (g[5] + 2 * b10[1]);
    g[8] = 1;
    g[9] = 1;
    g[10] = c0y[0];
    g[11] = c0y[1];
    g[12] = c0y[0] * c0y[0] + b10[0];
    g[13] = c0y[1] * c0y[1] + b10[1];
    g[14] = c0y[0] * (g[12] + 2 * b10[0]);
    g[15] = c0y[1] * (g[13] + 2 * b10[1]);
    //g[16] = w[0];
    //g[17] = w[1];
    g[18] = c0z[0] * g[16];
    g[19] = c0z[1] * g[17];
    g[20] = c0z[0] * g[18] + b10[0] * g[16];
    g[21] = c0z[1] * g[19] + b10[1] * g[17];
    g[22] = c0z[0] * g[20] + 2 * b10[0] * g[18];
    g[23] = c0z[1] * g[21] + 2 * b10[1] * g[19];
}

void CINTg0_2e_2d4d_unrolled(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    int type_ijkl = ((envs->li_ceil << 6) | (envs->lj_ceil << 4) |
        (envs->lk_ceil << 2) | (envs->ll_ceil));
    switch (type_ijkl) {
    case 0b00000000: _g0_2d4d_0000(g, bc, envs); return;
    case 0b00000001: _g0_2d4d_0001(g, bc, envs); return;
    case 0b00000010: _g0_2d4d_0002(g, bc, envs); return;
    case 0b00000011: _g0_2d4d_0003(g, bc, envs); return;
    case 0b00000100: _g0_2d4d_0010(g, bc, envs); return;
    case 0b00000101: _g0_2d4d_0011(g, bc, envs); return;
    case 0b00000110: _g0_2d4d_0012(g, bc, envs); return;
    case 0b00001000: _g0_2d4d_0020(g, bc, envs); return;
    case 0b00001001: _g0_2d4d_0021(g, bc, envs); return;
    case 0b00001100: _g0_2d4d_0030(g, bc, envs); return;
    case 0b00010000: _g0_2d4d_0100(g, bc, envs); return;
    case 0b00010001: _g0_2d4d_0101(g, bc, envs); return;
    case 0b00010010: _g0_2d4d_0102(g, bc, envs); return;
    case 0b00010100: _g0_2d4d_0110(g, bc, envs); return;
    case 0b00010101: _g0_2d4d_0111(g, bc, envs); return;
    case 0b00011000: _g0_2d4d_0120(g, bc, envs); return;
    case 0b00100000: _g0_2d4d_0200(g, bc, envs); return;
    case 0b00100001: _g0_2d4d_0201(g, bc, envs); return;
    case 0b00100100: _g0_2d4d_0210(g, bc, envs); return;
    case 0b00110000: _g0_2d4d_0300(g, bc, envs); return;
    case 0b01000000: _g0_2d4d_1000(g, bc, envs); return;
    case 0b01000001: _g0_2d4d_1001(g, bc, envs); return;
    case 0b01000010: _g0_2d4d_1002(g, bc, envs); return;
    case 0b01000100: _g0_2d4d_1010(g, bc, envs); return;
    case 0b01000101: _g0_2d4d_1011(g, bc, envs); return;
    case 0b01001000: _g0_2d4d_1020(g, bc, envs); return;
    case 0b01010000: _g0_2d4d_1100(g, bc, envs); return;
    case 0b01010001: _g0_2d4d_1101(g, bc, envs); return;
    case 0b01010100: _g0_2d4d_1110(g, bc, envs); return;
    case 0b01100000: _g0_2d4d_1200(g, bc, envs); return;
    case 0b10000000: _g0_2d4d_2000(g, bc, envs); return;
    case 0b10000001: _g0_2d4d_2001(g, bc, envs); return;
    case 0b10000100: _g0_2d4d_2010(g, bc, envs); return;
    case 0b10010000: _g0_2d4d_2100(g, bc, envs); return;
    case 0b11000000: _g0_2d4d_3000(g, bc, envs); return;
    }
    fprintf(stderr, "Dimension error for CINTg0_2e_lj2d4d: iklj = %d %d %d %d",
        (int)envs->li_ceil, (int)envs->lk_ceil,
        (int)envs->ll_ceil, (int)envs->lj_ceil);
}


static inline void _srg0_2d4d_0000(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    //g[4] = w[0];
    //g[5] = w[0];
}
static inline void _srg0_2d4d_0001(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    g[0] = 1;
    g[1] = 1;
    g[2] = cpx[0];
    g[3] = cpx[1];
    g[4] = 1;
    g[5] = 1;
    g[6] = cpy[0];
    g[7] = cpy[1];
    //g[8] = w[0];
    //g[9] = w[0];
    g[10] = cpz[0] * g[8];
    g[11] = cpz[1] * g[9];
}
static inline void _srg0_2d4d_0002(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[2];
    g[7] = cpx[3];
    g[8] = cpx[0] * cpx[0] + b01[0];
    g[9] = cpx[1] * cpx[1] + b01[1];
    g[10] = cpx[2] * cpx[2] + b01[2];
    g[11] = cpx[3] * cpx[3] + b01[3];
    g[12] = 1;
    g[13] = 1;
    g[14] = 1;
    g[15] = 1;
    g[16] = cpy[0];
    g[17] = cpy[1];
    g[18] = cpy[2];
    g[19] = cpy[3];
    g[20] = cpy[0] * cpy[0] + b01[0];
    g[21] = cpy[1] * cpy[1] + b01[1];
    g[22] = cpy[2] * cpy[2] + b01[2];
    g[23] = cpy[3] * cpy[3] + b01[3];
    //g[24] = w[0];
    //g[25] = w[0];
    //g[26] = w[1];
    //g[27] = w[1];
    g[28] = cpz[0] * g[24];
    g[29] = cpz[1] * g[25];
    g[30] = cpz[2] * g[26];
    g[31] = cpz[3] * g[27];
    g[32] = cpz[0] * g[28] + b01[0] * g[24];
    g[33] = cpz[1] * g[29] + b01[1] * g[25];
    g[34] = cpz[2] * g[30] + b01[2] * g[26];
    g[35] = cpz[3] * g[31] + b01[3] * g[27];
}
static inline void _srg0_2d4d_0003(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[2];
    g[7] = cpx[3];
    g[8] = cpx[0] * cpx[0] + b01[0];
    g[9] = cpx[1] * cpx[1] + b01[1];
    g[10] = cpx[2] * cpx[2] + b01[2];
    g[11] = cpx[3] * cpx[3] + b01[3];
    g[12] = cpx[0] * (g[8] + 2 * b01[0]);
    g[13] = cpx[1] * (g[9] + 2 * b01[1]);
    g[14] = cpx[2] * (g[10] + 2 * b01[2]);
    g[15] = cpx[3] * (g[11] + 2 * b01[3]);
    g[16] = 1;
    g[17] = 1;
    g[18] = 1;
    g[19] = 1;
    g[20] = cpy[0];
    g[21] = cpy[1];
    g[22] = cpy[2];
    g[23] = cpy[3];
    g[24] = cpy[0] * cpy[0] + b01[0];
    g[25] = cpy[1] * cpy[1] + b01[1];
    g[26] = cpy[2] * cpy[2] + b01[2];
    g[27] = cpy[3] * cpy[3] + b01[3];
    g[28] = cpy[0] * (g[24] + 2 * b01[0]);
    g[29] = cpy[1] * (g[25] + 2 * b01[1]);
    g[30] = cpy[2] * (g[26] + 2 * b01[2]);
    g[31] = cpy[3] * (g[27] + 2 * b01[3]);
    //g[32] = w[0];
    //g[33] = w[0];
    //g[34] = w[1];
    //g[35] = w[1];
    g[36] = cpz[0] * g[32];
    g[37] = cpz[1] * g[33];
    g[38] = cpz[2] * g[34];
    g[39] = cpz[3] * g[35];
    g[40] = cpz[0] * g[36] + b01[0] * g[32];
    g[41] = cpz[1] * g[37] + b01[1] * g[33];
    g[42] = cpz[2] * g[38] + b01[2] * g[34];
    g[43] = cpz[3] * g[39] + b01[3] * g[35];
    g[44] = cpz[0] * g[40] + 2 * b01[0] * g[36];
    g[45] = cpz[1] * g[41] + 2 * b01[1] * g[37];
    g[46] = cpz[2] * g[42] + 2 * b01[2] * g[38];
    g[47] = cpz[3] * g[43] + 2 * b01[3] * g[39];
}
static inline void _srg0_2d4d_0010(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    g[0] = 1;
    g[1] = 1;
    g[2] = cpx[0];
    g[3] = cpx[1];
    g[4] = 1;
    g[5] = 1;
    g[6] = cpy[0];
    g[7] = cpy[1];
    //g[8] = w[0];
    //g[9] = w[0];
    g[10] = cpz[0] * g[8];
    g[11] = cpz[1] * g[9];
}
static inline void _srg0_2d4d_0011(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b01 = bc->b01;
    double xkxl = envs->rkrl[0];
    double ykyl = envs->rkrl[1];
    double zkzl = envs->rkrl[2];
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[8] = cpx[0];
    g[9] = cpx[1];
    g[10] = cpx[2];
    g[11] = cpx[3];
    g[12] = cpx[0] * (xkxl + cpx[0]) + b01[0];
    g[13] = cpx[1] * (xkxl + cpx[1]) + b01[1];
    g[14] = cpx[2] * (xkxl + cpx[2]) + b01[2];
    g[15] = cpx[3] * (xkxl + cpx[3]) + b01[3];
    g[4] = xkxl + cpx[0];
    g[5] = xkxl + cpx[1];
    g[6] = xkxl + cpx[2];
    g[7] = xkxl + cpx[3];
    g[24] = 1;
    g[25] = 1;
    g[26] = 1;
    g[27] = 1;
    g[32] = cpy[0];
    g[33] = cpy[1];
    g[34] = cpy[2];
    g[35] = cpy[3];
    g[36] = cpy[0] * (ykyl + cpy[0]) + b01[0];
    g[37] = cpy[1] * (ykyl + cpy[1]) + b01[1];
    g[38] = cpy[2] * (ykyl + cpy[2]) + b01[2];
    g[39] = cpy[3] * (ykyl + cpy[3]) + b01[3];
    g[28] = ykyl + cpy[0];
    g[29] = ykyl + cpy[1];
    g[30] = ykyl + cpy[2];
    g[31] = ykyl + cpy[3];
    //g[48] = w[0];
    //g[49] = w[0];
    //g[50] = w[1];
    //g[51] = w[1];
    g[56] = cpz[0] * g[48];
    g[57] = cpz[1] * g[49];
    g[58] = cpz[2] * g[50];
    g[59] = cpz[3] * g[51];
    g[60] = g[56] * (zkzl + cpz[0]) + b01[0] * g[48];
    g[61] = g[57] * (zkzl + cpz[1]) + b01[1] * g[49];
    g[62] = g[58] * (zkzl + cpz[2]) + b01[2] * g[50];
    g[63] = g[59] * (zkzl + cpz[3]) + b01[3] * g[51];
    g[52] = g[48] * (zkzl + cpz[0]);
    g[53] = g[49] * (zkzl + cpz[1]);
    g[54] = g[50] * (zkzl + cpz[2]);
    g[55] = g[51] * (zkzl + cpz[3]);
}
static inline void _srg0_2d4d_0012(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b01 = bc->b01;
    double xkxl = envs->rkrl[0];
    double ykyl = envs->rkrl[1];
    double zkzl = envs->rkrl[2];
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[8] = cpx[0];
    g[9] = cpx[1];
    g[10] = cpx[2];
    g[11] = cpx[3];
    g[16] = cpx[0] * cpx[0] + b01[0];
    g[17] = cpx[1] * cpx[1] + b01[1];
    g[18] = cpx[2] * cpx[2] + b01[2];
    g[19] = cpx[3] * cpx[3] + b01[3];
    g[20] = g[16] * (xkxl + cpx[0]) + cpx[0] * 2 * b01[0];
    g[21] = g[17] * (xkxl + cpx[1]) + cpx[1] * 2 * b01[1];
    g[22] = g[18] * (xkxl + cpx[2]) + cpx[2] * 2 * b01[2];
    g[23] = g[19] * (xkxl + cpx[3]) + cpx[3] * 2 * b01[3];
    g[12] = cpx[0] * (xkxl + cpx[0]) + b01[0];
    g[13] = cpx[1] * (xkxl + cpx[1]) + b01[1];
    g[14] = cpx[2] * (xkxl + cpx[2]) + b01[2];
    g[15] = cpx[3] * (xkxl + cpx[3]) + b01[3];
    g[4] = xkxl + cpx[0];
    g[5] = xkxl + cpx[1];
    g[6] = xkxl + cpx[2];
    g[7] = xkxl + cpx[3];
    g[32] = 1;
    g[33] = 1;
    g[34] = 1;
    g[35] = 1;
    g[40] = cpy[0];
    g[41] = cpy[1];
    g[42] = cpy[2];
    g[43] = cpy[3];
    g[48] = cpy[0] * cpy[0] + b01[0];
    g[49] = cpy[1] * cpy[1] + b01[1];
    g[50] = cpy[2] * cpy[2] + b01[2];
    g[51] = cpy[3] * cpy[3] + b01[3];
    g[52] = g[48] * (ykyl + cpy[0]) + cpy[0] * 2 * b01[0];
    g[53] = g[49] * (ykyl + cpy[1]) + cpy[1] * 2 * b01[1];
    g[54] = g[50] * (ykyl + cpy[2]) + cpy[2] * 2 * b01[2];
    g[55] = g[51] * (ykyl + cpy[3]) + cpy[3] * 2 * b01[3];
    g[44] = cpy[0] * (ykyl + cpy[0]) + b01[0];
    g[45] = cpy[1] * (ykyl + cpy[1]) + b01[1];
    g[46] = cpy[2] * (ykyl + cpy[2]) + b01[2];
    g[47] = cpy[3] * (ykyl + cpy[3]) + b01[3];
    g[36] = ykyl + cpy[0];
    g[37] = ykyl + cpy[1];
    g[38] = ykyl + cpy[2];
    g[39] = ykyl + cpy[3];
    //g[64] = w[0];
    //g[65] = w[0];
    //g[66] = w[1];
    //g[67] = w[1];
    g[72] = cpz[0] * g[64];
    g[73] = cpz[1] * g[65];
    g[74] = cpz[2] * g[66];
    g[75] = cpz[3] * g[67];
    g[80] = cpz[0] * g[72] + b01[0] * g[64];
    g[81] = cpz[1] * g[73] + b01[1] * g[65];
    g[82] = cpz[2] * g[74] + b01[2] * g[66];
    g[83] = cpz[3] * g[75] + b01[3] * g[67];
    g[84] = g[80] * (zkzl + cpz[0]) + 2 * b01[0] * g[72];
    g[85] = g[81] * (zkzl + cpz[1]) + 2 * b01[1] * g[73];
    g[86] = g[82] * (zkzl + cpz[2]) + 2 * b01[2] * g[74];
    g[87] = g[83] * (zkzl + cpz[3]) + 2 * b01[3] * g[75];
    g[76] = g[72] * (zkzl + cpz[0]) + b01[0] * g[64];
    g[77] = g[73] * (zkzl + cpz[1]) + b01[1] * g[65];
    g[78] = g[74] * (zkzl + cpz[2]) + b01[2] * g[66];
    g[79] = g[75] * (zkzl + cpz[3]) + b01[3] * g[67];
    g[68] = g[64] * (zkzl + cpz[0]);
    g[69] = g[65] * (zkzl + cpz[1]);
    g[70] = g[66] * (zkzl + cpz[2]);
    g[71] = g[67] * (zkzl + cpz[3]);
}
static inline void _srg0_2d4d_0020(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[2];
    g[7] = cpx[3];
    g[8] = cpx[0] * cpx[0] + b01[0];
    g[9] = cpx[1] * cpx[1] + b01[1];
    g[10] = cpx[2] * cpx[2] + b01[2];
    g[11] = cpx[3] * cpx[3] + b01[3];
    g[12] = 1;
    g[13] = 1;
    g[14] = 1;
    g[15] = 1;
    g[16] = cpy[0];
    g[17] = cpy[1];
    g[18] = cpy[2];
    g[19] = cpy[3];
    g[20] = cpy[0] * cpy[0] + b01[0];
    g[21] = cpy[1] * cpy[1] + b01[1];
    g[22] = cpy[2] * cpy[2] + b01[2];
    g[23] = cpy[3] * cpy[3] + b01[3];
    //g[24] = w[0];
    //g[25] = w[0];
    //g[26] = w[1];
    //g[27] = w[1];
    g[28] = cpz[0] * g[24];
    g[29] = cpz[1] * g[25];
    g[30] = cpz[2] * g[26];
    g[31] = cpz[3] * g[27];
    g[32] = cpz[0] * g[28] + b01[0] * g[24];
    g[33] = cpz[1] * g[29] + b01[1] * g[25];
    g[34] = cpz[2] * g[30] + b01[2] * g[26];
    g[35] = cpz[3] * g[31] + b01[3] * g[27];
}
static inline void _srg0_2d4d_0021(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b01 = bc->b01;
    double xkxl = envs->rkrl[0];
    double ykyl = envs->rkrl[1];
    double zkzl = envs->rkrl[2];
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[2];
    g[7] = cpx[3];
    g[8] = cpx[0] * cpx[0] + b01[0];
    g[9] = cpx[1] * cpx[1] + b01[1];
    g[10] = cpx[2] * cpx[2] + b01[2];
    g[11] = cpx[3] * cpx[3] + b01[3];
    g[16] = xkxl + cpx[0];
    g[17] = xkxl + cpx[1];
    g[18] = xkxl + cpx[2];
    g[19] = xkxl + cpx[3];
    g[20] = cpx[0] * (xkxl + cpx[0]) + b01[0];
    g[21] = cpx[1] * (xkxl + cpx[1]) + b01[1];
    g[22] = cpx[2] * (xkxl + cpx[2]) + b01[2];
    g[23] = cpx[3] * (xkxl + cpx[3]) + b01[3];
    g[24] = g[8] * (xkxl + cpx[0]) + cpx[0] * 2 * b01[0];
    g[25] = g[9] * (xkxl + cpx[1]) + cpx[1] * 2 * b01[1];
    g[26] = g[10] * (xkxl + cpx[2]) + cpx[2] * 2 * b01[2];
    g[27] = g[11] * (xkxl + cpx[3]) + cpx[3] * 2 * b01[3];
    g[32] = 1;
    g[33] = 1;
    g[34] = 1;
    g[35] = 1;
    g[36] = cpy[0];
    g[37] = cpy[1];
    g[38] = cpy[2];
    g[39] = cpy[3];
    g[40] = cpy[0] * cpy[0] + b01[0];
    g[41] = cpy[1] * cpy[1] + b01[1];
    g[42] = cpy[2] * cpy[2] + b01[2];
    g[43] = cpy[3] * cpy[3] + b01[3];
    g[48] = ykyl + cpy[0];
    g[49] = ykyl + cpy[1];
    g[50] = ykyl + cpy[2];
    g[51] = ykyl + cpy[3];
    g[52] = cpy[0] * (ykyl + cpy[0]) + b01[0];
    g[53] = cpy[1] * (ykyl + cpy[1]) + b01[1];
    g[54] = cpy[2] * (ykyl + cpy[2]) + b01[2];
    g[55] = cpy[3] * (ykyl + cpy[3]) + b01[3];
    g[56] = g[40] * (ykyl + cpy[0]) + cpy[0] * 2 * b01[0];
    g[57] = g[41] * (ykyl + cpy[1]) + cpy[1] * 2 * b01[1];
    g[58] = g[42] * (ykyl + cpy[2]) + cpy[2] * 2 * b01[2];
    g[59] = g[43] * (ykyl + cpy[3]) + cpy[3] * 2 * b01[3];
    //g[64] = w[0];
    //g[65] = w[0];
    //g[66] = w[1];
    //g[67] = w[1];
    g[68] = cpz[0] * g[64];
    g[69] = cpz[1] * g[65];
    g[70] = cpz[2] * g[66];
    g[71] = cpz[3] * g[67];
    g[72] = cpz[0] * g[68] + b01[0] * g[64];
    g[73] = cpz[1] * g[69] + b01[1] * g[65];
    g[74] = cpz[2] * g[70] + b01[2] * g[66];
    g[75] = cpz[3] * g[71] + b01[3] * g[67];
    g[80] = g[64] * (zkzl + cpz[0]);
    g[81] = g[65] * (zkzl + cpz[1]);
    g[82] = g[66] * (zkzl + cpz[2]);
    g[83] = g[67] * (zkzl + cpz[3]);
    g[84] = g[68] * (zkzl + cpz[0]) + b01[0] * g[64];
    g[85] = g[69] * (zkzl + cpz[1]) + b01[1] * g[65];
    g[86] = g[70] * (zkzl + cpz[2]) + b01[2] * g[66];
    g[87] = g[71] * (zkzl + cpz[3]) + b01[3] * g[67];
    g[88] = g[72] * (zkzl + cpz[0]) + 2 * b01[0] * g[68];
    g[89] = g[73] * (zkzl + cpz[1]) + 2 * b01[1] * g[69];
    g[90] = g[74] * (zkzl + cpz[2]) + 2 * b01[2] * g[70];
    g[91] = g[75] * (zkzl + cpz[3]) + 2 * b01[3] * g[71];
}
static inline void _srg0_2d4d_0030(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[2];
    g[7] = cpx[3];
    g[8] = cpx[0] * cpx[0] + b01[0];
    g[9] = cpx[1] * cpx[1] + b01[1];
    g[10] = cpx[2] * cpx[2] + b01[2];
    g[11] = cpx[3] * cpx[3] + b01[3];
    g[12] = cpx[0] * (g[8] + 2 * b01[0]);
    g[13] = cpx[1] * (g[9] + 2 * b01[1]);
    g[14] = cpx[2] * (g[10] + 2 * b01[2]);
    g[15] = cpx[3] * (g[11] + 2 * b01[3]);
    g[16] = 1;
    g[17] = 1;
    g[18] = 1;
    g[19] = 1;
    g[20] = cpy[0];
    g[21] = cpy[1];
    g[22] = cpy[2];
    g[23] = cpy[3];
    g[24] = cpy[0] * cpy[0] + b01[0];
    g[25] = cpy[1] * cpy[1] + b01[1];
    g[26] = cpy[2] * cpy[2] + b01[2];
    g[27] = cpy[3] * cpy[3] + b01[3];
    g[28] = cpy[0] * (g[24] + 2 * b01[0]);
    g[29] = cpy[1] * (g[25] + 2 * b01[1]);
    g[30] = cpy[2] * (g[26] + 2 * b01[2]);
    g[31] = cpy[3] * (g[27] + 2 * b01[3]);
    //g[32] = w[0];
    //g[33] = w[0];
    //g[34] = w[1];
    //g[35] = w[1];
    g[36] = cpz[0] * g[32];
    g[37] = cpz[1] * g[33];
    g[38] = cpz[2] * g[34];
    g[39] = cpz[3] * g[35];
    g[40] = cpz[0] * g[36] + b01[0] * g[32];
    g[41] = cpz[1] * g[37] + b01[1] * g[33];
    g[42] = cpz[2] * g[38] + b01[2] * g[34];
    g[43] = cpz[3] * g[39] + b01[3] * g[35];
    g[44] = cpz[0] * g[40] + 2 * b01[0] * g[36];
    g[45] = cpz[1] * g[41] + 2 * b01[1] * g[37];
    g[46] = cpz[2] * g[42] + 2 * b01[2] * g[38];
    g[47] = cpz[3] * g[43] + 2 * b01[3] * g[39];
}
static inline void _srg0_2d4d_0100(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    g[0] = 1;
    g[1] = 1;
    g[2] = c0x[0];
    g[3] = c0x[1];
    g[4] = 1;
    g[5] = 1;
    g[6] = c0y[0];
    g[7] = c0y[1];
    //g[8] = w[0];
    //g[9] = w[0];
    g[10] = c0z[0] * g[8];
    g[11] = c0z[1] * g[9];
}
static inline void _srg0_2d4d_0101(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[2];
    g[7] = cpx[3];
    g[8] = c0x[0];
    g[9] = c0x[1];
    g[10] = c0x[2];
    g[11] = c0x[3];
    g[12] = cpx[0] * c0x[0] + b00[0];
    g[13] = cpx[1] * c0x[1] + b00[1];
    g[14] = cpx[2] * c0x[2] + b00[2];
    g[15] = cpx[3] * c0x[3] + b00[3];
    g[16] = 1;
    g[17] = 1;
    g[18] = 1;
    g[19] = 1;
    g[20] = cpy[0];
    g[21] = cpy[1];
    g[22] = cpy[2];
    g[23] = cpy[3];
    g[24] = c0y[0];
    g[25] = c0y[1];
    g[26] = c0y[2];
    g[27] = c0y[3];
    g[28] = cpy[0] * c0y[0] + b00[0];
    g[29] = cpy[1] * c0y[1] + b00[1];
    g[30] = cpy[2] * c0y[2] + b00[2];
    g[31] = cpy[3] * c0y[3] + b00[3];
    //g[32] = w[0];
    //g[33] = w[0];
    //g[34] = w[1];
    //g[35] = w[1];
    g[36] = cpz[0] * g[32];
    g[37] = cpz[1] * g[33];
    g[38] = cpz[2] * g[34];
    g[39] = cpz[3] * g[35];
    g[40] = c0z[0] * g[32];
    g[41] = c0z[1] * g[33];
    g[42] = c0z[2] * g[34];
    g[43] = c0z[3] * g[35];
    g[44] = cpz[0] * g[40] + b00[0] * g[32];
    g[45] = cpz[1] * g[41] + b00[1] * g[33];
    g[46] = cpz[2] * g[42] + b00[2] * g[34];
    g[47] = cpz[3] * g[43] + b00[3] * g[35];
}
static inline void _srg0_2d4d_0102(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    double* b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[2];
    g[7] = cpx[3];
    g[12] = c0x[0];
    g[13] = c0x[1];
    g[14] = c0x[2];
    g[15] = c0x[3];
    g[8] = cpx[0] * cpx[0] + b01[0];
    g[9] = cpx[1] * cpx[1] + b01[1];
    g[10] = cpx[2] * cpx[2] + b01[2];
    g[11] = cpx[3] * cpx[3] + b01[3];
    g[16] = cpx[0] * c0x[0] + b00[0];
    g[17] = cpx[1] * c0x[1] + b00[1];
    g[18] = cpx[2] * c0x[2] + b00[2];
    g[19] = cpx[3] * c0x[3] + b00[3];
    g[20] = cpx[0] * (g[16] + b00[0]) + b01[0] * c0x[0];
    g[21] = cpx[1] * (g[17] + b00[1]) + b01[1] * c0x[1];
    g[22] = cpx[2] * (g[18] + b00[2]) + b01[2] * c0x[2];
    g[23] = cpx[3] * (g[19] + b00[3]) + b01[3] * c0x[3];
    g[24] = 1;
    g[25] = 1;
    g[26] = 1;
    g[27] = 1;
    g[28] = cpy[0];
    g[29] = cpy[1];
    g[30] = cpy[2];
    g[31] = cpy[3];
    g[36] = c0y[0];
    g[37] = c0y[1];
    g[38] = c0y[2];
    g[39] = c0y[3];
    g[32] = cpy[0] * cpy[0] + b01[0];
    g[33] = cpy[1] * cpy[1] + b01[1];
    g[34] = cpy[2] * cpy[2] + b01[2];
    g[35] = cpy[3] * cpy[3] + b01[3];
    g[40] = cpy[0] * c0y[0] + b00[0];
    g[41] = cpy[1] * c0y[1] + b00[1];
    g[42] = cpy[2] * c0y[2] + b00[2];
    g[43] = cpy[3] * c0y[3] + b00[3];
    g[44] = cpy[0] * (g[40] + b00[0]) + b01[0] * c0y[0];
    g[45] = cpy[1] * (g[41] + b00[1]) + b01[1] * c0y[1];
    g[46] = cpy[2] * (g[42] + b00[2]) + b01[2] * c0y[2];
    g[47] = cpy[3] * (g[43] + b00[3]) + b01[3] * c0y[3];
    //g[48] = w[0];
    //g[49] = w[0];
    //g[50] = w[1];
    //g[51] = w[1];
    g[52] = cpz[0] * g[48];
    g[53] = cpz[1] * g[49];
    g[54] = cpz[2] * g[50];
    g[55] = cpz[3] * g[51];
    g[60] = c0z[0] * g[48];
    g[61] = c0z[1] * g[49];
    g[62] = c0z[2] * g[50];
    g[63] = c0z[3] * g[51];
    g[56] = cpz[0] * g[52] + b01[0] * g[48];
    g[57] = cpz[1] * g[53] + b01[1] * g[49];
    g[58] = cpz[2] * g[54] + b01[2] * g[50];
    g[59] = cpz[3] * g[55] + b01[3] * g[51];
    g[64] = cpz[0] * g[60] + b00[0] * g[48];
    g[65] = cpz[1] * g[61] + b00[1] * g[49];
    g[66] = cpz[2] * g[62] + b00[2] * g[50];
    g[67] = cpz[3] * g[63] + b00[3] * g[51];
    g[68] = cpz[0] * g[64] + b01[0] * g[60] + b00[0] * g[52];
    g[69] = cpz[1] * g[65] + b01[1] * g[61] + b00[1] * g[53];
    g[70] = cpz[2] * g[66] + b01[2] * g[62] + b00[2] * g[54];
    g[71] = cpz[3] * g[67] + b01[3] * g[63] + b00[3] * g[55];
}
static inline void _srg0_2d4d_0110(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[2];
    g[7] = cpx[3];
    g[8] = c0x[0];
    g[9] = c0x[1];
    g[10] = c0x[2];
    g[11] = c0x[3];
    g[12] = cpx[0] * c0x[0] + b00[0];
    g[13] = cpx[1] * c0x[1] + b00[1];
    g[14] = cpx[2] * c0x[2] + b00[2];
    g[15] = cpx[3] * c0x[3] + b00[3];
    g[16] = 1;
    g[17] = 1;
    g[18] = 1;
    g[19] = 1;
    g[20] = cpy[0];
    g[21] = cpy[1];
    g[22] = cpy[2];
    g[23] = cpy[3];
    g[24] = c0y[0];
    g[25] = c0y[1];
    g[26] = c0y[2];
    g[27] = c0y[3];
    g[28] = cpy[0] * c0y[0] + b00[0];
    g[29] = cpy[1] * c0y[1] + b00[1];
    g[30] = cpy[2] * c0y[2] + b00[2];
    g[31] = cpy[3] * c0y[3] + b00[3];
    //g[32] = w[0];
    //g[33] = w[0];
    //g[34] = w[1];
    //g[35] = w[1];
    g[36] = cpz[0] * g[32];
    g[37] = cpz[1] * g[33];
    g[38] = cpz[2] * g[34];
    g[39] = cpz[3] * g[35];
    g[40] = c0z[0] * g[32];
    g[41] = c0z[1] * g[33];
    g[42] = c0z[2] * g[34];
    g[43] = c0z[3] * g[35];
    g[44] = cpz[0] * g[40] + b00[0] * g[32];
    g[45] = cpz[1] * g[41] + b00[1] * g[33];
    g[46] = cpz[2] * g[42] + b00[2] * g[34];
    g[47] = cpz[3] * g[43] + b00[3] * g[35];
}
static inline void _srg0_2d4d_0111(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    double* b01 = bc->b01;
    double xkxl = envs->rkrl[0];
    double ykyl = envs->rkrl[1];
    double zkzl = envs->rkrl[2];
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[24] = c0x[0];
    g[25] = c0x[1];
    g[26] = c0x[2];
    g[27] = c0x[3];
    g[8] = cpx[0];
    g[9] = cpx[1];
    g[10] = cpx[2];
    g[11] = cpx[3];
    g[32] = cpx[0] * c0x[0] + b00[0];
    g[33] = cpx[1] * c0x[1] + b00[1];
    g[34] = cpx[2] * c0x[2] + b00[2];
    g[35] = cpx[3] * c0x[3] + b00[3];
    g[12] = cpx[0] * (xkxl + cpx[0]) + b01[0];
    g[13] = cpx[1] * (xkxl + cpx[1]) + b01[1];
    g[14] = cpx[2] * (xkxl + cpx[2]) + b01[2];
    g[15] = cpx[3] * (xkxl + cpx[3]) + b01[3];
    g[36] = g[32] * (xkxl + cpx[0]) + cpx[0] * b00[0] + b01[0] * c0x[0];
    g[37] = g[33] * (xkxl + cpx[1]) + cpx[1] * b00[1] + b01[1] * c0x[1];
    g[38] = g[34] * (xkxl + cpx[2]) + cpx[2] * b00[2] + b01[2] * c0x[2];
    g[39] = g[35] * (xkxl + cpx[3]) + cpx[3] * b00[3] + b01[3] * c0x[3];
    g[4] = xkxl + cpx[0];
    g[5] = xkxl + cpx[1];
    g[6] = xkxl + cpx[2];
    g[7] = xkxl + cpx[3];
    g[28] = c0x[0] * (xkxl + cpx[0]) + b00[0];
    g[29] = c0x[1] * (xkxl + cpx[1]) + b00[1];
    g[30] = c0x[2] * (xkxl + cpx[2]) + b00[2];
    g[31] = c0x[3] * (xkxl + cpx[3]) + b00[3];
    g[48] = 1;
    g[49] = 1;
    g[50] = 1;
    g[51] = 1;
    g[72] = c0y[0];
    g[73] = c0y[1];
    g[74] = c0y[2];
    g[75] = c0y[3];
    g[56] = cpy[0];
    g[57] = cpy[1];
    g[58] = cpy[2];
    g[59] = cpy[3];
    g[80] = cpy[0] * c0y[0] + b00[0];
    g[81] = cpy[1] * c0y[1] + b00[1];
    g[82] = cpy[2] * c0y[2] + b00[2];
    g[83] = cpy[3] * c0y[3] + b00[3];
    g[60] = cpy[0] * (ykyl + cpy[0]) + b01[0];
    g[61] = cpy[1] * (ykyl + cpy[1]) + b01[1];
    g[62] = cpy[2] * (ykyl + cpy[2]) + b01[2];
    g[63] = cpy[3] * (ykyl + cpy[3]) + b01[3];
    g[84] = g[80] * (ykyl + cpy[0]) + cpy[0] * b00[0] + b01[0] * c0y[0];
    g[85] = g[81] * (ykyl + cpy[1]) + cpy[1] * b00[1] + b01[1] * c0y[1];
    g[86] = g[82] * (ykyl + cpy[2]) + cpy[2] * b00[2] + b01[2] * c0y[2];
    g[87] = g[83] * (ykyl + cpy[3]) + cpy[3] * b00[3] + b01[3] * c0y[3];
    g[52] = ykyl + cpy[0];
    g[53] = ykyl + cpy[1];
    g[54] = ykyl + cpy[2];
    g[55] = ykyl + cpy[3];
    g[76] = c0y[0] * (ykyl + cpy[0]) + b00[0];
    g[77] = c0y[1] * (ykyl + cpy[1]) + b00[1];
    g[78] = c0y[2] * (ykyl + cpy[2]) + b00[2];
    g[79] = c0y[3] * (ykyl + cpy[3]) + b00[3];
    //g[96] = w[0];
    //g[97] = w[0];
    //g[98] = w[1];
    //g[99] = w[1];
    g[120] = c0z[0] * g[96];
    g[121] = c0z[1] * g[97];
    g[122] = c0z[2] * g[98];
    g[123] = c0z[3] * g[99];
    g[104] = cpz[0] * g[96];
    g[105] = cpz[1] * g[97];
    g[106] = cpz[2] * g[98];
    g[107] = cpz[3] * g[99];
    g[128] = cpz[0] * g[120] + b00[0] * g[96];
    g[129] = cpz[1] * g[121] + b00[1] * g[97];
    g[130] = cpz[2] * g[122] + b00[2] * g[98];
    g[131] = cpz[3] * g[123] + b00[3] * g[99];
    g[108] = g[104] * (zkzl + cpz[0]) + b01[0] * g[96];
    g[109] = g[105] * (zkzl + cpz[1]) + b01[1] * g[97];
    g[110] = g[106] * (zkzl + cpz[2]) + b01[2] * g[98];
    g[111] = g[107] * (zkzl + cpz[3]) + b01[3] * g[99];
    g[132] = g[128] * (zkzl + cpz[0]) + b01[0] * g[120] + b00[0] * g[104];
    g[133] = g[129] * (zkzl + cpz[1]) + b01[1] * g[121] + b00[1] * g[105];
    g[134] = g[130] * (zkzl + cpz[2]) + b01[2] * g[122] + b00[2] * g[106];
    g[135] = g[131] * (zkzl + cpz[3]) + b01[3] * g[123] + b00[3] * g[107];
    g[100] = g[96] * (zkzl + cpz[0]);
    g[101] = g[97] * (zkzl + cpz[1]);
    g[102] = g[98] * (zkzl + cpz[2]);
    g[103] = g[99] * (zkzl + cpz[3]);
    g[124] = g[120] * (zkzl + cpz[0]) + b00[0] * g[96];
    g[125] = g[121] * (zkzl + cpz[1]) + b00[1] * g[97];
    g[126] = g[122] * (zkzl + cpz[2]) + b00[2] * g[98];
    g[127] = g[123] * (zkzl + cpz[3]) + b00[3] * g[99];
}
static inline void _srg0_2d4d_0120(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    double* b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[2];
    g[7] = cpx[3];
    g[12] = c0x[0];
    g[13] = c0x[1];
    g[14] = c0x[2];
    g[15] = c0x[3];
    g[8] = cpx[0] * cpx[0] + b01[0];
    g[9] = cpx[1] * cpx[1] + b01[1];
    g[10] = cpx[2] * cpx[2] + b01[2];
    g[11] = cpx[3] * cpx[3] + b01[3];
    g[16] = cpx[0] * c0x[0] + b00[0];
    g[17] = cpx[1] * c0x[1] + b00[1];
    g[18] = cpx[2] * c0x[2] + b00[2];
    g[19] = cpx[3] * c0x[3] + b00[3];
    g[20] = cpx[0] * (g[16] + b00[0]) + b01[0] * c0x[0];
    g[21] = cpx[1] * (g[17] + b00[1]) + b01[1] * c0x[1];
    g[22] = cpx[2] * (g[18] + b00[2]) + b01[2] * c0x[2];
    g[23] = cpx[3] * (g[19] + b00[3]) + b01[3] * c0x[3];
    g[24] = 1;
    g[25] = 1;
    g[26] = 1;
    g[27] = 1;
    g[28] = cpy[0];
    g[29] = cpy[1];
    g[30] = cpy[2];
    g[31] = cpy[3];
    g[36] = c0y[0];
    g[37] = c0y[1];
    g[38] = c0y[2];
    g[39] = c0y[3];
    g[32] = cpy[0] * cpy[0] + b01[0];
    g[33] = cpy[1] * cpy[1] + b01[1];
    g[34] = cpy[2] * cpy[2] + b01[2];
    g[35] = cpy[3] * cpy[3] + b01[3];
    g[40] = cpy[0] * c0y[0] + b00[0];
    g[41] = cpy[1] * c0y[1] + b00[1];
    g[42] = cpy[2] * c0y[2] + b00[2];
    g[43] = cpy[3] * c0y[3] + b00[3];
    g[44] = cpy[0] * (g[40] + b00[0]) + b01[0] * c0y[0];
    g[45] = cpy[1] * (g[41] + b00[1]) + b01[1] * c0y[1];
    g[46] = cpy[2] * (g[42] + b00[2]) + b01[2] * c0y[2];
    g[47] = cpy[3] * (g[43] + b00[3]) + b01[3] * c0y[3];
    //g[48] = w[0];
    //g[49] = w[0];
    //g[50] = w[1];
    //g[51] = w[1];
    g[52] = cpz[0] * g[48];
    g[53] = cpz[1] * g[49];
    g[54] = cpz[2] * g[50];
    g[55] = cpz[3] * g[51];
    g[60] = c0z[0] * g[48];
    g[61] = c0z[1] * g[49];
    g[62] = c0z[2] * g[50];
    g[63] = c0z[3] * g[51];
    g[56] = cpz[0] * g[52] + b01[0] * g[48];
    g[57] = cpz[1] * g[53] + b01[1] * g[49];
    g[58] = cpz[2] * g[54] + b01[2] * g[50];
    g[59] = cpz[3] * g[55] + b01[3] * g[51];
    g[64] = cpz[0] * g[60] + b00[0] * g[48];
    g[65] = cpz[1] * g[61] + b00[1] * g[49];
    g[66] = cpz[2] * g[62] + b00[2] * g[50];
    g[67] = cpz[3] * g[63] + b00[3] * g[51];
    g[68] = cpz[0] * g[64] + b01[0] * g[60] + b00[0] * g[52];
    g[69] = cpz[1] * g[65] + b01[1] * g[61] + b00[1] * g[53];
    g[70] = cpz[2] * g[66] + b01[2] * g[62] + b00[2] * g[54];
    g[71] = cpz[3] * g[67] + b01[3] * g[63] + b00[3] * g[55];
}
static inline void _srg0_2d4d_0200(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = c0x[2];
    g[7] = c0x[3];
    g[8] = c0x[0] * c0x[0] + b10[0];
    g[9] = c0x[1] * c0x[1] + b10[1];
    g[10] = c0x[2] * c0x[2] + b10[2];
    g[11] = c0x[3] * c0x[3] + b10[3];
    g[12] = 1;
    g[13] = 1;
    g[14] = 1;
    g[15] = 1;
    g[16] = c0y[0];
    g[17] = c0y[1];
    g[18] = c0y[2];
    g[19] = c0y[3];
    g[20] = c0y[0] * c0y[0] + b10[0];
    g[21] = c0y[1] * c0y[1] + b10[1];
    g[22] = c0y[2] * c0y[2] + b10[2];
    g[23] = c0y[3] * c0y[3] + b10[3];
    //g[24] = w[0];
    //g[25] = w[0];
    //g[26] = w[1];
    //g[27] = w[1];
    g[28] = c0z[0] * g[24];
    g[29] = c0z[1] * g[25];
    g[30] = c0z[2] * g[26];
    g[31] = c0z[3] * g[27];
    g[32] = c0z[0] * g[28] + b10[0] * g[24];
    g[33] = c0z[1] * g[29] + b10[1] * g[25];
    g[34] = c0z[2] * g[30] + b10[2] * g[26];
    g[35] = c0z[3] * g[31] + b10[3] * g[27];
}
static inline void _srg0_2d4d_0201(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    double* b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[8] = c0x[0];
    g[9] = c0x[1];
    g[10] = c0x[2];
    g[11] = c0x[3];
    g[16] = c0x[0] * c0x[0] + b10[0];
    g[17] = c0x[1] * c0x[1] + b10[1];
    g[18] = c0x[2] * c0x[2] + b10[2];
    g[19] = c0x[3] * c0x[3] + b10[3];
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[2];
    g[7] = cpx[3];
    g[12] = cpx[0] * c0x[0] + b00[0];
    g[13] = cpx[1] * c0x[1] + b00[1];
    g[14] = cpx[2] * c0x[2] + b00[2];
    g[15] = cpx[3] * c0x[3] + b00[3];
    g[20] = c0x[0] * (g[12] + b00[0]) + b10[0] * cpx[0];
    g[21] = c0x[1] * (g[13] + b00[1]) + b10[1] * cpx[1];
    g[22] = c0x[2] * (g[14] + b00[2]) + b10[2] * cpx[2];
    g[23] = c0x[3] * (g[15] + b00[3]) + b10[3] * cpx[3];
    g[24] = 1;
    g[25] = 1;
    g[26] = 1;
    g[27] = 1;
    g[32] = c0y[0];
    g[33] = c0y[1];
    g[34] = c0y[2];
    g[35] = c0y[3];
    g[40] = c0y[0] * c0y[0] + b10[0];
    g[41] = c0y[1] * c0y[1] + b10[1];
    g[42] = c0y[2] * c0y[2] + b10[2];
    g[43] = c0y[3] * c0y[3] + b10[3];
    g[28] = cpy[0];
    g[29] = cpy[1];
    g[30] = cpy[2];
    g[31] = cpy[3];
    g[36] = cpy[0] * c0y[0] + b00[0];
    g[37] = cpy[1] * c0y[1] + b00[1];
    g[38] = cpy[2] * c0y[2] + b00[2];
    g[39] = cpy[3] * c0y[3] + b00[3];
    g[44] = c0y[0] * (g[36] + b00[0]) + b10[0] * cpy[0];
    g[45] = c0y[1] * (g[37] + b00[1]) + b10[1] * cpy[1];
    g[46] = c0y[2] * (g[38] + b00[2]) + b10[2] * cpy[2];
    g[47] = c0y[3] * (g[39] + b00[3]) + b10[3] * cpy[3];
    //g[48] = w[0];
    //g[49] = w[0];
    //g[50] = w[1];
    //g[51] = w[1];
    g[56] = c0z[0] * g[48];
    g[57] = c0z[1] * g[49];
    g[58] = c0z[2] * g[50];
    g[59] = c0z[3] * g[51];
    g[64] = c0z[0] * g[56] + b10[0] * g[48];
    g[65] = c0z[1] * g[57] + b10[1] * g[49];
    g[66] = c0z[2] * g[58] + b10[2] * g[50];
    g[67] = c0z[3] * g[59] + b10[3] * g[51];
    g[52] = cpz[0] * g[48];
    g[53] = cpz[1] * g[49];
    g[54] = cpz[2] * g[50];
    g[55] = cpz[3] * g[51];
    g[60] = cpz[0] * g[56] + b00[0] * g[48];
    g[61] = cpz[1] * g[57] + b00[1] * g[49];
    g[62] = cpz[2] * g[58] + b00[2] * g[50];
    g[63] = cpz[3] * g[59] + b00[3] * g[51];
    g[68] = c0z[0] * g[60] + b10[0] * g[52] + b00[0] * g[56];
    g[69] = c0z[1] * g[61] + b10[1] * g[53] + b00[1] * g[57];
    g[70] = c0z[2] * g[62] + b10[2] * g[54] + b00[2] * g[58];
    g[71] = c0z[3] * g[63] + b10[3] * g[55] + b00[3] * g[59];
}
static inline void _srg0_2d4d_0210(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    double* b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = cpx[0];
    g[5] = cpx[1];
    g[6] = cpx[2];
    g[7] = cpx[3];
    g[8] = c0x[0];
    g[9] = c0x[1];
    g[10] = c0x[2];
    g[11] = c0x[3];
    g[12] = cpx[0] * c0x[0] + b00[0];
    g[13] = cpx[1] * c0x[1] + b00[1];
    g[14] = cpx[2] * c0x[2] + b00[2];
    g[15] = cpx[3] * c0x[3] + b00[3];
    g[16] = c0x[0] * c0x[0] + b10[0];
    g[17] = c0x[1] * c0x[1] + b10[1];
    g[18] = c0x[2] * c0x[2] + b10[2];
    g[19] = c0x[3] * c0x[3] + b10[3];
    g[20] = c0x[0] * (g[12] + b00[0]) + b10[0] * cpx[0];
    g[21] = c0x[1] * (g[13] + b00[1]) + b10[1] * cpx[1];
    g[22] = c0x[2] * (g[14] + b00[2]) + b10[2] * cpx[2];
    g[23] = c0x[3] * (g[15] + b00[3]) + b10[3] * cpx[3];
    g[24] = 1;
    g[25] = 1;
    g[26] = 1;
    g[27] = 1;
    g[28] = cpy[0];
    g[29] = cpy[1];
    g[30] = cpy[2];
    g[31] = cpy[3];
    g[32] = c0y[0];
    g[33] = c0y[1];
    g[34] = c0y[2];
    g[35] = c0y[3];
    g[36] = cpy[0] * c0y[0] + b00[0];
    g[37] = cpy[1] * c0y[1] + b00[1];
    g[38] = cpy[2] * c0y[2] + b00[2];
    g[39] = cpy[3] * c0y[3] + b00[3];
    g[40] = c0y[0] * c0y[0] + b10[0];
    g[41] = c0y[1] * c0y[1] + b10[1];
    g[42] = c0y[2] * c0y[2] + b10[2];
    g[43] = c0y[3] * c0y[3] + b10[3];
    g[44] = c0y[0] * (g[36] + b00[0]) + b10[0] * cpy[0];
    g[45] = c0y[1] * (g[37] + b00[1]) + b10[1] * cpy[1];
    g[46] = c0y[2] * (g[38] + b00[2]) + b10[2] * cpy[2];
    g[47] = c0y[3] * (g[39] + b00[3]) + b10[3] * cpy[3];
    //g[48] = w[0];
    //g[49] = w[0];
    //g[50] = w[1];
    //g[51] = w[1];
    g[52] = cpz[0] * g[48];
    g[53] = cpz[1] * g[49];
    g[54] = cpz[2] * g[50];
    g[55] = cpz[3] * g[51];
    g[56] = c0z[0] * g[48];
    g[57] = c0z[1] * g[49];
    g[58] = c0z[2] * g[50];
    g[59] = c0z[3] * g[51];
    g[60] = cpz[0] * g[56] + b00[0] * g[48];
    g[61] = cpz[1] * g[57] + b00[1] * g[49];
    g[62] = cpz[2] * g[58] + b00[2] * g[50];
    g[63] = cpz[3] * g[59] + b00[3] * g[51];
    g[64] = c0z[0] * g[56] + b10[0] * g[48];
    g[65] = c0z[1] * g[57] + b10[1] * g[49];
    g[66] = c0z[2] * g[58] + b10[2] * g[50];
    g[67] = c0z[3] * g[59] + b10[3] * g[51];
    g[68] = c0z[0] * g[60] + b10[0] * g[52] + b00[0] * g[56];
    g[69] = c0z[1] * g[61] + b10[1] * g[53] + b00[1] * g[57];
    g[70] = c0z[2] * g[62] + b10[2] * g[54] + b00[2] * g[58];
    g[71] = c0z[3] * g[63] + b10[3] * g[55] + b00[3] * g[59];
}
static inline void _srg0_2d4d_0300(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = c0x[2];
    g[7] = c0x[3];
    g[8] = c0x[0] * c0x[0] + b10[0];
    g[9] = c0x[1] * c0x[1] + b10[1];
    g[10] = c0x[2] * c0x[2] + b10[2];
    g[11] = c0x[3] * c0x[3] + b10[3];
    g[12] = c0x[0] * (g[8] + 2 * b10[0]);
    g[13] = c0x[1] * (g[9] + 2 * b10[1]);
    g[14] = c0x[2] * (g[10] + 2 * b10[2]);
    g[15] = c0x[3] * (g[11] + 2 * b10[3]);
    g[16] = 1;
    g[17] = 1;
    g[18] = 1;
    g[19] = 1;
    g[20] = c0y[0];
    g[21] = c0y[1];
    g[22] = c0y[2];
    g[23] = c0y[3];
    g[24] = c0y[0] * c0y[0] + b10[0];
    g[25] = c0y[1] * c0y[1] + b10[1];
    g[26] = c0y[2] * c0y[2] + b10[2];
    g[27] = c0y[3] * c0y[3] + b10[3];
    g[28] = c0y[0] * (g[24] + 2 * b10[0]);
    g[29] = c0y[1] * (g[25] + 2 * b10[1]);
    g[30] = c0y[2] * (g[26] + 2 * b10[2]);
    g[31] = c0y[3] * (g[27] + 2 * b10[3]);
    //g[32] = w[0];
    //g[33] = w[0];
    //g[34] = w[1];
    //g[35] = w[1];
    g[36] = c0z[0] * g[32];
    g[37] = c0z[1] * g[33];
    g[38] = c0z[2] * g[34];
    g[39] = c0z[3] * g[35];
    g[40] = c0z[0] * g[36] + b10[0] * g[32];
    g[41] = c0z[1] * g[37] + b10[1] * g[33];
    g[42] = c0z[2] * g[38] + b10[2] * g[34];
    g[43] = c0z[3] * g[39] + b10[3] * g[35];
    g[44] = c0z[0] * g[40] + 2 * b10[0] * g[36];
    g[45] = c0z[1] * g[41] + 2 * b10[1] * g[37];
    g[46] = c0z[2] * g[42] + 2 * b10[2] * g[38];
    g[47] = c0z[3] * g[43] + 2 * b10[3] * g[39];
}
static inline void _srg0_2d4d_1000(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    g[0] = 1;
    g[1] = 1;
    g[2] = c0x[0];
    g[3] = c0x[1];
    g[4] = 1;
    g[5] = 1;
    g[6] = c0y[0];
    g[7] = c0y[1];
    //g[8] = w[0];
    //g[9] = w[0];
    g[10] = c0z[0] * g[8];
    g[11] = c0z[1] * g[9];
}
static inline void _srg0_2d4d_1001(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = c0x[2];
    g[7] = c0x[3];
    g[8] = cpx[0];
    g[9] = cpx[1];
    g[10] = cpx[2];
    g[11] = cpx[3];
    g[12] = cpx[0] * c0x[0] + b00[0];
    g[13] = cpx[1] * c0x[1] + b00[1];
    g[14] = cpx[2] * c0x[2] + b00[2];
    g[15] = cpx[3] * c0x[3] + b00[3];
    g[16] = 1;
    g[17] = 1;
    g[18] = 1;
    g[19] = 1;
    g[20] = c0y[0];
    g[21] = c0y[1];
    g[22] = c0y[2];
    g[23] = c0y[3];
    g[24] = cpy[0];
    g[25] = cpy[1];
    g[26] = cpy[2];
    g[27] = cpy[3];
    g[28] = cpy[0] * c0y[0] + b00[0];
    g[29] = cpy[1] * c0y[1] + b00[1];
    g[30] = cpy[2] * c0y[2] + b00[2];
    g[31] = cpy[3] * c0y[3] + b00[3];
    //g[32] = w[0];
    //g[33] = w[0];
    //g[34] = w[1];
    //g[35] = w[1];
    g[36] = c0z[0] * g[32];
    g[37] = c0z[1] * g[33];
    g[38] = c0z[2] * g[34];
    g[39] = c0z[3] * g[35];
    g[40] = cpz[0] * g[32];
    g[41] = cpz[1] * g[33];
    g[42] = cpz[2] * g[34];
    g[43] = cpz[3] * g[35];
    g[44] = cpz[0] * g[36] + b00[0] * g[32];
    g[45] = cpz[1] * g[37] + b00[1] * g[33];
    g[46] = cpz[2] * g[38] + b00[2] * g[34];
    g[47] = cpz[3] * g[39] + b00[3] * g[35];
}
static inline void _srg0_2d4d_1002(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    double* b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = c0x[2];
    g[7] = c0x[3];
    g[8] = cpx[0];
    g[9] = cpx[1];
    g[10] = cpx[2];
    g[11] = cpx[3];
    g[12] = cpx[0] * c0x[0] + b00[0];
    g[13] = cpx[1] * c0x[1] + b00[1];
    g[14] = cpx[2] * c0x[2] + b00[2];
    g[15] = cpx[3] * c0x[3] + b00[3];
    g[16] = cpx[0] * cpx[0] + b01[0];
    g[17] = cpx[1] * cpx[1] + b01[1];
    g[18] = cpx[2] * cpx[2] + b01[2];
    g[19] = cpx[3] * cpx[3] + b01[3];
    g[20] = cpx[0] * (g[12] + b00[0]) + b01[0] * c0x[0];
    g[21] = cpx[1] * (g[13] + b00[1]) + b01[1] * c0x[1];
    g[22] = cpx[2] * (g[14] + b00[2]) + b01[2] * c0x[2];
    g[23] = cpx[3] * (g[15] + b00[3]) + b01[3] * c0x[3];
    g[24] = 1;
    g[25] = 1;
    g[26] = 1;
    g[27] = 1;
    g[28] = c0y[0];
    g[29] = c0y[1];
    g[30] = c0y[2];
    g[31] = c0y[3];
    g[32] = cpy[0];
    g[33] = cpy[1];
    g[34] = cpy[2];
    g[35] = cpy[3];
    g[36] = cpy[0] * c0y[0] + b00[0];
    g[37] = cpy[1] * c0y[1] + b00[1];
    g[38] = cpy[2] * c0y[2] + b00[2];
    g[39] = cpy[3] * c0y[3] + b00[3];
    g[40] = cpy[0] * cpy[0] + b01[0];
    g[41] = cpy[1] * cpy[1] + b01[1];
    g[42] = cpy[2] * cpy[2] + b01[2];
    g[43] = cpy[3] * cpy[3] + b01[3];
    g[44] = cpy[0] * (g[36] + b00[0]) + b01[0] * c0y[0];
    g[45] = cpy[1] * (g[37] + b00[1]) + b01[1] * c0y[1];
    g[46] = cpy[2] * (g[38] + b00[2]) + b01[2] * c0y[2];
    g[47] = cpy[3] * (g[39] + b00[3]) + b01[3] * c0y[3];
    //g[48] = w[0];
    //g[49] = w[0];
    //g[50] = w[1];
    //g[51] = w[1];
    g[52] = c0z[0] * g[48];
    g[53] = c0z[1] * g[49];
    g[54] = c0z[2] * g[50];
    g[55] = c0z[3] * g[51];
    g[56] = cpz[0] * g[48];
    g[57] = cpz[1] * g[49];
    g[58] = cpz[2] * g[50];
    g[59] = cpz[3] * g[51];
    g[60] = cpz[0] * g[52] + b00[0] * g[48];
    g[61] = cpz[1] * g[53] + b00[1] * g[49];
    g[62] = cpz[2] * g[54] + b00[2] * g[50];
    g[63] = cpz[3] * g[55] + b00[3] * g[51];
    g[64] = cpz[0] * g[56] + b01[0] * g[48];
    g[65] = cpz[1] * g[57] + b01[1] * g[49];
    g[66] = cpz[2] * g[58] + b01[2] * g[50];
    g[67] = cpz[3] * g[59] + b01[3] * g[51];
    g[68] = cpz[0] * g[60] + b01[0] * g[52] + b00[0] * g[56];
    g[69] = cpz[1] * g[61] + b01[1] * g[53] + b00[1] * g[57];
    g[70] = cpz[2] * g[62] + b01[2] * g[54] + b00[2] * g[58];
    g[71] = cpz[3] * g[63] + b01[3] * g[55] + b00[3] * g[59];
}
static inline void _srg0_2d4d_1010(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = c0x[2];
    g[7] = c0x[3];
    g[8] = cpx[0];
    g[9] = cpx[1];
    g[10] = cpx[2];
    g[11] = cpx[3];
    g[12] = cpx[0] * c0x[0] + b00[0];
    g[13] = cpx[1] * c0x[1] + b00[1];
    g[14] = cpx[2] * c0x[2] + b00[2];
    g[15] = cpx[3] * c0x[3] + b00[3];
    g[16] = 1;
    g[17] = 1;
    g[18] = 1;
    g[19] = 1;
    g[20] = c0y[0];
    g[21] = c0y[1];
    g[22] = c0y[2];
    g[23] = c0y[3];
    g[24] = cpy[0];
    g[25] = cpy[1];
    g[26] = cpy[2];
    g[27] = cpy[3];
    g[28] = cpy[0] * c0y[0] + b00[0];
    g[29] = cpy[1] * c0y[1] + b00[1];
    g[30] = cpy[2] * c0y[2] + b00[2];
    g[31] = cpy[3] * c0y[3] + b00[3];
    //g[32] = w[0];
    //g[33] = w[0];
    //g[34] = w[1];
    //g[35] = w[1];
    g[36] = c0z[0] * g[32];
    g[37] = c0z[1] * g[33];
    g[38] = c0z[2] * g[34];
    g[39] = c0z[3] * g[35];
    g[40] = cpz[0] * g[32];
    g[41] = cpz[1] * g[33];
    g[42] = cpz[2] * g[34];
    g[43] = cpz[3] * g[35];
    g[44] = cpz[0] * g[36] + b00[0] * g[32];
    g[45] = cpz[1] * g[37] + b00[1] * g[33];
    g[46] = cpz[2] * g[38] + b00[2] * g[34];
    g[47] = cpz[3] * g[39] + b00[3] * g[35];
}
static inline void _srg0_2d4d_1011(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    double* b01 = bc->b01;
    double xkxl = envs->rkrl[0];
    double ykyl = envs->rkrl[1];
    double zkzl = envs->rkrl[2];
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = c0x[2];
    g[7] = c0x[3];
    g[16] = cpx[0];
    g[17] = cpx[1];
    g[18] = cpx[2];
    g[19] = cpx[3];
    g[20] = cpx[0] * c0x[0] + b00[0];
    g[21] = cpx[1] * c0x[1] + b00[1];
    g[22] = cpx[2] * c0x[2] + b00[2];
    g[23] = cpx[3] * c0x[3] + b00[3];
    g[24] = cpx[0] * (xkxl + cpx[0]) + b01[0];
    g[25] = cpx[1] * (xkxl + cpx[1]) + b01[1];
    g[26] = cpx[2] * (xkxl + cpx[2]) + b01[2];
    g[27] = cpx[3] * (xkxl + cpx[3]) + b01[3];
    g[28] = g[20] * (xkxl + cpx[0]) + cpx[0] * b00[0] + b01[0] * c0x[0];
    g[29] = g[21] * (xkxl + cpx[1]) + cpx[1] * b00[1] + b01[1] * c0x[1];
    g[30] = g[22] * (xkxl + cpx[2]) + cpx[2] * b00[2] + b01[2] * c0x[2];
    g[31] = g[23] * (xkxl + cpx[3]) + cpx[3] * b00[3] + b01[3] * c0x[3];
    g[8] = xkxl + cpx[0];
    g[9] = xkxl + cpx[1];
    g[10] = xkxl + cpx[2];
    g[11] = xkxl + cpx[3];
    g[12] = c0x[0] * (xkxl + cpx[0]) + b00[0];
    g[13] = c0x[1] * (xkxl + cpx[1]) + b00[1];
    g[14] = c0x[2] * (xkxl + cpx[2]) + b00[2];
    g[15] = c0x[3] * (xkxl + cpx[3]) + b00[3];
    g[48] = 1;
    g[49] = 1;
    g[50] = 1;
    g[51] = 1;
    g[52] = c0y[0];
    g[53] = c0y[1];
    g[54] = c0y[2];
    g[55] = c0y[3];
    g[64] = cpy[0];
    g[65] = cpy[1];
    g[66] = cpy[2];
    g[67] = cpy[3];
    g[68] = cpy[0] * c0y[0] + b00[0];
    g[69] = cpy[1] * c0y[1] + b00[1];
    g[70] = cpy[2] * c0y[2] + b00[2];
    g[71] = cpy[3] * c0y[3] + b00[3];
    g[72] = cpy[0] * (ykyl + cpy[0]) + b01[0];
    g[73] = cpy[1] * (ykyl + cpy[1]) + b01[1];
    g[74] = cpy[2] * (ykyl + cpy[2]) + b01[2];
    g[75] = cpy[3] * (ykyl + cpy[3]) + b01[3];
    g[76] = g[68] * (ykyl + cpy[0]) + cpy[0] * b00[0] + b01[0] * c0y[0];
    g[77] = g[69] * (ykyl + cpy[1]) + cpy[1] * b00[1] + b01[1] * c0y[1];
    g[78] = g[70] * (ykyl + cpy[2]) + cpy[2] * b00[2] + b01[2] * c0y[2];
    g[79] = g[71] * (ykyl + cpy[3]) + cpy[3] * b00[3] + b01[3] * c0y[3];
    g[56] = ykyl + cpy[0];
    g[57] = ykyl + cpy[1];
    g[58] = ykyl + cpy[2];
    g[59] = ykyl + cpy[3];
    g[60] = c0y[0] * (ykyl + cpy[0]) + b00[0];
    g[61] = c0y[1] * (ykyl + cpy[1]) + b00[1];
    g[62] = c0y[2] * (ykyl + cpy[2]) + b00[2];
    g[63] = c0y[3] * (ykyl + cpy[3]) + b00[3];
    //g[96] = w[0];
    //g[97] = w[0];
    //g[98] = w[1];
    //g[99] = w[1];
    g[100] = c0z[0] * g[96];
    g[101] = c0z[1] * g[97];
    g[102] = c0z[2] * g[98];
    g[103] = c0z[3] * g[99];
    g[112] = cpz[0] * g[96];
    g[113] = cpz[1] * g[97];
    g[114] = cpz[2] * g[98];
    g[115] = cpz[3] * g[99];
    g[116] = cpz[0] * g[100] + b00[0] * g[96];
    g[117] = cpz[1] * g[101] + b00[1] * g[97];
    g[118] = cpz[2] * g[102] + b00[2] * g[98];
    g[119] = cpz[3] * g[103] + b00[3] * g[99];
    g[120] = g[112] * (zkzl + cpz[0]) + b01[0] * g[96];
    g[121] = g[113] * (zkzl + cpz[1]) + b01[1] * g[97];
    g[122] = g[114] * (zkzl + cpz[2]) + b01[2] * g[98];
    g[123] = g[115] * (zkzl + cpz[3]) + b01[3] * g[99];
    g[124] = g[116] * (zkzl + cpz[0]) + b01[0] * g[100] + b00[0] * g[112];
    g[125] = g[117] * (zkzl + cpz[1]) + b01[1] * g[101] + b00[1] * g[113];
    g[126] = g[118] * (zkzl + cpz[2]) + b01[2] * g[102] + b00[2] * g[114];
    g[127] = g[119] * (zkzl + cpz[3]) + b01[3] * g[103] + b00[3] * g[115];
    g[104] = g[96] * (zkzl + cpz[0]);
    g[105] = g[97] * (zkzl + cpz[1]);
    g[106] = g[98] * (zkzl + cpz[2]);
    g[107] = g[99] * (zkzl + cpz[3]);
    g[108] = g[100] * (zkzl + cpz[0]) + b00[0] * g[96];
    g[109] = g[101] * (zkzl + cpz[1]) + b00[1] * g[97];
    g[110] = g[102] * (zkzl + cpz[2]) + b00[2] * g[98];
    g[111] = g[103] * (zkzl + cpz[3]) + b00[3] * g[99];
}
static inline void _srg0_2d4d_1020(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    double* b01 = bc->b01;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = c0x[2];
    g[7] = c0x[3];
    g[8] = cpx[0];
    g[9] = cpx[1];
    g[10] = cpx[2];
    g[11] = cpx[3];
    g[12] = cpx[0] * c0x[0] + b00[0];
    g[13] = cpx[1] * c0x[1] + b00[1];
    g[14] = cpx[2] * c0x[2] + b00[2];
    g[15] = cpx[3] * c0x[3] + b00[3];
    g[16] = cpx[0] * cpx[0] + b01[0];
    g[17] = cpx[1] * cpx[1] + b01[1];
    g[18] = cpx[2] * cpx[2] + b01[2];
    g[19] = cpx[3] * cpx[3] + b01[3];
    g[20] = cpx[0] * (g[12] + b00[0]) + b01[0] * c0x[0];
    g[21] = cpx[1] * (g[13] + b00[1]) + b01[1] * c0x[1];
    g[22] = cpx[2] * (g[14] + b00[2]) + b01[2] * c0x[2];
    g[23] = cpx[3] * (g[15] + b00[3]) + b01[3] * c0x[3];
    g[24] = 1;
    g[25] = 1;
    g[26] = 1;
    g[27] = 1;
    g[28] = c0y[0];
    g[29] = c0y[1];
    g[30] = c0y[2];
    g[31] = c0y[3];
    g[32] = cpy[0];
    g[33] = cpy[1];
    g[34] = cpy[2];
    g[35] = cpy[3];
    g[36] = cpy[0] * c0y[0] + b00[0];
    g[37] = cpy[1] * c0y[1] + b00[1];
    g[38] = cpy[2] * c0y[2] + b00[2];
    g[39] = cpy[3] * c0y[3] + b00[3];
    g[40] = cpy[0] * cpy[0] + b01[0];
    g[41] = cpy[1] * cpy[1] + b01[1];
    g[42] = cpy[2] * cpy[2] + b01[2];
    g[43] = cpy[3] * cpy[3] + b01[3];
    g[44] = cpy[0] * (g[36] + b00[0]) + b01[0] * c0y[0];
    g[45] = cpy[1] * (g[37] + b00[1]) + b01[1] * c0y[1];
    g[46] = cpy[2] * (g[38] + b00[2]) + b01[2] * c0y[2];
    g[47] = cpy[3] * (g[39] + b00[3]) + b01[3] * c0y[3];
    //g[48] = w[0];
    //g[49] = w[0];
    //g[50] = w[1];
    //g[51] = w[1];
    g[52] = c0z[0] * g[48];
    g[53] = c0z[1] * g[49];
    g[54] = c0z[2] * g[50];
    g[55] = c0z[3] * g[51];
    g[56] = cpz[0] * g[48];
    g[57] = cpz[1] * g[49];
    g[58] = cpz[2] * g[50];
    g[59] = cpz[3] * g[51];
    g[60] = cpz[0] * g[52] + b00[0] * g[48];
    g[61] = cpz[1] * g[53] + b00[1] * g[49];
    g[62] = cpz[2] * g[54] + b00[2] * g[50];
    g[63] = cpz[3] * g[55] + b00[3] * g[51];
    g[64] = cpz[0] * g[56] + b01[0] * g[48];
    g[65] = cpz[1] * g[57] + b01[1] * g[49];
    g[66] = cpz[2] * g[58] + b01[2] * g[50];
    g[67] = cpz[3] * g[59] + b01[3] * g[51];
    g[68] = cpz[0] * g[60] + b01[0] * g[52] + b00[0] * g[56];
    g[69] = cpz[1] * g[61] + b01[1] * g[53] + b00[1] * g[57];
    g[70] = cpz[2] * g[62] + b01[2] * g[54] + b00[2] * g[58];
    g[71] = cpz[3] * g[63] + b01[3] * g[55] + b00[3] * g[59];
}
static inline void _srg0_2d4d_1100(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* b10 = bc->b10;
    double xixj = envs->rirj[0];
    double yiyj = envs->rirj[1];
    double zizj = envs->rirj[2];
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[8] = c0x[0];
    g[9] = c0x[1];
    g[10] = c0x[2];
    g[11] = c0x[3];
    g[12] = c0x[0] * (xixj + c0x[0]) + b10[0];
    g[13] = c0x[1] * (xixj + c0x[1]) + b10[1];
    g[14] = c0x[2] * (xixj + c0x[2]) + b10[2];
    g[15] = c0x[3] * (xixj + c0x[3]) + b10[3];
    g[4] = xixj + c0x[0];
    g[5] = xixj + c0x[1];
    g[6] = xixj + c0x[2];
    g[7] = xixj + c0x[3];
    g[24] = 1;
    g[25] = 1;
    g[26] = 1;
    g[27] = 1;
    g[32] = c0y[0];
    g[33] = c0y[1];
    g[34] = c0y[2];
    g[35] = c0y[3];
    g[36] = c0y[0] * (yiyj + c0y[0]) + b10[0];
    g[37] = c0y[1] * (yiyj + c0y[1]) + b10[1];
    g[38] = c0y[2] * (yiyj + c0y[2]) + b10[2];
    g[39] = c0y[3] * (yiyj + c0y[3]) + b10[3];
    g[28] = yiyj + c0y[0];
    g[29] = yiyj + c0y[1];
    g[30] = yiyj + c0y[2];
    g[31] = yiyj + c0y[3];
    //g[48] = w[0];
    //g[49] = w[0];
    //g[50] = w[1];
    //g[51] = w[1];
    g[56] = c0z[0] * g[48];
    g[57] = c0z[1] * g[49];
    g[58] = c0z[2] * g[50];
    g[59] = c0z[3] * g[51];
    g[60] = g[56] * (zizj + c0z[0]) + b10[0] * g[48];
    g[61] = g[57] * (zizj + c0z[1]) + b10[1] * g[49];
    g[62] = g[58] * (zizj + c0z[2]) + b10[2] * g[50];
    g[63] = g[59] * (zizj + c0z[3]) + b10[3] * g[51];
    g[52] = g[48] * (zizj + c0z[0]);
    g[53] = g[49] * (zizj + c0z[1]);
    g[54] = g[50] * (zizj + c0z[2]);
    g[55] = g[51] * (zizj + c0z[3]);
}
static inline void _srg0_2d4d_1101(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    double* b10 = bc->b10;
    double xixj = envs->rirj[0];
    double yiyj = envs->rirj[1];
    double zizj = envs->rirj[2];
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[16] = c0x[0];
    g[17] = c0x[1];
    g[18] = c0x[2];
    g[19] = c0x[3];
    g[8] = cpx[0];
    g[9] = cpx[1];
    g[10] = cpx[2];
    g[11] = cpx[3];
    g[24] = cpx[0] * c0x[0] + b00[0];
    g[25] = cpx[1] * c0x[1] + b00[1];
    g[26] = cpx[2] * c0x[2] + b00[2];
    g[27] = cpx[3] * c0x[3] + b00[3];
    g[20] = c0x[0] * (xixj + c0x[0]) + b10[0];
    g[21] = c0x[1] * (xixj + c0x[1]) + b10[1];
    g[22] = c0x[2] * (xixj + c0x[2]) + b10[2];
    g[23] = c0x[3] * (xixj + c0x[3]) + b10[3];
    g[4] = xixj + c0x[0];
    g[5] = xixj + c0x[1];
    g[6] = xixj + c0x[2];
    g[7] = xixj + c0x[3];
    g[28] = g[24] * (xixj + c0x[0]) + c0x[0] * b00[0] + b10[0] * cpx[0];
    g[29] = g[25] * (xixj + c0x[1]) + c0x[1] * b00[1] + b10[1] * cpx[1];
    g[30] = g[26] * (xixj + c0x[2]) + c0x[2] * b00[2] + b10[2] * cpx[2];
    g[31] = g[27] * (xixj + c0x[3]) + c0x[3] * b00[3] + b10[3] * cpx[3];
    g[12] = cpx[0] * (xixj + c0x[0]) + b00[0];
    g[13] = cpx[1] * (xixj + c0x[1]) + b00[1];
    g[14] = cpx[2] * (xixj + c0x[2]) + b00[2];
    g[15] = cpx[3] * (xixj + c0x[3]) + b00[3];
    g[48] = 1;
    g[49] = 1;
    g[50] = 1;
    g[51] = 1;
    g[64] = c0y[0];
    g[65] = c0y[1];
    g[66] = c0y[2];
    g[67] = c0y[3];
    g[56] = cpy[0];
    g[57] = cpy[1];
    g[58] = cpy[2];
    g[59] = cpy[3];
    g[72] = cpy[0] * c0y[0] + b00[0];
    g[73] = cpy[1] * c0y[1] + b00[1];
    g[74] = cpy[2] * c0y[2] + b00[2];
    g[75] = cpy[3] * c0y[3] + b00[3];
    g[68] = c0y[0] * (yiyj + c0y[0]) + b10[0];
    g[69] = c0y[1] * (yiyj + c0y[1]) + b10[1];
    g[70] = c0y[2] * (yiyj + c0y[2]) + b10[2];
    g[71] = c0y[3] * (yiyj + c0y[3]) + b10[3];
    g[52] = yiyj + c0y[0];
    g[53] = yiyj + c0y[1];
    g[54] = yiyj + c0y[2];
    g[55] = yiyj + c0y[3];
    g[76] = g[72] * (yiyj + c0y[0]) + c0y[0] * b00[0] + b10[0] * cpy[0];
    g[77] = g[73] * (yiyj + c0y[1]) + c0y[1] * b00[1] + b10[1] * cpy[1];
    g[78] = g[74] * (yiyj + c0y[2]) + c0y[2] * b00[2] + b10[2] * cpy[2];
    g[79] = g[75] * (yiyj + c0y[3]) + c0y[3] * b00[3] + b10[3] * cpy[3];
    g[60] = cpy[0] * (yiyj + c0y[0]) + b00[0];
    g[61] = cpy[1] * (yiyj + c0y[1]) + b00[1];
    g[62] = cpy[2] * (yiyj + c0y[2]) + b00[2];
    g[63] = cpy[3] * (yiyj + c0y[3]) + b00[3];
    //g[96] = w[0];
    //g[97] = w[0];
    //g[98] = w[1];
    //g[99] = w[1];
    g[112] = c0z[0] * g[96];
    g[113] = c0z[1] * g[97];
    g[114] = c0z[2] * g[98];
    g[115] = c0z[3] * g[99];
    g[104] = cpz[0] * g[96];
    g[105] = cpz[1] * g[97];
    g[106] = cpz[2] * g[98];
    g[107] = cpz[3] * g[99];
    g[120] = cpz[0] * g[112] + b00[0] * g[96];
    g[121] = cpz[1] * g[113] + b00[1] * g[97];
    g[122] = cpz[2] * g[114] + b00[2] * g[98];
    g[123] = cpz[3] * g[115] + b00[3] * g[99];
    g[116] = g[112] * (zizj + c0z[0]) + b10[0] * g[96];
    g[117] = g[113] * (zizj + c0z[1]) + b10[1] * g[97];
    g[118] = g[114] * (zizj + c0z[2]) + b10[2] * g[98];
    g[119] = g[115] * (zizj + c0z[3]) + b10[3] * g[99];
    g[100] = g[96] * (zizj + c0z[0]);
    g[101] = g[97] * (zizj + c0z[1]);
    g[102] = g[98] * (zizj + c0z[2]);
    g[103] = g[99] * (zizj + c0z[3]);
    g[124] = g[120] * (zizj + c0z[0]) + b10[0] * g[104] + b00[0] * g[112];
    g[125] = g[121] * (zizj + c0z[1]) + b10[1] * g[105] + b00[1] * g[113];
    g[126] = g[122] * (zizj + c0z[2]) + b10[2] * g[106] + b00[2] * g[114];
    g[127] = g[123] * (zizj + c0z[3]) + b10[3] * g[107] + b00[3] * g[115];
    g[108] = zizj * g[104] + cpz[0] * g[112] + b00[0] * g[96];
    g[109] = zizj * g[105] + cpz[1] * g[113] + b00[1] * g[97];
    g[110] = zizj * g[106] + cpz[2] * g[114] + b00[2] * g[98];
    g[111] = zizj * g[107] + cpz[3] * g[115] + b00[3] * g[99];
}
static inline void _srg0_2d4d_1110(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    double* b10 = bc->b10;
    double xixj = envs->rirj[0];
    double yiyj = envs->rirj[1];
    double zizj = envs->rirj[2];
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[16] = c0x[0];
    g[17] = c0x[1];
    g[18] = c0x[2];
    g[19] = c0x[3];
    g[8] = cpx[0];
    g[9] = cpx[1];
    g[10] = cpx[2];
    g[11] = cpx[3];
    g[24] = cpx[0] * c0x[0] + b00[0];
    g[25] = cpx[1] * c0x[1] + b00[1];
    g[26] = cpx[2] * c0x[2] + b00[2];
    g[27] = cpx[3] * c0x[3] + b00[3];
    g[20] = c0x[0] * (xixj + c0x[0]) + b10[0];
    g[21] = c0x[1] * (xixj + c0x[1]) + b10[1];
    g[22] = c0x[2] * (xixj + c0x[2]) + b10[2];
    g[23] = c0x[3] * (xixj + c0x[3]) + b10[3];
    g[4] = xixj + c0x[0];
    g[5] = xixj + c0x[1];
    g[6] = xixj + c0x[2];
    g[7] = xixj + c0x[3];
    g[28] = g[24] * (xixj + c0x[0]) + c0x[0] * b00[0] + b10[0] * cpx[0];
    g[29] = g[25] * (xixj + c0x[1]) + c0x[1] * b00[1] + b10[1] * cpx[1];
    g[30] = g[26] * (xixj + c0x[2]) + c0x[2] * b00[2] + b10[2] * cpx[2];
    g[31] = g[27] * (xixj + c0x[3]) + c0x[3] * b00[3] + b10[3] * cpx[3];
    g[12] = cpx[0] * (xixj + c0x[0]) + b00[0];
    g[13] = cpx[1] * (xixj + c0x[1]) + b00[1];
    g[14] = cpx[2] * (xixj + c0x[2]) + b00[2];
    g[15] = cpx[3] * (xixj + c0x[3]) + b00[3];
    g[48] = 1;
    g[49] = 1;
    g[50] = 1;
    g[51] = 1;
    g[64] = c0y[0];
    g[65] = c0y[1];
    g[66] = c0y[2];
    g[67] = c0y[3];
    g[56] = cpy[0];
    g[57] = cpy[1];
    g[58] = cpy[2];
    g[59] = cpy[3];
    g[72] = cpy[0] * c0y[0] + b00[0];
    g[73] = cpy[1] * c0y[1] + b00[1];
    g[74] = cpy[2] * c0y[2] + b00[2];
    g[75] = cpy[3] * c0y[3] + b00[3];
    g[68] = c0y[0] * (yiyj + c0y[0]) + b10[0];
    g[69] = c0y[1] * (yiyj + c0y[1]) + b10[1];
    g[70] = c0y[2] * (yiyj + c0y[2]) + b10[2];
    g[71] = c0y[3] * (yiyj + c0y[3]) + b10[3];
    g[52] = yiyj + c0y[0];
    g[53] = yiyj + c0y[1];
    g[54] = yiyj + c0y[2];
    g[55] = yiyj + c0y[3];
    g[76] = g[72] * (yiyj + c0y[0]) + c0y[0] * b00[0] + b10[0] * cpy[0];
    g[77] = g[73] * (yiyj + c0y[1]) + c0y[1] * b00[1] + b10[1] * cpy[1];
    g[78] = g[74] * (yiyj + c0y[2]) + c0y[2] * b00[2] + b10[2] * cpy[2];
    g[79] = g[75] * (yiyj + c0y[3]) + c0y[3] * b00[3] + b10[3] * cpy[3];
    g[60] = cpy[0] * (yiyj + c0y[0]) + b00[0];
    g[61] = cpy[1] * (yiyj + c0y[1]) + b00[1];
    g[62] = cpy[2] * (yiyj + c0y[2]) + b00[2];
    g[63] = cpy[3] * (yiyj + c0y[3]) + b00[3];
    //g[96] = w[0];
    //g[97] = w[0];
    //g[98] = w[1];
    //g[99] = w[1];
    g[112] = c0z[0] * g[96];
    g[113] = c0z[1] * g[97];
    g[114] = c0z[2] * g[98];
    g[115] = c0z[3] * g[99];
    g[104] = cpz[0] * g[96];
    g[105] = cpz[1] * g[97];
    g[106] = cpz[2] * g[98];
    g[107] = cpz[3] * g[99];
    g[120] = cpz[0] * g[112] + b00[0] * g[96];
    g[121] = cpz[1] * g[113] + b00[1] * g[97];
    g[122] = cpz[2] * g[114] + b00[2] * g[98];
    g[123] = cpz[3] * g[115] + b00[3] * g[99];
    g[116] = g[112] * (zizj + c0z[0]) + b10[0] * g[96];
    g[117] = g[113] * (zizj + c0z[1]) + b10[1] * g[97];
    g[118] = g[114] * (zizj + c0z[2]) + b10[2] * g[98];
    g[119] = g[115] * (zizj + c0z[3]) + b10[3] * g[99];
    g[100] = g[96] * (zizj + c0z[0]);
    g[101] = g[97] * (zizj + c0z[1]);
    g[102] = g[98] * (zizj + c0z[2]);
    g[103] = g[99] * (zizj + c0z[3]);
    g[124] = g[120] * (zizj + c0z[0]) + b10[0] * g[104] + b00[0] * g[112];
    g[125] = g[121] * (zizj + c0z[1]) + b10[1] * g[105] + b00[1] * g[113];
    g[126] = g[122] * (zizj + c0z[2]) + b10[2] * g[106] + b00[2] * g[114];
    g[127] = g[123] * (zizj + c0z[3]) + b10[3] * g[107] + b00[3] * g[115];
    g[108] = zizj * g[104] + cpz[0] * g[112] + b00[0] * g[96];
    g[109] = zizj * g[105] + cpz[1] * g[113] + b00[1] * g[97];
    g[110] = zizj * g[106] + cpz[2] * g[114] + b00[2] * g[98];
    g[111] = zizj * g[107] + cpz[3] * g[115] + b00[3] * g[99];
}
static inline void _srg0_2d4d_1200(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* b10 = bc->b10;
    double xixj = envs->rirj[0];
    double yiyj = envs->rirj[1];
    double zizj = envs->rirj[2];
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[8] = c0x[0];
    g[9] = c0x[1];
    g[10] = c0x[2];
    g[11] = c0x[3];
    g[16] = c0x[0] * c0x[0] + b10[0];
    g[17] = c0x[1] * c0x[1] + b10[1];
    g[18] = c0x[2] * c0x[2] + b10[2];
    g[19] = c0x[3] * c0x[3] + b10[3];
    g[20] = g[16] * (xixj + c0x[0]) + c0x[0] * 2 * b10[0];
    g[21] = g[17] * (xixj + c0x[1]) + c0x[1] * 2 * b10[1];
    g[22] = g[18] * (xixj + c0x[2]) + c0x[2] * 2 * b10[2];
    g[23] = g[19] * (xixj + c0x[3]) + c0x[3] * 2 * b10[3];
    g[12] = c0x[0] * (xixj + c0x[0]) + b10[0];
    g[13] = c0x[1] * (xixj + c0x[1]) + b10[1];
    g[14] = c0x[2] * (xixj + c0x[2]) + b10[2];
    g[15] = c0x[3] * (xixj + c0x[3]) + b10[3];
    g[4] = xixj + c0x[0];
    g[5] = xixj + c0x[1];
    g[6] = xixj + c0x[2];
    g[7] = xixj + c0x[3];
    g[32] = 1;
    g[33] = 1;
    g[34] = 1;
    g[35] = 1;
    g[40] = c0y[0];
    g[41] = c0y[1];
    g[42] = c0y[2];
    g[43] = c0y[3];
    g[48] = c0y[0] * c0y[0] + b10[0];
    g[49] = c0y[1] * c0y[1] + b10[1];
    g[50] = c0y[2] * c0y[2] + b10[2];
    g[51] = c0y[3] * c0y[3] + b10[3];
    g[52] = g[48] * (yiyj + c0y[0]) + c0y[0] * 2 * b10[0];
    g[53] = g[49] * (yiyj + c0y[1]) + c0y[1] * 2 * b10[1];
    g[54] = g[50] * (yiyj + c0y[2]) + c0y[2] * 2 * b10[2];
    g[55] = g[51] * (yiyj + c0y[3]) + c0y[3] * 2 * b10[3];
    g[44] = c0y[0] * (yiyj + c0y[0]) + b10[0];
    g[45] = c0y[1] * (yiyj + c0y[1]) + b10[1];
    g[46] = c0y[2] * (yiyj + c0y[2]) + b10[2];
    g[47] = c0y[3] * (yiyj + c0y[3]) + b10[3];
    g[36] = yiyj + c0y[0];
    g[37] = yiyj + c0y[1];
    g[38] = yiyj + c0y[2];
    g[39] = yiyj + c0y[3];
    //g[64] = w[0];
    //g[65] = w[0];
    //g[66] = w[1];
    //g[67] = w[1];
    g[72] = c0z[0] * g[64];
    g[73] = c0z[1] * g[65];
    g[74] = c0z[2] * g[66];
    g[75] = c0z[3] * g[67];
    g[80] = c0z[0] * g[72] + b10[0] * g[64];
    g[81] = c0z[1] * g[73] + b10[1] * g[65];
    g[82] = c0z[2] * g[74] + b10[2] * g[66];
    g[83] = c0z[3] * g[75] + b10[3] * g[67];
    g[84] = g[80] * (zizj + c0z[0]) + 2 * b10[0] * g[72];
    g[85] = g[81] * (zizj + c0z[1]) + 2 * b10[1] * g[73];
    g[86] = g[82] * (zizj + c0z[2]) + 2 * b10[2] * g[74];
    g[87] = g[83] * (zizj + c0z[3]) + 2 * b10[3] * g[75];
    g[76] = g[72] * (zizj + c0z[0]) + b10[0] * g[64];
    g[77] = g[73] * (zizj + c0z[1]) + b10[1] * g[65];
    g[78] = g[74] * (zizj + c0z[2]) + b10[2] * g[66];
    g[79] = g[75] * (zizj + c0z[3]) + b10[3] * g[67];
    g[68] = g[64] * (zizj + c0z[0]);
    g[69] = g[65] * (zizj + c0z[1]);
    g[70] = g[66] * (zizj + c0z[2]);
    g[71] = g[67] * (zizj + c0z[3]);
}
static inline void _srg0_2d4d_2000(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = c0x[2];
    g[7] = c0x[3];
    g[8] = c0x[0] * c0x[0] + b10[0];
    g[9] = c0x[1] * c0x[1] + b10[1];
    g[10] = c0x[2] * c0x[2] + b10[2];
    g[11] = c0x[3] * c0x[3] + b10[3];
    g[12] = 1;
    g[13] = 1;
    g[14] = 1;
    g[15] = 1;
    g[16] = c0y[0];
    g[17] = c0y[1];
    g[18] = c0y[2];
    g[19] = c0y[3];
    g[20] = c0y[0] * c0y[0] + b10[0];
    g[21] = c0y[1] * c0y[1] + b10[1];
    g[22] = c0y[2] * c0y[2] + b10[2];
    g[23] = c0y[3] * c0y[3] + b10[3];
    //g[24] = w[0];
    //g[25] = w[0];
    //g[26] = w[1];
    //g[27] = w[1];
    g[28] = c0z[0] * g[24];
    g[29] = c0z[1] * g[25];
    g[30] = c0z[2] * g[26];
    g[31] = c0z[3] * g[27];
    g[32] = c0z[0] * g[28] + b10[0] * g[24];
    g[33] = c0z[1] * g[29] + b10[1] * g[25];
    g[34] = c0z[2] * g[30] + b10[2] * g[26];
    g[35] = c0z[3] * g[31] + b10[3] * g[27];
}
static inline void _srg0_2d4d_2001(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    double* b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = c0x[2];
    g[7] = c0x[3];
    g[8] = c0x[0] * c0x[0] + b10[0];
    g[9] = c0x[1] * c0x[1] + b10[1];
    g[10] = c0x[2] * c0x[2] + b10[2];
    g[11] = c0x[3] * c0x[3] + b10[3];
    g[12] = cpx[0];
    g[13] = cpx[1];
    g[14] = cpx[2];
    g[15] = cpx[3];
    g[16] = cpx[0] * c0x[0] + b00[0];
    g[17] = cpx[1] * c0x[1] + b00[1];
    g[18] = cpx[2] * c0x[2] + b00[2];
    g[19] = cpx[3] * c0x[3] + b00[3];
    g[20] = c0x[0] * (g[16] + b00[0]) + b10[0] * cpx[0];
    g[21] = c0x[1] * (g[17] + b00[1]) + b10[1] * cpx[1];
    g[22] = c0x[2] * (g[18] + b00[2]) + b10[2] * cpx[2];
    g[23] = c0x[3] * (g[19] + b00[3]) + b10[3] * cpx[3];
    g[24] = 1;
    g[25] = 1;
    g[26] = 1;
    g[27] = 1;
    g[28] = c0y[0];
    g[29] = c0y[1];
    g[30] = c0y[2];
    g[31] = c0y[3];
    g[32] = c0y[0] * c0y[0] + b10[0];
    g[33] = c0y[1] * c0y[1] + b10[1];
    g[34] = c0y[2] * c0y[2] + b10[2];
    g[35] = c0y[3] * c0y[3] + b10[3];
    g[36] = cpy[0];
    g[37] = cpy[1];
    g[38] = cpy[2];
    g[39] = cpy[3];
    g[40] = cpy[0] * c0y[0] + b00[0];
    g[41] = cpy[1] * c0y[1] + b00[1];
    g[42] = cpy[2] * c0y[2] + b00[2];
    g[43] = cpy[3] * c0y[3] + b00[3];
    g[44] = c0y[0] * (g[40] + b00[0]) + b10[0] * cpy[0];
    g[45] = c0y[1] * (g[41] + b00[1]) + b10[1] * cpy[1];
    g[46] = c0y[2] * (g[42] + b00[2]) + b10[2] * cpy[2];
    g[47] = c0y[3] * (g[43] + b00[3]) + b10[3] * cpy[3];
    //g[48] = w[0];
    //g[49] = w[0];
    //g[50] = w[1];
    //g[51] = w[1];
    g[52] = c0z[0] * g[48];
    g[53] = c0z[1] * g[49];
    g[54] = c0z[2] * g[50];
    g[55] = c0z[3] * g[51];
    g[56] = c0z[0] * g[52] + b10[0] * g[48];
    g[57] = c0z[1] * g[53] + b10[1] * g[49];
    g[58] = c0z[2] * g[54] + b10[2] * g[50];
    g[59] = c0z[3] * g[55] + b10[3] * g[51];
    g[60] = cpz[0] * g[48];
    g[61] = cpz[1] * g[49];
    g[62] = cpz[2] * g[50];
    g[63] = cpz[3] * g[51];
    g[64] = cpz[0] * g[52] + b00[0] * g[48];
    g[65] = cpz[1] * g[53] + b00[1] * g[49];
    g[66] = cpz[2] * g[54] + b00[2] * g[50];
    g[67] = cpz[3] * g[55] + b00[3] * g[51];
    g[68] = c0z[0] * g[64] + b10[0] * g[60] + b00[0] * g[52];
    g[69] = c0z[1] * g[65] + b10[1] * g[61] + b00[1] * g[53];
    g[70] = c0z[2] * g[66] + b10[2] * g[62] + b00[2] * g[54];
    g[71] = c0z[3] * g[67] + b10[3] * g[63] + b00[3] * g[55];
}
static inline void _srg0_2d4d_2010(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* cpx = bc->c0px;
    double* cpy = bc->c0py;
    double* cpz = bc->c0pz;
    double* b00 = bc->b00;
    double* b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = c0x[2];
    g[7] = c0x[3];
    g[8] = c0x[0] * c0x[0] + b10[0];
    g[9] = c0x[1] * c0x[1] + b10[1];
    g[10] = c0x[2] * c0x[2] + b10[2];
    g[11] = c0x[3] * c0x[3] + b10[3];
    g[12] = cpx[0];
    g[13] = cpx[1];
    g[14] = cpx[2];
    g[15] = cpx[3];
    g[16] = cpx[0] * c0x[0] + b00[0];
    g[17] = cpx[1] * c0x[1] + b00[1];
    g[18] = cpx[2] * c0x[2] + b00[2];
    g[19] = cpx[3] * c0x[3] + b00[3];
    g[20] = c0x[0] * (g[16] + b00[0]) + b10[0] * cpx[0];
    g[21] = c0x[1] * (g[17] + b00[1]) + b10[1] * cpx[1];
    g[22] = c0x[2] * (g[18] + b00[2]) + b10[2] * cpx[2];
    g[23] = c0x[3] * (g[19] + b00[3]) + b10[3] * cpx[3];
    g[24] = 1;
    g[25] = 1;
    g[26] = 1;
    g[27] = 1;
    g[28] = c0y[0];
    g[29] = c0y[1];
    g[30] = c0y[2];
    g[31] = c0y[3];
    g[32] = c0y[0] * c0y[0] + b10[0];
    g[33] = c0y[1] * c0y[1] + b10[1];
    g[34] = c0y[2] * c0y[2] + b10[2];
    g[35] = c0y[3] * c0y[3] + b10[3];
    g[36] = cpy[0];
    g[37] = cpy[1];
    g[38] = cpy[2];
    g[39] = cpy[3];
    g[40] = cpy[0] * c0y[0] + b00[0];
    g[41] = cpy[1] * c0y[1] + b00[1];
    g[42] = cpy[2] * c0y[2] + b00[2];
    g[43] = cpy[3] * c0y[3] + b00[3];
    g[44] = c0y[0] * (g[40] + b00[0]) + b10[0] * cpy[0];
    g[45] = c0y[1] * (g[41] + b00[1]) + b10[1] * cpy[1];
    g[46] = c0y[2] * (g[42] + b00[2]) + b10[2] * cpy[2];
    g[47] = c0y[3] * (g[43] + b00[3]) + b10[3] * cpy[3];
    //g[48] = w[0];
    //g[49] = w[0];
    //g[50] = w[1];
    //g[51] = w[1];
    g[52] = c0z[0] * g[48];
    g[53] = c0z[1] * g[49];
    g[54] = c0z[2] * g[50];
    g[55] = c0z[3] * g[51];
    g[56] = c0z[0] * g[52] + b10[0] * g[48];
    g[57] = c0z[1] * g[53] + b10[1] * g[49];
    g[58] = c0z[2] * g[54] + b10[2] * g[50];
    g[59] = c0z[3] * g[55] + b10[3] * g[51];
    g[60] = cpz[0] * g[48];
    g[61] = cpz[1] * g[49];
    g[62] = cpz[2] * g[50];
    g[63] = cpz[3] * g[51];
    g[64] = cpz[0] * g[52] + b00[0] * g[48];
    g[65] = cpz[1] * g[53] + b00[1] * g[49];
    g[66] = cpz[2] * g[54] + b00[2] * g[50];
    g[67] = cpz[3] * g[55] + b00[3] * g[51];
    g[68] = c0z[0] * g[64] + b10[0] * g[60] + b00[0] * g[52];
    g[69] = c0z[1] * g[65] + b10[1] * g[61] + b00[1] * g[53];
    g[70] = c0z[2] * g[66] + b10[2] * g[62] + b00[2] * g[54];
    g[71] = c0z[3] * g[67] + b10[3] * g[63] + b00[3] * g[55];
}
static inline void _srg0_2d4d_2100(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* b10 = bc->b10;
    double xixj = envs->rirj[0];
    double yiyj = envs->rirj[1];
    double zizj = envs->rirj[2];
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = c0x[2];
    g[7] = c0x[3];
    g[8] = c0x[0] * c0x[0] + b10[0];
    g[9] = c0x[1] * c0x[1] + b10[1];
    g[10] = c0x[2] * c0x[2] + b10[2];
    g[11] = c0x[3] * c0x[3] + b10[3];
    g[24] = g[8] * (xixj + c0x[0]) + c0x[0] * 2 * b10[0];
    g[25] = g[9] * (xixj + c0x[1]) + c0x[1] * 2 * b10[1];
    g[26] = g[10] * (xixj + c0x[2]) + c0x[2] * 2 * b10[2];
    g[27] = g[11] * (xixj + c0x[3]) + c0x[3] * 2 * b10[3];
    g[20] = c0x[0] * (xixj + c0x[0]) + b10[0];
    g[21] = c0x[1] * (xixj + c0x[1]) + b10[1];
    g[22] = c0x[2] * (xixj + c0x[2]) + b10[2];
    g[23] = c0x[3] * (xixj + c0x[3]) + b10[3];
    g[16] = xixj + c0x[0];
    g[17] = xixj + c0x[1];
    g[18] = xixj + c0x[2];
    g[19] = xixj + c0x[3];
    g[32] = 1;
    g[33] = 1;
    g[34] = 1;
    g[35] = 1;
    g[36] = c0y[0];
    g[37] = c0y[1];
    g[38] = c0y[2];
    g[39] = c0y[3];
    g[40] = c0y[0] * c0y[0] + b10[0];
    g[41] = c0y[1] * c0y[1] + b10[1];
    g[42] = c0y[2] * c0y[2] + b10[2];
    g[43] = c0y[3] * c0y[3] + b10[3];
    g[56] = g[40] * (yiyj + c0y[0]) + c0y[0] * 2 * b10[0];
    g[57] = g[41] * (yiyj + c0y[1]) + c0y[1] * 2 * b10[1];
    g[58] = g[42] * (yiyj + c0y[2]) + c0y[2] * 2 * b10[2];
    g[59] = g[43] * (yiyj + c0y[3]) + c0y[3] * 2 * b10[3];
    g[52] = c0y[0] * (yiyj + c0y[0]) + b10[0];
    g[53] = c0y[1] * (yiyj + c0y[1]) + b10[1];
    g[54] = c0y[2] * (yiyj + c0y[2]) + b10[2];
    g[55] = c0y[3] * (yiyj + c0y[3]) + b10[3];
    g[48] = yiyj + c0y[0];
    g[49] = yiyj + c0y[1];
    g[50] = yiyj + c0y[2];
    g[51] = yiyj + c0y[3];
    //g[64] = w[0];
    //g[65] = w[0];
    //g[66] = w[1];
    //g[67] = w[1];
    g[68] = c0z[0] * g[64];
    g[69] = c0z[1] * g[65];
    g[70] = c0z[2] * g[66];
    g[71] = c0z[3] * g[67];
    g[72] = c0z[0] * g[68] + b10[0] * g[64];
    g[73] = c0z[1] * g[69] + b10[1] * g[65];
    g[74] = c0z[2] * g[70] + b10[2] * g[66];
    g[75] = c0z[3] * g[71] + b10[3] * g[67];
    g[88] = g[72] * (zizj + c0z[0]) + 2 * b10[0] * g[68];
    g[89] = g[73] * (zizj + c0z[1]) + 2 * b10[1] * g[69];
    g[90] = g[74] * (zizj + c0z[2]) + 2 * b10[2] * g[70];
    g[91] = g[75] * (zizj + c0z[3]) + 2 * b10[3] * g[71];
    g[84] = g[68] * (zizj + c0z[0]) + b10[0] * g[64];
    g[85] = g[69] * (zizj + c0z[1]) + b10[1] * g[65];
    g[86] = g[70] * (zizj + c0z[2]) + b10[2] * g[66];
    g[87] = g[71] * (zizj + c0z[3]) + b10[3] * g[67];
    g[80] = g[64] * (zizj + c0z[0]);
    g[81] = g[65] * (zizj + c0z[1]);
    g[82] = g[66] * (zizj + c0z[2]);
    g[83] = g[67] * (zizj + c0z[3]);
}
static inline void _srg0_2d4d_3000(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    double* c0x = bc->c00x;
    double* c0y = bc->c00y;
    double* c0z = bc->c00z;
    double* b10 = bc->b10;
    g[0] = 1;
    g[1] = 1;
    g[2] = 1;
    g[3] = 1;
    g[4] = c0x[0];
    g[5] = c0x[1];
    g[6] = c0x[2];
    g[7] = c0x[3];
    g[8] = c0x[0] * c0x[0] + b10[0];
    g[9] = c0x[1] * c0x[1] + b10[1];
    g[10] = c0x[2] * c0x[2] + b10[2];
    g[11] = c0x[3] * c0x[3] + b10[3];
    g[12] = c0x[0] * (g[8] + 2 * b10[0]);
    g[13] = c0x[1] * (g[9] + 2 * b10[1]);
    g[14] = c0x[2] * (g[10] + 2 * b10[2]);
    g[15] = c0x[3] * (g[11] + 2 * b10[3]);
    g[16] = 1;
    g[17] = 1;
    g[18] = 1;
    g[19] = 1;
    g[20] = c0y[0];
    g[21] = c0y[1];
    g[22] = c0y[2];
    g[23] = c0y[3];
    g[24] = c0y[0] * c0y[0] + b10[0];
    g[25] = c0y[1] * c0y[1] + b10[1];
    g[26] = c0y[2] * c0y[2] + b10[2];
    g[27] = c0y[3] * c0y[3] + b10[3];
    g[28] = c0y[0] * (g[24] + 2 * b10[0]);
    g[29] = c0y[1] * (g[25] + 2 * b10[1]);
    g[30] = c0y[2] * (g[26] + 2 * b10[2]);
    g[31] = c0y[3] * (g[27] + 2 * b10[3]);
    //g[32] = w[0];
    //g[33] = w[0];
    //g[34] = w[1];
    //g[35] = w[1];
    g[36] = c0z[0] * g[32];
    g[37] = c0z[1] * g[33];
    g[38] = c0z[2] * g[34];
    g[39] = c0z[3] * g[35];
    g[40] = c0z[0] * g[36] + b10[0] * g[32];
    g[41] = c0z[1] * g[37] + b10[1] * g[33];
    g[42] = c0z[2] * g[38] + b10[2] * g[34];
    g[43] = c0z[3] * g[39] + b10[3] * g[35];
    g[44] = c0z[0] * g[40] + 2 * b10[0] * g[36];
    g[45] = c0z[1] * g[41] + 2 * b10[1] * g[37];
    g[46] = c0z[2] * g[42] + 2 * b10[2] * g[38];
    g[47] = c0z[3] * g[43] + 2 * b10[3] * g[39];
}


void CINTsrg0_2e_2d4d_unrolled(double* g, Rys2eT* bc, CINTEnvVars* envs)
{
    int type_ijkl = ((envs->li_ceil << 6) | (envs->lj_ceil << 4) |
        (envs->lk_ceil << 2) | (envs->ll_ceil));
    switch (type_ijkl) {
    case 0b00000000: _srg0_2d4d_0000(g, bc, envs); return;
    case 0b00000001: _srg0_2d4d_0001(g, bc, envs); return;
    case 0b00000010: _srg0_2d4d_0002(g, bc, envs); return;
    case 0b00000011: _srg0_2d4d_0003(g, bc, envs); return;
    case 0b00000100: _srg0_2d4d_0010(g, bc, envs); return;
    case 0b00000101: _srg0_2d4d_0011(g, bc, envs); return;
    case 0b00000110: _srg0_2d4d_0012(g, bc, envs); return;
    case 0b00001000: _srg0_2d4d_0020(g, bc, envs); return;
    case 0b00001001: _srg0_2d4d_0021(g, bc, envs); return;
    case 0b00001100: _srg0_2d4d_0030(g, bc, envs); return;
    case 0b00010000: _srg0_2d4d_0100(g, bc, envs); return;
    case 0b00010001: _srg0_2d4d_0101(g, bc, envs); return;
    case 0b00010010: _srg0_2d4d_0102(g, bc, envs); return;
    case 0b00010100: _srg0_2d4d_0110(g, bc, envs); return;
    case 0b00010101: _srg0_2d4d_0111(g, bc, envs); return;
    case 0b00011000: _srg0_2d4d_0120(g, bc, envs); return;
    case 0b00100000: _srg0_2d4d_0200(g, bc, envs); return;
    case 0b00100001: _srg0_2d4d_0201(g, bc, envs); return;
    case 0b00100100: _srg0_2d4d_0210(g, bc, envs); return;
    case 0b00110000: _srg0_2d4d_0300(g, bc, envs); return;
    case 0b01000000: _srg0_2d4d_1000(g, bc, envs); return;
    case 0b01000001: _srg0_2d4d_1001(g, bc, envs); return;
    case 0b01000010: _srg0_2d4d_1002(g, bc, envs); return;
    case 0b01000100: _srg0_2d4d_1010(g, bc, envs); return;
    case 0b01000101: _srg0_2d4d_1011(g, bc, envs); return;
    case 0b01001000: _srg0_2d4d_1020(g, bc, envs); return;
    case 0b01010000: _srg0_2d4d_1100(g, bc, envs); return;
    case 0b01010001: _srg0_2d4d_1101(g, bc, envs); return;
    case 0b01010100: _srg0_2d4d_1110(g, bc, envs); return;
    case 0b01100000: _srg0_2d4d_1200(g, bc, envs); return;
    case 0b10000000: _srg0_2d4d_2000(g, bc, envs); return;
    case 0b10000001: _srg0_2d4d_2001(g, bc, envs); return;
    case 0b10000100: _srg0_2d4d_2010(g, bc, envs); return;
    case 0b10010000: _srg0_2d4d_2100(g, bc, envs); return;
    case 0b11000000: _srg0_2d4d_3000(g, bc, envs); return;
    }
    fprintf(stderr, "Dimension error for CINTg0_2e_lj2d4d: iklj = %d %d %d %d",
        (int)envs->li_ceil, (int)envs->lk_ceil,
        (int)envs->ll_ceil, (int)envs->lj_ceil);
}




void CINTinit_int2c2e_EnvVars(CINTEnvVars* envs, int* ng, int* shls,
    int* atm, int natm, int* bas, int nbas, double* env)
{
    envs->natm = natm;
    envs->nbas = nbas;
    envs->atm = atm;
    envs->bas = bas;
    envs->env = env;
    envs->shls = shls;

    const int i_sh = shls[0];
    const int k_sh = shls[1];
    envs->i_l = bas(ANG_OF, i_sh);
    envs->j_l = 0;
    envs->k_l = bas(ANG_OF, k_sh);
    envs->l_l = 0;
    envs->x_ctr[0] = bas(NCTR_OF, i_sh);
    envs->x_ctr[1] = bas(NCTR_OF, k_sh);
    envs->x_ctr[2] = 1;
    envs->x_ctr[3] = 1;
    envs->nfi = (envs->i_l + 1) * (envs->i_l + 2) / 2;
    envs->nfj = 1;
    envs->nfk = (envs->k_l + 1) * (envs->k_l + 2) / 2;
    envs->nfl = 1;
    envs->nf = envs->nfi * envs->nfk;

    envs->ri = env + atm(PTR_COORD, bas(ATOM_OF, i_sh));
    envs->rk = env + atm(PTR_COORD, bas(ATOM_OF, k_sh));


    envs->common_factor = (constants::PI3) * 2 / constants::sqr_pi * CINTcommon_fac_sp(envs->i_l) * CINTcommon_fac_sp(envs->k_l);
    if (env[PTR_EXPCUTOFF] == 0) {
        envs->expcutoff = EXPCUTOFF;
    }
    else {
        envs->expcutoff = std::max((double)MIN_EXPCUTOFF, env[PTR_EXPCUTOFF]);
    }

    envs->gbits = ng[GSHIFT];
    envs->ncomp_e1 = ng[POS_E1];
    envs->ncomp_e2 = ng[POS_E2];
    envs->ncomp_tensor = ng[TENSOR];

    envs->li_ceil = envs->i_l + ng[IINC];
    envs->lj_ceil = 0;
    envs->lk_ceil = envs->k_l + ng[KINC];
    envs->ll_ceil = 0;
    int rys_order = (envs->li_ceil + envs->lk_ceil) / 2 + 1;
    int nrys_roots = rys_order;
    double omega = env[PTR_RANGE_OMEGA];
    if (omega < 0 && rys_order <= 3) {
        nrys_roots *= 2;
    }
    envs->rys_order = rys_order;
    envs->nrys_roots = nrys_roots;

    int dli = envs->li_ceil + 1;
    int dlk = envs->lk_ceil + 1;
    envs->g_stride_i = nrys_roots;
    envs->g_stride_k = nrys_roots * dli;
    envs->g_stride_l = envs->g_stride_k;
    envs->g_size = nrys_roots * dli * dlk;

    envs->aj[0] = 0;
    envs->al[0] = 0;
    envs->rij[0] = envs->ri[0];
    envs->rij[1] = envs->ri[1];
    envs->rij[2] = envs->ri[2];
    envs->rkl[0] = envs->rk[0];
    envs->rkl[1] = envs->rk[1];
    envs->rkl[2] = envs->rk[2];
    envs->g2d_ijmax = envs->g_stride_i;
    envs->g2d_klmax = envs->g_stride_k;
    envs->rkrl[0] = envs->rk[0];
    envs->rkrl[1] = envs->rk[1];
    envs->rkrl[2] = envs->rk[2];
    envs->rirj[0] = envs->ri[0];
    envs->rirj[1] = envs->ri[1];
    envs->rirj[2] = envs->ri[2];
    envs->rx_in_rklrx = envs->rk;
    envs->rx_in_rijrx = envs->ri;

    if (rys_order <= 2) {
        envs->f_g0_2d4d = &CINTg0_2e_2d4d_unrolled;
        if (rys_order != nrys_roots) {
            envs->f_g0_2d4d = &CINTsrg0_2e_2d4d_unrolled;
        }
    }
    else {
        envs->f_g0_2d4d = &CINTg0_2e_2d;
    }
    envs->f_g0_2e = &CINTg0_2e;

    // initialize j_l, j_ctr, nfj because they are used in c2s_sph_1e and
    // CINTg1e_index_xyz
    envs->j_l = envs->k_l;
    envs->nfj = envs->nfk;
    envs->g_stride_j = envs->g_stride_k;
}


/*
 * Note the 3c2e functions takes i,j,k parameters. But we initialize
 * ll_ceil, to reuse g2e_g02d function.
 */
void CINTinit_int3c2e_EnvVars(CINTEnvVars* envs, int* ng, int* shls,
    int* atm, int natm, int* bas, int nbas, double* env)
{
    envs->natm = natm;
    envs->nbas = nbas;
    envs->atm = atm;
    envs->bas = bas;
    envs->env = env;
    envs->shls = shls;

    const int i_sh = shls[0];
    const int j_sh = shls[1];
    const int k_sh = shls[2];
    envs->i_l = bas(ANG_OF, i_sh);
    envs->j_l = bas(ANG_OF, j_sh);
    envs->k_l = bas(ANG_OF, k_sh);
    envs->l_l = 0;
    envs->x_ctr[0] = bas(NCTR_OF, i_sh);
    envs->x_ctr[1] = bas(NCTR_OF, j_sh);
    envs->x_ctr[2] = bas(NCTR_OF, k_sh);
    envs->x_ctr[3] = 1;
    envs->nfi = (envs->i_l + 1) * (envs->i_l + 2) / 2;
    envs->nfj = (envs->j_l + 1) * (envs->j_l + 2) / 2;
    envs->nfk = (envs->k_l + 1) * (envs->k_l + 2) / 2;
    envs->nfl = 1;
    envs->nf = envs->nfi * envs->nfk * envs->nfj;

    envs->ri = env + atm(PTR_COORD, bas(ATOM_OF, i_sh));
    envs->rj = env + atm(PTR_COORD, bas(ATOM_OF, j_sh));
    envs->rk = env + atm(PTR_COORD, bas(ATOM_OF, k_sh));


    envs->common_factor = (constants::PI3) * 2 / constants::sqr_pi * CINTcommon_fac_sp(envs->i_l) * CINTcommon_fac_sp(envs->j_l) * CINTcommon_fac_sp(envs->k_l);
    if (env[PTR_EXPCUTOFF] == 0) {
        envs->expcutoff = EXPCUTOFF;
    }
    else {
        envs->expcutoff = std::max((double)MIN_EXPCUTOFF, env[PTR_EXPCUTOFF]);
    }

    envs->gbits = ng[GSHIFT];
    envs->ncomp_e1 = ng[POS_E1];
    envs->ncomp_e2 = ng[POS_E2];
    envs->ncomp_tensor = ng[TENSOR];

    envs->li_ceil = envs->i_l + ng[IINC];
    envs->lj_ceil = envs->j_l + ng[JINC];
    envs->lk_ceil = 0; // to reuse CINTg0_2e_lj2d4d
    envs->ll_ceil = envs->k_l + ng[KINC];
    int rys_order = (envs->li_ceil + envs->lj_ceil + envs->ll_ceil) / 2 + 1;
    int nrys_roots = rys_order;
    double omega = env[PTR_RANGE_OMEGA];
    if (omega < 0 && rys_order <= 3) {
        nrys_roots *= 2;
    }
    envs->rys_order = rys_order;
    envs->nrys_roots = nrys_roots;

    int dli, dlj, dlk;
    int ibase = envs->li_ceil > envs->lj_ceil;
    if (ibase) {
        dli = envs->li_ceil + envs->lj_ceil + 1;
        dlj = envs->lj_ceil + 1;
    }
    else {
        dli = envs->li_ceil + 1;
        dlj = envs->li_ceil + envs->lj_ceil + 1;
    }
    dlk = envs->ll_ceil + 1;

    envs->g_stride_i = nrys_roots;
    envs->g_stride_k = nrys_roots * dli;
    envs->g_stride_l = nrys_roots * dli;
    envs->g_stride_j = nrys_roots * dli * dlk;
    envs->g_size = nrys_roots * dli * dlk * dlj;

    envs->al[0] = 0;
    envs->rkl[0] = envs->rk[0];
    envs->rkl[1] = envs->rk[1];
    envs->rkl[2] = envs->rk[2];
    envs->g2d_klmax = envs->g_stride_k;
    envs->rkrl[0] = envs->rk[0];
    envs->rkrl[1] = envs->rk[1];
    envs->rkrl[2] = envs->rk[2];
    // in g0_2d rklrx = rkl - rx = 0 => rkl = rx
    envs->rx_in_rklrx = envs->rk;

    if (ibase) {
        envs->g2d_ijmax = envs->g_stride_i;
        envs->rx_in_rijrx = envs->ri;
        envs->rirj[0] = envs->ri[0] - envs->rj[0];
        envs->rirj[1] = envs->ri[1] - envs->rj[1];
        envs->rirj[2] = envs->ri[2] - envs->rj[2];
    }
    else {
        envs->g2d_ijmax = envs->g_stride_j;
        envs->rx_in_rijrx = envs->rj;
        envs->rirj[0] = envs->rj[0] - envs->ri[0];
        envs->rirj[1] = envs->rj[1] - envs->ri[1];
        envs->rirj[2] = envs->rj[2] - envs->ri[2];
    }

    if (rys_order <= 2) {
        envs->f_g0_2d4d = &CINTg0_2e_2d4d_unrolled;
        if (rys_order != nrys_roots) {
            envs->f_g0_2d4d = &CINTsrg0_2e_2d4d_unrolled;
        }
    }
    else if (ibase) {
        envs->f_g0_2d4d = &CINTg0_2e_il2d4d;
    }
    else {
        envs->f_g0_2d4d = &CINTg0_2e_lj2d4d;
    }
    envs->f_g0_2e = &CINTg0_2e;
}