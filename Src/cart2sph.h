/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * Cartisen GTO to spheric or spinor GTO transformation
 */
#pragma once
/*************************************************
 *
 * transform matrix
 *
 *************************************************/
#include <complex.h>
#include "integration_params.h"

double* a_bra_cart2spheric(double* gsph, int nket, double* gcart, int l);

double* a_ket_cart2spheric(double* gsph, double* gcart,
    int lds, int nbra, int l);

// transform s function from cartesian to spheric
double* s_bra_cart2spheric(double* gsph, int nket, double* gcart, int l);
double* s_ket_cart2spheric(double* gsph, double* gcart,
    int lds, int nbra, int l);

// transform p function from cartesian to spheric
double* p_bra_cart2spheric(double* gsph, int nket, double* gcart, int l);
double* p_ket_cart2spheric(double* gsph, double* gcart,
    int lds, int nbra, int l);

// transform d function from cartesian to spheric
double* d_bra_cart2spheric(double* gsph, int nket, double* gcart, int l);
double* d_ket_cart2spheric(double* gsph, double* gcart,
    int lds, int nbra, int l);

// transform f function from cartesian to spheric
double* f_bra_cart2spheric(double* gsph, int nket, double* gcart, int l);
double* f_ket_cart2spheric(double* gsph, double* gcart,
    int lds, int nbra, int l);

// transform g function from cartesian to spheric
double* g_bra_cart2spheric(double* gsph, int nket, double* gcart, int l);
double* g_ket_cart2spheric(double* gsph, double* gcart,
    int lds, int nbra, int l);

// [ca, cb] * [ 1+1j*z, y+1j*x]
//            [-y+1j*x, 1-1j*z]
void a_bra_cart2spinor_si(double* gspR, double* gspI,
    double* gx, double* gy, double* gz, double* g1,
    int nket, int kappa, int l);

void a_bra_cart2spinor_sf(double* gspR, double* gspI,
    double* gx, double* gy, double* gz, double* g1,
    int nket, int kappa, int l);

void a_bra1_cart2spinor_si(double* gspR, double* gspI,
    double* gx, double* gy, double* gz, double* g1,
    int ngrids, int nket, int kappa, int l);

void a_bra1_cart2spinor_sf(double* gspR, double* gspI,
    double* gx, double* gy, double* gz, double* g1,
    int ngrids, int nket, int kappa, int l);

void a_bra1_cart2spinor_zi(double* gspR, double* gspI,
    double* gx, double* gy, double* gz, double* g1,
    int ngrids, int nket, int kappa, int l);

void a_bra1_cart2spinor_zf(double* gspR, double* gspI,
    double* gx, double* gy, double* gz, double* g1,
    int ngrids, int nket, int kappa, int l);

void a_ket_cart2spinor(double* gspR, double* gspI,
    double* gcartR, double* gcartI,
    int nbra, int kappa, int l);

// with phase "i"
void a_iket_cart2spinor(double* gspR, double* gspI,
    double* gcartR, double* gcartI,
    int nbra, int kappa, int l);

void a_ket1_cart2spinor(double* gspR, double* gspI,
    double* gcartR, double* gcartI,
    int nbra, int counts, int kappa, int l);

// with phase "i"
void a_iket1_cart2spinor(double* gspR, double* gspI,
    double* gcartR, double* gcartI,
    int nbra, int counts, int kappa, int l);



/*
 * return the address of gemm results, for s,p function, results ==
 * input, so return the input address optimize
 */
inline double* (*c2s_bra_sph[])(double* gsph, int nket, double* gcart, int l) = {
        s_bra_cart2spheric,
        p_bra_cart2spheric,
        d_bra_cart2spheric,
        f_bra_cart2spheric,
        g_bra_cart2spheric,
        a_bra_cart2spheric,
        a_bra_cart2spheric,
        a_bra_cart2spheric,
        a_bra_cart2spheric,
        a_bra_cart2spheric,
        a_bra_cart2spheric,
        a_bra_cart2spheric,
        a_bra_cart2spheric,
        a_bra_cart2spheric,
        a_bra_cart2spheric,
        a_bra_cart2spheric,
};

inline double* (*c2s_ket_sph[])(double* gsph, double* gcart, int lds, int nbra, int l) = {
s_ket_cart2spheric,
p_ket_cart2spheric,
d_ket_cart2spheric,
f_ket_cart2spheric,
g_ket_cart2spheric,
a_ket_cart2spheric,
a_ket_cart2spheric,
a_ket_cart2spheric,
a_ket_cart2spheric,
a_ket_cart2spheric,
a_ket_cart2spheric,
a_ket_cart2spheric,
a_ket_cart2spheric,
a_ket_cart2spheric,
a_ket_cart2spheric,
a_ket_cart2spheric,
};

void c2s_sph_1e(double *opij, double *gctr, int *dims, CINTEnvVars *envs, double *cache);
void c2s_sph_2e1(double *fijkl, double *gctr, int *dims, CINTEnvVars *envs, double *cache);
//void c2s_sph_2e2();

void c2s_cart_1e(double *opij, double *gctr, int *dims, CINTEnvVars *envs, double *cache);
void c2s_cart_2e1(double *fijkl, double *gctr, int *dims, CINTEnvVars *envs, double *cache);
void c2s_cart_2e2();

//void c2s_sf_1e(double complex *opij, double *gctr, int *dims, CINTEnvVars *envs, double *cache);
//void c2s_sf_1ei(double complex *opij, double *gctr, int *dims, CINTEnvVars *envs, double *cache);
//
//void c2s_si_1e(double complex *opij, double *gctr, int *dims, CINTEnvVars *envs, double *cache);
//void c2s_si_1ei(double complex *opij, double *gctr, int *dims, CINTEnvVars *envs, double *cache);

void c2s_sph_1e_grids(double *out, double *gctr, int *dims, CINTEnvVars *envs, double *cache);
void c2s_cart_1e_grids(double *out, double *gctr, int *dims, CINTEnvVars *envs, double *cache);

//void c2s_sf_1e_grids(double complex *out, double *gctr, int *dims, CINTEnvVars *envs, double *cache);
//void c2s_sf_1e_gridsi(double complex *out, double *gctr, int *dims, CINTEnvVars *envs, double *cache);
//void c2s_si_1e_grids(double complex *out, double *gctr, int *dims, CINTEnvVars *envs, double *cache);
//void c2s_si_1e_gridsi(double complex *out, double *gctr, int *dims, CINTEnvVars *envs, double *cache);

void c2s_sf_2e1(double *opij, double *gctr, int *dims, CINTEnvVars *envs, double *cache);
void c2s_sf_2e1i(double *opij, double *gctr, int *dims, CINTEnvVars *envs, double *cache);

//void c2s_sf_2e2(double complex *fijkl, double *opij, int *dims, CINTEnvVars *envs, double *cache);
//void c2s_sf_2e2i(double complex *fijkl, double *opij, int *dims, CINTEnvVars *envs, double *cache);

void c2s_si_2e1(double *opij, double *gctr, int *dims, CINTEnvVars *envs, double *cache);
void c2s_si_2e1i(double *opij, double *gctr, int *dims, CINTEnvVars *envs, double *cache);

//void c2s_si_2e2(double complex *fijkl, double *opij, int *dims, CINTEnvVars *envs, double *cache);
//void c2s_si_2e2i(double complex *fijkl, double *opij, int *dims, CINTEnvVars *envs, double *cache);

void c2s_sph_3c2e1(double *fijkl, double *gctr, const int *dims, CINTEnvVars *envs, double *cache);
void c2s_cart_3c2e1(double *fijkl, double *gctr, int *dims, CINTEnvVars *envs, double *cache);
void c2s_sph_3c2e1_ssc(double *fijkl, double *gctr, int *dims, CINTEnvVars *envs, double *cache);

//void c2s_sf_3c2e1(double complex *opijk, double *gctr, int *dims, CINTEnvVars *envs, double *cache);
//void c2s_sf_3c2e1i(double complex *opijk, double *gctr, int *dims, CINTEnvVars *envs, double *cache);
//void c2s_si_3c2e1(double complex *opijk, double *gctr, int *dims, CINTEnvVars *envs, double *cache);
//void c2s_si_3c2e1i(double complex *opijk, double *gctr, int *dims, CINTEnvVars *envs, double *cache);
//void c2s_sf_3c2e1_ssc(double complex *opijk, double *gctr, int *dims, CINTEnvVars *envs, double *cache);
//void c2s_sf_3c2e1i_ssc(double complex *opijk, double *gctr, int *dims, CINTEnvVars *envs, double *cache);
//void c2s_si_3c2e1_ssc(double complex *opijk, double *gctr, int *dims, CINTEnvVars *envs, double *cache);
//void c2s_si_3c2e1i_ssc(double complex *opijk, double *gctr, int *dims, CINTEnvVars *envs, double *cache);

void c2s_sph_3c1e(double *fijkl, double *gctr, int *dims, CINTEnvVars *envs, double *cache);
void c2s_cart_3c1e(double *fijkl, double *gctr, int *dims, CINTEnvVars *envs, double *cache);

void c2s_dset0(double *out, const int *dims, int *counts);
//void c2s_zset0(double complex *out, int *dims, int *counts);
//void c2s_grids_dset0(double *out, int *dims, int *counts);
//void c2s_grids_zset0(double complex *out, int *dims, int *counts);

/*************************************************
 *
 * transform vectors
 *
 *************************************************/
//void c2s_sph_vec(double *sph, double *cart, int l, int nvec);

