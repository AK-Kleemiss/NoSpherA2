/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * Cartisen GTO to spheric or spinor GTO transformation
 */

/*************************************************
 *
 * transform matrix
 *
 *************************************************/
#include <complex.h>
#include "integration_params.h"


double* (*c2s_bra_sph[])(double* gsph, int nket, double* gcart, int l);
double* (*c2s_ket_sph[])(double* gsph, double* gcart, int lds, int nbra, int l);
double* (*c2s_ket_sph1[])(double* gsph, double* gcart, int lds, int nbra, int l);

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

