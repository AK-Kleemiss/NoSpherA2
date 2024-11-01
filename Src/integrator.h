#pragma once
#include "wfn_class.h"
#include "SALTED_predictor.h"

int density_fit(const WFN& wavy);

int lrys_schmidt(int nroots, double x, double lower, double* roots, double* weights);
int R_lsmit(long double* cs, long double* fmt_ints, int n);
int lrys_jacobi(int n, double x, double lower, double* roots, double* weights);
int lrys_laguerre(int n, double x, double lower, double* roots, double* weights);
void lnaive_jacobi_moments(int n, double t, double lower, long double* mus);

typedef struct {
    double rij[3];
    double eij;
    double cceij;
} PairData;

typedef struct {
    int** index_xyz_array; // LMAX1**4 pointers to index_xyz
    int** non0ctr;
    int** sortedidx;
    int nbas;
    double** log_max_coeff;
    PairData** pairdata;  // NULL indicates not-initialized, NO_VALUE can be skipped
} Opt;

typedef struct {
    double c00x[32];
    double c00y[32];
    double c00z[32];
    double c0px[32];
    double c0py[32];
    double c0pz[32];
    double  b01[32];
    double  b00[32];
    double  b10[32];
} Rys2eT;

struct Env {
    int* atm;
    int* bas;
    double* env;
    int* shls;
    int natm;
    int nbas;

    int i_l;
    int j_l;
    int k_l;
    int l_l;
    int nfi;  // number of cartesian components
    int nfj;
    // in int1e_grids, the grids_offset and the number of grids
    union { int nfk; int grids_offset; };
    union { int nfl; int ngrids; };
    int nf;  // = nfi*nfj*nfk*nfl;
    int rys_order; // = nrys_roots for regular ERIs. can be nrys_roots/2 for SR ERIs
    int x_ctr[4];

    int gbits;
    int ncomp_e1; // = 1 if spin free, = 4 when spin included, it
    int ncomp_e2; // corresponds to POSX,POSY,POSZ,POS1, see cint.h
    int ncomp_tensor; // e.g. = 3 for gradients

    /* values may diff based on the g0_2d4d algorithm */
    int li_ceil; // power of x, == i_l if nabla is involved, otherwise == i_l
    int lj_ceil;
    int lk_ceil;
    int ll_ceil;
    int g_stride_i; // nrys_roots * shift of (i++,k,l,j)
    int g_stride_k; // nrys_roots * shift of (i,k++,l,j)
    int g_stride_l; // nrys_roots * shift of (i,k,l++,j)
    int g_stride_j; // nrys_roots * shift of (i,k,l,j++)
    int nrys_roots;
    int g_size;  // ref to cint2e.c g = malloc(sizeof(double)*g_size)

    int g2d_ijmax;
    int g2d_klmax;
    double common_factor;
    double expcutoff;
    double rirj[3]; // diff by sign in different g0_2d4d algorithm
    double rkrl[3];
    double* rx_in_rijrx;
    double* rx_in_rklrx;

    double* ri;
    double* rj;
    double* rk;
    // in int2e or int3c2e, the coordinates of the fourth shell
    // in int1e_grids, the pointer for the grids coordinates
    union { double* rl; double* grids; };

    int(*f_g0_2e)(double*, double*, double*, double, Env*);
    void (*f_g0_2d4d)(double*, Rys2eT*, Env*);
    void (*f_gout)(double*, double*, int*, Env*, int);
    Opt* opt;

    /* values are assigned during calculation */
    int* idx;
    double ai[1];
    double aj[1];
    double ak[1];
    double al[1];
    double fac[1];
    double rij[3];
    double rkl[3];
};