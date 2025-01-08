#pragma once
#include "wfn_class.h"
#include "SALTED_predictor.h"

int density_fit(const WFN& wavy, const std::string auxname);

int lrys_schmidt(int nroots, double x, double lower, double* roots, double* weights);
int R_lsmit(long double* cs, long double* fmt_ints, int n);
int lrys_jacobi(int n, double x, double lower, double* roots, double* weights);
int lrys_laguerre(int n, double x, double lower, double* roots, double* weights);
void lnaive_jacobi_moments(int n, double t, double lower, long double* mus);


class Int_Params {
private:
    int ncen = 0;           // Number of centers (atoms)
    int wfn_origin = 0;
    std::vector<atom> atoms;

    std::map<int, vec2> basis_sets;  //Maps atomic number -> (vec(coefs), vec(expon))
    
    ivec _atm;              // Flattened ATM array
    ivec _bas;              // Flattened BAS array
    vec _env;

    int nao = 0;
    int nbas = 0;
    void calc_integration_parameters();
    void collect_basis_data();
 //   void collect_basis_data_from_gbw();
	//void collect_basis_data_internal();

    void populate_atm();
	void populate_bas();
	void populate_env();

    vec normalize_gto(vec coef, vec exp, int l);
    double gaussian_int(int n, double exp) {
        double n1 = (n + 1) * 0.5;
        return std::tgamma(n1) / (2.0 * pow(exp, n1));
    }
public:
    Int_Params();
	Int_Params(const WFN& wavy);
    Int_Params(WFN& wavy, const std::string auxname);
	Int_Params operator+(const Int_Params& other);

    // Getters
	int* get_ptr_atm() { return _atm.data(); };
	int* get_ptr_bas() { return _bas.data(); };
	double* get_ptr_env() { return _env.data(); };

    int get_nbas() { return nbas; };  // Number of basis segments
	int get_nao() { return nao; };    // Number of atomic orbitals  
    int get_natoms() { return ncen; };
};

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

int int2c2e_sph(double* out, int* dims, int* shls, int* atms,
    int natms, int* bas, int nbas, double* env, 
    Opt* opt, double* cache);
int int3c2e_sph(double* out, int* dims, int* shls, int* atm,
    int natm, int* bas, int nbas, double* env,
    Opt* opt, double* cache);
void GTOint2c(int (*intor)(double*, int*, int*, int*, int, int*, int, double*, Opt*, double*),
    double* mat, int comp, int hermi,
    int* shls_slice, int* ao_loc, Opt* opt,
    int* atm, int natm, int* bas, int nbas, double* env);
void GTOnr3c_fill_s1(int (*intor)(double*, int*, int*, int*, int, int*, int, double*, Opt*, double*), double* out, double* buf,
    int comp, int jobid,
    int* shls_slice, int* ao_loc, Opt* cintopt,
    int* atm, int natm, int* bas, int nbas, double* env);
void GTOnr3c_drv(int (*intor)(double*, int*, int*, int*, int, int*, int, double*, Opt*, double*),
    void (*fill)(int (*intor)(double*, int*, int*, int*, int, int*, int, double*, Opt*, double*),
        double*, double*, int, int, int*, int*, Opt*, int*, int, int*, int, double*),
    double* eri, int comp,
    int* shls_slice, int* ao_loc, Opt* cintopt,
    int* atm, int natm, int* bas, int nbas, double* env);

Opt int3c2e_optimizer(int* atm, int natm, int* bas, int nbas, double* env);
Opt int2c2e_optimizer(int* atm, int natm, int* bas, int nbas, double* env);
int fixed_density_fit_test();