#pragma once
#include "wfn_class.h"
#include "SALTED_predictor.h"

int density_fit(const WFN& wavy, const std::string auxname);

int lrys_schmidt(int nroots, double x, double lower, vec &roots, vec &weights);
int R_lsmit(long double* cs, long double* fmt_ints, int n);
int lrys_jacobi(int n, double x, double lower, vec& roots, vec& weights);
int lrys_laguerre(int n, double x, double lower, vec& roots, vec& weights);
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

	ivec get_atm() { return _atm; };
	ivec get_bas() { return _bas; };
	vec get_env() { return _env; };

    int get_nbas() { return nbas; };  // Number of basis segments
	int get_nao() { return nao; };    // Number of atomic orbitals  
    int get_natoms() { return ncen; };
};

typedef struct {
    std::array<double, 3> rij;
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
    ivec atm;
    ivec bas;
    vec env;
    ivec shls;
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
    vec rx_in_rijrx;
    vec rx_in_rklrx;

    //double* ri;
    vec ri = vec(3);
    vec rj = vec(3);
    vec rk = vec(3);
    // in int2e or int3c2e, the coordinates of the fourth shell
    // in int1e_grids, the pointer for the grids coordinates
    union { double* rl; double* grids; };

    int(*f_g0_2e)(double *, vec &, vec &, double, Env*);
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
    vec rkl = vec(3);
};

int int2c2e_sph(vec &out, ivec &dims, ivec &shls, ivec &atm,
    int natm, ivec &bas, int nbas, vec &env,
    Opt* opt, double* cache);
int int3c2e_sph(vec &out, ivec &dims, ivec &shls, ivec &atm,
    int natm, ivec &bas, int nbas, vec &env,
    Opt* opt, double* cache);
void GTOint2c(int (*intor)(vec &, ivec&, ivec&, ivec&, int, ivec&, int, vec&, Opt*, double*),
    vec &mat, int comp, int hermi,
    ivec& shls_slice, ivec& ao_loc, Opt* opt,
    ivec& atm, int natm, ivec& bas, int nbas, vec& env);
void GTOnr3c_fill_s1(int (*intor)(vec &, ivec&, ivec&, ivec&, int, ivec&, int, vec&, Opt*, double*), vec &out, double* buf,
    int comp, int jobid,
    ivec& shls_slice, ivec &ao_loc, Opt* cintopt,
    ivec &atm, int natm, ivec &bas, int nbas, vec &env);
void GTOnr3c_drv(int (*intor)(vec &, ivec &, ivec &, ivec &, int, ivec&, int, vec &, Opt*, double*),
    void (*fill)(int (*intor)(vec &, ivec &, ivec &, ivec &, int, ivec&, int, vec &, Opt*, double*),
        vec &, double*, int, int, ivec &, ivec &, Opt*, ivec&, int, ivec&, int, vec&),
    vec &eri, int comp,
   ivec &shls_slice, ivec &ao_loc, Opt* cintopt,
    ivec &atm, int natm, ivec &bas, int nbas, vec &env);

Opt int3c2e_optimizer(ivec &atm, int natm, ivec  &bas, int nbas, vec &env);
Opt int2c2e_optimizer(ivec &atm, int natm, ivec &bas, int nbas, vec &env);
int fixed_density_fit_test();