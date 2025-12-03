#pragma once
#include "convenience.h"
#include "constants.h"
#include "wfn_class.h"


// global parameters in env
// Overall cutoff for integral prescreening, value needs to be ~ln(threshold)
#define PTR_EXPCUTOFF           0
// R_C of (r-R_C) in dipole, GIAO operators
#define PTR_COMMON_ORIG         1
// R_O in 1/|r-R_O|
#define PTR_RINV_ORIG           4
// ZETA parameter for Gaussian charge distribution (Gaussian nuclear model)
#define PTR_RINV_ZETA           7
// omega parameter in range-separated coulomb operator
// LR interaction: erf(omega*r12)/r12 if omega > 0
// SR interaction: erfc(omega*r12)/r12 if omega < 0
#define PTR_RANGE_OMEGA         8
// Yukawa potential and Slater-type geminal e^{-zeta r}
#define PTR_F12_ZETA            9
// Gaussian type geminal e^{-zeta r^2}
#define PTR_GTG_ZETA            10
#define NGRIDS                  11
#define PTR_GRIDS               12
#define PTR_ENV_START           20


// slots of atm
#define CHARGE_OF       0
#define PTR_COORD       1
#define NUC_MOD_OF      2
#define PTR_ZETA        3
#define PTR_FRAC_CHARGE 4
#define RESERVE_ATMSLOT 5
#define ATM_SLOTS       6


// slots of bas
#define ATOM_OF         0
#define ANG_OF          1
#define NPRIM_OF        2
#define NCTR_OF         3
#define KAPPA_OF        4
#define PTR_EXP         5
#define PTR_COEFF       6
#define RESERVE_BASLOT  7
#define BAS_SLOTS       8

// slots of gout
#define POSX            0
#define POSY            1
#define POSZ            2
#define POS1            3

// For 2-electron integral with two spin operators
#define POSXX           0
#define POSYX           1
#define POSZX           2
#define POS1X           3
#define POSXY           4
#define POSYY           5
#define POSZY           6
#define POS1Y           7
#define POSXZ           8
#define POSYZ           9
#define POSZZ           10
#define POS1Z           11
#define POSX1           12
#define POSY1           13
#define POSZ1           14
#define POS11           15

// other boundaries
#define MXRYSROOTS      32 // > ANG_MAX*2+1 for 4c2e
#define ANG_MAX         15 // l = 0..15
#define LMAX1           16 // > ANG_MAX
#define CART_MAX        136 // > (ANG_MAX*(ANG_MAX+1)/2)
#define SHLS_MAX        1048576
#define NPRIM_MAX       64
#define NCTR_MAX        64

#define POINT_NUC       1
#define GAUSSIAN_NUC    2
#define FRAC_CHARGE_NUC 3

// ng[*]
#define IINC            0
#define JINC            1
#define KINC            2
#define LINC            3
#define GSHIFT          4
#define POS_E1          5
#define POS_E2          6
#define SLOT_RYS_ROOTS  6
#define TENSOR          7

#define EXPCUTOFF       60
#ifndef MIN_EXPCUTOFF
// ~ 1e-15
#define MIN_EXPCUTOFF   40
#endif

#define OF_CMPLX        2

#define GRID_BLKSIZE    104


#define bas(SLOT,I)     bas[BAS_SLOTS * (I) + (SLOT)]
#define atm(SLOT,I)     atm[ATM_SLOTS * (I) + (SLOT)]

#define DEF_GXYZ(type, G, GX, GY, GZ) \
        type *GX = G; \
        type *GY = G + envs->g_size; \
        type *GZ = G + envs->g_size * 2

#define MALLOC_INSTACK(var, n) \
        var = (double *)(((uintptr_t)cache + 7) & (~(uintptr_t)7)); \
        cache = (double *)(var + (n));

#define MALLOC_INSTACK_INT(var, n) \
        var = (int *)(((uintptr_t)cache + 7) & (~(uintptr_t)7)); \
        cache = (double *)(var + (n));

#define MALLOC_INSTACK_PAIRDATA(var, n) \
        var = (PairData *)(((uintptr_t)cache + 7) & (~(uintptr_t)7)); \
        cache = (double *)(var + (n));


#define BLKSIZE 8

//Temporary basis set format
struct LibCintBasis {
    vec coefficients;
    vec exponents;
    int env_idx;
    ivec shellcount;
    ivec shelltypes;
};


class Int_Params {
private:
    unsigned long long ncen = 0;           // Number of centers (atoms)
    int wfn_origin = 0;
    std::vector<atom> atoms;

    std::map<int, LibCintBasis> basis_sets;  //Maps atomic number -> (vec(coefs), vec(expon))

    ivec _atm;              // Flattened ATM array
    ivec _bas;              // Flattened BAS array
    vec _env;

    unsigned long long nao = 0;
    unsigned long long nbas = 0;
    void calc_integration_parameters();
    void collect_basis_data();
    //   void collect_basis_data_from_gbw();
       //void collect_basis_data_internal();

    void populate_atm();
    void populate_bas();
    void populate_env();

    double gaussian_int(int n, double exp) {
        double n1 = (n + 1) * 0.5;
        return std::tgamma(n1) / (2.0 * pow(exp, n1));
    }
public:
    Int_Params();
    Int_Params(const WFN& wavy);
    Int_Params operator+(const Int_Params& other);

    // Getters
    int* get_ptr_atm() { return _atm.data(); };
    int* get_ptr_bas() { return _bas.data(); };
    double* get_ptr_env() { return _env.data(); };

    ivec get_atm() { return _atm; };
    ivec get_bas() { return _bas; };
    vec get_env() { return _env; };

    unsigned long long get_nbas() { return nbas; };  // Number of basis segments
    unsigned long long get_nao() { return nao; };    // Number of atomic orbitals  
    unsigned long long get_natoms() { return ncen; };

    vec normalize_gto(vec coef, const vec& exp, const int l);

    void print_data(std::string name); //FOR DEBUG PURPOSES
};



struct PairData {
    double rij[3];
    double eij;
    double cceij;
};

struct CINTOpt {
    int** index_xyz_array; // LMAX1**4 pointers to index_xyz
    int** non0ctr;
    int** sortedidx;
    int nbas;
    double** log_max_coeff;
    PairData** pairdata;  // NULL indicates not-initialized, NO_VALUE can be skipped
};

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


struct CINTEnvVars {
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

    int (*f_g0_2e)(double*, double*, double*, double, CINTEnvVars*);
    void (*f_g0_2d4d)(double *, Rys2eT*, CINTEnvVars*);
    void (*f_gout)(double* , double* , int* , CINTEnvVars* , int );
    CINTOpt* opt;

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

ivec make_loc(ivec& bas, int nbas);


/*
 * to optimize memory copy in cart2sph.c, remove the common factor for s
 * and p function in cart2sph
 */
double CINTcommon_fac_sp(int l);
