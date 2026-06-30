#include "pch.h"
#include "convenience.h"

#include "libCintKernels.h"

#define BLKSIZE 8
/*
 * out[naoi,naoj,naok,comp] in F-order
 */
void GTOnr3c_fill_s1(
    int (*intor)(double*, int*, int*, int*, int, int*, int, double*, libcint::CINTOpt*, double*),
    double* out, double* buf, int comp, int jobid, int* shls_slice, int* ao_loc, libcint::CINTOpt* cintopt,
    int* atm, int natm, int* bas, int nbas, double* env)
{
    const int ish0 = shls_slice[0];
    const int ish1 = shls_slice[1];
    const int jsh0 = shls_slice[2];
    const int jsh1 = shls_slice[3];
    const int ksh0 = shls_slice[4];
    const int ksh1 = shls_slice[5];
    const int nksh = ksh1 - ksh0;

    const int ksh = jobid % nksh + ksh0;
    const int jstart = jobid / nksh * BLKSIZE + jsh0;
    const int jend = std::min(jstart + BLKSIZE, jsh1);
    if (jstart >= jend)
    {
        return;
    }

    const size_t naoi = ao_loc[ish1] - ao_loc[ish0];
    const size_t naoj = ao_loc[jsh1] - ao_loc[jsh0];
    const size_t naok = ao_loc[ksh1] - ao_loc[ksh0];
    int dims[] = { (int)naoi, (int)naoj, (int)naok };

    const int k0 = ao_loc[ksh] - ao_loc[ksh0];
    out += naoi * naoj * k0;

    int ish, jsh, i0, j0;
    int shls[3] = { 0, 0, ksh };

    for (jsh = jstart; jsh < jend; jsh++)
    {
        for (ish = ish0; ish < ish1; ish++)
        {
            shls[0] = ish;
            shls[1] = jsh;
            i0 = ao_loc[ish] - ao_loc[ish0];
            j0 = ao_loc[jsh] - ao_loc[jsh0];
            (*intor)(out + j0 * naoi + i0, dims, shls, atm, natm, bas, nbas, env,
                cintopt, buf);
        }
    }
}



int GTOmax_shell_dim(const int *ao_loc, const int *shls_slice, int ncenter)
{
    int i;
    int i0 = shls_slice[0];
    int i1 = shls_slice[1];
    int di = 0;
    for (i = 1; i < ncenter; i++)
    {
        i0 = std::min(i0, shls_slice[i * 2]);
        i1 = std::max(i1, shls_slice[i * 2 + 1]);
    }
    for (i = i0; i < i1; i++)
    {
        di = std::max(di, ao_loc[i + 1] - ao_loc[i]);
    }
    return di;
}

size_t GTOmax_cache_size(
    int (*intor)(double*, int*, int*, int*, int, int*, int, double*, libcint::CINTOpt*, double*),
    int* shls_slice, int ncenter,
    int* atm, int natm, int* bas, int nbas, double* env)
{
    int i;
    int i0 = shls_slice[0];
    int i1 = shls_slice[1];
    for (i = 1; i < ncenter; i++)
    {
        i0 = std::min(i0, shls_slice[i * 2]);
        i1 = std::max(i1, shls_slice[i * 2 + 1]);
    }
    int (*f)(double*, int*, int*, int*, int, int*, int, double*, libcint::CINTOpt*, double*) = (int (*)(double*, int*, int*, int*, int, int*, int, double*, libcint::CINTOpt*, double*))intor;
    int cache_size = 0;
    int n;
    int shls[4];
    for (i = i0; i < i1; i++)
    {
        shls[0] = i;
        shls[1] = i;
        shls[2] = i;
        shls[3] = i;
        n = (*f)(NULL, NULL, shls, atm, natm, bas, nbas, env, NULL, NULL);
        cache_size = std::max(cache_size, n);
    }
    return cache_size;
}


void GTOnr3c_drv(
    int (*intor)(double*, int*, int*, int*, int, int*, int, double*, libcint::CINTOpt*, double*),
    void (*fill)(int (*intor)(double*, int*, int*, int*, int, int*, int, double*, libcint::CINTOpt*, double*), double*, double*, int, int, int*, int*, libcint::CINTOpt*, int*, int, int*, int, double*),
    double* eri, int comp, int* shls_slice, int* ao_loc, libcint::CINTOpt* cintopt,
    int* atm, int natm, int* bas, int nbas, double* env)
{
    const int ish0 = shls_slice[0];
    const int ish1 = shls_slice[1];
    const int jsh0 = shls_slice[2];
    const int jsh1 = shls_slice[3];
    const int ksh0 = shls_slice[4];
    const int ksh1 = shls_slice[5];
    const int nish = ish1 - ish0;
    const int njsh = jsh1 - jsh0;
    const int nksh = ksh1 - ksh0;
    const int di = GTOmax_shell_dim(ao_loc, shls_slice, 3);
    const int cache_size = GTOmax_cache_size(intor, shls_slice, 3,
        atm, natm, bas, nbas, env);
    const int njobs = (std::max(nish, njsh) / BLKSIZE + 1) * nksh;

#pragma omp parallel
    {
        int jobid;
        double *buf = (double *)malloc(sizeof(double) * (di * di * di * comp + cache_size));
#pragma omp for nowait schedule(dynamic)
        for (jobid = 0; jobid < njobs; jobid++)
        {
            (*fill)(intor, eri, buf, comp, jobid, shls_slice, ao_loc,
                cintopt, atm, natm, bas, nbas, env);
        }
        free(buf);
    }
}


/*
 * mat(naoi,naoj,comp) in F-order
 */
void GTOint2c(int (*intor)(double*, int*, int*, int*, int, int*, int, double*, libcint::CINTOpt*, double*), double* mat, int comp, int hermi,
    int* shls_slice, int* ao_loc, libcint::CINTOpt* opt,
    int* atm, int natm, int* bas, int nbas, double* env)
{
    const int ish0 = shls_slice[0];
    const int ish1 = shls_slice[1];
    const int jsh0 = shls_slice[2];
    const int jsh1 = shls_slice[3];
    const int nish = ish1 - ish0;
    const int njsh = jsh1 - jsh0;
    const size_t naoi = ao_loc[ish1] - ao_loc[ish0];
    const size_t naoj = ao_loc[jsh1] - ao_loc[jsh0];
    const int cache_size = GTOmax_cache_size(intor, shls_slice, 2,
        atm, natm, bas, nbas, env);
#pragma omp parallel
    {
        int dims[] = { (int)naoi, (int)naoj };
        int ish, jsh, ij, i0, j0;
        int shls[2];
        double *cache = (double *)malloc(sizeof(double) * cache_size);
#pragma omp for schedule(dynamic, 4)
        for (ij = 0; ij < nish * njsh; ij++)
        {
            ish = ij / njsh;
            jsh = ij % njsh;

            ish += ish0;
            jsh += jsh0;
            shls[0] = ish;
            shls[1] = jsh;
            i0 = ao_loc[ish] - ao_loc[ish0];
            j0 = ao_loc[jsh] - ao_loc[jsh0];
            (*intor)(mat + j0 * naoi + i0, dims, shls,
                atm, natm, bas, nbas, env, opt, cache);
        }
        free(cache);
    }
}


enum COORDINATE_TYPE {
    SPH = 1,
    CART = 2
};

template<COORDINATE_TYPE CT>
ivec make_loc(ivec &bas, int nbas) {
    ivec dims(nbas, 0);
    // Calculate (2*l + 1) * nctr for spherical harmonics
    for (size_t i = 0; i < nbas; i++)
    {
        if constexpr (CT == COORDINATE_TYPE::CART) {
            int l = bas(ANG_OF, i);
            dims[i] = ((l + 1) * (l + 2)) / 2 * bas(NCTR_OF, i); // Number of cartesian functions for given l
            continue;
        }
        else if constexpr (CT == COORDINATE_TYPE::SPH) {
            dims[i] = (2 * bas(ANG_OF, i) + 1) * bas(NCTR_OF, i);
        }
    }

    // Create the ao_loc array
    ivec ao_loc(nbas + 1, 0);

    // Compute the cumulative sum
    std::partial_sum(dims.begin(), dims.end(), ao_loc.begin() + 1);

    return ao_loc;
}
template ivec make_loc<COORDINATE_TYPE::CART>(ivec &bas, int nbas);
template ivec make_loc<COORDINATE_TYPE::SPH>(ivec &bas, int nbas);

extern "C" {
    extern libcint::CINTOptimizerFunction int3c2e_optimizer;
    extern libcint::CINTIntegralFunction int3c2e_sph;
    extern libcint::CINTIntegralFunction int3c2e_cart;

    extern libcint::CINTOptimizerFunction int2c2e_optimizer;
    extern libcint::CINTIntegralFunction int2c2e_sph;
    extern libcint::CINTIntegralFunction int2c2e_cart;

    extern libcint::CINTOptimizerFunction int3c1e_optimizer;
    extern libcint::CINTIntegralFunction int3c1e_sph;
}

#define ADD_FUNCS_TO_KERNEL(Kernel, optimizer_func, driver_func, integral_func, coord_type) \
    void Kernel::optimizer(libcint::CINTOpt*& opt, \
        int* atm, int nat, int* bas, int nbas, double* env) { \
        optimizer_func(&opt, atm, nat, bas, nbas, env); \
    } \
    void Kernel::drv(double* out, int comp, int* shl_slice, int* aoloc, \
        libcint::CINTOpt* opt, int* atm, int nat, int* bas, int nbas, double* env) { \
        driver_func(integral_func, out, comp, 0, shl_slice, aoloc, opt, atm, nat, bas, nbas, env); \
    } \
    ivec Kernel::gen_loc(ivec& bas, int nbas) { \
        return make_loc<coord_type>(bas, nbas); \
    }

#define ADD_FUNCS_TO_KERNEL_NOOPT(Kernel, driver_func, integral_func, coord_type) \
    void Kernel::optimizer(libcint::CINTOpt*& opt, \
        int* atm, int nat, int* bas, int nbas, double* env) { \
        opt = nullptr; \
    } \
    void Kernel::drv(double* out, int comp, int* shl_slice, int* aoloc, \
        libcint::CINTOpt* opt, int* atm, int nat, int* bas, int nbas, double* env) { \
        driver_func(integral_func, out, comp, 0, shl_slice, aoloc, opt, atm, nat, bas, nbas, env); \
    } \
    ivec Kernel::gen_loc(ivec& bas, int nbas) { \
        return make_loc<coord_type>(bas, nbas); \
    }

#define ADD_FUNCS_TO_KERNEL_3C(Kernel, optimizer_func, driver_func, integral_func, fill_func, coord_type) \
    void Kernel::optimizer(libcint::CINTOpt*& opt, \
        int* atm, int nat, int* bas, int nbas, double* env) { \
        optimizer_func(&opt, atm, nat, bas, nbas, env); \
    } \
    void Kernel::drv(double* out, int comp, int* shl_slice, int* aoloc, \
        libcint::CINTOpt* opt, int* atm, int nat, int* bas, int nbas, double* env) { \
        driver_func(integral_func, fill_func, out, comp, shl_slice, aoloc, opt, atm, nat, bas, nbas, env); \
    } \
    ivec Kernel::gen_loc(ivec& bas, int nbas) { \
        return make_loc<coord_type>(bas, nbas); \
    }

ADD_FUNCS_TO_KERNEL(Coulomb2C_SPH, int2c2e_optimizer, GTOint2c, int2c2e_sph, COORDINATE_TYPE::SPH)
ADD_FUNCS_TO_KERNEL(Coulomb2C_CRT, int2c2e_optimizer, GTOint2c, int2c2e_cart, COORDINATE_TYPE::CART)
ADD_FUNCS_TO_KERNEL(Overlap2C_SPH, libcint::int1e_ovlp_optimizer, GTOint2c, libcint::int1e_ovlp_sph, COORDINATE_TYPE::SPH) //If this is suddenly wrong, use the NOOPT version
ADD_FUNCS_TO_KERNEL(Overlap2C_CRT, libcint::int1e_ovlp_optimizer, GTOint2c, libcint::int1e_ovlp_cart, COORDINATE_TYPE::CART)

ADD_FUNCS_TO_KERNEL_3C(Coulomb3C_SPH, int3c2e_optimizer, GTOnr3c_drv, int3c2e_sph, GTOnr3c_fill_s1, COORDINATE_TYPE::SPH)
ADD_FUNCS_TO_KERNEL_3C(Coulomb3C_CRT, int3c2e_optimizer, GTOnr3c_drv, int3c2e_cart, GTOnr3c_fill_s1, COORDINATE_TYPE::CART)
ADD_FUNCS_TO_KERNEL_3C(Overlap3C_SPH, int3c1e_optimizer, GTOnr3c_drv, int3c1e_sph, GTOnr3c_fill_s1, COORDINATE_TYPE::SPH) //If this is suddenly wrong, use the NOOPT version

#undef ADD_FUNCS_TO_KERNEL
#undef ADD_FUNCS_TO_KERNEL_NOOPT
#undef ADD_FUNCS_TO_KERNEL_3C



// 2 slots of int param[]
#define POS_E1   0
#define TENSOR   1

#define LMAX            ANG_MAX
#define SIMDD           8
// 128s42p21d12f8g6h4i3j
#define NCTR_CART       128
#define NPRIMAX         40
#define BLKSIZE         56
#define EXPCUTOFF       50  // 1e-22
#define NOTZERO(e)      (fabs(e)>1e-18)

#define ALIGN8_UP(buf) (void *)(((uintptr_t)buf + 7) & (-(uintptr_t)8))

typedef int (*FPtr_exp)(double* ectr, double* coord, double* alpha, double* coeff,
    int l, int nprim, int nctr, size_t ngrids, double fac);
typedef void (*FPtr_eval)(double* gto, double* ri, double* exps,
    double* coord, double* alpha, double* coeff,
    double* env, int l, int np, int nc,
    size_t nao, size_t ngrids, size_t blksize);

/*
 * to optimize memory copy in cart2sph.c, remove the common factor for s
 * and p function in cart2sph
 */
double CINTcommon_fac_sp(FINT l)
{
    switch (l) {
    case 0: return 0.282094791773878143;
    case 1: return 0.488602511902919921;
    default: return 1;
    }
}


int GTOshloc_by_atom(int* shloc, int* shls_slice, int* ao_loc, int* atm, int* bas)
{
    const int sh0 = shls_slice[0];
    const int sh1 = shls_slice[1];
    int ish, nshblk, lastatm;
    shloc[0] = sh0;
    nshblk = 1;
    lastatm = bas[BAS_SLOTS * sh0 + ATOM_OF];
    for (ish = sh0; ish < sh1; ish++) {
        if (lastatm != bas[BAS_SLOTS * ish + ATOM_OF]) {
            lastatm = bas[BAS_SLOTS * ish + ATOM_OF];
            shloc[nshblk] = ish;
            nshblk++;
        }
    }
    shloc[nshblk] = sh1;
    return nshblk;
}

/*
 * non0table[ngrids/blksize,natm] is the T/F table for ao values to
 * screen the ao evaluation for each shell
 */
void GTOeval_loop(void (*fiter)(FPtr_eval, FPtr_exp, double, size_t, size_t, size_t, int*, int*, int*, double*,double*, double*, uint8_t*, int*, int, int*, int, double*),
    FPtr_eval feval, FPtr_exp fexp, double fac,
    int ngrids, int param[],int* shls_slice, int* ao_loc,
    double* ao, double* coord, uint8_t* non0table,
    int* atm, int natm, int* bas, int nbas, double* env)
{
    ivec shloc_vec(shls_slice[1] - shls_slice[0] + 1);
	int* shloc = shloc_vec.data();
    const int nshblk = GTOshloc_by_atom(shloc, shls_slice, ao_loc, atm, bas);
    const int nblk = (ngrids + BLKSIZE - 1) / BLKSIZE;
    const size_t Ngrids = ngrids;

#pragma omp parallel
    {
        const int sh0 = shls_slice[0];
        const int sh1 = shls_slice[1];
        const size_t nao = ao_loc[sh1] - ao_loc[sh0];
        int ip, ib, k, iloc, ish;
        size_t aoff, bgrids;
        int ncart = NCTR_CART * param[TENSOR] * param[POS_E1];
        double* buf = (double *)std::malloc(sizeof(double) * BLKSIZE * (NPRIMAX * 2 + ncart + 1));
#pragma omp for schedule(dynamic, 4)
        for (k = 0; k < nblk * nshblk; k++) {
            iloc = k / nblk;
            ish = shloc[iloc];
            aoff = ao_loc[ish] - ao_loc[sh0];
            ib = k - iloc * nblk;
            ip = ib * BLKSIZE;
            bgrids = std::min(ngrids - ip, BLKSIZE);
            (*fiter)(feval, fexp, fac, nao, Ngrids, bgrids,
                param, shloc + iloc, ao_loc, buf, ao + aoff * Ngrids + ip,
                coord + ip, non0table + ib * nbas,
                atm, natm, bas, nbas, env);
        }
        free(buf);
    }
}

int GTOcontract_exp0(double* ectr, double* coord, double* alpha, double* coeff,
    int l, int nprim, int nctr, size_t ngrids, double fac)
{
    size_t i, j, k;
    double arr, eprim;
    double rr[BLKSIZE];
    double* gridx = coord;
    double* gridy = coord + BLKSIZE;
    double* gridz = coord + BLKSIZE * 2;

#if defined(_MSC_VER)
#pragma loop(ivdep)
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC ivdep
#endif
    for (i = 0; i < ngrids; i++) {
        rr[i] = gridx[i] * gridx[i] + gridy[i] * gridy[i] + gridz[i] * gridz[i];
    }

    for (i = 0; i < nctr * BLKSIZE; i++) {
        ectr[i] = 0;
    }
    for (j = 0; j < nprim; j++) {
        for (i = 0; i < ngrids; i++) {
            arr = alpha[j] * rr[i];
            eprim = exp(-arr) * fac;
            for (k = 0; k < nctr; k++) {
                ectr[k * BLKSIZE + i] += eprim * coeff[k * nprim + j];
            }
        }
    }
    return 1;
}

// grid2atm[atm_id,xyz,grid_id]
static void _fill_grid2atm(double* grid2atm, double* coord, size_t bgrids, size_t ngrids,
    int* atm, int natm, int* bas, int nbas, double* env)
{
    int atm_id;
    size_t ig;
    double* r_atm;
    for (atm_id = 0; atm_id < natm; atm_id++) {
        r_atm = env + atm[PTR_COORD + atm_id * ATM_SLOTS];
#if defined(_MSC_VER)
#pragma loop(ivdep)
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC ivdep
#endif
        for (ig = 0; ig < bgrids; ig++) {
            grid2atm[0 * BLKSIZE + ig] = coord[0 * ngrids + ig] - r_atm[0];
            grid2atm[1 * BLKSIZE + ig] = coord[1 * ngrids + ig] - r_atm[1];
            grid2atm[2 * BLKSIZE + ig] = coord[2 * ngrids + ig] - r_atm[2];
        }
        grid2atm += 3 * BLKSIZE;
    }
}

static void _dset0(double* out, size_t odim, size_t bgrids, int counts)
{
    size_t i, j;
    for (i = 0; i < counts; i++) {
        for (j = 0; j < bgrids; j++) {
            out[i * odim + j] = 0;
        }
    }
}

void GTOeval_sph_iter(FPtr_eval feval, FPtr_exp fexp, double fac,
    size_t nao, size_t ngrids, size_t bgrids,
    int param[], int* shls_slice, int* ao_loc, double* buf,
    double* ao, double* coord, uint8_t* non0table,
    int* atm, int natm, int* bas, int nbas, double* env)
{
    const int ncomp = param[TENSOR];
    const int sh0 = shls_slice[0];
    const int sh1 = shls_slice[1];
    const int atmstart = bas[sh0 * BAS_SLOTS + ATOM_OF];
    const int atmend = bas[(sh1 - 1) * BAS_SLOTS + ATOM_OF] + 1;
    const int atmcount = atmend - atmstart;
    int i, k, l, np, nc, atm_id, bas_id, deg, dcart, ao_id;
    size_t di;
    double fac1;
    double* p_exp, * pcoeff, * pcoord, * pcart, * ri, * pao;
    double* grid2atm = (double*)(((uintptr_t)buf + 7) & (~(uintptr_t)8)); // [atm_id,xyz,grid]
    double* eprim = grid2atm + atmcount * 3 * BLKSIZE;
    double* cart_gto = eprim + NPRIMAX * BLKSIZE * 2;

    _fill_grid2atm(grid2atm, coord, bgrids, ngrids,
        atm + atmstart * ATM_SLOTS, atmcount, bas, nbas, env);

    for (bas_id = sh0; bas_id < sh1; bas_id++) {
        np = bas[bas_id * BAS_SLOTS + NPRIM_OF];
        nc = bas[bas_id * BAS_SLOTS + NCTR_OF];
        l = bas[bas_id * BAS_SLOTS + ANG_OF];
        deg = l * 2 + 1;
        fac1 = fac * CINTcommon_fac_sp(l);
        p_exp = env + bas[bas_id * BAS_SLOTS + PTR_EXP];
        pcoeff = env + bas[bas_id * BAS_SLOTS + PTR_COEFF];
        atm_id = bas[bas_id * BAS_SLOTS + ATOM_OF];
        pcoord = grid2atm + (atm_id - atmstart) * 3 * BLKSIZE;
        ao_id = ao_loc[bas_id] - ao_loc[sh0];
        if (non0table[bas_id] &&
            (*fexp)(eprim, pcoord, p_exp, pcoeff, l, np, nc, bgrids, fac1)) {
            dcart = (l + 1) * (l + 2) / 2;
            di = nc * dcart;
            ri = env + atm[PTR_COORD + atm_id * ATM_SLOTS];
            if (l <= 1) { // s, p functions
                (*feval)(ao + ao_id * ngrids, ri, eprim, pcoord, p_exp, pcoeff,
                    env, l, np, nc, nao, ngrids, bgrids);
            }
            else {
                (*feval)(cart_gto, ri, eprim, pcoord, p_exp, pcoeff,
                    env, l, np, nc, di, bgrids, bgrids);
                pcart = cart_gto;
                for (i = 0; i < ncomp; i++) {
                    pao = ao + (i * nao + ao_id) * ngrids;
                    for (k = 0; k < nc; k++) {
                        libcint::CINTc2s_ket_sph1(pao, pcart,
                            ngrids, bgrids, l);
                        pao += deg * ngrids;
                        pcart += dcart * bgrids;
                    }
                }
            }
        }
        else {
            for (i = 0; i < ncomp; i++) {
                _dset0(ao + (i * nao + ao_id) * ngrids, ngrids, bgrids, nc * deg);
            }
        }
    }
}

// pre-contracted grid AO evaluator
// contracted factors = \sum c_{i} exp(-a_i*r_i**2)
void GTOshell_eval_grid_cart(double* gto, double* ri, double* exps,
    double* coord, double* alpha, double* coeff,
    double* env, int l, int np, int nc,
    size_t nao, size_t ngrids, size_t blksize)
{
    int lx, ly, lz;
    size_t i, k;
    double buf[(LMAX + 1) * 3 * BLKSIZE + 8];
    double* xpows = (double*)(((uintptr_t)buf + 7) & (~(uintptr_t)8));
    double* ypows = xpows + (LMAX + 1) * BLKSIZE;
    double* zpows = ypows + (LMAX + 1) * BLKSIZE;
    double* gridx = coord;
    double* gridy = coord + BLKSIZE;
    double* gridz = coord + BLKSIZE * 2;

    switch (l) {
    case 0:
        for (k = 0; k < nc; k++) {
            for (i = 0; i < blksize; i++) {
                gto[k * ngrids + i] = exps[k * BLKSIZE + i];
            }
        }
        break;
    case 1:
        for (k = 0; k < nc; k++) {
#if defined(_MSC_VER)
#pragma loop(ivdep)
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC ivdep
#endif
            for (i = 0; i < blksize; i++) {
                gto[i] = gridx[i] * exps[k * BLKSIZE + i];
                gto[1 * ngrids + i] = gridy[i] * exps[k * BLKSIZE + i];
                gto[2 * ngrids + i] = gridz[i] * exps[k * BLKSIZE + i];
            }
            gto += ngrids * 3;
        }
        break;
    case 2:
        for (k = 0; k < nc; k++) {
#if defined(_MSC_VER)
#pragma loop(ivdep)
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC ivdep
#endif
            for (i = 0; i < blksize; i++) {
                gto[i] = exps[k * BLKSIZE + i] * gridx[i] * gridx[i]; // xx
                gto[1 * ngrids + i] = exps[k * BLKSIZE + i] * gridx[i] * gridy[i]; // xy
                gto[2 * ngrids + i] = exps[k * BLKSIZE + i] * gridx[i] * gridz[i]; // xz
                gto[3 * ngrids + i] = exps[k * BLKSIZE + i] * gridy[i] * gridy[i]; // yy
                gto[4 * ngrids + i] = exps[k * BLKSIZE + i] * gridy[i] * gridz[i]; // yz
                gto[5 * ngrids + i] = exps[k * BLKSIZE + i] * gridz[i] * gridz[i]; // zz
            }
            gto += ngrids * 6;
        }
        break;
    case 3:
        for (k = 0; k < nc; k++) {
#if defined(_MSC_VER)
#pragma loop(ivdep)
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC ivdep
#endif
            for (i = 0; i < blksize; i++) {
                gto[i] = exps[k * BLKSIZE + i] * gridx[i] * gridx[i] * gridx[i]; // xxx
                gto[1 * ngrids + i] = exps[k * BLKSIZE + i] * gridx[i] * gridx[i] * gridy[i]; // xxy
                gto[2 * ngrids + i] = exps[k * BLKSIZE + i] * gridx[i] * gridx[i] * gridz[i]; // xxz
                gto[3 * ngrids + i] = exps[k * BLKSIZE + i] * gridx[i] * gridy[i] * gridy[i]; // xyy
                gto[4 * ngrids + i] = exps[k * BLKSIZE + i] * gridx[i] * gridy[i] * gridz[i]; // xyz
                gto[5 * ngrids + i] = exps[k * BLKSIZE + i] * gridx[i] * gridz[i] * gridz[i]; // xzz
                gto[6 * ngrids + i] = exps[k * BLKSIZE + i] * gridy[i] * gridy[i] * gridy[i]; // yyy
                gto[7 * ngrids + i] = exps[k * BLKSIZE + i] * gridy[i] * gridy[i] * gridz[i]; // yyz
                gto[8 * ngrids + i] = exps[k * BLKSIZE + i] * gridy[i] * gridz[i] * gridz[i]; // yzz
                gto[9 * ngrids + i] = exps[k * BLKSIZE + i] * gridz[i] * gridz[i] * gridz[i]; // zzz
            }
            gto += ngrids * 10;
        }
        break;
    default:
        for (k = 0; k < nc; k++) {
            for (i = 0; i < blksize; i++) {
                xpows[i] = 1;
                ypows[i] = 1;
                zpows[i] = 1;
            }
            for (lx = 1; lx <= l; lx++) {
#if defined(_MSC_VER)
#pragma loop(ivdep)
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC ivdep
#endif
                for (i = 0; i < blksize; i++) {
                    xpows[lx * BLKSIZE + i] = xpows[(lx - 1) * BLKSIZE + i] * gridx[i];
                    ypows[lx * BLKSIZE + i] = ypows[(lx - 1) * BLKSIZE + i] * gridy[i];
                    zpows[lx * BLKSIZE + i] = zpows[(lx - 1) * BLKSIZE + i] * gridz[i];
                }
            }
            for (lx = l; lx >= 0; lx--) {
                for (ly = l - lx; ly >= 0; ly--) {
                    lz = l - lx - ly;
#if defined(_MSC_VER)
#pragma loop(ivdep)
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC ivdep
#endif
                    for (i = 0; i < blksize; i++) {
                        gto[i] = xpows[lx * BLKSIZE + i]
                            * ypows[ly * BLKSIZE + i]
                            * zpows[lz * BLKSIZE + i] * exps[k * BLKSIZE + i];
                    }
                    gto += ngrids;
                }
            }
        }
    }
}

void GTOeval_sph_drv(FPtr_eval feval, FPtr_exp fexp, double fac, int ngrids,
    int param[], int* shls_slice, int* ao_loc,
    double* ao, double* coord, uint8_t* non0table,
    int* atm, int natm, int* bas, int nbas, double* env)
{
    GTOeval_loop(GTOeval_sph_iter, feval, fexp, fac, ngrids,
        param, shls_slice, ao_loc,
        ao, coord, non0table, atm, natm, bas, nbas, env);
}

void GTOval_sph(int ngrids, int* shls_slice, int* ao_loc,
    double* ao, double* coord, uint8_t* non0table,
    int* atm, int natm, int* bas, int nbas, double* env)
{
    int param[] = { 1, 1 };
    GTOeval_sph_drv(GTOshell_eval_grid_cart, GTOcontract_exp0, 1,
        ngrids, param, shls_slice, ao_loc,
        ao, coord, non0table, atm, natm, bas, nbas, env);
}