#pragma once
#include "convenience.h"

#include "libCintKernels.h"

#include "cint_funcs.h"

#define BLKSIZE 8
/*
 * out[naoi,naoj,naok,comp] in F-order
 */
void GTOnr3c_fill_s1(
    int (*intor)(double*, int*, int*, int*, int, int*, int, double*, CINTOpt*, double*),
    double* out, double* buf, int comp, int jobid, int* shls_slice, int* ao_loc, CINTOpt* cintopt,
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



int GTOmax_shell_dim(const int* ao_loc, const int* shls_slice, int ncenter)
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
    int (*intor)(double*, int*, int*, int*, int, int*, int, double*, CINTOpt*, double*),
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
    int (*f)(double*, int*, int*, int*, int, int*, int, double*, CINTOpt*, double*) = (int (*)(double*, int*, int*, int*, int, int*, int, double*, CINTOpt*, double*))intor;
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
    int (*intor)(double*, int*, int*, int*, int, int*, int, double*, CINTOpt*, double*),
    void (*fill)(int (*intor)(double*, int*, int*, int*, int, int*, int, double*, CINTOpt*, double*), double*, double*, int, int, int*, int*, CINTOpt*, int*, int, int*, int, double*),
    double* eri, int comp, int* shls_slice, int* ao_loc, CINTOpt* cintopt,
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
        double* buf = (double*)malloc(sizeof(double) * (di * di * di * comp + cache_size));
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
void GTOint2c(int (*intor)(double*, int*, int*, int*, int, int*, int, double*, CINTOpt*, double*), double* mat, int comp, int hermi,
    int* shls_slice, int* ao_loc, CINTOpt* opt,
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
        double* cache = (double*)malloc(sizeof(double) * cache_size);
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
ivec make_loc(ivec& bas, int nbas) {
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
template ivec make_loc<COORDINATE_TYPE::CART>(ivec& bas, int nbas);
template ivec make_loc<COORDINATE_TYPE::SPH>(ivec& bas, int nbas);

extern "C" {
    extern CINTOptimizerFunction int3c2e_optimizer;
    extern CINTIntegralFunction int3c2e_sph;

    extern CINTOptimizerFunction int2c2e_optimizer;
    extern CINTIntegralFunction int2c2e_sph;

    extern CINTOptimizerFunction int3c1e_optimizer;
    extern CINTIntegralFunction int3c1e_sph;
}

#define ADD_FUNCS_TO_KERNEL(Kernel, optimizer_func, driver_func, integral_func, coord_type) \
    void Kernel::optimizer(CINTOpt*& opt, \
        int* atm, int nat, int* bas, int nbas, double* env) { \
        optimizer_func(&opt, atm, nat, bas, nbas, env); \
    } \
    void Kernel::drv(double* out, int comp, int* shl_slice, int* aoloc, \
        CINTOpt* opt, int* atm, int nat, int* bas, int nbas, double* env) { \
        driver_func(integral_func, out, comp, 0, shl_slice, aoloc, opt, atm, nat, bas, nbas, env); \
    } \
    ivec Kernel::gen_loc(ivec& bas, int nbas) { \
        return make_loc<coord_type>(bas, nbas); \
    }

#define ADD_FUNCS_TO_KERNEL_NOOPT(Kernel, driver_func, integral_func, coord_type) \
    void Kernel::optimizer(CINTOpt*& opt, \
        int* atm, int nat, int* bas, int nbas, double* env) { \
        opt = nullptr; \
    } \
    void Kernel::drv(double* out, int comp, int* shl_slice, int* aoloc, \
        CINTOpt* opt, int* atm, int nat, int* bas, int nbas, double* env) { \
        driver_func(integral_func, out, comp, 0, shl_slice, aoloc, opt, atm, nat, bas, nbas, env); \
    } \
    ivec Kernel::gen_loc(ivec& bas, int nbas) { \
        return make_loc<coord_type>(bas, nbas); \
    }

#define ADD_FUNCS_TO_KERNEL_3C(Kernel, optimizer_func, driver_func, integral_func, fill_func, coord_type) \
    void Kernel::optimizer(CINTOpt*& opt, \
        int* atm, int nat, int* bas, int nbas, double* env) { \
        optimizer_func(&opt, atm, nat, bas, nbas, env); \
    } \
    void Kernel::drv(double* out, int comp, int* shl_slice, int* aoloc, \
        CINTOpt* opt, int* atm, int nat, int* bas, int nbas, double* env) { \
        driver_func(integral_func, fill_func, out, comp, shl_slice, aoloc, opt, atm, nat, bas, nbas, env); \
    } \
    ivec Kernel::gen_loc(ivec& bas, int nbas) { \
        return make_loc<coord_type>(bas, nbas); \
    }

ADD_FUNCS_TO_KERNEL(Coulomb2C_SPH, int2c2e_optimizer, GTOint2c, int2c2e_sph, COORDINATE_TYPE::SPH)
ADD_FUNCS_TO_KERNEL(Overlap2C_SPH, int1e_ovlp_optimizer, GTOint2c, int1e_ovlp_sph, COORDINATE_TYPE::SPH) //If this is suddenly wrong, use the NOOPT version
ADD_FUNCS_TO_KERNEL(Overlap2C_CRT, int1e_ovlp_optimizer, GTOint2c, int1e_ovlp_cart, COORDINATE_TYPE::CART)
ADD_FUNCS_TO_KERNEL_3C(Coulomb3C_SPH, int3c2e_optimizer, GTOnr3c_drv, int3c2e_sph, GTOnr3c_fill_s1, COORDINATE_TYPE::SPH)
ADD_FUNCS_TO_KERNEL_3C(Overlap3C_SPH, int3c1e_optimizer, GTOnr3c_drv, int3c1e_sph, GTOnr3c_fill_s1, COORDINATE_TYPE::SPH) //If this is suddenly wrong, use the NOOPT version

#undef ADD_FUNCS_TO_KERNEL
#undef ADD_FUNCS_TO_KERNEL_NOOPT
#undef ADD_FUNCS_TO_KERNEL_3C