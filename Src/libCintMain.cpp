#include "libCintMain.h"

#include "cart2sph.h"

#include "int_optimizer.h"
#include "int_g2e.h"

#include "nos_math.h"

#if defined(__APPLE__)
// On macOS we are using Accelerate for BLAS/LAPACK
#include <Accelerate/Accelerate.h>
#else
// Linux/Windows with oneMKL
#include <mkl.h>
#endif

#define gctrg gout
#define gctrm gctr
#define mempty empty
#define m_ctr n_comp

#define NOVALUE ((void *)0xffffffffffffffffuL)
#define COMMON_ENVS_AND_DECLARE                                                                 \
    int *shls = envs->shls;                                                                     \
    int *bas = envs->bas;                                                                       \
    double *env = envs->env;                                                                    \
    int i_sh = shls[0];                                                                         \
    int j_sh = shls[1];                                                                         \
    CINTOpt *opt = envs->opt;                                                                   \
    if (opt->pairdata != NULL &&                                                                \
        opt->pairdata[i_sh * opt->nbas + j_sh] == NOVALUE)                                      \
    {                                                                                           \
        return 0;                                                                               \
    }                                                                                           \
    int k_sh = shls[2];                                                                         \
    int i_ctr = envs->x_ctr[0];                                                                 \
    int j_ctr = envs->x_ctr[1];                                                                 \
    int k_ctr = envs->x_ctr[2];                                                                 \
    int i_prim = bas(NPRIM_OF, i_sh);                                                           \
    int j_prim = bas(NPRIM_OF, j_sh);                                                           \
    int k_prim = bas(NPRIM_OF, k_sh);                                                           \
    double *ai = env + bas(PTR_EXP, i_sh);                                                      \
    double *aj = env + bas(PTR_EXP, j_sh);                                                      \
    double *ak = env + bas(PTR_EXP, k_sh);                                                      \
    double *ci = env + bas(PTR_COEFF, i_sh);                                                    \
    double *cj = env + bas(PTR_COEFF, j_sh);                                                    \
    double *ck = env + bas(PTR_COEFF, k_sh);                                                    \
    double expcutoff = envs->expcutoff;                                                         \
    double rr_ij = SQUARE(envs->rirj);                                                          \
    PairData *pdata_base, *pdata_ij;                                                            \
    if (opt->pairdata != NULL)                                                                  \
    {                                                                                           \
        pdata_base = opt->pairdata[i_sh * opt->nbas + j_sh];                                    \
    }                                                                                           \
    else                                                                                        \
    {                                                                                           \
        double *log_maxci = opt->log_max_coeff[i_sh];                                           \
        double *log_maxcj = opt->log_max_coeff[j_sh];                                           \
        MALLOC_INSTACK_PAIRDATA(pdata_base, i_prim *j_prim);                                    \
        if (CINTset_pairdata(pdata_base, ai, aj, envs->ri, envs->rj,                            \
                             log_maxci, log_maxcj, envs->li_ceil, envs->lj_ceil,                \
                             i_prim, j_prim, rr_ij, expcutoff, env))                            \
        {                                                                                       \
            return 0;                                                                           \
        }                                                                                       \
    }                                                                                           \
    int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;                                           \
    size_t nf = envs->nf;                                                                       \
    double fac1i, fac1j, fac1k;                                                                 \
    int ip, jp, kp;                                                                             \
    int _empty[4] = {1, 1, 1, 1};                                                               \
    int *iempty = _empty + 0;                                                                   \
    int *jempty = _empty + 1;                                                                   \
    int *kempty = _empty + 2;                                                                   \
    int *gempty = _empty + 3;                                                                   \
    int *non0ctri = opt->non0ctr[i_sh];                                                         \
    int *non0ctrj = opt->non0ctr[j_sh];                                                         \
    int *non0idxi = opt->sortedidx[i_sh];                                                       \
    int *non0idxj = opt->sortedidx[j_sh];                                                       \
    int *non0ctrk, *non0idxk;                                                                   \
    MALLOC_INSTACK_INT(non0ctrk, k_prim + k_prim * k_ctr);                                      \
    non0idxk = non0ctrk + k_prim;                                                               \
    CINTOpt_non0coeff_byshell(non0idxk, non0ctrk, ck, k_prim, k_ctr);                           \
    double expij, cutoff;                                                                       \
    double *rij;                                                                                \
    double *rkl = envs->rkl;                                                                    \
    int *idx = opt->index_xyz_array[envs->i_l * LMAX1 * LMAX1 + envs->j_l * LMAX1 + envs->k_l]; \
    if (idx == NULL)                                                                            \
    {                                                                                           \
        MALLOC_INSTACK_INT(idx, nf * 3);                                                        \
        CINTg2e_index_xyz(idx, envs);                                                           \
    }

#define ADJUST_CUTOFF                                                                  \
    double omega = env[PTR_RANGE_OMEGA];                                               \
    if (omega < 0 && envs->rys_order > 1)                                              \
    {                                                                                  \
        double r_guess = 8.;                                                           \
        double omega2 = omega * omega;                                                 \
        int lij = envs->li_ceil + envs->lj_ceil;                                       \
        if (lij > 0)                                                                   \
        {                                                                              \
            double dist_ij = sqrt(rr_ij);                                              \
            double aij = ai[i_prim - 1] + aj[j_prim - 1];                              \
            double theta = omega2 / (omega2 + aij);                                    \
            expcutoff += lij * log(                                                    \
                                   (dist_ij + theta * r_guess + 1.) / (dist_ij + 1.)); \
        }                                                                              \
        if (envs->lk_ceil > 0)                                                         \
        {                                                                              \
            double theta = omega2 / (omega2 + ak[k_prim - 1]);                         \
            expcutoff += envs->lk_ceil * log(theta * r_guess + 1.);                    \
        }                                                                              \
    }

#define SET_RIJ                      \
    if (pdata_ij->cceij > expcutoff) \
    {                                \
        goto i_contracted;           \
    }                                \
    envs->ai[0] = ai[ip];            \
    expij = pdata_ij->eij;           \
    rij = pdata_ij->rij;

#define PAIRDATA_NON0IDX_SIZE(ps)        \
    int *bas = envs->bas;                \
    int *shls = envs->shls;              \
    int i_prim = bas(NPRIM_OF, shls[0]); \
    int j_prim = bas(NPRIM_OF, shls[1]); \
    int k_prim = bas(NPRIM_OF, shls[2]); \
    int ps = (i_prim * j_prim * 5 + i_prim * x_ctr[0] + j_prim * x_ctr[1] + k_prim * x_ctr[2] + (i_prim + j_prim) * 2 + k_prim + envs->nf * 3 + 16);

#define SQUARE(r) ((r)[0] * (r)[0] + (r)[1] * (r)[1] + (r)[2] * (r)[2])

#define ALIAS_ADDR_IF_EQUAL(x, y) \
    if (y##_ctr == 1)             \
    {                             \
        gctr##x = gctr##y;        \
        x##empty = y##empty;      \
    }                             \
    else                          \
    {                             \
        gctr##x = g1;             \
        g1 += len##x;             \
    }

#define PRIM2CTR0(ctrsymb, gp, ngp)                                           \
    if (ctrsymb##_ctr > 1)                                                    \
    {                                                                         \
        if (*ctrsymb##empty)                                                  \
        {                                                                     \
            CINTprim_to_ctr_0(gctr##ctrsymb, gp, c##ctrsymb + ctrsymb##p,     \
                              ngp, ctrsymb##_prim, ctrsymb##_ctr,             \
                              non0ctr##ctrsymb[ctrsymb##p],                   \
                              non0idx##ctrsymb + ctrsymb##p * ctrsymb##_ctr); \
        }                                                                     \
        else                                                                  \
        {                                                                     \
            CINTprim_to_ctr_1(gctr##ctrsymb, gp, c##ctrsymb + ctrsymb##p,     \
                              ngp, ctrsymb##_prim, ctrsymb##_ctr,             \
                              non0ctr##ctrsymb[ctrsymb##p],                   \
                              non0idx##ctrsymb + ctrsymb##p * ctrsymb##_ctr); \
        }                                                                     \
    }                                                                         \
    *ctrsymb##empty = 0

#define TRANSPOSE(a)                                  \
    if (*empty)                                       \
    {                                                 \
        CINTdmat_transpose(gctr, a, nf *nc, n_comp);  \
    }                                                 \
    else                                              \
    {                                                 \
        CINTdplus_transpose(gctr, a, nf *nc, n_comp); \
    }                                                 \
    *empty = 0;

void CINTprim_to_ctr_0(double *gc, double *gp, double *coeff, size_t nf,
                       int nprim, int nctr, int non0ctr, int *sortedidx)
{
    int i;
    size_t n;
    double c0;

    for (i = 0; i < nctr; i++)
    {
        c0 = coeff[nprim * i];
        for (n = 0; n < nf; n++)
        {
            gc[nf * i + n] = c0 * gp[n];
        }
    }
}

void CINTprim_to_ctr_1(double *gc, double *gp, double *coeff, size_t nf,
                       int nprim, int nctr, int non0ctr, int *sortedidx)
{
    int i, j;
    size_t n;
    double c0;

    for (i = 0; i < non0ctr; i++)
    {
        c0 = coeff[nprim * sortedidx[i]];
        j = sortedidx[i];
        for (n = 0; n < nf; n++)
        {
            gc[nf * j + n] += c0 * gp[n];
        }
    }
}

/*
 * a[m,n] -> a_t[n,m]
 */
void CINTdmat_transpose(double *a_t, double *a, int m, int n)
{
    int i, j, k;

    for (j = 0; j < n - 3; j += 4)
    {
#pragma ivdep
        for (i = 0; i < m; i++)
        {
            a_t[(j + 0) * m + i] = a[i * n + j + 0];
            a_t[(j + 1) * m + i] = a[i * n + j + 1];
            a_t[(j + 2) * m + i] = a[i * n + j + 2];
            a_t[(j + 3) * m + i] = a[i * n + j + 3];
        }
    }

    switch (n - j)
    {
    case 1:
#pragma ivdep
        for (i = 0; i < m; i++)
        {
            a_t[j * m + i] = a[i * n + j];
        }
        break;
    case 2:
#pragma ivdep
        for (i = 0; i < m; i++)
        {
            a_t[(j + 0) * m + i] = a[i * n + j + 0];
            a_t[(j + 1) * m + i] = a[i * n + j + 1];
        }
        break;
    case 3:
#pragma ivdep
        for (i = 0; i < m; i++)
        {
            a_t[(j + 0) * m + i] = a[i * n + j + 0];
            a_t[(j + 1) * m + i] = a[i * n + j + 1];
            a_t[(j + 2) * m + i] = a[i * n + j + 2];
        }
        break;
    }
}

/*
 * a_t[n,m] += a[m,n]
 */
void CINTdplus_transpose(double *a_t, double *a, int m, int n)
{
    int i, j, k;

    for (j = 0; j < n - 3; j += 4)
    {
#pragma ivdep
        for (i = 0; i < m; i++)
        {
            a_t[(j + 0) * m + i] += a[i * n + j + 0];
            a_t[(j + 1) * m + i] += a[i * n + j + 1];
            a_t[(j + 2) * m + i] += a[i * n + j + 2];
            a_t[(j + 3) * m + i] += a[i * n + j + 3];
        }
    }

    switch (n - j)
    {
    case 1:
#pragma ivdep
        for (i = 0; i < m; i++)
        {
            a_t[j * m + i] += a[i * n + j];
        }
        break;
    case 2:
#pragma ivdep
        for (i = 0; i < m; i++)
        {
            a_t[(j + 0) * m + i] += a[i * n + j + 0];
            a_t[(j + 1) * m + i] += a[i * n + j + 1];
        }
        break;
    case 3:
#pragma ivdep
        for (i = 0; i < m; i++)
        {
            a_t[(j + 0) * m + i] += a[i * n + j + 0];
            a_t[(j + 1) * m + i] += a[i * n + j + 1];
            a_t[(j + 2) * m + i] += a[i * n + j + 2];
        }
        break;
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
    int (*intor)(double *, int *, int *, int *, int, int *, int, double *, CINTOpt *, double *),
    int *shls_slice, int ncenter,
    int *atm, int natm, int *bas, int nbas, double *env)
{
    int i;
    int i0 = shls_slice[0];
    int i1 = shls_slice[1];
    for (i = 1; i < ncenter; i++)
    {
        i0 = std::min(i0, shls_slice[i * 2]);
        i1 = std::max(i1, shls_slice[i * 2 + 1]);
    }
    int (*f)(double *, int *, int *, int *, int, int *, int, double *, CINTOpt *, double *) = (int (*)(double *, int *, int *, int *, int, int *, int, double *, CINTOpt *, double *))intor;
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
    int (*intor)(double *, int *, int *, int *, int, int *, int, double *, CINTOpt *, double *),
    double *out, int comp, int *shls_slice, int *ao_loc, CINTOpt *cintopt,
    int *atm, int natm, int *bas, int nbas, double *env)
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
    const int naoi = ao_loc[ish1] - ao_loc[ish0];
    const int naoj = ao_loc[jsh1] - ao_loc[jsh0];
    const int naok = ao_loc[ksh1] - ao_loc[ksh0];

    const int di = GTOmax_shell_dim(ao_loc, shls_slice, 3);
    const int cache_size = GTOmax_cache_size(intor, shls_slice, 3,
                                             atm, natm, bas, nbas, env);
    const int njobs = (std::max(nish, njsh) / BLKSIZE + 1) * nksh;

    // #pragma omp parallel
    {
        int dims[3] = {naoi, naoj, naok};
        int shls[3] = {0, 0, 0};
        int ish, jsh, ksh, i0, j0, k0;
        double *cache = (double *)malloc(sizeof(double) * cache_size * di * di * di);
        // #pragma omp for nowait schedule(dynamic)
        for (ksh = ksh0; ksh < ksh1; ksh++)
        {
            for (jsh = jsh0; jsh < jsh1; jsh++)
            {
                for (ish = ish0; ish < ish1; ish++)
                {
                    shls[0] = ish;
                    shls[1] = jsh;
                    shls[2] = ksh;
                    i0 = ao_loc[ish] - ao_loc[ish0];
                    j0 = ao_loc[jsh] - ao_loc[jsh0];
                    k0 = ao_loc[ksh] - ao_loc[ksh0];
                    int offset = k0 * naoi * naoj + j0 * naoi + i0;
                    (*intor)(out + offset, dims, shls, atm, natm, bas, nbas, env,
                             cintopt, cache);
                }
            }
        }
        free(cache);
    }
}

void GTOnr3c_drv(
    int (*intor)(double *, int *, int *, int *, int, int *, int, double *, CINTOpt *, double *),
    void (*fill)(int (*intor)(double *, int *, int *, int *, int, int *, int, double *, CINTOpt *, double *), double *, double *, int, int, int *, int *, CINTOpt *, int *, int, int *, int, double *),
    double *eri, int comp, int *shls_slice, int *ao_loc, CINTOpt *cintopt,
    int *atm, int natm, int *bas, int nbas, double *env)
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
 * out[naoi,naoj,naok,comp] in F-order
 */
void GTOnr3c_fill_s1(
    int (*intor)(double *, int *, int *, int *, int, int *, int, double *, CINTOpt *, double *),
    double *out, double *buf, int comp, int jobid, int *shls_slice, int *ao_loc, CINTOpt *cintopt,
    int *atm, int natm, int *bas, int nbas, double *env)
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
    int dims[] = {(int)naoi, (int)naoj, (int)naok};

    const int k0 = ao_loc[ksh] - ao_loc[ksh0];
    out += naoi * naoj * k0;

    int ish, jsh, i0, j0;
    int shls[3] = {0, 0, ksh};

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

void CINTgout2e(double *gout,
                double *g,
                int *idx,
                CINTEnvVars *envs,
                int gout_empty)
{
    const int nf = envs->nf;
    const int nrys_roots = envs->nrys_roots;
    double s = 0.0;
    if (gout_empty)
    {
        // If gout is to be overwritten (empty), do a direct assignment
        for (int n = 0; n < nf; n++, idx += 3)
        {
            // Pointers to the relevant slices of g
            const double *gx = g + idx[0];
            const double *gy = g + idx[1];
            const double *gz = g + idx[2];

            // Inner loop: accumulate products
#pragma ivdep
            for (int i = 0; i < nrys_roots; i++)
            {
                s += gx[i] * gy[i] * gz[i];
            }
            gout[n] = s;
            s = 0.0;
        }
    }
    else
    {
        // If gout already has data, accumulate
        for (int n = 0; n < nf; n++, idx += 3)
        {
            const double *gx = g + idx[0];
            const double *gy = g + idx[1];
            const double *gz = g + idx[2];

#pragma ivdep
            for (int i = 0; i < nrys_roots; i++)
            {
                s += gx[i] * gy[i] * gz[i];
            }
            gout[n] += s;
            s = 0.0;
        }
    }
}

int CINT2c2e_loop(double *gctr, CINTEnvVars *envs, double *cache, int *empty)
{
    int *shls = envs->shls;
    int *bas = envs->bas;
    double *env = envs->env;
    int i_sh = shls[0];
    int k_sh = shls[1];
    int i_ctr = envs->x_ctr[0];
    int k_ctr = envs->x_ctr[1];
    int i_prim = bas(NPRIM_OF, i_sh);
    int k_prim = bas(NPRIM_OF, k_sh);
    double *ai = env + bas(PTR_EXP, i_sh);
    double *ak = env + bas(PTR_EXP, k_sh);
    double *ci = env + bas(PTR_COEFF, i_sh);
    double *ck = env + bas(PTR_COEFF, k_sh);
    double expcutoff = envs->expcutoff;
    double *ri = envs->ri;
    double *rk = envs->rk;
    int n_comp = envs->ncomp_tensor;
    double fac1i, fac1k;
    int ip, kp;
    int _empty[3] = {1, 1, 1};
    int *iempty = _empty + 0;
    int *kempty = _empty + 1;
    int *gempty = _empty + 2;
    int *non0ctri, *non0ctrk;
    int *non0idxi, *non0idxk;
    MALLOC_INSTACK_INT(non0ctri, i_prim + k_prim + i_prim * i_ctr + k_prim * k_ctr);
    non0ctrk = non0ctri + i_prim;
    non0idxi = non0ctrk + k_prim;
    non0idxk = non0idxi + i_prim * i_ctr;
    if (i_ctr > 1)
    {
        CINTOpt_non0coeff_byshell(non0idxi, non0ctri, ci, i_prim, i_ctr);
    }
    if (k_ctr > 1)
    {
        CINTOpt_non0coeff_byshell(non0idxk, non0ctrk, ck, k_prim, k_ctr);
    }

    int *idx = envs->opt->index_xyz_array[envs->i_l * LMAX1 + envs->k_l];

    size_t nf = envs->nf;
    const int nc = i_ctr * k_ctr;
    const int leng = envs->g_size * 3 * ((1 << envs->gbits) + 1);
    const int lenk = nf * nc * n_comp;    // gctrk
    const int leni = nf * i_ctr * n_comp; // gctri
    const int len0 = nf * n_comp;         // gout
    const int len = leng + lenk + leni + len0;
    double *g;
    MALLOC_INSTACK(g, len);
    double *g1 = g + leng;
    double *gout, *gctri, *gctrk;

    ALIAS_ADDR_IF_EQUAL(k, m);
    ALIAS_ADDR_IF_EQUAL(i, k);
    ALIAS_ADDR_IF_EQUAL(g, i);

    for (kp = 0; kp < k_prim; kp++)
    {
        envs->ak[0] = ak[kp];
        if (k_ctr == 1)
        {
            fac1k = envs->common_factor * ck[kp];
        }
        else
        {
            fac1k = envs->common_factor;
            *iempty = 1;
        }
        for (ip = 0; ip < i_prim; ip++)
        {
            envs->ai[0] = ai[ip];
            if (i_ctr == 1)
            {
                fac1i = fac1k * ci[ip];
            }
            else
            {
                fac1i = fac1k;
            }
            envs->fac[0] = fac1i;
            if ((*envs->f_g0_2e)(g, ri, rk, expcutoff, envs))
            {
                (*envs->f_gout)(gout, g, idx, envs, *gempty);
                PRIM2CTR0(i, gout, len0);
            }
        } // end loop i_prim
        if (!*iempty)
        {
            PRIM2CTR0(k, gctri, leni);
        }
    } // end loop k_prim

    if (n_comp > 1 && !*kempty)
    {
        TRANSPOSE(gctrk);
    }
    return !*empty;
}

int CINT2c2e_loop_nopt(double *gctr, CINTEnvVars *envs, double *cache, int *empty)
{
    int *shls = envs->shls;
    int *bas = envs->bas;
    double *env = envs->env;
    int i_sh = shls[0];
    int k_sh = shls[1];
    int i_ctr = envs->x_ctr[0];
    int k_ctr = envs->x_ctr[1];
    int i_prim = bas(NPRIM_OF, i_sh);
    int k_prim = bas(NPRIM_OF, k_sh);
    double *ai = env + bas(PTR_EXP, i_sh);
    double *ak = env + bas(PTR_EXP, k_sh);
    double *ci = env + bas(PTR_COEFF, i_sh);
    double *ck = env + bas(PTR_COEFF, k_sh);
    double expcutoff = envs->expcutoff;
    double *ri = envs->ri;
    double *rk = envs->rk;
    int n_comp = envs->ncomp_tensor;
    double fac1i, fac1k;
    int ip, kp;
    int _empty[3] = {1, 1, 1};
    int *iempty = _empty + 0;
    int *kempty = _empty + 1;
    int *gempty = _empty + 2;
    size_t nf = envs->nf;
    const int nc = i_ctr * k_ctr;
    const int leng = envs->g_size * 3 * ((1 << envs->gbits) + 1);
    const int lenk = nf * nc * n_comp;    // gctrk
    const int leni = nf * i_ctr * n_comp; // gctri
    const int len0 = nf * n_comp;         // gout
    const int len = leng + lenk + leni + len0;
    double *g;
    MALLOC_INSTACK(g, len);
    double *g1 = g + leng;
    double *gout, *gctri, *gctrk;

    ALIAS_ADDR_IF_EQUAL(k, m);
    ALIAS_ADDR_IF_EQUAL(i, k);
    ALIAS_ADDR_IF_EQUAL(g, i);

    int *idx;
    MALLOC_INSTACK_INT(idx, envs->nf * 3);
    CINTg1e_index_xyz(idx, envs);

    int *non0ctri, *non0ctrk;
    int *non0idxi, *non0idxk;
    MALLOC_INSTACK_INT(non0ctri, i_prim + k_prim + i_prim * i_ctr + k_prim * k_ctr);
    non0ctrk = non0ctri + i_prim;
    non0idxi = non0ctrk + k_prim;
    non0idxk = non0idxi + i_prim * i_ctr;
    if (i_ctr > 1)
    {
        CINTOpt_non0coeff_byshell(non0idxi, non0ctri, ci, i_prim, i_ctr);
    }
    if (k_ctr > 1)
    {
        CINTOpt_non0coeff_byshell(non0idxk, non0ctrk, ck, k_prim, k_ctr);
    }

    for (kp = 0; kp < k_prim; kp++)
    {
        envs->ak[0] = ak[kp];
        envs->al[0] = 0; // to use CINTg0_2e
        if (k_ctr == 1)
        {
            fac1k = envs->common_factor * ck[kp];
        }
        else
        {
            fac1k = envs->common_factor;
            *iempty = 1;
        }
        for (ip = 0; ip < i_prim; ip++)
        {
            envs->ai[0] = ai[ip];
            envs->aj[0] = 0;
            if (i_ctr == 1)
            {
                fac1i = fac1k * ci[ip];
            }
            else
            {
                fac1i = fac1k;
            }
            envs->fac[0] = fac1i;
            if ((*envs->f_g0_2e)(g, ri, rk, expcutoff, envs))
            {
                (*envs->f_gout)(gout, g, idx, envs, *gempty);
                PRIM2CTR0(i, gout, len0);
            }
        } // end loop i_prim
        if (!*iempty)
        {
            PRIM2CTR0(k, gctri, leni);
        }
    } // end loop k_prim

    if (n_comp > 1 && !*kempty)
    {
        TRANSPOSE(gctrk);
    }
    return !*empty;
}

int int1e_cache_size(CINTEnvVars *envs)
{
    int *shls = envs->shls;
    int *bas = envs->bas;
    int i_prim = bas(NPRIM_OF, shls[0]);
    int j_prim = bas(NPRIM_OF, shls[1]);
    int *x_ctr = envs->x_ctr;
    int nc = envs->nf * x_ctr[0] * x_ctr[1];
    int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
    int leng = envs->g_size * 3 * ((1 << envs->gbits) + 1);
    int lenj = envs->nf * nc * n_comp;
    int leni = envs->nf * x_ctr[0] * n_comp;
    int len0 = envs->nf * n_comp;
    int pdata_size = (i_prim * j_prim * 5 + i_prim * x_ctr[0] + j_prim * x_ctr[1] + (i_prim + j_prim) * 2 + envs->nf * 3);
    int cache_size = std::max(nc * n_comp + leng + lenj + leni + len0 + pdata_size,
                              nc * n_comp + envs->nf * 8 * OF_CMPLX);
    return cache_size;
}

int CINT2c2e_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
                 double *cache, void (*f_c2s)(double *, double *, int *, CINTEnvVars *, double *))
{
    int *x_ctr = envs->x_ctr;
    int nc = envs->nf * x_ctr[0] * x_ctr[1];
    int n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
    if (out == NULL)
    {
        return int1e_cache_size(envs);
    }
    double *stack = NULL;
    if (cache == NULL)
    {
        size_t cache_size = int1e_cache_size(envs);
        stack = (double *)malloc(sizeof(double) * cache_size);
        cache = stack;
    }
    double *gctr;
    MALLOC_INSTACK(gctr, nc * n_comp);

    int n;
    int empty = 1;
    if (opt != NULL)
    {
        envs->opt = opt;
        CINT2c2e_loop(gctr, envs, cache, &empty);
    }
    else
    {
        CINT2c2e_loop_nopt(gctr, envs, cache, &empty);
    }

    int counts[4];
    if (f_c2s == &c2s_sph_1e)
    {
        counts[0] = (envs->i_l * 2 + 1) * x_ctr[0];
        counts[1] = (envs->k_l * 2 + 1) * x_ctr[1];
    }
    else
    {
        counts[0] = envs->nfi * x_ctr[0];
        counts[1] = envs->nfk * x_ctr[1];
    }
    counts[2] = 1;
    counts[3] = 1;
    if (dims == NULL)
    {
        dims = counts;
    }
    int nout = dims[0] * dims[1];
    if (!empty)
    {
        for (n = 0; n < n_comp; n++)
        {
            (*f_c2s)(out + nout * n, gctr + nc * n, dims, envs, cache);
        }
    }
    else
    {
        for (n = 0; n < n_comp; n++)
        {
            c2s_dset0(out + nout * n, dims, counts);
        }
    }
    if (stack != NULL)
    {
        free(stack);
    }
    return !empty;
}

int int2c2e_sph(double *out, int *dims, int *shls, int *atm, int natm,
                int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
    int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
    CINTEnvVars envs;
    CINTinit_int2c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
    envs.f_gout = &CINTgout2e;
    return CINT2c2e_drv(out, dims, &envs, opt, cache, &c2s_sph_1e);
}

/*
 * mat(naoi,naoj,comp) in F-order
 */
void GTOint2c(int (*intor)(double *, int *, int *, int *, int, int *, int, double *, CINTOpt *, double *), double *mat, int comp, int hermi,
              int *shls_slice, int *ao_loc, CINTOpt *opt,
              int *atm, int natm, int *bas, int nbas, double *env)
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
        int dims[] = {(int)naoi, (int)naoj};
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

void computeEri2c(Int_Params &params, vec &eri2c)
{
    ivec bas = params.get_bas();
    ivec atm = params.get_atm();
    vec env = params.get_env();

    int nbas = params.get_nbas();
    int nat = params.get_natoms();

    ivec shl_slice = {0, nbas, 0, nbas};
    ivec aoloc = make_loc<COORDINATE_TYPE::SPH>(bas, nbas);

    int naoi = aoloc[shl_slice[1]] - aoloc[shl_slice[0]];
    int naoj = aoloc[shl_slice[3]] - aoloc[shl_slice[2]];

    CINTOpt *opty = nullptr;
    int2c2e_optimizer(&opty, atm.data(), nat, bas.data(), nbas, env.data());

    // Compute integrals
    vec res((size_t)naoi * (size_t)naoj, 0.0);
    eri2c.resize((size_t)naoi * (size_t)naoj, 0.0);
    GTOint2c(int2c2e_sph, res.data(), 1, 0, shl_slice.data(), aoloc.data(), opty, atm.data(), nat, bas.data(), nbas, env.data());

    // res is in fortran order, write the result in regular ordering
    for (int i = 0; i < naoi; i++)
    {
        for (int j = 0; j < naoj; j++)
        {
            eri2c[(size_t)j * (size_t)naoi + i] = res[(size_t)i * (size_t)naoj + j];
        }
    }
}

int CINT3c2e_loop_nopt(double *gctr, CINTEnvVars *envs, double *cache, int *empty)
{
    int *shls = envs->shls;
    int *bas = envs->bas;
    double *env = envs->env;
    int i_sh = shls[0];
    int j_sh = shls[1];
    int k_sh = shls[2];
    int i_ctr = envs->x_ctr[0];
    int j_ctr = envs->x_ctr[1];
    int k_ctr = envs->x_ctr[2];
    int i_prim = bas(NPRIM_OF, i_sh);
    int j_prim = bas(NPRIM_OF, j_sh);
    int k_prim = bas(NPRIM_OF, k_sh);
    // double *ri = envs->ri;
    // double *rj = envs->rj;
    double *ai = env + bas(PTR_EXP, i_sh);
    double *aj = env + bas(PTR_EXP, j_sh);
    double *ak = env + bas(PTR_EXP, k_sh);
    double *ci = env + bas(PTR_COEFF, i_sh);
    double *cj = env + bas(PTR_COEFF, j_sh);
    double *ck = env + bas(PTR_COEFF, k_sh);

    double expcutoff = envs->expcutoff;
    double rr_ij = SQUARE(envs->rirj);
    double *log_maxci, *log_maxcj;
    PairData *pdata_base, *pdata_ij;
    MALLOC_INSTACK(log_maxci, i_prim + j_prim);
    MALLOC_INSTACK_PAIRDATA(pdata_base, i_prim * j_prim);
    log_maxcj = log_maxci + i_prim;
    CINTOpt_log_max_pgto_coeff(log_maxci, ci, i_prim, i_ctr);
    CINTOpt_log_max_pgto_coeff(log_maxcj, cj, j_prim, j_ctr);
    if (CINTset_pairdata(pdata_base, ai, aj, envs->ri, envs->rj,
                         log_maxci, log_maxcj, envs->li_ceil, envs->lj_ceil,
                         i_prim, j_prim, rr_ij, expcutoff, env))
    {
        return 0;
    }

    int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
    size_t nf = envs->nf;
    double fac1i, fac1j, fac1k;
    int ip, jp, kp;
    int _empty[4] = {1, 1, 1, 1};
    int *iempty = _empty + 0;
    int *jempty = _empty + 1;
    int *kempty = _empty + 2;
    int *gempty = _empty + 3;
    /* COMMON_ENVS_AND_DECLARE end */

    double expij, cutoff;
    double *rij;
    double *rkl = envs->rk;
    double omega = env[PTR_RANGE_OMEGA];
    if (omega < 0 && envs->rys_order > 1)
    {
        double r_guess = 8.;
        double omega2 = omega * omega;
        int lij = envs->li_ceil + envs->lj_ceil;
        if (lij > 0)
        {
            double dist_ij = sqrt(rr_ij);
            double aij = ai[i_prim - 1] + aj[j_prim - 1];
            double theta = omega2 / (omega2 + aij);
            expcutoff += lij * std::log(
                                   (dist_ij + theta * r_guess + 1.) / (dist_ij + 1.));
        }
        if (envs->lk_ceil > 0)
        {
            double theta = omega2 / (omega2 + ak[k_prim - 1]);
            expcutoff += envs->lk_ceil * std::log(theta * r_guess + 1.);
        }
    }

    int *idx;
    MALLOC_INSTACK_INT(idx, nf * 3);
    CINTg2e_index_xyz(idx, envs);

    int *non0ctri, *non0ctrj, *non0ctrk;
    int *non0idxi, *non0idxj, *non0idxk;
    MALLOC_INSTACK_INT(non0ctri, i_prim + j_prim + k_prim + i_prim * i_ctr + j_prim * j_ctr + k_prim * k_ctr);
    non0ctrj = non0ctri + i_prim;
    non0ctrk = non0ctrj + j_prim;
    non0idxi = non0ctrk + k_prim;
    non0idxj = non0idxi + i_prim * i_ctr;
    non0idxk = non0idxj + j_prim * j_ctr;
    CINTOpt_non0coeff_byshell(non0idxi, non0ctri, ci, i_prim, i_ctr);
    CINTOpt_non0coeff_byshell(non0idxj, non0ctrj, cj, j_prim, j_ctr);
    CINTOpt_non0coeff_byshell(non0idxk, non0ctrk, ck, k_prim, k_ctr);

    int nc = i_ctr * j_ctr * k_ctr;
    // (irys,i,j,k,l,coord,0:1); +1 for nabla-r12
    size_t leng = (size_t)envs->g_size * 3 * ((size_t)(1 << envs->gbits) + 1);
    size_t lenk = nf * nc * n_comp;            // gctrk
    size_t lenj = nf * i_ctr * j_ctr * n_comp; // gctrj
    size_t leni = nf * i_ctr * n_comp;         // gctri
    size_t len0 = nf * n_comp;                 // gout
    size_t len = leng + lenk + lenj + leni + len0;
    double *g;
    MALLOC_INSTACK(g, len); // must be allocated last in this function
    double *g1 = g + leng;
    double *gout, *gctri, *gctrj, *gctrk;

    ALIAS_ADDR_IF_EQUAL(k, m);
    ALIAS_ADDR_IF_EQUAL(j, k);
    ALIAS_ADDR_IF_EQUAL(i, j);
    ALIAS_ADDR_IF_EQUAL(g, i);

    for (kp = 0; kp < k_prim; kp++)
    {
        envs->ak[0] = ak[kp];
        if (k_ctr == 1)
        {
            fac1k = envs->common_factor * ck[kp];
        }
        else
        {
            fac1k = envs->common_factor;
            *jempty = 1;
        }

        pdata_ij = pdata_base;
        for (jp = 0; jp < j_prim; jp++)
        {
            envs->aj[0] = aj[jp];
            if (j_ctr == 1)
            {
                fac1j = fac1k * cj[jp];
            }
            else
            {
                fac1j = fac1k;
                *iempty = 1;
            }
            for (ip = 0; ip < i_prim; ip++, pdata_ij++)
            {
                if (pdata_ij->cceij > expcutoff)
                {
                    goto i_contracted;
                }
                envs->ai[0] = ai[ip];
                expij = pdata_ij->eij;
                rij = pdata_ij->rij;
                cutoff = expcutoff - pdata_ij->cceij;
                if (i_ctr == 1)
                {
                    fac1i = fac1j * ci[ip] * expij;
                }
                else
                {
                    fac1i = fac1j * expij;
                }
                envs->fac[0] = fac1i;
                if ((*envs->f_g0_2e)(g, rij, rkl, cutoff, envs))
                {
                    (*envs->f_gout)(gout, g, idx, envs, *gempty);
                    PRIM2CTR0(i, gout, len0);
                }
            i_contracted:;
            } // end loop i_prim
            if (!*iempty)
            {
                PRIM2CTR0(j, gctri, leni);
            }
        } // end loop j_prim
        if (!*jempty)
        {
            PRIM2CTR0(k, gctrj, lenj);
        }
    } // end loop k_prim

    if (n_comp > 1 && !*kempty)
    {
        TRANSPOSE(gctrk);
    }
    return !*empty;
}

int CINT3c2e_loop(double *gctr, CINTEnvVars *envs, double *cache, int *empty)
{
    COMMON_ENVS_AND_DECLARE
    ADJUST_CUTOFF;
    int nc = i_ctr * j_ctr * k_ctr;
    // (irys,i,j,k,coord,0:1); +1 for nabla-r12
    size_t leng = (size_t)envs->g_size * 3 * (size_t)((1 << envs->gbits) + 1);
    size_t lenk = nf * nc * n_comp;            // gctrk
    size_t lenj = nf * i_ctr * j_ctr * n_comp; // gctrj
    size_t leni = nf * i_ctr * n_comp;         // gctri
    size_t len0 = nf * n_comp;                 // gout
    size_t len = leng + lenk + lenj + leni + len0;
    double *g;
    MALLOC_INSTACK(g, len);
    double *g1 = g + leng;
    double *gout, *gctri, *gctrj, *gctrk;

    ALIAS_ADDR_IF_EQUAL(k, m);
    ALIAS_ADDR_IF_EQUAL(j, k);
    ALIAS_ADDR_IF_EQUAL(i, j);
    ALIAS_ADDR_IF_EQUAL(g, i);

    for (kp = 0; kp < k_prim; kp++)
    {
        envs->ak[0] = ak[kp];
        if (k_ctr == 1)
        {
            fac1k = envs->common_factor * ck[kp];
        }
        else
        {
            fac1k = envs->common_factor;
            *jempty = 1;
        }
        pdata_ij = pdata_base;
        for (jp = 0; jp < j_prim; jp++)
        {
            envs->aj[0] = aj[jp];
            if (j_ctr == 1)
            {
                fac1j = fac1k * cj[jp];
            }
            else
            {
                fac1j = fac1k;
                *iempty = 1;
            }
            for (ip = 0; ip < i_prim; ip++, pdata_ij++)
            {
                SET_RIJ;
                cutoff = expcutoff - pdata_ij->cceij;
                if (i_ctr == 1)
                {
                    fac1i = fac1j * ci[ip] * expij;
                }
                else
                {
                    fac1i = fac1j * expij;
                }
                envs->fac[0] = fac1i;
                if ((*envs->f_g0_2e)(g, rij, rkl, cutoff, envs))
                {
                    (*envs->f_gout)(gout, g, idx, envs, *gempty);
                    PRIM2CTR0(i, gout, len0);
                }
            i_contracted:;
            } // end loop i_prim
            if (!*iempty)
            {
                PRIM2CTR0(j, gctri, leni);
            }
        } // end loop j_prim
        if (!*jempty)
        {
            PRIM2CTR0(k, gctrj, lenj);
        }
    } // end loop k_prim

    if (n_comp > 1 && !*kempty)
    {
        TRANSPOSE(gctrk);
    }
    return !*empty;
}

// i_ctr = n; j_ctr = k_ctr = 1;
int CINT3c2e_n11_loop(double *gctr, CINTEnvVars *envs, double *cache, int *empty)
{
    COMMON_ENVS_AND_DECLARE;
    ADJUST_CUTOFF;
    int nc = i_ctr;
    size_t leng = (size_t)envs->g_size * 3 * (size_t)((1 << envs->gbits) + 1);
    size_t leni = nf * i_ctr * n_comp; // gctri
    size_t len0 = nf * n_comp;         // gout
    size_t len = leng + leni + len0;
    double *g;
    MALLOC_INSTACK(g, len);
    double *g1 = g + leng;
    double *gout, *gctri;
    ALIAS_ADDR_IF_EQUAL(i, m);
    gout = g1;

    for (kp = 0; kp < k_prim; kp++)
    {
        envs->ak[0] = ak[kp];
        fac1k = envs->common_factor * ck[kp];
        pdata_ij = pdata_base;
        for (jp = 0; jp < j_prim; jp++)
        {
            envs->aj[0] = aj[jp];
            fac1j = fac1k * cj[jp];
            for (ip = 0; ip < i_prim; ip++, pdata_ij++)
            {
                SET_RIJ;
                cutoff = expcutoff - pdata_ij->cceij;
                fac1i = fac1j * expij;
                envs->fac[0] = fac1i;
                if ((*envs->f_g0_2e)(g, rij, rkl, cutoff, envs))
                {
                    (*envs->f_gout)(gout, g, idx, envs, 1);
                    PRIM2CTR0(i, gout, len0);
                }
            i_contracted:;
            } // end loop i_prim
        } // end loop j_prim
    } // end loop k_prim

    if (n_comp > 1 && !*iempty)
    {
        TRANSPOSE(gctri);
    }
    return !*empty;
}

// j_ctr = n; i_ctr = k_ctr = 1;
int CINT3c2e_1n1_loop(double *gctr, CINTEnvVars *envs, double *cache, int *empty)
{
    COMMON_ENVS_AND_DECLARE;
    ADJUST_CUTOFF;
    int nc = j_ctr;
    size_t leng = (size_t)envs->g_size * 3 * (size_t)((1 << envs->gbits) + 1);
    size_t lenj = nf * j_ctr * n_comp; // gctrj
    size_t len0 = nf * n_comp;         // gout
    size_t len = leng + lenj + len0;
    double *g;
    MALLOC_INSTACK(g, len);
    double *g1 = g + leng;
    double *gout, *gctrj;
    ALIAS_ADDR_IF_EQUAL(j, m);
    gout = g1;

    for (kp = 0; kp < k_prim; kp++)
    {
        envs->ak[0] = ak[kp];
        fac1k = envs->common_factor * ck[kp];
        pdata_ij = pdata_base;
        for (jp = 0; jp < j_prim; jp++)
        {
            envs->aj[0] = aj[jp];
            fac1j = fac1k;
            *iempty = 1;
            for (ip = 0; ip < i_prim; ip++, pdata_ij++)
            {
                SET_RIJ;
                cutoff = expcutoff - pdata_ij->cceij;
                fac1i = fac1j * ci[ip] * expij;
                envs->fac[0] = fac1i;
                if ((*envs->f_g0_2e)(g, rij, rkl, cutoff, envs))
                {
                    (*envs->f_gout)(gout, g, idx, envs, *iempty);
                    *iempty = 0;
                }
            i_contracted:;
            } // end loop i_prim
            if (!*iempty)
            {
                PRIM2CTR0(j, gout, len0);
            }
        } // end loop j_prim
    } // end loop k_prim

    if (n_comp > 1 && !*jempty)
    {
        TRANSPOSE(gctrj);
    }
    return !*empty;
}

// i_ctr = j_ctr = k_ctr = 1;
int CINT3c2e_111_loop(double *gctr, CINTEnvVars *envs, double *cache, int *empty)
{
    COMMON_ENVS_AND_DECLARE;
    ADJUST_CUTOFF;
    int nc = 1;
    size_t leng = (size_t)envs->g_size * 3 * (size_t)((1 << envs->gbits) + 1);
    size_t len0 = (size_t)envs->nf * (size_t)n_comp;
    size_t len = leng + len0;
    double *g;
    MALLOC_INSTACK(g, len);
    double *gout;
    if (n_comp == 1)
    {
        gout = gctr;
        gempty = empty;
    }
    else
    {
        gout = g + leng;
    }

    for (kp = 0; kp < k_prim; kp++)
    {
        envs->ak[0] = ak[kp];
        fac1k = envs->common_factor * ck[kp];

        pdata_ij = pdata_base;
        for (jp = 0; jp < j_prim; jp++)
        {
            envs->aj[0] = aj[jp];
            fac1j = fac1k * cj[jp];
            for (ip = 0; ip < i_prim; ip++, pdata_ij++)
            {
                SET_RIJ;
                cutoff = expcutoff - pdata_ij->cceij;
                fac1i = fac1j * ci[ip] * expij;
                envs->fac[0] = fac1i;
                if ((*envs->f_g0_2e)(g, rij, rkl, cutoff, envs))
                {
                    (*envs->f_gout)(gout, g, idx, envs, *gempty);
                    *gempty = 0;
                }
            i_contracted:;
            } // end loop i_prim
        } // end loop j_prim
    } // end loop k_prim

    if (n_comp > 1 && !*gempty)
    {
        TRANSPOSE(gout);
    }
    return !*empty;
}

static int (*CINTf_3c2e_loop[8])(double *, CINTEnvVars *, double *, int *) = {
    CINT3c2e_loop,
    CINT3c2e_loop,
    CINT3c2e_loop,
    CINT3c2e_n11_loop,
    CINT3c2e_loop,
    CINT3c2e_1n1_loop,
    CINT3c2e_loop,
    CINT3c2e_111_loop,
};

int CINT3c2e_drv(double *out, const int *dims, CINTEnvVars *envs, CINTOpt *opt,
                 double *cache, void (*f_e1_c2s)(double *, double *, const int *, CINTEnvVars *, double *), int is_ssc)
{
    int *x_ctr = envs->x_ctr;
    size_t nc = (size_t)envs->nf * (size_t)x_ctr[0] * (size_t)x_ctr[1] * (size_t)x_ctr[2];
    int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
    if (out == NULL)
    {
        PAIRDATA_NON0IDX_SIZE(pdata_size);
        int leng = envs->g_size * 3 * ((1 << envs->gbits) + 1);
        int len0 = envs->nf * n_comp;
        int cache_size = std::max(leng + len0 + nc * n_comp * 3 + pdata_size,
                                  nc * n_comp + envs->nf * 3);
        return cache_size;
    }
    double *stack = NULL;
    if (cache == NULL)
    {
        PAIRDATA_NON0IDX_SIZE(pdata_size);
        size_t leng = (size_t)envs->g_size * 3 * (size_t)((1 << envs->gbits) + 1);
        size_t len0 = (size_t)envs->nf * (size_t)n_comp;
        size_t cache_size = std::max(leng + len0 + nc * (size_t)n_comp * 3 + pdata_size,
                                     nc * (size_t)n_comp + (size_t)envs->nf * 3);
        stack = (double *)malloc(sizeof(double) * cache_size);
        cache = stack;
    }
    double *gctr;
    MALLOC_INSTACK(gctr, nc * n_comp);

    int n;
    int empty = 1;
    if (opt != NULL)
    {
        envs->opt = opt;
        n = ((envs->x_ctr[0] == 1) << 2) + ((envs->x_ctr[1] == 1) << 1) + (envs->x_ctr[2] == 1);
        CINTf_3c2e_loop[n](gctr, envs, cache, &empty);
    }
    else
    {
        CINT3c2e_loop_nopt(gctr, envs, cache, &empty);
    }

    int counts[4];
    if (f_e1_c2s == &c2s_sph_3c2e1)
    {
        counts[0] = (envs->i_l * 2 + 1) * x_ctr[0];
        counts[1] = (envs->j_l * 2 + 1) * x_ctr[1];
        if (is_ssc)
        {
            counts[2] = envs->nfk * x_ctr[2];
        }
        else
        {
            counts[2] = (envs->k_l * 2 + 1) * x_ctr[2];
        }
    }
    else
    {
        counts[0] = envs->nfi * x_ctr[0];
        counts[1] = envs->nfj * x_ctr[1];
        counts[2] = envs->nfk * x_ctr[2];
    }
    counts[3] = 1;
    if (dims == NULL)
    {
        dims = counts;
    }
    int nout = dims[0] * dims[1] * dims[2];
    if (!empty)
    {
        for (n = 0; n < n_comp; n++)
        {
            (*f_e1_c2s)(out + nout * n, gctr + nc * n, dims, envs, cache);
        }
    }
    else
    {
        for (n = 0; n < n_comp; n++)
        {
            c2s_dset0(out + nout * n, dims, counts);
        }
    }
    if (stack != NULL)
    {
        free(stack);
    }
    return !empty;
}

int int3c2e_sph(double *out, int *dims, int *shls, int *atm, int natm,
                int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
    int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
    CINTEnvVars envs;
    CINTinit_int3c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
    envs.f_gout = &CINTgout2e;
    return CINT3c2e_drv(out, dims, &envs, opt, cache, &c2s_sph_3c2e1, 0);
}

// Function to compute three-center two-electron integrals (eri3c)
void computeEri3c(Int_Params &param1,
                  Int_Params &param2,
                  vec &eri3c)
{
    int nQM = param1.get_nbas();
    int nAux = param2.get_nbas();

    Int_Params combined = param1 + param2;
    // combined.print_data("combined");

    ivec bas = combined.get_bas();
    ivec atm = combined.get_atm();
    vec env = combined.get_env();

    ivec shl_slice = {
        0,
        nQM,
        0,
        nQM,
        nQM,
        nQM + nAux,
    };

    int nat = combined.get_natoms();
    int nbas = combined.get_nbas();

    assert(shl_slice[1] <= nbas);
    assert(shl_slice[3] <= nbas);
    assert(shl_slice[5] <= nbas);

    ivec aoloc = make_loc<COORDINATE_TYPE::SPH>(bas, nbas);
    int naoi = aoloc[shl_slice[1]] - aoloc[shl_slice[0]];
    int naoj = aoloc[shl_slice[3]] - aoloc[shl_slice[2]];
    int naok = aoloc[shl_slice[5]] - aoloc[shl_slice[4]];

    CINTOpt *opty = nullptr;
    int3c2e_optimizer(&opty, atm.data(), nat, bas.data(), nbas, env.data());

    // Compute integrals
    vec res((size_t)naoi * (size_t)naoj * (size_t)naok, 0.0);
    eri3c.resize((size_t)naoi * (size_t)naoj * (size_t)naok, 0.0);

    GTOnr3c_drv(int3c2e_sph, GTOnr3c_fill_s1, res.data(), 1, shl_slice.data(), aoloc.data(), opty, atm.data(), nat, bas.data(), nbas, env.data());

    // FOR TESTING PURPOSES!!!!
    // GTOnr3c_drv(int3c2e_sph, res.data(), 1, shl_slice.data(), aoloc.data(), NULL, atm.data(), nat, bas.data(), nbas, env.data());

    // res is in fortran order, write the result in regular ordering
    for (int k = 0; k < naok; k++)
    {
        for (int j = 0; j < naoj; j++)
        {
            for (int i = 0; i < naoi; i++)
            {
                std::size_t idx_F = i + j * (size_t)naoi + k * ((size_t)naoi * (size_t)naoj);
                std::size_t idx_C = i * ((size_t)naoj * (size_t)naok) + j * (size_t)naok + k;
                eri3c[idx_C] = res[idx_F];
            }
        }
    }
}

int CINT3c1e_loop_nopt(double *gctr, CINTEnvVars *envs, double *cache, int *empty)
{
    int *shls = envs->shls;
    int *bas = envs->bas;
    double *env = envs->env;
    int i_sh = shls[0];
    int j_sh = shls[1];
    int k_sh = shls[2];
    int i_ctr = envs->x_ctr[0];
    int j_ctr = envs->x_ctr[1];
    int k_ctr = envs->x_ctr[2];
    int i_prim = bas(NPRIM_OF, i_sh);
    int j_prim = bas(NPRIM_OF, j_sh);
    int k_prim = bas(NPRIM_OF, k_sh);
    double *ri = envs->ri;
    double *rj = envs->rj;
    double *rk = envs->rk;
    double *ai = env + bas(PTR_EXP, i_sh);
    double *aj = env + bas(PTR_EXP, j_sh);
    double *ak = env + bas(PTR_EXP, k_sh);
    double *ci = env + bas(PTR_COEFF, i_sh);
    double *cj = env + bas(PTR_COEFF, j_sh);
    double *ck = env + bas(PTR_COEFF, k_sh);
    int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
    size_t nf = envs->nf;
    double fac1i, fac1j, fac1k;
    int ip, jp, kp;
    int _empty[4] = {1, 1, 1, 1};
    int *iempty = _empty + 0;
    int *jempty = _empty + 1;
    int *kempty = _empty + 2;
    int *gempty = _empty + 3;

    int *idx;
    MALLOC_INSTACK_INT(idx, envs->nf * 3);
    CINTg2e_index_xyz(idx, envs);

    int *non0ctri, *non0ctrj, *non0ctrk;
    int *non0idxi, *non0idxj, *non0idxk;
    MALLOC_INSTACK_INT(non0ctri, i_prim + j_prim + k_prim + i_prim * i_ctr + j_prim * j_ctr + k_prim * k_ctr);
    non0ctrj = non0ctri + i_prim;
    non0ctrk = non0ctrj + j_prim;
    non0idxi = non0ctrk + k_prim;
    non0idxj = non0idxi + i_prim * i_ctr;
    non0idxk = non0idxj + j_prim * j_ctr;
    CINTOpt_non0coeff_byshell(non0idxi, non0ctri, ci, i_prim, i_ctr);
    CINTOpt_non0coeff_byshell(non0idxj, non0ctrj, cj, j_prim, j_ctr);
    CINTOpt_non0coeff_byshell(non0idxk, non0ctrk, ck, k_prim, k_ctr);

    int nc = i_ctr * j_ctr * k_ctr;
    size_t leng = (size_t)envs->g_size * 3 * (size_t)((1 << envs->gbits) + 1);
    size_t lenk = (size_t)envs->nf * nc * n_comp;            // gctrk
    size_t lenj = (size_t)envs->nf * i_ctr * j_ctr * n_comp; // gctrj
    size_t leni = (size_t)envs->nf * i_ctr * n_comp;         // gctri
    size_t len0 = (size_t)envs->nf * n_comp;                 // gout
    size_t len = leng + lenk + lenj + leni + len0;
    double *g;
    MALLOC_INSTACK(g, len);
    double *g1 = g + leng;
    double *gout, *gctri, *gctrj, *gctrk;

    ALIAS_ADDR_IF_EQUAL(k, m);
    ALIAS_ADDR_IF_EQUAL(j, k);
    ALIAS_ADDR_IF_EQUAL(i, j);
    ALIAS_ADDR_IF_EQUAL(g, i);

    double eijk, dijk, aijk;
    double aiajrr, aiakrr, ajakrr;
    double rirk[3];
    double rjrk[3];
    rirk[0] = ri[0] - rk[0];
    rirk[1] = ri[1] - rk[1];
    rirk[2] = ri[2] - rk[2];
    rjrk[0] = rj[0] - rk[0];
    rjrk[1] = rj[1] - rk[1];
    rjrk[2] = rj[2] - rk[2];
    double rr_ij = SQUARE(envs->rirj);
    double rr_ik = SQUARE(rirk);
    double rr_jk = SQUARE(rjrk);

    for (kp = 0; kp < k_prim; kp++)
    {
        envs->ak[0] = ak[kp];
        if (k_ctr == 1)
        {
            fac1k = envs->common_factor * ck[kp];
        }
        else
        {
            fac1k = envs->common_factor;
            *jempty = 1;
        }

        for (jp = 0; jp < j_prim; jp++)
        {
            envs->aj[0] = aj[jp];
            if (j_ctr == 1)
            {
                fac1j = fac1k * cj[jp];
            }
            else
            {
                fac1j = fac1k;
                *iempty = 1;
            }
            ajakrr = aj[jp] * ak[kp] * rr_jk;
            for (ip = 0; ip < i_prim; ip++)
            {
                envs->ai[0] = ai[ip];
                aijk = ai[ip] + aj[jp] + ak[kp];
                aiakrr = ai[ip] * ak[kp] * rr_ik;
                aiajrr = ai[ip] * aj[jp] * rr_ij;
                eijk = (aiajrr + aiakrr + ajakrr) / aijk;
                if (eijk > EXPCUTOFF)
                {
                    continue;
                }

                if (i_ctr == 1)
                {
                    fac1i = fac1j * ci[ip] * exp(-eijk);
                }
                else
                {
                    fac1i = fac1j * exp(-eijk);
                }
                dijk = fac1i / (aijk * sqrt(aijk));
                envs->fac[0] = dijk;
                CINTg3c1e_ovlp(g, ai[ip], aj[jp], ak[kp], envs);
                (*envs->f_gout)(gout, g, idx, envs, *gempty);

                PRIM2CTR0(i, gout, len0);
            } // end loop i_prim
            if (!*iempty)
            {
                PRIM2CTR0(j, gctri, leni);
            }
        } // end loop j_prim
        if (!*jempty)
        {
            PRIM2CTR0(k, gctrj, lenj);
        }
    } // end loop k_prim

    if (n_comp > 1 && !*kempty)
    {
        TRANSPOSE(gctrk);
    }
    return !*empty;
}

int CINT3c1e_drv(double *out, int *dims, CINTEnvVars *envs, CINTOpt *opt,
                 double *cache, void (*f_e1_c2s)(double *, double *, const int *, CINTEnvVars *, double *), int int_type, int is_ssc)
{
    int *x_ctr = envs->x_ctr;
    int nc = envs->nf * x_ctr[0] * x_ctr[1] * x_ctr[2];
    int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
    if (out == NULL)
    {
        PAIRDATA_NON0IDX_SIZE(pdata_size);
        size_t leng = (size_t)envs->g_size * 3 * ((1 << envs->gbits) + 1);
        size_t len0 = (size_t)envs->nf * n_comp;
        int cache_size = std::max((double)(leng + len0 + (size_t)nc * n_comp * 4 + pdata_size), (double)(nc * n_comp + envs->nf * 3));
        return cache_size;
    }
    double *stack = NULL;
    if (cache == NULL)
    {
        PAIRDATA_NON0IDX_SIZE(pdata_size);
        size_t leng = (size_t)envs->g_size * 3 * ((1 << envs->gbits) + 1);
        size_t len0 = (size_t)envs->nf * n_comp;
        int cache_size = std::max((double)(leng + len0 + (size_t)nc * n_comp * 4 + pdata_size), (double)(nc * n_comp + envs->nf * 3));
        stack = (double *)malloc(sizeof(double) * cache_size);
        cache = stack;
    }
    double *gctr;
    int n;
    int empty = 1;

    if (int_type == 0)
    {
        MALLOC_INSTACK(gctr, nc * n_comp);
        CINT3c1e_loop_nopt(gctr, envs, cache, &empty);
    }
    else if (int_type == 1)
    {
        MALLOC_INSTACK(gctr, nc * n_comp);
        std::cout << "CINT3c1e_nuc_loop_nopt not Implementd!" << std::endl;
        exit(1);
        // CINT3c1e_nuc_loop_nopt(gctr, envs, 1, -1, cache, &empty);
    }
    else
    {
        std::cout << "CINT3c1e_drv: unknown int_type" << std::endl;
        exit(1);
        int *atm = envs->atm;
        //int i;
        //double fac;
        double *buf;
        MALLOC_INSTACK(gctr, nc * n_comp);
        MALLOC_INSTACK(buf, nc * n_comp);
        //for (n = 0; n < envs->natm; n++)
        //{
        //    if (atm(CHARGE_OF, n) != 0)
        //    {
        //        fac = -abs(atm(CHARGE_OF, n));
        //        // CINT3c1e_nuc_loop_nopt(buf, envs, fac, n, cache, &empty);
        //    }
        //}
        cache = buf; // release memory for c2s function
    }

    int counts[4];
    if (f_e1_c2s == &c2s_sph_3c1e || f_e1_c2s == &c2s_sph_3c2e1)
    {
        counts[0] = (envs->i_l * 2 + 1) * x_ctr[0];
        counts[1] = (envs->j_l * 2 + 1) * x_ctr[1];
        if (is_ssc)
        {
            counts[2] = envs->nfk * x_ctr[2];
        }
        else
        {
            counts[2] = (envs->k_l * 2 + 1) * x_ctr[2];
        }
    }
    else
    {
        counts[0] = envs->nfi * x_ctr[0];
        counts[1] = envs->nfj * x_ctr[1];
        counts[2] = envs->nfk * x_ctr[2];
    }
    counts[3] = 1;
    if (dims == NULL)
    {
        dims = counts;
    }
    int nout = dims[0] * dims[1] * dims[2];
    if (!empty)
    {
        for (n = 0; n < n_comp; n++)
        {
            (*f_e1_c2s)(out + nout * n, gctr + nc * n, dims, envs, cache);
        }
    }
    else
    {
        for (n = 0; n < n_comp; n++)
        {
            c2s_dset0(out + nout * n, dims, counts);
        }
    }
    if (stack != NULL)
    {
        free(stack);
    }
    return !empty;
}

void CINTgout1e(double *gout, double *g, int *idx, CINTEnvVars *envs, int empty)
{
    int nf = envs->nf;
    int n, ix, iy, iz;
    if (empty)
    {
        for (n = 0; n < nf; n++)
        {
            ix = idx[n * 3 + 0];
            iy = idx[n * 3 + 1];
            iz = idx[n * 3 + 2];
            gout[n] = g[ix] * g[iy] * g[iz];
        }
    }
    else
    {
        for (n = 0; n < nf; n++)
        {
            ix = idx[n * 3 + 0];
            iy = idx[n * 3 + 1];
            iz = idx[n * 3 + 2];
            gout[n] += g[ix] * g[iy] * g[iz];
        }
    }
}

int int3c1e_sph(double *out, int *dims, int *shls, int *atm, int natm,
                int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
    int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
    CINTEnvVars envs;
    CINTinit_int3c1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
    envs.f_gout = &CINTgout1e;
    return CINT3c1e_drv(out, dims, &envs, opt, cache, &c2s_sph_3c1e, 0, 0); // Second to last indice controls the type of integral (e.g. Overlap)
}

/*
 * 1e GTO integral basic loop for < i|j>, no 1/r
 */
int CINT1e_loop(double *gctr, CINTEnvVars *envs, double *cache, int int1e_type)
{
    int *shls = envs->shls;
    int *bas = envs->bas;
    double *env = envs->env;
    int i_sh = shls[0];
    int j_sh = shls[1];
    int i_ctr = envs->x_ctr[0];
    int j_ctr = envs->x_ctr[1];
    int i_prim = bas(NPRIM_OF, i_sh);
    int j_prim = bas(NPRIM_OF, j_sh);
    double *ai = env + bas(PTR_EXP, i_sh);
    double *aj = env + bas(PTR_EXP, j_sh);
    double *ci = env + bas(PTR_COEFF, i_sh);
    double *cj = env + bas(PTR_COEFF, j_sh);
    int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;

    double expcutoff = envs->expcutoff;
    double *log_maxci, *log_maxcj;
    PairData *pdata_base, *pdata_ij;
    MALLOC_INSTACK(log_maxci, i_prim + j_prim);
    MALLOC_INSTACK_PAIRDATA(pdata_base, i_prim * j_prim);
    log_maxcj = log_maxci + i_prim;
    CINTOpt_log_max_pgto_coeff(log_maxci, ci, i_prim, i_ctr);
    CINTOpt_log_max_pgto_coeff(log_maxcj, cj, j_prim, j_ctr);
    if (CINTset_pairdata(pdata_base, ai, aj, envs->ri, envs->rj,
                         log_maxci, log_maxcj, envs->li_ceil, envs->lj_ceil,
                         i_prim, j_prim, SQUARE(envs->rirj), expcutoff, env))
    {
        return 0;
    }

    double fac1i, fac1j, expij;
    int ip, jp;
    int empty[4] = {1, 1, 1, 1};
    int *gempty = empty + 0;
    int *iempty = empty + 1;
    int *jempty = empty + 2;
    double *rij;
    int *idx;
    MALLOC_INSTACK_INT(idx, envs->nf * 3);
    CINTg1e_index_xyz(idx, envs);

    int *non0ctri, *non0ctrj;
    int *non0idxi, *non0idxj;
    MALLOC_INSTACK_INT(non0ctri, i_prim + j_prim + i_prim * i_ctr + j_prim * j_ctr);
    non0ctrj = non0ctri + i_prim;
    non0idxi = non0ctrj + j_prim;
    non0idxj = non0idxi + i_prim * i_ctr;
    CINTOpt_non0coeff_byshell(non0idxi, non0ctri, ci, i_prim, i_ctr);
    CINTOpt_non0coeff_byshell(non0idxj, non0ctrj, cj, j_prim, j_ctr);

    const int nc = i_ctr * j_ctr;
    // (irys,i,j,k,l,coord,0:1); +1 for nabla-r12
    const int leng = envs->g_size * 3 * ((1 << envs->gbits) + 1);
    const int lenj = envs->nf * nc * n_comp;    // gctrj
    const int leni = envs->nf * i_ctr * n_comp; // gctri
    const int len0 = envs->nf * n_comp;         // gout
    const int len = leng + lenj + leni + len0;
    double *g, *gout, *gctri, *gctrj;
    MALLOC_INSTACK(g, len); // must be allocated last in this function
    double *g1 = g + leng;
    if (n_comp == 1)
    {
        gctrj = gctr;
    }
    else
    {
        gctrj = g1;
        g1 += lenj;
    }
    if (j_ctr == 1)
    {
        gctri = gctrj;
        iempty = jempty;
    }
    else
    {
        gctri = g1;
        g1 += leni;
    }
    if (i_ctr == 1)
    {
        gout = gctri;
        gempty = iempty;
    }
    else
    {
        gout = g1;
    }

    double common_factor = envs->common_factor * CINTcommon_fac_sp(envs->i_l) * CINTcommon_fac_sp(envs->j_l);

    pdata_ij = pdata_base;
    for (jp = 0; jp < j_prim; jp++)
    {
        envs->aj[0] = aj[jp];
        if (j_ctr == 1)
        {
            fac1j = common_factor * cj[jp];
        }
        else
        {
            fac1j = common_factor;
            *iempty = 1;
        }
        for (ip = 0; ip < i_prim; ip++, pdata_ij++)
        {
            if (pdata_ij->cceij > expcutoff)
            {
                continue;
            }
            envs->ai[0] = ai[ip];
            expij = pdata_ij->eij;
            rij = pdata_ij->rij;
            envs->rij[0] = rij[0];
            envs->rij[1] = rij[1];
            envs->rij[2] = rij[2];
            if (i_ctr == 1)
            {
                fac1i = fac1j * ci[ip] * expij;
            }
            else
            {
                fac1i = fac1j * expij;
            }
            envs->fac[0] = fac1i;
            make_g1e_gout(gout, g, idx, envs, *gempty, int1e_type);
            PRIM2CTR0(i, gout, (size_t)envs->nf * n_comp);
        }
        if (!*iempty)
        {
            PRIM2CTR0(j, gctri, (size_t)envs->nf * i_ctr * n_comp);
        }
    }

    if (n_comp > 1 && !*jempty)
    {
        CINTdmat_transpose(gctr, gctrj, (size_t)envs->nf * nc, n_comp);
    }
    return !*jempty;
}

/*
 * 1e integrals <i|O|j> without 1/r
 */
int CINT1e_drv(double *out, int *dims, CINTEnvVars *envs,
               double *cache, void (*f_c2s)(double *, double *, int *, CINTEnvVars *, double *), int int1e_type)
{
    if (out == NULL)
    {
        return int1e_cache_size(envs);
    }
    int *x_ctr = envs->x_ctr;
    int nc = envs->nf * x_ctr[0] * x_ctr[1];
    int n_comp = envs->ncomp_e1 * envs->ncomp_tensor;
    double *stack = NULL;
    if (cache == NULL)
    {
        size_t cache_size = int1e_cache_size(envs);
        stack = (double *)malloc(sizeof(double) * cache_size);
        cache = stack;
    }
    double *gctr;
    MALLOC_INSTACK(gctr, nc * n_comp);

    int has_value = CINT1e_loop(gctr, envs, cache, int1e_type);

    int counts[4];
    if (dims == NULL)
    {
        dims = counts;
    }
    if (f_c2s == &c2s_sph_1e)
    {
        counts[0] = (envs->i_l * 2 + 1) * x_ctr[0];
        counts[1] = (envs->j_l * 2 + 1) * x_ctr[1];
    }
    else if (f_c2s == &c2s_cart_1e)
    {
        counts[0] = envs->nfi * x_ctr[0];
        counts[1] = envs->nfj * x_ctr[1];
    }
    counts[2] = 1;
    counts[3] = 1;
    int nout = dims[0] * dims[1];
    int n;
    if (has_value)
    {
        for (n = 0; n < n_comp; n++)
        {
            (*f_c2s)(out + nout * n, gctr + nc * n, dims, envs, cache);
        }
    }
    else
    {
        for (n = 0; n < n_comp; n++)
        {
            c2s_dset0(out + nout * n, dims, counts);
        }
    }

    if (stack != NULL)
    {
        free(stack);
    }
    return has_value;
}

int int1e_ovlp_sph(double *out, int *dims, int *shls, int *atm, int natm,
                   int *bas, int nbas, double *env, CINTOpt *opt, double *cache)
{
    int ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
    CINTEnvVars envs;
    CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
    envs.f_gout = &CINTgout1e;
    return CINT1e_drv(out, dims, &envs, cache, &c2s_sph_1e, 0);
}

void compute2c_Overlap(Int_Params &params, vec &overlap_2c)
{
    ivec bas = params.get_bas();
    ivec atm = params.get_atm();
    vec env = params.get_env();

    int nbas = params.get_nbas();
    int nat = params.get_natoms();

    ivec shl_slice = {0, nbas, 0, nbas};
    ivec aoloc = make_loc<COORDINATE_TYPE::SPH>(bas, nbas);

    int naoi = aoloc[shl_slice[1]] - aoloc[shl_slice[0]];
    int naoj = aoloc[shl_slice[3]] - aoloc[shl_slice[2]];

    // Compute integrals
    vec res((size_t)naoi * naoj, 0.0);
    overlap_2c.resize((size_t)naoi * naoj, 0.0);
    GTOint2c(int1e_ovlp_sph, res.data(), 1, 0, shl_slice.data(), aoloc.data(), NULL, atm.data(), nat, bas.data(), nbas, env.data());

    // res is in fortran order, write the result in regular ordering
    for (int i = 0; i < naoi; i++)
    {
        for (int j = 0; j < naoj; j++)
        {
            overlap_2c[j * (size_t)naoi + i] = res[i * (size_t)naoj + j];
        }
    }
}

int int1e_ovlp_cart(double* out, int* dims, int* shls, int* atm, int natm,
    int* bas, int nbas, double* env, CINTOpt* opt, double* cache)
{
    int ng[] = { 0, 0, 0, 0, 0, 1, 1, 1 };
    CINTEnvVars envs;
    CINTinit_int1e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
    envs.f_gout = &CINTgout1e;
    return CINT1e_drv(out, dims, &envs, cache, &c2s_cart_1e, 0);
}

void compute2c_Overlap_Cart(Int_Params& params, vec& overlap_2c)
{
    ivec bas = params.get_bas();
    ivec atm = params.get_atm();
    vec env = params.get_env();

    int nbas = params.get_nbas();
    int nat = params.get_natoms();

    ivec shl_slice = { 0, nbas, 0, nbas };
    ivec aoloc = make_loc<COORDINATE_TYPE::CART>(bas, nbas);

    int naoi = aoloc[shl_slice[1]] - aoloc[shl_slice[0]];
    int naoj = aoloc[shl_slice[3]] - aoloc[shl_slice[2]];

    // Compute integrals
    vec res((size_t)naoi * naoj, 0.0);
    overlap_2c.resize((size_t)naoi * naoj, 0.0);
    GTOint2c(int1e_ovlp_cart, res.data(), 1, 0, shl_slice.data(), aoloc.data(), NULL, atm.data(), nat, bas.data(), nbas, env.data());

    // res is in fortran order, write the result in regular ordering
    for (int i = 0; i < naoi; i++)
    {
        for (int j = 0; j < naoj; j++)
        {
            overlap_2c[j * (size_t)naoi + i] = res[i * (size_t)naoj + j];
        }
    }
}

// Function to calculate the number of 3-center 2-electron integrals to compute at once based on the available memory
// naoi = number of basis functions in the first shell
// naoj = number of basis functions in the second shell
// aoloc = Running total of functions to computed based on the order of the basis functions
// nQM = number of basis functions in the QM basis
// nAux = number of basis functions in the auxiliary basis
// max_RAM = maximum available memory in MB
// Returns the number of functions to compute at once
ivec calc_3c_steps(const unsigned long long int naoi, const unsigned long long int naoj, const ivec aoloc, const int nQM, const int nAux, const double max_mem)
{
    ivec steps = {nQM};

    // First check if the maximum memory is enough to compute all integrals at once
    // Calculate minimum memory needed for all integrals
    unsigned long long int naok_end = (size_t)(aoloc[(size_t)nQM + nAux]) - aoloc[nQM];
    double min_mem = static_cast<double>(sizeof(double) * naoi * naoj * naok_end) * 1e-6 + 200; // Small buffer of 200MB for other memory usage
    if (min_mem < max_mem)
    {
        steps.push_back(nQM + nAux);
        return steps;
    }

    // Calculate maximum number of basis functions for every iteration to stay under the memory limit
    int current_step = 0;
    unsigned long long int naok_max = static_cast<unsigned long long int>(max_mem / ((static_cast<double>(sizeof(double) * naoi * naoj)) * 1e-6));
    for (int bas_i = 1; bas_i <= nAux; bas_i++)
    {
        unsigned long long int naok = aoloc[nQM + bas_i] - aoloc[steps[current_step]];
        if (naok > naok_max)
        {
            steps.push_back(nQM + bas_i - 1);
            current_step += 1;
        }
    }
    steps.push_back(nQM + nAux);
    return steps;
}

void computeRho_Coulomb(Int_Params &param1,
                        Int_Params &param2,
                        const dMatrix2 &dm,
                        vec &rho,
                        double max_mem)
{
    int nQM = param1.get_nbas();
    int nAux = param2.get_nbas();

    Int_Params combined = param1 + param2;
    // combined.print_data("combined");

    ivec bas = combined.get_bas();
    ivec atm = combined.get_atm();
    vec env = combined.get_env();

    int nat = combined.get_natoms();
    int nbas = combined.get_nbas();

    ivec aoloc = make_loc<COORDINATE_TYPE::SPH>(bas, nbas);

    rho.resize(param2.get_nao(), 0.0);

    unsigned long long int naoi = aoloc[nQM] - aoloc[0];
    unsigned long long int naoj = aoloc[nQM] - aoloc[0];
    unsigned long long tot_size = naoi * naoj;

    ivec steps = calc_3c_steps(naoi, naoj, aoloc, nQM, nAux, max_mem);

    // Calculate maximum vector size
    double min_mem = static_cast<double>(sizeof(double) * tot_size) / 1024 / 1024;
    unsigned long long int naok_max = static_cast<unsigned long long int>((max_mem - min_mem) / (static_cast<double>(min_mem))) + 5;
    if (naok_max > (aoloc[nQM + nAux] - aoloc[nQM]))
    {
        naok_max = aoloc[nQM + nAux] - aoloc[nQM];
    }

    std::cout << "Preallocating: " << static_cast<double>(sizeof(double) * tot_size * naok_max) / 1024 / 1024 << " MB" << std::endl;
    //vec res(tot_size * naok_max, 0.0);
    std::unique_ptr<double[]> res(new double[tot_size * naok_max]);

    CINTOpt *opty = nullptr;
    int3c2e_optimizer(&opty, atm.data(), nat, bas.data(), nbas, env.data());

    for (int step_idx = 1; step_idx < steps.size(); step_idx++)
    {
        ivec shl_slice = {
            0,
            nQM,
            0,
            nQM,
            steps[step_idx - 1],
            steps[step_idx]};

        unsigned int naok = aoloc[shl_slice[5]] - aoloc[shl_slice[4]];
        std::cout << "Memory needed for rho: " << static_cast<double>(sizeof(double) * tot_size * naok) / 1024 / 1024 << " MB" << std::endl;

        // Compute 3-center integrals
        GTOnr3c_drv(int3c2e_sph, GTOnr3c_fill_s1, res.get(), 1, shl_slice.data(), aoloc.data(), opty, atm.data(), nat, bas.data(), nbas, env.data());
        int idx_curr_rho = aoloc[steps[step_idx - 1]] - aoloc[nQM];
        // einsum('ijk,ij->k', res, dm, out=rho)
        // This is 100% pure magic, thanks ChatGPT
        cblas_dgemv(CblasRowMajor, 
                    CblasNoTrans,
                    naok,
                    tot_size,
                    1.0,
                    res.get(), tot_size,
                    dm.data(), 1,
                    0.0,
                    &rho[idx_curr_rho], 1);

        // Cheekyly not resetting the res vector, as it will be overwritten in the next iteration (Hopefully...)
        // std::fill(res.begin(), res.end(), 0); //This might be unnecessary, saves about 1sec for 15 steps  TEST THIS!!!!!
    }
}

void computeRho_Overlap(Int_Params &param1,
                        Int_Params &param2,
                        const dMatrix2 &dm,
                        vec &rho,
                        double max_mem)
{
    int nQM = param1.get_nbas();
    int nAux = param2.get_nbas();

    Int_Params combined = param1 + param2;
    // combined.print_data("combined");

    ivec bas = combined.get_bas();
    ivec atm = combined.get_atm();
    vec env = combined.get_env();

    int nat = combined.get_natoms();
    int nbas = combined.get_nbas();

    ivec aoloc = make_loc<COORDINATE_TYPE::SPH>(bas, nbas);

    rho.resize(param2.get_nao(), 0.0);

    unsigned long long int naoi = aoloc[nQM] - aoloc[0];
    unsigned long long int naoj = aoloc[nQM] - aoloc[0];

    ivec steps = calc_3c_steps(naoi, naoj, aoloc, nQM, nAux, max_mem);

    // Calculate maximum vector size
    double min_mem = static_cast<double>(sizeof(double) * naoi * naoj) / 1024 / 1024;
    unsigned long long int naok_max = static_cast<unsigned long long int>((max_mem - min_mem) / (static_cast<double>(min_mem))) + 5;
    if (naok_max > (aoloc[nQM + nAux] - aoloc[nQM]))
    {
        naok_max = aoloc[nQM + nAux] - aoloc[nQM];
    }

    std::cout << "Preallocating: " << static_cast<double>(sizeof(double) * naoi * naoj * naok_max) / 1024 / 1024 << " MB" << std::endl;
    vec res(naoi * naoj * naok_max, 0.0);

    for (int step_idx = 1; step_idx < steps.size(); step_idx++)
    {
        ivec shl_slice = {
            0,
            nQM,
            0,
            nQM,
            steps[step_idx - 1],
            steps[step_idx]};

        unsigned long long int naok = aoloc[shl_slice[5]] - aoloc[shl_slice[4]];
        std::cout << "Memory needed for rho: " << static_cast<double>(sizeof(double) * naoi * naoj * naok) / 1024 / 1024 << " MB" << std::endl;

        // Compute 3-center integrals
        GTOnr3c_drv(int3c1e_sph, GTOnr3c_fill_s1, res.data(), 1, shl_slice.data(), aoloc.data(), NULL, atm.data(), nat, bas.data(), nbas, env.data());

        int idx_curr_rho = aoloc[steps[step_idx - 1]] - aoloc[nQM];
        // einsum('ijk,ij->k', res, dm, out=rho)
        // This is 100% pure magic, thanks ChatGPT
        cblas_dgemv(CblasRowMajor, CblasNoTrans,
                    naok,
                    naoi * naoj,
                    1.0,
                    res.data(), naoi * naoj,
                    dm.data(), 1,
                    0.0,
                    &rho[idx_curr_rho], 1);
        // Cheekyly not resetting the res vector, as it will be overwritten in the next iteration (Hopefully...)
        // std::fill(res.begin(), res.end(), 0); //This might be unnecessary, saves about 1sec for 15 steps  TEST THIS!!!!!
    }
}
