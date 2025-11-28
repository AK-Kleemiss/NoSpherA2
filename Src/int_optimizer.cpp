
#include "int_optimizer.h"
#include "int_g2e.h"

#define NOVALUE                 ((void *)0xffffffffffffffffuL)
#define MAX_PGTO_FOR_PAIRDATA   2048

using namespace NoSpherA2;
int CINTset_pairdata(NoSpherA2::PairData* pairdata, double* ai, double* aj, double* ri, double* rj,
    double* log_maxci, double* log_maxcj,
    int li_ceil, int lj_ceil, int iprim, int jprim,
    double rr_ij, double expcutoff, double* env)
{
    int ip, jp, n;
    double aij, eij, cceij, wj;
    // Normally
    //    (aj*d/sqrt(aij)+1)^li * (ai*d/sqrt(aij)+1)^lj
    //    * pi^1.5/aij^{(li+lj+3)/2} * exp(-ai*aj/aij*rr_ij)
    // is a good approximation for overlap integrals.
    //    <~ (aj*d/aij+1/sqrt(aij))^li * (ai*d/aij+1/sqrt(aij))^lj * (pi/aij)^1.5
    //    <~ (d+1/sqrt(aij))^(li+lj) * (pi/aij)^1.5
    aij = ai[iprim - 1] + aj[jprim - 1];
    double log_rr_ij = 1.7 - 1.5 * std::log(aij);
    int lij = li_ceil + lj_ceil;
    if (lij > 0) {
        double dist_ij = sqrt(rr_ij);
        double omega = env[PTR_RANGE_OMEGA];
        if (omega < 0) {
            double r_guess = 8.;
            double omega2 = omega * omega;
            double theta = omega2 / (omega2 + aij);
            log_rr_ij += lij * std::log(dist_ij + theta * r_guess + 1.);
        }
        else {
            log_rr_ij += lij * std::log(dist_ij + 1.);
        }
    }
    NoSpherA2::PairData* pdata;

    int empty = 1;
    for (n = 0, jp = 0; jp < jprim; jp++) {
        for (ip = 0; ip < iprim; ip++, n++) {
            aij = 1 / (ai[ip] + aj[jp]);
            eij = rr_ij * ai[ip] * aj[jp] * aij;
            cceij = eij - log_rr_ij - log_maxci[ip] - log_maxcj[jp];
            pdata = pairdata + n;
            pdata->cceij = cceij;
            if (cceij < expcutoff) {
                empty = 0;
                wj = aj[jp] * aij;
                pdata->rij[0] = ri[0] + wj * (rj[0] - ri[0]);
                pdata->rij[1] = ri[1] + wj * (rj[1] - ri[1]);
                pdata->rij[2] = ri[2] + wj * (rj[2] - ri[2]);
                pdata->eij = exp(-eij);
            }
            else {
                pdata->rij[0] = 1e18;
                pdata->rij[1] = 1e18;
                pdata->rij[2] = 1e18;
                pdata->eij = 0;
            }
        }
    }
    return empty;
}


void CINTOpt_log_max_pgto_coeff(double* log_maxc, double* coeff, int nprim, int nctr)
{
    int i, ip;
    double maxc;
    for (ip = 0; ip < nprim; ip++) {
        maxc = 0;
        for (i = 0; i < nctr; i++) {
            maxc = std::max(maxc, fabs(coeff[i * nprim + ip]));
        }
        log_maxc[ip] = std::log(maxc);
    }
}

void CINTOpt_set_log_maxc(NoSpherA2::CINTOpt* opt, int* atm, int natm,
    int* bas, int nbas, double* env)
{
    int i, iprim, ictr;
    double* ci;
    size_t tot_prim = 0;
    for (i = 0; i < nbas; i++) {
        tot_prim += bas(NPRIM_OF, i);
    }
    if (tot_prim == 0) {
        return;
    }

    opt->log_max_coeff = (double**)malloc(sizeof(double*) * std::max(nbas, 1));
    double* plog_maxc = (double*)malloc(sizeof(double) * tot_prim);
    opt->log_max_coeff[0] = plog_maxc;
    for (i = 0; i < nbas; i++) {
        iprim = bas(NPRIM_OF, i);
        ictr = bas(NCTR_OF, i);
        ci = env + bas(PTR_COEFF, i);
        opt->log_max_coeff[i] = plog_maxc;
        CINTOpt_log_max_pgto_coeff(plog_maxc, ci, iprim, ictr);
        plog_maxc += iprim;
    }
}

void CINTOpt_non0coeff_byshell(int* sortedidx, int* non0ctr, double* ci,
    int iprim, int ictr)
{
    int ip, j, k, kp;
    //int zeroidx[ictr];
    //std::vector<int> zeroidx(ictr);
    int zeroidx[32];  //THIS IS A TEST, HOPE AND PRAY THAT zeroidx is not larger than 32!!!
    for (ip = 0; ip < iprim; ip++) {
        for (j = 0, k = 0, kp = 0; j < ictr; j++) {
            if (ci[iprim * j + ip] != 0) {
                sortedidx[k] = j;
                k++;
            }
            else {
                zeroidx[kp] = j;
                kp++;
            }
        }
        // Append the index of zero-coeff to sortedidx for function CINTprim_to_ctr_0
        for (j = 0; j < kp; j++) {
            sortedidx[k + j] = zeroidx[j];
        }
        non0ctr[ip] = k;
        sortedidx += ictr;
    }
}

void CINTOpt_set_non0coeff(NoSpherA2::CINTOpt* opt, int* atm, int natm,
    int* bas, int nbas, double* env)
{
    int i, iprim, ictr;
    double* ci;
    size_t tot_prim = 0;
    size_t tot_prim_ctr = 0;
    for (i = 0; i < nbas; i++) {
        tot_prim += bas(NPRIM_OF, i);
        tot_prim_ctr += (size_t)bas(NPRIM_OF, i) * bas(NCTR_OF, i);
    }
    if (tot_prim == 0) {
        return;
    }

    opt->non0ctr = (int**)malloc(sizeof(int*) * std::max(nbas, 1));
    opt->sortedidx = (int**)malloc(sizeof(int*) * std::max(nbas, 1));
    int* pnon0ctr = (int*)malloc(sizeof(int) * tot_prim);
    int* psortedidx = (int*)malloc(sizeof(int) * tot_prim_ctr);
    opt->non0ctr[0] = pnon0ctr;
    opt->sortedidx[0] = psortedidx;
    for (i = 0; i < nbas; i++) {
        iprim = bas(NPRIM_OF, i);
        ictr = bas(NCTR_OF, i);
        ci = env + bas(PTR_COEFF, i);
        opt->non0ctr[i] = pnon0ctr;
        opt->sortedidx[i] = psortedidx;
        CINTOpt_non0coeff_byshell(psortedidx, pnon0ctr, ci, iprim, ictr);
        pnon0ctr += iprim;
        psortedidx += iprim * ictr;
    }
}




// generate caller to CINTinit_2e_optimizer for each type of function
void CINTinit_2e_optimizer(NoSpherA2::CINTOpt** opt, int* atm, int natm,
    int* bas, int nbas, double* env)
{
    NoSpherA2::CINTOpt* opt0 = (NoSpherA2::CINTOpt*)malloc(sizeof(NoSpherA2::CINTOpt));
    opt0->index_xyz_array = NULL;
    opt0->non0ctr = NULL;
    opt0->sortedidx = NULL;
    opt0->nbas = nbas;
    opt0->log_max_coeff = NULL;
    opt0->pairdata = NULL;
    *opt = opt0;
}

void CINTOpt_setij(NoSpherA2::CINTOpt* opt, int* ng,
    int* atm, int natm, int* bas, int nbas, double* env)
{
    int i, j, ip, jp;
    int iprim, jprim, li, lj;
    double* ai, * aj, * ri, * rj;
    double expcutoff;
    if (env[PTR_EXPCUTOFF] == 0) {
        expcutoff = EXPCUTOFF;
    }
    else {
        expcutoff = std::max((double)MIN_EXPCUTOFF, env[PTR_EXPCUTOFF]);
    }

    if (opt->log_max_coeff == NULL) {
        CINTOpt_set_log_maxc(opt, atm, natm, bas, nbas, env);
    }
    double** log_max_coeff = opt->log_max_coeff;
    double* log_maxci, * log_maxcj;

    size_t tot_prim = 0;
    for (i = 0; i < nbas; i++) {
        tot_prim += bas(NPRIM_OF, i);
    }
    if (tot_prim == 0 || tot_prim > MAX_PGTO_FOR_PAIRDATA) {
        return;
    }
    opt->pairdata = (NoSpherA2::PairData**)malloc(sizeof(NoSpherA2::PairData*) * std::max(nbas * nbas, 1));
    NoSpherA2::PairData* pdata = (NoSpherA2::PairData*)malloc(sizeof(NoSpherA2::PairData) * tot_prim * tot_prim);
    opt->pairdata[0] = pdata;

    int ijkl_inc;
    if ((ng[IINC] + ng[JINC]) > (ng[KINC] + ng[LINC])) {
        ijkl_inc = ng[IINC] + ng[JINC];
    }
    else {
        ijkl_inc = ng[KINC] + ng[LINC];
    }

    int empty;
    double rr;
    NoSpherA2::PairData* pdata0;
    for (i = 0; i < nbas; i++) {
        ri = env + atm(PTR_COORD, bas(ATOM_OF, i));
        ai = env + bas(PTR_EXP, i);
        iprim = bas(NPRIM_OF, i);
        li = bas(ANG_OF, i);
        log_maxci = log_max_coeff[i];

        for (j = 0; j <= i; j++) {
            rj = env + atm(PTR_COORD, bas(ATOM_OF, j));
            aj = env + bas(PTR_EXP, j);
            jprim = bas(NPRIM_OF, j);
            lj = bas(ANG_OF, j);
            log_maxcj = log_max_coeff[j];
            rr = (ri[0] - rj[0]) * (ri[0] - rj[0])
                + (ri[1] - rj[1]) * (ri[1] - rj[1])
                + (ri[2] - rj[2]) * (ri[2] - rj[2]);

            empty = CINTset_pairdata(pdata, ai, aj, ri, rj, log_maxci, log_maxcj,
                li + ijkl_inc, lj, iprim, jprim, rr, expcutoff, env);
            if (i == 0 && j == 0) {
                opt->pairdata[0] = pdata;
                pdata += iprim * jprim;
            }
            else if (!empty) {
                opt->pairdata[i * nbas + j] = pdata;
                pdata += iprim * jprim;
                if (i != j) {
                    opt->pairdata[j * nbas + i] = pdata;
                    pdata0 = opt->pairdata[i * nbas + j];
                    // transpose pairdata
                    for (ip = 0; ip < iprim; ip++) {
                        for (jp = 0; jp < jprim; jp++, pdata++) {
                            memcpy(pdata, pdata0 + jp * iprim + ip,
                                sizeof(NoSpherA2::PairData));
                        }
                    }
                }
            }
            else {
                opt->pairdata[i * nbas + j] = (NoSpherA2::PairData*)NOVALUE;
                opt->pairdata[j * nbas + i] = (NoSpherA2::PairData*)NOVALUE;
            }
        }
    }
}

static int _make_fakebas(int* fakebas, int* bas, int nbas, double* env)
{
    int i;
    int max_l = 0;
    for (i = 0; i < nbas; i++) {
        max_l = std::max(max_l, bas(ANG_OF, i));
    }

    int fakenbas = max_l + 1;
    for (i = 0; i < BAS_SLOTS * fakenbas; i++) {
        fakebas[i] = 0;
    }
    // fakebas only initializes ANG_OF, since the others does not
    // affect index_xyz
    for (i = 0; i <= max_l; i++) {
        fakebas[BAS_SLOTS * i + ANG_OF] = i;
    }
    return max_l;
}

static int* _allocate_index_xyz(NoSpherA2::CINTOpt* opt, int max_l, int l_allow, int order)
{
    int i;
    int cumcart = (l_allow + 1) * (l_allow + 2) * (l_allow + 3) / 6;
    size_t ll = max_l + 1;
    size_t cc = cumcart;
    for (i = 1; i < order; i++) {
        ll *= LMAX1;
        cc *= cumcart;
    }
    int* buf = (int*)malloc(sizeof(int) * cc * 3);
    int** ppbuf = (int**)malloc(sizeof(int*) * ll);
    ppbuf[0] = buf;
    for (i = 1; i < ll; i++) {
        ppbuf[i] = NULL;
    }
    opt->index_xyz_array = ppbuf;
    return buf;
}

static void gen_idx(NoSpherA2::CINTOpt* opt, void (*finit)(NoSpherA2::CINTEnvVars*, int*, int*, int*, int, int*, int, double*), void (*findex_xyz)(int *,const NoSpherA2::CINTEnvVars *),
    int order, int l_allow, int* ng,
    int* atm, int natm, int* bas, int nbas, double* env)
{
    int i, j, k, l, ptr;
    int fakebas[BAS_SLOTS * LMAX1];
    int max_l = _make_fakebas(fakebas, bas, nbas, env);
    int fakenbas = max_l + 1;
    // index_xyz bufsize may blow up for large max_l
    l_allow = std::min(max_l, l_allow);
    int* buf = _allocate_index_xyz(opt, max_l, l_allow, order);

    NoSpherA2::CINTEnvVars envs;
    int shls[4] = { 0, };
    if (order == 2) {
        for (i = 0; i <= l_allow; i++) {
            for (j = 0; j <= l_allow; j++) {
                shls[0] = i; shls[1] = j;
                (*finit)(&envs, ng, shls, atm, natm, fakebas, fakenbas, env);
                ptr = i * LMAX1 + j;
                opt->index_xyz_array[ptr] = buf;
                (*findex_xyz)(buf, &envs);
                buf += envs.nf * 3;
            }
        }

    }
    else if (order == 3) {
        for (i = 0; i <= l_allow; i++) {
            for (j = 0; j <= l_allow; j++) {
                for (k = 0; k <= l_allow; k++) {
                    shls[0] = i; shls[1] = j; shls[2] = k;
                    (*finit)(&envs, ng, shls, atm, natm, fakebas, fakenbas, env);
                    ptr = i * LMAX1 * LMAX1 + j * LMAX1 + k;
                    opt->index_xyz_array[ptr] = buf;
                    (*findex_xyz)(buf, &envs);
                    buf += envs.nf * 3;
                }
            }
        }

    }
    else {
        for (i = 0; i <= l_allow; i++) {
            for (j = 0; j <= l_allow; j++) {
                for (k = 0; k <= l_allow; k++) {
                    for (l = 0; l <= l_allow; l++) {
                        shls[0] = i; shls[1] = j; shls[2] = k; shls[3] = l;
                        (*finit)(&envs, ng, shls, atm, natm, fakebas, fakenbas, env);
                        ptr = i * LMAX1 * LMAX1 * LMAX1
                            + j * LMAX1 * LMAX1
                            + k * LMAX1
                            + l;
                        opt->index_xyz_array[ptr] = buf;
                        (*findex_xyz)(buf, &envs);
                        buf += envs.nf * 3;
                    }
                }
            }
        }
    }
}


void CINTall_2c2e_optimizer(NoSpherA2::CINTOpt** opt, int* ng,
    int* atm, int natm, int* bas, int nbas, double* env)
{
    CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
    CINTOpt_set_log_maxc(*opt, atm, natm, bas, nbas, env);
    CINTOpt_set_non0coeff(*opt, atm, natm, bas, nbas, env);
    gen_idx(*opt, &CINTinit_int2c2e_EnvVars, &CINTg1e_index_xyz,
        2, ANG_MAX, ng, atm, natm, bas, nbas, env);
}

void int2c2e_optimizer(NoSpherA2::CINTOpt** opt, int* atm, int natm,
    int* bas, int nbas, double* env)
{
    int ng[] = { 0, 0, 0, 0, 0, 1, 1, 1 };
    CINTall_2c2e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}


void CINTall_3c2e_optimizer(NoSpherA2::CINTOpt** opt, int* ng,
    int* atm, int natm, int* bas, int nbas, double* env)
{
    CINTinit_2e_optimizer(opt, atm, natm, bas, nbas, env);
    CINTOpt_setij(*opt, ng, atm, natm, bas, nbas, env);
    CINTOpt_set_non0coeff(*opt, atm, natm, bas, nbas, env);
    gen_idx(*opt, &CINTinit_int3c2e_EnvVars, &CINTg2e_index_xyz,
        3, 12, ng, atm, natm, bas, nbas, env);
}

void int3c2e_optimizer(NoSpherA2::CINTOpt** opt, int* atm, int natm,
    int* bas, int nbas, double* env)
{
    int ng[] = { 0, 0, 0, 0, 0, 1, 1, 1 };
    CINTall_3c2e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}
