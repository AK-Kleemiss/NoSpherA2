#pragma once
#include "integration_params.h"

using namespace NoSpherA2;
int CINTset_pairdata(NoSpherA2::PairData* pairdata, double* ai, double* aj, double* ri, double* rj,

    double* log_maxci, double* log_maxcj,
    int li_ceil, int lj_ceil, int iprim, int jprim,
    double rr_ij, double expcutoff, double* env);

void CINTOpt_non0coeff_byshell(int* sortedidx, int* non0ctr, double* ci,
    int iprim, int ictr);

void int3c2e_optimizer(NoSpherA2::CINTOpt** opt, int* atm, int natm,
    int* bas, int nbas, double* env);

void int2c2e_optimizer(NoSpherA2::CINTOpt** opt, int* atm, int natm,
    int* bas, int nbas, double* env);

void CINTOpt_log_max_pgto_coeff(double* log_maxc, double* coeff, int nprim, int nctr);
