#pragma once
#include "integration_params.h"


void CINTinit_int2c2e_EnvVars(CINTEnvVars* envs, int* ng, int* shls,
    int* atm, int natm, int* bas, int nbas, double* env);

void CINTinit_int3c2e_EnvVars(CINTEnvVars* envs, int* ng, int* shls,
    int* atm, int natm, int* bas, int nbas, double* env);

void CINTg1e_index_xyz(int* idx, const CINTEnvVars* envs);

void CINTg2e_index_xyz(int* idx, const CINTEnvVars* envs);