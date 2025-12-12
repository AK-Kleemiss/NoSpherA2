#pragma once
#include "integration_params.h"

void CINTinit_int1e_EnvVars(NoSpherA2::CINTEnvVars* envs, int* ng, int* shls,
                            int* atm, int natm,
                            int* bas, int nbas, double* env);

void CINTinit_int3c1e_EnvVars(NoSpherA2::CINTEnvVars* envs, int* ng, int* shls,
    int* atm, int natm, int* bas, int nbas, double* env);

void CINTinit_int2c2e_EnvVars(NoSpherA2::CINTEnvVars* envs, int* ng, int* shls,
    int* atm, int natm, int* bas, int nbas, double* env);

void CINTinit_int3c2e_EnvVars(NoSpherA2::CINTEnvVars* envs, int* ng, int* shls,
    int* atm, int natm, int* bas, int nbas, double* env);

void CINTg1e_index_xyz(int* idx, const NoSpherA2::CINTEnvVars* envs);

void CINTg2e_index_xyz(int* idx, const NoSpherA2::CINTEnvVars* envs);

void CINTg3c1e_ovlp(double* g, double ai, double aj, double ak,
    const NoSpherA2::CINTEnvVars* envs);

void make_g1e_gout(double* gout, double* g, int* idx,
    NoSpherA2::CINTEnvVars* envs, int empty, int int1e_type);
