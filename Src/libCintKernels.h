#pragma once
#include "cint.h"

#define DECLARE_CINT_KERNEL(NAME, NEEDS_OPT) \
    struct NAME { \
        static constexpr bool NeedsOpt = NEEDS_OPT; \
        static void optimizer (CINTOpt*& opt, \
        int* atm, int nat, int* bas, int nbas, double* env); \
        static void drv(double* out, int comp, int* shl_slice, int* aoloc, \
            CINTOpt* opt, int* atm, int nat, int* bas, int nbas, double* env); \
        static ivec gen_loc(ivec& bas, int nbas); \
    };

/*
Declares a libCint Kernel struct with optimizer and drv functions.
With functions:
- optimizer: initializes the CINTOpt structure for the kernel.
- drv: computes the integrals for the kernel.
- gen_loc: generates the ao_loc array for the kernel.
*/

DECLARE_CINT_KERNEL(Coulomb2C_SPH, true);
DECLARE_CINT_KERNEL(Overlap2C_SPH, false);
DECLARE_CINT_KERNEL(Overlap2C_CRT, false);

DECLARE_CINT_KERNEL(Coulomb3C_SPH, true);
DECLARE_CINT_KERNEL(Overlap3C_SPH, false);

#undef DECLARE_CINT_KERNEL