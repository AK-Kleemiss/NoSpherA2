#pragma once

//This is really stupid... but at least it works.
//Although this will probably break complex support in libCint if ever needed!
#if !defined(_MSC_VER)
#define _FAKE_MSC_VER
#define _MSC_VER 1930   // any reasonable MSVC version
#endif

typedef struct _Dcomplex {
    double _Val[2];   // [0]=real, [1]=imag
} _Dcomplex;
#define _SILENCE_CXX17_C_HEADER_DEPRECATION_WARNING
#pragma warning(push)
#pragma warning(disable:4996)

extern "C" {
    #include "cint.h"
}

#pragma warning(pop)
#ifdef _FAKE_MSC_VER
#undef _MSC_VER
#undef _FAKE_MSC_VER
#endif


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
DECLARE_CINT_KERNEL(Coulomb2C_CRT, true);
DECLARE_CINT_KERNEL(Overlap2C_SPH, false);
DECLARE_CINT_KERNEL(Overlap2C_CRT, false);
DECLARE_CINT_KERNEL(Coulomb3C_SPH, true);
DECLARE_CINT_KERNEL(Coulomb3C_CRT, true);
DECLARE_CINT_KERNEL(Overlap3C_SPH, false);

#undef DECLARE_CINT_KERNEL