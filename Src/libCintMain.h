#pragma once

#include "convenience.h"
#include "integration_params.h"


//DEPRICATED::Function to compute electron repulsion integrals
void computeEri3c(Int_Params& param1,
    Int_Params& param2,
    vec& eri3c);

struct Coulomb2C {
    static constexpr bool NeedsOpt = true;
    static void optimizer(CINTOpt*& opt,
        int* atm, int nat, int* bas, int nbas, double* env);
    static void drv(double* out, int comp, int* shl_slice, int* aoloc,
        CINTOpt* opt, int* atm, int nat, int* bas, int nbas, double* env);
};

struct Overlap2C {
    static constexpr bool NeedsOpt = false;
    static void optimizer(CINTOpt*& opt,
        const int*, int, const int*, int, const double*);
    static void drv(double* out, int comp, int* shl_slice, int* aoloc,
        CINTOpt* opt, int* atm, int nat, int* bas, int nbas, double* env);
};

/**
 * @brief Computes 2-center integrals using the specified kernel type.
 * 
 * This function calculates 2-center integrals (overlap or Coulomb) based on the
 * integration parameters provided. The specific integral type is determined by
 * the kernel implementation used internally.
 * 
 * @param params Integration parameters containing basis set information, 
 *               atomic coordinates, and computational settings
 * @param ret Output vector to store the computed 2-center integral values
 * 
 * @note The function automatically handles optimization setup and memory 
 *       management based on the kernel requirements (NeedsOpt flag)
 * @see Coulomb2C, Overlap2C for supported 2-center integral types
 */
template <typename Kernel>
void compute2C(Int_Params& params, vec& ret);


struct Coulomb3C {
    static constexpr bool NeedsOpt = true;

    static void optimizer(CINTOpt*& opt,
        int* atm, int nat, int* bas, int nbas, double* env);

    static void drv(double* out, int comp, int* shl_slice, int* aoloc,
        CINTOpt* opt, int* atm, int nat, int* bas, int nbas, double* env);
};

struct Overlap3C {
    static constexpr bool NeedsOpt = false;

    static void optimizer(CINTOpt*& opt,
        const int*, int, const int*, int, const double*);

        static void drv(double* out, int comp, int* shl_slice, int* aoloc,
            CINTOpt* opt, int* atm, int nat, int* bas, int nbas, double* env);
};

// Function to compute Overlap integrals
void compute2c_Overlap_Cart(Int_Params& param1,
    vec& overlap_2c);

/*
* Function to compute fitting coefficients for the density matrix using the desired metric indicated by the Kernel template parameter
* This function combines the calculation of 3-center 2-electron integrals and the fitting coefficients
*
* SUM^K (ij | k) dm_k = c_k
*
* As the matrix (ij | k) can become very large, the computation is performed in steps to avoid memory issues
* Depending on the memory available, chunks for only specifc k are computed at once
*
 * param1: Int_Params object for the QM basis
 * param2: Int_Params object for the auxiliary basis
 * dm: Density matrix
 * rho: Resulting density matrix
 * max_mem: Maximum available memory in MB
 */
template <typename Kernel>
void computeRho(Int_Params& param1,
    Int_Params& param2,
    const dMatrix2& dm,
    vec& rho,
    double max_mem);
