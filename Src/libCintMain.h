#pragma once

#include "convenience.h"

//Include both, as they are both required when calling the functions declared here
#include "integration_params.h"
#include "libCintKernels.h"

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
 */
template <typename Kernel>
void computeRho(
    const Int_Params& normal_basis,
    const Int_Params& aux_basis,
    const dMatrix2& dm,
    vec& rho);


//DEPRICATED::Function to compute electron repulsion integrals
template <typename Kernel>
void computeEri3c(Int_Params& param1,
    Int_Params& param2,
    vec& eri3c);

template <typename Kernel>
void compute3C(Int_Params& param1,
    Int_Params& param2,
    vec& eri3c);