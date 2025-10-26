#pragma once

#include "convenience.h"
#include "integration_params.h"


//Functions to compute electron repulsion integrals
void computeEri3c(Int_Params& param1,
    Int_Params& param2,
    vec& eri3c);

void computeEri2c(Int_Params& param1,
    vec& eri2c);


// Function to compute Overlap integrals
void compute2c_Overlap(Int_Params& param1,
    vec& overlap_2c);

// Function to compute Overlap integrals
void compute2c_Overlap_Cart(Int_Params& param1,
    vec& overlap_2c);

/*
* Function to compute fitting coefficients for the density matrix using the Coulomb metric
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
void computeRho_Coulomb(Int_Params& param1,
    Int_Params& param2,
    const dMatrix2& dm,
    vec& rho,
    double max_mem);

/*
* Function to compute fitting coefficients for the density matrix using the Overlap Metric
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
void computeRho_Overlap(Int_Params& param1,
    Int_Params& param2,
    const dMatrix2& dm,
    vec& rho,
    double max_mem);