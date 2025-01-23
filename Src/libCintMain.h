#pragma once

#include "convenience.h"
#include "integration_params.h"


void computeEri3c(Int_Params& param1,
    Int_Params& param2,
    vec& eri3c);

// Function to compute three-center two-electron integrals (eri3c)
void computeEri2c(Int_Params& param1,
    vec& eri2c);


//max_RAM in MB
void computeRho(Int_Params& param1,
    Int_Params& param2,
    vec2& dm,
	vec& rho,
    int max_RAM);