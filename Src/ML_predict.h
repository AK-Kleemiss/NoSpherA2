#pragma once
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <set>
#include <algorithm>


#include "convenience.h"
#include "npy.h"
#include "properties.h"
#include "JKFit.h"
//#include "equicomb.h"

#ifdef _WIN32
#include "H5Cpp.h"
#include "rascaline.hpp"
#include "metatensor.h"
#include <H5public.h>
#include <H5Dpublic.h>
#include <H5Opublic.h>
#endif

#include <complex>
#include <cmath>
#include <omp.h>

vec predict(WFN wavy, std::string model_folder);


