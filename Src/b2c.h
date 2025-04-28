#pragma once
#include "atoms.h"
#include <vector>

class WFN;
class cube;

bool b2c(const cube* cub, const std::vector<atom> &atoms, bool debug, bool bcp);

#include "wfn_class.h"
#include "cube.h"
