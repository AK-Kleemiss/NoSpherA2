#pragma once
#include "atoms.h"
#include "cube.h"
#include <vector>

class WFN;
class cube;

bool b2c(const cube* cub, const std::vector<atom> &atoms, bool debug, bool bcp);
cubei topological_cube_analysis(const cube* cub, const std::vector<atom>& atoms, bool debug, bool bcp, double value_floor = 0.0, double grad_epsilon = 1e-12, double assignment_radius = -1.0);
vec integrate_values_in_basins(const cube *cub, const cubei *basin_cube, bool debug);

#include "wfn_class.h"
