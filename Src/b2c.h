#pragma once
#include "atoms.h"
#include "cube.h"
#include <vector>

class WFN;
class cube;

struct cubepoint {
    int x;
    int y;
    int z;
    double value;
};

bool b2c(const cube* cub, const std::vector<atom> &atoms, bool debug, bool bcp);
std::pair<cubei, std::vector<d4>> topological_cube_analysis(const cube* cub, const std::vector<atom>& atoms, bool debug, bool bcp, double value_floor = 0.0, double grad_epsilon = 1e-12, double assignment_radius = -1.0);
vec integrate_values_in_basins(const cube *cub, const cubei *basin_cube, svec &basin_label, bool debug);
svec assign_labels_to_basins(const std::vector<d4> &Maxima, const std::vector<atom> &atoms, bool debug, int type_switch = 0);

#include "wfn_class.h"
