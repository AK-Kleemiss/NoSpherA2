#pragma once
#include "atoms.h"
#include "cube.h"
#include <array>
#include <string>
#include <vector>

class WFN;
class cube;

struct cubepoint {
    int x;
    int y;
    int z;
    double value;
};

struct critical_point_seed {
    i3 grid_index;
    d3 position;
    double value;
    double gradient_norm;
    bool is_nuclear_seed = false;
    int nucleus_index = -1;
};

struct critical_point {
    i3 grid_index;
    d3 seed_position;
    d3 position;
    d3 gradient;
    d3 hessian_eigenvalues;
    std::array<d3, 3> hessian_eigenvectors;
    double seed_value;
    double density;
    double gradient_norm;
    double laplacian;
    double ellipticity;
    double virial_field;
    double kinetic_lagrangian;
    double kinetic_hamiltonian;
    double lagrangian_density;
    std::string type;
    int negative_eigenvalues;
    int positive_eigenvalues;
    int zero_eigenvalues;
    int iterations;
    bool converged;
};

bool b2c(const cube* cub, const std::vector<atom> &atoms, bool debug, bool bcp);
std::pair<cubei, std::vector<d4>> topological_cube_analysis(const cube* cub, const std::vector<atom>& atoms, bool debug, bool bcp, double value_floor = 0.0, double grad_epsilon = 1e-12, double assignment_radius = -1.0);
std::vector<critical_point_seed> find_cube_critical_point_seeds(const cube* cub, bool debug, double value_floor = -1.0, double gradient_epsilon = -1.0);
std::vector<critical_point> refine_cube_critical_points(const cube* cub, const WFN& wavy, const std::vector<critical_point_seed>& seeds, bool debug, double value_floor = -1.0, double gradient_tolerance = 1e-8, double step_tolerance = 1e-6, int max_iterations = 32);
std::vector<critical_point> analyze_cube_critical_points(const cube* cub, const WFN& wavy, bool debug, double value_floor = -1.0, double gradient_epsilon = -1.0, double gradient_tolerance = 1e-8, double step_tolerance = 1e-6, int max_iterations = 32);
vec integrate_values_in_basins(const cube *cub, const cubei *basin_cube, svec &basin_label, bool debug);
svec assign_labels_to_basins(const std::vector<d4> &Maxima, const std::vector<atom> &atoms, bool debug, int type_switch = 0);

#include "wfn_class.h"
