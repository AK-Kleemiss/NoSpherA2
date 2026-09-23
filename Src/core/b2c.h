#pragma once
#include "atoms.h"
#include "cube.h"
#include <array>
#include <functional>
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
//The density the basin code follows off the grid. The wavefunction's orbitals unless one of
//these is given, which -ri_fit and -SALTED do with the fitted density: one loop over the
//auxiliary functions for rho and its gradient instead of a sum over orbitals
struct density_field {
	std::function<double(const d3&)> rho;
	std::function<void(const d3&, d3&)> grad;
};
//core_density and core_gradient, when given, are added to the wavefunction's density and
//gradient wherever the field is followed: the spherical core an ECP took out
std::pair<cubei, std::vector<d4>> topological_cube_analysis(const cube* cub, const std::vector<atom>& atoms, bool debug, bool bcp, double value_floor = 0.0, double grad_epsilon = 1e-12, double assignment_radius = -1.0, double merge_persistence = 5e-3, const std::vector<d3>* seeds = nullptr, const WFN* field_wfn = nullptr, const std::function<double(const d3&)>* core_density = nullptr, const std::function<void(const d3&, d3&)>* core_gradient = nullptr, const density_field* field = nullptr);
//Radius holding the ELI-D maxima of an atom's core shells, by period; the outermost shell the
//element keeps beneath its valence peaks at about 0.7 bohr for the first transition row
double core_shell_radius(const int Z);
//Every basin whose maximum lies within an atom's core radius becomes that atom's one core
//basin, as DGrid's ELIDcore does; returns the number of basins merged away
int unify_core_basins(cubei& basin_cube, std::vector<d4>& maxima, const std::vector<atom>& atoms);
//Atomic overlap matrices S^b_ij = int_b phi_i phi_j, taken on the same quadrature points and
//with the same basin assignment as the populations, so a basin's trace is its population by
//construction. One packed lower triangle per basin over the occupied MOs
struct basin_overlaps {
	int nmo = 0;   //MOs in the triangles
	ivec mo_index; //their indices in the WFN, size nmo
	vec2 S;        //S[basin][packed(i,j)], basins in the order of the maxima
	static size_t packed(const int i, const int j) { return i >= j ? (size_t)i * (i + 1) / 2 + j : (size_t)j * (j + 1) / 2 + i; }
	size_t triangle() const { return (size_t)nmo * (nmo + 1) / 2; }
	double at(const int b, const int i, const int j) const { return S[b][packed(i, j)]; }
};
//delta(A,B) = 2 m sum_{ij, same spin} n_i n_j S^A_ij S^B_ij, lambda(A) the same with B = A,
//over spin-orbital occupations n in [0,1]; m = 2 for a restricted wavefunction, whose single
//set of MOs stands for both spins, 1 when alpha and beta are listed separately
struct delocalization_result {
	vec lambda;                            //localization index per basin
	vec population;                        //m sum_i n_i S^A_ii, the trace population
	std::vector<std::array<int, 2>> pairs; //basin pairs, first < second
	vec di;                                //delta for each pair, same order
	double identity_error = 0.0;           //max |sum_A S^A_ij - delta_ij|: the integration's own error
};
delocalization_result delocalization_indices(const WFN& wavy, const basin_overlaps& ovl);
void report_delocalization(const WFN& wavy, const basin_overlaps& ovl, const svec& labels, std::ostream& log, const double threshold = 0.01);
//ovl, when given, is filled with the basin overlap matrices. Only meaningful for the orbital
//density (field == nullptr): it is the orbitals that are being partitioned
vec integrate_basins_on_atomic_grids(const cube* cub, const cubei* basin_cube, const std::vector<d4>& maxima, const WFN& wavy, const int accuracy, const bool eli_field, vec& volumes, double& outside, const std::function<double(const d3&)>* core_density = nullptr, const std::function<void(const d3&, d3&)>* core_gradient = nullptr, const int grid_boost = 1, const density_field* field = nullptr, basin_overlaps* ovl = nullptr);
std::vector<critical_point_seed> find_cube_critical_point_seeds(const cube* cub, bool debug, double value_floor = -1.0, double gradient_epsilon = -1.0);
std::vector<critical_point> refine_cube_critical_points(const cube* cub, const WFN& wavy, const std::vector<critical_point_seed>& seeds, bool debug, double value_floor = -1.0, double gradient_tolerance = 1e-8, double step_tolerance = 1e-6, int max_iterations = 32);
std::vector<critical_point> analyze_cube_critical_points(const cube* cub, const WFN& wavy, bool debug, double value_floor = -1.0, double gradient_epsilon = -1.0, double gradient_tolerance = 1e-8, double step_tolerance = 1e-6, int max_iterations = 32);
vec integrate_values_in_basins(const cube *cub, const cubei *basin_cube, svec &basin_label, bool debug);
svec assign_labels_to_basins(const std::vector<d4> &Maxima, const std::vector<atom> &atoms, bool debug, int type_switch = 0);

#include "wfn_class.h"
