#pragma once
#include "atoms.h"
#include "cube.h"
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

class WFN;
class cube;

constexpr double basin_density_cutoff = 1e-4;

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
//rho, its Hessian eigensystem, the CP type and V/G/K/L at an already refined point
critical_point evaluate_critical_point(const critical_point_seed& seed, const d3& position, const WFN& wavy, int iterations, bool converged);
//The density the basin code follows off the grid. The wavefunction's orbitals unless one of
//these is given, which -ri_fit and -SALTED do with the fitted density: one loop over the
//auxiliary functions for rho and its gradient instead of a sum over orbitals
struct density_field {
	std::function<double(const d3&)> rho;
	std::function<void(const d3&, d3&)> grad;
};
//core_density/core_gradient, the spherical core an ECP removed, are added to the followed density and gradient
std::pair<cubei, std::vector<d4>> topological_cube_analysis(const cube* cub, const std::vector<atom>& atoms, bool debug, bool bcp, double value_floor = 0.0, double grad_epsilon = 1e-12, double assignment_radius = -1.0, double merge_persistence = 5e-3, const std::vector<d3>* seeds = nullptr, const WFN* field_wfn = nullptr, const std::function<double(const d3&)>* core_density = nullptr, const std::function<void(const d3&, d3&)>* core_gradient = nullptr, const density_field* field = nullptr);
//A field's value at a point, gradient filled in: all the streaming analysis asks, so it serves rho and ELI-D alike
using scalar_field = std::function<double(const d3&, d3&)>;
//Newton-Raphson from p onto the nearest critical point; a maximum needs a vanishing gradient and a negative definite
//Hessian (central differences of the analytic gradient). This, not a population threshold, separates a non-nuclear
//attractor from grid debris. On true, p is the converged point.
bool converge_to_maximum(const scalar_field& field, d3& p, double step_limit = 0.5, int max_iterations = 60, double gradient_tolerance = 1e-6);
//Every nucleus plus each non-nuclear attractor found by the CP search that survives converge_to_maximum; no voxel adds to it
std::vector<d4> streaming_density_attractors(const WFN& wavy, const std::vector<critical_point>& critical_points, const std::function<double(const d3&)>* core_density = nullptr, const std::function<void(const d3&, d3&)>* core_gradient = nullptr, bool debug = false);
//A spin member of the ELI family in place of ELI-D: value and gradient at p; aux, when given, receives rho, rho_alpha,
//rho_beta and rho_s * ELI-q_s (WFN::computeELISpinGrad)
using eli_spin_field = std::function<void(const d3&, double&, d3&, double*)>;
//ELI-D maxima without a cube: ascent on computeELIGrad from atom-centred seed shells inside rho >= basin_density_cutoff,
//ends within 0.1 bohr merged, sorted by value, highest first. Core and shattered shells are left to unify_core_basins/unify_shell_basins
std::vector<d4> analytic_eli_maxima(const WFN& wavy, bool debug = false, const eli_spin_field* eli = nullptr);
//A spin field can rise to the rho = basin_density_cutoff isosurface; its attractor is then the field's maximum on that
//surface, reached by sliding along it. analytic_eli_maxima keeps such attractors and the basin integration ends
//trajectories that leave the domain with this.
bool bounded_eli_ascent(const WFN& wavy, const eli_spin_field& eli, d3& p, double& f);
//Surface attractors form a ring or sphere of equal maxima around an axial or spherical open shell, as does a spin field's
//degenerate shell just inside it; they lie beyond unify_shell_basins' reach, linked only over a flat arc. One basin per
//such set; basin_map as for unify_core_basins.
int unify_boundary_basins(std::vector<d4>& maxima, const WFN& wavy, const eli_spin_field& eli, ivec* basin_map = nullptr);
//Radius holding the ELI-D maxima of an atom's core shells, by period (outermost core shell about 0.7 bohr for the first transition row)
double core_shell_radius(const int Z);
//The ELI-D core: core_shell_radius, except for a d-block metal (and K, Ca) whose outer core shell stays outside it
//and keeps its own basins, labelled "shell"
double eli_core_radius(const int Z);
//Electrons the eli_core_radius core holds: the closed shells beneath the valence s,p shell, or beneath a metal's kept outer core shell
int eli_core_electrons(const int Z);
//Every basin whose maximum lies within an atom's eli_core_radius becomes that atom's one core basin (DGrid's ELIDcore);
//returns the number merged away. basin_map, when given, is sized maxima.size() + 1 and maps each input maximum to its
//1-based basin: a streaming walk must still reach every core shell's maximum, while the report shows one core basin per atom.
int unify_core_basins(cubei& basin_cube, std::vector<d4>& maxima, const std::vector<atom>& atoms, ivec* basin_map = nullptr);

//Fold a shattered shell (maxima close together and near-degenerate in value) into one basin each. max_dist is in bohr,
//not voxels: the defect worsens as the grid is refined. With atoms, maxima in a metal's outer core shell (see
//eli_core_radius) merge only within max_dist / 2.
int unify_shell_basins(cubei& basin_cube, std::vector<d4>& maxima, ivec* basin_map = nullptr, double max_dist = 1.2, double rel_tol = 0.05, const std::vector<atom>* atoms = nullptr, const WFN* wavy = nullptr, const eli_spin_field* eli = nullptr);
//Atomic overlap matrices S^b_ij = int_b phi_i phi_j on the populations' points and basin assignment, so a basin's
//trace is its population. One packed lower triangle per basin over the occupied MOs
struct basin_overlaps {
	int nmo = 0;   //MOs in the triangles
	ivec mo_index; //their indices in the WFN, size nmo
	vec2 S;        //S[basin][packed(i,j)], basins in the order of the maxima
	vec outside;   //overlap beyond the density isosurface
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
	vec outside_half;                      //delta(A,outside)/2
	double identity_error = 0.0;           //max |sum_A S^A_ij + S^outside_ij - delta_ij|
};
delocalization_result delocalization_indices(const WFN& wavy, const basin_overlaps& ovl);
void report_delocalization(const WFN& wavy, const basin_overlaps& ovl, const svec& labels, std::ostream& log, const double threshold = 0.01);
//Beta sphere: a radius around an attractor that no ascent leaves, so points inside are assigned without climbing and
//a trajectory entering it stops. Streaming quadrature only, on by default; -no_beta_spheres turns them off.
void beta_spheres_set_enabled(const bool on);
bool beta_spheres_enabled();
//The fraction (0.9 by default, at most 1, NOS_BETA_MARGIN) of the smallest safe radius over 302 sampled directions at
//which the sphere is drawn; it trades the dominant quadrature stage against population accuracy.
double basin_beta_margin();
//Angle-adaptive RK2 step: step_at's step is a floor, a longer one has to be earned from the midpoint gradient.
//Off by default (-adaptive_step) because the longer steps cost population accuracy.
void basin_adaptive_step_set_enabled(const bool on);
bool basin_adaptive_step_enabled();
//The grown step's knobs as the run uses them; enabling the growth re-reads NOS_ADP_CAP, NOS_ADP_GROW, NOS_ADP_KEEP and NOS_ADP_REACH
void basin_adaptive_step_knobs(double &cap, double &grow, double &keep, double &reach);
//Grown-step counters since the last reset: steps, longer steps proposed, turned back by the angle test, dropped again
//because the field stopped rising. Only a run with the growth enabled moves them; fell_back <= proposed is asserted.
void basin_adaptive_step_counters(long long &steps, long long &proposed, long long &turned_back, long long &fell_back);
void basin_adaptive_step_counters_reset();
//-basin_timing: the wall clock of every analysis stage; off by default because the golden files capture this log.
//-basin_step <f>: scales the streaming trajectory step, a fraction of the voxel; 1.0 is the default.
void basin_step_scale_set(const double f);
double basin_step_scale();
//Density trajectories of the last integration that stopped rising on a gradient >= 1e-2 e/bohr^4. Should be zero: such
//a point is a step too long for the local curvature, so the step controller shrinks below its floor rather than handing
//the point to the nearest attractor.
long long basin_stalls_on_a_slope();
void basin_timing_set_enabled(const bool on);
bool basin_timing_enabled();
struct basin_stage_timer {
	std::chrono::steady_clock::time_point t = std::chrono::steady_clock::now();
	//Seconds since the last lap or construction, printed only under -basin_timing.
	void lap(const std::string &what);
};
//Which basin the streaming climbs through a voxel ended in, so a later climb entering a block they all agreed on can stop.
//Fixed capacity, open addressing, no deletion, lock-free: one word per voxel packs the key (3 x 14-bit indices), the
//1-based basin (13 bits), a visit count saturating at 255 and a sticky conflict bit. Results depend on climb order;
//NOS_BASIN_MEMO=0 turns it off, NOS_BASIN_MEMO_MB caps the table (256 MB by default).
struct basin_memo {
	basin_memo(const d3 &lo, const d3 &hi, double voxel, size_t megabytes);
	//The voxel of p, 0 outside the 16382-voxel box centred on lo..hi
	uint64_t key(const d3 &p) const;
	//One clean climb's voxels, consecutive repeats already dropped, all ending in basin label
	void write(const std::vector<uint64_t> &path, int label);
	//The basin of p's voxel when settled: visited at least twice and its 3x3x3 block populated, conflict-free and of one basin; 0 otherwise
	int settled(const d3 &p) const;
	size_t entries() const { return used.load(std::memory_order_relaxed); }
	size_t capacity() const { return table.size(); }
private:
	static constexpr int bits = 14;
	static constexpr uint64_t key_mask = (1ull << 3 * bits) - 1, label_mask = (1ull << 13) - 1;
	static constexpr uint64_t one = 1ull << 55, conflict = 1ull << 63;
	d3 origin;
	double inv;
	std::vector<std::atomic<uint64_t>> table;
	size_t mask, limit;
	std::atomic<size_t> used{ 0 };
	uint64_t find(uint64_t k) const;
	void add(uint64_t k, uint64_t label);
};
//ovl, when given, receives the basin overlap matrices; only for the orbital density (field == nullptr).
//cub or basin_cube null: streaming quadrature, every point climbs the analytic field to one of the maxima; step and
//arrival radius then assume a 0.1 A grid. maximum_basin, when given, is each maximum's 1-based basin (as
//unify_core_basins reports it); the result then has one entry per basin.
//eli (with eli_field): the ELI member whose boundaries are followed instead of ELI-D. spin_pop (with eli): per basin
//{N_alpha, N_beta, int rho_s * ELI-q_s} from the populations' points and weights.
vec integrate_basins_on_atomic_grids(const cube* cub, const cubei* basin_cube, const std::vector<d4>& maxima, const WFN& wavy, const int accuracy, const bool eli_field, vec& volumes, double& outside, const std::function<double(const d3&)>* core_density = nullptr, const std::function<void(const d3&, d3&)>* core_gradient = nullptr, const int grid_boost = 1, const density_field* field = nullptr, basin_overlaps* ovl = nullptr, const ivec* maximum_basin = nullptr, const eli_spin_field* eli = nullptr, vec2* spin_pop = nullptr);
std::vector<critical_point_seed> find_cube_critical_point_seeds(const cube* cub, bool debug, double value_floor = -1.0, double gradient_epsilon = -1.0);
std::vector<critical_point> refine_cube_critical_points(const cube* cub, const WFN& wavy, const std::vector<critical_point_seed>& seeds, bool debug, double value_floor = -1.0, double gradient_tolerance = 1e-8, double step_tolerance = 1e-6, int max_iterations = 32);
std::vector<critical_point> analyze_cube_critical_points(const cube* cub, const WFN& wavy, bool debug, double value_floor = -1.0, double gradient_epsilon = -1.0, double gradient_tolerance = 1e-8, double step_tolerance = 1e-6, int max_iterations = 32);
vec integrate_values_in_basins(const cube *cub, const cubei *basin_cube, svec &basin_label, bool debug);
svec assign_labels_to_basins(const std::vector<d4> &Maxima, const std::vector<atom> &atoms, bool debug, int type_switch = 0);

#include "wfn_class.h"
