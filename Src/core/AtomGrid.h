#pragma once

#include <vector>
#include <iosfwd>
#include "density_source.h"

class AtomGrid
{
public:
	AtomGrid(const double radial_precision,
		const int min_num_angular_points,
		const int max_num_angular_points,
		const int proton_charge,
		const double alpha_max,
		const int max_l_quantum_number,
		const double alpha_min[],
		std::ostream &file,
		const double radial_step_scale = 1.0,
		const bool generic_rotation = false);

	~AtomGrid();

	int get_num_grid_points() const;

	int get_num_radial_grid_points() const;

	void get_grid(const int num_centers,
		const int center_index,
		const double *x_coordinates_bohr,
		const double *y_coordinates_bohr,
		const double *z_coordinates_bohr,
		const int *proton_charges,
		double grid_x_bohr[],
		double grid_y_bohr[],
		double grid_z_bohr[],
		double grid_aw[],
		double grid_becke_w[],
		double grid_TFVC_w[],
		const vec &chi) const;

	void get_radial_grid(double grid_r_bohr[], double grid_w[]) const;
	void get_radial_distances(double grid_r_bohr[]) const;
	void get_radial_grid_omp(double grid_r_bohr[], double grid_w[]) const;
	void get_radial_distances_omp(double grid_r_bohr[]) const;

	double *get_gridx_ptr(void) { return atom_grid_x_bohr_.data(); };
	double *get_gridy_ptr(void) { return atom_grid_y_bohr_.data(); };
	double *get_gridz_ptr(void) { return atom_grid_z_bohr_.data(); };
	double *get_gridw_ptr(void) { return atom_grid_w_.data(); };

	double get_gridx(const int &i) { return atom_grid_x_bohr_[i]; };
	double get_gridy(const int &i) { return atom_grid_y_bohr_[i]; };
	double get_gridz(const int &i) { return atom_grid_z_bohr_[i]; };
	const vec &get_gridx() const { return atom_grid_x_bohr_; };
	const vec &get_gridy() const { return atom_grid_y_bohr_; };
	const vec &get_gridz() const { return atom_grid_z_bohr_; };
	//shell s holds points [shell_start[s], shell_start[s + 1]), one radius and one Lebedev rule each
	const ivec &get_shell_start() const { return shell_start_; };

private:

	vec atom_grid_x_bohr_;
	vec atom_grid_y_bohr_;
	vec atom_grid_z_bohr_;
	vec atom_grid_w_;
	ivec shell_start_;

	int num_radial_grid_points_;

	vec radial_atom_grid_r_bohr_;
	vec radial_atom_grid_w_;
};

//Where each kept density point sits in its atom's prototype grid, so the Fourier-Bessel transform takes shells and
//Lebedev rules from the prototype instead of guessing them from the points
struct sf_shell_layout
{
	struct prototype { vec x, y, z; ivec shell_start; }; //points about the nucleus, shells as AtomGrid::get_shell_start
	std::vector<prototype> protos;
	ivec proto; //[density atom] its prototype, -1 none
	ivec2 idx;  //[density atom][kept point] index in that prototype
};

//The line density comes from density when set (a fitted density on an aux-only wfn), else from the orbitals of wfn
vec make_chi(const WFN &wfn, int samples = 50, bool refine = true, bool debug = false, const DensityBatch &density = {});

std::array<double, 2> get_integration_weights(const int &num_centers,
	const int *proton_charges,
	const double *x_coordinates_bohr,
	const double *y_coordinates_bohr,
	const double *z_coordinates_bohr,
	const int &center_index,
	const double &x,
	const double &y,
	const double &z,
	vec &pa_b,
	vec &pa_tv,
	const vec &chi,
	const double *dist_ab = nullptr); //interatomic distances, num_centers x num_centers, or computed per call

const double get_r_inner(const double &max_error, const double &alpha_inner);

double get_r_outer(const double &max_error,
	const double &alpha_outer,
	const int &l,
	const double &guess);

double get_h(const double &max_error, const int &l, const double &guess);

std::vector<std::pair<vec, vec>> make_MBIS_vectors(
	const WFN &wavy,
	const vec3 &grid,
	const ivec &num_grid_points,
	const bool debug = false,
	std::ostream &file = std::cout);

std::vector<std::pair<vec2, vec>> make_EMBIS_tensors(
	const WFN &wavy,
	const vec3 &grid,
	const ivec &num_grid_points,
	const bool debug = false,
	const std::vector<std::pair<vec, vec>> MBIS_vectors = {},
	std::ostream &file = std::cout);