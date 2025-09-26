#pragma once

#include <vector>

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
             std::ostream& file);

    ~AtomGrid();

    int get_num_grid_points() const;

    int get_num_radial_grid_points() const;

    void get_grid(const int num_centers,
                  const int center_index,
                  const double* x_coordinates_bohr,
                  const double* y_coordinates_bohr,
                  const double* z_coordinates_bohr,
                  const int* proton_charges,
                  double grid_x_bohr[],
                  double grid_y_bohr[],
                  double grid_z_bohr[],
                  double grid_aw[],
                  double grid_becke_w[],
                  double grid_TFVC_w[],
                  const WFN& wfn,
                  vec2& chi,
                  bool debug = false) const;

    void get_radial_grid(double grid_r_bohr[], double grid_w[]) const;
    void get_radial_distances(double grid_r_bohr[]) const;
    void get_radial_grid_omp(double grid_r_bohr[], double grid_w[]) const;
    void get_radial_distances_omp(double grid_r_bohr[]) const;

    double* get_gridx_ptr(void) { return atom_grid_x_bohr_.data(); };
    double* get_gridy_ptr(void) { return atom_grid_y_bohr_.data(); };
    double* get_gridz_ptr(void) { return atom_grid_z_bohr_.data(); };
    double* get_gridw_ptr(void) { return atom_grid_w_.data(); };

    double get_gridx(const int& i) { return atom_grid_x_bohr_[i]; };
    double get_gridy(const int& i) { return atom_grid_y_bohr_[i]; };
    double get_gridz(const int& i) { return atom_grid_z_bohr_[i]; };

  private:

    vec atom_grid_x_bohr_;
    vec atom_grid_y_bohr_;
    vec atom_grid_z_bohr_;
    vec atom_grid_w_;

    int num_radial_grid_points_;

    vec radial_atom_grid_r_bohr_;
    vec radial_atom_grid_w_;
};

std::pair<double,double> get_integration_weights(const int& num_centers,
    const int* proton_charges,
    const double* x_coordinates_bohr,
    const double* y_coordinates_bohr,
    const double* z_coordinates_bohr,
    const int& center_index,
    const double& x,
    const double& y,
    const double& z,
    vec& pa_b,
    vec& pa_tv, 
    const vec2& chi);

const double get_r_inner(const double& max_error, const double& alpha_inner);

double get_r_outer(const double& max_error,
    const double& alpha_outer,
    const int& l,
    const double& guess);

double get_h(const double& max_error, const int& l, const double& guess);
