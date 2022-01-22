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
             const double alpha_min[]);

    AtomGrid(const double radial_precision,
        const int min_num_angular_points,
        const int max_num_angular_points,
        const int proton_charge,
        const double alpha_max,
        const int max_l_quantum_number,
        const double alpha_min[],
        const bool debug,
        std::ofstream& file);

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
        double grid_mw[]) const;

    void get_grid(const int num_centers,
                  const int center_index,
                  const double* x_coordinates_bohr,
                  const double* y_coordinates_bohr,
                  const double* z_coordinates_bohr,
                  const int* proton_charges,
                  double grid_x_bohr[],
                  double grid_y_bohr[],
                  double grid_z_bohr[],
                  double grid_w[]) const;
    void get_atom_grid_omp(
        double grid_x_bohr[],
        double grid_y_bohr[],
        double grid_z_bohr[],
        double grid_w[]) const;

    void get_grid_omp(const int num_centers,
        const int center_index,
        const double* x_coordinates_bohr,
        const double* y_coordinates_bohr,
        const double* z_coordinates_bohr,
        const int* proton_charges,
        double grid_x_bohr[],
        double grid_y_bohr[],
        double grid_z_bohr[],
        double grid_w[]) const;

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

    std::vector<double> atom_grid_x_bohr_;
    std::vector<double> atom_grid_y_bohr_;
    std::vector<double> atom_grid_z_bohr_;
    std::vector<double> atom_grid_w_;

    std::size_t num_radial_grid_points_;

    std::vector<double> radial_atom_grid_r_bohr_;
    std::vector<double> radial_atom_grid_w_;
};

double get_becke_w(const int num_centers,
    const int proton_charges[],
    const double x_coordinates_bohr[],
    const double y_coordinates_bohr[],
    const double z_coordinates_bohr[],
    const int center_index,
    const double x,
    const double y,
    const double z);

double get_r_inner(const double max_error, const double alpha_inner);

double get_r_outer(const double max_error,
    const double alpha_outer,
    const int l,
    const double guess);

double get_h(const double max_error, const int l, const double guess);
