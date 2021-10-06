#include <cmath>
#include <limits>
#include <fstream>

#include "AtomGrid.h"
#include "convenience.h"
#include "sphere_lebedev_rule.h"

#ifdef _WIN32
    #define NOMIMAX
    #include <algorithm>
    #include <io.h>
#endif
//Constants for later use
const int hardness = 3;
const double cutoff = 1.0e-20;
const double PI = 3.14159265358979323846;
const double C0 = 4.0 * std::sqrt(2.0) * PI;
const double PI2 = pow(PI, 2);
const double TG32 = tgamma(3.0 / 2.0);
const int max_LT = 33;
const int MAG = 5810;

//                       3,     5     7,    9,    11,   13,   15,   17
//                      19,    21
int lebedev_table[33] = {6,    14,   26,   38,   50,   74,   86,   110,
                         146,  170,  194,  230,  266,  302,  350,  434,
                         590,  770,  974,  1202, 1454, 1730, 2030, 2354,
                         2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810};

constexpr int get_closest_num_angular(int n) {
    int m = 0;

    for (int i = 0; i < max_LT; i++) {
        m = lebedev_table[i];
        if (m >= n)
            return m;
    }

    std::cout << "Input n too high";
    return -1;
};

constexpr int get_angular_order(int n) {
    for (int i = 0; i < max_LT; i++) {
        if (lebedev_table[i] == n)
            return i;
    }

    std::cout << "Not found in get_angular_offset";
    return -1;
};

AtomGrid::AtomGrid(const double radial_precision,
                   const int min_num_angular_points,
                   const int max_num_angular_points,
                   const int proton_charge,
                   const double alpha_max,
                   const int max_l_quantum_number,
                   const double alpha_min[])
{
    int min_num_angular_points_closest = get_closest_num_angular(min_num_angular_points);
    int max_num_angular_points_closest = get_closest_num_angular(max_num_angular_points);

    double *angular_x = new double[max_LT * MAG];
    double *angular_y = new double[max_LT * MAG];
    double *angular_z = new double[max_LT * MAG];
    double *angular_w = new double[max_LT * MAG];

    for (int i = get_angular_order(min_num_angular_points_closest); i < get_angular_order(max_num_angular_points_closest) + 1; i++) {
        int angular_off = i * MAG;
        ld_by_order(lebedev_table[i],
                    &angular_x[angular_off],
                    &angular_y[angular_off],
                    &angular_z[angular_off],
                    &angular_w[angular_off]);
    }

    // obtain radial parameters
    double r_inner = get_r_inner(radial_precision, alpha_max * 2.0);
	double h = (std::numeric_limits<double>::max)();
    double r_outer = 0.0;
    for (int l = 0; l <= max_l_quantum_number; l++) {
        if (alpha_min[l] > 0.0) {
            r_outer = std::max(r_outer,
                         get_r_outer(radial_precision,
                                     alpha_min[l],
                                     l,
                                     4.0 * bragg_angstrom[proton_charge])
                      );
            h = std::min(h, get_h(radial_precision, l, 0.1 * (r_outer - r_inner)));
        }
    }

    num_radial_grid_points_ = 0;

    double rb = bragg_angstrom[proton_charge] / (5.0 * 0.529177249);
    double c = r_inner / (exp(h) - 1.0);
    int num_radial = int(log(1.0 + (r_outer / c)) / h);
    for (int irad = 0; irad < num_radial; irad++) {
        double radial_r = c * (exp(double(irad + 1.0) * h) - 1.0);
        double radial_w = (radial_r + c) * radial_r * radial_r * h;

        radial_atom_grid_r_bohr_.push_back(radial_r);
        radial_atom_grid_w_.push_back(radial_w);
        num_radial_grid_points_++;

        int num_angular = max_num_angular_points_closest;
        if (radial_r < rb) {
            num_angular = static_cast<int>(max_num_angular_points_closest * (radial_r / rb));
            num_angular = get_closest_num_angular(num_angular);
            if (num_angular < min_num_angular_points_closest)
                num_angular = min_num_angular_points_closest;
        }

        int angular_off = get_angular_order(num_angular) * MAG;

        for (int iang = 0; iang < num_angular; iang++) {
            atom_grid_x_bohr_.push_back(angular_x[angular_off + iang] * radial_r);
            atom_grid_y_bohr_.push_back(angular_y[angular_off + iang] * radial_r);
            atom_grid_z_bohr_.push_back(angular_z[angular_off + iang] * radial_r);

            atom_grid_w_.push_back(4.0 * PI * angular_w[angular_off + iang] * radial_w);
        }
    }

    delete[] angular_x;
    delete[] angular_y;
    delete[] angular_z;
    delete[] angular_w;
}

AtomGrid::AtomGrid(const double radial_precision,
                   const int min_num_angular_points,
                   const int max_num_angular_points,
                   const int proton_charge,
                   const double alpha_max,
                   const int max_l_quantum_number,
                   const double alpha_min[],
                   const bool debug,
                   std::ofstream& file,
                   const bool slater)
 {
    using namespace std;
     int min_num_angular_points_closest =
         get_closest_num_angular(min_num_angular_points);
     int max_num_angular_points_closest =
         get_closest_num_angular(max_num_angular_points);

     double* angular_x = new double[max_LT * MAG];
     double* angular_y = new double[max_LT * MAG];
     double* angular_z = new double[max_LT * MAG];
     double* angular_w = new double[max_LT * MAG];
     int angular_off;

     for (int i = get_angular_order(min_num_angular_points_closest); i < get_angular_order(max_num_angular_points_closest) + 1; i++) {
         angular_off = i * MAG;
         ld_by_order(lebedev_table[i],
             &angular_x[angular_off],
             &angular_y[angular_off],
             &angular_z[angular_off],
             &angular_w[angular_off]);
     }

     // radial parameters
     double r_inner = get_r_inner(radial_precision, alpha_max * 2.0); // factor 2.0 to match DIRAC
     if (slater)
         r_inner *= 1E-5;
     double h = (std::numeric_limits<double>::max)();
     double r_outer = 0.0;

     for (int l = 0; l <= max_l_quantum_number; l++) {
         if (alpha_min[l] > 0.0) {
             //if (debug) 
             //    file << "ATOM GRID: " 
             //       << "l= " << l 
             //       << " r_inner: " << r_inner 
             //       << " alpha_min: " << alpha_min[l] << endl;
             r_outer = std::max(r_outer, 
                get_r_outer(radial_precision, 
                        alpha_min[l], 
                        l, 
                        4.0 * bragg_angstrom[proton_charge])
                        );
             h = std::min(h, get_h(radial_precision, 
                                    l, 
                                    0.1 * (r_outer - r_inner))
                          );
         }
     }

     if(debug) 
         file << "ATOM GRID: " 
            << "r_inner: " << r_inner 
            << " h: " << h 
            << " r_outer: " << r_outer << endl;

     num_radial_grid_points_ = 0;

     double rb = bragg_angstrom[proton_charge] / (5.0 * 0.529177249);
     double c = r_inner / (exp(h) - 1.0);
     int num_radial = int(log(1.0 + (r_outer / c)) / h);

     if (debug) 
         file << "ATOM GRID: " 
            << "rb: " << rb 
            << " c: " << c 
            << " num_radial: " << num_radial << endl;

     for (int irad = 0; irad < num_radial; irad++) {
         double radial_r = c * (exp(double(irad + 1.0) * h) - 1.0);
         double radial_w;
         radial_w = (radial_r + c) * radial_r * radial_r * h;

         radial_atom_grid_r_bohr_.push_back(radial_r);
         radial_atom_grid_w_.push_back(radial_w);
         num_radial_grid_points_++;

         int num_angular = max_num_angular_points_closest;
         if (radial_r < rb) {
             num_angular = static_cast<int>(max_num_angular_points_closest *
                 (radial_r / rb));
             num_angular = get_closest_num_angular(num_angular);
             if (num_angular < min_num_angular_points_closest)
                 num_angular = min_num_angular_points_closest;
         }

         int angular_off = get_angular_order(num_angular) * MAG;
         int start = atom_grid_x_bohr_.size();
         atom_grid_x_bohr_.resize(start + num_angular);
         atom_grid_y_bohr_.resize(start + num_angular);
         atom_grid_z_bohr_.resize(start + num_angular);
         atom_grid_w_.resize(start + num_angular);
#pragma omp parallel for
         for (int iang = start; iang < start + num_angular; iang++) {
             atom_grid_x_bohr_[iang] = angular_x[angular_off + iang - start] * radial_r;
             atom_grid_y_bohr_[iang] = angular_y[angular_off + iang - start] * radial_r;
             atom_grid_z_bohr_[iang] = angular_z[angular_off + iang - start] * radial_r;

             atom_grid_w_[iang] = 4.0 * PI * angular_w[angular_off + iang - start] * radial_w;
         }
     }

     delete[] angular_x;
     delete[] angular_y;
     delete[] angular_z;
     delete[] angular_w;
 }

AtomGrid::~AtomGrid() {}

int AtomGrid::get_num_grid_points() const { return atom_grid_x_bohr_.size(); }

int AtomGrid::get_num_radial_grid_points() const {
    return num_radial_grid_points_;
}

void AtomGrid::get_grid(const int num_centers,
    const int center_index,
    const double* x_coordinates_bohr,
    const double* y_coordinates_bohr,
    const double* z_coordinates_bohr,
    const int* proton_charges,
    double grid_x_bohr[],
    double grid_y_bohr[],
    double grid_z_bohr[],
    double grid_aw[],
    double grid_mw[]) const
{
    if (num_centers > 1)
#pragma omp parallel for schedule(dynamic)
        for (int ipoint = 0; ipoint < get_num_grid_points(); ipoint++) {
            grid_x_bohr[ipoint] = atom_grid_x_bohr_[ipoint] + x_coordinates_bohr[center_index];
            grid_y_bohr[ipoint] = atom_grid_y_bohr_[ipoint] + y_coordinates_bohr[center_index];
            grid_z_bohr[ipoint] = atom_grid_z_bohr_[ipoint] + z_coordinates_bohr[center_index];
            grid_mw[ipoint] = atom_grid_w_[ipoint] * get_becke_w(num_centers,
                proton_charges,
                x_coordinates_bohr,
                y_coordinates_bohr,
                z_coordinates_bohr,
                center_index,
                grid_x_bohr[ipoint],
                grid_y_bohr[ipoint],
                grid_z_bohr[ipoint]);
            grid_aw[ipoint] = atom_grid_w_[ipoint];
        }
    else
#pragma omp parallel for schedule(dynamic)
        for (int ipoint = 0; ipoint < get_num_grid_points(); ipoint++) {
            grid_x_bohr[ipoint] = atom_grid_x_bohr_[ipoint] + x_coordinates_bohr[center_index];
            grid_y_bohr[ipoint] = atom_grid_y_bohr_[ipoint] + y_coordinates_bohr[center_index];
            grid_z_bohr[ipoint] = atom_grid_z_bohr_[ipoint] + z_coordinates_bohr[center_index];
            grid_mw[ipoint] = atom_grid_w_[ipoint];
            grid_aw[ipoint] = atom_grid_w_[ipoint];
        }
}

void AtomGrid::get_grid(const int num_centers,
                        const int center_index,
                        const double* x_coordinates_bohr,
                        const double* y_coordinates_bohr,
                        const double* z_coordinates_bohr,
                        const int* proton_charges,
                        double grid_x_bohr[],
                        double grid_y_bohr[],
                        double grid_z_bohr[],
                        double grid_w[]) const
{
    if (num_centers > 1)
        for (size_t ipoint = 0; ipoint < get_num_grid_points(); ipoint++) {
            grid_x_bohr[ipoint] = atom_grid_x_bohr_[ipoint] + x_coordinates_bohr[center_index];
            grid_y_bohr[ipoint] = atom_grid_y_bohr_[ipoint] + y_coordinates_bohr[center_index];
            grid_z_bohr[ipoint] = atom_grid_z_bohr_[ipoint] + z_coordinates_bohr[center_index];
            grid_w[ipoint] = atom_grid_w_[ipoint] * get_becke_w(num_centers,
                proton_charges,
                x_coordinates_bohr,
                y_coordinates_bohr,
                z_coordinates_bohr,
                center_index,
                grid_x_bohr[ipoint],
                grid_y_bohr[ipoint],
                grid_z_bohr[ipoint]);
        }
    else
        for (size_t ipoint = 0; ipoint < get_num_grid_points(); ipoint++) {
            grid_x_bohr[ipoint] = atom_grid_x_bohr_[ipoint] + x_coordinates_bohr[center_index];
            grid_y_bohr[ipoint] = atom_grid_y_bohr_[ipoint] + y_coordinates_bohr[center_index];
            grid_z_bohr[ipoint] = atom_grid_z_bohr_[ipoint] + z_coordinates_bohr[center_index];
            grid_w[ipoint] = atom_grid_w_[ipoint];
        }
}

void AtomGrid::get_grid_omp(const int num_centers,
    const int center_index,
    const double* x_coordinates_bohr,
    const double* y_coordinates_bohr,
    const double* z_coordinates_bohr,
    const int* proton_charges,
    double grid_x_bohr[],
    double grid_y_bohr[],
    double grid_z_bohr[],
    double grid_w[]) const
{
    if (num_centers > 1)
#pragma omp parallel for
        for (int ipoint = 0; ipoint < get_num_grid_points(); ipoint++) {
            grid_x_bohr[ipoint] = atom_grid_x_bohr_[ipoint] + x_coordinates_bohr[center_index];
            grid_y_bohr[ipoint] = atom_grid_y_bohr_[ipoint] + y_coordinates_bohr[center_index];
            grid_z_bohr[ipoint] = atom_grid_z_bohr_[ipoint] + z_coordinates_bohr[center_index];
            grid_w[ipoint] = atom_grid_w_[ipoint] * get_becke_w(num_centers,
                proton_charges,
                x_coordinates_bohr,
                y_coordinates_bohr,
                z_coordinates_bohr,
                center_index,
                grid_x_bohr[ipoint],
                grid_y_bohr[ipoint],
                grid_z_bohr[ipoint]);
        }

    else
#pragma omp parallel for
        for (int ipoint = 0; ipoint < get_num_grid_points(); ipoint++) {
            grid_x_bohr[ipoint] = atom_grid_x_bohr_[ipoint] + x_coordinates_bohr[center_index];
            grid_y_bohr[ipoint] = atom_grid_y_bohr_[ipoint] + y_coordinates_bohr[center_index];
            grid_z_bohr[ipoint] = atom_grid_z_bohr_[ipoint] + z_coordinates_bohr[center_index];
            grid_w[ipoint] = atom_grid_w_[ipoint];
        }
    
}

void AtomGrid::get_atom_grid_omp(
    double grid_x_bohr[],
    double grid_y_bohr[],
    double grid_z_bohr[],
    double grid_w[]) const
{
#pragma omp parallel for
        for (int ipoint = 0; ipoint < get_num_grid_points(); ipoint++) {
            grid_x_bohr[ipoint] = atom_grid_x_bohr_[ipoint];
            grid_y_bohr[ipoint] = atom_grid_y_bohr_[ipoint];
            grid_z_bohr[ipoint] = atom_grid_z_bohr_[ipoint];
            grid_w[ipoint] = atom_grid_w_[ipoint];
        }

}

void AtomGrid::get_radial_grid(double grid_r_bohr[], double grid_w[]) const {
    for (size_t ipoint = 0; ipoint < num_radial_grid_points_; ipoint++) {
        grid_r_bohr[ipoint] = radial_atom_grid_r_bohr_[ipoint];
        grid_w[ipoint] = radial_atom_grid_w_[ipoint];
    }
}

void AtomGrid::get_radial_distances(double grid_r_bohr[]) const {
    for (size_t ipoint = 0; ipoint < num_radial_grid_points_; ipoint++)
        grid_r_bohr[ipoint] = radial_atom_grid_r_bohr_[ipoint];
}

void AtomGrid::get_radial_grid_omp(double grid_r_bohr[], double grid_w[]) const {
#pragma omp parallel for
    for (int ipoint = 0; ipoint < num_radial_grid_points_; ipoint++) {
        grid_r_bohr[ipoint] = radial_atom_grid_r_bohr_[ipoint];
        grid_w[ipoint] = radial_atom_grid_w_[ipoint];
    }
}

void AtomGrid::get_radial_distances_omp(double grid_r_bohr[]) const {
#pragma omp parallel for
    for (int ipoint = 0; ipoint < num_radial_grid_points_; ipoint++)
        grid_r_bohr[ipoint] = radial_atom_grid_r_bohr_[ipoint];
}

// JCP 88, 2547 (1988), eq. 20
inline double f3(const double x)
{
    double f = x;
    for (int i = 0; i < hardness; i++)
        f *= (1.5 - 0.5 * f * f);
    return f;
}

// JCP 88, 2547 (1988)
double get_becke_w(const int num_centers,
    const int proton_charges[],
    const double x_coordinates_bohr[],
    const double y_coordinates_bohr[],
    const double z_coordinates_bohr[],
    const int center_index,
    const double x,
    const double y,
    const double z)
{
    double R_a, R_b;
    double u_ab, a_ab, mu_ab, nu_ab;
    double f, chi;
    double dist_a, dist_b, dist_ab;
    double vx, vy, vz;

    std::vector<double> pa;
    pa.resize(num_centers);

    for (int a = 0; a < num_centers; a++)
        pa[a] = 1.0;

    for (int a = 0; a < num_centers; a++) {
        vx = x_coordinates_bohr[a] - x;
        vy = y_coordinates_bohr[a] - y;
        vz = z_coordinates_bohr[a] - z;
        dist_a = vx * vx + vy * vy + vz * vz;
        dist_a = std::sqrt(dist_a);

        if (a != center_index && dist_a > 15) {
            pa[a] = 0.0;
            continue;
        }

        R_a = bragg_angstrom[proton_charges[a]];

        for (int b = 0; b < a; b++) {
            vx = x_coordinates_bohr[b] - x;
            vy = y_coordinates_bohr[b] - y;
            vz = z_coordinates_bohr[b] - z;
            dist_b = vx * vx + vy * vy + vz * vz;
            dist_b = std::sqrt(dist_b);

            R_b = bragg_angstrom[proton_charges[b]];

            vx = x_coordinates_bohr[b] - x_coordinates_bohr[a];
            vy = y_coordinates_bohr[b] - y_coordinates_bohr[a];
            vz = z_coordinates_bohr[b] - z_coordinates_bohr[a];
            dist_ab = vx * vx + vy * vy + vz * vz;
            dist_ab = std::sqrt(dist_ab);

            // JCP 88, 2547 (1988), eq. 11
            mu_ab = (dist_a - dist_b) / dist_ab;

            if (std::abs(R_a - R_b) > cutoff) {
                chi = R_a / R_b;
                u_ab = (chi - 1) / (chi + 1);
                a_ab = u_ab / (u_ab * u_ab - 1.0);

                // JCP 88, 2547 (1988), eq. A3
                if (a_ab > 0.5)
                    a_ab = 0.5;
                if (a_ab < -0.5)
                    a_ab = -0.5;

                nu_ab = mu_ab + a_ab * (1.0 - mu_ab * mu_ab);
            }
            else
                nu_ab = mu_ab;

            f = f3(nu_ab);

            if (std::abs(1.0 - f) < cutoff)
                pa[a] = 0.0;
            else {
                if (pa[a] > 1E-250 || pa[a] < -1E-250)
                    pa[a] *= 0.5 * (1.0 - f);
                else
                    pa[a] = 0.0;
                if (pa[b] > 1E-250 || pa[b] < -1E-250)
                    pa[b] *= 0.5 * (1.0 + f);
                else
                    pa[b] = 0.0;
            }
        }
    }

    double w = 0.0;
    for (int a = 0; a < num_centers; a++)
        w += pa[a];

    double res = 1.0;
    if (std::abs(w) > cutoff)
        res = pa[center_index] / w;

    pa.clear();

    return res;
}

// TCA 106, 178 (2001), eq. 25
// we evaluate r_inner for s functions
double get_r_inner(const double max_error, const double alpha_inner)
{
    int m = 0;
    double d = 1.9;

    double r = d - std::log(1.0 / max_error);
    r = r * 2.0 / (m + 3.0);
    r = std::exp(r) / (alpha_inner);
    r = std::sqrt(r);

    return r;
}

// TCA 106, 178 (2001), eq. 19
double get_r_outer(const double max_error,
    const double alpha_outer,
    const int l,
    const double guess)
{
    const int m = 2 * l;
    double r = guess;
    double r_old = 1.0e50;
    double c, a, e;
    double step = 0.5;
    double sign, sign_old;
    double f = 1.0e50;

    (f > max_error) ? (sign = 1.0) : (sign = -1.0);

    while (std::abs(r_old - r) > cutoff)
    {
        c = tgamma((m + 3.0) / 2.0);
        a = std::pow(alpha_outer * r * r, (m + 1.0) / 2.0);
        e = std::exp(-alpha_outer * r * r);
        f = c * a * e;

        sign_old = sign;
        (f > max_error) ? (sign = 1.0) : (sign = -1.0);
        if (r < 0.0)
            sign = 1.0;
        if (sign != sign_old)
            step *= 0.1;

        r_old = r;
        r += sign * step;
    }

    return r;
}

// TCA 106, 178 (2001), eqs. 17 and 18
double get_h(const double max_error, const int l, const double guess)
{
    const int m = 2 * l;
    double h = guess;
    double h_old = h * 2;
    double step = 0.1 * guess;
    double sign = -1.0, sign_old, f, pl, rd0, e0;
    const double cm = TG32 / tgamma((m + 3.0) / 2.0);

    while (std::abs(h_old - h) > cutoff)
    {
        e0 = std::exp(-PI2 / (2.0 * h));
        pl = std::pow(PI / h, l);
        rd0 = C0 / h * e0;
        f = cm * pl * rd0;

        sign_old = sign;
        (f > max_error) ? (sign = -1.0) : (sign = 1.0);
        if (h < 0.0)
            sign = 1.0;
        if (sign != sign_old)
            step *= 0.1;

        h_old = h;
        h += sign * step;
        //if (h < 0.007) h = 0.007;
    }

    return h;
}
