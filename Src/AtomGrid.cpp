#include "pch.h"
#include "convenience.h"
#include "AtomGrid.h"
#include "sphere_lebedev_rule.h"
#include "constants.h"
#include "wfn_class.h"
#include "spherical_density.h"

#ifdef _WIN32
#include <algorithm>
#include <io.h>
#endif

// Helper to describe extrema along a line between two atoms
struct DensityExtremum
{
    enum class Type { Minimum, Maximum } type;
    d3 pos{ 0, 0, 0 };     // Cartesian position (Bohr)
    double density{ 0 };           // Electron density at extremum
    double distA{ 0 };             // Distance to atom A (Bohr)
    double distB{ 0 };             // Distance to atom B (Bohr)
    double t{ 0 };                 // Parameter along AB in [0,1]
};

// Find all local minima and maxima of the electron density along the straight
// line segment between atomA and atomB. The function samples uniformly and
// detects zero-crossings of the first derivative, with optional quadratic
// interpolation refinement around the detected extremum.
static std::vector<DensityExtremum> find_line_density_extrema(
    const WFN& wfn,
    int atomA,
    int atomB,
    int samples = 400,
    bool refine = true)
{
    std::vector<DensityExtremum> result;
    if (atomA < 0 || atomB < 0 || atomA == atomB) return result;

    // Clamp minimum samples
    samples = std::max(samples, 5);

    // Atom positions (Bohr)
    const d3 pos_a = wfn.get_atom_pos((unsigned)atomA);
    const d3 pos_b = wfn.get_atom_pos((unsigned)atomB);

    const d3 dx = { pos_b[0] - pos_a[0],
       pos_b[1] - pos_a[1],
       pos_b[2] - pos_a[2] };

    const double ab_len = array_length(dx);
    if (ab_len <= 0) return result;

    // Sampling t in (0,1); avoid endpoints to reduce numerical issues near nuclei
    const double eps = 1e-6;
    const double t0 = eps, t1 = 1.0 - eps;
    const double dt = (t1 - t0) / (samples - 1);

    vec dens(samples);
    vec ts(samples);
    Thakkar* spherical_temp1 = NULL;
    Thakkar* spherical_temp2 = NULL;
    Spherical_Gaussian_Density* a1 = NULL;
    Spherical_Gaussian_Density* a2 = NULL;
    if (wfn.get_ECP_mode() != 0) {
        spherical_temp1 = new Thakkar(wfn.get_atom_charge(atomA), wfn.get_ECP_mode());
        spherical_temp2 = new Thakkar(wfn.get_atom_charge(atomB), wfn.get_ECP_mode());
        a1 = new Spherical_Gaussian_Density(wfn.get_atom_charge(atomA), wfn.get_ECP_mode());
        a2 = new Spherical_Gaussian_Density(wfn.get_atom_charge(atomB), wfn.get_ECP_mode());
    }
#pragma omp parallel for
    for (int i = 0; i < samples; i++) {
        const double t = t0 + i * dt;
        ts[i] = t;
        const d3 Pos = { pos_a[0] + dx[0] * t,
          pos_a[1] + dx[1] * t,
          pos_a[2] + dx[2] * t };
        dens[i] = wfn.compute_dens(Pos);
        if (wfn.get_ECP_mode() != 0) {
            dens[i] += spherical_temp1->get_core_density(ab_len * t, wfn.get_atom_ECP_electrons(atomA));
            dens[i] += spherical_temp2->get_core_density(ab_len * (1.0 - t), wfn.get_atom_ECP_electrons(atomB));
            dens[i] += a1->get_radial_density(ab_len * t);
            dens[i] += a2->get_radial_density(ab_len * (1.0 - t));
        }
    }

    if (wfn.get_ECP_mode() != 0) {
        delete spherical_temp1;
        delete spherical_temp2;
        delete a1;
        delete a2;
    }

    // Detect extrema via sign changes of the discrete derivative
    for (int i = 1; i < samples - 1; i++) {
        const double d1 = dens[i] - dens[i - 1];
        const double d2 = dens[i + 1] - dens[i];
        if (d1 == 0 && d2 == 0) continue; // flat region

        bool isMax = (d1 > 0 && d2 < 0);
        bool isMin = (d1 < 0 && d2 > 0);
        if (!isMax && !isMin) continue;

        double t_ext = ts[i];
        double f_ext = dens[i];

        // Optional parabolic interpolation using points (i-1, i, i+1)
        if (refine) {
            const double y0 = dens[i - 1];
            const double y1 = dens[i];
            const double y2 = dens[i + 1];
            // Vertex of parabola fitted to equally spaced samples
            const double denom = (y0 - 2.0 * y1 + y2);
            if (std::abs(denom) > 1e-20) {
                const double delta = 0.5 * (y0 - y2) / denom; // in units of dt around center
                // Clamp delta to [-1,1] to remain within the bracket
                const double clamped = std::max(-1.0, std::min(1.0, delta));
                t_ext = ts[i] + clamped * dt;
                const d3 Pos = { pos_a[0] + dx[0] * t_ext,
                  pos_a[1] + dx[1] * t_ext,
                  pos_a[2] + dx[2] * t_ext };
                f_ext = wfn.compute_dens(Pos);
            }
        }

        DensityExtremum e;
        e.type = isMax ? DensityExtremum::Type::Maximum : DensityExtremum::Type::Minimum;
        e.t = t_ext;
        e.pos = { pos_a[0] + dx[0] * t_ext,
          pos_a[1] + dx[1] * t_ext,
          pos_a[2] + dx[2] * t_ext };
        e.density = f_ext;
        e.distA = t_ext * ab_len;
        e.distB = (1.0 - t_ext) * ab_len;
        result.push_back(e);
    }

    return result;
}

AtomGrid::AtomGrid(const double radial_precision,
    const int min_num_angular_points,
    const int max_num_angular_points,
    const int proton_charge,
    const double alpha_max,
    const int max_l_quantum_number,
    const double alpha_min[],
    std::ostream& file)
{
    using namespace std;
    const int min_num_angular_points_closest =
        constants::get_closest_num_angular(min_num_angular_points);
    const int max_num_angular_points_closest =
        constants::get_closest_num_angular(max_num_angular_points);
    err_checkf(min_num_angular_points_closest != -1 && max_num_angular_points_closest != -1, "No valid value for angular number found!", file);

    vec angular_x(constants::max_LT * constants::MAG, 0.0);
    vec angular_y(constants::max_LT * constants::MAG, 0.0);
    vec angular_z(constants::max_LT * constants::MAG, 0.0);
    vec angular_w(constants::max_LT * constants::MAG, 0.0);
    int angular_off;
    lebedev_sphere ls;

    for (int i = constants::get_angular_order(min_num_angular_points_closest); i < constants::get_angular_order(max_num_angular_points_closest) + 1; i++) {
        angular_off = i * constants::MAG;
        ls.ld_by_order(constants::lebedev_table[i],
            angular_x.data() + angular_off,
            angular_y.data() + angular_off,
            angular_z.data() + angular_off,
            angular_w.data() + angular_off);
    }

    // radial parameters
    const double r_inner = get_r_inner(radial_precision, alpha_max * 2.0); // factor 2.0 to match DIRAC
    double h = (std::numeric_limits<double>::max)();
    double r_outer = 0.0;

    for (int l = 0; l <= max_l_quantum_number; l++) {
        if (alpha_min[l] > 0.0) {
            //if (debug) 
            //    file << "ATOM GRID: " 
            //       << "l= " << l 
            //       << " r_inner: " << r_inner 
            //       << " alpha_min: " << alpha_min[l] << endl;
            r_outer = (((r_outer) > (get_r_outer(radial_precision, alpha_min[l], l, 4.0 * constants::bragg_angstrom[proton_charge]))) ? (r_outer) : (get_r_outer(radial_precision, alpha_min[l], l, 4.0 * constants::bragg_angstrom[proton_charge])));
            h = (((h) < (get_h(radial_precision, l, 0.1 * (r_outer - r_inner)))) ? (h) : (get_h(radial_precision, l, 0.1 * (r_outer - r_inner))));
        }
    }

    //if (debug)
    //  file << "ATOM GRID: "
    //  << "r_inner: " << r_inner
    //  << " h: " << h
    //  << " r_outer: " << r_outer << endl;

    num_radial_grid_points_ = 0;

    const double rb = constants::bragg_angstrom[proton_charge] / (5.0E10 * constants::a0);
    const double c = r_inner / (exp(h) - 1.0);
    const int num_radial = int(log(1.0 + (r_outer / c)) / h);
    radial_atom_grid_r_bohr_.reserve(num_radial);
    radial_atom_grid_w_.reserve(num_radial);

    //if (debug)
    //  file << "ATOM GRID: "
    //  << "rb: " << rb
    //  << " c: " << c
    //  << " num_radial: " << num_radial << endl;

    for (int irad = 0; irad < num_radial; irad++) {
        double radial_r = c * (exp(static_cast<double>(irad + 1.0) * h) - 1.0);
        double radial_w = (radial_r + c) * radial_r * radial_r * h;

        radial_atom_grid_r_bohr_.emplace_back(radial_r);
        radial_atom_grid_w_.emplace_back(radial_w);
        num_radial_grid_points_++;

        int num_angular = max_num_angular_points_closest;
        if (radial_r < rb) {
            num_angular = static_cast<int>(max_num_angular_points_closest *
                (radial_r / rb));
            num_angular = constants::get_closest_num_angular(num_angular);
            err_checkf(num_angular != -1, "No valid value for angular number found!", file);
            if (num_angular < min_num_angular_points_closest)
                num_angular = min_num_angular_points_closest;
        }

        angular_off = constants::get_angular_order(num_angular) * constants::MAG;
        err_checkf(angular_off != -constants::MAG, "Invalid angular order!", file);
        const int start = (int)atom_grid_x_bohr_.size();
        angular_off -= start;
        const int size = start + num_angular;
        int p = 0;
        atom_grid_x_bohr_.resize(size);
        atom_grid_y_bohr_.resize(size);
        atom_grid_z_bohr_.resize(size);
        atom_grid_w_.resize(size);
#pragma omp parallel for private(p)
        for (int iang = start; iang < size; iang++) {
            p = angular_off + iang;
            atom_grid_x_bohr_[iang] = angular_x[p] * radial_r;
            atom_grid_y_bohr_[iang] = angular_y[p] * radial_r;
            atom_grid_z_bohr_[iang] = angular_z[p] * radial_r;

            atom_grid_w_[iang] = constants::FOUR_PI * angular_w[p] * radial_w;
        }
    }
}

AtomGrid::~AtomGrid() {}

int AtomGrid::get_num_grid_points() const { return (int)atom_grid_x_bohr_.size(); }

int AtomGrid::get_num_radial_grid_points() const { return num_radial_grid_points_; }

vec make_chi(const WFN& wfn, int samples, bool refine, bool debug) {
    const int ncen = wfn.get_ncen();
    vec chi(ncen * ncen, 0.0);
    std::vector<std::vector<bool>> neighbours(ncen, bvec(ncen, true));
    double rijx2, rijy2, rijz2, xdist, disth;
    for (int a = 0; a < wfn.get_ncen(); a++) {
        for (int b = a + 1; b < wfn.get_ncen(); b++) {
            // Evaluating midpoint between a pair of atoms
            rijx2 = (wfn.get_atom_coordinate(b, 0) + wfn.get_atom_coordinate(a, 0)) / 2.0;
            rijy2 = (wfn.get_atom_coordinate(b, 1) + wfn.get_atom_coordinate(a, 1)) / 2.0;
            rijz2 = (wfn.get_atom_coordinate(b, 2) + wfn.get_atom_coordinate(a, 2)) / 2.0;
            disth = vec_length({ wfn.get_atom_coordinate(b, 0) - wfn.get_atom_coordinate(a, 0),
                                 wfn.get_atom_coordinate(b, 1) - wfn.get_atom_coordinate(a, 1),
                                 wfn.get_atom_coordinate(b, 2) - wfn.get_atom_coordinate(a, 2) }) / 2.0; // Half inter-atomic distance
            //std::cout << a << " " << b << " disth: " << disth << std::endl;
            // Check if any other atom is considered "neighbour"
            for (int ii = 0; ii < wfn.get_ncen(); ii++) {
                if (ii != a && ii != b) {
                    // Calculate distance between ii and the midpoint of atom pair AB
                    xdist = vec_length({ rijx2 - wfn.get_atom_coordinate(ii, 0),
                                         rijy2 - wfn.get_atom_coordinate(ii, 1),
                                         rijz2 - wfn.get_atom_coordinate(ii, 2) });
                    //std::cout << "ii: " << ii << " xdist: " << xdist;
                    // Condition to satisfy for neightbouring: xdist lower than disth
                    if (xdist < disth) {
                        // Set chi to 1.0 temporarily
                        //std::cout << " Entering" << std::endl;
                        neighbours[a][b] = false;
                        neighbours[b][a] = false;
                        break;
                    }
                    //std::cout << std::endl;
                }
            }
        }
    }
    for (int a = 0; a < wfn.get_ncen(); a++) {
        rijx2 = wfn.get_atom_coordinate(a, 0); //abuse old variable
        rijy2 = wfn.get_atom_coordinate(a, 1);
        rijz2 = wfn.get_atom_coordinate(a, 2);
#pragma omp parallel for
        for (int b = a + 1; b < wfn.get_ncen(); b++) {
            if (neighbours[a][b]) {
                if (vec_length({ wfn.get_atom_coordinate(b, 0) - rijx2,
                                 wfn.get_atom_coordinate(b, 1) - rijy2,
                                 wfn.get_atom_coordinate(b, 2) - rijz2 }) > constants::far_away) {
                    // If they are neighbours but far apart we will consider them equal radius
                    chi[a * ncen + b] = 1.0;
                    chi[b * ncen + a] = 1.0;
                }
                else {
                    // If they are neighbours and close enough we do the line search for topology evaluation!
                    auto extrema = find_line_density_extrema(wfn, a, b, samples, refine);
                    size_t use_extr = 0;
                    if (wfn.get_atom_ECP_electrons(a) != 0 || wfn.get_atom_ECP_electrons(b) != 0) {
                        double closeness = 1.0;
                        for (int i = 0; i < extrema.size(); i++) {
                            if (abs(extrema[i].t - 0.5) < closeness && extrema[i].type == DensityExtremum::Type::Minimum) {
                                use_extr = i;
                                closeness = abs(extrema[i].t - 0.5);
                            }
                        }
                    }
                    else {
                        if (extrema.size() == 1 && extrema[0].type == DensityExtremum::Type::Minimum) {
                            use_extr = 0;
                        }
                        else if (extrema.size() % 2 != 0) {
                            int maxmin = 0;
                            for (auto ex : extrema) {
                                if (ex.type == DensityExtremum::Type::Maximum)
                                    maxmin++;
                                else
                                    maxmin--;
                            }
                            if (maxmin > 0) {
                                double closeness = 1.0;
                                for (int i = 0; i < extrema.size(); i++) {
                                    if (abs(extrema[i].t - 0.5) < closeness && extrema[i].type == DensityExtremum::Type::Minimum) {
                                        use_extr = i;
                                        closeness = extrema[i].t - 0.5;
                                    }
                                }
                            }
                            else {
                                if (extrema[1].type == DensityExtremum::Type::Minimum)
                                    use_extr = 1;
                                else if (extrema[1].type == DensityExtremum::Type::Maximum) {
                                    if (debug)
                                        std::cout << "Strange case with 3 extrema but middle is not a minimum. Using Maximum." << std::endl;
                                    use_extr = 1;
                                }
                            }
                        }
                        else if (extrema.size() == 2) {
                            //"Hydrogen case"
                            if (extrema[0].type == DensityExtremum::Type::Minimum && extrema[1].type == DensityExtremum::Type::Maximum) {
                                use_extr = 0;
                            }
                            else if (extrema[0].type == DensityExtremum::Type::Maximum && extrema[1].type == DensityExtremum::Type::Minimum) {
                                use_extr = 1;
                            }
                            else {
                                std::cout << "WARNING: Unexpected types of extrema found between atoms " << a << " and " << b << ". Setting chi to 1.0." << std::endl;
                                chi[a * ncen + b] = 1.0;
                                chi[b * ncen + a] = 1.0;
                                continue;
                            }
                        }
                        else {
                            std::cout << "WARNING: Unexpected number of extrema (" << extrema.size() << ") found between atoms " << a << " and " << b << ". Setting chi to 1.0." << std::endl;
                            chi[a * ncen + b] = 1.0;
                            chi[b * ncen + a] = 1.0;
                            continue;
                        }
                    }
                    // And save the "chi" in the same matrix
                    if (extrema[use_extr].distB == 0.0) {
                        std::cout << "WARNING: Division by zero detected for chi between atoms " << a << " and " << b << ". Setting chi to 1.0." << std::endl;
                        chi[a * ncen + b] = 1.0;
                        chi[b * ncen + a] = 1.0;
                    }
                    else {
                        chi[a * ncen + b] = extrema[use_extr].distA / extrema[use_extr].distB;
                        chi[b * ncen + a] = 1.0 / chi[a * ncen + b];
                    }
                    // For if one have a very polarized bond (bcp super displaced towards one of the atoms
                    double ri = extrema[use_extr].distA / (extrema[use_extr].distA + extrema[use_extr].distB);
                    if (ri >= 0.90 || ri <= 0.10) {
                        std::cout << "WARNING: Large bond polarization detected between atoms " << a << " and " << b
                            << " (ri = " << ri << ", ri >= 0.90 or ri <= 0.10). Please verify TFVC results." << std::endl;
                    }
                    // Some printing on the ratio between pair of atoms, for if strange results (or against chemical intuition) are obtained one can see if it is TFVC problem
                    // This is helpful to see if the "bcp" evaluation is broken, or one face a nasty scenario that is not yet included in the code...
                    if (debug)
                        std::cout << "TFVC -- Atom pair: " << a << ", " << b << ", Ratio: " << chi[a * ncen + b] << std::endl;
                }
            }
            else {
                chi[a * ncen + b] = 1.0;
                chi[b * ncen + a] = 1.0;
            }
        }
    }

    return chi;
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
    double grid_becke_w[],
    double grid_TFVC_w[],
    const WFN& wfn,
    vec& chi,
    bool debug) const
{

    if (num_centers > 1) {
        if (chi.size() == 0)
            chi = make_chi(wfn, 40, true, debug);
        const int np = get_num_grid_points();
#pragma omp parallel
        {
            vec pa_b(num_centers);
            vec pa_tv(num_centers);
            std::array<double, 2> result_weights;
            double temp;
#pragma omp for
            for (int ipoint = 0; ipoint < np; ipoint++) {
                grid_x_bohr[ipoint] = atom_grid_x_bohr_[ipoint] + x_coordinates_bohr[center_index];
                grid_y_bohr[ipoint] = atom_grid_y_bohr_[ipoint] + y_coordinates_bohr[center_index];
                grid_z_bohr[ipoint] = atom_grid_z_bohr_[ipoint] + z_coordinates_bohr[center_index];
                temp = atom_grid_w_[ipoint];
                result_weights = get_integration_weights(
                    num_centers,
                    proton_charges,
                    x_coordinates_bohr,
                    y_coordinates_bohr,
                    z_coordinates_bohr,
                    center_index,
                    grid_x_bohr[ipoint],
                    grid_y_bohr[ipoint],
                    grid_z_bohr[ipoint],
                    pa_b,
                    pa_tv,
                    chi
                );
                grid_becke_w[ipoint] = temp * result_weights[0];
                grid_TFVC_w[ipoint] = temp * result_weights[1];
                grid_aw[ipoint] = temp;
            }
        }
    }
    else
#pragma omp parallel for schedule(dynamic)
        for (int ipoint = 0; ipoint < get_num_grid_points(); ipoint++) {
            grid_x_bohr[ipoint] = atom_grid_x_bohr_[ipoint] + x_coordinates_bohr[center_index];
            grid_y_bohr[ipoint] = atom_grid_y_bohr_[ipoint] + y_coordinates_bohr[center_index];
            grid_z_bohr[ipoint] = atom_grid_z_bohr_[ipoint] + z_coordinates_bohr[center_index];
            grid_becke_w[ipoint] = atom_grid_w_[ipoint];
            grid_TFVC_w[ipoint] = atom_grid_w_[ipoint];
            grid_aw[ipoint] = atom_grid_w_[ipoint];
        }
}

void AtomGrid::get_radial_grid(double grid_r_bohr[], double grid_w[]) const
{
    for (size_t ipoint = 0; ipoint < num_radial_grid_points_; ipoint++) {
        grid_r_bohr[ipoint] = radial_atom_grid_r_bohr_[ipoint];
        grid_w[ipoint] = radial_atom_grid_w_[ipoint];
    }
}

void AtomGrid::get_radial_distances(double grid_r_bohr[]) const
{
    for (size_t ipoint = 0; ipoint < num_radial_grid_points_; ipoint++)
        grid_r_bohr[ipoint] = radial_atom_grid_r_bohr_[ipoint];
}

void AtomGrid::get_radial_grid_omp(double grid_r_bohr[], double grid_w[]) const
{
#pragma omp parallel for
    for (int ipoint = 0; ipoint < num_radial_grid_points_; ipoint++) {
        grid_r_bohr[ipoint] = radial_atom_grid_r_bohr_[ipoint];
        grid_w[ipoint] = radial_atom_grid_w_[ipoint];
    }
}

void AtomGrid::get_radial_distances_omp(double grid_r_bohr[]) const
{
#pragma omp parallel for
    for (int ipoint = 0; ipoint < num_radial_grid_points_; ipoint++)
        grid_r_bohr[ipoint] = radial_atom_grid_r_bohr_[ipoint];
}

// JCP 88, 2547 (1988), eq. 20
constexpr double f3(const double& x)
{
    double f = x;
    f *= (1.5 - 0.5 * f * f); // First iteration
    f *= (1.5 - 0.5 * f * f); // Second iteration
    f *= (1.5 - 0.5 * f * f); // Third iteration
    return f;
}

// JCP 88, 2547 (1988), eq. 20
constexpr double f4(const double& x)
{
    double f = x;
    f *= (1.5 - 0.5 * f * f); // First iteration
    f *= (1.5 - 0.5 * f * f); // Second iteration
    f *= (1.5 - 0.5 * f * f); // Third iteration
    f *= (1.5 - 0.5 * f * f); // Fourth iteration
    return f;
}

constexpr double f(const double& x)
{
    double f = x;
    for (int i = 0; i < constants::hardness; i++)
        f *= (1.5 - 0.5 * f * f);
    return f;
}

// JCP 139, 071103 (2013) for TFVC
// JCP 88, 2547 (1988) for Becke
std::array<double, 2> get_integration_weights(const int& num_centers,
    const int* proton_charges,
    const double* x_coordinates_bohr,
    const double* y_coordinates_bohr,
    const double* z_coordinates_bohr,
    const int& center_index,
    const double& x,
    const double& y,
    const double& z,
    std::vector<double>& pa_b,
    std::vector<double>& pa_tv,
    const vec& chi)
{
    double mu_ab, nu_ab, f, dist_ab;
    double dist_a, dist_b;
    double vx, vy, vz;
    double R_a, R_b, chi_becke, u_ab, chi_mod;
    const double* chi_off;

    for (int a = 0; a < num_centers; a++) {
        pa_b[a] = 1.0;
        pa_tv[a] = 1.0;
    }

    for (int a = 0; a < num_centers; a++) {
        vx = x_coordinates_bohr[a] - x;
        vy = y_coordinates_bohr[a] - y;
        vz = z_coordinates_bohr[a] - z;
        dist_a = std::sqrt(vx * vx + vy * vy + vz * vz);

        double& pa_b_a = pa_b[a];
        double& pa_tv_a = pa_tv[a];

        if (dist_a > constants::far_away) {
            pa_b_a = 0.0;
            pa_tv_a = 0.0;
            continue;
        }

        R_a = constants::bragg_angstrom[proton_charges[a]];
        chi_off = chi.data() + a * num_centers;

        for (int b = a + 1; b < num_centers; b++) {
            double& pa_b_b = pa_b[b];
            double& pa_tv_b = pa_tv[b];

            vx = x_coordinates_bohr[b] - x_coordinates_bohr[a];
            vy = y_coordinates_bohr[b] - y_coordinates_bohr[a];
            vz = z_coordinates_bohr[b] - z_coordinates_bohr[a];
            dist_ab = std::sqrt(vx * vx + vy * vy + vz * vz);

            vx = x_coordinates_bohr[b] - x;
            vy = y_coordinates_bohr[b] - y;
            vz = z_coordinates_bohr[b] - z;
            dist_b = std::sqrt(vx * vx + vy * vy + vz * vz);
            R_b = constants::bragg_angstrom[proton_charges[b]];

            // JCP 139, 071103 (2013), eq. 7
            // JCP 88, 2547 (1988), eq. 11
            mu_ab = (dist_a - dist_b) / dist_ab;

            chi_mod = *(chi_off + b) * (1.0 - mu_ab);

            nu_ab = 1.0 + mu_ab;

            nu_ab = (nu_ab - chi_mod) / (nu_ab + chi_mod);

            f = 1.0 - f4(nu_ab);

            if (std::abs(f) < constants::cutoff)
                pa_tv_a = 0.0;
            else {
                //Reduce numerical jittering
                if (pa_tv_a > 1E-250 || pa_tv_a < -1E-250)
                    pa_tv_a *= 0.5 * f;
                else
                    pa_tv_a = 0.0;
                if (pa_tv_b > 1E-250 || pa_tv_b < -1E-250)
                    pa_tv_b *= 0.5 * (2.0 - f);
                else
                    pa_tv_b = 0.0;
            }


            if (std::abs(R_a - R_b) > constants::cutoff) {
                chi_becke = R_a / R_b;
                u_ab = (chi_becke - 1.0) / (chi_becke + 1.0);
                u_ab = u_ab / (u_ab * u_ab - 1.0);

                // JCP 88, 2547 (1988), eq. A3
                if (u_ab > 0.5)
                    u_ab = 0.5;
                else if (u_ab < -0.5)
                    u_ab = -0.5;

                nu_ab = mu_ab + u_ab * (1.0 - mu_ab * mu_ab);
            }
            else
                nu_ab = mu_ab;

            f = 1.0 - f3(nu_ab);

            if (std::abs(f) < constants::cutoff)
                pa_b_a = 0.0;
            else {
                if (pa_b_a > 1E-250 || pa_b_a < -1E-250)
                    pa_b_a *= 0.5 * f;
                else
                    pa_b_a = 0.0;
                if (pa_b_b > 1E-250 || pa_b_b < -1E-250)
                    pa_b_b *= 0.5 * (2.0 - f);
                else
                    pa_b_b = 0.0;
            }
        }
    }

    double w_becke = 0.0, w_tfvc = 0.0;
    for (int a = 0; a < num_centers; a++) {
        w_becke += pa_b[a];
        w_tfvc += pa_tv[a];
    }

    return { std::abs(w_becke) > constants::cutoff ? pa_b[center_index] / w_becke : 1.0, std::abs(w_tfvc) > constants::cutoff ? pa_tv[center_index] / w_tfvc : 1.0 };
}

// TCA 106, 178 (2001), eq. 25
// we evaluate r_inner for s functions
const double get_r_inner(const double& max_error, const double& alpha_inner)
{
    double d = 1.9;

    double r = (d - log(1.0 / max_error)) * 2. / 3.;
    r = exp(r) / (alpha_inner);
    r = std::sqrt(r);

    return r;
}

std::pair<vec, vec> get_shsig_shpop(const int& atom_type) {
    vec shsig(constants::MBIS_function[atom_type], 0.0);
    vec shpop(constants::MBIS_function[atom_type], 0.0);

    switch (constants::MBIS_function[atom_type]) {
    case 0:
        err_not_impl_f("MBIS does not work with ghost atoms. Please dont do that to me!", std::cout);
        break;
    case 1:
        shsig[0] = 1.0 / (2 * atom_type);
        shpop[0] = atom_type;
        break;
    case 2:
        shsig[0] = 1.0 / (2 * atom_type);
        shpop[0] = 2;
        shsig[1] = 0.5;
        shpop[1] = atom_type - 2;
        break;
    case 3:
        shsig[0] = 1.0 / (2 * atom_type);
        shpop[0] = 2;
        shsig[1] = 1.0 / (2 * pow(atom_type, 0.5));
        shpop[1] = 8;
        shsig[2] = 0.5;
        shpop[2] = atom_type - 10;
        break;
    case 4:
        shsig[0] = 1.0 / (2 * atom_type);
        shpop[0] = 2;
        shsig[1] = 1.0 / (2 * pow(atom_type, constants::c_23));
        shpop[1] = 8;
        shsig[2] = 1.0 / (2 * pow(atom_type, constants::c_13));
        shpop[2] = 8;
        shsig[3] = 0.5;
        shpop[3] = atom_type - 18;
        break;
    case 5:
        shsig[0] = 1.0 / (2 * atom_type);
        shpop[0] = 2;
        shsig[1] = 1.0 / (2 * pow(atom_type, 0.25));
        shpop[1] = 8;
        shsig[2] = 1.0 / (2 * pow(atom_type, 0.5));
        shpop[2] = 8;
        shsig[3] = 1.0 / (2 * pow(atom_type, 0.75));
        shpop[3] = 18;
        shsig[4] = 0.5;
        shpop[4] = atom_type - 36;
        break;
    case 6:
        shsig[0] = 1.0 / (2 * atom_type);
        shpop[0] = 2;
        shsig[1] = 1.0 / (2 * pow(atom_type, 0.2));
        shpop[1] = 8;
        shsig[2] = 1.0 / (2 * pow(atom_type, 0.4));
        shpop[2] = 8;
        shsig[3] = 1.0 / (2 * pow(atom_type, 0.6));
        shpop[3] = 18;
        shsig[4] = 1.0 / (2 * pow(atom_type, 0.8));
        shpop[4] = 18;
        shsig[5] = 0.5;
        shpop[5] = atom_type - 54;
        break;
    default:
        err_not_impl_f("MBIS does not work with more than 6 shells. Please dont do that to me!", std::cout);
        break;
    }

    return std::make_pair(shsig, shpop);
}

std::pair<vec2, vec> get_shalpha_shpop(const int& atom_type) {
    // Order here is 11, 12, 13, 21, 22, 23, 31, 32, 33
    vec2 shalpha(constants::MBIS_function[atom_type], vec(6, 0.0));
    vec shpop(constants::MBIS_function[atom_type], 0.0);

    switch (constants::MBIS_function[atom_type]) {
    case 0:
        err_not_impl_f("MBIS does not work with ghost atoms. Please dont do that to me!", std::cout);
        break;
    case 1:
        shalpha[0][0] = pow(2.0 * atom_type, 2);
        shalpha[0][3] = pow(2.0 * atom_type, 2);
        shalpha[0][5] = pow(2.0 * atom_type, 2);
        shpop[0] = atom_type;
        break;
    case 2:
        shalpha[0][0] = pow(2.0 * atom_type, 2);
        shalpha[0][3] = pow(2.0 * atom_type, 2);
        shalpha[0][5] = pow(2.0 * atom_type, 2);
        shpop[0] = 2;
        shalpha[1][0] = 4.0;
        shalpha[1][3] = 4.0;
        shalpha[1][5] = 4.0;
        shpop[1] = atom_type - 2;
        break;
    case 3:
        shalpha[0][0] = pow(2.0 * atom_type, 2);
        shalpha[0][3] = pow(2.0 * atom_type, 2);
        shalpha[0][5] = pow(2.0 * atom_type, 2);
        shpop[0] = 2;
        shalpha[1][0] = pow(2.0 * sqrt(atom_type), 2);
        shalpha[1][3] = pow(2.0 * sqrt(atom_type), 2);
        shalpha[1][5] = pow(2.0 * sqrt(atom_type), 2);
        shpop[1] = 8;
        shalpha[2][0] = 4.0;
        shalpha[2][3] = 4.0;
        shalpha[2][5] = 4.0;
        shpop[2] = atom_type - 10;
        break;
    case 4:
        shalpha[0][0] = pow(2.0 * atom_type, 2);
        shalpha[0][3] = pow(2.0 * atom_type, 2);
        shalpha[0][5] = pow(2.0 * atom_type, 2);
        shpop[0] = 2;
        shalpha[1][0] = pow(2 * pow(atom_type, constants::c_23), 2);
        shalpha[1][3] = pow(2 * pow(atom_type, constants::c_23), 2);
        shalpha[1][5] = pow(2 * pow(atom_type, constants::c_23), 2);
        shpop[1] = 8;
        shalpha[2][0] = pow(2 * pow(atom_type, constants::c_13), 2);
        shalpha[2][3] = pow(2 * pow(atom_type, constants::c_13), 2);
        shalpha[2][5] = pow(2 * pow(atom_type, constants::c_13), 2);
        shpop[2] = 8;
        shalpha[3][0] = 4.0;
        shalpha[3][3] = 4.0;
        shalpha[3][5] = 4.0;
        shpop[3] = atom_type - 18;
        break;
    case 5:
        shalpha[0][0] = pow(2.0 * atom_type, 2);
        shalpha[0][3] = pow(2.0 * atom_type, 2);
        shalpha[0][5] = pow(2.0 * atom_type, 2);
        shpop[0] = 2;
        shalpha[1][0] = pow(2 * pow(atom_type, 0.25), 2);
        shalpha[1][3] = pow(2 * pow(atom_type, 0.25), 2);
        shalpha[1][5] = pow(2 * pow(atom_type, 0.25), 2);
        shpop[1] = 8;
        shalpha[2][0] = pow(2 * pow(atom_type, 0.5), 2);
        shalpha[2][3] = pow(2 * pow(atom_type, 0.5), 2);
        shalpha[2][5] = pow(2 * pow(atom_type, 0.5), 2);
        shpop[2] = 8;
        shalpha[3][0] = pow(2 * pow(atom_type, 0.75), 2);
        shalpha[3][3] = pow(2 * pow(atom_type, 0.75), 2);
        shalpha[3][5] = pow(2 * pow(atom_type, 0.75), 2);
        shpop[3] = 18;
        shalpha[4][0] = 4.0;
        shalpha[4][3] = 4.0;
        shalpha[4][5] = 4.0;
        shpop[4] = atom_type - 36;
        break;
    case 6:
        shalpha[0][0] = pow(2.0 * atom_type, 2);
        shalpha[0][3] = pow(2.0 * atom_type, 2);
        shalpha[0][5] = pow(2.0 * atom_type, 2);
        shpop[0] = 2;
        shalpha[1][0] = pow(2 * pow(atom_type, 0.2), 2);
        shalpha[1][3] = pow(2 * pow(atom_type, 0.2), 2);
        shalpha[1][5] = pow(2 * pow(atom_type, 0.2), 2);
        shpop[1] = 8;
        shalpha[2][0] = pow(2 * pow(atom_type, 0.4), 2);
        shalpha[2][3] = pow(2 * pow(atom_type, 0.4), 2);
        shalpha[2][5] = pow(2 * pow(atom_type, 0.4), 2);
        shpop[2] = 8;
        shalpha[3][0] = pow(2 * pow(atom_type, 0.6), 2);
        shalpha[3][3] = pow(2 * pow(atom_type, 0.6), 2);
        shalpha[3][5] = pow(2 * pow(atom_type, 0.6), 2);
        shpop[3] = 18;
        shalpha[4][0] = pow(2 * pow(atom_type, 0.8), 2);
        shalpha[4][3] = pow(2 * pow(atom_type, 0.8), 2);
        shalpha[4][5] = pow(2 * pow(atom_type, 0.8), 2);
        shpop[4] = 18;
        shalpha[5][0] = 4.0;
        shalpha[5][3] = 4.0;
        shalpha[5][5] = 4.0;
        shpop[5] = atom_type - 54;
        break;
    default:
        err_not_impl_f("MBIS does not work with more than 6 shells. Please dont do that to me!", std::cout);
        break;
    }

    return std::make_pair(shalpha, shpop);
}

// grids here are the total_grid in the SF calculator functions, just as a not to myself and future me
//Implemented according to the modified Multiwfn version by FJR and Anker
std::vector<std::pair<vec, vec>> make_MBIS_vectors(
    const WFN& wavy,
    const vec3& grid,
    const ivec& num_grid_points,
    const bool debug)
{
    using sp_vec = std::vector<std::pair<vec, vec>>;
    auto atoms = wavy.get_atoms();
    vec charges(atoms.size(), 0.0);
    vec last_charges(atoms.size(), 0.0);
    vec zeros_6(6, 0.0);
    sp_vec sig_pop_vector;
    for (auto& atom : atoms) {
        sig_pop_vector.emplace_back(get_shsig_shpop(atom.get_charge()));
    }
    const double crit = 0.001;
    sp_vec copy_of_input = sig_pop_vector;
    double varmax = 0, varsig = 0;
    // Cache atom coordinates and charges to avoid repeated function calls
    const int ncen = wavy.get_ncen();
    vec atom_coords(ncen * 3);
    ivec atom_charges(ncen);
    ivec nshell_cache(ncen);
    for (int j = 0; j < ncen; j++) {
        atom_coords[j * 3 + 0] = atoms[j].get_coordinate(0);
        atom_coords[j * 3 + 1] = atoms[j].get_coordinate(1);
        atom_coords[j * 3 + 2] = atoms[j].get_coordinate(2);
        atom_charges[j] = wavy.get_atom_charge(j);
        nshell_cache[j] = constants::MBIS_function[atom_charges[j]];
    }
    for (size_t it = 0; it < 2000; it++) {
        for (int j = 0; j < wavy.get_ncen(); j++) {
            std::fill(sig_pop_vector[j].first.begin(), sig_pop_vector[j].first.end(), 0.0);
            std::fill(sig_pop_vector[j].second.begin(), sig_pop_vector[j].second.end(), 0.0);
        }
        it == 0 ? std::cout << "Starting MBIS iterations..." << std::endl : std::cout << "MBIS iteration: " << it << " max change: " << varmax << std::endl;
        varmax = 0.0, varsig = 0.0;
        for (int i = 0; i < wavy.get_ncen(); i++) {
            const int end = num_grid_points[i];
            //Assuming 3 is the quadrature weight and 7 is the electron density 
            const double* b_weight = grid[i][5].data();
            const double* dens = grid[i][7].data();
            //This assumes GridIndex enum being X = 0, Y = 1, Z = 2
            const double* gx = grid[i][0].data();
            const double* gy = grid[i][1].data();
            const double* gz = grid[i][2].data();

#pragma omp parallel
            {
                vec rho0shell(ncen * 6, 0.0);
                double tmp = 0.0, density = 0.0, rho0 = 0.0, temp_res = 0.0, r0s = 0.0, sigval, dist_sq;
                int j, shell, nshell;
                double dx[3];
                sp_vec local = sig_pop_vector;
                for (j = 0; j < wavy.get_ncen(); j++) {
                    std::fill(local[j].first.begin(), local[j].first.end(), 0.0);
                    std::fill(local[j].second.begin(), local[j].second.end(), 0.0);
                }

#pragma omp for schedule(dynamic)
                for (int point = 0; point < end; point++) {
                    rho0 = 0.0;
                    vec dists(ncen, 0.0);
                    std::fill(rho0shell.begin(), rho0shell.end(), 0.0);
                    for (j = 0; j < wavy.get_ncen(); j++) {
                        //This assumes GridIndex enum being X = 0, Y = 1, Z = 2
                        dx[0] = gx[point] - atom_coords[j * 3 + 0];
                        dx[1] = gy[point] - atom_coords[j * 3 + 1];
                        dx[2] = gz[point] - atom_coords[j * 3 + 2];
                        // Check distance squared first to avoid sqrt for far away points
                        dist_sq = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
                        if (dist_sq > constants::far_away_sq)
                            continue;

                        dists[j] = std::sqrt(dist_sq);
                        nshell = nshell_cache[j];

                        for (shell = 0; shell < nshell; shell++) {
                            sigval = copy_of_input[j].first[shell];
                            tmp = copy_of_input[j].second[shell] * constants::INV_EIGHT_PI * pow(sigval, -3) * exp(-dists[j] / sigval);
                            if (abs(tmp) < 1e-20)
                                continue;
                            rho0shell[j * 6 + shell] = tmp;
                            rho0 += tmp;
                        }
                    }
                    density = b_weight[point] * dens[point];
                    for (j = 0; j < wavy.get_ncen(); j++) {
                        nshell = nshell_cache[j];
                        for (shell = 0; shell < nshell; shell++) {
                            r0s = rho0shell[j * 6 + shell];
                            if (r0s == 0)
                                continue;
                            temp_res = density * r0s / rho0;
                            local[j].first[shell] += temp_res * dists[j];
                            local[j].second[shell] += temp_res;
                        }
                    }
                }

                for (j = 0; j < wavy.get_ncen(); j++) {
                    nshell = nshell_cache[j];
                    for (shell = 0; shell < nshell; shell++) {
#pragma omp atomic
                        sig_pop_vector[j].first[shell] += local[j].first[shell];
#pragma omp atomic
                        sig_pop_vector[j].second[shell] += local[j].second[shell];
                    }
                }
            }
        }
        //back to the cycle main loop, we updated sig and pop based on information loss :)

        for (int i = 0; i < wavy.get_ncen(); i++) {
            for (int shell = 0; shell < constants::MBIS_function[wavy.get_atom_charge(i)]; shell++) {
                if (sig_pop_vector[i].second[shell] != 0)
                    sig_pop_vector[i].first[shell] /= 3.0 * sig_pop_vector[i].second[shell];
                varsig = std::max(varsig, std::abs(sig_pop_vector[i].first[shell] - copy_of_input[i].first[shell]));
            }
            charges[i] = wavy.get_atom_charge(i) - vec_sum(sig_pop_vector[i].second) + wavy.get_atom_ECP_electrons(i);
            varmax = std::max(varmax, std::abs(charges[i] - last_charges[i]));
            if (debug)
                std::cout << "Atom " << std::setw(3) << i << " charge: " << charges[i] << std::endl;
        }
        if (varmax < crit || varsig < crit) {
            std::cout << "MBIS converged after " << it << " iterations with max charge change: " << varmax << " and max sig change: " << varsig << std::endl;
            std::cout << "Promolecular charges:\n";
            for (int i = 0; i < wavy.get_ncen(); i++) {
                std::cout << "Atom " << std::setw(3) << i << ": " << charges[i] << "\n";
            }
            return sig_pop_vector;
        }
        copy_of_input = sig_pop_vector;
        last_charges = charges;

    }
    std::cout << "MBIS NOT converged after " << 200 << " iterations with max charge change: " << varmax << " and max sig change: " << varsig << std::endl << "Returning last iteration results." << std::endl << "BE CAREFUL WITH THESE RESULTS, THEY MIGHT NOT BE RELIABLE!" << std::endl;
    return sig_pop_vector;
}

// grids here are the total_grid in the SF calculator functions, just as a not to myself and future me
//Implemented according to the modified Multiwfn version by FJR and Anker
std::vector<std::pair<vec2, vec>> make_EMBIS_tensors(
    const WFN& wavy,
    const vec3& grid,
    const ivec& num_grid_points,
    const bool debug)
{
    using sp_vec = std::vector<std::pair<vec2, vec>>;
    const double crit = 0.001;
    auto atoms = wavy.get_atoms();
    vec charges(atoms.size(), 0.0);
    vec last_charges(atoms.size(), 0.0);
    double varmax = 0, varsig = 0;
    sp_vec sig_pop_vector;
    vec zeros_6(6, 0.0);
    for (auto& atom : atoms) {
        sig_pop_vector.emplace_back(get_shalpha_shpop(atom.get_charge()));
    }
    sp_vec copy_of_input = sig_pop_vector;

    const int ncen = wavy.get_ncen();
    vec atom_coords(ncen * 3);
    ivec atom_charges(ncen);
    ivec nshell_cache(ncen);
    for (int j = 0; j < ncen; j++) {
        atom_coords[j * 3 + 0] = atoms[j].get_coordinate(0);
        atom_coords[j * 3 + 1] = atoms[j].get_coordinate(1);
        atom_coords[j * 3 + 2] = atoms[j].get_coordinate(2);
        atom_charges[j] = wavy.get_atom_charge(j);
        nshell_cache[j] = constants::MBIS_function[atom_charges[j]];
    }

    for (size_t it = 0; it < 2000; it++) {
        for (int j = 0; j < ncen; j++) {
            std::fill(sig_pop_vector[j].first.begin(), sig_pop_vector[j].first.end(), zeros_6);
            std::fill(sig_pop_vector[j].second.begin(), sig_pop_vector[j].second.end(), 0.0);
        }
        it == 0 ? std::cout << "Starting EMBIS iterations..." << std::endl : std::cout << "EMBIS iteration: " << std::setw(4) << it << " max charge/alpha change: " << varsig << "/" << varmax << std::endl;
        varmax = 0.0, varsig = 0.0;
        for (int i = 0; i < ncen; i++) {
            const int end = num_grid_points[i];
            //Assuming 5 is the becke weight and 7 is the electron density 
            const double* b_weight = grid[i][5].data();
            const double* dens = grid[i][7].data();
            //This assumes GridIndex enum being X = 0, Y = 1, Z = 2
            const double* gx = grid[i][0].data();
            const double* gy = grid[i][1].data();
            const double* gz = grid[i][2].data();
#pragma omp parallel
            {
                vec rho0shell(ncen * 6, 0.0);
                double tmp = 0.0, density = 0.0, rho0 = 0.0, temp_res = 0.0, r0s = 0.0, g, det;
                int j, shell, nshell, ind;
                double* alpha;
                vec dx(ncen * 3);
                double* d_local;
                sp_vec local = sig_pop_vector;
                vec g_cache(6 * ncen, 0.0);
                double d_cache[6] = { 0.0 };
                std::pair<vec2, vec>* coi;

                // Precompute determinants for all atoms and shells - they don't change per point
                vec det_pop_cache(ncen * 6, 0.0);
                for (int k = 0; k < ncen; k++) {
                    nshell = nshell_cache[k];
                    coi = &copy_of_input[k];
                    for (shell = 0; shell < nshell; shell++) {
                        alpha = coi->first[shell].data();
                        det = sqrt(alpha[0] * alpha[3] * alpha[5] -
                            alpha[0] * alpha[4] * alpha[4] -
                            alpha[3] * alpha[2] * alpha[2] -
                            alpha[5] * alpha[1] * alpha[1] +
                            2 * alpha[1] * alpha[2] * alpha[4]);
                        det_pop_cache[k * 6 + shell] = coi->second[shell] * constants::INV_EIGHT_PI * det;
                    }
                }

                for (j = 0; j < ncen; j++) {
                    std::fill(local[j].first.begin(), local[j].first.end(), zeros_6);
                    std::fill(local[j].second.begin(), local[j].second.end(), 0.0);
                }
#pragma omp for schedule(dynamic)
                for (int point = 0; point < end; point++) {
                    rho0 = 0.0;
                    std::fill(rho0shell.begin(), rho0shell.end(), 0.0);
                    density = b_weight[point] * dens[point];
                    for (j = 0; j < ncen; j++) {
                        ind = j * 3;
                        d_local = dx.data() + ind;
                        d_local[0] = gx[point] - atom_coords[ind + 0];
                        d_local[1] = gy[point] - atom_coords[ind + 1];
                        d_local[2] = gz[point] - atom_coords[ind + 2];
                        d_cache[0] = d_local[0] * d_local[0];
                        d_cache[3] = d_local[1] * d_local[1];
                        d_cache[5] = d_local[2] * d_local[2];
                        // Check distance squared first to avoid sqrt for far away points
                        if (d_cache[0] + d_cache[3] + d_cache[5] > constants::far_away_sq)
                            continue;

                        d_cache[1] = d_local[0] * d_local[1];
                        d_cache[2] = d_local[0] * d_local[2];
                        d_cache[4] = d_local[1] * d_local[2];
                        nshell = nshell_cache[j];
                        coi = &copy_of_input[j];
                        ind *= 2;
                        for (shell = 0; shell < nshell; shell++) {
                            alpha = coi->first[shell].data();
                            // Order here is 11, 12, 13, 21, 22, 23, 31, 32, 33
                            g = sqrt(alpha[0] * d_cache[0] +
                                alpha[3] * d_cache[3] +
                                alpha[5] * d_cache[5] +
                                2 * alpha[1] * d_cache[1] +
                                2 * alpha[2] * d_cache[2] +
                                2 * alpha[4] * d_cache[4]);
                            g_cache[ind + shell] = g;
                            if (g > 42) // avoid vanishingly small contributions due to exp(-g) and potential overflow in exp(g)
                                continue;
                            tmp = det_pop_cache[ind + shell] * exp(-g);
                            if (abs(tmp) < constants::cutoff)
                                continue;
                            rho0shell[ind + shell] = tmp;
                            rho0 += tmp;
                        }
                    } //We first need to finish building total rho0 before we can calculate the contributions to the sigmas and populations
                    for (j = 0; j < ncen; j++) {
                        d_local = dx.data() + j * 3;
                        nshell = nshell_cache[j];
                        d_cache[0] = d_local[0] * d_local[0];
                        d_cache[1] = d_local[0] * d_local[1];
                        d_cache[2] = d_local[0] * d_local[2];
                        d_cache[3] = d_local[1] * d_local[1];
                        d_cache[4] = d_local[1] * d_local[2];
                        d_cache[5] = d_local[2] * d_local[2];
                        coi = &local[j];
                        ind = j * 6;
                        for (shell = 0; shell < nshell; shell++) {
                            r0s = rho0shell[ind + shell];
                            if (r0s <= constants::cutoff)
                                continue;
                            //alpha = copy_of_input[j].first[shell].data();
                            g = 1.0 / g_cache[ind + shell];
                            temp_res = density * r0s / rho0;
                            alpha = coi->first[shell].data();
                            // Order here is 11, 12, 13, 21, 22, 23, 31, 32, 33
                            alpha[0] += temp_res * d_cache[0] * g;
                            alpha[1] += temp_res * d_cache[1] * g;
                            alpha[2] += temp_res * d_cache[2] * g;
                            alpha[3] += temp_res * d_cache[3] * g;
                            alpha[4] += temp_res * d_cache[4] * g;
                            alpha[5] += temp_res * d_cache[5] * g;
                            coi->second[shell] += temp_res;
                        }
                    }
                }

                for (j = 0; j < ncen; j++) {
                    nshell = nshell_cache[j];
                    for (shell = 0; shell < nshell; shell++) {
                        for (int n = 0; n < 6; n++) {
#pragma omp atomic
                            sig_pop_vector[j].first[shell][n] += local[j].first[shell][n];
                        }
#pragma omp atomic
                        sig_pop_vector[j].second[shell] += local[j].second[shell];
                    }
                }
            }
        }
        //back to the cycle main loop, we updated sig and pop based on information loss :)

        for (int i = 0; i < ncen; i++) {
            for (int shell = 0; shell < constants::MBIS_function[wavy.get_atom_charge(i)]; shell++) {
                double det = 1 / (sig_pop_vector[i].first[shell][0] * sig_pop_vector[i].first[shell][3] * sig_pop_vector[i].first[shell][5] -
                    sig_pop_vector[i].first[shell][0] * sig_pop_vector[i].first[shell][4] * sig_pop_vector[i].first[shell][4] -
                    sig_pop_vector[i].first[shell][3] * sig_pop_vector[i].first[shell][2] * sig_pop_vector[i].first[shell][2] -
                    sig_pop_vector[i].first[shell][5] * sig_pop_vector[i].first[shell][1] * sig_pop_vector[i].first[shell][1] +
                    2 * sig_pop_vector[i].first[shell][1] * sig_pop_vector[i].first[shell][2] * sig_pop_vector[i].first[shell][4]);
                if (det < 1E-14)
                    continue;
                auto sigmas = sig_pop_vector[i].first[shell];
                //                0   1   2   X   3   4   X   X   5
                // Order here is 11, 12, 13, 21, 22, 23, 31, 32, 33
                sig_pop_vector[i].first[shell][0] = (sigmas[3] * sigmas[5] - sigmas[4] * sigmas[4]) * det;
                sig_pop_vector[i].first[shell][1] = (sigmas[2] * sigmas[4] - sigmas[1] * sigmas[5]) * det;
                sig_pop_vector[i].first[shell][2] = (sigmas[1] * sigmas[4] - sigmas[2] * sigmas[3]) * det;
                sig_pop_vector[i].first[shell][3] = (sigmas[0] * sigmas[5] - sigmas[2] * sigmas[2]) * det;
                sig_pop_vector[i].first[shell][4] = (sigmas[2] * sigmas[1] - sigmas[0] * sigmas[4]) * det;
                sig_pop_vector[i].first[shell][5] = (sigmas[0] * sigmas[3] - sigmas[1] * sigmas[1]) * det;

                if (sig_pop_vector[i].second[shell] != 0)
                    for (int n = 0; n < 6; n++)
                        sig_pop_vector[i].first[shell][n] *= sig_pop_vector[i].second[shell];
                for (int n = 0; n < 6; n++)
                    varsig = std::max(varsig, std::abs(sig_pop_vector[i].first[shell][n] - copy_of_input[i].first[shell][n]));
            }
            charges[i] = wavy.get_atom_charge(i) - vec_sum(sig_pop_vector[i].second) + wavy.get_atom_ECP_electrons(i);
            varmax = std::max(varmax, std::abs(charges[i] - last_charges[i]));
            if (debug)
                std::cout << "Atom " << std::setw(3) << i << " charge: " << charges[i] << std::endl;
        }
        if (varmax < crit || varsig < crit) {
            std::cout << "EMBIS converged after " << it << " iterations with max charge change: " << varmax << " and max sig change: " << varsig << std::endl;
            std::cout << "Promolecular charges:\n";
            for (int i = 0; i < wavy.get_ncen(); i++) {
                std::cout << "Atom " << std::setw(3) << i << ": " << charges[i] << "\n";
            }
            return sig_pop_vector;
        }
        copy_of_input = sig_pop_vector;
        last_charges = charges;

    }
    std::cout << "EMBIS NOT converged after " << 200 << " iterations with max charge change: " << varmax << " and max sig change: " << varsig << std::endl << "Returning last iteration results." << std::endl << "BE CAREFUL WITH THESE RESULTS, THEY MIGHT NOT BE RELIABLE!" << std::endl;
    return sig_pop_vector;
}

#pragma GCC optimize("-fno-fast-math") //gcc fails to converge in this function when using fast math, so not use it
// TCA 106, 178 (2001), eq. 19
double get_r_outer(const double& max_error,
    const double& alpha_outer,
    const int& l,
    const double& guess)
{
    const int m = 2 * l;
    double r = guess;
    double r_old = 1.0e50;
    double c, a, e;
    double step = 0.5;
    double sign, sign_old;
    double f = 1.0e50;

    (f > max_error) ? (sign = 1.0) : (sign = -1.0);

    while (std::abs(r_old - r) > constants::cutoff)
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
double get_h(const double& max_error, const int& l, const double& guess)
{
    const int m = 2 * l;
    double h = guess;
    double h_old = h * 2;
    double step = 0.1 * guess;
    double sign = -1.0, sign_old, f, pl, rd0, e0;
    const double cm = constants::TG32 / tgamma((m + 3.0) / 2.0);

    while (std::abs(h_old - h) > constants::cutoff)
    {
        e0 = std::exp(-constants::PI2 / (2.0 * h));
        pl = std::pow(constants::PI / h, l);
        rd0 = constants::C0 / h * e0;
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
#pragma GCC optimize("fast-math")
