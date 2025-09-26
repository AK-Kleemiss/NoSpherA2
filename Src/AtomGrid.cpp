#include "pch.h"
#include "convenience.h"
#include "AtomGrid.h"
#include "sphere_lebedev_rule.h"
#include "constants.h"
#include "wfn_class.h"

#ifdef _WIN32
#include <algorithm>
#include <io.h>
#endif

// Helper to describe extrema along a line between two atoms
struct DensityExtremum
{
    enum class Type { Minimum, Maximum } type;
    double x{0}, y{0}, z{0};     // Cartesian position (Bohr)
    double density{0};           // Electron density at extremum
    double distA{0};             // Distance to atom A (Bohr)
    double distB{0};             // Distance to atom B (Bohr)
    double t{0};                 // Parameter along AB in [0,1]
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
    const double ax = wfn.get_atom_coordinate((unsigned)atomA, 0);
    const double ay = wfn.get_atom_coordinate((unsigned)atomA, 1);
    const double az = wfn.get_atom_coordinate((unsigned)atomA, 2);
    const double bx = wfn.get_atom_coordinate((unsigned)atomB, 0);
    const double by = wfn.get_atom_coordinate((unsigned)atomB, 1);
    const double bz = wfn.get_atom_coordinate((unsigned)atomB, 2);

    const double dx = bx - ax;
    const double dy = by - ay;
    const double dz = bz - az;

    const double ab_len = std::sqrt(dx*dx + dy*dy + dz*dz);
    if (ab_len <= 0) return result;

    // Sampling t in (0,1); avoid endpoints to reduce numerical issues near nuclei
    const double eps = 1e-6;
    const double t0 = eps, t1 = 1.0 - eps;
    const double dt = (t1 - t0) / (samples - 1);

    vec dens(samples);
    vec ts(samples);
#pragma omp parallel for
    for (int i = 0; i < samples; ++i) {
        const double t = t0 + i * dt;
        ts[i] = t;
        const double x = ax + dx * t;
        const double y = ay + dy * t;
        const double z = az + dz * t;
        dens[i] = wfn.compute_dens(x, y, z);
    }

    // Detect extrema via sign changes of the discrete derivative
    for (int i = 1; i < samples - 1; ++i) {
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
                const double x = ax + dx * t_ext;
                const double y = ay + dy * t_ext;
                const double z = az + dz * t_ext;
                f_ext = wfn.compute_dens(x, y, z);
            }
        }

        DensityExtremum e;
        e.type = isMax ? DensityExtremum::Type::Maximum : DensityExtremum::Type::Minimum;
        e.t = t_ext;
        e.x = ax + dx * t_ext;
        e.y = ay + dy * t_ext;
        e.z = az + dz * t_ext;
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
    const int start = (int) atom_grid_x_bohr_.size();
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

int AtomGrid::get_num_grid_points() const { return (int) atom_grid_x_bohr_.size(); }

int AtomGrid::get_num_radial_grid_points() const { return num_radial_grid_points_; }

vec2 make_chi(const WFN& wfn, int samples = 50, bool refine = true, bool debug = false) {
	vec2 chi(wfn.get_ncen(), vec(wfn.get_ncen(), 0.0));
    std::vector<std::vector<bool>> neighbours(wfn.get_ncen(), bvec(wfn.get_ncen(), true));
    double rijx2, rijy2, rijz2, xdist, disth;
    for (int a = 0; a < wfn.get_ncen(); a++){
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
                    chi[a][b] = 1.0;
                    chi[b][a] = 1.0;
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
                                chi[a][b] = 1.0;
                                chi[b][a] = 1.0;
                                continue;
                            }
                        }
                        else {
                            std::cout << "WARNING: Unexpected number of extrema (" << extrema.size() << ") found between atoms " << a << " and " << b << ". Setting chi to 1.0." << std::endl;
                            chi[a][b] = 1.0;
                            chi[b][a] = 1.0;
                            continue;
                        }
                    }
                    // And save the "chi" in the same matrix
                    chi[a][b] = extrema[use_extr].distA / extrema[use_extr].distB;
                    chi[b][a] = 1.0 / chi[a][b];
                    // For if one have a very polarized bond (bcp super displaced towards one of the atoms
                    double ri = extrema[use_extr].distA / (extrema[use_extr].distA + extrema[use_extr].distB);
                    if (ri >= 0.90 || ri <= 0.10) {
                        std::cout << "WARNING, LARGE POLARIZATION. PLEASE CHECK" << std::endl;
                    }
                    // Some printing on the ratio between pair of atoms, for if strange results (or against chemical intuition) are obtained one can see if it is TFVC problem
                    // This is helpful to see if the "bcp" evaluation is broken, or one face a nasty scenario that is not yet included in the code...
                    if(debug)
                        std::cout << "TFVC -- Atom pair: " << a << ", " << b << ", Ratio: " << chi[a][b] << std::endl;
                }
            }
            else {
                chi[a][b] = 1.0;
                chi[b][a] = 1.0;
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
    vec2& chi,
    bool debug) const
{
    //vec pa(num_centers);
    //double x = -1.836008, y = -0.292080, z = -0.861742;
    //double weight1 = get_tfvc_w(num_centers,
    //    proton_charges,
    //    x_coordinates_bohr,
    //    y_coordinates_bohr,
    //    z_coordinates_bohr,
    //    0,
    //    x,
    //    y,
	//	z, pa, chi);
    //double weight2 = get_tfvc_w(num_centers,
    //    proton_charges,
    //    x_coordinates_bohr,
    //    y_coordinates_bohr,
    //    z_coordinates_bohr,
    //    1,
    //    x,
    //    y,
	//	z, pa, chi);

    if (num_centers > 1) {
        if (chi.size() == 0)
            chi = make_chi(wfn, 20, true, debug);
#pragma omp parallel
        {
            vec pa_b(num_centers);
            vec pa_tv(num_centers);
            std::pair<double, double> result_weights;
            double temp;
#pragma omp for
            for (int ipoint = 0; ipoint < get_num_grid_points(); ipoint++) {
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
			    grid_becke_w[ipoint] = temp * result_weights.first;
				grid_TFVC_w[ipoint] = temp * result_weights.second;
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
std::pair<double, double> get_integration_weights(const int& num_centers,
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
    const vec2& chi)
{
    double mu_ab, nu_ab, f, dist_ab;
    double dist_a, dist_b;
    double vx, vy, vz;
	double R_a, R_b, chi_becke, u_ab;

    for (int a = 0; a < num_centers; a++) {
        pa_b[a] = 1.0;
		pa_tv[a] = 1.0;
    }

    for (int a = 0; a < num_centers; a++) {
        vx = x_coordinates_bohr[a] - x;
        vy = y_coordinates_bohr[a] - y;
        vz = z_coordinates_bohr[a] - z;
        dist_a = std::sqrt(vx * vx + vy * vy + vz * vz);

        if (dist_a > constants::far_away) {
            pa_b[a] = 0.0;
            pa_tv[a] = 0.0;
            continue;
        }

        R_a = constants::bragg_angstrom[proton_charges[a]];

        for (int b = a + 1; b < num_centers; b++) {
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

            nu_ab = (1.0 + mu_ab - chi[a][b] * (1.0 - mu_ab)) / (1.0 + mu_ab + chi[a][b] * (1.0 - mu_ab));

            f = f4(nu_ab);

            if (std::abs(1.0 - f) < constants::cutoff)
                pa_tv[a] = 0.0;
            else {
				//Reduce numerical jittering
                if (pa_tv[a] > 1E-250 || pa_tv[a] < -1E-250)
                    pa_tv[a] *= 0.5 * (1.0 - f);
                else
                    pa_tv[a] = 0.0;
                if (pa_tv[b] > 1E-250 || pa_tv[b] < -1E-250)
                    pa_tv[b] *= 0.5 * (1.0 + f);
                else
                    pa_tv[b] = 0.0;
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

            f = f3(nu_ab);

            if (std::abs(1.0 - f) < constants::cutoff)
                pa_b[a] = 0.0;
            else {
                if (pa_b[a] > 1E-250 || pa_b[a] < -1E-250)
                    pa_b[a] *= 0.5 * (1.0 - f);
                else
                    pa_b[a] = 0.0;
                if (pa_b[b] > 1E-250 || pa_b[b] < -1E-250)
                    pa_b[b] *= 0.5 * (1.0 + f);
                else
                    pa_b[b] = 0.0;
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

  double r = (d - log(1.0 / max_error)) * 2./3.;
  r = exp(r) / (alpha_inner);
  r = std::sqrt(r);

  return r;
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
