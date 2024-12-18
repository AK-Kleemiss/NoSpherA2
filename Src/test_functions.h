#pragma once
#include "spherical_density.h"
#include "scattering_factors.h"
#include "AtomGrid.h"
#include "convenience.h"
#include "npy.h"
#include "properties.h"
#include "JKFit.h"
#include "SALTED_utilities.h"
#include "sphere_lebedev_rule.h"
#if has_RAS
#include "math.h"
#include "rascaline.hpp"
#include "metatensor.h"
#endif
#ifdef __APPLE__
#include "TargetConditionals.h"
#endif

void thakkar_d_test(options &opt)
{
    using namespace std;
    Thakkar Os(76), Ca(20), C(6), O(8), H(1), P(15);
    Thakkar_Cation C_cat(6), O_cat(8), P_cat(15), Ca_cat(20);
    Thakkar_Anion C_an(6), O_an(8), H_an(1), P_an(15);
    double k_value = 0.0;
    if (!opt.electron_diffraction)
    {
        ofstream result("thakkar.dat", ios::out);
        for (double i = 0.001; i <= 4.0; i += 0.001)
        {
            k_value = constants::bohr2ang(constants::FOUR_PI * i);
            result << showpos << setw(6) << setprecision(3) << fixed << i;
            result << showpos << setw(16) << setprecision(8) << scientific << H.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << C.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << O.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << P.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << Ca.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << Os.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << C_cat.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << O_cat.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << P_cat.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << Ca_cat.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << H_an.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << C_an.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << O_an.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << P_an.get_form_factor(k_value);
            result << endl;
        }
        result.flush();
        result.close();
    }
    else
    {
        ofstream result("thakkar_ED.dat", ios::out);
        for (double i = 0.001; i <= 4.0; i += 0.001)
        {
            k_value = constants::bohr2ang(constants::FOUR_PI * i);
            result << showpos << setw(6) << setprecision(3) << fixed << i;
            complex<double> temp = H.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(1, temp, i, 0).real();
            temp = C.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(6, temp, i, 0).real();
            temp = O.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(8, temp, i, 0).real();
            temp = P.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(15, temp, i, 0).real();
            temp = Ca.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(20, temp, i, 0).real();
            temp = Os.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(76, temp, i, 0).real();
            temp = 0.0;
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(1, temp, i, 0).real();
            temp = C_cat.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(6, temp, i, 0).real();
            temp = O_cat.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(8, temp, i, 0).real();
            temp = P_cat.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(15, temp, i, 0).real();
            temp = Ca_cat.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(20, temp, i, 0).real();
            temp = H_an.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(1, temp, i, 0).real();
            temp = C_an.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(6, temp, i, 0).real();
            temp = O_an.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(8, temp, i, 0).real();
            temp = P_an.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(15, temp, i, 0).real();
            temp = C.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(6, temp, i, 1).real();
            temp = O.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(8, temp, i, 1).real();
            temp = C.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(6, temp, i, -1).real();
            temp = O.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(7, temp, i, -1).real();
            result << endl;
        }
        result.flush();
        result.close();
    }
}

void sfac_scan(options &opt, std::ostream &log_file)
{
    using namespace std;
    std::vector<WFN> wavy;
    auto t = new WFN(1);
    wavy.push_back(*t);
    delete t;
    wavy[0].read_known_wavefunction_format(opt.wfn, std::cout, opt.debug);
    Thakkar O(wavy[0].atoms[0].charge);
    Thakkar_Cation O_cat(wavy[0].atoms[0].charge);
    Thakkar_Anion O_an(wavy[0].atoms[0].charge);
    err_checkf(wavy[0].get_ncen() != 0, "No Atoms in the wavefunction, this will not work!! ABORTING!!", std::cout);
    err_checkf(exists(opt.cif), "CIF does not exists!", std::cout);
    // err_checkf(exists(asym_cif), "Asym/Wfn CIF does not exists!", file);

    std::vector<time_point> time_points;
    std::vector<std::string> time_descriptions;

    cell unit_cell(opt.cif, std::cout, opt.debug);
    ifstream cif_input(opt.cif.c_str(), std::ios::in);
    ivec atom_type_list;
    ivec asym_atom_to_type_list;
    ivec asym_atom_list;
    bvec constant_atoms;
    bvec needs_grid(wavy[0].get_ncen(), false);
    svec known_atoms;

    auto labels = read_atoms_from_CIF(cif_input,
                                      opt.groups[0],
                                      unit_cell,
                                      wavy[0],
                                      known_atoms,
                                      atom_type_list,
                                      asym_atom_to_type_list,
                                      asym_atom_list,
                                      needs_grid,
                                      std::cout,
                                      constant_atoms,
                                      opt.debug);

    cif_input.close();
    vec2 d1, d2, d3, dens;

    make_hirshfeld_grids(opt.pbc,
                         opt.accuracy,
                         unit_cell,
                         wavy[0],
                         atom_type_list,
                         asym_atom_list,
                         needs_grid,
                         d1, d2, d3, dens,
                         labels,
                         std::cout,
                         time_points,
                         time_descriptions,
                         opt.debug,
                         opt.no_date);

    std::cout << "finished partitioning" << endl;
    const int size = 4000;
    const int phi_size = 50;
    const int theta_size = 50;
    const double phi_step = 360.0 / phi_size * constants::PI_180;
    const double theta_step = 180.0 / phi_size * constants::PI_180;

    // This bit is basically the substitute for make_k_pts, where we sample the whole sphere
    //  by iterating over both spherical angles by a fixed step defined above
    vec2 k_pt;
    k_pt.resize(4);
#pragma omp parallel for
    for (int i = 0; i < 4; i++)
        k_pt[i].resize(size * phi_size * theta_size, 0.0);

        // int null = 0;
#pragma omp parallel for
    for (int ref = 1; ref <= size; ref++)
    {
        for (int p = 0; p < phi_size; p++)
        {
            for (int _t = 0; _t < theta_size; _t++)
            {
                int ind = _t + (p + (ref - 1) * phi_size) * theta_size;
                double k_length = constants::bohr2ang(constants::FOUR_PI * ref / size * opt.d_sfac_scan);
                k_pt[0][ind] = k_length * sin(_t * theta_step) * cos(p * phi_step);
                k_pt[1][ind] = k_length * sin(_t * theta_step) * sin(p * phi_step);
                k_pt[2][ind] = k_length * cos(_t * theta_step);
                k_pt[3][ind] = k_length;
            }
        }
    }
    // below is a strip of Calc_SF without the file IO or progress bar
    cvec2 sf;

    const int imax = (int)dens.size();
    const int smax = (int)k_pt[0].size();
    int pmax = (int)dens[0].size();
    std::cout << "Done with making k_pt " << smax << " " << imax << " " << pmax << endl;
    sf.resize(imax);
#pragma omp parallel for
    for (int i = 0; i < imax; i++)
        sf[i].resize(k_pt[0].size());
    double *dens_local, *d1_local, *d2_local, *d3_local;
    complex<double> *sf_local;
    const double *k1_local = k_pt[0].data();
    const double *k2_local = k_pt[1].data();
    const double *k3_local = k_pt[2].data();
    double work, rho;
    ProgressBar * progress = new ProgressBar(smax, 50, "=", " ", "Calculating Scattering factors");
    for (int i = 0; i < imax; i++)
    {
        pmax = (int)dens[i].size();
        dens_local = dens[i].data();
        d1_local = d1[i].data();
        d2_local = d2[i].data();
        d3_local = d3[i].data();
        sf_local = sf[i].data();
#pragma omp parallel for private(work, rho)
        for (int s = 0; s < smax; s++)
        {
            for (int p = pmax - 1; p >= 0; p--)
            {
                rho = dens_local[p];
                work = k1_local[s] * d1_local[p] + k2_local[s] * d2_local[p] + k3_local[s] * d3_local[p];
                sf_local[s] += complex<double>(rho * cos(work), rho * sin(work));
            }
            progress->update();
        }
    }
    delete (progress);
    if (true)
    { // Change if oyu do not want X-ray
        ofstream result("sfacs.dat", ios::out);
        log_file << "Writing X-ray sfacs...";
        log_file.flush();
        // Now we just need to write the result to a file, together with the spherical results and separated for valence and core
        for (int i = 0; i < k_pt[0].size(); i++)
        {
            result << showpos << setw(8) << setprecision(5) << fixed << constants::ang2bohr(k_pt[3][i] / constants::FOUR_PI);
            result << showpos << setw(16) << setprecision(8) << scientific << O.get_form_factor((k_pt[3][i]));
            result << showpos << setw(16) << setprecision(8) << scientific << O_an.get_form_factor((k_pt[3][i]));
            result << showpos << setw(16) << setprecision(8) << scientific << O_cat.get_form_factor((k_pt[3][i]));
            result << showpos << setw(16) << setprecision(8) << scientific << sqrt(pow(sf[0][i].real(), 2) + pow(sf[0][i].imag(), 2));
            result << showpos << setw(16) << setprecision(8) << scientific << O.get_custom_form_factor((k_pt[3][i]), 1, 0, 0, 0, 0, 0, 0, 0);
            result << showpos << setw(16) << setprecision(8) << scientific << O.get_custom_form_factor((k_pt[3][i]), 2, 1, 0, 0, 1, 0, 0, 0);
            result << "\n";
        }
        log_file << " ... done!" << endl;
        result.flush();
        result.close();
    }
    if (true)
    { // change if you do not want ED sfacs
        log_file << "Writing ED sfacs...";
        log_file.flush();
        ofstream result = ofstream("sfacs_ED.dat", ios::out);
        const double fact = 0.023934;
        double h2;
        for (int s = 0; s < k_pt[0].size(); s++)
        {
            h2 = pow(constants::ang2bohr(k_pt[3][s] / constants::FOUR_PI), 2);
            sf[0][s] = std::complex<double>(fact * (wavy[0].get_atom_charge(0) - sf[0][s].real()) / h2, -fact * sf[0][s].imag() / h2);

            result << showpos << setw(8) << setprecision(5) << fixed << constants::ang2bohr(k_pt[3][s] / constants::FOUR_PI);
            double temp = fact * (wavy[0].get_atom_charge(0) - O.get_form_factor(k_pt[3][s])) / h2;
            result << showpos << setw(16) << setprecision(8) << scientific << temp;
            temp = fact * (wavy[0].get_atom_charge(0) - O_an.get_form_factor(k_pt[3][s])) / h2;
            result << showpos << setw(16) << setprecision(8) << scientific << temp;
            temp = fact * (wavy[0].get_atom_charge(0) - O_cat.get_form_factor(k_pt[3][s])) / h2;
            result << showpos << setw(16) << setprecision(8) << scientific << temp;

            result << showpos << setw(16) << setprecision(8) << scientific << sqrt(pow(sf[0][s].real(), 2) + pow(sf[0][s].imag(), 2));

            temp = fact * (2 - O.get_custom_form_factor(k_pt[3][s], 1, 0, 0, 0, 0, 0, 0, 0)) / h2;
            result << showpos << setw(16) << setprecision(8) << scientific << temp;
            temp = fact * (6 - O.get_custom_form_factor(k_pt[3][s], 2, 1, 0, 0, 1, 0, 0, 0)) / h2;
            result << showpos << setw(16) << setprecision(8) << scientific << temp;
            result << "\n";
        }
        result.flush();
        result.close();
        log_file << " ... done!" << endl;
    }
}

template <typename T>
std::vector<T> geomspace(T start, T stop, int num)
{
    std::vector<T> result;
    T delta = (log(stop) - log(start)) / (num - 1);
    for (int i = 0; i < num; i++)
    {
        T exponent = log(start) + i * delta;
        result.push_back(exp(exponent));
    }
    return result;
}

double calc_spherically_averaged_at_r(const WFN &wavy,
                                      const double &r,
                                      const double rel_precision = 1E-4,
                                      const double start_angle = 5.0,
                                      const bool print = false)
{
    double old_result = 1E100;
    double new_result = 0.9E100;
    double angle = start_angle;
    double ratio = abs(old_result / new_result - 1.0);
    vec2 d;
    vec _phi(wavy.get_nmo(), 0.0);
    d.resize(16);
    for (int i = 0; i < 16; i++)
        d[i].resize(1, 0.0);

    while (ratio > rel_precision)
    {
        const double da = constants::PI / angle;
        const double da2 = da * da;
        // Make angular grids
        double x0, y0, z0, st0;
        double v = 0;
        for (double phi = 0; phi <= constants::TWO_PI; phi += da)
        {
            const double cp = cos(phi);
            const double sp = sin(phi);
            for (double theta = da; theta <= constants::PI; theta += da)
            {
                const double st = sin(theta);
                st0 = st * da2;
                x0 = st * cp, y0 = st * sp, z0 = cos(theta);
                v += wavy.compute_dens(r * x0, r * y0, r * z0, d, _phi) * st0;
            }
        }
        if (v != 0.0)
        {
            old_result = new_result;
            new_result = v / constants::FOUR_PI;
            ratio = abs(old_result / new_result - 1.0);
            if ((new_result < 1E-15 && angle > 36)) {
#pragma omp critical
                std::cout << "Aborted due to too small density at r = " << r << " and angle = " << angle << std::endl;
                break;
            }
            if (new_result > 1E90) {
#pragma omp critical
                std::cout << "Aborted due too large density value at r = " << r << std::endl;
                break;
            }
            if (angle > 360) {
#pragma omp critical
                std::cout << "Aborted due to too large angle = " << angle << std::endl;
                break;
            }
        }
        else {
            new_result = 0.0;
#pragma omp critical
            std::cout << "ZERO result at r = " << r << std::endl;
            break;
        }
        angle *= 1.2;
    }
    if (print)
        std::cout << "Done with " << std::setw(10) << std::setprecision(5) << r << " " << (ratio > rel_precision) << std::setw(14) << angle << std::endl;
    return new_result;
}

double calc_grid_averaged_at_r(const WFN& wavy,
    const double& r,
    const int& min_angular = 60,
    const int& max_angular = 1000,
    const bool print = false)
{

    const int min_num_angular_points_closest =
        constants::get_closest_num_angular(min_angular);
    const int max_num_angular_points_closest =
        constants::get_closest_num_angular(max_angular);
    err_checkf(min_num_angular_points_closest != -1 && max_num_angular_points_closest != -1, "No valid value for angular number found!", std::cout);

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

    const double rb = constants::bragg_angstrom[wavy.atoms[0].charge] / (5.0E10 * constants::a0);

    int num_angular = max_num_angular_points_closest;
    if (r < rb) {
        num_angular = static_cast<int>(max_num_angular_points_closest *
            (r / rb));
        num_angular = constants::get_closest_num_angular(num_angular);
        err_checkf(num_angular != -1, "No valid value for angular number found!", std::cout);
        if (num_angular < min_num_angular_points_closest)
            num_angular = min_num_angular_points_closest;
    }

    angular_off = constants::get_angular_order(num_angular) * constants::MAG;
    err_checkf(angular_off != -constants::MAG, "Invalid angular order!", std::cout);
    const int start = 0;
    vec2 d;
    vec _phi(wavy.get_nmo(), 0.0);
    d.resize(16);
    for (int i = 0; i < 16; i++)
        d[i].resize(1, 0.0);
    angular_off -= start;
    const int size = start + num_angular;
    int p = 0;
    double dens = 0.0;
    for (int iang = start; iang < size; iang++) {
        p = angular_off + iang;
        double x = angular_x[p] * r + wavy.atoms[0].x;
        double y = angular_y[p] * r + wavy.atoms[0].y;
        double z = angular_z[p] * r + wavy.atoms[0].z;
        dens += wavy.compute_dens(x, y, z, d, _phi) * constants::FOUR_PI * angular_w[p];
    }
    if (print)
        std::cout << "Done with " << std::setw(10) << std::setprecision(5) << r << std::endl;
    return dens;
}
/// not like above //////////////////////////////////////////////////////////////////////////
double calc_hirsh_grid_averaged_at_r(const WFN& wavy,
    const int i, /// now counts for several atoms
    const double& r,
    const int& min_angular = 60,
    const int& max_angular = 1000,
    const bool print = false)
{

    const int min_num_angular_points_closest =
        constants::get_closest_num_angular(min_angular);
    const int max_num_angular_points_closest =
        constants::get_closest_num_angular(max_angular);
    err_checkf(min_num_angular_points_closest != -1 && max_num_angular_points_closest != -1, "No valid value for angular number found!", std::cout);

    vec angular_x(constants::max_LT * constants::MAG, 0.0);
    vec angular_y(constants::max_LT * constants::MAG, 0.0);
    vec angular_z(constants::max_LT * constants::MAG, 0.0);
    vec angular_w(constants::max_LT * constants::MAG, 0.0);
    int angular_off;
    lebedev_sphere ls;

    for (int j = constants::get_angular_order(min_num_angular_points_closest); j < constants::get_angular_order(max_num_angular_points_closest) + 1; j++) {
        angular_off = j * constants::MAG;
        ls.ld_by_order(constants::lebedev_table[j],
            angular_x.data() + angular_off,
            angular_y.data() + angular_off,
            angular_z.data() + angular_off,
            angular_w.data() + angular_off);
    }

    const double rb = constants::bragg_angstrom[wavy.atoms[i].charge] / (5.0E10 * constants::a0);

    int num_angular = max_num_angular_points_closest;
    if (r < rb) {
        num_angular = static_cast<int>(max_num_angular_points_closest *
            (r / rb));
        num_angular = constants::get_closest_num_angular(num_angular);
        err_checkf(num_angular != -1, "No valid value for angular number found!", std::cout);
        if (num_angular < min_num_angular_points_closest)
            num_angular = min_num_angular_points_closest;
    }

    angular_off = constants::get_angular_order(num_angular) * constants::MAG;
    err_checkf(angular_off != -constants::MAG, "Invalid angular order!", std::cout);
    const int start = 0;
    angular_off -= start;
    const int size = start + num_angular;
    int p = 0;
    double dens = 0.0;
	Thakkar A(wavy.atoms[i].charge); /// for the atom i, a Thakkar object is created
#pragma omp parallel
    {
        vec2 d(16);
        vec _phi(wavy.get_nmo(), 0.0);
        for (int j = 0; j < 16; j++)
            d[j].resize(wavy.get_ncen(), 0.0);
#pragma omp for reduction(+:dens)
        for (int iang = start; iang < size; iang++) {
            p = angular_off + iang;
            const double x = angular_x[p] * r + wavy.atoms[i].x; /// moving the lebedev points to the atom i
            const double y = angular_y[p] * r + wavy.atoms[i].y;
            const double z = angular_z[p] * r + wavy.atoms[i].z;
            const std::array<double, 3> d_ = { angular_x[p] * r, angular_y[p] * r, angular_z[p] * r }; /// as the function get_radial_density needs a distance, the distance is calculated. The atom A is in the origin, the wavy.atoms[i].x is 0 
            const double dist = array_length(d_);
            const double rho_a = A.get_radial_density(dist); /// the radial density of the atom i is calculated
            double rho_all = rho_a; /// the molecular density based on pro-atoms
            for (int atom = 0; atom < wavy.atoms.size(); atom++) { /// a new for loop is started, which calculates the sum of rhos, with respect to the lebedev points
                if (atom == i) continue; /// if the atom is the same as the atom i, the loop is skipped
                const std::array<double, 3> d_atom = { x - wavy.atoms[atom].x, y - wavy.atoms[atom].y, z - wavy.atoms[atom].z };
                const double dist_atom = array_length(d_atom); /// is like d_, but for another atom 
                Thakkar thakkar_atom(wavy.atoms[atom].charge);
                rho_all += thakkar_atom.get_radial_density(dist_atom);
            }
            const double hirsh_weight = rho_a / rho_all;
            dens += wavy.compute_dens(x, y, z, d, _phi) * hirsh_weight * constants::FOUR_PI * angular_w[p];
        }
    }
    if (print)
        std::cout << "Done with " << std::setw(10) << std::setprecision(5) << r << std::endl;
    return dens;
}

void bondwise_laplacian_plots(std::string& wfn_name) {
    char cwd[1024];
#ifdef _WIN32
    if (_getcwd(cwd, sizeof(cwd)) != NULL)
#else
    if (getcwd(cwd, sizeof(cwd)) != NULL)
#endif
        std::cout << "Current working dir: " << cwd << std::endl;
    WFN wavy(9);
    wavy.read_known_wavefunction_format(wfn_name, std::cout);

    err_checkf(wavy.get_ncen() != 0, "No Atoms in the wavefunction, this will not work!! ABORTING!!", std::cout);

    int points = 1001;

    for (int i = 0; i < wavy.get_ncen(); i++) {
        for (int j = i+1; j < wavy.get_ncen(); j++) {
            std::string path = cwd;
            double distance = sqrt(pow(wavy.atoms[i].x - wavy.atoms[j].x, 2) + pow(wavy.atoms[i].y - wavy.atoms[j].y, 2) + pow(wavy.atoms[i].z - wavy.atoms[j].z, 2));
            double svdW = constants::ang2bohr(constants::covalent_radii[wavy.atoms[i].charge] + constants::covalent_radii[wavy.atoms[j].charge]);
            if (distance < 1.35 * svdW)
            {
                std::cout << "Bond between " << i << " ("<< wavy.atoms[i].charge << ") and " << j << " (" << wavy.atoms[j].charge << ") with distance " << distance << " and svdW " << svdW << std::endl;
                const vec bond_vec = { (wavy.atoms[j].x - wavy.atoms[i].x)/points, (wavy.atoms[j].y - wavy.atoms[i].y)/points, (wavy.atoms[j].z - wavy.atoms[i].z)/points };
                double dr = distance / points;
                vec lapl(points, 0.0);
                vec pos = { wavy.atoms[i].x, wavy.atoms[i].y, wavy.atoms[i].z };
#pragma omp parallel for
                for (int k = 0; k < points; k++) {
                    std::array<double, 3> t_pos = { pos[0], pos[1], pos[2] };
                    t_pos[0] += k*bond_vec[0];
                    t_pos[1] += k*bond_vec[1];
                    t_pos[2] += k*bond_vec[2];
                    lapl[k] = wavy.computeLap(t_pos);
                }
                std::string outname(wfn_name + "_bondwise_laplacian_" + std::to_string(i) + "_" + std::to_string(j) + ".dat");
                join_path(path, outname);
                std::ofstream result(path, std::ios::out);
                for (int k = 0; k < points; k++) {
                    result << std::setw(10) << std::scientific << std::setprecision(6) << dr*k << " " << std::setw(10) << std::scientific << std::setprecision(6) << lapl[k] << std::endl;
                }
                result.flush();
                result.close();
            }
        }
    }
}

void calc_partition_densities() {
    using namespace std;
    WFN Hartree_Fock(9);
    WFN DFT(9);
    Hartree_Fock.read_known_wavefunction_format("HF.gbw", std::cout, false);
    DFT.read_known_wavefunction_format("DFT.gbw", std::cout, false);

    vec Pos_C = { Hartree_Fock.atoms[0].x, Hartree_Fock.atoms[0].y, Hartree_Fock.atoms[0].z };
    vec Pos_H1 = { Hartree_Fock.atoms[1].x, Hartree_Fock.atoms[1].y, Hartree_Fock.atoms[1].z };
    vec Pos_H2 = { Hartree_Fock.atoms[2].x, Hartree_Fock.atoms[2].y, Hartree_Fock.atoms[2].z };
    vec Pos_H3 = { Hartree_Fock.atoms[3].x, Hartree_Fock.atoms[3].y, Hartree_Fock.atoms[3].z };
    vec Pos_H4 = { Hartree_Fock.atoms[4].x, Hartree_Fock.atoms[4].y, Hartree_Fock.atoms[4].z };
    vec x_coords = {
        Hartree_Fock.atoms[0].x,
        Hartree_Fock.atoms[1].x,
        Hartree_Fock.atoms[2].x,
        Hartree_Fock.atoms[3].x,
        Hartree_Fock.atoms[4].x
    };
    vec y_coords = {
        Hartree_Fock.atoms[0].y,
        Hartree_Fock.atoms[1].y,
        Hartree_Fock.atoms[2].y,
        Hartree_Fock.atoms[3].y,
        Hartree_Fock.atoms[4].y
    };
    vec z_coords = {
        Hartree_Fock.atoms[0].z,
        Hartree_Fock.atoms[1].z,
        Hartree_Fock.atoms[2].z,
        Hartree_Fock.atoms[3].z,
        Hartree_Fock.atoms[4].z
    };
    ivec charges{ 6,1,1,1,1 };
    const double fac = 0.01;
    const int min = -100, max = 500;
    const int size = -min + max + 1;
    vec C_dens (size, 0.0), H_dens(size, 0.0), total_dens(size, 0.0);
    vec HF_densities(size, 0.0);
    vec DFT_densities(size, 0.0);
    vec B_weights_C(size, 0.0), B_weights_H(size, 0.0);
    Thakkar C(6);
    Thakkar H(1);
    
    ofstream result("densities.dat", ios::out);
    ProgressBar pb(size, 100, "=", "", "Calculating densities");
//#pragma omp parallel 
//    {
        vec2 d(16);
        for (int i = 0; i < 16; i++) d[i].resize(DFT.get_ncen());
        vec phi(DFT.get_nmo(), 0.0);
        double x = 0;
        vec pa(5);
//#pragma omp for
        for (int i = 0; i < size; i++) {
            x = (i+min) * fac;
            HF_densities[i] = Hartree_Fock.compute_dens(x, 0, 0, d, phi);
            DFT_densities[i] = DFT.compute_dens(x, 0, 0, d, phi);
            double temp = abs(x - Pos_C[0]);
            C_dens[i] = C.get_radial_density(temp);
            total_dens[i] = C_dens[i];
            temp = abs(x - Pos_H1[0]);
            H_dens[i] = H.get_radial_density(temp);
            total_dens[i] += H_dens[i];
            temp = sqrt(pow(x - Pos_H2[0], 2) + pow(Pos_H2[1], 2) + pow(Pos_H2[2], 2));
            total_dens[i] += H.get_radial_density(temp);
            temp = sqrt(pow(x - Pos_H3[0], 2) + pow(Pos_H3[1], 2) + pow(Pos_H3[2], 2));
            total_dens[i] += H.get_radial_density(temp);
            temp = sqrt(pow(x - Pos_H4[0], 2) + pow(Pos_H4[1], 2) + pow(Pos_H4[2], 2));
            total_dens[i] += H.get_radial_density(temp);
            B_weights_C[i] = get_becke_w(5, charges.data(), x_coords.data(), y_coords.data(), z_coords.data(), 0, x, 0, 0, pa);
            B_weights_H[i] = get_becke_w(5, charges.data(), x_coords.data(), y_coords.data(), z_coords.data(), 1, x, 0, 0, pa);
            pb.update();
        }
//    }
    for (int i = 0; i < size; i++) {
        result << setw(10) << setprecision(4) << scientific << (i+min) * fac
            << setw(16) << setprecision(8) << scientific << HF_densities[i]
            << setw(16) << setprecision(8) << scientific << DFT_densities[i]
            << setw(16) << setprecision(8) << scientific << C_dens[i]
            << setw(16) << setprecision(8) << scientific << H_dens[i]
            << setw(16) << setprecision(8) << scientific << total_dens[i] 
            << setw(16) << setprecision(8) << scientific << B_weights_C[i]
            << setw(16) << setprecision(8) << scientific << B_weights_H[i] << endl;
    }
    result.flush();
    result.close();
    exit(0);
};

cdouble calc_spherically_averaged_at_k(vec2 &d1, vec2 &d2, vec2 &d3, vec2 &dens,
                                       const double &k,
                                       const double rel_precision = 1E-4,
                                       const double start_angle = 5.0,
                                       const bool print = false)
{
    cdouble old_result = 1E100;
    cdouble new_result = 0.9E100;
    double angle = start_angle;
    double ratio = abs(old_result / new_result - 1.0);
    double phi_ratio = 2.0;
    double rho, work;
    double *d1_local, *d2_local, *d3_local, *dens_local;
    int pmax;

    while ((ratio > rel_precision || phi_ratio > rel_precision) && angle / constants::PI < 1E5)
    {
        const double da = constants::PI / angle;
        const double da2 = da * da;
        // Make angular grids
        double x0, y0, z0, st0;
        cdouble sf_local;
        for (double phi = 0; phi <= constants::TWO_PI; phi += da)
        {
            const double cp = cos(phi);
            const double sp = sin(phi);
            for (double theta = da; theta <= constants::PI; theta += da)
            {
                const double st = sin(theta);
                st0 = st * da2;
                x0 = st * cp, y0 = st * sp, z0 = cos(theta);

                pmax = (int)dens[0].size();
                dens_local = dens[0].data();
                d1_local = d1[0].data();
                d2_local = d2[0].data();
                d3_local = d3[0].data();
                for (int p = pmax - 1; p >= 0; p--)
                {
                    rho = dens_local[p];
                    work = k * x0 * d1_local[p] + k * y0 * d2_local[p] + k * z0 * d3_local[p];
                    sf_local += cdouble(rho * cos(work), rho * sin(work)) * st0;
                }
            }
        }
        if (sf_local != 0.0)
        {
            old_result = new_result;
            new_result = sf_local / cdouble(constants::FOUR_PI);
            phi_ratio = std::abs(std::arg(old_result) - std::arg(new_result));
            // if(print)
            //     std::cout << "Phi ratio: " << phi_ratio << " " << std::arg(old_result) << " " << std::arg(new_result) << std::endl;
            ratio = std::abs(std::abs(old_result) / std::abs(new_result) - 1.0);
            // if(print)
            //     std::cout << "Ratio:     " << ratio << " " << std::abs(old_result) << " " << std::abs(new_result) << std::endl;
        }
        angle *= 1.2;
    }
    if (print)
        std::cout << "Done with " << std::setw(10) << std::setprecision(5) << k << " " << (ratio > rel_precision) << " " << angle << std::endl;
    return new_result;
}

void spherically_averaged_density(options &opt, const ivec val_els_alpha, const ivec val_els_beta)
{
    std::cout << "Reading wavefunction" << std::endl;
    using namespace std;
    WFN wavy(9);
    wavy.set_charge(opt.charge);
    wavy.set_multi(opt.mult);
    wavy.read_known_wavefunction_format(opt.wfn, cout, opt.debug);
    cout << "Number of MOs before: " << wavy.get_nmo() << endl;
    wavy.delete_unoccupied_MOs();
    cout << "Number of occupied MOs before: " << wavy.get_nmo() << endl;
    bvec MOs_to_delete(wavy.get_nmo(), false);
    int deleted = 0;
    if (val_els_alpha.size() > 0)
    {
        // Delete core orbitals
        int offset = wavy.get_MO_op_count(0);
        for (int i = offset - 1; i >= 0; i--)
            // only delete if i is not an element of val_els
            if (find(val_els_alpha.begin(), val_els_alpha.end(), i) == val_els_alpha.end())
            {
                cout << "Deleting from Alpha: " << i << endl;
                wavy.delete_MO(i);
                MOs_to_delete[i] = true;
                deleted++;
            }
        offset = wavy.get_MO_op_count(0);
        for (int i = wavy.get_nmo() - 1; i >= offset; i--)
            if (find(val_els_beta.begin(), val_els_beta.end(), i - offset) == val_els_beta.end())
            {
                cout << "Deleting from Beta: " << i - offset << endl;
                wavy.delete_MO(i);
                MOs_to_delete[i + deleted] = true;
                deleted++;
            }
    }
    cout << "MOs deleted: " << deleted << endl;
    cout << "MO map:" << endl;
    for (int i = 0; i < MOs_to_delete.size(); i++)
        cout << i << " " << MOs_to_delete[i] << endl;
    cout << "Number of MOs after: " << wavy.get_nmo() << endl;
    cout << "\n\nEnergies after:" << endl;
    for (int i = 0; i < wavy.get_nmo(); i++)
        cout << wavy.get_MO_energy(i) << " " << wavy.get_MO_occ(i) << endl;

    // Make radial grids on logarithmic scale
    const double dr = 0.00008;
    const long long int upper_r = 250000;
    // Calcualte density on angular grid at each point of radial grid, average and integrate
    long double tot_int2 = 0;
    vec radial_dens2(upper_r, 0.0);
    ProgressBar * progress = new ProgressBar(upper_r, 85, "=", " ", "Calculating Densities");
    cout << endl;

#pragma omp parallel for reduction(+ : tot_int2) num_threads(opt.threads)
    for (long long int _r = 1; _r < upper_r; _r++)
    {
        double r = _r * dr;
        radial_dens2[_r] = calc_grid_averaged_at_r(wavy, r, 1200, 5800, opt.debug);
        if (_r >= 1) {
            tot_int2 += radial_dens2[_r] * r * r * (r - (_r - 1) * dr);
        }
        progress->update();
    }
    delete (progress);
    cout << "Start writing the file" << endl;
    string el = constants::atnr2letter(wavy.get_atom_charge(0));
    ofstream out(el + ".dat", ios::out);
    out << "Total Integral: " << setw(18) << scientific << setprecision(10) << tot_int2 << "\n";
    for (int i = 0; i < upper_r; i++)
        out << setw(24) << scientific << setprecision(15) << i * dr << setw(24) << scientific << setprecision(16) << radial_dens2[i] / constants::FOUR_PI << "\n";
    out.flush();
    out.close();
}

void add_ECP_contribution_test(const ivec &asym_atom_list,
                               const WFN &wave,
                               cvec2 &sf,
                               vec2 &k_pt)
{
    double k = 1.0;
    // Using a Thakkar core density
    std::vector<Thakkar> temp;
    for (int i = 0; i < asym_atom_list.size(); i++)
    {
        temp.push_back(Thakkar(wave.atoms[asym_atom_list[i]].charge));
    }

#pragma omp parallel for private(k)
    for (int s = 0; s < sf[0].size(); s++)
    {
        k = k_pt[3][s];
        for (int i = 0; i < asym_atom_list.size(); i++)
        {
            if (wave.atoms[asym_atom_list[i]].ECP_electrons != 0)
                sf[i][s] += temp[i].get_core_form_factor(k, wave.atoms[asym_atom_list[i]].ECP_electrons);
        }
    }
}

void spherical_harmonic_test()
{
    const double phi = 0.3, theta = 0.4;
    for (int lam = 0; lam <= 5; lam++)
    {
        for (int m = -lam; m <= lam; m++)
        {
            vec d{std::sin(theta) * std::cos(phi),
                  std::sin(theta) * std::sin(phi),
                  std::cos(theta), 1.0, 1.0};
            std::cout << constants::spherical_harmonic(lam, m, d.data()) << " ";
        }
        std::cout << "\n";
    }
};


// Only valid for one atom positioned at 0,0,0
double compute_MO_spherical_orig(double x, double y, double z, double expon, double coef, int type)
{
    int l[3]{0, 0, 0};
    double ex = 0;
    double temp = 0;

    // x, y, z and dsqd
    vec d{
        x,
        y,
        z,
        x * x + y * y + z * z};

    // if (iat != atom) continue;
    constants::type2vector(type, l);
    temp = -expon * d[3];
    if (temp < -46.0517) // corresponds to cutoff of ex ~< 1E-20
        return 0.;
    ex = exp(temp);
    for (int k = 0; k < 3; k++)
    {
        if (l[k] == 0)
            continue;
        else if (l[k] == 1)
            ex *= d[k];
        else if (l[k] == 2)
            ex *= d[k] * d[k];
        else if (l[k] == 3)
            ex *= pow(d[k], 3);
        else if (l[k] == 4)
            ex *= pow(d[k], 4);
        else if (l[k] == 5)
            ex *= pow(d[k], 5);
    }
    return coef * ex; // build MO values at this point
}

// Only valid for one atom positioned at 0,0,0
double compute_MO_spherical_new(double x, double y, double z, double expon, double coef, int type, int m = 1)
{
    double result = 0.0;
    int l[3]{0, 0, 0};
    int type_local = type;
    double N = 1;
    // convert "type" to l
    if (type == 1)
    {
        type_local = 0;
        N = 1 / sqrt(1 / (constants::FOUR_PI));
    }
    else if (type > 1)
    {
        type_local = 1;
        N = 1 / sqrt(3 / (constants::FOUR_PI));
    }

    // r, theta, phi
    vec spher = constants::cartesian_to_spherical(x, y, z);
    double R = exp(-expon * spher[0] * spher[0]) * pow(spher[0], type_local);
    double Y = constants::real_spherical(type_local, m, spher[1], spher[2]) * N;
    result = coef * R * Y;

    double test = compute_MO_spherical_orig(x, y, z, expon, coef, type);
    std::cout << "New: " << result << " Old: " << test << " diff: " << result - test << " fact: " << result / test << std::endl;
    return result;
}

void Calc_MO_spherical_harmonics(
    cube &CubeMO,
    WFN &wavy,
    std::ostream &file,
    bool nodate)
{
    using namespace std;
    time_point start = get_time();
    ProgressBar *progress = NULL;
    if (!nodate)
        progress = new ProgressBar(CubeMO.get_size(0), 60, "=", " ", "Calculating MOs");

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < CubeMO.get_size(0); i++)
    {
        for (int j = 0; j < CubeMO.get_size(1); j++)
            for (int k = 0; k < CubeMO.get_size(2); k++)
            {

                double PosGrid[3]{
                    i * CubeMO.get_vector(0, 0) + j * CubeMO.get_vector(0, 1) + k * CubeMO.get_vector(0, 2) + CubeMO.get_origin(0),
                    i * CubeMO.get_vector(1, 0) + j * CubeMO.get_vector(1, 1) + k * CubeMO.get_vector(1, 2) + CubeMO.get_origin(1),
                    i * CubeMO.get_vector(2, 0) + j * CubeMO.get_vector(2, 1) + k * CubeMO.get_vector(2, 2) + CubeMO.get_origin(2)};

                CubeMO.set_value(i, j, k, compute_MO_spherical_new(PosGrid[0], PosGrid[1], PosGrid[2], wavy.get_exponent(0), wavy.get_MO_coef(0, 0), wavy.get_type(0)));
            }
        if (!nodate)
            progress->update();
    }
    if (!nodate)
    {
        delete (progress);

        time_point end = get_time();
        if (get_sec(start, end) < 60)
            file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) << " s" << endl;
        else if (get_sec(start, end) < 3600)
            file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 60 << " m " << get_sec(start, end) % 60 << " s" << endl;
        else
            file << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 3600 << " h " << (get_sec(start, end) % 3600) / 60 << " m" << endl;
    }
};


void calc_rho_cube(WFN &dummy)
{
    using namespace std;
    double MinMax[6]{0, 0, 0, 0, 0, 0};
    int steps[3]{0, 0, 0};
    readxyzMinMax_fromWFN(dummy, MinMax, steps, 3., 0.025, true);
    cube CubeRho(steps[0], steps[1], steps[2], dummy.get_ncen(), true);
    dummy.delete_unoccupied_MOs();
    CubeRho.give_parent_wfn(dummy);
    cout << "Starting work..." << endl;

    for (int i = 0; i < 3; i++)
    {
        CubeRho.set_origin(i, MinMax[i]);
        CubeRho.set_vector(i, i, (MinMax[i + 3] - MinMax[i]) / steps[i]);
    }
    CubeRho.set_comment1("Calculated density using NoSpherA2");
    CubeRho.set_comment2("from " + dummy.get_path());
    CubeRho.path = get_basename_without_ending(dummy.get_path()) + "_rho.cube";

    time_point start = get_time();
    const int s1 = CubeRho.get_size(0), s2 = CubeRho.get_size(1), s3 = CubeRho.get_size(2), total_size = s1 * s2 * s3;
    ;
    cout << "Lets go into the loop! There is " << total_size << " points" << endl;
    ProgressBar* progress = new ProgressBar(total_size, 50, "=", " ", "Calculating Rho");

    vec v1{
        CubeRho.get_vector(0, 0),
        CubeRho.get_vector(1, 0),
        CubeRho.get_vector(2, 0)},
        v2{
            CubeRho.get_vector(0, 1),
            CubeRho.get_vector(1, 1),
            CubeRho.get_vector(2, 1)},
        v3{
            CubeRho.get_vector(0, 2),
            CubeRho.get_vector(1, 2),
            CubeRho.get_vector(2, 2)},
        orig{
            CubeRho.get_origin(0),
            CubeRho.get_origin(1),
            CubeRho.get_origin(2)};

#pragma omp parallel
    {
        vec2 d;
        vec phi(dummy.get_nmo(), 0.0);
        d.resize(16);
        for (int i = 0; i < 16; i++)
            d[i].resize(dummy.get_ncen(), 0.0);
#pragma omp for schedule(dynamic)
        for (int index = 0; index < total_size; index++)
        {
            int i = index / (s2 * s3);
            int j = (index / s3) % s2;
            int k = index % s3;

            vec PosGrid{
                i * v1[0] + j * v2[0] + k * v3[0] + orig[0],
                i * v1[1] + j * v2[1] + k * v3[1] + orig[1],
                i * v1[2] + j * v2[2] + k * v3[2] + orig[2]};

            CubeRho.set_value(i, j, k, dummy.compute_dens(PosGrid[0], PosGrid[1], PosGrid[2], d, phi));
            progress->update();
        }
    }
    delete (progress);

    using namespace std;
    time_point end = get_time();
    if (get_sec(start, end) < 60)
        std::cout << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) << " s" << endl;
    else if (get_sec(start, end) < 3600)
        std::cout << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 60 << " m " << get_sec(start, end) % 60 << " s" << endl;
    else
        std::cout << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 3600 << " h " << (get_sec(start, end) % 3600) / 60 << " m" << endl;
    CubeRho.calc_dv();
    std::cout << "Number of electrons: " << std::fixed << std::setprecision(4) << CubeRho.sum() << std::endl;

    std::string fn(get_basename_without_ending(dummy.get_path()) + "_rho.cube");
    CubeRho.write_file(fn, false);
};

double s_value(double *d)
{
    primitive p_0(0, 0, 3.5, 1.0);
    int m = 0;
    return pow(gaussian_radial(p_0, d[3]) * constants::spherical_harmonic(p_0.type, m, d), 2);
}

void cube_from_coef_npy(std::string &coef_fn, std::string &xyzfile)
{
    std::vector<unsigned long> shape{};
    bool fortran_order;
    vec data{};

    npy::LoadArrayFromNumpy(coef_fn, shape, fortran_order, data);
    WFN dummy(7);
    dummy.read_xyz(xyzfile, std::cout);

    // const std::vector<std::vector<primitive>> basis(TZVP_JKfit.begin(), TZVP_JKfit.end());

    const int nr_coefs = load_basis_into_WFN(dummy, def2_qzvppd_rifit);
    std::cout << data.size() << " vs. " << nr_coefs << " ceofficients" << std::endl;
    calc_cube_ML(data, dummy);

    for (int i = 0; i < dummy.get_ncen(); i++)
        calc_cube_ML(data, dummy, i);
}

void test_xtb_molden(options &opt, std::ostream &log_file)
{
    using namespace std;
    for (int i = 0; i < 1; i++)
    {
        WFN wavy(8);
        wavy.read_molden("Co2.molden", cout, true);
        opt.cif = "Co2.cif";
        opt.dmin = 0.5;
        cout << "STARTING CALC" << endl;
        calculate_scattering_factors_HF(
            opt,
            wavy,
            log_file);
    }
}

void test_core_dens_corrected(double &precisison, int ncpus = 4, std::string ele = "Au", ivec val_els_alpha = {}, ivec val_els_beta = {})
{
    if (ncpus == -1)
        ncpus = 4;
    using namespace std;
#ifdef _OPENMP
    omp_set_num_threads(ncpus);
#endif
    vec2 res(6);
    for (int i = 0; i < res.size(); i++)
        res[i].resize(10000, 0.0);

    string dat = "core_dens_" + ele + ".dat";
    int el_nr = constants::get_Z_from_label(ele.c_str()) + 1;
    Thakkar T_Au(el_nr);

    WFN ECP_way_Au(9);
    string def2 = ele + "_def2TZVP.gbw";
    ECP_way_Au.read_gbw(def2, cout, false, false);
    ECP_way_Au.delete_unoccupied_MOs();
    ECP_way_Au.set_has_ECPs(true, true, 1);

    string jorge = ele + "_jorge.gbw";
    WFN wavy_full_Au(9);
    wavy_full_Au.read_gbw(jorge, cout, false);
    wavy_full_Au.delete_unoccupied_MOs();
    WFN wavy_val_Au(9);
    wavy_val_Au.read_gbw(jorge, cout, false);
    wavy_val_Au.delete_unoccupied_MOs();
    cout << "Number of occupied MOs before: " << wavy_val_Au.get_nmo() << endl;
    bvec MOs_to_delete(wavy_val_Au.get_nmo(), false);
    int deleted = 0;
    if (val_els_alpha.size() > 0)
    {
        // Delete core orbitals
        int offset = wavy_val_Au.get_MO_op_count(0);
        for (int i = offset - 1; i >= 0; i--)
            // only delete if i is not an element of val_els
            if (find(val_els_alpha.begin(), val_els_alpha.end(), i) == val_els_alpha.end())
            {
                cout << "Deleting from Alpha: " << i << endl;
                wavy_val_Au.delete_MO(i);
                MOs_to_delete[i] = true;
                deleted++;
            }
        offset = wavy_val_Au.get_MO_op_count(0);
        for (int i = wavy_val_Au.get_nmo() - 1; i >= offset; i--)
            if (find(val_els_beta.begin(), val_els_beta.end(), i - offset) == val_els_beta.end())
            {
                cout << "Deleting from Beta: " << i - offset << endl;
                wavy_val_Au.delete_MO(i);
                MOs_to_delete[i + deleted] = true;
            }
    }
    cout << "MOs deleted: " << deleted << endl;
    cout << "MO map:" << endl;
    for (int i = 0; i < MOs_to_delete.size(); i++)
        cout << i << " " << MOs_to_delete[i] << endl;
    cout << "Number of MOs after: " << wavy_val_Au.get_nmo() << endl;
    cout << "\n\nEnergies / Occu after:" << endl;
    for (int i = 0; i < wavy_val_Au.get_nmo(); i++)
        cout << wavy_val_Au.get_MO_energy(i) << " / " << wavy_val_Au.get_MO_occ(i) << endl;

    time_point start = get_time();

    const int upper = static_cast<int>(res[0].size());
    ProgressBar* progress = new ProgressBar(upper, 60, "=", " ", "Calculating Densities");
    
#pragma omp parallel for schedule(dynamic)
    for (int i = 1; i < upper; i++)
    {
        double sr = i * 0.001;
        res[0][i] = sr;
        res[1][i] = T_Au.get_core_density(sr, ECP_way_Au.atoms[0].ECP_electrons);
        res[2][i] = T_Au.get_radial_density(sr);
        res[3][i] = calc_spherically_averaged_at_r(ECP_way_Au, sr, precisison, 60);
        res[4][i] = calc_spherically_averaged_at_r(wavy_full_Au, sr, precisison, 60);
        res[5][i] = calc_spherically_averaged_at_r(wavy_val_Au, sr, precisison, 60);
        progress->update();
    }
    delete (progress);
    time_point end = get_time();
    cout << "Time taken: " << round(get_sec(start, end) / 60) << " m " << get_sec(start, end) % 60 << " s " << get_msec(start, end) << " ms" << endl;
    ofstream dat_out(dat, ios::out);
    dat_out << scientific << setprecision(12) << setw(20);
    for (int i = 0; i < res[0].size(); i++)
    {
        for (int j = 0; j < 6; j++)
        {
            auto t = res[j][i];
            dat_out << t;
            dat_out << " ";
        }
        dat_out << "\n";
    }
    dat_out << flush;
    dat_out.close();
}

void test_core_sfac_corrected(double &precisison, int ncpus = 4, std::string ele = "Au", ivec val_els_alpha = {}, ivec val_els_beta = {})
{
    if (ncpus == -1)
        ncpus = 4;
    using namespace std;
#ifdef _OPENMP
    omp_set_num_threads(ncpus);
#endif
    cvec2 res(6);
    for (int i = 0; i < res.size(); i++)
        res[i].resize(1000, 0.0);

    string dat = "core_sfac_" + ele + ".dat";
    const int el_nr = constants::get_Z_from_label(ele.c_str()) + 1;
    Thakkar T_Au(el_nr);

    WFN ECP_way_Au(9);
    ECP_way_Au.read_gbw(ele + "_def2TZVP.gbw", cout, false, false);
    ECP_way_Au.delete_unoccupied_MOs();
    ECP_way_Au.set_has_ECPs(true, true, 1);

    WFN wavy_full_Au(9);
    wavy_full_Au.read_gbw(ele + "_jorge.gbw", cout, false);
    wavy_full_Au.delete_unoccupied_MOs();
    WFN wavy_val_Au(9);
    wavy_val_Au.read_gbw(ele + "_jorge.gbw", cout, false);
    wavy_val_Au.delete_unoccupied_MOs();
    cout << "Number of occupied MOs before: " << wavy_val_Au.get_nmo() << endl;
    bvec MOs_to_delete(wavy_val_Au.get_nmo(), false);
    int deleted = 0;
    if (val_els_alpha.size() > 0)
    {
        // Delete core orbitals
        int offset = wavy_val_Au.get_MO_op_count(0);
        for (int i = offset - 1; i >= 0; i--)
            // only delete if i is not an element of val_els
            if (find(val_els_alpha.begin(), val_els_alpha.end(), i) == val_els_alpha.end())
            {
                cout << "Deleting from Alpha: " << i << endl;
                wavy_val_Au.delete_MO(i);
                MOs_to_delete[i] = true;
                deleted++;
            }
        offset = wavy_val_Au.get_MO_op_count(0);
        for (int i = wavy_val_Au.get_nmo() - 1; i >= offset; i--)
            if (find(val_els_beta.begin(), val_els_beta.end(), i - offset) == val_els_beta.end())
            {
                cout << "Deleting from Beta: " << i - offset << endl;
                wavy_val_Au.delete_MO(i);
                MOs_to_delete[i + deleted] = true;
            }
    }
    cout << "MOs deleted: " << deleted << endl;
    cout << "MO map:" << endl;
    for (int i = 0; i < MOs_to_delete.size(); i++)
        cout << i << " " << MOs_to_delete[i] << endl;
    cout << "Number of MOs after: " << wavy_val_Au.get_nmo() << endl;
    cout << "\n\nEnergies / Occu after:" << endl;
    for (int i = 0; i < wavy_val_Au.get_nmo(); i++)
        cout << wavy_val_Au.get_MO_energy(i) << " / " << wavy_val_Au.get_MO_occ(i) << endl;

    std::vector<time_point> time_points;
    std::vector<std::string> time_descriptions;

    ivec atom_type_list({el_nr});
    ivec asym_atom_list({0});
    bvec needs_grid({true});
    vec2 d1_ECP, d2_ECP, d3_ECP, dens_ECP;
    vec2 d1_all, d2_all, d3_all, dens_all;
    vec2 d1_val, d2_val, d3_val, dens_val;
    svec labels({ele});

    auto temp_cell = cell();
    make_hirshfeld_grids(0,
                         3,
                         temp_cell,
                         ECP_way_Au,
                         atom_type_list,
                         asym_atom_list,
                         needs_grid,
                         d1_ECP, d2_ECP, d3_ECP, dens_ECP,
                         labels,
                         cout,
                         time_points,
                         time_descriptions);
    make_hirshfeld_grids(0,
                         3,
                         temp_cell,
                         wavy_full_Au,
                         atom_type_list,
                         asym_atom_list,
                         needs_grid,
                         d1_all, d2_all, d3_all, dens_all,
                         labels,
                         cout,
                         time_points,
                         time_descriptions);
    make_hirshfeld_grids(0,
                         3,
                         temp_cell,
                         wavy_val_Au,
                         atom_type_list,
                         asym_atom_list,
                         needs_grid,
                         d1_val, d2_val, d3_val, dens_val,
                         labels,
                         cout,
                         time_points,
                         time_descriptions);

    const int upper = static_cast<int>(res[0].size());
    ProgressBar* progress = new ProgressBar(upper, 60, "=", " ", "Calculating SFACs");
    
#pragma omp parallel for schedule(dynamic)
    for (int i = 1; i < upper; i++)
    {
        double sr = i * 0.01;
        res[0][i] = sr;
        res[1][i] = T_Au.get_core_form_factor(sr, ECP_way_Au.atoms[0].ECP_electrons);
        res[2][i] = T_Au.get_form_factor(sr);
        res[3][i] = calc_spherically_averaged_at_k(d1_ECP, d2_ECP, d3_ECP, dens_ECP, sr, precisison, 150.0);
        res[4][i] = calc_spherically_averaged_at_k(d1_all, d2_all, d3_all, dens_all, sr, precisison, 150.0);
        res[5][i] = calc_spherically_averaged_at_k(d1_val, d2_val, d3_val, dens_val, sr, precisison, 150.0);
        progress->update();
    }
    delete (progress);

    ofstream dat_out("core_sfac_" + ele + ".dat", ios::out);
    time_point end = get_time();
    cout << "Time taken: " << get_sec(time_points.front(), time_points.back()) << " s " << get_msec(time_points.front(), time_points.back()) << " ms" << endl;
    dat_out << scientific << setprecision(12) << setw(20);
    double t = 0;
    for (int i = 0; i < res[0].size(); i++)
    {
        for (int j = 0; j < 6; j++)
        {
            t = res[j][i].real();
            dat_out << t;
            dat_out << " ";
        }
        dat_out << "\n";
    }
    dat_out << flush;
    dat_out.close();
}

double calc_pot_by_integral(vec3 &grid, const double &r, const double &cube_dist, const double &dr)
{
    double res = 0;
    const double dr3 = dr * dr * dr;
    double d = 0;
    for (double x = -cube_dist / dr; x <= cube_dist / dr; x++)
    {
        for (double y = -cube_dist / dr; y <= cube_dist / dr; y++)
        {
            for (double z = -cube_dist / dr; z <= cube_dist / dr; z++)
            {
                d = sqrt((x - r) * (x - r) + y * y + z * z);
                int offset = static_cast<int>(cube_dist / dr);
                res += grid[static_cast<int>(x) + offset][static_cast<int>(y) + offset][static_cast<int>(z) + offset] / d * dr3;
            }
        }
    }
    return res / constants::FOUR_PI;
}

void test_openblas()
{
#if has_RAS
    SALTEDPredictor SP;
    SP.load_BLAS();
    std::cout << "Running Openblas test" << std::endl;
    _test_openblas();
    SP.unload_BLAS();
    exit(0);
#else
    err_not_impl_f("Openblas not included in Executable!", std::cout);
#endif
}

void test_analytical_fourier()
{
    // Generate grid and k_pts
    vec2 kpts;
    for (int i = 1; i < 100; i++)
    {
        kpts.push_back({0.01 * i, 0, 0});                   // X
        kpts.push_back({0, 0.01 * i, 0});                   // Y
        kpts.push_back({0, 0, 0.01 * i});                   // Z
        kpts.push_back({0.01 * i, 0.01 * i, 0});            // XY
        kpts.push_back({0.01 * i, 0, 0.01 * i});            // XZ
        kpts.push_back({0, 0.01 * i, 0.01 * i});            // YZ
        kpts.push_back({0.01 * i * 2, 0.01 * i, 0.01 * i}); // XYZ
        kpts.push_back({-0.01 * i, 0, 0});
        kpts.push_back({0, -0.01 * i, 0});
        kpts.push_back({0, 0, -0.01 * i});
        kpts.push_back({-0.01 * i, -0.01 * i, 0});
        kpts.push_back({-0.01 * i, 0, -0.01 * i});
        kpts.push_back({0, -0.01 * i, -0.01 * i});
        kpts.push_back({-0.01 * i * 2, -0.01 * i, -0.01 * i});
    }
    vec2 grid;
    grid.resize(5); // x, y, z, dens, atomic_weight

    double alpha_min[] = {0.5};
    AtomGrid griddy(1E-25,
                    350,
                    770,
                    1,
                    4.5,
                    1,
                    alpha_min,
                    std::cout);

    double pos[] = {0};
    for (int i = 0; i < grid.size(); i++)
    {
        grid[i].resize(griddy.get_num_grid_points(), 0.0);
    }
    griddy.get_atomic_grid(0, pos, pos, pos, grid[0].data(), grid[1].data(), grid[2].data(), grid[4].data());

    // Initialize the vectors sf_A and sf_N
    cvec2 sf_A, sf_N;
    sf_A.resize(1);
    sf_N.resize(1);
    sf_A[0].resize(kpts.size(), 0.0);
    sf_N[0].resize(kpts.size(), 0.0);

    // Conditions for the Wavefunction
    const double c_exp = 2.0;
    double vals[] = {1.0};

    // double MinMax[6]{ 0, 0, 0, 0, 0, 0 };
    // int steps[3]{ 0, 0, 0 };
    // readxyzMinMax_fromWFN(wavy, MinMax, steps, 3., 0.025, true);
    // cube CubeMO(steps[0], steps[1], steps[2], 1, true);
    // CubeMO.give_parent_wfn(wavy);
    // for (int i = 0; i < 3; i++)
    //{
    //     CubeMO.set_origin(i, MinMax[i]);
    //     CubeMO.set_vector(i, i, (MinMax[i + 3] - MinMax[i]) / steps[i]);
    // }
    // Calc_MO_spherical_harmonics(CubeMO, wavy, 18, std::cout, false);
    // CubeMO.path = "MO.cube";
    // CubeMO.write_file(true);

    bool correct = true;

    cdouble max_diff, diff;
    for (int type = 6; type < 7; type++)
    {
        std::cout << "Testing l = " << type << std::endl;
        vec coefs(type * 2 + 1);

        // Initialize the Wavefunction
        WFN wavy(0);
        wavy.push_back_MO(0, 1.0, -13);
        wavy.push_back_atom("H", 0, 0, 0, 1);
        wavy.atoms[0].push_back_basis_set(c_exp, vals[0], type, 0);
        primitive p(1, type, c_exp, vals[0]);

        for (int l = 0; l < type * 2 + 1; l++)
        {
            for (int i = 0; i < coefs.size(); i++)
            {
                coefs[i] = 0.0;
            }
            max_diff = 0.0;
            coefs[l] = 1.0;
            // print content of coefs
            std::cout << "Coefs: ";
            for (int i = 0; i < coefs.size(); i++)
            {
                std::cout << coefs[i] << " ";
            }
            std::cout << "    |  m: " << l - type << std::endl;

            for (int i = 0; i < grid[0].size(); i++)
            {
                // grid[3][i] = wavy.compute_dens(grid[0][i], grid[1][i], grid[2][i]);
                grid[3][i] = calc_density_ML(grid[0][i], grid[1][i], grid[2][i], coefs, wavy.atoms);
            }

            // Empty the vectors sf:A nad sf_N
            for (int i = 0; i < kpts.size(); i++)
            {
                sf_A[0][i] = 0.0;
                sf_N[0][i] = 0.0;
            }

#pragma omp parallel for
            for (int i = 0; i < kpts.size(); i++)
            {
                double k_pt_local[3];
                for (int d = 0; d < 3; d++)
                {
                    k_pt_local[d] = kpts[i][d] * 2 * constants::PI;
                }
                sf_A[0][i] = sfac_bessel(p, k_pt_local, coefs);
                for (int _p = 0; _p < grid[0].size(); _p++)
                {
                    double work = 2 * constants::PI * (kpts[i][0] * grid[0][_p] + kpts[i][1] * grid[1][_p] + kpts[i][2] * grid[2][_p]);
                    sf_N[0][i] += std::polar(grid[3][_p] * grid[4][_p], work);
                }
                diff = abs(sf_A[0][i] - sf_N[0][i]);
                if (abs(diff) > abs(max_diff))
                {
                    max_diff = diff;
                }
                if (abs(diff) > 2.1E-4)
                {
                    correct = false;
                }
            }
            if (!correct)
            {
                std::cout << "Error at m: " << l - type << "   Max diff: " << std::setprecision(6) << max_diff << std::endl;
                break;
            }
            else
            {
                std::cout << "m: " << l - type << " passed!" << " Max diff: " << std::setprecision(6) << max_diff << std::endl;
            }
        }
        if (!correct)
            break;
        std::cout << "l = " << type << " passed!\n"
                  << std::endl;
    }
    if (!correct)
    {
        using namespace std;
        ofstream result("sfacs.dat", ios::out);
        for (int i = 0; i < kpts.size(); i++)
        {
            result << setw(8) << setprecision(2) << fixed << kpts[i][0];
            result << setw(8) << setprecision(2) << fixed << kpts[i][1];
            result << setw(8) << setprecision(2) << fixed << kpts[i][2];
            result << setw(16) << setprecision(8) << scientific << sf_A[0][i].real();
            result << setw(16) << setprecision(8) << scientific << sf_A[0][i].imag();
            result << setw(16) << setprecision(8) << scientific << sf_N[0][i].real();
            result << setw(16) << setprecision(8) << scientific << sf_N[0][i].imag();
            result << setw(16) << setprecision(8) << scientific << abs(sf_A[0][i] - sf_N[0][i]);
            result << setw(35) << setprecision(8) << scientific << sf_A[0][i] / sf_N[0][i];
            result << "\n";
        }
        result.flush();
        result.close();
    }
    else
    {
        std::cout << "All tests passed!" << std::endl;
    }
};

void draw_orbital(const int lambda, const int m, const double resulution = 0.025, const double radius = 3.5)
{
    if (m > lambda || m < -lambda)
    {
        std::cout << "m must be between -l and l" << std::endl;
        return;
    }

    // Initialize the Wavefunction
    WFN wavy(0);
    wavy.push_back_MO(0, 1.0, -13);
    wavy.push_back_atom("H", 0, 0, 0, 1);
    wavy.atoms[0].push_back_basis_set(1.0, 1.0, lambda, 0);
    double MinMax[6]{0, 0, 0, 0, 0, 0};
    int steps[3]{0, 0, 0};
    readxyzMinMax_fromWFN(wavy, MinMax, steps, radius, resulution, true);
    cube CubeMO(steps[0], steps[1], steps[2], 1, true);
    CubeMO.give_parent_wfn(wavy);
    for (int i = 0; i < 3; i++)
    {
        CubeMO.set_origin(i, MinMax[i]);
        CubeMO.set_vector(i, i, (MinMax[i + 3] - MinMax[i]) / steps[i]);
    }
#ifdef _OPENMP
    omp_lock_t l;
    omp_init_lock(&l);
#endif
    ProgressBar pb(CubeMO.get_size(0));

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < CubeMO.get_size(0); i++)
    {
        for (int j = 0; j < CubeMO.get_size(1); j++)
        {
            for (int k = 0; k < CubeMO.get_size(2); k++)
            {
                double PosGrid[3]{
                    i * CubeMO.get_vector(0, 0) + j * CubeMO.get_vector(0, 1) + k * CubeMO.get_vector(0, 2) + CubeMO.get_origin(0),
                    i * CubeMO.get_vector(1, 0) + j * CubeMO.get_vector(1, 1) + k * CubeMO.get_vector(1, 2) + CubeMO.get_origin(1),
                    i * CubeMO.get_vector(2, 0) + j * CubeMO.get_vector(2, 1) + k * CubeMO.get_vector(2, 2) + CubeMO.get_origin(2)};

                const double leng = sqrt(PosGrid[0] * PosGrid[0] + PosGrid[1] * PosGrid[1] + PosGrid[2] * PosGrid[2]);

                const std::pair<double, double> spherical = constants::norm_cartesian_to_spherical(PosGrid[0] / leng, PosGrid[1] / leng, PosGrid[2] / leng);
                const double radial = std::exp(-leng * leng) * std::pow(leng, lambda);                                  // Radial part of the wavefunction with b==1
                const double result = radial * constants::real_spherical(lambda, m, spherical.first, spherical.second); // Coefficients are 1
                CubeMO.set_value(i, j, k, result);
            }
        }
#ifdef _OPENMP
        omp_set_lock(&l);
#endif
        pb.update(std::cout);
#ifdef _OPENMP
        omp_unset_lock(&l);
#endif
    }
#ifdef _OPENMP
    omp_destroy_lock(&l);
#endif
    CubeMO.path = "Oribital-lam" + std::to_string(lambda) + "-m-" + std::to_string(m) + ".cube";
    CubeMO.write_file();
}

// const double dlm_function(const unsigned int& l, const int& m, const double& theta, const double& phi) {
//     double result = (double)NAN;
//     const double x = cos(theta);
//     const double s = -sin(theta);
//     const int m_ = abs(m);
//     // the formula for a real spherical harmonic is:
//     // for positive m:
//     // (-1)^m cos(m * phi)
//     // for negative m:
//     // (-1)^m sin(-m * phi)
//     const double p((m_) % 2 == 1 ? (m < 0 ? -sin(m_ * phi) : -cos(m_ * phi)) : (m < 0 ? sin(m_ * phi) : cos(m_ * phi)));
//     switch (l) {
//     case 0: switch (m_) {
//     case 0: result = 1.; break;
//     }; break;
//     case 1: switch (m_) {
//     case 0: result = x; break;
//     case 1: result = p * s; break;
//     }; break;
//     case 2: switch (m_) {
//     case 0: result = 0.5 * (3. * x * x - 1.); break;
//     case 1: result = p * 3. * x * s; break;
//     case 2: result = p * 3. * s * s; break;
//     }; break;
//     case 3: switch (m_) {
//     case 0: result = 0.5 * x * (5. * x * x - 3.); break;
//     case 1: result = p * 1.5 * (5. * x * x - 1.) * s; break;
//     case 2: result = p * 15. * x * s * s; break;
//     case 3: result = p * 15. * s * s * s; break;
//     }; break;
//     case 4: switch (m_) {
//     case 0: result = 0.125 * (35. * x * x * x * x - 30. * x * x + 3.); break;
//     case 1: result = p * 2.5 * (7. * x * x * x - 3. * x) * s; break;
//     case 2: result = p * 7.5 * (7. * x * x - 1.) * s * s; break;
//     case 3: result = p * 105. * x * s * s * s; break;
//     case 4: result = p * 105. * s * s * s * s; break;
//     }; break;
//     case 5: switch (m_) {
//     case 0: result = 0.125 * x * (63. * x * x * x * x - 70. * x * x + 15.); break;
//     case 1: result = p * 1.875 * (21. * x * x * x * x - 14. * x * x + 1.) * s; break;
//     case 2: result = p * 52.5 * x * (3. * x * x - 1.) * s * s; break;
//     case 3: result = p * 52.5 * (9. * x * x - 1.) * s * s * s; break;
//     case 4: result = p * 945. * x * s * s * s * s; break;
//     case 5: result = p * 945. * s * s * s * s * s; break;
//     }; break;
//     case 6: switch (m_) {
//     case 0: result = 0.0625 * (231. * x * x * x * x * x * x - 315. * x * x * x * x + 105. * x * x - 5.); break;
//     case 1: result = p * 2.625 * x * (33. * x * x * x * x - 30. * x * x + 5.) * s; break;
//     case 2: result = p * 13.125 * (33. * x * x * x * x - 18. * x * x + 1.) * s * s; break;
//     case 3: result = p * 157.5 * x * (11. * x * x - 3.) * s * s * s; break;
//     case 4: result = p * 472.5 * (11. * x * x - 1.) * s * s * s * s; break;
//     case 5: result = p * 10395. * x * s * s * s * s * s; break;
//     case 6: result = p * 10395. * s * s * s * s * s * s; break;
//     }; break;
//     case 7: switch (m_) {
//     case 0: result = 0.0625 * x * (429. * x * x * x * x * x * x - 693. * x * x * x * x + 315. * x * x - 35.); break;
//     case 1: result = p * 0.4375 * (429. * x * x * x * x * x * x - 495. * x * x * x * x + 135. * x * x - 5.) * s; break;
//     case 2: result = p * 7.875 * x * (143. * x * x * x * x - 110. * x * x + 15.) * s * s; break;
//     case 3: result = p * 39.375 * (143. * x * x * x * x - 66. * x * x + 3.) * s * s * s; break;
//     case 4: result = p * 1732.5 * x * (13. * x * x - 3.) * s * s * s * s; break;
//     case 5: result = p * 5197.5 * (13. * x * x - 1.) * s * s * s * s * s; break;
//     case 6: result = p * 135135. * x * s * s * s * s * s * s; break;
//     case 7: result = p * 135135. * s * s * s * s * s * s * s; break;
//     }; break;
//     case 8: switch (m_) {
//     case 0: result = 0.0078125 * (6435. * x * x * x * x * x * x * x * x - 12012. * x * x * x * x * x * x + 6930. * x * x * x * x - 1260. * x * x + 35.); break;
//     case 1: result = p * 0.5625 * x * (715. * x * x * x * x * x * x - 1001. * x * x * x * x + 385. * x * x - 35.) * s; break;
//     case 2: result = p * 19.6875 * (143. * x * x * x * x * x * x - 143. * x * x * x * x + 33. * x * x - 1.) * s * s; break;
//     case 3: result = p * 433.125 * x * (39. * x * x * x * x - 26. * x * x + 3.) * s * s * s; break;
//     case 4: result = p * 1299.375 * (65. * x * x * x * x - 26. * x * x + 1.) * s * s * s * s; break;
//     case 5: result = p * 67567.5 * x * (5. * x * x - 1.) * s * s * s * s * s; break;
//     case 6: result = p * 67567.5 * (15. * x * x - 1.) * s * s * s * s * s * s; break;
//     case 7: result = p * 2027025. * x * s * s * s * s * s * s * s; break;
//     case 8: result = p * 2027025. * s * s * s * s * s * s * s * s; break;
//     }; break;
//     };
//     return result;
// }

// TESTING THE FOURIER BESSEL TRANSFORM for WAVE FUNCTIONS
// cdouble S_n_recursion(int n, double H, double b);
// cdouble C_n_recursion(int n, double H, double b);
// For the case J_l(H) = int_0^inf j_l(Hr) * R_l(r)^2 * r^2 dr    |   Wave Functions!!
// cdouble S_0(double H, double b) {
//	using namespace std::complex_literals;
//     double two_32 = pow(2.0, 1.5);
//	return -(cerf((1.0i * H) / (two_32 * sqrt(b))) * constants::sqr_pi * 1.0i * exp(-H * H / (8. * b))) / (two_32 * sqrt(b));
// }
// double C_0(double H, double b) {
//	return  (constants::sqr_pi * exp(-H * H / (8. * b))) / (pow(2., 1.5) * sqrt(b));
// }
////Following 1/(4b) * ((n-1)C_(n-2) - H*S_(n-1)) = C_n
// cdouble C_n_recursion(int n, double H, double b) {
//     using namespace std::complex_literals;
//     if (n == 0) {
//         return C_0(H, b);
//     }
//     else if (n == 1) {
//         return (1. / (4. * b)) * (1. + S_0(H, b));
//     }
//     else {
//         return (1. / (4. * b)) * ((n - 1.) * C_n_recursion(n - 2, H, b) - H * S_n_recursion(n - 1, H, b));
//     }
// }
////\int_{ 0 }^ {\infty} r^ n\sin(Hr) \cdot e^ { -2br ^ 2 } dr & = \frac{ 1 }{4b}\left((n - 1)S_{ n - 2 } + HC_{ n - 1 }\right) = S_n
// cdouble S_n_recursion(int n, double H, double b) {
//     using namespace std::complex_literals;
//     if (n == 0) {
//         return S_0(H, b);
//     }
//     else if (n == 1) {
//         return (1. / (4. * b)) * H * C_0(H, b);
//     }
//     else {
//         return (1. / (4. * b)) * ((n - 1.) * S_n_recursion(n - 2, H, b) + H * C_n_recursion(n - 1, H, b));
//     }
// }
//// This function yields the fourier bessel transform of the radial integral of a gaussian density function (compare equation 1.2.7.9 in 10.1107/97809553602060000759),a ssuming that H = 2 \pi S
// cdouble fourier_bessel_integral(
//     primitive& p,
//     const double H
//)
//{
//     using namespace std::complex_literals;
//     const int l = p.type;
//     const double b = p.exp;
//     double N;
//     //double N = pow(
//     //    pow(2, 7 + 4 * l) * pow(b, 3 + 2 * l) / constants::PI / pow(doublefactorial(2 * l + 1), 2),
//     //    0.25);
//     //pow(8 * pow(b, 3) / pow(constants::PI, 3), 0.25);
//     //double N = p.norm_const;
//     //return N * (pow(H, l * 2) * constants::sqr_pi * exp(-H * H / (4 * b))) / (pow(2, l + 2) * pow(b, l + 1.5));
//     if (l == 0)
//     {
//         N = 1.;
//         //return N * N * (pow(H, l * 2) * constants::sqr_pi * exp(-H * H / (8 * b))) / (pow(2, l + 3.5) * pow(b, l + 1.5));
//         return  ((N * N) / (4 * b)) * C_n_recursion(0, H, b);
//     }
//     else if (l == 1)
//     {
//         N = 1.; //p.norm_const;
//         return ((N * N) / (H * H)) * (S_n_recursion(2, H, b) - H * C_n_recursion(3, H, b));
//     }
// }
// cdouble sfac_bessel(
//     const primitive& p,
//     const vec& k_point,
//     const vec& coefs
//)
//{
//     using namespace std::complex_literals;
//     vec local_k = k_point;
//     double leng = sqrt(local_k[0] * local_k[0] + local_k[1] * local_k[1] + local_k[2] * local_k[2]);
//     double H = 2 * constants::PI * leng;
//     //normalize the spherical harmonics k_point
//     for (int i = 0; i < 3; i++)
//         local_k[i] /= leng;
//
//     vec spherical = cartesian_to_spherical(local_k[0], local_k[1], local_k[2]);
//      for (int m = -p.type; m <= p.type; m++) {
//            cdouble angular = real_spherical(p.type, m, spherical[1], spherical[2]);
//            result += constants::FOUR_PI * pow(1.0i, p.type) * radial * angular * coefs[m + p.type];
//   cdouble radial = fourier_bessel_integral(p, H);
//   if (p.type == 0)
//   {
//       //constants::FOUR_PI^2 weil in angular^2 ein Faktor 1/4pi drin ist
//       return  constants::FOUR_PI * constants::FOUR_PI * pow(1.0i, p.type) * radial * p.coefficient * p.coefficient * abs(angular * angular);
//   }
//   else if (p.type == 1) {
////double angular = dlm_function(p.type, m+1, spherical[1], spherical[2]);
//      double angular2 = constants::spherical_harmonic(p.type, m, local_k.data());
//      return pow(0.75*constants::PI, 2) * constants::FOUR_PI * pow(1.0i, p.type) * radial * p.coefficient * p.coefficient * abs(angular * angular);
//  }
//}