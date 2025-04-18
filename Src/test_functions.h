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
#include "integrator.h"
#include "nos_math.h"
#include "rascaline.hpp"
#include "metatensor.h"
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
    auto t = new WFN(opt.wfn);
    wavy.push_back(*t);
    delete t;
    Thakkar O(wavy[0].get_atom_charge(0));
    Thakkar_Cation O_cat(wavy[0].get_atom_charge(0));
    Thakkar_Anion O_an(wavy[0].get_atom_charge(0));
    err_checkf(wavy[0].get_ncen() != 0, "No Atoms in the wavefunction, this will not work!! ABORTING!!", std::cout);
    err_checkf(exists(opt.cif), "CIF does not exists!", std::cout);
    // err_checkf(exists(asym_cif), "Asym/Wfn CIF does not exists!", file);

    std::vector<_time_point> time_points;
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
    ProgressBar *progress = new ProgressBar(smax, 50, "=", " ", "Calculating Scattering factors");
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
            if ((new_result < 1E-15 && angle > 36))
            {
#pragma omp critical
                std::cout << "Aborted due to too small density at r = " << r << " and angle = " << angle << std::endl;
                break;
            }
            if (new_result > 1E90)
            {
#pragma omp critical
                std::cout << "Aborted due too large density value at r = " << r << std::endl;
                break;
            }
            if (angle > 360)
            {
#pragma omp critical
                std::cout << "Aborted due to too large angle = " << angle << std::endl;
                break;
            }
        }
        else
        {
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

double calc_grid_averaged_at_r(const WFN &wavy,
                               const double &r,
                               const int &min_angular = 60,
                               const int &max_angular = 1000)
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

    for (int i = constants::get_angular_order(min_num_angular_points_closest); i < constants::get_angular_order(max_num_angular_points_closest) + 1; i++)
    {
        angular_off = i * constants::MAG;
        ls.ld_by_order(constants::lebedev_table[i],
                       angular_x.data() + angular_off,
                       angular_y.data() + angular_off,
                       angular_z.data() + angular_off,
                       angular_w.data() + angular_off);
    }

    const double rb = constants::bragg_angstrom[wavy.get_atom_charge(0)] / (5.0E10 * constants::a0);

    int num_angular = max_num_angular_points_closest;
    if (r < rb)
    {
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
    for (int iang = start; iang < size; iang++)
    {
        p = angular_off + iang;
        double x = angular_x[p] * r + wavy.get_atom_coordinate(0, 0);
        double y = angular_y[p] * r + wavy.get_atom_coordinate(0, 1);
        double z = angular_z[p] * r + wavy.get_atom_coordinate(0, 2);
        dens += wavy.compute_dens(x, y, z, d, _phi) * constants::FOUR_PI * angular_w[p];
    }
    // if (print)
    //     std::cout << "Done with " << std::setw(10) << std::setprecision(5) << r << std::endl;
    return dens;
}
/// not like above //////////////////////////////////////////////////////////////////////////-spherical_aver_hirsh in convience.cpp////////
double calc_hirsh_grid_averaged_at_r(const WFN &wavy,
                                     const int i, /// now counts for several atoms
                                     const double &r,
                                     const int &min_angular = 60,
                                     const int &max_angular = 1000,
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

    for (int j = constants::get_angular_order(min_num_angular_points_closest); j < constants::get_angular_order(max_num_angular_points_closest) + 1; j++)
    {
        angular_off = j * constants::MAG;
        ls.ld_by_order(constants::lebedev_table[j],
                       angular_x.data() + angular_off,
                       angular_y.data() + angular_off,
                       angular_z.data() + angular_off,
                       angular_w.data() + angular_off);
    }

    const double rb = constants::bragg_angstrom[wavy.get_atom_charge(i)] / (5.0E10 * constants::a0);

    int num_angular = max_num_angular_points_closest;
    if (r < rb)
    {
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
    Thakkar A(wavy.get_atom_charge(i)); /// for the atom i, a Thakkar object is created
#pragma omp parallel
    {
        vec2 d(16);
        vec _phi(wavy.get_nmo(), 0.0);
        for (int j = 0; j < 16; j++)
            d[j].resize(wavy.get_ncen(), 0.0);
#pragma omp for reduction(+ : dens)
        for (int iang = start; iang < size; iang++)
        {
            p = angular_off + iang;
            const double x = angular_x[p] * r + wavy.get_atom_coordinate(i, 0); /// moving the lebedev points to the atom i
            const double y = angular_y[p] * r + wavy.get_atom_coordinate(i, 1);
            const double z = angular_z[p] * r + wavy.get_atom_coordinate(i, 2);
            const std::array<double, 3> d_ = {angular_x[p] * r, angular_y[p] * r, angular_z[p] * r}; /// as the function get_radial_density needs a distance, the distance is calculated. The atom A is in the origin, the wavy.get_atom_coordinate(i,0) is 0
            const double dist = array_length(d_);
            const double rho_a = A.get_radial_density(dist); /// the radial density of the atom i is calculated
            double rho_all = rho_a;                          /// the molecular density based on pro-atoms
            for (int atom = 0; atom < wavy.get_ncen(); atom++)
            { /// a new for loop is started, which calculates the sum of rhos, with respect to the lebedev points
                if (atom == i)
                    continue; /// if the atom is the same as the atom i, the loop is skipped
                const std::array<double, 3> d_atom = {x - wavy.get_atom_coordinate(atom, 0), y - wavy.get_atom_coordinate(atom, 1), z - wavy.get_atom_coordinate(atom, 2)};
                const double dist_atom = array_length(d_atom); /// is like d_, but for another atom
                Thakkar thakkar_atom(wavy.get_atom_charge(atom));
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
double calc_grid_averaged_at_r_from_cube(const cube &cuby,
                                         const double &r,
                                         const int &min_angular = 60,
                                         const int &max_angular = 1000,
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

    for (int j = constants::get_angular_order(min_num_angular_points_closest); j < constants::get_angular_order(max_num_angular_points_closest) + 1; j++)
    {
        angular_off = j * constants::MAG;
        ls.ld_by_order(constants::lebedev_table[j],
                       angular_x.data() + angular_off,
                       angular_y.data() + angular_off,
                       angular_z.data() + angular_off,
                       angular_w.data() + angular_off);
    }

    const double rb = constants::bragg_angstrom[cuby.get_parent_wfn_atoms()[0].get_charge()] / (5.0E10 * constants::a0);

    int num_angular = max_num_angular_points_closest;
    if (r < rb)
    {
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

#pragma omp parallel
    {
        vec2 d(16);
        for (int j = 0; j < 16; j++)
            d[j].resize(cuby.get_na(), 0.0);
#pragma omp for reduction(+ : dens)
        for (int iang = start; iang < size; iang++)
        {
            p = angular_off + iang;
            const double x = angular_x[p] * r + cuby.get_parent_wfn_atoms()[0].get_coordinate(0); /// moving the lebedev points to the atom i
            const double y = angular_y[p] * r + cuby.get_parent_wfn_atoms()[0].get_coordinate(1);
            const double z = angular_z[p] * r + cuby.get_parent_wfn_atoms()[0].get_coordinate(2);
            dens += cuby.get_interpolated_value(x, y, z) * constants::FOUR_PI * angular_w[p];
        }
    }
    if (print)
        std::cout << "Done with " << std::setw(10) << std::setprecision(5) << r << std::endl;
    return dens;
}

double calc_fukui_averaged_at_r(const WFN &wavy1,
                                const WFN &wavy2,
                                const double &r,
                                const int &min_angular = 60,
                                const int &max_angular = 1000,
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

    for (int j = constants::get_angular_order(min_num_angular_points_closest); j < constants::get_angular_order(max_num_angular_points_closest) + 1; j++)
    {
        angular_off = j * constants::MAG;
        ls.ld_by_order(constants::lebedev_table[j],
                       angular_x.data() + angular_off,
                       angular_y.data() + angular_off,
                       angular_z.data() + angular_off,
                       angular_w.data() + angular_off);
    }

    const double rb = constants::bragg_angstrom[wavy1.get_atom_charge(0)] / (5.0E10 * constants::a0);

    int num_angular = max_num_angular_points_closest;
    if (r < rb)
    {
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

#pragma omp parallel
    {
        vec2 d(16);
        vec _phi(std::max(wavy1.get_nmo(), wavy2.get_nmo()), 0.0);
        for (int j = 0; j < 16; j++)
            d[j].resize(wavy1.get_ncen(), 0.0);
#pragma omp for reduction(+ : dens)
        for (int iang = start; iang < size; iang++)
        {
            p = angular_off + iang;
            const double x = angular_x[p] * r + wavy1.get_atom_coordinate(0, 0); /// moving the lebedev points to the atom i
            const double y = angular_y[p] * r + wavy1.get_atom_coordinate(0, 1);
            const double z = angular_z[p] * r + wavy1.get_atom_coordinate(0, 2);
            dens += (wavy1.compute_dens(x, y, z, d, _phi) - wavy2.compute_dens(x, y, z, d, _phi)) * constants::FOUR_PI * angular_w[p];
        }
    }
    if (print)
        std::cout << "Done with " << std::setw(10) << std::setprecision(5) << r << std::endl;
    return dens;
}

double subtract_dens_from_gbw(std::filesystem::path &wfn_name_1,
                              std::filesystem::path &wfn_name_2, const double &r, const double &resol)
{
    using namespace std;

    WFN wavy1(wfn_name_1);
    WFN wavy2(wfn_name_2); //  reading the WFNs

    double MinMax1[6];
    int steps1[3];
    readxyzMinMax_fromWFN(wavy1, MinMax1, steps1, r, resol, true);

    double MinMax2[6];
    int steps2[3];
    readxyzMinMax_fromWFN(wavy2, MinMax2, steps2, r, resol, true); // and getting the MinMax and steps

    double MinMax[6]{100, 100, 100, -100, -100, -100};
    int steps[3]{0, 0, 0};
    for (int i = 0; i < 3; i++)
    {
        if (MinMax1[i] < MinMax[i])
            MinMax[i] = MinMax1[i];
        if (MinMax1[i + 3] > MinMax[i + 3])
            MinMax[i + 3] = MinMax1[i + 3];
    }
    for (int i = 0; i < 3; i++)
    {
        if (MinMax2[i] < MinMax[i])
            MinMax[i] = MinMax2[i];
        if (MinMax2[i + 3] > MinMax[i + 3])
            MinMax[i + 3] = MinMax2[i + 3];
        steps[i] = (int)ceil(constants::bohr2ang(MinMax[i + 3] - MinMax[i]) / resol);
    }

    cube dens1(steps[0], steps[1], steps[2], wavy1.get_ncen(), true);
    cube dens2(steps[0], steps[1], steps[2], wavy2.get_ncen(), true);
    dens1.give_parent_wfn(wavy1);
    dens2.give_parent_wfn(wavy2);

    for (int i = 0; i < 3; i++)
    {
        dens1.set_origin(i, MinMax[i]);
        dens1.set_vector(i, i, (MinMax[i + 3] - MinMax[i]) / steps[i]);
        dens2.set_origin(i, MinMax[i]);
        dens2.set_vector(i, i, (MinMax[i + 3] - MinMax[i]) / steps[i]);
    }

    //   dens1.path = get_basename_without_ending(wavy1.get_path()) + "_rho.cube";
    // dens2.path = get_basename_without_ending(wavy2.get_path()) + "_rho.cube";
    dens1.set_path(wavy1.get_path().parent_path() / (wavy1.get_path().stem().string() + "_rho.cube"));
    dens2.set_path(wavy2.get_path().parent_path() / (wavy2.get_path().stem().string() + "_rho.cube"));

    Calc_Rho(dens1, wavy1, r, std::cout, false);
    Calc_Rho(dens2, wavy2, r, std::cout, false);
    dens1.write_file(true);
    dens2.write_file(true);

    cube difference = dens1 - dens2;

    for (int i = 0; i < 3; i++)
    {
        difference.set_origin(i, MinMax[i]);
        difference.set_vector(i, i, (MinMax[i + 3] - MinMax[i]) / steps[i]);
    }
    difference.give_parent_wfn(wavy1);
    difference.write_file("difference.cube", false);
    difference.calc_dv();
    dens1.calc_dv();
    dens2.calc_dv();
    std::cout << dens1.sum() << std::endl;
    std::cout << dens2.sum() << std::endl;
    std::cout << difference.sum() << std::endl;
    return 0;
}

void bondwise_laplacian_plots(std::filesystem::path &wfn_name)
{
    char cwd[1024];
#ifdef _WIN32
    if (_getcwd(cwd, sizeof(cwd)) != NULL)
#else
    if (getcwd(cwd, sizeof(cwd)) != NULL)
#endif
        std::cout << "Current working dir: " << cwd << std::endl;
    WFN wavy(wfn_name);

    err_checkf(wavy.get_ncen() != 0, "No Atoms in the wavefunction, this will not work!! ABORTING!!", std::cout);

    int points = 1001;

    for (int i = 0; i < wavy.get_ncen(); i++)
    {
        for (int j = i + 1; j < wavy.get_ncen(); j++)
        {
            std::filesystem::path path = cwd;
            double distance = sqrt(pow(wavy.get_atom_coordinate(i, 0) - wavy.get_atom_coordinate(j, 0), 2) + pow(wavy.get_atom_coordinate(i, 1) - wavy.get_atom_coordinate(j, 1), 2) + pow(wavy.get_atom_coordinate(i, 2) - wavy.get_atom_coordinate(j, 2), 2));
            double svdW = constants::ang2bohr(constants::covalent_radii[wavy.get_atom_charge(i)] + constants::covalent_radii[wavy.get_atom_charge(j)]);
            if (distance < 1.35 * svdW)
            {
                std::cout << "Bond between " << i << " (" << wavy.get_atom_charge(i) << ") and " << j << " (" << wavy.get_atom_charge(j) << ") with distance " << distance << " and svdW " << svdW << std::endl;
                const vec bond_vec = {(wavy.get_atom_coordinate(j, 0) - wavy.get_atom_coordinate(i, 0)) / points, (wavy.get_atom_coordinate(j, 1) - wavy.get_atom_coordinate(i, 1)) / points, (wavy.get_atom_coordinate(j, 2) - wavy.get_atom_coordinate(i, 2)) / points};
                double dr = distance / points;
                vec lapl(points, 0.0);
                const vec pos = {wavy.get_atom_coordinate(i, 0), wavy.get_atom_coordinate(i, 1), wavy.get_atom_coordinate(i, 2)};
#pragma omp parallel for
                for (int k = 0; k < points; k++)
                {
                    std::array<double, 3> t_pos = {pos[0], pos[1], pos[2]};
                    t_pos[0] += k * bond_vec[0];
                    t_pos[1] += k * bond_vec[1];
                    t_pos[2] += k * bond_vec[2];
                    lapl[k] = wavy.computeLap(t_pos);
                }
                std::filesystem::path outname(wfn_name.string() + "_bondwise_laplacian_" + std::to_string(i) + "_" + std::to_string(j) + ".dat");
                path = path / std::filesystem::path(outname);
                std::ofstream result(path, std::ios::out);
                for (int k = 0; k < points; k++)
                {
                    result << std::setw(10) << std::scientific << std::setprecision(6) << dr * k << " " << std::setw(10) << std::scientific << std::setprecision(6) << lapl[k] << std::endl;
                }
                result.flush();
                result.close();
            }
            else
            {
                std::cout << "No bond between " << i << " and " << j << " with distance " << distance << " and svdW " << svdW << std::endl;
            }
        }
    }
}

void calc_partition_densities()
{
    using namespace std;
    WFN Hartree_Fock("HF.gbw");
    WFN DFT("DFT.gbw");

    vec Pos_C = {Hartree_Fock.get_atom_coordinate(0, 0), Hartree_Fock.get_atom_coordinate(0, 1), Hartree_Fock.get_atom_coordinate(0, 2)};
    vec Pos_H1 = {Hartree_Fock.get_atom_coordinate(1, 0), Hartree_Fock.get_atom_coordinate(1, 1), Hartree_Fock.get_atom_coordinate(1, 2)};
    vec Pos_H2 = {Hartree_Fock.get_atom_coordinate(2, 0), Hartree_Fock.get_atom_coordinate(2, 1), Hartree_Fock.get_atom_coordinate(2, 2)};
    vec Pos_H3 = {Hartree_Fock.get_atom_coordinate(3, 0), Hartree_Fock.get_atom_coordinate(3, 1), Hartree_Fock.get_atom_coordinate(3, 2)};
    vec Pos_H4 = {Hartree_Fock.get_atom_coordinate(4, 0), Hartree_Fock.get_atom_coordinate(4, 1), Hartree_Fock.get_atom_coordinate(4, 2)};
    vec x_coords = {
        Hartree_Fock.get_atom_coordinate(0, 0),
        Hartree_Fock.get_atom_coordinate(1, 0),
        Hartree_Fock.get_atom_coordinate(2, 0),
        Hartree_Fock.get_atom_coordinate(3, 0),
        Hartree_Fock.get_atom_coordinate(4, 0)};
    vec y_coords = {
        Hartree_Fock.get_atom_coordinate(0, 1),
        Hartree_Fock.get_atom_coordinate(1, 1),
        Hartree_Fock.get_atom_coordinate(2, 1),
        Hartree_Fock.get_atom_coordinate(3, 1),
        Hartree_Fock.get_atom_coordinate(4, 1)};
    vec z_coords = {
        Hartree_Fock.get_atom_coordinate(0, 2),
        Hartree_Fock.get_atom_coordinate(1, 2),
        Hartree_Fock.get_atom_coordinate(2, 2),
        Hartree_Fock.get_atom_coordinate(3, 2),
        Hartree_Fock.get_atom_coordinate(4, 2)};
    ivec charges{6, 1, 1, 1, 1};
    const double fac = 0.01;
    const int min = -100, max = 500;
    const int size = -min + max + 1;
    vec C_dens(size, 0.0), H_dens(size, 0.0), total_dens(size, 0.0);
    vec HF_densities(size, 0.0);
    vec DFT_densities(size, 0.0);
    vec B_weights_C(size, 0.0), B_weights_H(size, 0.0);
    Thakkar C(6);
    Thakkar H(1);

    ofstream result("densities.dat", ios::out);
    ProgressBar pb(size, 100, "=", "", "Calculating densities");
    // #pragma omp parallel
    //     {
    vec2 d(16);
    for (int i = 0; i < 16; i++)
        d[i].resize(DFT.get_ncen());
    vec phi(DFT.get_nmo(), 0.0);
    double x = 0;
    vec pa(5);
    // #pragma omp for
    for (int i = 0; i < size; i++)
    {
        x = (i + min) * fac;
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
    for (int i = 0; i < size; i++)
    {
        result << setw(10) << setprecision(4) << scientific << (i + min) * fac
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
    WFN wavy(opt.wfn, opt.charge, opt.mult);
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
    ProgressBar *progress = new ProgressBar(upper_r, 85, "=", " ", "Calculating Densities");
    cout << endl;

#pragma omp parallel for reduction(+ : tot_int2) num_threads(opt.threads)
    for (long long int _r = 1; _r < upper_r; _r++)
    {
        double r = _r * dr;
        radial_dens2[_r] = calc_grid_averaged_at_r(wavy, r, 1200, 5800);
        if (_r >= 1)
        {
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
        temp.push_back(Thakkar(wave.get_atom_charge(asym_atom_list[i])));
    }

#pragma omp parallel for private(k)
    for (int s = 0; s < sf[0].size(); s++)
    {
        k = k_pt[3][s];
        for (int i = 0; i < asym_atom_list.size(); i++)
        {
            if (wave.get_atom_ECP_electrons(asym_atom_list[i]) != 0)
                sf[i][s] += temp[i].get_core_form_factor(k, wave.get_atom_ECP_electrons(asym_atom_list[i]));
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
    _time_point start = get_time();
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

        _time_point end = get_time();
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
    CubeRho.set_comment2("from " + dummy.get_path().string());
    CubeRho.set_path((dummy.get_path().parent_path() / dummy.get_path().stem()).string() + "_rho.cube");

    _time_point start = get_time();
    const int s1 = CubeRho.get_size(0), s2 = CubeRho.get_size(1), s3 = CubeRho.get_size(2), total_size = s1 * s2 * s3;
    ;
    cout << "Lets go into the loop! There is " << total_size << " points" << endl;
    ProgressBar *progress = new ProgressBar(total_size, 50, "=", " ", "Calculating Rho");

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
    _time_point end = get_time();
    if (get_sec(start, end) < 60)
        std::cout << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) << " s" << endl;
    else if (get_sec(start, end) < 3600)
        std::cout << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 60 << " m " << get_sec(start, end) % 60 << " s" << endl;
    else
        std::cout << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 3600 << " h " << (get_sec(start, end) % 3600) / 60 << " m" << endl;
    CubeRho.calc_dv();
    std::cout << "Number of electrons: " << std::fixed << std::setprecision(4) << CubeRho.sum() << std::endl;

    std::filesystem::path fn((dummy.get_path().parent_path() / dummy.get_path().stem()).string() + "_rho.cube");
    CubeRho.write_file(fn, false);
};

void cube_from_coef_npy(std::string &coef_fn, std::string &xyzfile)
{
    std::vector<unsigned long> shape{};
    bool fortran_order;
    vec data{};

    npy::LoadArrayFromNumpy(coef_fn, shape, fortran_order, data);
    WFN dummy(xyzfile);

    // const std::vector<std::vector<primitive>> basis(TZVP_JKfit.begin(), TZVP_JKfit.end());

    const int nr_coefs = load_basis_into_WFN(dummy, BasisSetLibrary().get_basis_set("cc-pvqz-jkfit"));
    std::cout << data.size() << " vs. " << nr_coefs << " ceofficients" << std::endl;
    cube res = calc_cube_ML(data, dummy);
    res.set_path((dummy.get_path().parent_path() / dummy.get_path().stem()).string() + "_PySCF_COEFS_rho.cube");
    res.write_file(true);
    // for (int i = 0; i < dummy.get_ncen(); i++)
    //     calc_cube_ML(data, dummy, i);
}

void test_xtb_molden(options &opt, std::ostream &log_file)
{
    using namespace std;
    for (int i = 0; i < 1; i++)
    {
        std::vector<WFN> wavy;
        wavy.emplace_back("Co2.molden");
        opt.cif = "Co2.cif";
        opt.dmin = 0.5;
        cout << "STARTING CALC" << endl;
        svec empty = {};
        calculate_scattering_factors(
            opt,
            wavy,
            log_file,
            empty,
            0);
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
    string def2 = ele + "_def2TZVP.gbw";
    WFN ECP_way_Au(def2);
    ECP_way_Au.delete_unoccupied_MOs();
    ECP_way_Au.set_has_ECPs(true, true, 1);

    string jorge = ele + "_jorge.gbw";
    WFN wavy_full_Au(jorge);
    wavy_full_Au.delete_unoccupied_MOs();
    WFN wavy_val_Au(jorge);
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

    _time_point start = get_time();

    const int upper = static_cast<int>(res[0].size());
    ProgressBar *progress = new ProgressBar(upper, 60, "=", " ", "Calculating Densities");

#pragma omp parallel for schedule(dynamic)
    for (int i = 1; i < upper; i++)
    {
        double sr = i * 0.001;
        res[0][i] = sr;
        res[1][i] = T_Au.get_core_density(sr, ECP_way_Au.get_atom_ECP_electrons(0));
        res[2][i] = T_Au.get_radial_density(sr);
        res[3][i] = calc_spherically_averaged_at_r(ECP_way_Au, sr, precisison, 60);
        res[4][i] = calc_spherically_averaged_at_r(wavy_full_Au, sr, precisison, 60);
        res[5][i] = calc_spherically_averaged_at_r(wavy_val_Au, sr, precisison, 60);
        progress->update();
    }
    delete (progress);
    _time_point end = get_time();
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

    WFN ECP_way_Au(ele + "_def2TZVP.gbw");
    ECP_way_Au.delete_unoccupied_MOs();
    ECP_way_Au.set_has_ECPs(true, true, 1);

    WFN wavy_full_Au(ele + "_jorge.gbw");
    wavy_full_Au.delete_unoccupied_MOs();
    WFN wavy_val_Au(ele + "_jorge.gbw");
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

    std::vector<_time_point> time_points;
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
    ProgressBar *progress = new ProgressBar(upper, 60, "=", " ", "Calculating SFACs");

#pragma omp parallel for schedule(dynamic)
    for (int i = 1; i < upper; i++)
    {
        double sr = i * 0.01;
        res[0][i] = sr;
        res[1][i] = T_Au.get_core_form_factor(sr, ECP_way_Au.get_atom_ECP_electrons(0));
        res[2][i] = T_Au.get_form_factor(sr);
        res[3][i] = calc_spherically_averaged_at_k(d1_ECP, d2_ECP, d3_ECP, dens_ECP, sr, precisison, 150.0);
        res[4][i] = calc_spherically_averaged_at_k(d1_all, d2_all, d3_all, dens_all, sr, precisison, 150.0);
        res[5][i] = calc_spherically_averaged_at_k(d1_val, d2_val, d3_val, dens_val, sr, precisison, 150.0);
        progress->update();
    }
    delete (progress);

    ofstream dat_out("core_sfac_" + ele + ".dat", ios::out);
    _time_point end = get_time();
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
    set_BLAS_threads(4);
    std::cout << "Running Openblas test" << std::endl;
    _test_openblas();
}

void test_analytical_fourier(bool full)
{
    // Generate grid and k_pts
    vec2 kpts;

    for (int i = 1; i < 1000; i++) {
        //Generate random k-points with values between -1 and 1
		kpts.push_back({ (double)rand() / RAND_MAX * 2 - 1, (double)rand() / RAND_MAX * 2 - 1, (double)rand() / RAND_MAX * 2 - 1 });
    }
    vec2 grid;
    grid.resize(5); // x, y, z, dens, atomic_weight

    // Conditions for the Wavefunction
    const double c_exp = 2.0;
    double vals[] = { 1.0 };
    unsigned int max_l = 6;
    double radial_res = 1E-17;
    if (full) {
        max_l = 8;
        radial_res = 1E-25;
    }

    double alpha_min[] = {0.5};
    AtomGrid griddy(radial_res,
        350,
        5000,
        1,
        c_exp,
        max_l,
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
    bool all_correct = true; //Veryfiy if all m for one l are correct break if one failed
    bool correct = true; //Verify if the current m is correct

    cdouble max_diff, diff;
    for (unsigned int type = 0; type <= max_l; type++)
    {
        std::cout << "Testing l = " << type << std::endl;
        vec coefs(type * 2 + 1);

        // Initialize the Wavefunction
        WFN wavy(0);
        wavy.push_back_MO(0, 1.0, -13);
        wavy.push_back_atom("H", 0, 0, 0, 1);
        wavy.push_back_atom_basis_set(0, c_exp, vals[0], type, 0);
        primitive p(1, type, c_exp, vals[0]);


        for (unsigned int l = 0; l < type * 2 + 1; l++)
        {
            int m = static_cast<int>(l) - static_cast<int>(type);
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
            std::cout << "    |  m: " << m << std::endl;

            for (int i = 0; i < grid[0].size(); i++)
            {
                //grid[3][i] = wavy.compute_dens(grid[0][i], grid[1][i], grid[2][i]);
                grid[3][i] = calc_density_ML(grid[0][i], grid[1][i], grid[2][i], coefs, wavy.get_atoms());
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
				double k_pt_local[4] = { kpts[i][0] * 2 * constants::PI , kpts[i][1] * 2 * constants::PI , kpts[i][2] * 2 * constants::PI , 0.0 };
				k_pt_local[3] = sqrt(k_pt_local[0] * k_pt_local[0] + k_pt_local[1] * k_pt_local[1] + k_pt_local[2] * k_pt_local[2]);
				for (int d = 0; d < 3; d++) k_pt_local[d] /= k_pt_local[3];

                sf_A[0][i] = sfac_bessel(p, k_pt_local, coefs);
				//sf_A[0][i] = sfac_bessel(p, k_pt_local, ri_coefs);
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
                if (abs(diff) > 2E-5)
            {
                all_correct = false;
                    correct = false;
            }
            }
            if (!correct)
            {
                std::cout << "Error at m: " << m << "   Max diff: " << std::setprecision(6) << max_diff << std::endl;
                correct = true;
            }
            else
            {
                std::cout << "m: " << m << " passed!" << std::endl;//"   Max diff: " << std::setprecision(6) << max_diff << std::endl;
            }
        }
        if (!all_correct)
            break;
        std::cout << "l = " << type << " passed!\n" << std::endl;
    }
    if (!all_correct)
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
        std::cout << "Error in the calculations!" << std::endl;
        exit(1);
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
    wavy.push_back_atom_basis_set(0, 1.0, 1.0, lambda, 0);
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
    CubeMO.set_path("Oribital-lam" + std::to_string(lambda) + "-m-" + std::to_string(m) + ".cube");
    CubeMO.write_file();
};

// Function to calculate a cubefile containing the electron density of a given wavefunction using the RI-Fit and the regular basis set
// Also calculate the difference between the two densities
void gen_CUBE_for_RI(WFN wavy, const std::string aux_basis, const options *opt)
{
    using namespace std;
    vec ri_coefs = density_fit(wavy, aux_basis, (*opt).mem, 'C');
    ri_coefs = { -0.08681301814153446,0.6893612840685148,1.851248983035524,0.7200605057888373,2.533894657455927,0.06769927976810412,0.16486727327428144,0.10214660297188578,-0.43587807175273924,0.6716291007095919,0.1535971503230409,0.20214410134459035,1.833057564987861e-05,-0.03171100856597612,0.0007015833908988174,-0.0007310482318135529,-0.0007614062703634593,0.005968004468057129,0.0037140369644389705,0.0037517475933302566,0.003619491878756544,-0.008629677384914885,-0.008639374812152208,0.01275789017798526,0.013069500960470134,0.012272174192901655,-0.012041396738739657,-0.014737982245753286,-0.010456261157422214,-0.026127834999407715,0.02366277836449113,0.0059551850182221165,-0.02015832878377889,-0.029930488471320703,-0.0052848515192608775,-0.031017181709458143,-0.005842558941586559,-0.0190190089404415,-0.00010210903489001871,0.01283490562235719,-0.0039905933051441606,0.004468321986025498,-0.0022428396874364822,0.00015253304343987873,-0.0007916457410083985,-0.0007366171809580118,-0.001198288841071443,-0.000130916275132039,0.00015102670476952503,-0.019217137291393495,0.0005331701092988829,0.0027301474783096415,-0.007835301670683595,0.00700305753482342,-0.02965984945947426,-0.001347954356702975,0.0013017283357612333,-0.012094916901129163,0.010115778910707755,-0.042700331598951806,-0.009124941225256936,0.015530859072978355,-0.013575110353053515,0.007537274596419995,-0.0025812545151315575,0.006741793915856968,0.0010414160788195228,-0.002825888932429283,0.00531220816215558,0.003369849953836014,0.0014275825084176554,-0.0018606440890580737,0.000607057628134346,0.000891316375502051,0.006274992256412367,4.3160139482851816e-05,-0.006899670545902283,-0.008991869659576348,0.0012530462176821636,-0.0024768499135500513,-0.005493239680105929,0.01509270547893319,-0.005861503749366906,-0.015510795696269536,-0.021981877111553822,0.009438884114965832,-0.006321113586864953,-0.01223483443392395,-0.0014637222674403682,0.0063536149311118045,-0.0032187962806403004,0.0019211369740930416,0.0027254473139287187,-0.0013536396857086286,-0.005458656450159162,-0.005588315818205726,-0.001473184454485633,-0.0005544395398563846,0.09487229864209928,0.15381051632863044,0.012830923086788766,-0.0014376208558240606,0.006097596464243993,0.024271830357534906,0.0017320775415754985,-0.001201077004920023,-0.0049096220092887865,-0.0002736975710647482,0.002376528661194932,0.0053307744841738,-0.0006080509121370887,-0.0004522693154487152,-0.0002334650701154084,-0.0006192927742801253,-0.001164986721310418,0.000507642752614726,-0.0003168266570422083,-0.08507280328032563,0.6839962914235038,1.8654586548033918,0.6993420926649301,2.572072200600347,0.006886635173806106,0.25958667410039205,-0.03861600616155372,-0.23694541875220396,0.42706148862166277,0.4080478884993377,-0.004942103756574089,-0.016749585144505103,-0.0634844130284931,-0.0005340437198496752,-4.876030614452072e-05,6.819393852316126e-05,0.0031751569511941904,-0.0001969086563999575,-0.00021371806360071928,-0.0059841400457550774,-0.00040606423094525346,0.0006317064485312727,0.009663685287401238,6.915124079394795e-05,-0.00046460340077295146,-0.005890375582218066,-0.0028984386593294675,-0.0011474257227171792,0.011640458186605033,0.004948500930925754,0.0016597601024863089,-0.020890857330992756,-0.0022838927036026957,-0.003168997879223405,-0.017998941633342747,-0.010479417735768905,-0.0066788923656819574,0.023716808109736812,-0.004497557571525183,-0.0027514526715679837,0.01381355471182715,-0.005412986914857128,-0.0022240461597631384,0.0002932073628091074,-0.0004542351372461123,0.0003745876408836843,-3.957654019871617e-05,-0.00032272287762341194,0.001315613675601378,-0.013014619559867489,0.007962561842640622,-0.0009936632811907173,-0.004089257979649182,0.0051981097193574246,-0.012096213322501416,0.011080931184629612,0.00033897843501305805,-0.009885548320681954,0.01922298178529268,0.016971520270784437,-0.0003505238713526819,0.013129192137913618,-0.008990935278376425,-0.001689290975127314,-0.01247119778819468,0.0008545058235718797,-0.0022203313664695616,0.006422010104477395,-0.007134686635851543,-0.008594237414273054,2.7419701699097637e-05,-0.0050996658058192765,0.004664000155280782,-0.006109915880358499,-0.0017669668811111908,0.001717041127668536,0.0009127353581353208,-0.002469669723916244,-0.004469268355650732,0.002997295033119306,-0.022354945656802575,-0.011841001755244177,0.00783530161164578,0.006922025274806692,-0.005704435027895609,-0.014408209024117621,0.01104645088952426,-0.0013691430311808007,0.0015352228087613436,0.001027378723012383,-0.0008601832689122663,0.000606376519964327,0.0007483617296510767,-0.0008710408884889968,-0.0004306867719941204,-0.0004815384826829635,-0.23576269958254878,2.202311861110941,2.222868495566912,2.771367221727892,2.799102967824014,-0.8850711392946219,0.8287188006889831,-0.9480777589725496,1.0823017640867234,0.9825704463863301,0.08487317450882244,0.1745199499782286,-0.01813806549357264,0.04995127885016351,0.0003979893666664279,-0.0030964318388569984,-0.0021807424709330233,0.00505964365498927,-0.03019702126603878,-0.01993278178604891,-0.001493867727953757,0.003866701653914171,0.0014814727083645827,0.007434547816450484,-0.03887353379574684,-0.02428298506754608,-0.019126347831524372,0.11632162527042784,0.07611204533275799,-0.008372660125075267,0.04978560006218384,0.0359899789108796,0.002505843450649891,-0.001991276573246951,-0.0044741639144046256,0.004137802150847714,-0.03161164127996354,-0.02018471860987474,0.0056274326879405005,-0.01313728279052992,-0.008762608548449662,-0.0006834826471512694,-0.00021429758923554185,-0.00016219342978880593,-0.0009508387480425727,0.0007251144468927528,0.001974658805141642,-0.0011950480183036745,-0.003641949139270999,0.003438162184316552,-0.004344258830499491,-0.005773828731011663,0.0040287410866097235,0.01283114377683168,0.010035761510890879,-0.020709436369163,-0.010279703894438696,0.010317546553159003,0.03607296016290538,0.018533554232721662,-0.02642166145708479,-0.02666763636207603,0.02097096640115567,0.06689077033400395,0.0004964653231775063,0.012890032713600365,-0.007178145229527855,0.0023622027703213513,0.004881265399811345,-0.003416240772509798,0.010831866862312009,-0.0003133764987691933,-0.001999139743346828,-0.008376604865276699,0.003991033458084582,0.001781335859586958,-0.001785409906398532,0.0015179027339469234,0.0006225934638857192,0.004592179342260691,-0.0015736506934121341,0.006080326970283076,0.003316267858606847,-0.003763599476248201,0.0033681188004270313,0.0012569724037651306,0.008869581302461947,-0.0031143247884653858,0.0011727838829784358,-0.002375582218047913,-0.0009859017173609712,-7.755504889626345e-05,-0.0006539673078125416,-0.0002662043149570598,-0.000955919332770109,0.0015259105747063664,0.0016286751154723127,-0.2358697324295431,2.202948461581667,2.222796879809394,2.7727907853688207,2.7976685007236473,-0.8792429546659035,0.8209901081072289,-0.9396300770058097,1.070040759399461,1.0027266678075029,0.04698137575373166,0.1772747327211374,-0.03418498042454093,0.0552600045279181,0.002888120943299002,0.0019461047024776704,0.001546469869114806,0.02680872185357136,0.01925525757011164,0.0155419008502233,-0.002571740129737749,-0.0027700434034530425,-0.0024441815158268727,0.033638185836739394,0.025196371030673105,0.020614786554832537,-0.10339476857495064,-0.07448097283144794,-0.06027550485091356,-0.04362528359669436,-0.029215857486079442,-0.022837425976730105,-0.008518568755786174,-0.007315510628689127,-0.006854204349627159,0.02216713030298856,0.01515139634119377,0.011970069114849548,0.013728604498036813,0.006968393062584815,0.005310330176786095,-3.846212082776595e-05,0.000659293920476244,0.00047105654522560454,0.003364947416377935,-0.0020088316787933444,0.000617228610922546,0.00235676166830499,0.0004322924565514646,-0.012431880033646255,0.00493645566973835,-0.0012345919167692536,-0.008961019929316657,-0.002118576048520116,-0.038048768052159974,0.0031383334888067426,0.00119599575118248,-0.02873312327270275,-0.009049864525398462,-0.06697699758558023,0.020254882386430795,-0.004125018909333606,-0.04909002948728123,-0.013112789420948736,-0.007091662351841624,0.011281624476934649,-0.005186337280421391,-0.005246646832436704,-0.00032927327452974436,0.008267317488779425,0.003146603331191929,-0.003092961762781206,0.006304437545834044,0.0037798134610006607,0.00418518066540268,0.003575849874270269,-0.0007783383717300498,-0.0020039735598586946,0.00045158398439702947,0.002405416710379833,-0.0017274353974220929,0.006759767970527159,0.005446344153857734,-0.0012351600246597544,-0.0036999809014398805,3.195313251329335e-05,0.003374368167295069,-0.0026661415661430594,0.0011575843889027422,0.0020789374756033604,-9.330405191712378e-05,-0.0004946492881508596,-0.00023040728755126323,-0.0007702277868727982,0.0011636778478614346,6.821758072489958e-05,-0.00215301266333831,-0.17331665336915641,1.3988877234582546,2.02516348605469,1.6440206313894816,2.8518008588611146,-0.49067323261541096,0.520757445189782,-0.38997208785557447,0.039349065646602545,0.8733671438086376,0.2438867288314144,0.14728240118906388,-0.0081349796689332,-0.007779133190614983,0.0011077082976357832,0.0013224128711906392,0.00010534684361874152,-0.0002391311583964918,-0.0007475045316903224,-0.0003592512307975999,0.0050337260151094385,0.006103326362341144,0.0006003069744326626,-0.0011635627949530642,-0.002533181336054609,-0.00080176916414727,-0.004737786584166688,-0.003920163460449928,0.0002665645379408945,-0.00155773169940268,-0.002320547234392584,0.003319210070425585,-0.023198542337826582,-0.02133654018107026,-0.006971281604361585,0.014451618126808733,0.009155726831047324,0.0016628253712651381,-0.0037215591408672893,0.001899741102370645,0.00018547763321123066,-0.00025547142697132625,0.0001804062680123568,6.637451304330613e-05,-0.0004783660732812326,-9.506285510518686e-05,0.00016261598174775139,-0.0001081737841194563,0.00016136517202786386,0.0011956954446577734,0.0002689332108363106,-0.0003972070790463037,0.00026453149906050183,-0.0004150546935525095,0.006991257814496737,0.0026355075824655216,-0.002231603986226554,0.0014896571578302694,-0.002223954123473409,0.006321746656429517,0.0022019483244732443,-0.0002491649558290266,0.001085720447116653,-0.0037240752071263923,0.01930604852263946,0.005472284413816955,-0.008367037957597159,0.004735826010362845,-0.0033641127229612645,0.00038035494924048424,-0.0018450527129051592,8.590202680151614e-05,0.0011782091261263621,-0.002776787347402735,-0.0076062521881949405,0.00020480758388964477,0.008958346804262723,0.007503434851942676,-0.003062317219167238,0.0021607961057562145,0.0056384088360996515,-0.01551708013206894,-0.0031092753589507856,0.01689216938796152,0.014583399424938403,-0.0009494815772025582,0.0036119420272047406,0.015895164098664417,-0.0021150047197340703,0.0032119705753079473,-0.003297680483558907,-0.0010956300935772688,0.002601374960403305,-0.0022693537868436198,-0.0022371696892624347,-0.0031873872128468663,-0.0042632700289505325,-0.08556272338803589,0.6844817635255568,1.8561067691665558,0.7059201091418343,2.551209566542173,0.03426424112917056,0.22202066925733663,0.013232228261252596,-0.3061468338904918,0.5010554784950112,0.2939270850818008,0.16930608224622712,0.04024508884780622,-0.015754253550975362,0.00012707300360273113,-0.0004418886460862467,0.0003828751834631668,-0.0007046316918607121,0.0023460630819548546,-0.002019342260642881,0.0015150721617499406,-0.00503142196519773,0.004439766179048681,-0.0022882824515278823,0.007396779000749745,-0.0065949749526948535,0.002226328163671824,-0.005328799877720361,0.005419587250048863,-0.000657108301179072,-0.0012009438958287273,-0.000549049786905344,2.2573127921779505e-05,0.0008665513671100548,0.0004091797254828577,0.003010332555768091,-0.016549867574585164,0.009801660638468393,-0.001598036260991517,0.010370278352212678,-0.003911904165560951,0.0005685695312529943,-0.00014380139187023343,0.00016312590441365285,-2.341761051828933e-05,-0.000958478685671724,-0.00016235987528953152,0.0002401005974446749,-0.0007763570156713491,0.002803667973849678,0.003812459688064305,0.00015936628502453515,-0.00043334850446181543,0.002024795697811651,0.0048557588164199585,0.003488320967421477,-0.00025818011768914993,-7.002462013167456e-05,0.0006602525559310197,0.00521953752883282,0.004577718373457078,-0.0006523158095866523,0.00016336877647889065,0.002111336779725526,-0.00023717116813769,0.000831590178043768,-0.00014876378072572552,-0.00016266235140618927,-0.0010808630145911053,-0.00020394683376852562,-0.0006918783955250131,0.0003976686365013385,9.400656854317948e-05,-0.0011007740197459454,-0.0087629053318891,0.001491286363141169,0.00778009024955056,0.009960745639449367,-0.0007747782201967543,0.0014606643063885636,0.00596832750277871,-0.018421865052583843,0.0035611791332697865,0.016252237989364923,0.02032449662004464,-0.0022648212402482633,0.004084239971331266,0.012699080563633823,-0.00039120740624274536,0.005409022354173752,-0.002123336292177247,0.0027384103758144244,0.0008636750179564531,-1.9610294704658308e-05,-0.004405585474417817,-0.0040574876457534145,-0.0009316803593656504,-0.00045465284170991373,0.09743914963025628,0.16183911620839286,0.015107860336351848,-0.024486565318824673,0.009977638354850653,0.005977828544432221,0.004519586094323767,-0.000787107771108477,-0.0017631965999769606,-0.003628704328792993,0.0007810810762371919,-0.002768423792591578,-0.002348012720731554,0.0037935111578550635,0.0007142608625752695,-0.00044939788947488435,0.0006583653980633328,0.0009140265572475312,-0.0013491440885929507,-0.0005309143266242666,0.08833481306900115,0.15077619909966927,0.014201053989535752,0.018711255057062625,0.01884876327840628,0.004390927369336767,-0.003447952046203477,-0.0017021626888461643,-0.001512196535515893,0.005397855505380924,0.0012438079355942034,-0.002982499122937034,0.0014307135131040602,7.872048332065014e-05,-0.001266716601332511,-0.0007556341369046814,0.0006902106315738091,-0.000582082532797278,-0.0006901399688169125,-0.0005590818818568679,0.0942997250133335,0.15867969118864278,0.01137592445382254,0.00039322893556583636,-0.006746260389568438,-0.02692436311432345,-0.00010461334797769311,0.001384546086023641,0.004773235236543802,-5.695188868693325e-05,0.002892736572464779,0.005750395616486729,-0.00024051803066603642,-0.0005390391975041489,6.044820083181168e-05,-0.000824334079363417,-0.00139978906632613,0.00012000982681702477,-0.00011628883948937296,-0.001360949187811544,0.07308024667112523,0.13087490515575834,0.010976872419609944,0.002677993169803291,-0.009499006265877791,-0.02853228449321657,-0.0012554692891001639,-0.0005462367932003389,-0.00043229506989730417,-0.0002727303159876394,0.004222648012265378,0.007203959162991366,-0.0014454484085305824,-0.0005597961170759408,-0.00017385064906058644,0.00016308883339285593,-0.00010581965477153756,0.0003480666927293428,0.0004385548839448275,-0.001641179550740793,0.07294006556330769,0.13009397116607688,0.0077728777656136125,0.005002247855770243,-0.023138951529334394,0.019682423448294954,-0.0010242197775228194,-0.0007339217505112518,-0.0008254678549794377,-0.0019122087856253973,-0.007267292514753519,0.0009786035688137873,0.0017953669184210528,-0.004088152180591257,0.00032017658230294174,-0.00012795859500914632,-0.00033326143826153135,-0.0006361645821989993,-0.00036344474248125615,-0.0015102748357514418,0.07801999464840029,0.13779930602145593,0.011664116405709102,-0.028672607642381262,0.011648973543755967,0.004271876194759714,0.0003524172718833207,-0.0007044006247829684,-0.0011244756601173679,-0.005273859564007602,0.0008255949099622277,-0.004199986900312371,-0.0018651625936797814,0.005266283620739999,0.0005495163420748253,-0.0002667822066026791,0.0001605751970627145,0.0006145494185696306,-4.362552940910637e-05};
    //for (int i = 0; i < ri_coefs.size(); i++) {
    //    std::cout << ri_coefs[i] << "    " << ri_coefs2[i] << "   " << ri_coefs[i] - ri_coefs2[i] << std::endl;
    //}

    WFN wavy_aux(0);
    wavy_aux.set_atoms(wavy.get_atoms());
    wavy_aux.set_ncen(wavy.get_ncen());
    wavy_aux.delete_basis_set();
    load_basis_into_WFN(wavy_aux, BasisSetLibrary().get_basis_set(aux_basis));


    vec3 grid;
    ivec pointy = fuckery(wavy, grid, 1.0);
#pragma omp parallel
    {
        vec2 d_temp(16);
        for (int i = 0; i < 16; i++)
        {
            d_temp[i].resize(wavy.get_ncen(), 0.0);
        }
        vec phi_temp(wavy.get_nmo(), 0.0);

        for (int a = 0; a < wavy.get_ncen(); a++) {
#pragma omp for
            for (int i = 0; i < pointy[a]; i++)
            {
                grid[a][4][i] = wavy.compute_dens(
                    grid[a][0][i],
                    grid[a][1][i],
                    grid[a][2][i],
                    d_temp,
                    phi_temp) * grid[a][5][i];
                //grid[a][6][i] = calc_density_ML(
                //    grid[a][0][i],
                //    grid[a][1][i],
                //    grid[a][2][i],
                //    ri_coefs,
                //    wavy_aux.get_atoms(),
                //    a) * grid[a][3][i];
                grid[a][6][i] = calc_density_ML(
                    grid[a][0][i],
                    grid[a][1][i],
                    grid[a][2][i],
                    ri_coefs,
                    wavy_aux.get_atoms()) * grid[a][5][i];
            }
        }
        for (int i = 0; i < 16; i++)
            shrink_vector<double>(d_temp[i]);
        shrink_vector<vec>(d_temp);
        shrink_vector<double>(phi_temp);
    }
    vec elecs_DFT(wavy.get_ncen(), 0.0);
    vec elecs_RI(wavy.get_ncen(), 0.0);
    for (int a = 0; a < wavy.get_ncen(); a++)
    {
        elecs_DFT[a] = vec_sum(grid[a][4]);
        elecs_RI[a] = vec_sum(grid[a][6]);
    }
    vec analytic_RI = calc_atomic_density(wavy_aux.get_atoms(), ri_coefs);

    std::cout << "Table of Charges in electrons\n"
        << "       Atom  DFT     RI      Diff" << std::endl;
    for (int i = 0; i < wavy.get_ncen(); i++)
    {
        std::cout << std::setw(10) << wavy.get_atom_label(i)
            << std::fixed << std::setw(10) << std::setprecision(4) << elecs_DFT[i]
            << std::fixed << std::setw(10) << std::setprecision(4) << elecs_RI[i]
                << std::fixed << std::setw(10) << std::setprecision(4) << analytic_RI[i]
                << std::fixed << std::setw(10) << std::setprecision(4) << elecs_DFT[i] - elecs_RI[i]
                    << std::endl;
    }

    std::cout << "SUM OF CHARGES: " << std::endl;
    std::cout << "DFT: " << std::fixed << std::setprecision(4) << vec_sum(elecs_DFT) << std::endl;
    std::cout << "RI:  " << std::fixed << std::setprecision(4) << vec_sum(elecs_RI) << std::endl;
    std::cout << "Analytic:  " << std::fixed << std::setprecision(4) << vec_sum(analytic_RI) << std::endl;
    







	//const double radius = 3.5;
	//const double resolution = 0.02;

 //   std::cout << "-------------------------------------DENSITY USING ORCA GBW-------------------------------------" << std::endl;
 //   double MinMax[6]{0, 0, 0, 0, 0, 0};
 //   int steps[3]{0, 0, 0};
 //   readxyzMinMax_fromWFN(wavy, MinMax, steps, radius, resolution, true);
 //   cube cube_normal(steps[0], steps[1], steps[2], wavy.get_ncen(), true);
 //   cube cube_RI_FIT(steps[0], steps[1], steps[2], wavy.get_ncen(), true);


 //   std::filesystem::path normal_cube_path = std::filesystem::path(wavy.get_path().stem().string() + "_normal_rho.cube");
 //   cube_normal.set_path(normal_cube_path);
 //   // Check if the cube file already exists, if so read it
 //   for (int i = 0; i < 3; i++)
 //   {
 //       cube_normal.set_origin(i, MinMax[i]);
	//	cube_RI_FIT.set_origin(i, MinMax[i]);
 //       cube_normal.set_vector(i, i, (MinMax[i + 3] - MinMax[i]) / steps[i]);
	//	cube_RI_FIT.set_vector(i, i, (MinMax[i + 3] - MinMax[i]) / steps[i]);
 //   }
 //   if (std::filesystem::exists(normal_cube_path))
 //   {
 //       //cube_normal.read_file(true, true);
 //   }
 //   else
 //   {
 //       cube_normal.give_parent_wfn(wavy);
 //       cube_normal.set_comment1("Calculated density using NoSpherA2");
 //       cube_normal.set_comment2("from " + wavy.get_path().string());
 //       Calc_Rho(cube_normal, wavy, radius, std::cout, false);
 //       cube_normal.write_file(true);
 //   }
 //   cube_normal.calc_dv();
 //   std::cout << "Number of electrons: " << std::fixed << std::setprecision(4) << cube_normal.sum() << std::endl;

 //   std::cout << "-------------------------------------RI Fit-------------------------------------" << std::endl;
 //   vec ri_coefs = density_fit(wavy, aux_basis, (*opt).mem, 'C');
 //   std::cout << "Calculating RI-FIT cube" << std::endl;

 //   WFN wavy_aux(0);
 //   wavy_aux.set_atoms(wavy.get_atoms());
 //   wavy_aux.set_ncen(wavy.get_ncen());
 //   wavy_aux.delete_basis_set();
 //   load_basis_into_WFN(wavy_aux, BasisSetLibrary().get_basis_set(aux_basis));

 //   
 //   cube_RI_FIT.give_parent_wfn(wavy_aux);
 //   //calc_cube_ML(ri_coefs, wavy_aux, -1, cube_RI_FIT);
 //   //cube_RI_FIT.set_path(std::filesystem::path(wavy.get_path().stem().string() + "_RI_FIT_rho.cube"));
 //   //cube_RI_FIT.write_file(true);

 //   std::cout << "-------------------------------------DIFFERENCE CUBE-------------------------------------" << std::endl;
 //   //cube cube_diff = cube_normal - cube_RI_FIT;
 //   //cube_diff.give_parent_wfn(wavy);
 //   //cube_diff.set_comment1("Difference between RI-FIT and normal density");
 //   //cube_diff.set_comment2("from " + wavy.get_path().string());
 //   //cube_diff.set_path(std::filesystem::path(wavy.get_path().stem().string() + "_diff.cube"));
 //   //cube_diff.write_file(true);
 //   //cube_diff.calc_dv();
 //   //std::cout << "Number of electrons in difference cube: " << std::fixed << std::setprecision(5) << cube_diff.sum() << std::endl;
 //   //std::cout << "RRS of difference cube: " << std::fixed << std::setprecision(5) << cube_normal.rrs(cube_RI_FIT) << std::endl;
 //   

 //   std::cout << "------------------------------------RI Fit analytical integral------------------------------------" << std::endl;
 //   vec atom_elecs = calc_atomic_density(wavy_aux.get_atoms(), ri_coefs);

 //   vec atom_elecs2(wavy_aux.get_ncen());
 //   for (int i = 0; i < wavy_aux.get_ncen(); i++) {
 //       cube_RI_FIT.set_zero();
 //       calc_cube_ML(ri_coefs, wavy_aux, i, cube_RI_FIT);
 //       atom_elecs2[i] = cube_RI_FIT.sum();
 //       cube_RI_FIT.write_file(true);
 //   }


 //   std::cout << "Table of Charges in electrons\n"
 //       << "       Atom  RI-Ana    NElec   NElec_calc" << std::endl;

 //   for (int i = 0; i < wavy_aux.get_ncen(); i++)
 //   {
 //       std::cout << std::setw(10) << wavy_aux.get_atom_label(i)
 //           << std::fixed << std::setw(10) << std::setprecision(3) << wavy_aux.get_atom_charge(i) - atom_elecs[i];
 //       std::cout << " " << std::setw(4) << wavy_aux.get_atom_charge(i) << " " << std::fixed << std::setw(10) << std::setprecision(10) << atom_elecs[i];
 //       std::cout << " " << std::setw(10) << std::setprecision(10) << atom_elecs2[i];
 //       std::cout << std::endl;
 //   }

 //   std::cout << "Done!" << std::endl;
}

void test_reading_SALTED_binary_file() {
    std::filesystem::path path("Model/model.salted");
    SALTED_BINARY_FILE file = SALTED_BINARY_FILE(path, true);
    Config config;
    file.populate_config(config);
    std::unordered_map<int, std::vector<int64_t>> fps = file.read_fps();
	std::unordered_map<std::string, vec> averages = file.read_averages();
    std::unordered_map<int, vec> wigners = file.read_wigners();
    vec weights = file.read_weights();
	std::unordered_map<std::string, dMatrix2> feats = file.read_features();
    std::unordered_map<std::string, dMatrix2> proj = file.read_projectors();

	// TEST if both configs are the same
	std::cout << "Comparing configs" << std::endl;
	std::cout << "Average:" << config.average << std::endl;
	std::cout << "Field:" << config.field << std::endl;
	std::cout << "Sparsify:" << config.sparsify << std::endl;
	std::cout << "Ncut:" << config.ncut << std::endl;
	std::cout << "Ntrain:" << config.Ntrain << std::endl;
	std::cout << "Menv:" << config.Menv << std::endl;
	std::cout << "trainfrac:" << config.trainfrac << std::endl;
    std::cout << "Rcut1:" << config.rcut1 << std::endl;
	std::cout << "Rcut2:" << config.rcut2 << std::endl;
	std::cout << "nang1:" << config.nang1 << std::endl;
	std::cout << "nang2:" << config.nang2 << std::endl;
	std::cout << "sig1:" << config.sig1 << std::endl;
	std::cout << "sig2:" << config.sig2 << std::endl;
	std::cout << "zeta:" << config.zeta << std::endl;
	std::cout << "neighspe size:" << config.neighspe1.size() << std::endl;
    for (int i = 0; i < config.neighspe1.size(); i++)
    {
		std::cout << "neighspe1[" << i << "]:" << config.neighspe1[i] << std::endl;
	}
	std::cout << "neighspe2 size:" << config.neighspe2.size() << std::endl;
	for (int i = 0; i < config.neighspe2.size(); i++)
	{
		std::cout << "neighspe2[" << i << "]:"  << config.neighspe2[i] << std::endl;
	}
	std::cout << "dfBasis:" << config.dfbasis << std::endl;

	std::cout << "Comparing wigners" << std::endl;
	for (int i = 0; i < wigners.size(); i++)
	{
		for (int j = 0; j < wigners[i].size(); j+=10)
		{
			std::cout << wigners[i][j] << " ";
		}
		std::cout << std::endl;
	}

	std::cout << "Comparing FPS" << std::endl;
	for (int i = 0; i < fps.size(); i++)
	{
		for (int j = 0; j < fps[i].size(); j+=10)
		{
			std::cout << fps[i][j] << " ";
		}
        std::cout << std::endl;
	}

	std::cout << "All tests passed!" << std::endl;
}


// Convert a row-major matrix to column-major format
std::vector<double> rowToColMajor(const std::vector<double> &rowMajorMatrix, int rows, int cols)
{
    std::vector<double> colMajorMatrix(rows * cols);

    for (int r = 0; r < rows; ++r)
    {
        for (int c = 0; c < cols; ++c)
        {
            colMajorMatrix[c * rows + r] = rowMajorMatrix[r * cols + c];
        }
    }
    return colMajorMatrix;
}

void test_NNLS()
{
    // Define the matrix A in column-major order (Fortran-style storage)
    std::vector<double> A_in = {
        54.88135039, 71.51893664, 60.27633761, 54.4883183, 42.36547993,
        64.58941131, 43.75872113, 89.17730008, 96.36627605, 38.34415188,
        79.17250381, 52.88949198, 56.80445611, 92.55966383, 7.10360582,
        8.71292997, 2.02183974, 83.26198455, 77.81567509, 87.00121482,
        97.86183422, 79.91585642, 46.14793623, 78.05291763, 11.82744259,
        63.99210213, 14.33532874, 94.4668917, 52.18483218, 41.466194,
        26.45556121, 77.42336894, 45.61503322, 56.84339489, 1.87898004,
        61.76354971, 61.20957227, 61.69339969, 94.37480785, 68.18202991,
        35.95079006, 43.70319538, 69.76311959, 6.02254716, 66.67667154,
        67.06378696, 21.03825611, 12.89262977, 31.54283509, 36.37107709};
     //Dimensions
    const int m = 10; // Number of rows
    const int n = 5;  // Number of columns
     //Right-hand side vector B
    std::vector<double> B_in = {57.01967704, 43.86015135, 98.83738381, 10.20448107, 20.88767561, 16.13095179, 65.31083255, 25.32916025, 46.63107729, 24.4425592};

    dMatrix1 B(B_in.size());
    std::copy(B_in.begin(), B_in.end(), B.data());
    dMatrix2 A(m, n);
	std::copy(A_in.begin(), A_in.end(), A.data());
    NNLSResult res = nnls(A, B);

	//double sum = std::accumulate(res.x.begin(), res.x.end(), 0.0);
	//std::cout << "Sum of coefficients: " << sum << std::endl;

    std::cout << "NNLS solution: " << std::endl;
    std::cout << "Mode: " << res.status << std::endl;
    std::cout << "X: ";
    for (int i = 0; i < res.x.size(); i++)
    {
        std::cout << res.x[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "Residual: " << res.rnorm << std::endl;

    // Check vs correct solution
    std::vector<double> correct = {0.00370583, 0.58600469, 0.14755389, 0.04786599, 0.0};
    for (int i = 0; i < res.x.size(); i++)
    {
        if (std::abs(res.x[i] - correct[i]) > 1E-8)
        {
            std::cout << "NNLS solution is incorrect!" << std::endl;
            exit(1);
        }
    }
    std::cout << "NNLS solution is correct!" << std::endl;
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