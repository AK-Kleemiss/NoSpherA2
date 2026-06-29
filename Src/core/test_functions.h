#pragma once
#include "spherical_density.h"
#include "scattering_factors.h"
#include "AtomGrid.h"
#include "convenience.h"
#include "npy.h"
#include "properties.h"
#include "basis_set.h"
#include "SALTED_utilities.h"
#include "sphere_lebedev_rule.h"
#include "integrator.h"
#include "nos_math.h"
#include "featomic.hpp"
#include "metatensor.h"
#include "GridManager.h"
#ifdef __APPLE__
#include "TargetConditionals.h"
#endif

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
    d.resize(1);
    for (int i = 0; i < 1; i++)
        d[i].resize(16, 0.0);
    angular_off -= start;
    const int size = start + num_angular;
    int p = 0;
    double dens = 0.0;
    for (int iang = start; iang < size; iang++)
    {
        p = angular_off + iang;
        const d3 Pos = { angular_x[p] * r + wavy.get_atom_coordinate(0, 0), angular_y[p] * r + wavy.get_atom_coordinate(0, 1), angular_z[p] * r + wavy.get_atom_coordinate(0, 2) };
        dens += wavy.compute_dens(Pos, d, _phi) * constants::FOUR_PI * angular_w[p];
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
        vec2 d(wavy.get_ncen());
        vec _phi(wavy.get_nmo(), 0.0);
        for (int j = 0; j < wavy.get_ncen(); j++)
            d[j].resize(16, 0.0);
#pragma omp for reduction(+ : dens)
        for (int iang = start; iang < size; iang++)
        {
            p = angular_off + iang;
            const d3 Pos = { angular_x[p] * r + wavy.get_atom_coordinate(i, 0), /// moving the lebedev points to the atom i
            angular_y[p] * r + wavy.get_atom_coordinate(i, 1),
            angular_z[p] * r + wavy.get_atom_coordinate(i, 2) };
            const d3 d_ = { angular_x[p] * r, angular_y[p] * r, angular_z[p] * r }; /// as the function get_radial_density needs a distance, the distance is calculated. The atom A is in the origin, the wavy.get_atom_coordinate(i,0) is 0
            const double dist = array_length(d_);
            const double rho_a = A.get_radial_density(dist); /// the radial density of the atom i is calculated
            double rho_all = rho_a;                          /// the molecular density based on pro-atoms
            for (int atom = 0; atom < wavy.get_ncen(); atom++)
            { /// a new for loop is started, which calculates the sum of rhos, with respect to the lebedev points
                if (atom == i)
                    continue; /// if the atom is the same as the atom i, the loop is skipped
                const d3 d_atom = { Pos[0] - wavy.get_atom_coordinate(atom, 0), Pos[1] - wavy.get_atom_coordinate(atom, 1), Pos[2] - wavy.get_atom_coordinate(atom, 2) };
                const double dist_atom = array_length(d_atom); /// is like d_, but for another atom
                Thakkar thakkar_atom(wavy.get_atom_charge(atom));
                rho_all += thakkar_atom.get_radial_density(dist_atom);
            }
            const double hirsh_weight = rho_a / rho_all;
            dens += wavy.compute_dens(Pos, d, _phi) * hirsh_weight * constants::FOUR_PI * angular_w[p];
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
        vec2 d(wavy1.get_ncen());
        vec _phi(std::max(wavy1.get_nmo(), wavy2.get_nmo()), 0.0);
        for (int j = 0; j < wavy1.get_ncen(); j++)
            d[j].resize(16, 0.0);
#pragma omp for reduction(+ : dens)
        for (int iang = start; iang < size; iang++)
        {
            p = angular_off + iang;
            const d3 Pos = { angular_x[p] * r + wavy1.get_atom_coordinate(0, 0),
                                                 angular_y[p] * r + wavy1.get_atom_coordinate(0, 1),
                                                 angular_z[p] * r + wavy1.get_atom_coordinate(0, 2) };
            dens += (wavy1.compute_dens(Pos, d, _phi) - wavy2.compute_dens(Pos, d, _phi)) * constants::FOUR_PI * angular_w[p];
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

    properties_options opts;
    opts.radius = r;
    opts.resolution = resol;
    readxyzMinMax_fromWFN(wavy1, opts, true);

    properties_options opts2;
    opts2.radius = r;
    opts2.resolution = resol;
    readxyzMinMax_fromWFN(wavy2, opts2, true); // and getting the MinMax and steps

    double MinMax[6]{ 100, 100, 100, -100, -100, -100 };
    int steps[3]{ 0, 0, 0 };
    for (int i = 0; i < 3; i++)
    {
        if (opts.MinMax[i] < MinMax[i])
            MinMax[i] = opts.MinMax[i];
        if (opts.MinMax[i + 3] > MinMax[i + 3])
            MinMax[i + 3] = opts.MinMax[i + 3];
    }
    for (int i = 0; i < 3; i++)
    {
        if (opts2.MinMax[i] < MinMax[i])
            MinMax[i] = opts2.MinMax[i];
        if (opts2.MinMax[i + 3] > MinMax[i + 3])
            MinMax[i + 3] = opts2.MinMax[i + 3];
        steps[i] = (int)ceil(constants::bohr2ang(MinMax[i + 3] - MinMax[i]) / resol);
    }

    cube dens1({ steps[0], steps[1], steps[2] }, wavy1.get_ncen(), true);
    cube dens2({ steps[0], steps[1], steps[2] }, wavy2.get_ncen(), true);
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

void calc_partition_densities()
{
    using namespace std;
    WFN Hartree_Fock("HF.gbw");
    WFN DFT("DFT.gbw");
    Hartree_Fock.delete_unoccupied_MOs();
    DFT.delete_unoccupied_MOs();

    vec Pos_C = { Hartree_Fock.get_atom_coordinate(0, 0), Hartree_Fock.get_atom_coordinate(0, 1), Hartree_Fock.get_atom_coordinate(0, 2) };
    vec Pos_O1 = { Hartree_Fock.get_atom_coordinate(1, 0), Hartree_Fock.get_atom_coordinate(1, 1), Hartree_Fock.get_atom_coordinate(1, 2) };
    vec Pos_H1 = { Hartree_Fock.get_atom_coordinate(2, 0), Hartree_Fock.get_atom_coordinate(2, 1), Hartree_Fock.get_atom_coordinate(2, 2) };
    vec Pos_O3 = { Hartree_Fock.get_atom_coordinate(3, 0), Hartree_Fock.get_atom_coordinate(3, 1), Hartree_Fock.get_atom_coordinate(3, 2) };
    vec x_coords = {
        Hartree_Fock.get_atom_coordinate(0, 0),
        Hartree_Fock.get_atom_coordinate(1, 0),
        Hartree_Fock.get_atom_coordinate(2, 0),
        Hartree_Fock.get_atom_coordinate(3, 0) };
    vec y_coords = {
        Hartree_Fock.get_atom_coordinate(0, 1),
        Hartree_Fock.get_atom_coordinate(1, 1),
        Hartree_Fock.get_atom_coordinate(2, 1),
        Hartree_Fock.get_atom_coordinate(3, 1) };
    vec z_coords = {
        Hartree_Fock.get_atom_coordinate(0, 2),
        Hartree_Fock.get_atom_coordinate(1, 2),
        Hartree_Fock.get_atom_coordinate(2, 2),
        Hartree_Fock.get_atom_coordinate(3, 2) };
    ivec charges{ 6, 8, 1, 8 };
    const int min = -250, max = 700;
    const int size = -min + max + 1;
    vec C_dens(size, 0.0), H_dens(size, 0.0), total_dens(size, 0.0);
    vec HF_densities(size, 0.0);
    vec DFT_densities(size, 0.0);
    vec B_weights_C(size, 0.0), B_weights_H(size, 0.0);
    vec TFVC_weights_C(size, 0.0), TFVC_weights_H(size, 0.0);
    vec TFVC_weights_C_DFT(size, 0.0), TFVC_weights_H_DFT(size, 0.0);
    Thakkar C(6);
    Thakkar H(1);
    Thakkar O(8);
    const double dx = (Hartree_Fock.get_atom_coordinate(2, 0) - Hartree_Fock.get_atom_coordinate(0, 0)) / (min + max);
    const double dy = (Hartree_Fock.get_atom_coordinate(2, 1) - Hartree_Fock.get_atom_coordinate(0, 1)) / (min + max);
    const double dz = (Hartree_Fock.get_atom_coordinate(2, 2) - Hartree_Fock.get_atom_coordinate(0, 2)) / (min + max);

    std::cout << "Positions:\nH:" << Hartree_Fock.get_atom_coordinate(2, 0) << " " << Hartree_Fock.get_atom_coordinate(2, 1) << " " << Hartree_Fock.get_atom_coordinate(2, 2) << "\n";
    std::cout << "C:" << Hartree_Fock.get_atom_coordinate(0, 0) << " " << Hartree_Fock.get_atom_coordinate(0, 1) << " " << Hartree_Fock.get_atom_coordinate(0, 2) << "\n";
    std::cout << "dr:" << dx << " " << dy << " " << dz << "\n";



    ProgressBar *pb = new ProgressBar(size, 100, "=", "", "Calculating densities");

#pragma omp parallel 
    {
        vec2 d(DFT.get_ncen());
        for (int i = 0; i < DFT.get_ncen(); i++)
            d[i].resize(16);
        vec phi(DFT.get_nmo(), 0.0);
        vec pa_b(5);
        vec pa_tv(5);
        vec chi_HF = make_chi(Hartree_Fock, 200);
        vec chi_DFT = make_chi(DFT, 200);
#pragma omp for
        for (int i = 0; i < size; i++)
        {
            const double x = Hartree_Fock.get_atom_coordinate(0, 0) + (i + min) * dx;
            const double y = Hartree_Fock.get_atom_coordinate(0, 1) + (i + min) * dy;
            const double z = Hartree_Fock.get_atom_coordinate(0, 2) + (i + min) * dz;
            HF_densities[i] = Hartree_Fock.compute_dens({ x, y, z }, d, phi);
            DFT_densities[i] = DFT.compute_dens({ x, y, z }, d, phi);
            double temp = abs(x - Pos_C[0]);
            C_dens[i] = C.get_radial_density(temp);
            total_dens[i] = C_dens[i];
            temp = abs(x - Pos_H1[0]);
            H_dens[i] = H.get_radial_density(temp);
            total_dens[i] += H_dens[i];
            temp = sqrt(pow(x - Pos_O1[0], 2) + pow(y - Pos_O1[1], 2) + pow(z - Pos_O1[2], 2));
            total_dens[i] += O.get_radial_density(temp);
            temp = sqrt(pow(x - Pos_O3[0], 2) + pow(y - Pos_O3[1], 2) + pow(z - Pos_O3[2], 2));
            total_dens[i] += O.get_radial_density(temp);
            auto res = get_integration_weights(4, charges.data(), x_coords.data(), y_coords.data(), z_coords.data(), 0, x, y, z, pa_b, pa_tv, chi_HF);
            B_weights_C[i] = res[0];
            TFVC_weights_C[i] = res[1];
            res = get_integration_weights(4, charges.data(), x_coords.data(), y_coords.data(), z_coords.data(), 2, x, y, z, pa_b, pa_tv, chi_HF);
            B_weights_H[i] = res[0];
            TFVC_weights_H[i] = res[1];
            res = get_integration_weights(4, charges.data(), x_coords.data(), y_coords.data(), z_coords.data(), 0, x, y, z, pa_b, pa_tv, chi_DFT);
            TFVC_weights_C_DFT[i] = res[1];
            res = get_integration_weights(4, charges.data(), x_coords.data(), y_coords.data(), z_coords.data(), 2, x, y, z, pa_b, pa_tv, chi_DFT);
            TFVC_weights_H_DFT[i] = res[1];
            pb->update();
        }
    }
    delete pb;

    GridConfiguration config;
    config.all_charges = true;
    GridManager grid_HF(config);
    GridManager grid_DFT(config);

    grid_HF.setup3DGridsForMolecule(Hartree_Fock, { 0, 1, 2, 3 });
    auto res_HF = grid_HF.calculatePartitionedCharges(Hartree_Fock);

    grid_DFT.setup3DGridsForMolecule(DFT, { 0, 1, 2, 3 });
    auto res_DFT = grid_DFT.calculatePartitionedCharges(DFT);

    std::cout << "\n\nHartree fock charges: " << std::endl;
    grid_HF.printChargeTable({ "C1", "O1", "H1", "O3" }, Hartree_Fock, { 0, 1, 2, 3 }, std::cout, res_HF);
    std::cout << "\n\n\nDFT charges: " << std::endl;
    grid_DFT.printChargeTable({ "C1", "O1", "H1", "O3" }, DFT, { 0, 1, 2, 3 }, std::cout, res_DFT);

    ofstream result("densities.dat", ios::out);
    result << setw(10) << "x" << setw(16) << "HF_density" << setw(16) << "DFT_density" << setw(16) << "C_density" << setw(16) << "H_density" << setw(16) << "total_density"
        << setw(16) << "B_weight_C" << setw(16) << "B_weight_H"
        << setw(16) << "TFVC_weight_C_HF" << setw(16) << "TFVC_weight_H_HF"
        << setw(16) << "TFVC_weight_C_DFT" << setw(16) << "TFVC_weight_H_DFT"
        << endl;
    for (int i = 0; i < size; i++)
    {
        result << setw(10) << setprecision(4) << scientific << (i + min) * dx
            << setw(16) << setprecision(8) << scientific << HF_densities[i]
            << setw(16) << setprecision(8) << scientific << DFT_densities[i]
            << setw(16) << setprecision(8) << scientific << C_dens[i]
            << setw(16) << setprecision(8) << scientific << H_dens[i]
            << setw(16) << setprecision(8) << scientific << total_dens[i]
            << setw(16) << setprecision(8) << scientific << B_weights_C[i]
            << setw(16) << setprecision(8) << scientific << B_weights_H[i]
            << setw(16) << setprecision(8) << scientific << TFVC_weights_C[i]
            << setw(16) << setprecision(8) << scientific << TFVC_weights_H[i]
            << setw(16) << setprecision(8) << scientific << TFVC_weights_C_DFT[i]
            << setw(16) << setprecision(8) << scientific << TFVC_weights_H_DFT[i]
            << endl;
    }
    result.flush();
    result.close();
    exit(0);
};

cdouble calc_spherically_averaged_at_k(vec2 &d1, vec2 &d2, vec2 &d3, vec2 &dens,
    const double &k)
{
    double rho, work;
    double *d1_local, *d2_local, *d3_local, *dens_local;
    const int max_num_angular_points_closest = constants::get_closest_num_angular(5810);

    int num_angular = max_num_angular_points_closest;
    vec angular_x(num_angular, 0.0);
    vec angular_y(num_angular, 0.0);
    vec angular_z(num_angular, 0.0);
    vec angular_w(num_angular, 0.0);
    lebedev_sphere ls;
    ls.ld_by_order(max_num_angular_points_closest,
        angular_x.data(),
        angular_y.data(),
        angular_z.data(),
        angular_w.data());


    dens_local = dens[0].data();
    d1_local = d1[0].data();
    d2_local = d2[0].data();
    d3_local = d3[0].data();
    const int dsize = dens[0].size();
    double re = 0, im = 0, c, si, k1, k2, k3;
#pragma omp parallel for reduction(+ : re, im)
    for (int p = 0; p < num_angular; p++)
    {
        k1 = k * angular_x[p];
        k2 = k * angular_y[p];
        k3 = k * angular_z[p];
        for (int d = 0; d < dsize; d++) {
            rho = dens_local[d];
            work = k1 * d1_local[d] + k2 * d2_local[d] + k3 * d3_local[d];
            c = cos(work);
            si = sin(work);
            re += rho * c * angular_w[p];
            im += rho * si * angular_w[p];
        }
    }
    return { re, im };
}

void spherically_averaged_density(options &opt, const ivec val_els_alpha, const ivec val_els_beta)
{
    std::cout << "Reading wavefunction" << std::endl;
    using namespace std;
    WFN wavy(opt.wfn, opt.charge, opt.mult);
    std::cout << "Number of MOs before: " << wavy.get_nmo() << endl;
    wavy.delete_unoccupied_MOs();
    std::cout << "Number of occupied MOs before: " << wavy.get_nmo() << endl;
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
                std::cout << "Deleting from Alpha: " << i << endl;
                wavy.delete_MO(i);
                MOs_to_delete[i] = true;
                deleted++;
            }
        offset = wavy.get_MO_op_count(0);
        for (int i = wavy.get_nmo() - 1; i >= offset; i--)
            if (find(val_els_beta.begin(), val_els_beta.end(), i - offset) == val_els_beta.end())
            {
                std::cout << "Deleting from Beta: " << i - offset << endl;
                wavy.delete_MO(i);
                MOs_to_delete[i + deleted] = true;
                deleted++;
            }
    }
    std::cout << "MOs deleted: " << deleted << endl;
    std::cout << "MO map:" << endl;
    for (int i = 0; i < MOs_to_delete.size(); i++)
        std::cout << i << " " << MOs_to_delete[i] << endl;
    std::cout << "Number of MOs after: " << wavy.get_nmo() << endl;
    std::cout << "\n\nEnergies after:" << endl;
    for (int i = 0; i < wavy.get_nmo(); i++)
        std::cout << wavy.get_MO_energy(i) << " " << wavy.get_MO_occ(i) << endl;

    // Make radial grids on logarithmic scale
    const double dr = 0.00008;
    const long long int upper_r = 250000;
    // Calcualte density on angular grid at each point of radial grid, average and integrate
    long double tot_int2 = 0;
    vec radial_dens2(upper_r, 0.0);
    ProgressBar *progress = new ProgressBar(upper_r, 85, "=", " ", "Calculating Densities");
    std::cout << endl;

#pragma omp parallel for reduction(+ : tot_int2) num_threads(opt.threads)
    for (long long int _r = 1; _r < upper_r; _r++)
    {
        double r = _r * dr;
        radial_dens2[_r] = calc_grid_averaged_at_r(wavy, r, 1200, 5800);
        if (_r >= 1)
        {
            tot_int2 += (long double)radial_dens2[_r] * r * r * (r - (_r - 1) * dr);
        }
        progress->update();
    }
    delete (progress);
    std::cout << "Start writing the file" << endl;
    string el = constants::atnr2letter(wavy.get_atom_charge(0));
    ofstream out(el + ".dat", ios::out);
    out << "Total Integral: " << setw(18) << scientific << setprecision(10) << tot_int2 << "\n";
    for (int i = 0; i < upper_r; i++)
        out << setw(24) << scientific << setprecision(15) << i * dr << setw(24) << scientific << setprecision(16) << radial_dens2[i] / constants::FOUR_PI << "\n";
    out.flush();
    out.close();
}

void spherical_harmonic_test()
{
    const double phi = 0.3, theta = 0.4;
    for (int lam = 0; lam <= 5; lam++)
    {
        for (int m = -lam; m <= lam; m++)
        {
            vec d{ std::sin(theta) * std::cos(phi),
                  std::sin(theta) * std::sin(phi),
                  std::cos(theta), 1.0, 1.0 };
            std::cout << constants::spherical_harmonic(lam, m, d.data()) << " ";
        }
        std::cout << "\n";
    }
};

// Only valid for one atom positioned at 0,0,0
double compute_MO_spherical_orig(double x, double y, double z, double expon, double coef, int type)
{
    int l[3]{ 0, 0, 0 };
    double ex = 0;
    double temp = 0;

    // x, y, z and dsqd
    vec d{
        x,
        y,
        z,
        x * x + y * y + z * z };

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
    int l[3]{ 0, 0, 0 };
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
                    i * CubeMO.get_vector(2, 0) + j * CubeMO.get_vector(2, 1) + k * CubeMO.get_vector(2, 2) + CubeMO.get_origin(2) };

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

void cube_from_coef_npy(vec &coefs, WFN &wave)
{
    cube res = calc_cube_ML(coefs, wave);
    res.set_path((wave.get_path().parent_path() / wave.get_path().stem()).string() + "_PySCF_COEFS_rho.cube");
    res.write_file(true);
    // for (int i = 0; i < dummy.get_ncen(); i++)
    //     calc_cube_ML(data, dummy, i);
}


void get1DGridData(WFN &wavy, std::vector<std::shared_ptr<BasisSet>> &aux_basis, const int atom_idx_1, const int atom_idx_2, const int gridpoints = 1000, const double padding = 2.0) {
    GridConfiguration config;
    config.accuracy = 0;
    config.partition_type = PartitionType::TFVC;
    config.pbc = 0;
    config.debug = false;
    GridManager grid_manager(config);

    grid_manager.setup1DGridsForMolecule(wavy, atom_idx_1, atom_idx_2, gridpoints, padding);
    const vec3 atomic_grids_local = grid_manager.getGridData().atomic_grids;


    WFN wavy_aux = generate_aux_wfn(wavy, aux_basis);
    DensityFitting::CONFIG RI_config;
    RI_config.analyze_quality = true;
    vec ri_coefs_u = DensityFitting::density_fit(wavy, wavy_aux, RI_config);

    RI_config.restrain_type = DensityFitting::RESTRAINT_TYPE::SIMPLE_AND_TIK;
    RI_config.charge_scheme = DensityFitting::CHARGE_SCHEME::TFVC;
    RI_config.restraint_strength = 2.0E-4;
    RI_config.tikhonov_lambda = 1E-6;
    vec ri_coefs_tfvc = DensityFitting::density_fit(wavy, wavy_aux, RI_config);

    vec2 dens_hirsh(2);
    vec2 dens_becke(2);
    vec2 dens_tfvc(2);
    vec2 dens_RI_u(2);
    vec2 dens_RI_tfvc(2);
    ivec atom_indices = { atom_idx_1, atom_idx_2 };
    for (int g = 0; g < 2; ++g) {
        dens_hirsh[g].resize(gridpoints), dens_becke[g].resize(gridpoints), dens_tfvc[g].resize(gridpoints), dens_RI_u[g].resize(gridpoints), dens_RI_tfvc[g].resize(gridpoints);
        for (int p = 0; p < gridpoints; ++p) {
            dens_hirsh[g][p] = atomic_grids_local[g][GridData::WFN_DENSITY][p] * atomic_grids_local[g][GridData::HIRSH_WEIGHT][p];
            dens_becke[g][p] = atomic_grids_local[g][GridData::WFN_DENSITY][p] * atomic_grids_local[g][GridData::BECKE_WEIGHT][p];
            dens_tfvc[g][p] = atomic_grids_local[g][GridData::WFN_DENSITY][p] * atomic_grids_local[g][GridData::TFVC_WEIGHT][p];
            dens_RI_u[g][p] = calc_density_ML(
                atomic_grids_local[g][GridData::X][p],
                atomic_grids_local[g][GridData::Y][p],
                atomic_grids_local[g][GridData::Z][p],
                ri_coefs_u,
                wavy_aux.get_atoms(),
                atom_indices[g]);
            dens_RI_tfvc[g][p] = calc_density_ML(
                atomic_grids_local[g][GridData::X][p],
                atomic_grids_local[g][GridData::Y][p],
                atomic_grids_local[g][GridData::Z][p],
                ri_coefs_tfvc,
                wavy_aux.get_atoms(),
                atom_indices[g]);
        }
    }

    std::vector<std::pair<std::string, vec>> data_atom1(
        {
           {"hirsh", dens_hirsh[0]},
           {"becke", dens_becke[0]},
           {"tfvc", dens_tfvc[0]},
           {"ri_u", dens_RI_u[0]},
           { "ri_tfvc", dens_RI_tfvc[0] }
        }
    );
    std::vector<std::pair<std::string, vec>> data_atom2(
        {
            {"hirsh", dens_hirsh[1]},
            {"becke", dens_becke[1]},
            {"tfvc", dens_tfvc[1]},
            {"ri_u", dens_RI_u[1]},
            { "ri_tfvc", dens_RI_tfvc[1] }
        }
    );

    //std::string atom_specifier = std::format("{}{}_{}{}", wavy.get_atom_label(atom_idx_1), atom_idx_1, wavy.get_atom_label(atom_idx_2), atom_idx_2);
    std::string atom_specifier = wavy.get_atom_label(atom_idx_1) + std::to_string(atom_idx_1) + "_" + wavy.get_atom_label(atom_idx_2) + std::to_string(atom_idx_2);
    vec2 grid_points = { atomic_grids_local[0][GridData::X], atomic_grids_local[0][GridData::Y], atomic_grids_local[0][GridData::Z] };
    grid_manager.writeSimpleGrid("WFN_density_" + atom_specifier + ".dat", grid_points, { {"WFN_density",atomic_grids_local[0][GridData::WFN_DENSITY]} });
    grid_manager.writeSimpleGrid("Density_" + wavy.get_atom_label(atom_idx_1) + std::to_string(atom_idx_1) + ".dat", grid_points, data_atom1);
    grid_manager.writeSimpleGrid("Density_" + wavy.get_atom_label(atom_idx_2) + std::to_string(atom_idx_2) + ".dat", grid_points, data_atom2);
}

void draw_orbital(const int lambda, const int m, const double resulution = 0.025, const double radius = 3.5)
{
    if (m > lambda || m < -lambda)
    {
        std::cout << "m must be between -l and l" << std::endl;
        return;
    }

    // Initialize the Wavefunction
    WFN wavy(e_origin::NOT_YET_DEFINED);
    wavy.push_back_MO(0, 1.0, -13);
    wavy.push_back_atom("H", 0, 0, 0, 1);
    wavy.push_back_atom_basis_set(0, 1.0, 1.0, lambda, 0);
    properties_options opts;
    opts.radius = radius;
    opts.resolution = resulution;
    readxyzMinMax_fromWFN(wavy, opts, true);
    cube CubeMO(opts.NbSteps, 1, true);
    CubeMO.give_parent_wfn(wavy);
    for (int i = 0; i < 3; i++)
    {
        CubeMO.set_origin(i, opts.MinMax[i]);
        CubeMO.set_vector(i, i, (opts.MinMax[i + 3] - opts.MinMax[i]) / opts.NbSteps[i]);
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
                    i * CubeMO.get_vector(2, 0) + j * CubeMO.get_vector(2, 1) + k * CubeMO.get_vector(2, 2) + CubeMO.get_origin(2) };

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
    //std::shared_ptr<BasisSet> aux_basis_set = std::make_shared<BasisSet>(wavy);
    std::vector<std::shared_ptr<BasisSet>> aux_basis_set{ BasisSetLibrary::get_basis_set("def2-universal-jkfit") };

    WFN wavy_aux = generate_aux_wfn(wavy, aux_basis_set);

    DensityFitting::CONFIG RI_config;
    RI_config.restrain_type = DensityFitting::RESTRAINT_TYPE::SIMPLE_AND_TIK;
    RI_config.analyze_quality = opt->debug;

    vec ri_coefs = density_fit(wavy, wavy_aux, RI_config);


    bvec needs_grid(wavy.get_ncen(), false);
    needs_grid[0] = true;
    GridConfiguration conf;
    conf.partition_type = PartitionType::Hirshfeld;
    conf.accuracy = 4;
    GridManager grid(conf);
    vec2 d1, d2, d3, dens;

    ivec asym_atom_list(1, 0);
    grid.setup3DGridsForMolecule(wavy, asym_atom_list, needs_grid);
    grid.getDensityVectors(wavy, asym_atom_list, d1, d2, d3, dens);
    vec2 ML_grid(wavy.get_ncen());
    auto grid_data = grid.getGridData();
    int size;
#pragma omp parallel
    {
        vec2 d_temp(wavy.get_ncen());
        for (int i = 0; i < wavy.get_ncen(); i++)
        {
            d_temp[i].resize(16, 0.0);
        }
        vec phi_temp(wavy.get_nmo(), 0.0);

        for (int a = 0; a < wavy.get_ncen(); a++) {
            size = grid.getNumPointsForAtom(a);
            ML_grid[a].resize(size, 0.0);
#pragma omp for
            for (int i = 0; i < size; i++)
            {
                ML_grid[a][i] = calc_density_ML(
                    d1[a][i],
                    d2[a][i],
                    d3[a][i],
                    ri_coefs,
                    wavy_aux.get_atoms()) * grid_data.atomic_grids[a][GridData::GridIndex::WEIGHT][i];
            }
        }
        for (int i = 0; i < 16; i++)
            shrink_vector<double>(d_temp[i]);
        shrink_vector<vec>(d_temp);
        shrink_vector<double>(phi_temp);
    }
    vec elecs_DFT(wavy.get_ncen(), 0.0);
    auto charge_results = grid.calculatePartitionedCharges(wavy);
    vec elecs_RI(wavy.get_ncen(), 0.0);
    for (int a = 0; a < wavy.get_ncen(); a++)
    {
        elecs_DFT[a] = charge_results.atom_charges[a][PartitionResults::CHARGE_ORDER::S_HIRSH];
        elecs_RI[a] = vec_sum(ML_grid[a]);
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

}

// Convert a row-major matrix to column-major format
std::vector<double> rowToColMajor(const std::vector<double> &rowMajorMatrix, int rows, int cols)
{
    std::vector<double> colMajorMatrix((size_t)rows * cols);

    for (int r = 0; r < rows; ++r)
    {
        for (int c = 0; c < cols; ++c)
        {
            colMajorMatrix[(size_t)c * rows + r] = rowMajorMatrix[(size_t)r * cols + c];
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
        67.06378696, 21.03825611, 12.89262977, 31.54283509, 36.37107709 };
    //Dimensions
    const int m = 10; // Number of rows
    const int n = 5;  // Number of columns
    //Right-hand side vector B
    std::vector<double> B_in = { 57.01967704, 43.86015135, 98.83738381, 10.20448107, 20.88767561, 16.13095179, 65.31083255, 25.32916025, 46.63107729, 24.4425592 };

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
    std::vector<double> correct = { 0.00370583, 0.58600469, 0.14755389, 0.04786599, 0.0 };
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