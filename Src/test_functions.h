#pragma once
#include "spherical_density.h"
#include "scattering_factors.h"
#include "AtomGrid.h"
#include "convenience.h"
#include "npy.h"
#include "properties.h"
#include "JKFit.h"
#include "SALTED_utilities.h"
#include "cerf.h"
#if has_RAS
#include "SALTED_math.h"
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

void test_cerf() {
    using namespace std;
    cerf(cdouble(1.0, 1.0));
    ofstream name("test.dat", ios::out);
    for (int i = 0; i < 1000; i++) {
        name << setprecision(8) << scientific << setw(16) << i * 0.01;
        name << setprecision(8) << scientific << setw(16) << cerf(cdouble(i * 0.01, 0.0)).real();
        name << setprecision(8) << scientific << setw(16) << cerf(cdouble(i * 0.01, 0.0)).imag();
        name << setprecision(8) << scientific << setw(16) << cerf(cdouble(i * 0.01, 1.0)).real();
        name << setprecision(8) << scientific << setw(16) << cerf(cdouble(i * 0.01, 1.0)).imag();
        name << setprecision(8) << scientific << setw(16) << cerf(cdouble(0.0, i * 0.01)).real();
        name << setprecision(8) << scientific << setw(16) << cerf(cdouble(0.0, i * 0.01)).imag();
        name << setprecision(8) << scientific << setw(16) << cerf(cdouble(1.0, i * 0.01)).real();
        name << setprecision(8) << scientific << setw(16) << cerf(cdouble(1.0, i * 0.01)).imag();
        name << endl;
    }
    name.close();
    exit(0);
}

void test_density_cubes(options &opt, std::ostream &log_file)
{
    using namespace std;
    vector<WFN> wavy(10);
    // ScF2+ test file against ORCA calcualted cubes
    log_file << "====================ScF2+ Test===============================" << endl;
    wavy[0].read_known_wavefunction_format("test.molden", log_file);
    cube Rho(100, 100, 100, wavy[0].get_ncen(), true);
    Rho.set_origin(0, -7), Rho.set_origin(1, -7), Rho.set_origin(2, -7);
    Rho.set_vector(0, 0, 0.141414);
    Rho.set_vector(1, 1, 0.141414);
    Rho.set_vector(2, 2, 0.141414);
    Rho.path = get_basename_without_ending(wavy[0].get_path()) + "_rho.cube";
    log_file << "Calcualting Rho...";
    Calc_Rho(Rho, wavy[0], opt.threads, 7.0, log_file);
    log_file << " ...done!" << endl;
    // Rho.write_file(wavy[0], true);
    const double test_molden = Rho.sum();
    log_file << "Number of electrons in the cube: " << setprecision(4) << fixed << test_molden << endl;
    cube Rho2("test.eldens.cube", true, wavy[0], log_file);
    const double test_ref = Rho2.sum();
    log_file << "Number of electrons in the reference cube: " << setprecision(4) << fixed << test_ref << endl;
    for (int i = 0; i < 1; i++)
    {
        cube MO(100, 100, 100, wavy[0].get_ncen(), true);
        MO.set_origin(0, -7), MO.set_origin(1, -7), MO.set_origin(2, -7);
        MO.set_vector(0, 0, 0.141414);
        MO.set_vector(1, 1, 0.141414);
        MO.set_vector(2, 2, 0.141414);
        MO.path = get_basename_without_ending(wavy[0].get_path()) + "_MO_" + to_string(i) + ".cube";
        log_file << "Calcualting MO " + to_string(i) + "...";
        Calc_MO(MO, i, wavy[0], opt.threads, 4.0, log_file);
        log_file << " ...done!" << endl;
        // MO.write_file(wavy[0], true);
        string name("test.mo" + to_string(i) + "a.cube");
        cube MO2(name, true, wavy[0], log_file);
        log_file << "sum in the cube: " << setprecision(4) << fixed << MO.sum() << endl;
        log_file << "sum in the reference cube: " << setprecision(4) << fixed << MO2.sum() << endl;
    }

    // F- ion calculations
    log_file << "====================F Test===============================" << endl;
    wavy[1].read_known_wavefunction_format("F_full.molden", log_file);
    wavy[1].write_wfn("F_conv.wfn", false, true);
    cube Rho_2(71, 71, 71, wavy[1].get_ncen(), true);
    Rho_2.set_origin(0, -7), Rho_2.set_origin(1, -7), Rho_2.set_origin(2, -7);
    Rho_2.set_vector(0, 0, 0.2);
    Rho_2.set_vector(1, 1, 0.2);
    Rho_2.set_vector(2, 2, 0.2);
    Rho_2.path = get_basename_without_ending(wavy[1].get_path()) + "_rho.cube";
    log_file << "Calcualting Rho...";
    Calc_Rho(Rho_2, wavy[1], opt.threads, 7.0, log_file);
    log_file << " ...done!" << endl;
    // Rho_2.write_file(wavy[1], true);
    const double F_molden = Rho_2.sum();
    log_file << "Number of electrons in the cube: " << setprecision(4) << fixed << F_molden << endl;
    wavy[4].read_known_wavefunction_format("f_ref.wfx", log_file);
    cube Rho_3(71, 71, 71, wavy[4].get_ncen(), true);
    Rho_3.set_origin(0, -7), Rho_3.set_origin(1, -7), Rho_3.set_origin(2, -7);
    Rho_3.set_vector(0, 0, 0.2);
    Rho_3.set_vector(1, 1, 0.2);
    Rho_3.set_vector(2, 2, 0.2);
    Rho_3.path = get_basename_without_ending(wavy[4].get_path()) + "_rho.cube";
    log_file << "Calcualting Rho...";
    Calc_Rho_spherical_harmonics(Rho_3, wavy[4], opt.threads, log_file);
    log_file << " ...done!" << endl;
    // Rho_3.write_file(wavy[1], true);
    const double F_ref = Rho_3.sum();
    log_file << "Number of electrons in the reference cube: " << setprecision(4) << fixed << F_ref << endl;

    // Ce conevrsion for test of f type
    log_file << "====================Ce Test===============================" << endl;
    wavy[2].read_known_wavefunction_format("Ce_full.molden", log_file);
    wavy[2].write_wfn("Ce_conv.wfn", false, true);
    cube Rho_4(141, 141, 141, wavy[2].get_ncen(), true);
    Rho_4.set_origin(0, -7), Rho_4.set_origin(1, -7), Rho_4.set_origin(2, -7);
    Rho_4.set_vector(0, 0, 0.1);
    Rho_4.set_vector(1, 1, 0.1);
    Rho_4.set_vector(2, 2, 0.1);
    Rho_4.path = get_basename_without_ending(wavy[2].get_path()) + "_rho.cube";
    log_file << "Calcualting Rho...";
    Calc_Rho(Rho_4, wavy[2], opt.threads, 7.0, log_file);
    log_file << " ...done!" << endl;
    // Rho_4.write_file(wavy[2], true);
    const double Ce_molden = Rho_4.sum();
    log_file << "Number of electrons in the cube: " << setprecision(4) << fixed << Ce_molden << endl;

    Rho_4.set_zero();
    wavy[5].read_known_wavefunction_format("Ce_full.wfn", log_file);
    Calc_Rho(Rho_4, wavy[5], opt.threads, 7.0, log_file);
    Rho_4.path = get_basename_without_ending(wavy[5].get_path()) + "_reference_rho.cube";
    // Rho_4.write_file(wavy[5], true);
    const double Ce_ref = Rho_4.sum();
    log_file << "Number of electrons in the ref cube: " << setprecision(4) << fixed << Ce_ref << endl;

    // Sc conversion
    log_file << "====================Sc Test===============================" << endl;
    wavy[3].read_known_wavefunction_format("Sc_full.molden", log_file);
    wavy[3].write_wfn("Sc_conv.wfn", false, true);

    // Lu g-type test
    log_file << "====================Lu Test===============================" << endl;
    wavy[6].read_known_wavefunction_format("Lu_jorge.molden", log_file);
    wavy[6].write_wfn("Lu_conv.wfn", false, true);
    cube Rho4(141, 141, 141, wavy[6].get_ncen(), true);
    Rho4.set_origin(0, -7), Rho4.set_origin(1, -7), Rho4.set_origin(2, -7);
    Rho4.set_vector(0, 0, 0.1);
    Rho4.set_vector(1, 1, 0.1);
    Rho4.set_vector(2, 2, 0.1);
    Calc_Rho(Rho4, wavy[6], opt.threads, 7.0, log_file);
    Rho4.path = get_basename_without_ending(wavy[6].get_path()) + "_rho.cube";
    // Rho4.write_file(wavy[6], true);
    const double Lu_molden = Rho4.sum();
    log_file << "Number of electrons in cube: " << setprecision(4) << fixed << Lu_molden << endl;

    wavy[7].read_known_wavefunction_format("Lu_jorge.wfn", log_file);
    cube Rho_5(141, 141, 141, wavy[7].get_ncen(), true);
    Rho_5.set_origin(0, -7), Rho_5.set_origin(1, -7), Rho_5.set_origin(2, -7);
    Rho_5.set_vector(0, 0, 0.1);
    Rho_5.set_vector(1, 1, 0.1);
    Rho_5.set_vector(2, 2, 0.1);
    Calc_Rho(Rho_5, wavy[7], opt.threads, 7.0, log_file);
    Rho_5.path = get_basename_without_ending(wavy[7].get_path()) + "_reference_rho.cube";
    // Rho_5.write_file(wavy[7], true);
    const double Lu_wfn = Rho_5.sum();
    Rho_5 -= Rho4;
    Rho_5.path = get_basename_without_ending(wavy[7].get_path()) + "_diff_rho.cube";
    // Rho_5.write_file(wavy[7], true);
    log_file << "Number of electrons in the ref cube: " << setprecision(4) << fixed << Lu_wfn << endl;
    log_file << "Number of electrons in the diff cube: " << setprecision(4) << fixed << Rho_5.diff_sum() << endl;

    wavy[8].read_known_wavefunction_format("Lu_def2.molden", log_file);
    wavy[8].write_wfn("Lu_def2_conv.wfn", false, true);
    wavy[8].set_has_ECPs(true);
    cube Rho5(141, 141, 141, wavy[8].get_ncen(), true);
    Rho5.set_origin(0, -7), Rho5.set_origin(1, -7), Rho5.set_origin(2, -7);
    Rho5.set_vector(0, 0, 0.1);
    Rho5.set_vector(1, 1, 0.1);
    Rho5.set_vector(2, 2, 0.1);
    Calc_Rho(Rho5, wavy[8], opt.threads, 7.0, log_file);
    Rho5.path = get_basename_without_ending(wavy[8].get_path()) + "_rho.cube";
    // Rho4.write_file(wavy[6], true);
    const double Lu_def2 = Rho5.sum();
    log_file << "Number of electrons in cube: " << setprecision(4) << fixed << Lu_def2 << endl;

    err_checkf(abs(test_molden - test_ref) < 0.1, "Difference in test too big!", log_file);
    err_checkf(abs(F_molden - F_ref) < 0.1, "Difference in F too big!", log_file);
    err_checkf(abs(Ce_molden - Ce_ref) < 0.1, "Difference in Ce too big!", log_file);
    err_checkf(abs(Lu_molden - Lu_wfn) < 0.1, "Difference in Lu too big!", log_file);
    log_file << "All tests successfull!" << endl;
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
    const int step = max((int)floor(smax / 20), 1);
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
    progress_bar *progress = new progress_bar{std::cout, 60u, "Calculating scattering factors"};
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
            if (s != 0 && s % step == 0)
                progress->write(s / static_cast<double>(smax));
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

    while (ratio > rel_precision && angle / constants::PI < 1E4)
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
        }
        angle *= 1.2;
    }
    if (print)
        std::cout << "Done with " << std::setw(10) << std::setprecision(5) << r << " " << (ratio > rel_precision) << std::setw(14) << angle << std::endl;
    return new_result;
}

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
    using namespace std;
    WFN wavy(9);
    wavy.read_known_wavefunction_format(opt.wfn, cout, opt.debug);
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
    long double tot_int = 0;
    vec radial_dens(upper_r, 0.0);
    progress_bar *progress = new progress_bar{cout, 60u, "Calculating densities"};
    const long long int step = max(static_cast<long long int>(floor(upper_r / 20)), 1LL);

#pragma omp parallel for reduction(+ : tot_int) num_threads(opt.threads) schedule(dynamic)
    for (long long int _r = 1; _r < upper_r; _r++)
    {
        double r = _r * dr;
        double v = calc_spherically_averaged_at_r(wavy, r, 1E-4, 20);
        radial_dens[_r] = v;
        if (_r >= 1)
            tot_int += v * r * r * (r - (_r - 1) * dr);
        if (_r != 0 && _r % step == 0)
            progress->write(_r / static_cast<double>(upper_r));
    }
    delete (progress);
    cout << "Start writing the file" << endl;
    string el = constants::atnr2letter(wavy.get_atom_charge(0));
    ofstream out(el + ".dat", ios::out);
    out << "Total Integral: " << setw(18) << scientific << setprecision(10) << tot_int * constants::FOUR_PI << "\n";
    for (int i = 0; i < upper_r; i++)
        out << setw(24) << scientific << setprecision(15) << i * dr << setw(24) << scientific << setprecision(16) << radial_dens[i] << "\n";
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

void sfac_scan_ECP(options &opt, std::ostream &log_file)
{
#ifdef _OPENMP
    omp_set_num_threads(opt.threads);
#endif
    using namespace std;
    std::vector<WFN> wavy;
    wavy.emplace_back(1);
    wavy[0].read_known_wavefunction_format(opt.wfn, std::cout, opt.debug);
    wavy.emplace_back(9);
    wavy[1].read_known_wavefunction_format("Au_def2.gbw", std::cout, opt.debug);
    wavy.emplace_back(9);
    wavy[2].read_known_wavefunction_format("Au_alle.gbw", std::cout, opt.debug);
    Thakkar Au(wavy[0].atoms[0].charge);
    if (opt.ECP)
    {
        if (opt.ECP_mode == 2)
        {

            wavy[0].set_has_ECPs(true, true, true);
            wavy[1].set_has_ECPs(true, true, true);
        }
        else
        {
            wavy[0].set_has_ECPs(true);
            wavy[1].set_has_ECPs(true);
        }
    }
    err_checkf(wavy[0].get_ncen() != 0, "No Atoms in the wavefunction, this will not work!! ABORTING!!", std::cout);
    err_checkf(exists(opt.cif), "CIF does not exists!", std::cout);

    std::vector<time_point> time_points;
    std::vector<std::string> time_descriptions;

    cell unit_cell(opt.cif, log_file, opt.debug);
    ifstream cif_input(opt.cif.c_str(), std::ios::in);
    ivec atom_type_list;
    ivec asym_atom_to_type_list;
    ivec asym_atom_list;
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
                                      log_file,
                                      opt.debug);

    cif_input.close();
    vec2 d1, d2, d3, dens;

    make_hirshfeld_grids(opt.pbc,
                         4,
                         unit_cell,
                         wavy[0],
                         atom_type_list,
                         asym_atom_list,
                         needs_grid,
                         d1, d2, d3, dens,
                         labels,
                         log_file,
                         time_points,
                         time_descriptions,
                         opt.debug,
                         opt.no_date);

    vec2 d1_def2, d2_def2, d3_def2, dens_def2;
    ivec atl_def2{79};
    ivec aal_def2{0};
    bvec ng_def2(1, true);

    make_hirshfeld_grids(opt.pbc,
                         4,
                         unit_cell,
                         wavy[1],
                         atl_def2,
                         aal_def2,
                         ng_def2,
                         d1_def2, d2_def2, d3_def2, dens_def2,
                         labels,
                         log_file,
                         time_points,
                         time_descriptions,
                         opt.debug,
                         opt.no_date);

    vec2 d1_all, d2_all, d3_all, dens_all;
    ivec atl_all{79};
    ivec aal_all{0};
    bvec ng_all(1, true);

    make_hirshfeld_grids(opt.pbc,
                         4,
                         unit_cell,
                         wavy[2],
                         atl_all,
                         aal_all,
                         ng_all,
                         d1_all, d2_all, d3_all, dens_all,
                         labels,
                         log_file,
                         time_points,
                         time_descriptions,
                         opt.debug,
                         opt.no_date);

    WFN wavy_all_val(9);
    wavy_all_val.read_known_wavefunction_format("Au_alle.gbw", std::cout, opt.debug);
    for (int i = 0; i < 23; i++)
        wavy_all_val.delete_MO(0);
    for (int i = 0; i < 7; i++)
        wavy_all_val.delete_MO(1);
    vec2 d1_all_val, d2_all_val, d3_all_val, dens_all_val;

    make_hirshfeld_grids(opt.pbc,
                         4,
                         unit_cell,
                         wavy_all_val,
                         atl_all,
                         aal_all,
                         ng_all,
                         d1_all_val, d2_all_val, d3_all_val, dens_all_val,
                         labels,
                         log_file,
                         time_points,
                         time_descriptions,
                         opt.debug,
                         opt.no_date);

    WFN wavy_ZORA(9);
    wavy_ZORA.read_known_wavefunction_format("Au_alle_ZORA.gbw", std::cout, opt.debug);
    vec2 d1_ZORA, d2_ZORA, d3_ZORA, dens_ZORA;

    make_hirshfeld_grids(opt.pbc,
                         4,
                         unit_cell,
                         wavy_ZORA,
                         atl_all,
                         aal_all,
                         ng_all,
                         d1_ZORA, d2_ZORA, d3_ZORA, dens_ZORA,
                         labels,
                         log_file,
                         time_points,
                         time_descriptions,
                         opt.debug,
                         opt.no_date);

    WFN wavy_ZORA_val(9);
    wavy_ZORA_val.read_known_wavefunction_format("Au_alle_ZORA.gbw", std::cout, opt.debug);
    for (int i = 0; i < 23; i++)
        wavy_ZORA_val.delete_MO(0);
    for (int i = 0; i < 7; i++)
        wavy_ZORA_val.delete_MO(1);
    vec2 d1_ZORA_val, d2_ZORA_val, d3_ZORA_val, dens_ZORA_val;

    make_hirshfeld_grids(opt.pbc,
                         4,
                         unit_cell,
                         wavy_ZORA_val,
                         atl_all,
                         aal_all,
                         ng_all,
                         d1_ZORA_val, d2_ZORA_val, d3_ZORA_val, dens_ZORA_val,
                         labels,
                         log_file,
                         time_points,
                         time_descriptions,
                         opt.debug,
                         opt.no_date);

    WFN wavy_x2c(9);
    wavy_x2c.read_known_wavefunction_format("Au_alle_x2c.gbw", std::cout, opt.debug);
    vec2 d1_x2c, d2_x2c, d3_x2c, dens_x2c;

    make_hirshfeld_grids(opt.pbc,
                         4,
                         unit_cell,
                         wavy_x2c,
                         atl_all,
                         aal_all,
                         ng_all,
                         d1_x2c, d2_x2c, d3_x2c, dens_x2c,
                         labels,
                         log_file,
                         time_points,
                         time_descriptions,
                         opt.debug,
                         opt.no_date);

    WFN wavy_x2c_val(9);
    wavy_x2c_val.read_known_wavefunction_format("Au_alle_ZORA.gbw", std::cout, opt.debug);
    for (int i = 0; i < 23; i++)
        wavy_x2c_val.delete_MO(0);
    for (int i = 0; i < 7; i++)
        wavy_x2c_val.delete_MO(1);
    vec2 d1_x2c_val, d2_x2c_val, d3_x2c_val, dens_x2c_val;

    make_hirshfeld_grids(opt.pbc,
                         4,
                         unit_cell,
                         wavy_x2c_val,
                         atl_all,
                         aal_all,
                         ng_all,
                         d1_x2c_val, d2_x2c_val, d3_x2c_val, dens_x2c_val,
                         labels,
                         log_file,
                         time_points,
                         time_descriptions,
                         opt.debug,
                         opt.no_date);

    std::cout << "finished partitioning" << endl;
    const int size = 2000;
    const int phi_size = 30;
    const int theta_size = 30;
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
#pragma omp parallel for schedule(dynamic)
    for (int ref = 1; ref <= size; ref++)
    {
        for (int p = 0; p < phi_size; p++)
        {
            for (int t = 0; t < theta_size; t++)
            {
                int ind = t + (p + (ref - 1) * phi_size) * theta_size;
                double k_length = constants::bohr2ang(constants::FOUR_PI * ref / size * opt.d_sfac_scan);
                k_pt[0][ind] = k_length * sin(t * theta_step) * cos(p * phi_step);
                k_pt[1][ind] = k_length * sin(t * theta_step) * sin(p * phi_step);
                k_pt[2][ind] = k_length * cos(t * theta_step);
                k_pt[3][ind] = k_length;
            }
        }
    }
    // below is a strip of Calc_SF without the file IO or progress bar
    cvec2 sf;
    cvec2 sf_def2;
    cvec2 sf_all;
    cvec2 sf_all_val;
    cvec2 sf_ZORA;
    cvec2 sf_ZORA_val;
    cvec2 sf_x2c;
    cvec2 sf_x2c_val;
    const int smax = (int)k_pt[0].size();
    int pmax = (int)dens[0].size();
    std::cout << "Done with making k_pt " << smax << " " << pmax << endl;
    sf.resize(1);
    sf_def2.resize(1);
    sf_all.resize(1);
    sf_all_val.resize(1);
    sf_ZORA.resize(1);
    sf_ZORA_val.resize(1);
    sf_x2c.resize(1);
    sf_x2c_val.resize(1);
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < 1; i++)
    {
        sf[i].resize(k_pt[0].size());
        sf_def2[i].resize(k_pt[0].size());
        sf_all[i].resize(k_pt[0].size());
        sf_all_val[i].resize(k_pt[0].size());
        sf_ZORA[i].resize(k_pt[0].size());
        sf_ZORA_val[i].resize(k_pt[0].size());
        sf_x2c[i].resize(k_pt[0].size());
        sf_x2c_val[i].resize(k_pt[0].size());
    }
    double *dens_local, *d1_local, *d2_local, *d3_local;
    complex<double> *sf_local;
    const double *k1_local = k_pt[0].data();
    const double *k2_local = k_pt[1].data();
    const double *k3_local = k_pt[2].data();
    double work, rho;
    progress_bar *progress = new progress_bar{std::cout, 60u, "Calculating scattering factors"};
    for (int i = 0; i < 1; i++)
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
#ifdef __APPLE__
#if TARGET_OS_MAC
                if (rho < 0)
                {
                    rho = -rho;
                    work += M_PI;
                }
#endif
#endif
                sf_local[s] += polar(rho, work);
            }
        }
        log_file << "Done with HAR SFs" << endl;
        pmax = (int)dens_def2[i].size();
        dens_local = dens_def2[i].data();
        d1_local = d1_def2[i].data();
        d2_local = d2_def2[i].data();
        d3_local = d3_def2[i].data();
        sf_local = sf_def2[i].data();
#pragma omp parallel for private(work, rho)
        for (int s = 0; s < smax; s++)
        {
            for (int p = pmax - 1; p >= 0; p--)
            {
                rho = dens_local[p];
                work = k1_local[s] * d1_local[p] + k2_local[s] * d2_local[p] + k3_local[s] * d3_local[p];
#ifdef __APPLE__
#if TARGET_OS_MAC
                if (rho < 0)
                {
                    rho = -rho;
                    work += M_PI;
                }
#endif
#endif
                sf_local[s] += polar(rho, work);
            }
        }
        log_file << "Done with def2 SFs" << endl;
        pmax = (int)dens_all[i].size();
        dens_local = dens_all[i].data();
        d1_local = d1_all[i].data();
        d2_local = d2_all[i].data();
        d3_local = d3_all[i].data();
        sf_local = sf_all[i].data();
#pragma omp parallel for private(work, rho)
        for (int s = 0; s < smax; s++)
        {
            for (int p = pmax - 1; p >= 0; p--)
            {
                rho = dens_local[p];
                work = k1_local[s] * d1_local[p] + k2_local[s] * d2_local[p] + k3_local[s] * d3_local[p];
#ifdef __APPLE__
#if TARGET_OS_MAC
                if (rho < 0)
                {
                    rho = -rho;
                    work += M_PI;
                }
#endif
#endif
                sf_local[s] += polar(rho, work);
            }
        }
        log_file << "Done with Jorge SFs" << endl;
        pmax = (int)dens_all_val[i].size();
        dens_local = dens_all_val[i].data();
        d1_local = d1_all_val[i].data();
        d2_local = d2_all_val[i].data();
        d3_local = d3_all_val[i].data();
        sf_local = sf_all_val[i].data();
#pragma omp parallel for private(work, rho)
        for (int s = 0; s < smax; s++)
        {
            for (int p = pmax - 1; p >= 0; p--)
            {
                rho = dens_local[p];
                work = k1_local[s] * d1_local[p] + k2_local[s] * d2_local[p] + k3_local[s] * d3_local[p];
#ifdef __APPLE__
#if TARGET_OS_MAC
                if (rho < 0)
                {
                    rho = -rho;
                    work += M_PI;
                }
#endif
#endif
                sf_local[s] += polar(rho, work);
            }
        }
        log_file << "Done with Jorge Valence SFs" << endl;
        pmax = (int)dens_ZORA[i].size();
        dens_local = dens_ZORA[i].data();
        d1_local = d1_ZORA[i].data();
        d2_local = d2_ZORA[i].data();
        d3_local = d3_ZORA[i].data();
        sf_local = sf_ZORA[i].data();
#pragma omp parallel for private(work, rho)
        for (int s = 0; s < smax; s++)
        {
            for (int p = pmax - 1; p >= 0; p--)
            {
                rho = dens_local[p];
                work = k1_local[s] * d1_local[p] + k2_local[s] * d2_local[p] + k3_local[s] * d3_local[p];
#ifdef __APPLE__
#if TARGET_OS_MAC
                if (rho < 0)
                {
                    rho = -rho;
                    work += M_PI;
                }
#endif
#endif
                sf_local[s] += polar(rho, work);
            }
        }
        log_file << "Done with ZORA-Jorge SFs" << endl;
        pmax = (int)dens_ZORA_val[i].size();
        dens_local = dens_ZORA_val[i].data();
        d1_local = d1_ZORA_val[i].data();
        d2_local = d2_ZORA_val[i].data();
        d3_local = d3_ZORA_val[i].data();
        sf_local = sf_ZORA_val[i].data();
#pragma omp parallel for private(work, rho)
        for (int s = 0; s < smax; s++)
        {
            for (int p = pmax - 1; p >= 0; p--)
            {
                rho = dens_local[p];
                work = k1_local[s] * d1_local[p] + k2_local[s] * d2_local[p] + k3_local[s] * d3_local[p];
#ifdef __APPLE__
#if TARGET_OS_MAC
                if (rho < 0)
                {
                    rho = -rho;
                    work += M_PI;
                }
#endif
#endif
                sf_local[s] += polar(rho, work);
            }
        }
        log_file << "Done with ZORA-Jorge Valence SFs" << endl;
        pmax = (int)dens_x2c[i].size();
        dens_local = dens_x2c[i].data();
        d1_local = d1_x2c[i].data();
        d2_local = d2_x2c[i].data();
        d3_local = d3_x2c[i].data();
        sf_local = sf_x2c[i].data();
#pragma omp parallel for private(work, rho)
        for (int s = 0; s < smax; s++)
        {
            for (int p = pmax - 1; p >= 0; p--)
            {
                rho = dens_local[p];
                work = k1_local[s] * d1_local[p] + k2_local[s] * d2_local[p] + k3_local[s] * d3_local[p];
#ifdef __APPLE__
#if TARGET_OS_MAC
                if (rho < 0)
                {
                    rho = -rho;
                    work += M_PI;
                }
#endif
#endif
                sf_local[s] += polar(rho, work);
            }
        }
        log_file << "Done with X2C SFs" << endl;
        pmax = (int)dens_x2c_val[i].size();
        dens_local = dens_x2c_val[i].data();
        d1_local = d1_x2c_val[i].data();
        d2_local = d2_x2c_val[i].data();
        d3_local = d3_x2c_val[i].data();
        sf_local = sf_x2c_val[i].data();
#pragma omp parallel for private(work, rho)
        for (int s = 0; s < smax; s++)
        {
            for (int p = pmax - 1; p >= 0; p--)
            {
                rho = dens_local[p];
                work = k1_local[s] * d1_local[p] + k2_local[s] * d2_local[p] + k3_local[s] * d3_local[p];
#ifdef __APPLE__
#if TARGET_OS_MAC
                if (rho < 0)
                {
                    rho = -rho;
                    work += M_PI;
                }
#endif
#endif
                sf_local[s] += polar(rho, work);
            }
        }
        log_file << "Done with X2C Valence SFs" << endl;
    }
    delete (progress);
    log_file << "adding ECP contribution" << endl;
    auto sf2 = sf;
    auto sf2_def2 = sf_def2;
    add_ECP_contribution_test(
        asym_atom_list,
        wavy[0],
        sf2,
        k_pt);
    add_ECP_contribution_test(
        aal_def2,
        wavy[1],
        sf2_def2,
        k_pt);
    log_file << "done adding ECP contribution" << endl;
    log_file << "Calculating sfacs..." << endl;
    vec thakkar_sfac(k_pt[0].size());
    vec thakkar_core_sfac(k_pt[0].size());
    vec sf1_sfac(k_pt[0].size());         // Has the HAR density without ECP
    vec sf2_sfac(k_pt[0].size());         // Has the HAR density with ECP
    vec sf_def2_sfac(k_pt[0].size());     // Has the def2 density without ECP
    vec sf2_def2_sfac(k_pt[0].size());    // Has the def2 density with ECP
    vec sf_all_sfac(k_pt[0].size());      // Has the all electron density
    vec sf_all_val_sfac(k_pt[0].size());  // Has the valence electron density off the all electron wavefunction
    vec sf_ZORA_sfac(k_pt[0].size());     // Has the all electron density using ZORA
    vec sf_ZORA_val_sfac(k_pt[0].size()); // Has the valence electron density off the all electron wavefunction using ZORA
    vec sf_x2c_sfac(k_pt[0].size());      // Has the all electron density using X2C
    vec sf_x2c_val_sfac(k_pt[0].size());  // Has the valence electron density off the all electron wavefunction using X2C
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < k_pt[0].size(); i++)
    {
        sf1_sfac[i] = sqrt(pow(sf[0][i].real(), 2) + pow(sf[0][i].imag(), 2));
        sf2_sfac[i] = sqrt(pow(sf2[0][i].real(), 2) + pow(sf2[0][i].imag(), 2));
        sf_def2_sfac[i] = sqrt(pow(sf_def2[0][i].real(), 2) + pow(sf_def2[0][i].imag(), 2));
        sf2_def2_sfac[i] = sqrt(pow(sf2_def2[0][i].real(), 2) + pow(sf2_def2[0][i].imag(), 2));
        sf_all_sfac[i] = sqrt(pow(sf_all[0][i].real(), 2) + pow(sf_all[0][i].imag(), 2));
        sf_all_val_sfac[i] = sqrt(pow(sf_all_val[0][i].real(), 2) + pow(sf_all_val[0][i].imag(), 2));
        sf_ZORA_sfac[i] = sqrt(pow(sf_ZORA[0][i].real(), 2) + pow(sf_ZORA[0][i].imag(), 2));
        sf_ZORA_val_sfac[i] = sqrt(pow(sf_ZORA_val[0][i].real(), 2) + pow(sf_ZORA_val[0][i].imag(), 2));
        sf_x2c_sfac[i] = sqrt(pow(sf_x2c[0][i].real(), 2) + pow(sf_x2c[0][i].imag(), 2));
        sf_x2c_val_sfac[i] = sqrt(pow(sf_x2c_val[0][i].real(), 2) + pow(sf_x2c_val[0][i].imag(), 2));
        if (sf[0][i].real() < 0)
            sf1_sfac[i] = -sf1_sfac[i];
        if (sf_def2[0][i].real() < 0)
            sf_def2_sfac[i] = -sf_def2_sfac[i];
        if (sf_all_val[0][i].real() < 0)
            sf_all_val_sfac[i] = -sf_all_val_sfac[i];
        if (sf_ZORA_val[0][i].real() < 0)
            sf_ZORA_val_sfac[i] = -sf_ZORA_val_sfac[i];
        if (sf_x2c_val[0][i].real() < 0)
            sf_x2c_val_sfac[i] = -sf_x2c_val_sfac[i];
        thakkar_sfac[i] = Au.get_form_factor(k_pt[3][i]);
        thakkar_core_sfac[i] = Au.get_core_form_factor(k_pt[3][i], 60);
    }
    if (true)
    { // Change if you do not want X-ray
        ofstream result("sfacs.dat", ios::out);
        log_file << "Writing X-ray sfacs...";
        log_file.flush();
        // Now we just need to write the result to a file, together with the spherical results and separated for valence and core
        for (int i = 0; i < k_pt[0].size(); i++)
        {
            result << showpos << setw(8) << setprecision(5) << fixed << constants::ang2bohr(k_pt[3][i] / constants::FOUR_PI);
            result << showpos << setw(16) << setprecision(8) << scientific << thakkar_sfac[i];
            result << showpos << setw(16) << setprecision(8) << scientific << thakkar_core_sfac[i];
            result << showpos << setw(16) << setprecision(8) << scientific << sf1_sfac[i];
            result << showpos << setw(16) << setprecision(8) << scientific << sf2_sfac[i];
            result << showpos << setw(16) << setprecision(8) << scientific << sf_def2_sfac[i];
            result << showpos << setw(16) << setprecision(8) << scientific << sf2_def2_sfac[i];
            result << showpos << setw(16) << setprecision(8) << scientific << sf_all_sfac[i];
            result << showpos << setw(16) << setprecision(8) << scientific << sf_all_val_sfac[i];
            result << showpos << setw(16) << setprecision(8) << scientific << sf_ZORA_sfac[i];
            result << showpos << setw(16) << setprecision(8) << scientific << sf_ZORA_val_sfac[i];
            result << showpos << setw(16) << setprecision(8) << scientific << sf_x2c_sfac[i];
            result << showpos << setw(16) << setprecision(8) << scientific << sf_x2c_val_sfac[i];
            result << "\n";
        }
        log_file << " ... done!" << endl;
        result.flush();
        result.close();
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

void calc_cube_ML(vec data, WFN &dummy, const int &exp_coef, const int atom = -1)
{
    double MinMax[6]{0, 0, 0, 0, 0, 0};
    int steps[3]{0, 0, 0};
    readxyzMinMax_fromWFN(dummy, MinMax, steps, 2.5, 0.1);
    cube CubeRho(steps[0], steps[1], steps[2], dummy.get_ncen(), true);
    CubeRho.give_parent_wfn(dummy);

    for (int i = 0; i < 3; i++)
    {
        CubeRho.set_origin(i, MinMax[i]);
        CubeRho.set_vector(i, i, (MinMax[i + 3] - MinMax[i]) / steps[i]);
    }
    CubeRho.set_comment1("Calculated density using NoSpherA2 from ML Data");
    CubeRho.set_comment2("from " + dummy.get_path());
    CubeRho.path = get_basename_without_ending(dummy.get_path()) + "_rho.cube";

    time_point start = get_time();

    progress_bar *progress = new progress_bar{std::cout, 50u, "Calculating Values"};
    const int step = (int)std::max(floor(CubeRho.get_size(0) / 20.0), 1.0);
    if (atom != -1)
        std::cout << "Calculation for atom " << atom << std::endl;

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < CubeRho.get_size(0); i++)
    {
        for (int j = 0; j < CubeRho.get_size(1); j++)
            for (int k = 0; k < CubeRho.get_size(2); k++)
            {

                vec PosGrid{
                    i * CubeRho.get_vector(0, 0) + j * CubeRho.get_vector(0, 1) + k * CubeRho.get_vector(0, 2) + CubeRho.get_origin(0),
                    i * CubeRho.get_vector(1, 0) + j * CubeRho.get_vector(1, 1) + k * CubeRho.get_vector(1, 2) + CubeRho.get_origin(1),
                    i * CubeRho.get_vector(2, 0) + j * CubeRho.get_vector(2, 1) + k * CubeRho.get_vector(2, 2) + CubeRho.get_origin(2)};

                if (atom == -1)
                    CubeRho.set_value(i, j, k, calc_density_ML(PosGrid[0], PosGrid[1], PosGrid[2], data, dummy.atoms));
                else
                    CubeRho.set_value(i, j, k, calc_density_ML(PosGrid[0], PosGrid[1], PosGrid[2], data, dummy.atoms, atom));
            }
        if (i != 0 && i % step == 0)
            progress->write(i / static_cast<double>(CubeRho.get_size(0)));
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
    std::cout << "Number of electrons: " << std::fixed << std::setprecision(4) << CubeRho.sum() << std::endl;
    if (atom == -1)
    {
        CubeRho.write_file(true);
    }
    else
    {
        std::string fn(get_basename_without_ending(dummy.get_path()) + "_rho_" + std::to_string(atom) + ".cube");
        CubeRho.write_file(fn, false);
    }
};


const double dlm_function(const unsigned int& l, const int& m, const double& theta, const double& phi) {
    double result = (double)NAN;
    const double x = cos(theta);
    const double s = -sin(theta);
    const int m_ = abs(m);
    // the formula for a real spherical harmonic is:
    // for positive m:
    // (-1)^m cos(m * phi)
    // for negative m:
    // (-1)^m sin(-m * phi)
    const double p((m_) % 2 == 1 ? (m < 0 ? -sin(m_ * phi) : -cos(m_ * phi)) : (m < 0 ? sin(m_ * phi) : cos(m_ * phi)));
    switch (l) {
    case 0: switch (m_) {
    case 0: result = 1.; break;
    }; break;
    case 1: switch (m_) {
    case 0: result = x; break;
    case 1: result = p * s; break;
    }; break;
    case 2: switch (m_) {
    case 0: result = 0.5 * (3. * x * x - 1.); break;
    case 1: result = p * 3. * x * s; break;
    case 2: result = p * 3. * s * s; break;
    }; break;
    case 3: switch (m_) {
    case 0: result = 0.5 * x * (5. * x * x - 3.); break;
    case 1: result = p * 1.5 * (5. * x * x - 1.) * s; break;
    case 2: result = p * 15. * x * s * s; break;
    case 3: result = p * 15. * s * s * s; break;
    }; break;
    case 4: switch (m_) {
    case 0: result = 0.125 * (35. * x * x * x * x - 30. * x * x + 3.); break;
    case 1: result = p * 2.5 * (7. * x * x * x - 3. * x) * s; break;
    case 2: result = p * 7.5 * (7. * x * x - 1.) * s * s; break;
    case 3: result = p * 105. * x * s * s * s; break;
    case 4: result = p * 105. * s * s * s * s; break;
    }; break;
    case 5: switch (m_) {
    case 0: result = 0.125 * x * (63. * x * x * x * x - 70. * x * x + 15.); break;
    case 1: result = p * 1.875 * (21. * x * x * x * x - 14. * x * x + 1.) * s; break;
    case 2: result = p * 52.5 * x * (3. * x * x - 1.) * s * s; break;
    case 3: result = p * 52.5 * (9. * x * x - 1.) * s * s * s; break;
    case 4: result = p * 945. * x * s * s * s * s; break;
    case 5: result = p * 945. * s * s * s * s * s; break;
    }; break;
    case 6: switch (m_) {
    case 0: result = 0.0625 * (231. * x * x * x * x * x * x - 315. * x * x * x * x + 105. * x * x - 5.); break;
    case 1: result = p * 2.625 * x * (33. * x * x * x * x - 30. * x * x + 5.) * s; break;
    case 2: result = p * 13.125 * (33. * x * x * x * x - 18. * x * x + 1.) * s * s; break;
    case 3: result = p * 157.5 * x * (11. * x * x - 3.) * s * s * s; break;
    case 4: result = p * 472.5 * (11. * x * x - 1.) * s * s * s * s; break;
    case 5: result = p * 10395. * x * s * s * s * s * s; break;
    case 6: result = p * 10395. * s * s * s * s * s * s; break;
    }; break;
    case 7: switch (m_) {
    case 0: result = 0.0625 * x * (429. * x * x * x * x * x * x - 693. * x * x * x * x + 315. * x * x - 35.); break;
    case 1: result = p * 0.4375 * (429. * x * x * x * x * x * x - 495. * x * x * x * x + 135. * x * x - 5.) * s; break;
    case 2: result = p * 7.875 * x * (143. * x * x * x * x - 110. * x * x + 15.) * s * s; break;
    case 3: result = p * 39.375 * (143. * x * x * x * x - 66. * x * x + 3.) * s * s * s; break;
    case 4: result = p * 1732.5 * x * (13. * x * x - 3.) * s * s * s * s; break;
    case 5: result = p * 5197.5 * (13. * x * x - 1.) * s * s * s * s * s; break;
    case 6: result = p * 135135. * x * s * s * s * s * s * s; break;
    case 7: result = p * 135135. * s * s * s * s * s * s * s; break;
    }; break;
    case 8: switch (m_) {
    case 0: result = 0.0078125 * (6435. * x * x * x * x * x * x * x * x - 12012. * x * x * x * x * x * x + 6930. * x * x * x * x - 1260. * x * x + 35.); break;
    case 1: result = p * 0.5625 * x * (715. * x * x * x * x * x * x - 1001. * x * x * x * x + 385. * x * x - 35.) * s; break;
    case 2: result = p * 19.6875 * (143. * x * x * x * x * x * x - 143. * x * x * x * x + 33. * x * x - 1.) * s * s; break;
    case 3: result = p * 433.125 * x * (39. * x * x * x * x - 26. * x * x + 3.) * s * s * s; break;
    case 4: result = p * 1299.375 * (65. * x * x * x * x - 26. * x * x + 1.) * s * s * s * s; break;
    case 5: result = p * 67567.5 * x * (5. * x * x - 1.) * s * s * s * s * s; break;
    case 6: result = p * 67567.5 * (15. * x * x - 1.) * s * s * s * s * s * s; break;
    case 7: result = p * 2027025. * x * s * s * s * s * s * s * s; break;
    case 8: result = p * 2027025. * s * s * s * s * s * s * s * s; break;
    }; break;
    };
    return result;
}


//cdouble S_n_recursion(int n, double H, double b);
//cdouble C_n_recursion(int n, double H, double b);
//For the case J_l(H) = int_0^inf j_l(Hr) * R_l(r)^2 * r^2 dr    |   Wave Functions!!

//cdouble S_0(double H, double b) {
//	using namespace std::complex_literals;
//    double two_32 = pow(2.0, 1.5);
//	return -(cerf((1.0i * H) / (two_32 * sqrt(b))) * constants::sqr_pi * 1.0i * exp(-H * H / (8. * b))) / (two_32 * sqrt(b));
//}
//double C_0(double H, double b) {
//	return  (constants::sqr_pi * exp(-H * H / (8. * b))) / (pow(2., 1.5) * sqrt(b));
//}
////Following 1/(4b) * ((n-1)C_(n-2) - H*S_(n-1)) = C_n
//cdouble C_n_recursion(int n, double H, double b) {
//    using namespace std::complex_literals;
//    if (n == 0) {
//        return C_0(H, b);
//    }
//    else if (n == 1) {
//        return (1. / (4. * b)) * (1. + S_0(H, b));
//    }
//    else {
//        return (1. / (4. * b)) * ((n - 1.) * C_n_recursion(n - 2, H, b) - H * S_n_recursion(n - 1, H, b));
//    }
//}
////\int_{ 0 }^ {\infty} r^ n\sin(Hr) \cdot e^ { -2br ^ 2 } dr & = \frac{ 1 }{4b}\left((n - 1)S_{ n - 2 } + HC_{ n - 1 }\right) = S_n
//cdouble S_n_recursion(int n, double H, double b) {
//    using namespace std::complex_literals;
//    if (n == 0) {
//        return S_0(H, b);
//    }
//    else if (n == 1) {
//        return (1. / (4. * b)) * H * C_0(H, b);
//    }
//    else {
//        return (1. / (4. * b)) * ((n - 1.) * S_n_recursion(n - 2, H, b) + H * C_n_recursion(n - 1, H, b));
//    }
//}
//// This function yields the fourier bessel transform of the radial integral of a gaussian density function (compare equation 1.2.7.9 in 10.1107/97809553602060000759),a ssuming that H = 2 \pi S
//cdouble fourier_bessel_integral(
//    primitive& p,
//    const double H
//)
//{
//    using namespace std::complex_literals;
//    const int l = p.type;
//    const double b = p.exp;
//    double N;
//    //double N = pow(
//    //    pow(2, 7 + 4 * l) * pow(b, 3 + 2 * l) / constants::PI / pow(doublefactorial(2 * l + 1), 2),
//    //    0.25);
//    //pow(8 * pow(b, 3) / pow(constants::PI, 3), 0.25);
//    //double N = p.norm_const;
//    //return N * (pow(H, l * 2) * constants::sqr_pi * exp(-H * H / (4 * b))) / (pow(2, l + 2) * pow(b, l + 1.5));
//    if (l == 0)
//    {
//        N = 1.;
//        //return N * N * (pow(H, l * 2) * constants::sqr_pi * exp(-H * H / (8 * b))) / (pow(2, l + 3.5) * pow(b, l + 1.5));
//        return  ((N * N) / (4 * b)) * C_n_recursion(0, H, b);
//    }
//    else if (l == 1)
//    {
//        N = 1.; //p.norm_const;
//        return ((N * N) / (H * H)) * (S_n_recursion(2, H, b) - H * C_n_recursion(3, H, b));
//    }
//}

//For the case J_l(H) = int_0^inf j_l(Hr) * R_l(r) * r^2 dr    |   Densities!!


// This function yields the fourier bessel transform of the radial integral of a gaussian density function (compare equation 1.2.7.9 in 10.1107/97809553602060000759),a ssuming that H = 2 \pi S
double fourier_bessel_integral(
    primitive& p,
    const double H
)
{
    using namespace std::complex_literals;
    const int l = p.type;
    const double b = p.exp;
    //double N = p.norm_const;
    double N = pow(
        pow(2, 7 + 4 * l) * pow(b, 3 + 2 * l) / constants::PI / pow(doublefactorial(2 * l + 1), 2),
        0.25);

    return N * (pow(H, l) * constants::sqr_pi * exp(-H * H / (4 * b))) / (pow(2, l + 2) * pow(b, l + 1.5));
}

//Returns the spherical coordinates of a given cartesian vector
//the output is a vector with the (radius, theta and phi)
vec cartesian_to_spherical(double x, double y, double z) {
    double theta, phi;
    double r = sqrt(x * x + y * y + z * z);
    if (r == 0) {
        theta = 0;
        phi = 0;
    }
    else {
        theta = acos(z / r);
        phi = atan2(y, x);
    }
    return { r, theta, phi };
}

//https://en.wikipedia.org/wiki/Associated_Legendre_polynomials
double legendre_polynomial(int l, int m, double x) {
    if (l == 0) {
        return 1.;
    }
    else if (l == 1) {
        if (m == 0) {
            return x;
        }
        else if (m == 1) {
            return sqrt(1 - x * x);
        }
        else if (m == -1) {
            return -0.5 * sqrt(1. - x * x);
        }
        else err_checkf(false, "This is impossible!", std::cout);
    }
    else if (l == 2) {
        if (m == 0) {
            return 0.5 * (3 * x * x - 1.);
        }
        else if (m == 1) {
            return -3. * x * sqrt(1. - x * x);
        }
        else if (m == 2) {
            return 3. * (1. - x * x);
        }
        else if (m == -1) {
            return 0.5 * x * sqrt(1. - x * x);
        }
        else if (m == -2) {
            return 0.125 * (1. - x * x);
        }
        else err_checkf(false, "This is impossible!", std::cout);
    }
    else if (l == 3) {
        if (m == 0) {
            return 0.5 * x * (5. * x * x - 3.);
        }
        else if (m == 1) {
            return 1.5 * (1 - 5 * x * x) * sqrt((1 - x * x));
        }
        else if (m == 2) {
            return 15. * x * (1 - x * x);
        }
        else if (m == 3) {
            return -15. * sqrt(1 - x * x) * (1 - x * x);
        }
        else if (m == -1) {
            return -0.125 * (1 - 5 * x * x) * sqrt((1 - x * x));
        }
        else if (m == -2) {
            return 0.125 * x * (1 - x * x);
        }
        else if (m == -3) {
            return (15. / 720.) * sqrt(1 - x * x) * (1 - x * x);
        }
        else err_checkf(false, "This is impossible!", std::cout);
    }
    else err_not_impl_f("This Legendre Polynomial doesn't yet exist in the code!", std::cout);
}

//Original implementation after P. Coppens DOI: 10.1107/97809553602060000759 Eq. 1.2.7.2b
//I omitted the abs(m) in the factorial as most other sources do not include it
double real_spherical(int l, int m, double theta, double phi) {
    if (m > 0) {
        return sqrt(((2 * l + 1) / constants::TWO_PI) * double(constants::ft[l - (m)]) / double(constants::ft[l + (m)])) * legendre_polynomial(l, m, cos(theta)) * cos(m * phi);
    }
    else if (m < 0) {
        return sqrt(((2 * l + 1) / constants::TWO_PI) * double(constants::ft[l - (m)]) / double(constants::ft[l + (m)])) * legendre_polynomial(l, m, cos(theta)) * sin(m * phi);
    }
    else
        return sqrt((2 * l + 1) / constants::FOUR_PI) * legendre_polynomial(l, m, cos(theta));
}

cdouble sfac_bessel(
    primitive& p,
    vec& k_point,
    const vec& coefs
)
{
    using namespace std::complex_literals;
    vec local_k = k_point;
    double leng = sqrt(local_k[0] * local_k[0] + local_k[1] * local_k[1] + local_k[2] * local_k[2]);
    double H = 2 * constants::PI * leng;
    //normalize the spherical harmonics k_point
    for (int i = 0; i < 3; i++)
        local_k[i] /= leng;

    vec spherical = cartesian_to_spherical(local_k[0], local_k[1], local_k[2]);


    // IMPLEMENTATION FOR DENSITY AS INPUT
    double radial = fourier_bessel_integral(p, H) * p.coefficient;
    cdouble result(0.0, 0.0);
    for (int m = -p.type; m <= p.type; m++) {
        cdouble angular = real_spherical(p.type, m, spherical[1], spherical[2]);
        result += constants::FOUR_PI * pow(1.0i, p.type) * radial * angular * coefs[m + p.type];
    }
    return result;

    //// IMPLEMENTATION FOR WAVE FUNCTIONS AS INPUT
    //  cdouble radial = fourier_bessel_integral(p, H);
    //  if (p.type == 0)
    //  {
    //      //constants::FOUR_PI^2 weil in angular^2 ein Faktor 1/4pi drin ist
    //      return  constants::FOUR_PI * constants::FOUR_PI * pow(1.0i, p.type) * radial * p.coefficient * p.coefficient * abs(angular * angular);
    //  }
    //  else if (p.type == 1) {
          ////double angular = dlm_function(p.type, m+1, spherical[1], spherical[2]);
    //      double angular2 = constants::spherical_harmonic(p.type, m, local_k.data());
    //      return pow(0.75*constants::PI, 2) * constants::FOUR_PI * pow(1.0i, p.type) * radial * p.coefficient * p.coefficient * abs(angular * angular);
    //  }
}

//Only valid for one atom positioned at 0,0,0
double compute_MO_spherical_orig(double x, double y, double z, double expon, double coef, int type) {
    double result = 0.0;
    int iat = 0;
    int l[3]{ 0, 0, 0 };
    double ex = 0;
    double temp = 0;

    // x, y, z and dsqd
    vec d{
        x,
        y,
        z,
        x * x + y * y + z * z
    };


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

//Only valid for one atom positioned at 0,0,0
double compute_MO_spherical_new(double x, double y, double z, double expon, double coef, int type, int m = 1) {
    double result = 0.0;
    int iat = 0;
    int l[3]{ 0, 0, 0 };
    double ex = 0;
    double temp = 0;
    int type_local = type;
    double N = 1;
    //convert "type" to l
    if (type == 1) {
        type_local = 0;
        N =1/ sqrt(1 / (constants::FOUR_PI));
    }
    else if (type > 1 && type < 5) {
        type_local = 1;
        N = 1/ sqrt(3 / (constants::FOUR_PI));
    }
    
    // r, theta, phi
    vec spher = cartesian_to_spherical(x, y, z);
    double R = exp(-expon * spher[0] * spher[0]) * pow(spher[0], type_local);
    double Y = real_spherical(type_local, m, spher[1],spher[2]) * N;
    result = coef * R * Y;
    double test = compute_MO_spherical_orig(x, y, z, expon, coef, type);
    std::cout << "New: " << result << " Old: " << test << " diff: " << result - test << " fact: " << result / test << std::endl;
    return result;
}


void Calc_MO_spherical_harmonics(
    cube& CubeMO,
    WFN& wavy,
    int cpus,
    std::ostream& file,
    bool nodate)
{
    using namespace std;
#ifdef _OPENMP
    if (cpus != -1)
    {
        if (cpus > 1)
            omp_set_nested(1);
    }
#endif

    time_point start = get_time();

    progress_bar* progress = NULL;
    if (!nodate)
        progress = new progress_bar{ file, 50u, "Calculating Values" };
    const int step = (int)max(floor(CubeMO.get_size(0) / 20.0), 1.0);

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

                CubeMO.set_value(i, j, k, compute_MO_spherical_new(PosGrid[0], PosGrid[1], PosGrid[2], wavy.get_exponent(0), wavy.get_MO_coef(0,0), wavy.get_type(0)));
            }
        if (i != 0 && i % step == 0 && !nodate)
            progress->write((i) / static_cast<double>(CubeMO.get_size(0)));
    }
    if (!nodate) {
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

void test_analytical_fourier() {
    WFN wavy(0);
    wavy.push_back_MO(0, 1.0, -13);
    wavy.push_back_atom("H", 0, 0, 0, 1);
    const double c_exp = 2.0;
    double vals[] = { 1.0 };
    vec coefs = { 0.0, 0.0, 0.0 , 0.0 , 1.0};
    int type = 2;
    //wavy.add_primitive(1, type + 1, c_exp, vals);
    wavy.atoms[0].push_back_basis_set(c_exp, vals[0],  type, 0);
    //double MinMax[6]{ 0, 0, 0, 0, 0, 0 };
    //int steps[3]{ 0, 0, 0 };
    //readxyzMinMax_fromWFN(wavy, MinMax, steps, 3., 0.025, true);
    //cube CubeMO(steps[0], steps[1], steps[2], 1, true);
    //CubeMO.give_parent_wfn(wavy);
    //for (int i = 0; i < 3; i++)
    //{
    //    CubeMO.set_origin(i, MinMax[i]);
    //    CubeMO.set_vector(i, i, (MinMax[i + 3] - MinMax[i]) / steps[i]);
    //}
    //Calc_MO_spherical_harmonics(CubeMO, wavy, 18, std::cout, false);
    //CubeMO.path = "MO.cube";
    //CubeMO.write_file(true);

    vec2 kpts;
    for (int i = 1; i < 100; i++) {
        kpts.push_back({ 0.01 * i, 0, 0 });
        kpts.push_back({ 0, 0.01 * i, 0 });
        kpts.push_back({ 0,0,0.01 * i });
        kpts.push_back({ -0.01 * i, 0, 0 });
        kpts.push_back({ 0, -0.01 * i, 0 });
        kpts.push_back({ 0, 0, -0.01 * i });
    }
    cvec2 sf_A, sf_N;
    sf_A.resize(1);
    sf_N.resize(1);
    sf_A[0].resize(kpts.size(), 0.0);
    sf_N[0].resize(kpts.size(), 0.0);
    vec2 grid;
    grid.resize(5); //x, y, z, dens, atomic_weight
    
    double alpha_min[] = { 0.5 };
    AtomGrid griddy(1E-40,
        1120,
        1350,
        1,
        3.6,
        1,
        alpha_min,
        std::cout);

    double pos[] = { 0 };
    for (int i = 0; i < grid.size(); i++) {
        grid[i].resize(griddy.get_num_grid_points(), 0.0);
    }
    griddy.get_atomic_grid(0, pos,pos, pos, grid[0].data(), grid[1].data(), grid[2].data(), grid[4].data());

    for (int i = 0; i < grid[0].size(); i++) {
        //grid[3][i] = wavy.compute_dens(grid[0][i], grid[1][i], grid[2][i]);
        grid[3][i] = calc_density_ML(grid[0][i], grid[1][i], grid[2][i], coefs, wavy.atoms);
    }

    primitive p(1, type, c_exp, vals[0]);
#pragma omp parallel for
    for (int i = 0; i < kpts.size(); i++) {
        sf_A[0][i] = sfac_bessel(p, kpts[i], coefs);
        for (int p = 0; p < grid[0].size(); p++) {
            double work = 2 * constants::PI * (kpts[i][0] * grid[0][p] + kpts[i][1] * grid[1][p] + kpts[i][2] * grid[2][p]);
            sf_N[0][i] += std::polar(grid[3][p] * grid[4][p],work);
        }
	}
    using namespace std;
	ofstream result("sfacs.dat", ios::out);
	for (int i = 0; i < kpts.size(); i++) {
        vec local_k = kpts[i];
        double leng = sqrt(local_k[0] * local_k[0] + local_k[1] * local_k[1] + local_k[2] * local_k[2]);
        //normalize the spherical harmonics k_point
        for (int i = 0; i < 3; i++)
            local_k[i] /= leng;

        //double angular = dlm_function(p.type, m, acos(k_point[2]), atan2(k_point[1], k_point[0]));
        //double angular_m1 = constants::spherical_harmonic(p.type, -1, local_k.data());
        //double angular_1 = constants::spherical_harmonic(p.type, 1, local_k.data());
        //double angular_0 = constants::spherical_harmonic(p.type, 0, local_k.data());
		result << setw(8) << setprecision(2) << fixed << kpts[i][0];
        result << setw(8) << setprecision(2) << fixed << kpts[i][1];
        result << setw(8) << setprecision(2) << fixed << kpts[i][2];
		result << setw(16) << setprecision(8) << scientific << sf_A[0][i].real();
		result << setw(16) << setprecision(8) << scientific << sf_A[0][i].imag();
		result << setw(16) << setprecision(8) << scientific << sf_N[0][i].real();
		result << setw(16) << setprecision(8) << scientific << sf_N[0][i].imag();
        result << setw(16) << setprecision(8) << scientific << abs(sf_A[0][i] - sf_N[0][i]);
        result << setw(35) << setprecision(8) << scientific << sf_A[0][i] / sf_N[0][i];
        //result << setw(16) << setprecision(8) << scientific << abs(angular_m1 * angular_m1);
        //result << setw(16) << setprecision(8) << scientific << abs(angular_0 * angular_0);
        //result << setw(16) << setprecision(8) << scientific << abs(angular_1 * angular_1);
		result << "\n";
	}
	result.flush();
	result.close();
        
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
    progress_bar *progress = new progress_bar{std::cout, 50u, "Calculating Values"};
    const int step = (int)std::max(floor(total_size / 20.0), 1.0);

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
            if (index != 0 && index % step == 0)
                progress->write(index / static_cast<double>(total_size));
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

double numerical_3d_integral(std::function<double(double *)> f, double stepsize = 0.001, double r_max = 10.0)
{
    double tot_int = 0;
    double da = stepsize * constants::PI;
    double dr = stepsize;
    int upper = static_cast<int>(constants::TWO_PI / da);
    long long int total_calcs = upper * (long long int)(constants::PI / da * r_max / dr);
    std::cout << "Integrating over " << total_calcs << " points" << std::endl;

    progress_bar *progress = new progress_bar{std::cout, 50u, "Integrating"};
    const long long int step = (long long int)std::max(floor(upper / 20.0), 1.0);

#pragma omp parallel for reduction(+ : tot_int) schedule(dynamic)
    for (int phic = 0; phic <= upper; phic++)
    {
        double phi = phic * da;
        double cp = std::cos(phi);
        double sp = std::sin(phi);
        for (double theta = 0; theta <= constants::PI; theta += da)
        {
            double st = std::sin(theta);
            double ct = std::cos(theta);
            double x0 = st * cp, y0 = st * sp;
            double ang_mom = dr * da * da * st;
            for (double r = r_max; r >= 0; r -= dr)
            {
                double d[4]{r * x0,
                            r * y0,
                            r * ct,
                            r};
                tot_int += f(d) * r * r * ang_mom;
            }
        }
        if (phic != 0 && phic % step == 0)
            progress->write(static_cast<double>(phic) / static_cast<double>(upper));
    }

    std::cout << "\ntotal integral : " << std::fixed << std::setprecision(4) << tot_int << std::endl;
    return tot_int;
}

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
    calc_cube_ML(data, dummy, nr_coefs);

    for (int i = 0; i < dummy.get_ncen(); i++)
        calc_cube_ML(data, dummy, nr_coefs, i);
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

void test_core_dens_corrected(double& precisison, int ncpus = 4, std::string ele = "Au", ivec val_els_alpha = {}, ivec val_els_beta = {})
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

    progress_bar *progress = new progress_bar{cout, 60u, "Calculating densities"};
    const int upper = static_cast<int>(res[0].size());
    const long long int step = max(static_cast<long long int>(floor(upper / 20)), 1LL);
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
        if (i != 0 && i % step == 0)
            progress->write(i / static_cast<double>(upper));
    }
    delete (progress);
    time_point end = get_time();
    cout << "Time taken: " << round(get_sec(start, end)/60) << " m " << get_sec(start, end) % 60 << " s " << get_msec(start, end) << " ms" << endl;
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
    ECP_way_Au.read_gbw(ele+"_def2TZVP.gbw", cout, false, false);
    ECP_way_Au.delete_unoccupied_MOs();
    ECP_way_Au.set_has_ECPs(true, true, 1);

    WFN wavy_full_Au(9);
    wavy_full_Au.read_gbw(ele+"_jorge.gbw", cout, false);
    wavy_full_Au.delete_unoccupied_MOs();
    WFN wavy_val_Au(9);
    wavy_val_Au.read_gbw(ele+"_jorge.gbw", cout, false);
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

    progress_bar *progress = new progress_bar{cout, 60u, "Calculating SFACs"};
    const int upper = static_cast<int>(res[0].size());
    const long long int step = max(static_cast<long long int>(floor(upper / 20)), 1LL);
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
        if (i != 0 && i % step == 0)
            progress->write(i / static_cast<double>(upper));
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

void test_openblas() {
#if has_RAS
    std::cout << "Running Openblas test" << std::endl;
    _test_openblas();
    exit(0);
#else
    err_not_impl_f("Openblas not included in Executable!", std::cout);
#endif
}
