#pragma once
#include "spherical_density.h"
#include "scattering_factors.h"
#include "convenience.h"
#include "npy.h"
#include "properties.h"

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
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single<int>(1, temp, i).real();
            temp = C.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single<int>(6, temp, i).real();
            temp = O.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single<int>(8, temp, i).real();
            temp = P.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single<int>(15, temp, i).real();
            temp = Ca.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single<int>(20, temp, i).real();
            temp = Os.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single<int>(76, temp, i).real();
            temp = 0.0;
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single<int>(1, temp, i).real();
            temp = C_cat.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single<int>(6, temp, i).real();
            temp = O_cat.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single<int>(8, temp, i).real();
            temp = P_cat.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single<int>(15, temp, i).real();
            temp = Ca_cat.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single<int>(20, temp, i).real();
            temp = H_an.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single<int>(1, temp, i).real();
            temp = C_an.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single<int>(6, temp, i).real();
            temp = O_an.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single<int>(8, temp, i).real();
            temp = P_an.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single<int>(15, temp, i).real();
            temp = C.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single<int>(6, temp, i, 1).real();
            temp = O.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single<int>(8, temp, i, 1).real();
            temp = C.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single<int>(6, temp, i, -1).real();
            temp = O.get_form_factor(k_value);
            result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single<int>(7, temp, i, -1).real();
            result << endl;
        }
        result.flush();
        result.close();
    }
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

    time_point start = get_time();
    time_point end_becke, end_prototypes, end_spherical, end_prune, end_aspherical;

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
                         start,
                         end_becke,
                         end_prototypes,
                         end_spherical,
                         end_prune,
                         end_aspherical,
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

void test_timing()
{
    using namespace std;
    cout << NoSpherA2_message() << endl;
    options opt;
    opt.dmin = 0.5;
    opt.combined_tsc_calc_files.push_back("olex2\\Wfn_job\\Part_1\\thpp.wfx");
    opt.combined_tsc_calc_files.push_back("olex2\\Wfn_job\\Part_2\\thpp.wfx");
    opt.accuracy = 2;
    ivec t1 = {0, 1};
    ivec t2 = {0, 2};
    opt.groups.pop_back();
    opt.groups.push_back(t1);
    opt.groups.push_back(t2);
    opt.cif = "thpp.cif";
    opt.combined_tsc_calc = true;
    opt.binary_tsc = true;

    int max_threads = 48;
    vector<long long int> time(max_threads);

    for (int ncpus = 1; ncpus < 49; ncpus++)
    {
        auto start = std::chrono::high_resolution_clock::now();
        opt.threads = ncpus;
#ifdef _OPENMP
        omp_set_num_threads(ncpus);
        // print the number of threads used
        std::cout << "Using " << opt.threads << " threads" << std::endl;
#pragma omp parallel
        {
#pragma omp single
            cout << omp_get_num_threads() << endl;
        }
#else
        std::cout << "You are running an OmpenMP Benchamrk without parallelization" << std::endl;
        std::cout << "This will not give you any useful information" << std::endl;
#endif
        std::vector<WFN> wavy;
        err_checkf(opt.hkl != "" || opt.dmin != 99.0, "No hkl specified and no dmin value given", std::cout);
        if (opt.combined_tsc_calc)
            err_checkf(opt.cif != "", "No cif specified", std::cout);
        // First make sure all files exist
        for (int i = 0; i < opt.combined_tsc_calc_files.size(); i++)
        {
            err_checkf(exists(opt.combined_tsc_calc_files[i]), "Specified file for combined calculation doesn't exist! " + opt.combined_tsc_calc_files[i], std::cout);
            if (opt.cif_based_combined_tsc_calc)
                err_checkf(exists(opt.combined_tsc_calc_cifs[i]), "Specified file for combined calculation doesn't exist! " + opt.combined_tsc_calc_cifs[i], std::cout);
        }
        wavy.resize(opt.combined_tsc_calc_files.size());
        for (int i = 0; i < opt.combined_tsc_calc_files.size(); i++)
        {
            std::cout << "Reading: " << setw(44) << opt.combined_tsc_calc_files[i] << flush;
            wavy[i].read_known_wavefunction_format(opt.combined_tsc_calc_files[i], std::cout);
            std::cout << " done!\nNumber of atoms in Wavefunction file: " << wavy[i].get_ncen() << " Number of MOs: " << wavy[i].get_nmo() << endl;
        }

        svec known_scatterer;
        vec2 known_kpts;
        tsc_block<int, cdouble> result;
        for (int i = 0; i < opt.combined_tsc_calc_files.size(); i++)
        {
            known_scatterer = result.get_scatterers();
            if (wavy[i].get_origin() != 7)
            {
                result.append(calculate_scattering_factors_MTC(
                                  opt,
                                  wavy,
                                  std::cout,
                                  known_scatterer,
                                  i,
                                  &known_kpts),
                              std::cout);
            }
            else
            {
                result.append(MTC_thakkar_sfac(
                                  opt,
                                  std::cout,
                                  known_scatterer,
                                  wavy,
                                  i),
                              std::cout);
            }
        }

        known_scatterer = result.get_scatterers();
        std::cout << "Final number of atoms in .tsc file: " << known_scatterer.size() << endl;
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        time[ncpus - 1] = duration.count();
    }

    // Print out results of timing
    for (int i = 0; i < max_threads; i++)
    {
        std::cout << "Time for " << i + 1 << " threads: " << time[i] << " microseconds" << std::endl;
    }

    exit(0);
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

void spherically_averaged_density(options &opt, const ivec val_els_alpha, const ivec val_els_beta)
{
    using namespace std;
    WFN wavy(9);
    wavy.read_known_wavefunction_format(opt.wfn, cout, opt.debug);
    wavy.delete_unoccupied_MOs();
    cout << "Number of MOs before: " << wavy.get_nmo() << endl;
    bvec MOs_to_delete(wavy.get_nmo(), false);
    int deleted = 0;
    if (val_els_alpha.size() > 0)
    {
        // Delete core orbitals
        int offset = wavy.get_MO_op_count(0);
        for (int i = offset - 1; i >= 0; i--)
            // only delete if i is not an element of core-size
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
    const long double da = constants::PI / 360.0;
    const int upper = static_cast<int>(constants::TWO_PI / da);
    // Make angular grids
    vector<long double> x0, y0, z0, st0;
    for (long double phi = 0; phi <= constants::TWO_PI; phi += da)
    {
        const long double cp = cos(phi);
        const long double sp = sin(phi);
        for (long double theta = 0; theta <= constants::PI; theta += da)
        {
            const long double st = sin(theta);
            st0.push_back(st * da * da);
            const long double ct = cos(theta);
            x0.push_back(st * cp), y0.push_back(st * sp), z0.push_back(ct);
        }
    }
    const long long int ang_size = x0.size();
    cout << "Averaging over " << ang_size << " angular points" << endl;

    // Make radial grids on logarithmic scale
    vec radial_dist = geomspace<double>(1E-7, 15, 20000);
    const long long int upper_r = radial_dist.size();
    const long long int total_points = static_cast<long long int>(upper_r * ang_size);
    cout << "Calculating " << total_points << " points" << endl;
    // Calcualte density on angular grid at each point of radial grid, average and integrate
    long double tot_int = 0;
    vec radial_dens(upper_r, 0.0);

#pragma omp parallel num_threads(opt.threads)
    {
        vec2 d;
        vec phi(wavy.get_nmo(), 0.0);
        d.resize(16);
        for (int i = 0; i < 16; i++)
            d[i].resize(1, 0.0);
#pragma omp for reduction(+ : tot_int)
        for (long long int _r = 0; _r < upper_r; _r++)
        {
            double r = radial_dist[_r];
            long double v = 0;
            for (long long int i = 0; i < ang_size; i++)
            {
                v += wavy.compute_dens(r * x0[i], r * y0[i], r * z0[i], d, phi) * st0[i];
            }
            if (_r >= 1)
                tot_int += v * r * r * (r - radial_dist[_r - 1]);
            else
                tot_int += v * r * r * (r);
            radial_dens[_r] = v / 4.0 / constants::PI;
        }
    }
    cout << "Start writing the file" << endl;
    string el = atnr2letter(wavy.get_atom_charge(0));
    ofstream out(el + ".dat", ios::out);
    out << "Total Integral: " << setw(18) << scientific << setprecision(10) << tot_int << "\n";
    for (int i = 0; i < upper_r; i++)
        out << setw(24) << scientific << setprecision(15) << radial_dist[i] << setw(24) << scientific << setprecision(16) << radial_dens[i] << "\n";
    out.flush();
    out.close();
}

void add_ECP_contribution_test(const ivec &asym_atom_list,
                               const WFN &wave,
                               cvec2 &sf,
                               vec2 &k_pt,
                               std::ostream &file,
                               const bool debug)
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

    time_point start = get_time(), end_becke, end_prototypes, end_spherical, end_prune, end_aspherical;

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
                         start,
                         end_becke,
                         end_prototypes,
                         end_spherical,
                         end_prune,
                         end_aspherical,
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
                         start,
                         end_becke,
                         end_prototypes,
                         end_spherical,
                         end_prune,
                         end_aspherical,
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
                         start,
                         end_becke,
                         end_prototypes,
                         end_spherical,
                         end_prune,
                         end_aspherical,
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
                         start,
                         end_becke,
                         end_prototypes,
                         end_spherical,
                         end_prune,
                         end_aspherical,
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
                         start,
                         end_becke,
                         end_prototypes,
                         end_spherical,
                         end_prune,
                         end_aspherical,
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
                         start,
                         end_becke,
                         end_prototypes,
                         end_spherical,
                         end_prune,
                         end_aspherical,
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
                         start,
                         end_becke,
                         end_prototypes,
                         end_spherical,
                         end_prune,
                         end_aspherical,
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
                         start,
                         end_becke,
                         end_prototypes,
                         end_spherical,
                         end_prune,
                         end_aspherical,
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
        k_pt,
        log_file,
        opt.debug);
    add_ECP_contribution_test(
        aal_def2,
        wavy[1],
        sf2_def2,
        k_pt,
        log_file,
        opt.debug);
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
            std::cout << spherical_harmonic(lam, m, d.data()) << " ";
        }
        std::cout << "\n";
    }
};

void calc_cube_ML(vec data, WFN &dummy, int &exp_coef, int atom = -1)
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
                    CubeRho.set_value(i, j, k, calc_density_ML(PosGrid[0], PosGrid[1], PosGrid[2], data, dummy.atoms, exp_coef));
                else
                    CubeRho.set_value(i, j, k, calc_density_ML(PosGrid[0], PosGrid[1], PosGrid[2], data, dummy.atoms, exp_coef, atom));
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

void calc_rho_cube(WFN& dummy)
{
    using namespace std;
    double MinMax[6]{ 0, 0, 0, 0, 0, 0 };
    int steps[3]{ 0, 0, 0 };
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
    const int s1 = CubeRho.get_size(0), s2 = CubeRho.get_size(1), s3 = CubeRho.get_size(2), total_size = s1 * s2 * s3;;
    cout << "Lets go into the loop! There is " << total_size << " points" << endl;
    progress_bar* progress = new progress_bar{ std::cout, 50u, "Calculating Values" };
    const int step = (int)std::max(floor(total_size / 20.0), 1.0);
    
    vec v1{
    CubeRho.get_vector(0, 0),
    CubeRho.get_vector(1, 0),
    CubeRho.get_vector(2, 0)
    }, v2{
    CubeRho.get_vector(0, 1),
    CubeRho.get_vector(1, 1),
    CubeRho.get_vector(2, 1)
    }, v3{
    CubeRho.get_vector(0, 2),
    CubeRho.get_vector(1, 2),
    CubeRho.get_vector(2, 2)
    }, orig{
    CubeRho.get_origin(0),
    CubeRho.get_origin(1),
    CubeRho.get_origin(2)
    };


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
                i * v1[2] + j * v2[2] + k * v3[2] + orig[2] };


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
    return pow(gaussian_radial(p_0, d[3]) * spherical_harmonic(p_0.type, m, d), 2);
}

double one(double *d)
{
    return 1;
}

void cube_from_coef_npy(std::string &coef_fn, std::string &xyzfile)
{
    std::vector<unsigned long> shape{};
    bool fortran_order;
    vec data{};

    npy::LoadArrayFromNumpy(coef_fn, shape, fortran_order, data);
    WFN dummy(7);
    dummy.read_xyz(xyzfile, std::cout);

    int nr_coefs = load_basis_into_WFN(dummy, TZVP_JKfit);
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

void test_core_dens()
{
    using namespace std;
    ofstream out("core_dens.log", ios::out);
    vec2 res(12);
    for (int i = 0; i < 12; i++)
        res[i].resize(1000000, 0.0);

    ofstream dat_out("core_dens_Li.dat", ios::out);

    Thakkar T_Li(3);

    WFN ECP_way_LI(9);
    ECP_way_LI.read_gbw("Li_ECP.gbw", out, true, true);

    WFN wavy3_LI(9);
    wavy3_LI.read_gbw("Li_QZVP.gbw", out, false);
    wavy3_LI.delete_unoccupied_MOs();
    WFN wavy4_LI(9);
    wavy4_LI.read_gbw("Li_QZVP.gbw", cout, false);
    wavy4_LI.delete_unoccupied_MOs();
    wavy4_LI.delete_MO(1); // delete the 2nd MO
#pragma omp parallel for
    for (int i = 0; i < res[0].size(); i++)
    {
        double sr = i * 0.00001;
        res[0][i] = sr;
        res[1][i] = T_Li.get_core_density(sr, 2);
        res[2][i] = T_Li.get_radial_density(sr);
        res[3][i] = ECP_way_LI.compute_dens(sr, 0, 0);
        res[4][i] = wavy3_LI.compute_dens(sr, 0, 0);
        res[5][i] = wavy4_LI.compute_dens(sr, 0, 0);
    }
    dat_out << fixed;
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
    exit(0);

    Thakkar T_Rb(37);
    string base = "def2-ECP";
    Gaussian_Atom G_Rb(37, base);
    WFN ECP_way(9);
    ECP_way.read_gbw("Atomic_densities\\atom_atom37.gbw", out, true, true);

    WFN multi_ecp_wavy(9);
    multi_ecp_wavy.read_gbw("Rb4.gbw", out, true, true);

    WFN wavy(9);
    wavy.read_gbw("Rb.gbw", out, false);
    wavy.delete_unoccupied_MOs();
    WFN wavy2(9);
    wavy2.read_gbw("Rb.gbw", cout, false);
    wavy2.delete_unoccupied_MOs();
    wavy2.pop_back_MO(); // delete the last 5 MOs: 5s, 3* 4p & 4s
    wavy2.pop_back_MO();
    wavy2.pop_back_MO();
    wavy2.pop_back_MO();
    wavy2.pop_back_MO();

    WFN wavy3(9);
    wavy3.read_gbw("Rb_QZVP.gbw", out, false);
    wavy3.delete_unoccupied_MOs();
    WFN wavy4(9);
    wavy4.read_gbw("Rb_QZVP.gbw", cout, false);
    wavy4.delete_unoccupied_MOs();
    wavy4.pop_back_MO(); // delete the last 5 MOs: 5s, 3* 4p & 4s
    wavy4.pop_back_MO();
    wavy4.pop_back_MO();
    wavy4.pop_back_MO();
    wavy4.pop_back_MO();

#pragma omp parallel for
    for (int i = 0; i < res[0].size(); i++)
    {
        double r, sr;
        r = i * 0.001;
        sr = r * 0.01;
        res[0][i] = r;
        res[1][i] = T_Rb.get_core_form_factor(r, 28);
        res[2][i] = G_Rb.get_core_form_factor(r, 28);
        res[3][i] = T_Rb.get_core_density(sr, 28);
        res[4][i] = G_Rb.get_radial_density(sr);
        res[5][i] = T_Rb.get_radial_density(sr);
        res[6][i] = wavy.compute_dens(sr, 0, 0);
        res[7][i] = wavy2.compute_dens(sr, 0, 0);
        res[8][i] = sr;
        res[9][i] = ECP_way.compute_dens(sr, 0, 0);
        res[10][i] = wavy3.compute_dens(sr, 0, 0);
        res[11][i] = wavy4.compute_dens(sr, 0, 0);
    }
    dat_out << fixed;
    for (int i = 0; i < res[0].size(); i++)
    {
        for (int j = 0; j < res.size(); j++)
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
    return res / 4 / constants::PI;
}

void test_esp_dens()
{
    using namespace std;
    Thakkar T_Li(3);
    ofstream out("core_pot.log", ios::out);
    WFN ECP_way(9);
    ECP_way.read_gbw("Li_QZVP.gbw", out, true, false);
    ECP_way.delete_unoccupied_MOs();

    WFN ECP_way_core(9);
    ECP_way_core.read_gbw("Li_QZVP.gbw", out, true, false);
    ECP_way_core.delete_unoccupied_MOs();
    ECP_way_core.delete_MO(1);
    const int points = 100000;
    const double len = 8;
    const double dr = len / points;

    vec2 res(10);
    for (int i = 0; i < 10; i++)
        res[i].resize(points, 0.0);

    ofstream dat_out("core_pot.dat", ios::out);
    vec2 d2;
    d2.resize(ECP_way.get_ncen());
    for (int i = 0; i < ECP_way.get_ncen(); i++)
    {
        d2[i].resize(ECP_way.get_ncen());
        for (int j = 0; j < ECP_way.get_ncen(); j++)
        {
            if (i == j)
            {
                d2[i][j] = 0;
                continue;
            }
            d2[i][j] = pow(ECP_way.atoms[i].x - ECP_way.atoms[j].x, 2) + pow(ECP_way.atoms[i].y - ECP_way.atoms[j].y, 2) + pow(ECP_way.atoms[i].z - ECP_way.atoms[j].z, 2);
        }
    }
    /*
    const double incr = pow(1.005, max(1, 4));
    const double lincr = log(incr);
    const double min_dist = 0.0000001;
    double current = 1;
    double dist = min_dist;
    vec2 radial_density(1);
    vec2 radial_density_core(1);
    vec2 radial_dist(1);
    const double cube_dist = 5.0;
    while (dist < 2.3 * cube_dist)
    {
        radial_dist[0].push_back(dist);
        current = T_Li.get_radial_density(dist);
        radial_density[0].push_back(current);
        current = T_Li.get_core_density(dist, 2);
        radial_density_core[0].push_back(current);
        dist *= incr;
    }
    vec3 dens_grid;
    double dr = 0.02;
    const double dr3 = dr * dr * dr;
    const int grid_size = (2*cube_dist / dr) + 1;
    dens_grid.resize(grid_size);
    int points = (int)cube_dist / dr;
#pragma omp parallel for
    for (int x = -points; x <= points; x++)
    {
        dens_grid[x + points].resize(grid_size);
        for (int y = -points; y <= points; y++)
        {
            dens_grid[x + points][y + points].resize(grid_size, 0.0);
        }
    }
    double NULLY = 0;
#pragma omp parallel for
    for (int x = -points; x <= points; x++)
    {
        for (int y = -points; y <= points; y++)
        {
            for (int z = -points; z <= points; z++) {
                dens_grid[x+points][y + points][z + points] = linear_interpolate_spherical_density(radial_density[0], radial_dist[0], sqrt(x * x + y * y + z * z)*dr, lincr, 0.0000001);
            }
        }
    }
    out << "Done with radial densities" << endl;
    */

    progress_bar *progress = new progress_bar{out, 50u, "Calculating Values", '-', 0.01};
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < res[0].size(); i++)
    {
        double sr = i * dr;
        double pos[3] = {sr, 0, 0};
        res[0][i] = sr;
        res[1][i] = T_Li.get_core_density(sr, 2);
        res[2][i] = T_Li.get_radial_density(sr);
        res[3][i] = ECP_way.compute_dens(sr, 0, 0);
        res[4][i] = ECP_way.computeESP_noCore(pos, d2);
        // res[5][i] = calc_pot_by_integral(dens_grid, sr, cube_dist, dr);
        // res[6][i] = calc_pot_by_integral(dens_grid, sr, cube_dist, dr);
        res[7][i] = ECP_way_core.compute_dens(sr, 0, 0);
        res[8][i] = ECP_way_core.computeESP_noCore(pos, d2);
        if (i != 0 && i % static_cast<int>(points / 100) == 0)
            progress->write(i / static_cast<double>(res[0].size()));
    }
    dat_out << scientific << setprecision(14);
    for (int i = 0; i < res[0].size(); i++)
    {
        for (int j = 0; j < 9; j++)
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
