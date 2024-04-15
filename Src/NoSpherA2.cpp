#define WIN32_LEAN_AND_MEAN
#include "tsc_block.h"
#include "convenience.h"
#include "fchk.h"
#include "cube.h"
#include "scattering_factors.h"
#include "properties.h"

using namespace std;

int main(int argc, char **argv)
{
    std::cout << NoSpherA2_message();
    char cwd[1024];
#ifdef _WIN32
    if (_getcwd(cwd, sizeof(cwd)) != NULL)
#else
    if (getcwd(cwd, sizeof(cwd)) != NULL)
#endif
    {
        std::cout << "Current working directory: " << cwd << std::endl;
    }
    else
    {
        std::cerr << "getcwd() error" << std::endl;
        return 1;
    }
    ofstream log_file("NoSpherA2.log", ios::out);
    auto coutbuf = std::cout.rdbuf(log_file.rdbuf()); // save and redirect
    vector<WFN> wavy;
    options opt(argc, argv, log_file);
    opt.digest_options();
    if (opt.threads != -1)
    {
#ifdef _OPENMP
        omp_set_num_threads(opt.threads);
#endif
    }
    log_file << NoSpherA2_message();
    if (!opt.no_date)
    {
        log_file << build_date();
    }
    log_file.flush();
    // Perform fractal dimensional analysis and quit
    if (opt.fract)
    {
        WFN wav(6);
        cube residual(opt.fract_name, true, wav, std::cout, opt.debug);
        residual.fractal_dimension(0.01);
        log_file.flush();
        std::cout.rdbuf(coutbuf); // reset to standard output again
        std::cout << "Finished!" << endl;
        return 0;
    }
    // Perform calcualtion of difference between two wavefunctions using the resolution, radius, wfn and wfn2 keywords. wfn2 keaword is provided by density-difference flag
    if (opt.wfn2.length() != 0)
    {
        WFN wav(9), wav2(9);
        wav.read_known_wavefunction_format(opt.wfn, log_file, opt.debug), wav2.read_known_wavefunction_format(opt.wfn2, log_file, opt.debug);
        if (opt.debug)
            cout << opt.wfn << opt.wfn2 << endl;
        wav.delete_unoccupied_MOs();
        wav2.delete_unoccupied_MOs();
        readxyzMinMax_fromWFN(wav, opt.MinMax, opt.NbSteps, opt.radius, opt.resolution);
        cube Rho1(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wav.get_ncen(), true);
        cube Rho2(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wav.get_ncen(), true);
        Rho1.give_parent_wfn(wav);
        Rho2.give_parent_wfn(wav2);
        double len[3]{0, 0, 0};
        for (int i = 0; i < 3; i++)
        {
            len[i] = (opt.MinMax[3 + i] - opt.MinMax[i]) / opt.NbSteps[i];
        }
        for (int i = 0; i < 3; i++)
        {
            Rho1.set_origin(i, opt.MinMax[i]);
            Rho2.set_origin(i, opt.MinMax[i]);
            Rho1.set_vector(i, i, len[i]);
            Rho2.set_vector(i, i, len[i]);
        }
        Rho1.set_comment1("Calculated density using NoSpherA2");
        Rho1.set_comment2("from " + wav.get_path());
        Rho2.set_comment1("Calculated density using NoSpherA2");
        Rho2.set_comment2("from " + wav2.get_path());
        Rho1.path = get_basename_without_ending(wav.get_path()) + "_rho.cube";
        Rho2.path = get_basename_without_ending(wav2.get_path()) + "_rho.cube";
        Calc_Rho_no_trans(Rho1, wav, opt.threads, opt.radius, log_file);
        Calc_Rho_no_trans(Rho2, wav2, opt.threads, opt.radius, log_file);
        cube Rho_diff(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wav.get_ncen(), true);
#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < Rho1.get_size(0); i++)
        {
            for (int j = 0; j < Rho1.get_size(1); j++)
                for (int k = 0; k < Rho1.get_size(2); k++)
                {
                    Rho_diff.set_value(i, j, k, Rho1.get_value(i, j, k) - Rho2.get_value(i, j, k));
                }
        }
        for (int i = 0; i < 3; i++)
        {
            Rho_diff.set_origin(i, opt.MinMax[i]);
            Rho_diff.set_vector(i, i, len[i]);
        }
        Rho_diff.give_parent_wfn(wav);
        cout << "RSR between the two cubes: " << setw(16) << scientific << setprecision(16) << Rho1.rrs(Rho2) << endl;
        cout << "Ne of shifted electrons: " << Rho_diff.diff_sum() << endl;
        cout << "Writing cube 1..." << flush;
        Rho1.write_file(Rho1.path, false);
        cout << " ... done!\nWriting cube 2..." << flush;
        Rho2.write_file(Rho2.path, false);
        Rho_diff.path = get_basename_without_ending(wav2.get_path()) + "_diff.cube";
        cout << " ... done\nWriting difference..." << flush;
        Rho_diff.write_file(Rho_diff.path, false);
        cout << " ... done :)" << endl;
        cout << "Bye Bye!" << endl;
        return 0;
    }
    // Performs MTC and CMTC calcualtions, that is multiple wfns with either one or multiple cifs and 1 common hkl.
    if (opt.cif_based_combined_tsc_calc || opt.combined_tsc_calc)
    {
        err_checkf(opt.hkl != "" || opt.dmin != 99.0, "No hkl specified and no dmin value given", log_file);
        if (opt.combined_tsc_calc)
            err_checkf(opt.cif != "", "No cif specified", log_file);
        // First make sure all files exist
        if (opt.cif_based_combined_tsc_calc)
        {
            err_checkf(opt.combined_tsc_calc_files.size() == opt.combined_tsc_calc_cifs.size(), "Unequal number of CIFs and WFNs impossible!", log_file);
            err_checkf(opt.combined_tsc_calc_mult.size() == opt.combined_tsc_calc_cifs.size(), "Unequal number of CIFs and mults impossible!", log_file);
            err_checkf(opt.combined_tsc_calc_charge.size() == opt.combined_tsc_calc_cifs.size(), "Unequal number of CIFs and charges impossible!", log_file);
        }
        for (int i = 0; i < opt.combined_tsc_calc_files.size(); i++)
        {
            err_checkf(exists(opt.combined_tsc_calc_files[i]), "Specified file for combined calculation doesn't exist! " + opt.combined_tsc_calc_files[i], log_file);
            if (opt.cif_based_combined_tsc_calc)
                err_checkf(exists(opt.combined_tsc_calc_cifs[i]), "Specified file for combined calculation doesn't exist! " + opt.combined_tsc_calc_cifs[i], log_file);
        }
        wavy.resize(opt.combined_tsc_calc_files.size());
        for (int i = 0; i < opt.combined_tsc_calc_files.size(); i++)
        {
            log_file << "Reading: " << setw(44) << opt.combined_tsc_calc_files[i] << flush;
            if (opt.debug)
            {
                log_file << "m: " << opt.combined_tsc_calc_mult[i] << endl;
                log_file << "c: " << opt.combined_tsc_calc_charge[i] << "\n";
            }
            wavy[i].set_multi(opt.combined_tsc_calc_mult[i]);
            wavy[i].set_charge(opt.combined_tsc_calc_charge[i]);
            wavy[i].read_known_wavefunction_format(opt.combined_tsc_calc_files[i], log_file);
            if (opt.ECP)
            {
                wavy[i].set_has_ECPs(true, true, opt.ECP_mode);
            }
            if (opt.set_ECPs)
            {
                wavy[i].set_ECPs(opt.ECP_nrs, opt.ECP_elcounts);
            }
            log_file << " done!\nNumber of atoms in Wavefunction file: " << wavy[i].get_ncen() << " Number of MOs: " << wavy[i].get_nmo() << endl;
        }

        vector<string> known_scatterer;
        vector<vec> known_kpts;
        tsc_block<int, cdouble> result;
        for (int i = 0; i < opt.combined_tsc_calc_files.size(); i++)
        {
            known_scatterer = result.get_scatterers();
            if (wavy[i].get_origin() != 7)
            {
                result.append(calculate_scattering_factors_MTC(
                                  opt,
                                  wavy,
                                  log_file,
                                  known_scatterer,
                                  i,
                                  &known_kpts),
                              log_file);
            }
            else
            {
                result.append(MTC_thakkar_sfac(
                                  opt,
                                  log_file,
                                  known_scatterer,
                                  wavy,
                                  i),
                              log_file);
            }
        }

        known_scatterer = result.get_scatterers();
        log_file << "Final number of atoms in .tsc file: " << known_scatterer.size() << endl;
        time_point start = get_time();
        log_file << "Writing tsc file... " << flush;
        if (opt.binary_tsc)
            result.write_tscb_file();
        if (opt.old_tsc)
        {
            result.write_tsc_file(opt.cif);
        }
        log_file << " ... done!" << endl;
        if (!opt.no_date)
        {
            time_point end_write = get_time();
            if (get_sec(start, end_write) < 60)
                log_file << "Writing Time: " << fixed << setprecision(0) << get_sec(start, end_write) << " s\n";
            else if (get_sec(start, end_write) < 3600)
                log_file << "Writing Time: " << fixed << setprecision(0) << floor(get_sec(start, end_write) / 60) << " m " << get_sec(start, end_write) % 60 << " s\n";
            else
                log_file << "Writing Time: " << fixed << setprecision(0) << floor(get_sec(start, end_write) / 3600) << " h " << (get_sec(start, end_write) % 3600) / 60 << " m\n";
            log_file << endl;
        }
        log_file.flush();
        std::cout.rdbuf(coutbuf); // reset to standard output again
        std::cout << "Finished!" << endl;
        return 0;
    }
    // Performs the Thakkar IAM
    if (opt.iam_switch)
    {
        if (opt.debug)
        {
            log_file << "I am doing a Thakkar IAM!" << endl;
        }
        err_checkf(opt.xyz_file != "", "No xyz specified", log_file);
        err_checkf(exists(opt.xyz_file), "xyz doesn't exist", log_file);
        auto t = new WFN(0);
        wavy.push_back(*t);
        delete t;
        wavy[0].read_known_wavefunction_format(opt.xyz_file, log_file, opt.debug);

        if (opt.electron_diffraction && opt.debug)
            log_file << "Making Electron diffraction scattering factors, be carefull what you are doing!" << endl;
        if (opt.debug)
            log_file << "Entering scattering Factor Calculation!" << endl;
        err_checkf(thakkar_sfac(opt, log_file, wavy[0]), "Error during SF Calculation!", log_file);
        log_file.flush();
        std::cout.rdbuf(coutbuf); // reset to standard output again
        std::cout << "Finished!" << endl;
        return 0;
    }
    // This one has conversion to fchk and calculation of one single tsc file
    if (opt.wfn != "" && !opt.calc && !opt.gbw2wfn && opt.d_sfac_scan == 0.0)
    {
        auto t = new WFN(0);
        wavy.push_back(*t);
        delete t;
        wavy[0].set_method(opt.method);
        wavy[0].set_multi(opt.mult);
        wavy[0].set_charge(opt.charge);
        if (opt.debug)
            log_file << "method/mult/charge: " << opt.method << " " << opt.mult << " " << opt.charge << endl;
        log_file << "Reading: " << setw(44) << opt.wfn << flush;
        wavy[0].read_known_wavefunction_format(opt.wfn, log_file, opt.debug);
        if (opt.ECP)
        {
            wavy[0].set_has_ECPs(true, true, opt.ECP_mode);
        }
        if (opt.set_ECPs)
        {
            log_file << "Adding ECPs" << endl;
            wavy[0].set_ECPs(opt.ECP_nrs, opt.ECP_elcounts);
        }
        log_file << " done!\nNumber of atoms in Wavefunction file: " << wavy[0].get_ncen() << " Number of MOs: " << wavy[0].get_nmo() << endl;

        if (opt.basis_set != "" || opt.fchk != "")
        {
            // Make a fchk out of the wfn/wfx file
            join_path(opt.basis_set_path, opt.basis_set);
            if (opt.debug)
                log_file << "Checking for " << opt.basis_set_path << " " << exists(opt.basis_set_path) << endl;
            // err_checkf(exists(opt.basis_set_path), "Basis set file does not exist!", log_file);
            wavy[0].set_basis_set_name(opt.basis_set_path);

            string outputname;
            if (opt.fchk != "")
                outputname = opt.fchk;
            else
            {
                outputname = wavy[0].get_path();
                if (opt.debug)
                    log_file << "Loaded path..." << endl;
                size_t where(0);
                int len = 4;
                if (wavy[0].get_origin() == 2)
                    where = outputname.find(".wfn");
                else if (wavy[0].get_origin() == 4)
                    where = outputname.find(".ffn");
                else if (wavy[0].get_origin() == 6)
                    where = outputname.find(".wfx");
                else if (wavy[0].get_origin() == 8)
                    where = outputname.find(".molden"), len = 7;
                else if (wavy[0].get_origin() == 9)
                    where = outputname.find(".gbw");
                if (where >= outputname.length() && where != string::npos)
                    err_checkf(false, "Cannot make output file name!", log_file);
                else
                    outputname.erase(where, len);
            }
            wavy[0].assign_charge(wavy[0].calculate_charge());
            if (opt.mult == 0)
                err_checkf(wavy[0].guess_multiplicity(log_file), "Error guessing multiplicity", log_file);
            else
                wavy[0].assign_multi(opt.mult);
            free_fchk(log_file, outputname, "", wavy[0], opt.debug, true);
        }
        if (opt.cif != "" || opt.hkl != "")
        {
            // Calculate tsc file from given files
            if (opt.debug)
                log_file << "Entering scattering Factor Calculation!" << endl;
            if (opt.electron_diffraction)
                log_file << "Making Electron diffraction scattering factors, be carefull what you are doing!" << endl;
            if (wavy[0].get_origin() != 7)
                err_checkf(calculate_scattering_factors_HF(
                               opt,
                               wavy[0],
                               log_file),
                           "Error during SF Calcualtion", log_file);
            else
                err_checkf(thakkar_sfac(
                               opt,
                               log_file,
                               wavy[0]),
                           "Error during SF Calcualtion", log_file);
        }
        log_file.flush();
        std::cout.rdbuf(coutbuf); // reset to standard output again
        std::cout << "Finished!" << endl;
        if (opt.write_CIF)
            wavy[0].write_wfn_CIF(opt.wfn + ".cif");
        // log_file.close();
        return 0;
    }
    // Contains all calculations of properties and cubes
    if (opt.calc)
    {
        properties_calculation(opt);
        log_file.flush();
        std::cout.rdbuf(coutbuf); // reset to standard output again
        std::cout << "Finished!" << endl;
        return 0;
    }
    // Converts gbw file to wfn file and leaves
    if (opt.gbw2wfn)
    {
        err_checkf(opt.wfn != "", "No Wavefunction given!", log_file);
        auto t = new WFN(9);
        wavy.push_back(*t);
        delete t;
        wavy[0].read_known_wavefunction_format(opt.wfn, log_file);
        wavy[0].write_wfn("converted.wfn", false, false);
        log_file.flush();
        std::cout.rdbuf(coutbuf); // reset to standard output again
        std::cout << "Finished!" << endl;
        if (opt.write_CIF)
            wavy[0].write_wfn_CIF(opt.wfn + ".cif");
        return 0;
    }
    if (opt.coef_file != "")
    {
        if (opt.cif != "" && opt.xyz_file != "")
        {
            auto t = new WFN(7);
            wavy.push_back(*t);
            delete t;
            // Read the xyzfile
            wavy[0].read_xyz(opt.xyz_file, std::cout, opt.debug);
            // Fill WFN wil the primitives of the JKFit basis (currently hardcoded)
            int nr_coefs = load_basis_into_WFN(wavy[0], TZVP_JKfit);
            // Run generation of tsc file
            if (opt.SALTED_BECKE || opt.SALTED)
                err_checkf(calculate_scattering_factors_RI_No_H(
                               opt,
                               wavy[0],
                               log_file,
                               nr_coefs),
                           "Error during ML-SF Calcualtion", log_file);
            else
                err_checkf(calculate_scattering_factors_RI(
                               opt,
                               wavy[0],
                               log_file,
                               nr_coefs),
                           "Error during ML-SF Calcualtion", log_file);
        }
        else
        {
            // #include "test_functions.h"
            //       cube_from_coef_npy(opt.coef_file, opt.xyz_file);
        }
        std::cout.rdbuf(coutbuf); // reset to standard output again
        std::cout << "Finished!" << endl;
        return 0;
    }
    std::cout << NoSpherA2_message();
    if (!opt.no_date)
        std::cout << build_date();
    std::cout << "Did not understand the task to perform!\n"
              << help_message() << endl;
    log_file.flush();
    return 0;
}