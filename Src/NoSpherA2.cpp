#define WIN32_LEAN_AND_MEAN
#include "tsc_block.h"
#include "convenience.h"
#include "JKFit.h"
#include "fchk.h"
#include "cube.h"
#include "scattering_factors.h"
#include "properties.h"
#include "isosurface.h"

bool myGlobalBool = false;
#ifdef _WIN32
bool has_BLAS = false;
#else
bool has_BLAS = true;
#endif
void* BLAS_pointer = nullptr;

int main(int argc, char **argv)
{
    using namespace std;
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
    auto _coutbuf = std::cout.rdbuf(log_file.rdbuf()); // save and redirect
    options opt(argc, argv, log_file);
    opt.digest_options();
    vector<WFN> wavy;
    
    if (opt.threads != -1)
    {
#ifdef _OPENMP
        omp_set_num_threads(opt.threads);
#endif
    }
    log_file << NoSpherA2_message(opt.no_date);
    if (!opt.no_date)
    {
        log_file << build_date;
    }
    log_file.flush();
    // Perform fractal dimensional analysis and quit
    if (opt.fract)
    {
        wavy.push_back(WFN(6));
        cube residual(opt.fract_name, true, wavy[0], std::cout, opt.debug);
        residual.fractal_dimension(0.01);
        log_file.flush();
        std::cout.rdbuf(coutbuf); // reset to standard output again
        std::cout << "Finished!" << endl;
        return 0;
    }
    //Perform Hirshfeld surface based on input and quit
    if (opt.hirshfeld_surface != "") {
        if (opt.radius < 2.5) {
            std::cout << "Resetting Radius to at least 2.5!" << endl;
            opt.radius = 2.5;
        }
        wavy.push_back(WFN(opt.hirshfeld_surface, opt.debug));
        wavy.push_back(WFN(opt.hirshfeld_surface2, opt.debug));
        readxyzMinMax_fromWFN(wavy[0], opt.MinMax, opt.NbSteps, opt.radius, opt.resolution);
        cube Hirshfeld_grid(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy[0].get_ncen(), true);
        cube Hirshfeld_grid2(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy[1].get_ncen(), true);
        Hirshfeld_grid.give_parent_wfn(wavy[0]);
        Hirshfeld_grid2.give_parent_wfn(wavy[1]);
        double len[3]{ 0, 0, 0 };
        for (int i = 0; i < 3; i++)
        {
            len[i] = (opt.MinMax[3 + i] - opt.MinMax[i]) / opt.NbSteps[i];
        }
        for (int i = 0; i < 3; i++)
        {
            Hirshfeld_grid.set_origin(i, opt.MinMax[i]);
            Hirshfeld_grid2.set_origin(i, opt.MinMax[i]);
            Hirshfeld_grid.set_vector(i, i, len[i]);
            Hirshfeld_grid2.set_vector(i, i, len[i]);
        }
        Hirshfeld_grid.set_comment1("Calculated density using NoSpherA2");
        Hirshfeld_grid.set_comment2("from " + wavy[0].get_path().string());
        Hirshfeld_grid2.set_comment1("Calculated density using NoSpherA2");
        Hirshfeld_grid2.set_comment2("from " + wavy[1].get_path().string());
        Calc_Spherical_Dens(Hirshfeld_grid, wavy[0], opt.radius, log_file, false);
        Calc_Spherical_Dens(Hirshfeld_grid2, wavy[1], opt.radius, log_file, false);
        cube Total_Dens = Hirshfeld_grid + Hirshfeld_grid2;
        Total_Dens.give_parent_wfn(wavy[0]);
        cube Hirshfeld_weight = Hirshfeld_grid / Total_Dens;
        Hirshfeld_weight.give_parent_wfn(wavy[0]);
        std::array<std::array<int, 3>, 3> Colourcode;
        
        Colourcode[0] = { 255, 0, 0 };
        Colourcode[1] = { 255, 255, 255 };
        Colourcode[2] = { 0, 0, 255 };

        std::vector<Triangle> triangles_i = marchingCubes(Hirshfeld_weight, 0.5, 1);
        std::cout << "Found " << triangles_i.size() << " triangles!" << endl;
        auto triangles_e = triangles_i;
        double area = 0.0;
        double volume = 0.0;
        double low_lim_di = 1E7;
        double high_lim_di = 0.0;
        double low_lim_de = 1E7;
        double high_lim_de = 0.0;
        ofstream fingerprint_file("Hirshfeld_fingerprint.dat");
        fingerprint_file << "d_i\td_e" << endl;
#pragma omp parallel for reduction(+:area, volume)
        for (int i = 0; i < triangles_i.size(); i++) {
            area += triangles_i[i].calc_area();
            volume += triangles_i[i].calc_inner_volume();
            std::array<double, 3> pos = triangles_i[i].calc_center();
            double d_i = calc_d_i(pos, wavy[0]);
            double d_e = calc_d_i(pos, wavy[1]);
#pragma omp critical
            {
                if (d_i < low_lim_di)
                    low_lim_di = d_i;
                if (d_i > high_lim_di)
                    high_lim_di = d_i;
                if (d_e < low_lim_de)
                    low_lim_de = d_e;
                if (d_e > high_lim_de)
                    high_lim_de = d_e;
                fingerprint_file << d_i << "\t" << d_e << "\n";
            }
        }
        fingerprint_file.flush();
        fingerprint_file.close();
        std::cout << "d_i is scaled from " << low_lim_di << " to " << high_lim_di * 0.9 << endl;
        std::cout << "d_e is scaled from " << low_lim_de << " to " << high_lim_de * 0.9 << endl;
#pragma omp parallel for
        for (int i=0; i<triangles_i.size(); i++)
        {
            triangles_i[i].colour = get_colour(triangles_i[i], calc_d_i, wavy[0], Colourcode, low_lim_di, high_lim_di*0.9);
            triangles_e[i].colour = get_colour(triangles_e[i], calc_d_i, wavy[1], Colourcode, low_lim_de, high_lim_de*0.9);
        }
        std::cout << "Total area: " << area << endl;
        std::cout << "Total volume: " << volume << endl;
        writeColourObj("Hirshfeld_surface_i.obj", triangles_i);
        writeColourObj("Hirshfeld_surface_e.obj", triangles_e);
        std::cout.rdbuf(coutbuf); // reset to standard output again
        std::cout << "Finished!" << endl;
        return 0;
    }
    // Perform calcualtion of difference between two wavefunctions using the resolution, radius, wfn and wfn2 keywords. wfn2 keaword is provided by density-difference flag
    if (!opt.wfn2.empty())
    {
        wavy.push_back(WFN(opt.wfn, opt.debug));
        wavy.push_back(WFN(opt.wfn2, opt.debug));
        if (opt.debug)
            cout << opt.wfn << opt.wfn2 << endl;
        wavy[0].delete_unoccupied_MOs();
        wavy[1].delete_unoccupied_MOs();
        readxyzMinMax_fromWFN(wavy[0], opt.MinMax, opt.NbSteps, opt.radius, opt.resolution);
        cube Rho1(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy[0].get_ncen(), true);
        cube Rho2(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy[0].get_ncen(), true);
        Rho1.give_parent_wfn(wavy[0]);
        Rho2.give_parent_wfn(wavy[1]);
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
        Rho1.set_comment2("from " + wavy[0].get_path().string());
        Rho2.set_comment1("Calculated density using NoSpherA2");
        Rho2.set_comment2("from " + wavy[1].get_path().string());
        Rho1.set_path(std::filesystem::path(wavy[0].get_path().stem().string() + "_rho.cube"));
        Rho2.set_path(std::filesystem::path(wavy[1].get_path().stem().string() + "_rho.cube"));
        Calc_Rho(Rho1, wavy[0], opt.radius, log_file, false);
        Calc_Rho(Rho2, wavy[1], opt.radius, log_file, false);
        cube Rho_diff = Rho1 - Rho2;
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
        Rho_diff.give_parent_wfn(wavy[0]);
        cout << "RSR between the two cubes: " << setw(16) << scientific << setprecision(16) << Rho1.rrs(Rho2) << endl;
        cout << "Ne of shifted electrons: " << Rho_diff.diff_sum() << endl;
        cout << "Writing cube 1..." << flush;
        Rho1.write_file(Rho1.get_path(), false);
        cout << " ... done!\nWriting cube 2..." << flush;
        Rho2.write_file(Rho2.get_path(), false);
        Rho_diff.set_path(wavy[1].get_path().stem().string() + "_diff.cube");
        cout << " ... done\nWriting difference..." << flush;
        Rho_diff.write_file(Rho_diff.get_path(), false);
        cout << " ... done :)" << endl;
        cout << "Bye Bye!" << endl;
        return 0;
    }
    if (opt.pol_wfns.size() != 0) {
        polarizabilities(opt, log_file);
        exit(0);
    }
    // Performs MTC and CMTC calcualtions, that is multiple wfns with either one or multiple cifs and 1 common hkl.
    if (opt.cif_based_combined_tsc_calc || opt.combined_tsc_calc)
    {
        err_checkf(opt.hkl != "" || opt.dmin != 99.0 || opt.hkl_min_max[0][0] != -100, "No hkl specified and no dmin value given", log_file);
        if (opt.combined_tsc_calc)
            err_checkf(opt.cif != "", "No cif specified", log_file);
        // First make sure all files exist
        if (opt.cif_based_combined_tsc_calc)
        {
            err_checkf(opt.combined_tsc_calc_files.size() == opt.combined_tsc_calc_cifs.size(), "Unequal number of CIFs and WFNs impossible!", log_file);
        }
        err_checkf(opt.combined_tsc_calc_mult.size() == opt.combined_tsc_calc_files.size(), "Unequal number of WFNs and mults impossible!", log_file);
        err_checkf(opt.combined_tsc_calc_charge.size() == opt.combined_tsc_calc_files.size(), "Unequal number of WFNs and charges impossible!", log_file);

        for (int i = 0; i < opt.combined_tsc_calc_files.size(); i++)
        {
            err_checkf(std::filesystem::exists(opt.combined_tsc_calc_files[i]), "Specified file for combined calculation doesn't exist! " + opt.combined_tsc_calc_files[i].string(), log_file);
            if (opt.cif_based_combined_tsc_calc)
                err_checkf(std::filesystem::exists(opt.combined_tsc_calc_cifs[i]), "Specified file for combined calculation doesn't exist! " + opt.combined_tsc_calc_cifs[i].string(), log_file);
        }
        wavy.resize(opt.combined_tsc_calc_files.size());
        for (int i = 0; i < opt.combined_tsc_calc_files.size(); i++)
        {
            log_file << "Reading: " << setw(44) << opt.combined_tsc_calc_files[i] << flush;
            if (opt.debug)
            {
                log_file << "\nm: " << opt.combined_tsc_calc_mult[i] << endl;
                log_file << "c: " << opt.combined_tsc_calc_charge[i] << "\n";
            }
            wavy[i].set_multi(opt.combined_tsc_calc_mult[i]);
            wavy[i].set_charge(opt.combined_tsc_calc_charge[i]);
            wavy[i].read_known_wavefunction_format(opt.combined_tsc_calc_files[i], log_file, opt.debug);
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

        svec known_scatterer;
        vec2 known_kpts;
        tsc_block<int, cdouble> result;
        for (int i = 0; i < opt.combined_tsc_calc_files.size(); i++)
        {
            known_scatterer = result.get_scatterers();
            if (wavy[i].get_origin() != 7 && !opt.SALTED)
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
            else if (opt.SALTED) {
#if has_RAS == 1
                // Fill WFN wil the primitives of the JKFit basis (currently hardcoded)
                // const std::vector<std::vector<primitive>> basis(QZVP_JKfit.begin(), QZVP_JKfit.end());
#ifdef _WIN32
                check_OpenBLAS_DLL(opt.debug);
#endif
                BasisSetLibrary basis_library;
                SALTEDPredictor* temp_pred = new SALTEDPredictor(wavy[i], opt);
                string df_basis_name = temp_pred->get_dfbasis_name();
                filesystem::path h5file = temp_pred->get_h5_filename();
                log_file << "Using " << h5file << " for the prediction" << endl;
                load_basis_into_WFN(temp_pred->wavy, basis_library.get_basis_set(df_basis_name));

                if (opt.debug)
                    log_file << "Entering scattering ML Factor Calculation with H part!" << endl;
                result.append(calculate_scattering_factors_MTC_SALTED(
                    opt,
					*temp_pred,
                    log_file,
                    known_scatterer,
                    i,
                    &known_kpts), log_file);
                delete temp_pred;
#else
                log_file << "SALTED is not available in this build!" << endl;
                exit(-1);
#endif
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
        wavy.emplace_back(opt.xyz_file, opt.debug);

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
        log_file << "Reading: " << setw(44) << opt.wfn << flush;
        wavy.emplace_back(opt.wfn, opt.charge, opt.mult, opt.debug);
        wavy[0].set_method(opt.method);
        wavy[0].set_multi(opt.mult);
        wavy[0].set_charge(opt.charge);
        if (opt.debug)
            log_file << "method/mult/charge: " << opt.method << " " << opt.mult << " " << opt.charge << endl;
        
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

        // this one is for generation of an fchk file
        if (opt.basis_set != "" || opt.fchk != "")
        {
            // Make a fchk out of the wfn/wfx file
            filesystem::path tmp = opt.basis_set_path / opt.basis_set;
            if (opt.debug)
                log_file << "Checking for " << opt.basis_set_path << " " << exists(opt.basis_set_path) << endl;
            // err_checkf(exists(opt.basis_set_path), "Basis set file does not exist!", log_file);
            wavy[0].set_basis_set_name(tmp.string());

            std::filesystem::path outputname;
            if (opt.fchk != "")
                outputname = opt.fchk;
            else
            {
                outputname = wavy[0].get_path();
                outputname.replace_extension(".fchk");
            }
            wavy[0].assign_charge(wavy[0].calculate_charge());
            if (opt.mult == 0)
                err_checkf(wavy[0].guess_multiplicity(log_file), "Error guessing multiplicity", log_file);
            free_fchk(log_file, outputname, "", wavy[0], opt.debug, true);
        }

        // This one will calcualte a single tsc/tscb file form a single wfn
        if (opt.cif != "" || opt.hkl != "")
        {
            if (!opt.SALTED)
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
            else
            {
#if has_RAS == 1
                // Fill WFN wil the primitives of the JKFit basis (currently hardcoded)
                // const std::vector<std::vector<primitive>> basis(QZVP_JKfit.begin(), QZVP_JKfit.end());
#ifdef _WIN32
                check_OpenBLAS_DLL(opt.debug);
#endif
                BasisSetLibrary basis_library;
                SALTEDPredictor* temp_pred = new SALTEDPredictor(wavy[0], opt);
                string df_basis_name = temp_pred->get_dfbasis_name();
                filesystem::path h5file = temp_pred->get_h5_filename();
                log_file << "Using " << h5file << " for the prediction" << endl;
                load_basis_into_WFN(temp_pred->wavy, basis_library.get_basis_set(df_basis_name));

                if (opt.debug)
                    log_file << "Entering scattering ML Factor Calculation with H part!" << endl;
                err_checkf(calculate_scattering_factors_ML(
                    opt,
                    *temp_pred,
                    log_file),
                    "Error during ML-SF Calcualtion", log_file);
                delete temp_pred;
#else
							log_file << "SALTED is not available in this build!" << endl;
              exit(-1);
#endif
            }
        }
        log_file.flush();
        std::cout.rdbuf(_coutbuf); // reset to standard output again
        std::cout << "Finished!" << endl;
        if (opt.write_CIF)
            wavy[0].write_wfn_CIF(opt.wfn.replace_extension(".cif"));
        // log_file.close();
        return 0;
    }
    // Contains all calculations of properties and cubes
    if (opt.calc)
    {
        properties_calculation(opt);
        log_file.flush();
        std::cout.rdbuf(_coutbuf); // reset to standard output again
        std::cout << "Finished!" << endl;
        return 0;
    }
    // Converts gbw file to wfn file and leaves
    if (opt.gbw2wfn)
    {
        err_checkf(opt.wfn != "", "No Wavefunction given!", log_file);
        wavy.emplace_back(opt.wfn, opt.debug);
        wavy[0].write_wfn("converted.wfn", false, false);
        log_file.flush();
        std::cout.rdbuf(_coutbuf); // reset to standard output again
        std::cout << "Finished!" << endl;
        if (opt.write_CIF)
            wavy[0].write_wfn_CIF(opt.wfn.replace_extension(".cif"));
        return 0;
    }
    std::cout << NoSpherA2_message(opt.no_date);
    if (!opt.no_date)
        std::cout << build_date;
    std::cout << "Did not understand the task to perform!\n"
              << help_message << endl;
    log_file.flush();
    return 0;
}