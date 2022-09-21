
#include "convenience.h"
#include "fchk.h"
#include "cube.h"
#include "basis_set.h"
#include "structure_factors.h"
#include "properties.h"
#include "cell.h"

using namespace std;
//bool expert = false;

int main(int argc, char** argv)
{
  ofstream log_file("NoSpherA2.log", ios::out);
  auto coutbuf = std::cout.rdbuf(log_file.rdbuf()); //save and redirect
  vector<WFN> wavy;
  options opt(argc, argv);
  if (opt.debug)
    cout << "argc:" << argc << endl;
  if (opt.threads != -1) {
    omp_set_num_threads(opt.threads);
    omp_set_dynamic(0);
  }
  //See if user wants help
  if (argc > 1) {
    if (string(argv[1]).find("--h") != string::npos) {
      cout << NoSpherA2_message(opt.no_date) << help_message() << endl;
      return 0;
    }
    if (string(argv[1]).find("-help") != string::npos) {
      cout << NoSpherA2_message(opt.no_date) << help_message() << endl;
      return 0;
    }
  }
  //Lets print what was the command line, for debugging
  if (opt.debug) {
    ofstream log("NoSpherA2.log", ios::out);
    log << " Recap of input:\n";
    for (int i = 0; i < argc; i++) {
      //cout << argv[i] << endl;
      log << argv[i] << "\n";
    }
    log.flush();
    log.close();
  }
  if (opt.fract) {
    WFN wavy(6);
    cube residual(opt.fract_name, true, wavy, cout, opt.debug);
    residual.fractal_dimension(0.01);
    return 0;
  }
  if (opt.combined_tsc_calc) {
    log_file << NoSpherA2_message(opt.no_date);
    log_file.flush();
    err_checkf(opt.hkl != "", "No hkl specified", log_file);
    err_checkf(opt.cif != "", "No cif specified", log_file);
    //First make sure all files exist
    for (int i = 0; i < opt.combined_tsc_calc_files.size(); i++)
      err_checkf(exists(opt.combined_tsc_calc_files[i]), "Specified file for combined calculation doesn't exist! " + opt.combined_tsc_calc_files[i], log_file);

    wavy.resize(opt.combined_tsc_calc_files.size());
    for (int i = 0; i < opt.combined_tsc_calc_files.size(); i++) {
      log_file << "Reading: " << setw(44) << opt.combined_tsc_calc_files[i] << flush;
      wavy[i].read_known_wavefunction_format(opt.combined_tsc_calc_files[i], log_file);
      if (opt.ECP) {
        wavy[i].set_has_ECPs(true);
      }
      log_file << " done!" << endl << "Number of atoms in Wavefunction file: " << wavy[i].get_ncen() << " Number of MOs: " << wavy[i].get_nmo() << endl;
    }

    vector<string> known_scatterer;
    vector<vector<int>> known_indices;
    tsc_block result;
    for (int i = 0; i < opt.combined_tsc_calc_files.size(); i++) {
      known_scatterer = result.get_scatterers();
      known_indices = result.get_index_vector();
      if (wavy[i].get_origin() != 7) {
        result.append(calculate_structure_factors_MTC(
          opt,
          wavy,
          log_file,
          known_scatterer,
          known_indices,
          i
        ), log_file);
      }
      else {
        result.append(MTC_thakkar_sfac(
          opt,
          log_file,
          known_scatterer,
          known_indices,
          wavy,
          i
        ), log_file);
      }
    }

    if (opt.debug) {
      for (int i = 0; i < opt.combined_tsc_calc_files.size(); i++) {
        log_file << opt.combined_tsc_calc_files[i] << " ";
        for (int j = 0; j < opt.groups[i].size(); j++) log_file << opt.groups[i][j] << " ";
        log_file << endl;
      }
    }
    known_scatterer = result.get_scatterers();
    log_file << "Final number of atoms in .tsc file: " << known_scatterer.size() << endl;
#ifdef _WIN64
    time_t start = time(NULL);
    time_t end_write;
#endif
    log_file << "Writing tsc file... " << flush;
    if (opt.binary_tsc)
      result.write_tscb_file();
    if (opt.old_tsc) {
      result.write_tsc_file(opt.cif);
    }
    log_file << " ... done!" << endl;
#ifdef _WIN64
    end_write = time(NULL);

    if (end_write - start < 60) log_file << "Writing Time: " << fixed << setprecision(0) << end_write - start << " s\n";
    else if (end_write - start < 3600) log_file << "Writing Time: " << fixed << setprecision(0) << floor((end_write - start) / 60) << " m " << (end_write - start) % 60 << " s\n";
    else log_file << "Writing Time: " << fixed << setprecision(0) << floor((end_write - start) / 3600) << " h " << ((end_write - start) % 3600) / 60 << " m\n";
    log_file << endl;
#endif
    log_file.close();
    return 0;
  }
  if (opt.cif_based_combined_tsc_calc) {
    log_file << NoSpherA2_message(opt.no_date);
    log_file.flush();
    err_checkf(opt.hkl != "", "No hkl specified", log_file);
    //First make sure all files exist
    err_checkf(opt.combined_tsc_calc_files.size() == opt.combined_tsc_calc_cifs.size(), "Unequal number of CIFs and WFNs impossible!", log_file);
    for (int i = 0; i < opt.combined_tsc_calc_files.size(); i++) {
      err_checkf(exists(opt.combined_tsc_calc_files[i]), "Specified file for combined calculation doesn't exist! " + opt.combined_tsc_calc_files[i], log_file);
      err_checkf(exists(opt.combined_tsc_calc_cifs[i]), "Specified file for combined calculation doesn't exist! " + opt.combined_tsc_calc_cifs[i], log_file);
    }

    wavy.resize(opt.combined_tsc_calc_files.size());
    for (int i = 0; i < opt.combined_tsc_calc_files.size(); i++) {
      log_file << "Reading: " << setw(44) << opt.combined_tsc_calc_files[i] << flush;
      wavy[i].read_known_wavefunction_format(opt.combined_tsc_calc_files[i], log_file);
      if (opt.ECP) {
        wavy[i].set_has_ECPs(true);
      }
      log_file << " done!" << endl << "Number of atoms in Wavefunction file: " << wavy[i].get_ncen() << " Number of MOs: " << wavy[i].get_nmo() << endl;
    }

    vector<string> known_scatterer;
    vector<vector<int>> known_indices;
    tsc_block result;
    for (int i = 0; i < opt.combined_tsc_calc_files.size(); i++) {
      known_scatterer = result.get_scatterers();
      known_indices = result.get_index_vector();
      if (wavy[i].get_origin() != 7) {
        result.append(calculate_structure_factors_MTC(
          opt,
          wavy,
          log_file,
          known_scatterer,
          known_indices,
          i
        ), log_file);
      }
      else {
        result.append(MTC_thakkar_sfac(
          opt,
          log_file,
          known_scatterer,
          known_indices,
          wavy,
          i
        ), log_file);
      }
    }

    if (opt.debug) {
      for (int i = 1; i < opt.combined_tsc_calc_files.size(); i++) {
        log_file << opt.combined_tsc_calc_files[i] << " ";
        for (int j = 0; j < opt.groups[i].size(); j++) log_file << opt.groups[i][j] << " ";
        log_file << endl;
      }
    }
    known_scatterer = result.get_scatterers();
    log_file << "Final number of atoms in .tsc file: " << known_scatterer.size() << endl;
#ifdef _WIN64
    time_t start = time(NULL);
    time_t end_write;
#endif
    log_file << "Writing tsc file... " << flush;
    if (opt.binary_tsc)
      result.write_tscb_file();
    if (opt.old_tsc) {
      result.write_tsc_file(opt.cif);
    }
    log_file << " ... done!" << endl;
#ifdef _WIN64
    end_write = time(NULL);

    if (end_write - start < 60) log_file << "Writing Time: " << fixed << setprecision(0) << end_write - start << " s\n";
    else if (end_write - start < 3600) log_file << "Writing Time: " << fixed << setprecision(0) << floor((end_write - start) / 60) << " m " << (end_write - start) % 60 << " s\n";
    else log_file << "Writing Time: " << fixed << setprecision(0) << floor((end_write - start) / 3600) << " h " << ((end_write - start) % 3600) / 60 << " m\n";
    log_file << endl;
#endif
    log_file.close();
    return 0;
  }
  if (opt.iam_switch) {
    ofstream log_file("NoSpherA2.log", ios::out);
    if (opt.debug)
      for (int i = 0; i < argc; i++)
        log_file << argv[i] << endl;

    log_file << NoSpherA2_message(opt.no_date);
    log_file.flush();
    if (argc > 1) {
      if (string(argv[1]).find("--help") != string::npos) {
        log_file << help_message() << endl;
        log_file.close();
        return 0;
      }
    }

    log_file.flush();
    if (opt.debug)
      log_file << "status: " << opt.cif << "&" << opt.hkl << "&" << opt.xyz_file << endl;
    err_checkf(opt.xyz_file != "", "No xyz specified", log_file);
    err_checkf(exists(opt.xyz_file), "xyz doesn't exist", log_file);
    wavy.push_back(WFN(0));
    wavy[0].read_known_wavefunction_format(opt.xyz_file, log_file, opt.debug);

    if (opt.electron_diffraction && opt.debug)
      log_file << "Making Electron diffraction scattering factors, be carefull what you are doing!" << endl;
    if (opt.debug)
      log_file << "Entering Structure Factor Calculation!" << endl;
    err_checkf(thakkar_sfac(opt, log_file, wavy[0]), "Error during SF Calculation!", log_file);
    return 0;
  }
  //This one has conversion to fchk and calculation of one single tsc file
  if (opt.hkl != "" && opt.wfn != "") {
    log_file << NoSpherA2_message(opt.no_date);
    //Lets print what was the command line, for debugging
    if (opt.debug)
      for (int i = 0; i < argc; i++)
        log_file << argv[i] << endl;
    if (opt.debug) {
      log_file << "status:" << opt.wfn << "&" << opt.fchk << "&" << opt.basis_set << "&" << opt.basis_set_path << "&" << opt.cif << "&" << opt.hkl << "&" << opt.groups.size();
      if (opt.groups.size() != 0) log_file << "&" << opt.groups[0].size();
      log_file << endl;
    }
    wavy.push_back(WFN(0));
    log_file << "Reading: " << setw(44) << opt.wfn << flush;
    wavy[0].read_known_wavefunction_format(opt.wfn, log_file, opt.debug);
    wavy[0].set_method(opt.method);
    if (opt.ECP) {
      wavy[0].set_has_ECPs(true);
    }
    log_file << " done!" << endl << "Number of atoms in Wavefunction file: " << wavy[0].get_ncen() << " Number of MOs: " << wavy[0].get_nmo() << endl;

    if (opt.basis_set != "" || opt.fchk != "") {
      //Make a fchk out of the wfn/wfx file
      join_path(opt.basis_set_path, opt.basis_set);
      err_checkf(exists(opt.basis_set_path), "Basis set file does not exist!", log_file);
      wavy[0].set_basis_set_name(opt.basis_set_path);

      string outputname;
      if (opt.fchk != "")
        outputname = opt.fchk;
      else {
        outputname = wavy[0].get_path();
        if (opt.debug) log_file << "Loaded path..." << endl;
        size_t where;
        if (wavy[0].get_origin() == 2) where = outputname.find("wfn");
        else if (wavy[0].get_origin() == 4) where = outputname.find("ffn");
        else if (wavy[0].get_origin() == 4) where = outputname.find(".wfx");
        if (where >= outputname.length() && where != string::npos) err_checkf(false, "Cannot make output file name!", log_file);
        else outputname.erase(where, 3);
      }
      wavy[0].assign_charge(wavy[0].calculate_charge());
      if (opt.mult == 0) err_checkf(wavy[0].guess_multiplicity(log_file), "Error guessing multiplicity", log_file);
      else wavy[0].assign_multi(opt.mult);
      free_fchk(log_file, outputname, "", wavy[0], opt.debug, true);
    }
    if (opt.cif != "" || opt.hkl != "") {
      // Calculate tsc file from given files
      if (opt.debug)
        log_file << "Entering Structure Factor Calculation!" << endl;
      if (opt.electron_diffraction)
        log_file << "Making Electron diffraction scattering factors, be carefull what you are doing!" << endl;
      if (wavy[0].get_origin() != 7)
        err_checkf(calculate_structure_factors_HF(
          opt,
          wavy[0],
          log_file), "Error during SF Calcualtion", log_file);
      else
        err_checkf(thakkar_sfac(
          opt,
          log_file,
          wavy[0]), "Error during SF Calcualtion", log_file);
    }
    return 0;
  }
  if (opt.calc) {
    ofstream log2("NoSpherA2_cube.log", ios::out);
    auto coutbuf = std::cout.rdbuf(log2.rdbuf()); //save and redirect
    if (opt.debug)
      for (int i = 0; i < argc; i++)
        log2 << argv[i] << endl;
    log2 << NoSpherA2_message(opt.no_date);
    log2.flush();

    err_checkf(opt.wfn != "", "Error, no wfn file specified!", log2);
    wavy.push_back(WFN(0));
    wavy[0].read_known_wavefunction_format(opt.wfn, log2, opt.debug);
    if (opt.debug)
      log2 << "Starting calculation of properties" << endl;
    if (opt.all_mos)
      for (int mo = 0; mo < wavy[0].get_nmo(); mo++)
        opt.MOs.push_back(mo);
    if (opt.debug)
      log2 << "Size of MOs: " << opt.MOs.size() << endl;

    vector < vector < double > > cell_matrix;
    cell_matrix.resize(3);
    for (int i = 0; i < 3; i++)
      cell_matrix[i].resize(3);
    if (opt.debug)
      log2 << opt.cif << " " << opt.resolution << " " << opt.radius << endl;
    readxyzMinMax_fromCIF(opt.cif, opt.MinMax, opt.NbSteps, cell_matrix, opt.resolution, log2, opt.debug);
    if (opt.debug) {
      log2 << "Resolution: " << opt.resolution << endl;
      log2 << "MinMax:" << endl;
      for (int i = 0; i < 6; i++)
        log2 << setw(14) << scientific << opt.MinMax[i];
      log2 << endl;
      log2 << "Steps:" << endl;
      for (int i = 0; i < 3; i++)
        log2 << setw(14) << scientific << opt.NbSteps[i];
      log2 << endl;
      log2 << "Cell Matrix:" << endl;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
          log2 << setw(14) << scientific << cell_matrix[i][j];
        log2 << endl;
      }
    }
    cube Rho  (opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy[0].get_ncen(), true);
    cube RDG(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy[0].get_ncen(), opt.rdg);
    cube Elf(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy[0].get_ncen(), opt.elf);
    cube Eli(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy[0].get_ncen(), opt.eli);
    cube Lap(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy[0].get_ncen(), opt.lap);
    cube ESP(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy[0].get_ncen(), opt.esp);
    cube MO(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy[0].get_ncen(), true);
    cube HDEF(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy[0].get_ncen(), opt.hdef);
    cube DEF(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy[0].get_ncen(), opt.def);
    cube Hirsh(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy[0].get_ncen(), opt.hirsh);

    Rho.give_parent_wfn(wavy[0]);
    RDG.give_parent_wfn(wavy[0]);
    Elf.give_parent_wfn(wavy[0]);
    Eli.give_parent_wfn(wavy[0]);
    Lap.give_parent_wfn(wavy[0]);
    ESP.give_parent_wfn(wavy[0]);
    MO.give_parent_wfn(wavy[0]);
    HDEF.give_parent_wfn(wavy[0]);
    DEF.give_parent_wfn(wavy[0]);
    Hirsh.give_parent_wfn(wavy[0]);

    for (int i = 0; i < 3; i++) {
      Rho.set_origin  (i, opt.MinMax[i]);
      RDG.set_origin  (i, opt.MinMax[i]);
      Elf.set_origin  (i, opt.MinMax[i]);
      Eli.set_origin  (i, opt.MinMax[i]);
      Lap.set_origin  (i, opt.MinMax[i]);
      ESP.set_origin  (i, opt.MinMax[i]);
      MO.set_origin   (i, opt.MinMax[i]);
      HDEF.set_origin (i, opt.MinMax[i]);
      DEF.set_origin  (i, opt.MinMax[i]);
      Hirsh.set_origin(i, opt.MinMax[i]);
      for (int j = 0; j < 3; j++) {
        Rho.set_vector(i, j, cell_matrix[i][j]);
        RDG.set_vector(i, j, cell_matrix[i][j]);
        Elf.set_vector(i, j, cell_matrix[i][j]);
        Eli.set_vector(i, j, cell_matrix[i][j]);
        Lap.set_vector(i, j, cell_matrix[i][j]);
        ESP.set_vector(i, j, cell_matrix[i][j]);
        MO.set_vector(i, j, cell_matrix[i][j]);
        HDEF.set_vector(i, j, cell_matrix[i][j]);
        DEF.set_vector(i, j, cell_matrix[i][j]);
        Hirsh.set_vector(i, j, cell_matrix[i][j]);
      }
    }
    if (opt.debug)
      log2 << "Origins etc are set up" << endl;
    Rho.set_comment1("Calculated density using NoSpherA2");
    RDG.set_comment1("Calculated reduced density gradient using NoSpherA2");
    Elf.set_comment1("Calculated electron localization function using NoSpherA2");
    Eli.set_comment1("Calculated same-spin electron localizability indicator using NoSpherA2");
    Lap.set_comment1("Calculated laplacian of electron density using NoSpherA2");
    ESP.set_comment1("Calculated electrostatic potential using NoSpherA2");
    MO.set_comment1("Calcualted MO values using NoSpherA2");
    HDEF.set_comment1("Calculated Atomic Hirshfeld deformation density values using NoSpherA2");
    DEF.set_comment1("Calculated static deformation density values using NoSpherA2");
    Hirsh.set_comment1("Calculated Hirshfeld atom density values using NoSpherA2");
    Rho.set_comment2("from " + wavy[0].get_path());
    RDG.set_comment2("from " + wavy[0].get_path());
    Elf.set_comment2("from " + wavy[0].get_path());
    Eli.set_comment2("from " + wavy[0].get_path());
    Lap.set_comment2("from " + wavy[0].get_path());
    ESP.set_comment2("from " + wavy[0].get_path());
    MO.set_comment2("from" + wavy[0].get_path());
    HDEF.set_comment2("from" + wavy[0].get_path());
    DEF.set_comment2("from" + wavy[0].get_path());
    Hirsh.set_comment2("from" + wavy[0].get_path());
    Rho.path = get_basename_without_ending(wavy[0].get_path()) + "_rho.cube";
    RDG.path = get_basename_without_ending(wavy[0].get_path()) + "_rdg.cube";
    Elf.path = get_basename_without_ending(wavy[0].get_path()) + "_elf.cube";
    Eli.path = get_basename_without_ending(wavy[0].get_path()) + "_eli.cube";
    Lap.path = get_basename_without_ending(wavy[0].get_path()) + "_lap.cube";
    ESP.path = get_basename_without_ending(wavy[0].get_path()) + "_esp.cube";
    DEF.path = get_basename_without_ending(wavy[0].get_path()) + "_def.cube";
    Hirsh.path = get_basename_without_ending(wavy[0].get_path()) + "_hirsh.cube";

    if (opt.debug) {
      log2 << "Status: " << opt.hdef << opt.def << opt.hirsh << opt.lap << opt.eli << opt.elf << opt.rdg << opt.esp << endl;
      log2 << "Everything is set up; starting calculation..." << endl;
    }

    log2 << "Calculating for " << fixed << setprecision(0) << opt.NbSteps[0] * opt.NbSteps[1] * opt.NbSteps[2] << " Gridpoints." << endl;

    if (opt.MOs.size() != 0)
      for (int i = 0; i < opt.MOs.size(); i++) {
        log2 << "Calcualting MO: " << opt.MOs[i] << endl;
        MO.set_zero();
        MO.path = get_basename_without_ending(wavy[0].get_path()) + "_MO_" + to_string(opt.MOs[i]) + ".cube";
        Calc_MO(MO, opt.MOs[i], wavy[0], opt.threads, opt.radius, log2);
        MO.write_file(true);
      }

    if (opt.hdef || opt.def || opt.hirsh) {
      log2 << "Calcualting Rho...";
      Calc_Rho(Rho, wavy[0], opt.threads, opt.radius, log2);
      log2 << " ...done!" << endl;
      cube temp(opt.NbSteps[0], opt.NbSteps[1], opt.NbSteps[2], wavy[0].get_ncen(), opt.hdef || opt.hirsh);
      for (int i = 0; i < 3; i++) {
        temp.set_origin(i, opt.MinMax[i]);
        for (int j = 0; j < 3; j++)
          temp.set_vector(i, j, cell_matrix[i][j]);
      }
      if (opt.hdef || opt.hirsh) {
        log2 << "Calcualting spherical Rho...";
        Calc_Spherical_Dens(temp, wavy[0], opt.threads, opt.radius, log2);
        log2 << " ...done!" << endl;
      }

      if (opt.def) {
        log2 << "Calculating static deformation density...";
        if (opt.hdef)
          Calc_Static_Def(DEF, Rho, temp, wavy[0], opt.threads, opt.radius, log2);
        else
          Calc_Static_Def(DEF, Rho, wavy[0], opt.threads, opt.radius, log2);
        log2 << " ...done!" << endl;
      }

      if (opt.hdef) {
        for (int a = 0; a < wavy[0].get_ncen(); a++) {
          log2 << "Calcualting Hirshfeld deformation density for atom: " << a << endl;
          HDEF.path = get_basename_without_ending(wavy[0].get_path()) + "_HDEF_" + to_string(a) + ".cube";
          Calc_Hirshfeld(HDEF, Rho, temp, wavy[0], opt.threads, opt.radius, a, log2);
          HDEF.write_file(true);
          HDEF.set_zero();
        }
      }

      if (opt.hirsh)
      {
        log2 << "Calcualting Hirshfeld density for atom: " << opt.hirsh_number << endl;
        Calc_Hirshfeld_atom(Hirsh, Rho, temp, wavy[0], opt.threads, opt.radius, opt.hirsh_number, log2);
        log2 << "..done!" << endl;
      }
    }

    if (opt.lap || opt.eli || opt.elf || opt.rdg || opt.esp)
      Calc_Prop(Rho, RDG, Elf, Eli, Lap, ESP, wavy[0], opt.threads, opt.radius, log2, opt.no_date);

    log2 << "Writing cubes to Disk..." << flush;
    if (opt.rdg) {
      Rho.path = get_basename_without_ending(wavy[0].get_path()) + "_signed_rho.cube";
      Rho.write_file(true);
      Rho.path = get_basename_without_ending(wavy[0].get_path()) + "_rho.cube";
      Rho.write_file(true, true);
    }
    else if (opt.lap || opt.eli || opt.elf || opt.esp) Rho.write_file(true);
    if (opt.rdg) RDG.write_file(true);
    if (opt.lap) Lap.write_file(true);
    if (opt.elf) Elf.write_file(true);
    if (opt.eli) Eli.write_file(true);
    if (opt.def) {
      DEF.write_file(true);
      Rho.write_file(true);
    }
    if (opt.hirsh) Hirsh.write_file(true);

    log2 << " done!" << endl;

    if (opt.esp) {
      log2 << "Calculating ESP..." << flush;
      WFN temp = wavy[0];
      temp.delete_unoccupied_MOs();
      Calc_ESP(ESP, temp, opt.threads, opt.radius, opt.no_date, log2);
      log2 << "Writing cube to Disk..." << flush;
      ESP.write_file(true);
      log2 << "  done!" << endl;
    }
    return 0;
  }
  if (opt.density_test_cube) {
    log_file << NoSpherA2_message(opt.no_date);
    log_file.flush();
    wavy.resize(10);
    //ScF2+ test file against ORCA calcualted cubes
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
    //Rho.write_file(wavy[0], true);
    const double test_molden = Rho.sum();
    log_file << "Number of electrons in the cube: " << setprecision(4) << fixed << test_molden << endl;
    cube Rho2("test.eldens.cube", true, wavy[0], log_file);
    const double test_ref = Rho2.sum();
    log_file << "Number of electrons in the reference cube: " << setprecision(4) << fixed << test_ref << endl;
    for (int i = 0; i < 1; i++) {
      cube MO(100, 100, 100, wavy[0].get_ncen(), true);
      MO.set_origin(0, -7), MO.set_origin(1, -7), MO.set_origin(2, -7);
      MO.set_vector(0, 0, 0.141414);
      MO.set_vector(1, 1, 0.141414);
      MO.set_vector(2, 2, 0.141414);
      MO.path = get_basename_without_ending(wavy[0].get_path()) + "_MO_" + to_string(i) + ".cube";
      log_file << "Calcualting MO " + to_string(i) + "...";
      Calc_MO(MO, i, wavy[0], opt.threads, 4.0, log_file);
      log_file << " ...done!" << endl;
      //MO.write_file(wavy[0], true);
      string name("test.mo" + to_string(i) + "a.cube");
      cube MO2(name, true, wavy[0], log_file);
      log_file << "sum in the cube: " << setprecision(4) << fixed << MO.sum() << endl;
      log_file << "sum in the reference cube: " << setprecision(4) << fixed << MO2.sum() << endl;
    }

    //F- ion calculations
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
    //Rho_2.write_file(wavy[1], true);
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
    //Rho_3.write_file(wavy[1], true);
    const double F_ref = Rho_3.sum();
    log_file << "Number of electrons in the reference cube: " << setprecision(4) << fixed << F_ref << endl;

    //Ce conevrsion for test of f type
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
    //Rho_4.write_file(wavy[2], true);
    const double Ce_molden = Rho_4.sum();
    log_file << "Number of electrons in the cube: " << setprecision(4) << fixed << Ce_molden << endl;

    Rho_4.set_zero();
    wavy[5].read_known_wavefunction_format("Ce_full.wfn", log_file);
    Calc_Rho(Rho_4, wavy[5], opt.threads, 7.0, log_file);
    Rho_4.path = get_basename_without_ending(wavy[5].get_path()) + "_reference_rho.cube";
    //Rho_4.write_file(wavy[5], true);
    const double Ce_ref = Rho_4.sum();
    log_file << "Number of electrons in the ref cube: " << setprecision(4) << fixed << Ce_ref << endl;

    //Sc conversion
    log_file << "====================Sc Test===============================" << endl;
    wavy[3].read_known_wavefunction_format("Sc_full.molden", log_file);
    wavy[3].write_wfn("Sc_conv.wfn", false, true);

    //Lu g-type test
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
    //Rho4.write_file(wavy[6], true);
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
    //Rho_5.write_file(wavy[7], true);
    const double Lu_wfn = Rho_5.sum();
    Rho_5 -= Rho4;
    Rho_5.path = get_basename_without_ending(wavy[7].get_path()) + "_diff_rho.cube";
    //Rho_5.write_file(wavy[7], true);
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
    //Rho4.write_file(wavy[6], true);
    const double Lu_def2 = Rho5.sum();
    log_file << "Number of electrons in cube: " << setprecision(4) << fixed << Lu_def2 << endl;

    err_checkf(abs(test_molden - test_ref) < 0.1, "Difference in test too big!", log_file);
    err_checkf(abs(F_molden - F_ref) < 0.1, "Difference in F too big!", log_file);
    err_checkf(abs(Ce_molden - Ce_ref) < 0.1, "Difference in Ce too big!", log_file);
    err_checkf(abs(Lu_molden - Lu_wfn) < 0.1, "Difference in Lu too big!", log_file);
    log_file << "All tests successfull!" << endl;
    return 0;
  }
  if (opt.gbw2wfn) {
    err_checkf(opt.wfn != "", "No Wavefunction given!", log_file);
    wavy.push_back(WFN(9));
    wavy[0].read_known_wavefunction_format(opt.wfn, log_file);
    wavy[0].write_wfn("converted.wfn", false, true);
    return 0;
  }
  cout << NoSpherA2_message(opt.no_date) << "Did not understand the task to perform!\n" << help_message() << endl;
  return 0;
}