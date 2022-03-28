
#include "convenience.h"
#include "fchk.h"
#include "cube.h"
#include "basis_set.h"
#include "structure_factors.h"
#include "properties.h"
#include "cell.h"

using namespace std;
bool debug_main = false;
bool debug_all = false;
bool expert = false;

int main(int argc, char** argv)
{
  if (debug_main)
    cout << "argc:" << argc << endl;
  vector<WFN> wavy;
  string wavename;
  string gaussian_path;
  string turbomole_path;
  string basis_set_path;
  int ncpus = 0;
  double mem = 0.0;
  int mult = 0;
  string wfn("");
  string fchk("");
  string basis_set("");
  string hkl("");
  string cif("");
  string method("rhf");
  string xyz_file("");
  string temp;
  string fract_name("");
  vector<string> combined_tsc_calc_files;
  vector<string> combined_tsc_calc_cifs;
  int accuracy = 2;
  int threads = -1;
  int pbc = 0;
  double resolution = 0.1;
  double radius = 2.0;
  bool becke = false;
  bool electron_diffraction = false;
  bool calc = false;
  bool eli = false;
  bool esp = false;
  bool elf = false;
  bool lap = false;
  bool rdg = false;
  bool hdef = false;
  bool def = false;
  bool fract = false;
  bool hirsh = false;
  bool Olex2_1_3_switch = false;
  bool iam_switch = false;
  bool read_k_pts = false;
  bool save_k_pts = false;
  bool combined_tsc_calc = false;
  bool binary_tsc = true;
  bool cif_based_combined_tsc_calc = false;
  int hirsh_number = 0;
  double MinMax[6];
  double NbSteps[3];
  vector <int> MOs;
  vector < vector <int> > groups;
  vector < vector <double> > twin_law;
  vector < vector <int> > combined_tsc_groups;
  bool all_mos = false;
  bool test = false;
  groups.resize(1);
  for (int i = 0; i < argc; i++) {
    temp = argv[i];
    if (temp.find("-") > 0) continue;
    if (temp.find("-wfn") < 1)
      wfn = argv[i + 1];
    else if (temp.find("-fchk") < 1)
      fchk = argv[i + 1];
    else if (temp.find("-b") < 1)
      basis_set = argv[i + 1];
    else if (temp.find("-d") < 1)
      basis_set_path = argv[i + 1];
    else if (temp.find("-hkl") < 1)
      hkl = argv[i + 1];
    else if (temp.find("-cif") < 1)
      cif = argv[i + 1];
    else if (temp.find("-acc") < 1)
      accuracy = stoi(argv[i + 1]);
    else if (temp.find("-mult") < 1)
      mult = stoi(argv[i + 1]);
    else if (temp.find("-method") < 1)
      method = argv[i + 1];
    else if (temp.find("-cpus") < 1)
      threads = stoi(argv[i + 1]);
    else if (temp.find("-pbc") < 1)
      pbc = stoi(argv[i + 1]);
    else if (temp.find("-ED") < 1)
      electron_diffraction = true;
    else if (temp.find("-Olex2_1_3") < 1)
      Olex2_1_3_switch = true;
    else if (temp.find("-v2") < 1)
      cout << "Turning on verbose mode 2!" << endl, debug_all = debug_main = true;
    else if (temp.find("-v") < 1)
      cout << "Turning on verbose mode!" << endl, debug_main = true;
    else if (temp.find("-test") < 1)
      cout << "Running in test mode!" << endl, test = true;
    else if (temp.find("-eli") < 1)
      calc = eli = true;
    else if (temp.find("-elf") < 1)
      calc = elf = true;
    else if (temp.find("-lap") < 1)
      calc = lap = true;
    else if (temp.find("-esp") < 1)
      calc = esp = true;
    else if (temp.find("-rdg") < 1)
      calc = rdg = true;
    else if (temp.find("-hirsh") < 1)
      calc = hirsh = true, hirsh_number = stoi(argv[i + 1]);
    else if (temp.find("-resolution") < 1)
      resolution = stod(argv[i + 1]);
    else if (temp.find("-radius") < 1)
      radius = stod(argv[i + 1]);
    else if (temp.find("-MO") < 1) {
      if (string(argv[i + 1]) != "all")
        MOs.push_back(stoi(argv[i + 1]));
      else
        all_mos = true;
      calc = true;
    }
    else if (temp.find("-fractal") < 1)
      fract = true, fract_name = argv[i + 1];
    else if (temp.find("-HDEF") < 1)
      hdef = calc = true;
    else if (temp.find("-def") < 1)
      def = calc = true;
    else if (temp.find("-skpts") < 1)
      save_k_pts = true;
    else if (temp.find("-rkpts") < 1)
      read_k_pts = true;
    else if (temp.find("-group") < 1) {
      int n = 1;
      while (i + n < argc && string(argv[i + n]).find("-") == string::npos) {
        int group;
        if (argv[i + 1][0] == '+')
          group = -stoi(argv[i + n]);
        else
          group = stoi(argv[i + n]);
        groups[0].push_back(group), n++;
      }
      i += n;
    }
    else if (temp.find("-twin") < 1) {
      twin_law.resize(twin_law.size() + 1);
      twin_law[twin_law.size() - 1].resize(9);
      for (int twl = 0; twl < 9; twl++)
        twin_law[twin_law.size() - 1][twl] = stod(argv[i + 1 + twl]);
      if (debug_main) {
        cout << "twin_law: ";
        for (int twl = 0; twl < 9; twl++)
          cout << setw(7) << setprecision(2) << twin_law[twin_law.size() - 1][twl];
        cout << endl;
      }
      i += 9;
    }
    else if (temp.find("-IAM") != string::npos)
      iam_switch = true;
    else if (temp.find("-xyz") != string::npos)
      xyz_file = argv[i + 1];
    else if (temp.find("-merge") != string::npos) {
      vector<string> filenames;
      int n = 1;
      while (i + n < argc && string(argv[i + n]).find("-") > 0) {
        filenames.push_back(argv[i + n]);
        n++;
      }
      merge_tscs("combine", filenames, debug_all);
      return 0;
    }
    else if (temp.find("-merge_nocheck") != string::npos) {
      vector<string> filenames;
      int n = 1;
      while (i + n < argc && string(argv[i + n]).find("-") > 0) {
        filenames.push_back(argv[i + n]);
        n++;
      }
      merge_tscs_without_checks("combine", filenames, debug_all);
      return 0;
    }
    else if (temp.find("-mtc") != string::npos) {
      combined_tsc_calc = true;
      int n = 1;
      string delimiter = ",";
      groups.pop_back();
      while (i + n < argc && string(argv[i + n]).find("-") > 0) {
        combined_tsc_calc_files.push_back(argv[i + n]);
        n++;
        const string temp = argv[i + n];
        groups.push_back(split_string<int>(temp, delimiter));
        n++;
      }
    }
    else if (temp.find("-cmtc") != string::npos) {
      cif_based_combined_tsc_calc = true;
      int n = 1;
      string delimiter = ",";
      groups.pop_back();
      while (i + n < argc && string(argv[i + n]).find("-") > 0) {
        combined_tsc_calc_files.push_back(argv[i + n]);
        n++;
        combined_tsc_calc_cifs.push_back(argv[i + n]);
        n++;
        const string temp = argv[i + n];
        groups.push_back(split_string<int>(temp, delimiter));
        n++;
      }
    }
  }
  if (threads != -1) {
    omp_set_num_threads(threads);
    omp_set_dynamic(0);
  }
  if (argc > 1) {
    if (string(argv[1]).find("--help") != string::npos) {
      cout << NoSpherA2_message() << help_message() << endl;
      return 0;
    }
    if (string(argv[1]).find("--h") != string::npos) {
      cout << NoSpherA2_message() << help_message() << endl;
      return 0;
    }
    if (string(argv[1]).find("-help") != string::npos) {
      cout << NoSpherA2_message() << help_message() << endl;
      return 0;
    }
  }
  //Lets print what was the command line, for debugging
  if (debug_main) {
    ofstream log("NoSpherA2.log", ios::out);
    for (int i = 0; i < argc; i++) {
      //cout << argv[i] << endl;
      log << argv[i] << endl;
    }
    log.close();
  }
  if (fract) {
    cube residual(fract_name, true, debug_all);
    residual.fractal_dimension(0.01);
    return 0;
  }
  if (combined_tsc_calc) {
    ofstream log_file("NoSpherA2.log", ios::out);
    log_file << NoSpherA2_message();
    log_file.flush();
    err_checkf(hkl != "", "No hkl specified", log_file);
    err_checkf(exists(hkl), "hkl doesn't exist", log_file);
    err_checkf(cif != "", "No cif specified", log_file);
    err_checkf(exists(cif), "CIF doesn't exist", log_file);
    //Make sure we have more than 2 files...
    //err_checkf(combined_tsc_calc_files.size() > 1, "Need at least 2 wfn files", log_file);
    //First make sure all files exist
    for (int i = 0; i < combined_tsc_calc_files.size(); i++)
      err_checkf(exists(combined_tsc_calc_files[i]), "Specified file for combined calculation doesn't exist! " + combined_tsc_calc_files[i], log_file);

    wavy.resize(combined_tsc_calc_files.size());
    for (int i = 0; i < combined_tsc_calc_files.size(); i++)
      wavy[i].read_known_wavefunction_format(combined_tsc_calc_files[i]);

    vector<string> known_scatterer;
    tsc_block result;
    for (int i = 0; i < combined_tsc_calc_files.size(); i++) {
      known_scatterer = result.get_scatterers();
      if (wavy[i].get_origin() != 7) {
        result.append(calculate_structure_factors_MTC(
          hkl,
          cif,
          wavy[i],
          debug_main,
          accuracy,
          log_file,
          groups[i],
          twin_law,
          known_scatterer,
          threads,
          electron_diffraction,
          0,
          i == 0,
          i != 0
        ), log_file);
      }
      else {
        result.append(MTC_thakkar_sfac(
          hkl,
          cif,
          debug_main,
          log_file,
          groups[i],
          twin_law,
          known_scatterer,
          wavy[i],
          threads,
          electron_diffraction,
          i == 0,
          i != 0
        ), log_file);
      }
    }

    if (debug_main) {
      for (int i = 0; i < combined_tsc_calc_files.size(); i++) {
        log_file << combined_tsc_calc_files[i] << " ";
        for (int j = 0; j < groups[i].size(); j++) log_file << groups[i][j] << " ";
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
    if (binary_tsc)
      result.write_tscb_file();
    result.write_tsc_file(cif);
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
  if (cif_based_combined_tsc_calc) {
    ofstream log_file("NoSpherA2.log", ios::out);
    log_file << NoSpherA2_message();
    log_file.flush();
    err_checkf(hkl != "", "No hkl specified", log_file);
    err_checkf(exists(hkl), "hkl doesn't exist", log_file);
    //Make sure we have more than 2 files...
    //err_checkf(combined_tsc_calc_files.size() > 1, "Need at least 2 wfn files", log_file);
    //First make sure all files exist
    for (int i = 0; i < combined_tsc_calc_files.size(); i++)
      err_checkf(exists(combined_tsc_calc_files[i]), "Specified file for combined calculation doesn't exist! " + combined_tsc_calc_files[i], log_file);
    for (int i = 0; i < combined_tsc_calc_cifs.size(); i++)
      err_checkf(exists(combined_tsc_calc_cifs[i]), "Specified file for combined calculation doesn't exist! " + combined_tsc_calc_cifs[i], log_file);

    wavy.resize(combined_tsc_calc_files.size());
    for (int i = 0; i < combined_tsc_calc_files.size(); i++) {
      log_file << "Reading: " << setw(44) << combined_tsc_calc_files[i] << flush;
      wavy[i].read_known_wavefunction_format(combined_tsc_calc_files[i]);
      log_file << " done!" << endl << "Number of atoms in Wavefunction file: " << wavy[i].get_ncen() << " Number of MOs: " << wavy[i].get_nmo() << endl;
    }

    vector<string> known_scatterer;
    tsc_block result;
    for (int i = 0; i < combined_tsc_calc_files.size(); i++) {
      known_scatterer = result.get_scatterers();
      if (wavy[i].get_origin() != 7) {
        result.append(calculate_structure_factors_MTC(
          hkl,
          combined_tsc_calc_cifs[i],
          wavy[i],
          debug_main,
          accuracy,
          log_file,
          groups[i],
          twin_law,
          known_scatterer,
          threads,
          electron_diffraction,
          0,
          i == 0,
          i != 0
        ), log_file);
      }
      else {
        result.append(MTC_thakkar_sfac(
          hkl,
          combined_tsc_calc_cifs[i],
          debug_main,
          log_file,
          groups[i],
          twin_law,
          known_scatterer,
          wavy[i],
          threads,
          electron_diffraction,
          i == 0,
          i != 0
        ), log_file);
      }
    }

    if (debug_main) {
      for (int i = 1; i < combined_tsc_calc_files.size(); i++) {
        log_file << combined_tsc_calc_files[i] << " ";
        for (int j = 0; j < groups[i].size(); j++) log_file << groups[i][j] << " ";
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
    if (binary_tsc)
      result.write_tscb_file();
    result.write_tsc_file(cif);
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
  if (iam_switch) {
    ofstream log_file("NoSpherA2.log", ios::out);
    if (debug_main)
      for (int i = 0; i < argc; i++)
        log_file << argv[i] << endl;

    log_file << NoSpherA2_message();
    log_file.flush();
    if (argc > 1)
      if (string(argv[1]).find("--help") != string::npos) {
        log_file << help_message() << endl;
        log_file.close();
        return 0;
      }

    log_file.flush();
    if (debug_main)
      log_file << "status: " << cif << "&" << hkl << "&" << xyz_file << endl;
    err_checkf(exists(cif), "CIF doesn't exist", log_file);
    err_checkf(xyz_file != "", "No xyz specified", log_file);
    err_checkf(exists(xyz_file), "xyz doesn't exist", log_file);
    wavy.push_back(WFN(0));
    wavy[0].read_known_wavefunction_format(xyz_file, log_file, debug_main);

    if (electron_diffraction && debug_main)
      log_file << "Making Electron diffraction scattering factors, be carefull what you are doing!" << endl;
    if (cif != "" || hkl != "") {
      if (debug_main)
        log_file << "Entering Structure Factor Calculation!" << endl;
      err_checkf(thakkar_sfac(hkl, cif, debug_main, log_file, groups[0], twin_law, wavy[0], threads, electron_diffraction), "Error during SF Calculation!", log_file);
    }
    return 0;
  }
  else if (hkl != "" || basis_set != "" || fchk != "") {
    ofstream log_file("NoSpherA2.log", ios::out);

    log_file << NoSpherA2_message();
    //Lets print what was the command line, for debugging
    if (debug_main)
      for (int i = 0; i < argc; i++)
        log_file << argv[i] << endl;
    err_checkf(argc >= 4, "Not enough arguments given, at least provide -wfn <FILENAME>.wfn/.wfx -b <basis_set>", log_file);
    if (debug_main) {
      log_file << "status:" << wfn << "&" << fchk << "&" << basis_set << "&" << basis_set_path << "&" << cif << "&" << hkl << "&" << groups.size();
      if (groups.size() != 0) log_file << "&" << groups[0].size();
      log_file << endl;
    }
    err_checkf(wfn != "", "No wfn specified", log_file);
    wavy.push_back(WFN(0));
    log_file << "Reading: " << setw(44) << wfn << flush;
    wavy[0].read_known_wavefunction_format(wfn, log_file, debug_main);
    wavy[0].set_method(method);
    log_file << " done!" << endl << "Number of atoms in Wavefunction file: " << wavy[0].get_ncen() << " Number of MOs: " << wavy[0].get_nmo() << endl;
    if (electron_diffraction)
      log_file << "Making Electron diffraction scattering factors, be carefull what you are doing!" << endl;

    if (basis_set != "" || fchk != "") {
      //Make and fchk out of the wfn/wfx file
      join_path(basis_set_path, basis_set);
      err_checkf(exists(basis_set_path), "Basis set file does not exist!", log_file);
      wavy[0].set_basis_set_name(basis_set_path);

      string outputname;
      if (fchk != "")
        outputname = fchk;
      else {
        outputname = wavy[0].get_path();
        if (debug_main) log_file << "Loaded path..." << endl;
        size_t where;
        if (wavy[0].get_origin() == 2) where = outputname.find("wfn");
        else if (wavy[0].get_origin() == 4) where = outputname.find("ffn");
        else if (wavy[0].get_origin() == 4) where = outputname.find(".wfx");
        if (where >= outputname.length() && where != string::npos) err_checkf(false, "Cannot make output file name!", log_file);
        else outputname.erase(where, 3);
      }
      wavy[0].assign_charge(wavy[0].calculate_charge());
      if (mult == 0) err_checkf(wavy[0].guess_multiplicity(log_file, expert), "Error guessing multiplicity", log_file);
      else wavy[0].assign_multi(mult);
      free_fchk(log_file, outputname, "", wavy[0], debug_main, true);
    }
    if (cif != "" || hkl != "") {
      // Calculate tsc fiel from given files
      err_checkf(exists(hkl), "Hkl file doesn't exist!", log_file);
      err_checkf(exists(cif), "CIF doesn't exist!", log_file);
      if (debug_main)
        log_file << "Entering Structure Factor Calculation!" << endl;
      if (wavy[0].get_origin() != 7)
        err_checkf(calculate_structure_factors_HF(
          hkl,
          cif,
          wavy[0],
          debug_main,
          accuracy,
          log_file,
          groups[0],
          twin_law,
          threads,
          electron_diffraction,
          pbc,
          Olex2_1_3_switch,
          save_k_pts,
          read_k_pts), "Error during SF Calcualtion", log_file);
      else
        err_checkf(thakkar_sfac(
          hkl,
          cif,
          debug_main,
          log_file,
          groups[0],
          twin_law,
          wavy[0],
          threads,
          electron_diffraction,
          save_k_pts,
          read_k_pts), "Error during SF Calcualtion", log_file);
    }
    return 0;
  }
  if (calc) {
    ofstream log2("NoSpherA2_cube.log", ios::out);
    if (debug_main)
      for (int i = 0; i < argc; i++)
        log2 << argv[i] << endl;
    log2 << NoSpherA2_message();
    log2.flush();

    err_checkf(wfn != "", "Error, no wfn file specified!", log2);
    wavy.push_back(WFN(0));
    wavy[0].read_known_wavefunction_format(wfn, log2, debug_main);
    if (debug_main)
      log2 << "Starting calculation of properties" << endl;
    if (all_mos)
      for (int mo = 0; mo < wavy[0].get_nmo(); mo++)
        MOs.push_back(mo);
    if (debug_main)
      log2 << "Size of MOs: " << MOs.size() << endl;

    vector < vector < double > > cell_matrix;
    cell_matrix.resize(3);
    for (int i = 0; i < 3; i++)
      cell_matrix[i].resize(3);
    if (debug_main)
      log2 << cif << " " << resolution << " " << radius << endl;
    readxyzMinMax_fromCIF(cif, MinMax, NbSteps, cell_matrix, resolution, log2, debug_main);
    if (debug_main) {
      log2 << "Resolution: " << resolution << endl;
      log2 << "MinMax:" << endl;
      for (int i = 0; i < 6; i++)
        log2 << setw(14) << scientific << MinMax[i];
      log2 << endl;
      log2 << "Steps:" << endl;
      for (int i = 0; i < 3; i++)
        log2 << setw(14) << scientific << NbSteps[i];
      log2 << endl;
      log2 << "Cell Matrix:" << endl;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
          log2 << setw(14) << scientific << cell_matrix[i][j];
        log2 << endl;
      }
    }
    cube Rho(NbSteps[0], NbSteps[1], NbSteps[2], wavy[0].get_ncen(), true);
    cube RDG(NbSteps[0], NbSteps[1], NbSteps[2], wavy[0].get_ncen(), rdg);
    cube Elf(NbSteps[0], NbSteps[1], NbSteps[2], wavy[0].get_ncen(), elf);
    cube Eli(NbSteps[0], NbSteps[1], NbSteps[2], wavy[0].get_ncen(), eli);
    cube Lap(NbSteps[0], NbSteps[1], NbSteps[2], wavy[0].get_ncen(), lap);
    cube ESP(NbSteps[0], NbSteps[1], NbSteps[2], wavy[0].get_ncen(), esp);
    cube MO(NbSteps[0], NbSteps[1], NbSteps[2], wavy[0].get_ncen(), true);
    cube HDEF(NbSteps[0], NbSteps[1], NbSteps[2], wavy[0].get_ncen(), hdef);
    cube DEF(NbSteps[0], NbSteps[1], NbSteps[2], wavy[0].get_ncen(), def);
    cube Hirsh(NbSteps[0], NbSteps[1], NbSteps[2], wavy[0].get_ncen(), hirsh);

    for (int i = 0; i < 3; i++) {
      Rho.set_origin(i, MinMax[i]);
      RDG.set_origin(i, MinMax[i]);
      Elf.set_origin(i, MinMax[i]);
      Eli.set_origin(i, MinMax[i]);
      Lap.set_origin(i, MinMax[i]);
      ESP.set_origin(i, MinMax[i]);
      MO.set_origin(i, MinMax[i]);
      HDEF.set_origin(i, MinMax[i]);
      DEF.set_origin(i, MinMax[i]);
      Hirsh.set_origin(i, MinMax[i]);
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
    if (debug_main)
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

    if (debug_main)
      log2 << "Everything is set up; starting calculation..." << endl;

    log2 << "Calculating for " << fixed << setprecision(0) << NbSteps[0] * NbSteps[1] * NbSteps[2] << " Gridpoints." << endl;

    if (MOs.size() != 0)
      for (int i = 0; i < MOs.size(); i++) {
        log2 << "Calcualting MO: " << MOs[i] << endl;
        MO.path = get_basename_without_ending(wavy[0].get_path()) + "_MO_" + to_string(MOs[i]) + ".cube";
        Calc_MO(MO, MOs[i], wavy[0], ncpus, radius, log2);
        MO.write_file(wavy[0], true);
      }

    if (hdef || def || hirsh) {
      log2 << "Calcualting Rho...";
      Calc_Rho(Rho, wavy[0], ncpus, radius, log2);
      log2 << " ...done!" << endl;
      cube temp(NbSteps[0], NbSteps[1], NbSteps[2], wavy[0].get_ncen(), hdef || hirsh);
      for (int i = 0; i < 3; i++) {
        temp.set_origin(i, MinMax[i]);
        for (int j = 0; j < 3; j++)
          temp.set_vector(i, j, cell_matrix[i][j]);
      }
      if (hdef || hirsh) {
        log2 << "Calcualting spherical Rho...";
        Calc_Spherical_Dens(temp, wavy[0], ncpus, radius, log2);
        log2 << " ...done!" << endl;
      }

      if (def) {
        log2 << "Calculating static deformation density...";
        if (hdef)
          Calc_Static_Def(DEF, Rho, temp, wavy[0], ncpus, radius, log2);
        else
          Calc_Static_Def(DEF, Rho, wavy[0], ncpus, radius, log2);
        log2 << " ...done!" << endl;
      }

      if (hdef)
        for (int a = 0; a < wavy[0].get_ncen(); a++) {
          log2 << "Calcualting Hirshfeld deformation density for atom: " << a << endl;
          HDEF.path = get_basename_without_ending(wavy[0].get_path()) + "_HDEF_" + to_string(a) + ".cube";
          Calc_Hirshfeld(HDEF, Rho, temp, wavy[0], ncpus, radius, a, log2);
          HDEF.write_file(wavy[0], true);
          HDEF.set_zero();
        }

      if (hirsh)
      {
        log2 << "Calcualting Hirshfeld density for atom: " << hirsh_number << endl;
        Calc_Hirshfeld_atom(Hirsh, Rho, temp, wavy[0], ncpus, radius, hirsh_number, log2);
        log2 << "..done!" << endl;
      }
    }

    if (lap || eli || elf || rdg || esp)
      Calc_Prop(Rho, RDG, Elf, Eli, Lap, ESP, wavy[0], ncpus, radius, log2, test);

    log2 << "Writing cubes to Disk..." << flush;
    if (rdg) {
      Rho.path = get_basename_without_ending(wavy[0].get_path()) + "_signed_rho.cube";
      Rho.write_file(wavy[0], true);
      Rho.path = get_basename_without_ending(wavy[0].get_path()) + "_rho.cube";
      Rho.write_file(wavy[0], true, true);
    }
    else if (lap || eli || elf || esp) Rho.write_file(wavy[0], true);
    if (rdg) RDG.write_file(wavy[0], true);
    if (lap) Lap.write_file(wavy[0], true);
    if (elf) Elf.write_file(wavy[0], true);
    if (eli) Eli.write_file(wavy[0], true);
    if (def) {
      DEF.write_file(wavy[0], true);
      Rho.write_file(wavy[0], true);
    }
    if (hirsh) Hirsh.write_file(wavy[0], true);

    log2 << " done!" << endl;

    if (esp) {
      log2 << "Calculating ESP..." << flush;
      Calc_ESP(ESP, wavy[0], ncpus, radius, log2);
      log2 << "Writing cube to Disk..." << flush;
      ESP.write_file(wavy[0], true);
      log2 << "  done!" << endl;
    }
    return 0;
  }
  cout << NoSpherA2_message() << "Did not understand the task to perform!\n" << help_message() << endl;
  return 0;
}