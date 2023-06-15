#pragma once
#include "spherical_density.h"
#include "convenience.h"
#include "npy.h"


void thakkar_d_test(options& opt) {
  using namespace std;
  Thakkar Os(76), Ca(20), C(6), O(8), H(1), P(15);
  Thakkar_Cation C_cat(6), O_cat(8), P_cat(15), Ca_cat(20);
  Thakkar_Anion C_an(6), O_an(8), H_an(1), P_an(15);
  double k_value = 0.0;
  if (!opt.electron_diffraction) {
    ofstream result("thakkar.dat", ios::out);
    for (double i = 0.001; i <= 4.0; i += 0.001) {
      k_value = bohr2ang(FOUR_PI * i);
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
  else {
    ofstream result("thakkar_ED.dat", ios::out);
    for (double i = 0.001; i <= 4.0; i += 0.001) {
      k_value = bohr2ang(FOUR_PI * i);
      result << showpos << setw(6) << setprecision(3) << fixed << i;
      complex<double> temp = H.get_form_factor(k_value);
      result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(1, temp, i).real();
      temp = C.get_form_factor(k_value);
      result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(6, temp, i).real();
      temp = O.get_form_factor(k_value);
      result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(8, temp, i).real();
      temp = P.get_form_factor(k_value);
      result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(15, temp, i).real();
      temp = Ca.get_form_factor(k_value);
      result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(20, temp, i).real();
      temp = Os.get_form_factor(k_value);
      result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(76, temp, i).real();
      temp = 0.0;
      result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(1, temp, i).real();
      temp = C_cat.get_form_factor(k_value);
      result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(6, temp, i).real();
      temp = O_cat.get_form_factor(k_value);
      result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(8, temp, i).real();
      temp = P_cat.get_form_factor(k_value);
      result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(15, temp, i).real();
      temp = Ca_cat.get_form_factor(k_value);
      result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(20, temp, i).real();
      temp = H_an.get_form_factor(k_value);
      result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(1, temp, i).real();
      temp = C_an.get_form_factor(k_value);
      result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(6, temp, i).real();
      temp = O_an.get_form_factor(k_value);
      result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(8, temp, i).real();
      temp = P_an.get_form_factor(k_value);
      result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(15, temp, i).real();
      temp = C.get_form_factor(k_value);
      result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(6, 1, temp, i).real();
      temp = O.get_form_factor(k_value);
      result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(8, 1, temp, i).real();
      temp = C.get_form_factor(k_value);
      result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(6, -1, temp, i).real();
      temp = O.get_form_factor(k_value);
      result << showpos << setw(16) << setprecision(8) << scientific << convert_to_ED_single(7, -1, temp, i).real();
      result << endl;
    }
    result.flush();
    result.close();
  }
}

void test_density_cubes(options& opt, std::ofstream& log_file) {
  using namespace std;
  vector<WFN> wavy(10);
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
}

void sfac_scan(options& opt, std::ofstream& log_file) {
  using namespace std;
  std::vector<WFN> wavy;
  wavy.push_back(WFN(1));
  wavy[0].read_known_wavefunction_format(opt.wfn, std::cout, opt.debug);
  Thakkar O(wavy[0].atoms[0].charge);
  Thakkar_Cation O_cat(wavy[0].atoms[0].charge);
  Thakkar_Anion O_an(wavy[0].atoms[0].charge);
  err_checkf(wavy[0].get_ncen() != 0, "No Atoms in the wavefunction, this will not work!! ABORTING!!", std::cout);
  err_checkf(exists(opt.cif), "CIF does not exists!", std::cout);
  //err_checkf(exists(asym_cif), "Asym/Wfn CIF does not exists!", file);
  if (opt.threads != -1) {
    omp_set_num_threads(opt.threads);
    omp_set_dynamic(0);
  }

#ifdef _WIN64
  time_t start = time(NULL);
  time_t end_becke, end_prototypes, end_spherical, end_prune, end_aspherical;
#else
  struct timeval t1, t2;

  gettimeofday(&t1, 0);
#endif

  cell unit_cell(opt.cif, std::cout, opt.debug);
  ifstream cif_input(opt.cif.c_str(), std::ios::in);
  vector <int> atom_type_list;
  vector <int> asym_atom_to_type_list;
  vector <int> asym_atom_list;
  vector <bool> needs_grid(wavy[0].get_ncen(), false);
  vector<string> known_atoms;

  read_atoms_from_CIF(cif_input,
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
  vector<vector<double>> d1, d2, d3, dens;

  make_hirshfeld_grids(opt.pbc,
    opt.accuracy,
    unit_cell,
    wavy[0],
    atom_type_list,
    asym_atom_list,
    needs_grid,
    d1, d2, d3, dens,
    std::cout,
#ifdef _WIN64
    start,
    end_becke,
    end_prototypes,
    end_spherical,
    end_prune,
    end_aspherical,
#else
    t1,
    t2,
#endif
    opt.debug,
    opt.no_date);


  std::cout << "finished partitioning" << endl;
  const int size = 4000;
  const int phi_size = 50;
  const int theta_size = 50;
  const double phi_step = 360.0 / phi_size * PI_180;
  const double theta_step = 180.0 / phi_size * PI_180;

  //This bit is basically the substitute for make_k_pts, where we sample the whole sphere 
  // by iterating over both spherical angles by a fixed step defined above
  vector<vector<double>> k_pt;
  k_pt.resize(4);
#pragma omp parallel for
  for (int i = 0; i < 4; i++)
    k_pt[i].resize(size * phi_size * theta_size, 0.0);

  //int null = 0;
#pragma omp parallel for
  for (int ref = 1; ref <= size; ref++) {
    for (int p = 0; p < phi_size; p++) {
      for (int t = 0; t < theta_size; t++) {
        int ind = t + (p + (ref - 1) * phi_size) * theta_size;
        double k_length = bohr2ang(FOUR_PI * ref / size * opt.sfac_scan);
        k_pt[0][ind] = k_length * sin(t * theta_step) * cos(p * phi_step);
        k_pt[1][ind] = k_length * sin(t * theta_step) * sin(p * phi_step);
        k_pt[2][ind] = k_length * cos(t * theta_step);
        k_pt[3][ind] = k_length;
      }
    }
  }
  // below is a strip of Calc_SF without the file IO or progress bar
  vector<vector<complex<double>>> sf;

  const int imax = (int)dens.size();
  const int smax = (int)k_pt[0].size();
  int pmax = (int)dens[0].size();
  const int step = max((int)floor(smax / 20), 1);
  std::cout << "Done with making k_pt " << smax << " " << imax << " " << pmax << endl;
  sf.resize(imax);
#pragma omp parallel for
  for (int i = 0; i < imax; i++)
    sf[i].resize(k_pt[0].size());
  double* dens_local, * d1_local, * d2_local, * d3_local;
  complex<double>* sf_local;
  const double* k1_local = k_pt[0].data();
  const double* k2_local = k_pt[1].data();
  const double* k3_local = k_pt[2].data();
  double work, rho;
  progress_bar* progress = new progress_bar{ std::cout, 60u, "Calculating scattering factors" };
  for (int i = 0; i < imax; i++) {
    pmax = (int)dens[i].size();
    dens_local = dens[i].data();
    d1_local = d1[i].data();
    d2_local = d2[i].data();
    d3_local = d3[i].data();
    sf_local = sf[i].data();
#pragma omp parallel for private(work,rho)
    for (int s = 0; s < smax; s++) {
      for (int p = pmax - 1; p >= 0; p--) {
        rho = dens_local[p];
        work = k1_local[s] * d1_local[p] + k2_local[s] * d2_local[p] + k3_local[s] * d3_local[p];
        sf_local[s] += complex<double>(rho * cos(work), rho * sin(work));
      }
      if (s != 0 && s % step == 0)
        progress->write(s / double(smax));
    }
  }
  delete(progress);
  if (true) { //Change if oyu do not want X-ray
    ofstream result("sfacs.dat", ios::out);
    log_file << "Writing X-ray sfacs...";
    log_file.flush();
    //Now we just need to write the result to a file, together with the spherical results and separated for valence and core
    for (int i = 0; i < k_pt[0].size(); i++) {
      result << showpos << setw(8) << setprecision(5) << fixed << ang2bohr(k_pt[3][i] / FOUR_PI);
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
  if (true) { //change if you do not want ED sfacs
    log_file << "Writing ED sfacs...";
    log_file.flush();
    ofstream result = ofstream("sfacs_ED.dat", ios::out);
    const double fact = 0.023934;
    double h2;
    for (int s = 0; s < k_pt[0].size(); s++) {
      h2 = pow(ang2bohr(k_pt[3][s] / FOUR_PI), 2);
      sf[0][s] = std::complex<double>(fact * (wavy[0].get_atom_charge(0) - sf[0][s].real()) / h2, -fact * sf[0][s].imag() / h2);

      result << showpos << setw(8) << setprecision(5) << fixed << ang2bohr(k_pt[3][s] / FOUR_PI);
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
  log_file.close();
}

void spherical_harmonic_test() {
  const double phi = 0.3, theta = 0.4;
  for (int lam = 0; lam <= 5; lam++) {
    for (int m = -lam; m <= lam; m++) {
      std::vector<double> d {std::sin(theta)*std::cos(phi),
                             std::sin(theta)*std::sin(phi),
                             std::cos(theta), 1.0, 1.0};
      std::cout << spherical_harmonic(lam, m, d.data()) << " ";
    }
    std::cout << "\n";
  }
};

const double gaussian_radial(primitive& p, double& r) {
  return pow(r, p.type) * std::exp(-p.exp * r * r) * p.norm_const;
}

double calc_density_ML(double& x, 
                       double& y, 
                       double& z, 
                       vec& coefficients,
                       std::vector<atom>& atoms,
                       const int& exp_coefs) {
  double dens = 0, radial = 0;
  int coef_counter = 0;
  int e = 0, size = 0;
  for (int a = 0; a < atoms.size(); a++) {
    size = (int)atoms[a].basis_set.size();
    basis_set_entry* bf;
    double d[4]{
    x - atoms[a].x,
    y - atoms[a].y,
    z - atoms[a].z, 0.0 };
    //store r in last element
    d[3] = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    //normalize distances for spherical harmonic
    for(e = 0; e<3; e++)
      d[e] /= d[3];
    for (e = 0; e < size; e++) {
      bf = &atoms[a].basis_set[e];
      primitive p(a, bf->type, bf->exponent, bf->coefficient);
      radial = gaussian_radial(p, d[3]);
      for (int m = -p.type; m <= p.type; m++) {
        dens += coefficients[coef_counter + m + p.type] * radial * spherical_harmonic(p.type, m, d);
      }
      coef_counter += (2 * p.type + 1);
    }
  }
  err_checkf(coef_counter == exp_coefs, "WRONG NUMBER OF COEFFICIENTS! " + std::to_string(coef_counter) + " vs. " + std::to_string(exp_coefs), std::cout);
  return dens;
}

std::vector<std::vector<primitive>> TZVP_JKfit(
  {
    { //center  type  exp         coef
      {0,       0, 9.5302493327,  1.0},
      {0,       0, 1.9174506246,  1.0},
      {0,       0, 0.68424049142, 1.0},
      {0,       0, 0.28413255710, 1.0},
      {0,       1, 2.9133232035,  1.0},
      {0,       1, 1.2621205398,  1.0},
      {0,       1, 0.50199775874, 1.0},
      {0,       2, 2.3135329149,  1.0},
      {0,       2, 0.71290724024, 1.0},
      {0,       3, 1.6565726132 , 1.0},
    },//H
    {},//He
    {},//Li
    {},//Be
    {},//B
    {
      {0,        0, 1113.9867719, 1.0},
      {0,        0, 369.16234180, 1.0},
      {0,        0, 121.79275232, 1.0},
      {0,        0, 48.127114540, 1.0},
      {0,        0, 20.365074004, 1.0},
      {0,        0, 8.0883596856, 1.0},
      {0,        0, 2.5068656570, 1.0},
      {0,        0, 1.2438537380, 1.0},
      {0,        0,0.48449899601, 1.0},
      {0,        0,0.19185160296, 1.0},
      {0,        1, 102.99176249, 1.0},
      {0,        1, 28.132594009, 1.0},
      {0,        1, 9.8364318173, 1.0},
      {0,        1, 3.3490544980, 1.0},
      {0,        1, 1.4947618613, 1.0},
      {0,        1,0.57690108899, 1.0},
      {0,        1,0.20320063291, 1.0},
      {0,        2, 10.594068356, 1.0},
      {0,        2, 3.5997195366, 1.0},
      {0,        2, 1.3355691094, 1.0},
      {0,        2,0.51949764954, 1.0},
      {0,        2,0.19954125200, 1.0},
      {0,        3,1.194866338369, 1.0},
      {0,        3, .415866338369, 1.0},
      {0,        4, .858866338369, 1.0}
    },//C
    {
      {0,        0, 1102.8622453, 1.0},
      {0,        0, 370.98041153, 1.0},
      {0,        0, 136.73555938, 1.0},
      {0,        0, 50.755871924, 1.0},
      {0,        0, 20.535656095, 1.0},
      {0,        0, 7.8318737184, 1.0},
      {0,        0, 3.4784063855, 1.0},
      {0,        0, 1.4552856603, 1.0},
      {0,        0,0.63068989071, 1.0},
      {0,        0,0.27276596483, 1.0},
      {0,        1, 93.540954073, 1.0},
      {0,        1, 29.524019527, 1.0},
      {0,        1, 10.917502987, 1.0},
      {0,        1, 4.3449288991, 1.0},
      {0,        1, 1.8216912640, 1.0},
      {0,        1,0.75792424494, 1.0},
      {0,        1,0.28241469033, 1.0},
      {0,        2, 16.419378926, 1.0},
      {0,        2, 5.0104049385, 1.0},
      {0,        2, 1.9793971884, 1.0},
      {0,        2,0.78495771518, 1.0},
      {0,        2,0.28954065963, 1.0},
      {0,        3,1.79354239843, 1.0},
      {0,        3, .60854239843, 1.0},
      {0,        4,1.23254239843, 1.0}
    },//N
    {
      {0,        0, 1517.8667506, 1.0},
      {0,        0, 489.67952008, 1.0},
      {0,        0, 176.72118665, 1.0},
      {0,        0, 63.792233137, 1.0},
      {0,        0, 25.366499130, 1.0},
      {0,        0, 9.9135491200, 1.0},
      {0,        0, 4.4645306584, 1.0},
      {0,        0, 1.8017743661, 1.0},
      {0,        0,0.80789710965, 1.0},
      {0,        0,0.33864326862, 1.0},
      {0,        1, 120.16030921, 1.0},
      {0,        1, 34.409622474, 1.0},
      {0,        1, 12.581148610, 1.0},
      {0,        1, 5.0663824249, 1.0},
      {0,        1, 2.0346927092, 1.0},
      {0,        1,0.86092967212, 1.0},
      {0,        1,0.36681356726, 1.0},
      {0,        2, 19.043062805, 1.0},
      {0,        2, 5.8060381104, 1.0},
      {0,        2, 2.1891841580, 1.0},
      {0,        2,0.87794613558, 1.0},
      {0,        2,0.35623646700, 1.0},
      {0,        3,2.493914788135, 1.0},
      {0,        3, .824914788135, 1.0},
      {0,        4,1.607914788135, 1.0}
    },//O
    {},//F
    {} //Ne
  }
);

std::vector<std::vector<primitive>> QZVP_JKfit(
  { 
    { //center  type  exp         coef
      {0,       0, 9.5302493327,  1.0},
      {0,       0, 1.9174506246,  1.0},
      {0,       0, 0.68424049142, 1.0},
      {0,       0, 0.28413255710, 1.0},
      {0,       1, 2.9133232035,  1.0},
      {0,       1, 1.2621205398,  1.0},
      {0,       1, 0.50199775874, 1.0},
      {0,       2, 2.8832083931,  1.0},
      {0,       2, 1.2801701725 , 1.0},
      {0,       2, 0.52511317770, 1.0},
      {0,       3, 2.7489448439 , 1.0},
      {0,       3, 1.1900885456 , 1.0},
      {0,       4, 1.4752662714 , 1.0}
    },//H
    {},//He
    {},//Li
    {},//Be
    {},//B
    {
      {0,        0, 1113.9867719, 1.0},
      {0,        0, 369.16234180, 1.0},
      {0,        0, 121.79275232, 1.0},
      {0,        0, 48.127114540, 1.0},
      {0,        0, 20.365074004, 1.0},
      {0,        0, 8.0883596856, 1.0},
      {0,        0, 2.5068656570, 1.0},
      {0,        0, 1.2438537380, 1.0},
      {0,        0,0.48449899601, 1.0},
      {0,        0,0.19185160296, 1.0},
      {0,        1, 102.99176249, 1.0},
      {0,        1, 28.132594009, 1.0},
      {0,        1, 9.8364318173, 1.0},
      {0,        1, 3.3490544980, 1.0},
      {0,        1, 1.4947618613, 1.0},
      {0,        1,0.57690108899, 1.0},
      {0,        1,0.20320063291, 1.0},
      {0,        2, 10.594068356, 1.0},
      {0,        2, 3.5997195366, 1.0},
      {0,        2, 1.3355691094, 1.0},
      {0,        2,0.51949764954, 1.0},
      {0,        2,0.19954125200, 1.0},
      {0,        3,1.95390000000, 1.0},
      {0,        3, .75490000000, 1.0},
      {0,        3, .33390000000, 1.0},
      {0,        4,1.52490000000, 1.0},
      {0,        4, .59090000000, 1.0},
      {0,        5,1.11690000000, 1.0}
    },//C
    {
      {0,        0, 1102.8622453, 1.0},
      {0,        0, 370.98041153, 1.0},
      {0,        0, 136.73555938, 1.0},
      {0,        0, 50.755871924, 1.0},
      {0,        0, 20.535656095, 1.0},
      {0,        0, 7.8318737184, 1.0},
      {0,        0, 3.4784063855, 1.0},
      {0,        0, 1.4552856603, 1.0},
      {0,        0,0.63068989071, 1.0},
      {0,        0,0.27276596483, 1.0},
      {0,        1, 93.540954073, 1.0},
      {0,        1, 29.524019527, 1.0},
      {0,        1, 10.917502987, 1.0},
      {0,        1, 4.3449288991, 1.0},
      {0,        1, 1.8216912640, 1.0},
      {0,        1,0.75792424494, 1.0},
      {0,        1,0.28241469033, 1.0},
      {0,        2, 16.419378926, 1.0},
      {0,        2, 5.0104049385, 1.0},
      {0,        2, 1.9793971884, 1.0},
      {0,        2,0.78495771518, 1.0},
      {0,        2,0.28954065963, 1.0},
      {0,        3,2.98600000000, 1.0},
      {0,        3,1.11700000000, 1.0},
      {0,        3, .48400000000, 1.0},
      {0,        4,2.17600000000, 1.0},
      {0,        4, .83400000000, 1.0},
      {0,        5,1.57600000000, 1.0}
    },//N
    {
      {0,        0, 1517.8667506, 1.0},
      {0,        0, 489.67952008, 1.0},
      {0,        0, 176.72118665, 1.0},
      {0,        0, 63.792233137, 1.0},
      {0,        0, 25.366499130, 1.0},
      {0,        0, 9.9135491200, 1.0},
      {0,        0, 4.4645306584, 1.0},
      {0,        0, 1.8017743661, 1.0},
      {0,        0,0.80789710965, 1.0},
      {0,        0,0.33864326862, 1.0},
      {0,        1, 120.16030921, 1.0},
      {0,        1, 34.409622474, 1.0},
      {0,        1, 12.581148610, 1.0},
      {0,        1, 5.0663824249, 1.0},
      {0,        1, 2.0346927092, 1.0},
      {0,        1,0.86092967212, 1.0},
      {0,        1,0.36681356726, 1.0},
      {0,        2, 19.043062805, 1.0},
      {0,        2, 5.8060381104, 1.0},
      {0,        2, 2.1891841580, 1.0},
      {0,        2,0.87794613558, 1.0},
      {0,        2,0.35623646700, 1.0},
      {0,        3,3.96585000000, 1.0},
      {0,        3,1.49085000000, 1.0},
      {0,        3, .63485000000, 1.0},
      {0,        4,2.85685000000, 1.0},
      {0,        4,1.04985000000, 1.0},
      {0,        5,2.03685000000, 1.0}
    },//O
    {},//F
    {} //Ne
  }
);

void calc_cube(vec data, WFN& dummy, int& exp_coef) {
  double MinMax[6]{ 0,0,0,0,0,0 };
  int steps[3]{ 0,0,0 };
  readxyzMinMax_fromWFN(dummy, MinMax, steps, 2.5, 0.1);
  cube CubeRho(steps[0], steps[1], steps[2], dummy.get_ncen(), true);
  CubeRho.give_parent_wfn(dummy);

  for (int i = 0; i < 3; i++) {
    CubeRho.set_origin(i, MinMax[i]);
    CubeRho.set_vector(i, i, (MinMax[i + 3] - MinMax[i]) / steps[i]);
  }
  CubeRho.set_comment1("Calculated density using NoSpherA2 from ML Data");
  CubeRho.set_comment2("from " + dummy.get_path());
  CubeRho.path = get_basename_without_ending(dummy.get_path()) + "_rho.cube";

  time_t start;
  std::time(&start);

  progress_bar* progress = new progress_bar{ std::cout, 50u, "Calculating Values" };
  const int step = (int)std::max(floor(CubeRho.get_size(0) / 20.0), 1.0);

#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < CubeRho.get_size(0); i++) {
    for (int j = 0; j < CubeRho.get_size(1); j++)
      for (int k = 0; k < CubeRho.get_size(2); k++) {

        vec PosGrid{
          i * CubeRho.get_vector(0, 0) + j * CubeRho.get_vector(0, 1) + k * CubeRho.get_vector(0, 2) + CubeRho.get_origin(0),
          i * CubeRho.get_vector(1, 0) + j * CubeRho.get_vector(1, 1) + k * CubeRho.get_vector(1, 2) + CubeRho.get_origin(1),
          i * CubeRho.get_vector(2, 0) + j * CubeRho.get_vector(2, 1) + k * CubeRho.get_vector(2, 2) + CubeRho.get_origin(2)
        };

        CubeRho.set_value(i, j, k, calc_density_ML(PosGrid[0], PosGrid[1], PosGrid[2], data, dummy.atoms, exp_coef));
      }
    if (i != 0 && i % step == 0)
      progress->write(i / double(CubeRho.get_size(0)));
  }
  delete(progress);

  time_t end;
  std::time(&end);
  if (difftime(end, start) < 60) std::cout << "Time to calculate Values: " << std::fixed << std::setprecision(0) << difftime(end, start) << " s" << std::endl;
  else if (difftime(end, start) < 3600) std::cout << "Time to calculate Values: " << std::fixed << std::setprecision(0) << std::floor(difftime(end, start) / 60) << " m " << int(std::floor(difftime(end, start))) % 60 << " s" << std::endl;
  else std::cout << "Time to calculate Values: " << std::fixed << std::setprecision(0) << std::floor(difftime(end, start) / 3600) << " h " << (int(std::floor(difftime(end, start))) % 3600) / 60 << " m" << std::endl;
  std::cout << "Number of electrons: " << std::fixed << std::setprecision(4) << CubeRho.sum() << std::endl;
  CubeRho.write_file(true);
};

double numerical_3d_integral(std::function<double(double*)> f, double stepsize = 0.001, double r_max = 10.0) {
  double tot_int = 0;
  double da = stepsize * PI;
  double dr = stepsize;
  int upper = TWO_PI / da;
  long long int total_calcs = upper * PI / da * r_max / dr;
  std::cout << "Integrating over " << total_calcs << " points" << std::endl;

  progress_bar* progress = new progress_bar{ std::cout, 50u, "Integrating" };
  const long long int step = (long long int)std::max(floor(upper / 20.0), 1.0);

#pragma omp parallel for reduction(+:tot_int) schedule(dynamic)
  for (int phic = 0; phic <= upper; phic++) {
    double phi = phic * da;
    double cp = std::cos(phi);
    double sp = std::sin(phi);
    for (double theta = 0; theta <= PI; theta += da) {
      double st = std::sin(theta);
      double ct = std::cos(theta);
      double x0 = st * cp, y0 = st * sp;
      double ang_mom = dr * da * da * st;
      for (double r = r_max; r >= 0; r -= dr) {
        double d[4] {r* x0,
          r* y0,
          r* ct,
          r};
        tot_int += f(d) * r * r * ang_mom;
      }
    }
    if (phic != 0 && phic % step == 0)
      progress->write((double)phic / (double)upper);
  }

  std::cout << "\ntotal integral : " << std::fixed << std::setprecision(4) << tot_int << std::endl;
}

double s_value(double* d) {
  primitive p_0(0, 0, 3.5, 1.0);
  int m = 0;
  return pow(gaussian_radial(p_0, d[3]) * spherical_harmonic(p_0.type, m, d),2);
}

double one(double* d) {
  return 1;
}

double gaussian(double* d) {
  double g = exp(-3.5 * d[3] * d[3]);
  return g;
}

int load_basis_into_WFN(WFN& wavy, auto b) {
  int nr_coefs = 0;
  for (int i = 0; i < wavy.atoms.size(); i++) {
    int current_charge = wavy.atoms[i].charge - 1;
    primitive* basis = b[current_charge].data();
    int size = (int)b[current_charge].size();
    for (int e = 0; e < size; e++) {
      wavy.atoms[i].push_back_basis_set(basis[e].exp, 1.0, basis[e].type, e);
      nr_coefs += 2 * basis[e].type + 1;
    }
  }
  return nr_coefs;
}

void ML_test() {
  using namespace std;
  //cout << "c_1_4p " << c_1_4p << endl;
  //cout << "1/c_1_4p " << 1/ c_1_4p << endl;
  //cout << "4pi " << FOUR_PI << endl;
  //cout << "4/3 pi " << 4.0 / 3.0 * PI << endl;
  //expected value integrating to r=1 
  //for f=one is 4/3 pi
  //for f=gaussian: 0.789257
  //for f=s_value = 1?
  //numerical_3d_integral(s_value,0.002,10.0);

  vector<unsigned long> shape {};
  bool fortran_order;
  vec data{};

  string path {"coefficients_conf0.npy"};
  npy::LoadArrayFromNumpy(path, shape, fortran_order, data);

  WFN dummy(7);
  dummy.read_xyz("water_monomer.xyz", std::cout);

  int nr_coefs = load_basis_into_WFN(dummy, QZVP_JKfit);

  cout << data.size() << " vs. " << nr_coefs << " ceofficients" << endl;
  calc_cube(data, dummy, nr_coefs);

  path = "prediction_conf0.npy";
  npy::LoadArrayFromNumpy(path, shape, fortran_order, data);

  WFN dummy2(7);
  dummy2.read_xyz("water.xyz", std::cout);

  nr_coefs = load_basis_into_WFN(dummy2, QZVP_JKfit);;
  
  cout << data.size() << " vs. " << nr_coefs << " ceofficients" << endl;
  calc_cube(data, dummy2, nr_coefs);

  path = "alanine2.npy";
  npy::LoadArrayFromNumpy(path, shape, fortran_order, data);

  WFN dummy3(7);
  dummy3.read_xyz("alanine.xyz", std::cout);

  nr_coefs = load_basis_into_WFN(dummy3, TZVP_JKfit);;

  cout << data.size() << " vs. " << nr_coefs << " ceofficients" << endl;
  calc_cube(data, dummy3, nr_coefs);

  cout << "Done :)" << endl;
}

