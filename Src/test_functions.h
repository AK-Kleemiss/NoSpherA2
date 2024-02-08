#pragma once
#include "spherical_density.h"
#include "structure_factors.h"
#include "convenience.h"
#include "npy.h"
#include "properties.h"


void thakkar_d_test(options& opt) {
	using namespace std;
	Thakkar Os(76), Ca(20), C(6), O(8), H(1), P(15);
	Thakkar_Cation C_cat(6), O_cat(8), P_cat(15), Ca_cat(20);
	Thakkar_Anion C_an(6), O_an(8), H_an(1), P_an(15);
	double k_value = 0.0;
	if (!opt.electron_diffraction) {
		ofstream result("thakkar.dat", ios::out);
		for (double i = 0.001; i <= 4.0; i += 0.001) {
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
	else {
		ofstream result("thakkar_ED.dat", ios::out);
		for (double i = 0.001; i <= 4.0; i += 0.001) {
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

void test_density_cubes(options& opt, std::ostream& log_file) {
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

void sfac_scan(options& opt, std::ostream& log_file) {
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
	//err_checkf(exists(asym_cif), "Asym/Wfn CIF does not exists!", file);

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
	vector<vec> d1, d2, d3, dens;

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
	const double phi_step = 360.0 / phi_size * constants::PI_180;
	const double theta_step = 180.0 / phi_size * constants::PI_180;

	//This bit is basically the substitute for make_k_pts, where we sample the whole sphere 
	// by iterating over both spherical angles by a fixed step defined above
	vector<vec> k_pt;
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
				double k_length = constants::bohr2ang(constants::FOUR_PI * ref / size * opt.d_sfac_scan);
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
	if (true) { //change if you do not want ED sfacs
		log_file << "Writing ED sfacs...";
		log_file.flush();
		ofstream result = ofstream("sfacs_ED.dat", ios::out);
		const double fact = 0.023934;
		double h2;
		for (int s = 0; s < k_pt[0].size(); s++) {
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

void spherical_harmonic_test() {
	const double phi = 0.3, theta = 0.4;
	for (int lam = 0; lam <= 5; lam++) {
		for (int m = -lam; m <= lam; m++) {
			vec d{ std::sin(theta) * std::cos(phi),
				   std::sin(theta) * std::sin(phi),
				   std::cos(theta), 1.0, 1.0 };
			std::cout << spherical_harmonic(lam, m, d.data()) << " ";
		}
		std::cout << "\n";
	}
};

void calc_cube(vec data, WFN& dummy, int& exp_coef, int atom = -1) {
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
	if (atom != -1)
		std::cout << "Calculation for atom " << atom << std::endl;

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < CubeRho.get_size(0); i++) {
		for (int j = 0; j < CubeRho.get_size(1); j++)
			for (int k = 0; k < CubeRho.get_size(2); k++) {

				vec PosGrid{
					i * CubeRho.get_vector(0, 0) + j * CubeRho.get_vector(0, 1) + k * CubeRho.get_vector(0, 2) + CubeRho.get_origin(0),
					i * CubeRho.get_vector(1, 0) + j * CubeRho.get_vector(1, 1) + k * CubeRho.get_vector(1, 2) + CubeRho.get_origin(1),
					i * CubeRho.get_vector(2, 0) + j * CubeRho.get_vector(2, 1) + k * CubeRho.get_vector(2, 2) + CubeRho.get_origin(2)
				};

				if (atom == -1)
					CubeRho.set_value(i, j, k, calc_density_ML(PosGrid[0], PosGrid[1], PosGrid[2], data, dummy.atoms, exp_coef));
				else
					CubeRho.set_value(i, j, k, calc_density_ML(PosGrid[0], PosGrid[1], PosGrid[2], data, dummy.atoms, exp_coef, atom));
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
	if (atom == -1) {
		CubeRho.write_file(true);
	}
	else {
		std::string fn(get_basename_without_ending(dummy.get_path()) + "_rho_" + std::to_string(atom) + ".cube");
		CubeRho.write_file(fn, false);
	}
};

double numerical_3d_integral(std::function<double(double*)> f, double stepsize = 0.001, double r_max = 10.0) {
	double tot_int = 0;
	double da = stepsize * constants::PI;
	double dr = stepsize;
	int upper = int(constants::TWO_PI / da);
	long long int total_calcs = upper * (long long int)(constants::PI / da * r_max / dr);
	std::cout << "Integrating over " << total_calcs << " points" << std::endl;

	progress_bar* progress = new progress_bar{ std::cout, 50u, "Integrating" };
	const long long int step = (long long int)std::max(floor(upper / 20.0), 1.0);

#pragma omp parallel for reduction(+:tot_int) schedule(dynamic)
	for (int phic = 0; phic <= upper; phic++) {
		double phi = phic * da;
		double cp = std::cos(phi);
		double sp = std::sin(phi);
		for (double theta = 0; theta <= constants::PI; theta += da) {
			double st = std::sin(theta);
			double ct = std::cos(theta);
			double x0 = st * cp, y0 = st * sp;
			double ang_mom = dr * da * da * st;
			for (double r = r_max; r >= 0; r -= dr) {
				double d[4]{ r * x0,
					r * y0,
					r * ct,
					r };
				tot_int += f(d) * r * r * ang_mom;
			}
		}
		if (phic != 0 && phic % step == 0)
			progress->write((double)phic / (double)upper);
	}

	std::cout << "\ntotal integral : " << std::fixed << std::setprecision(4) << tot_int << std::endl;
	return tot_int;
}

double s_value(double* d) {
	primitive p_0(0, 0, 3.5, 1.0);
	int m = 0;
	return pow(gaussian_radial(p_0, d[3]) * spherical_harmonic(p_0.type, m, d), 2);
}

double one(double* d) {
	return 1;
}

double gaussian(double* d) {
	double g = exp(-3.5 * d[3] * d[3]);
	return g;
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

	vector<unsigned long> shape{};
	bool fortran_order;
	vec data{};

	string path{ "coefficients_conf0.npy" };
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

void cube_from_coef_npy(std::string& coef_fn, std::string& xyzfile) {
	std::vector<unsigned long> shape{};
	bool fortran_order;
	vec data{};

	npy::LoadArrayFromNumpy(coef_fn, shape, fortran_order, data);
	WFN dummy(7); dummy.read_xyz(xyzfile, std::cout);

	int nr_coefs = load_basis_into_WFN(dummy, TZVP_JKfit);
	std::cout << data.size() << " vs. " << nr_coefs << " ceofficients" << std::endl;
	calc_cube(data, dummy, nr_coefs);

	for (int i = 0; i < dummy.get_ncen(); i++)
		calc_cube(data, dummy, nr_coefs, i);
}

void test_xtb_molden(options& opt, std::ostream& log_file) {
	using namespace std;
	for (int i = 0; i < 1; i++) {
		WFN wavy(8);
		wavy.read_molden("Co2.molden", cout, true);
		opt.cif = "Co2.cif";
		opt.dmin = 0.5;
		cout << "STARTING CALC" << endl;
		calculate_structure_factors_HF(
			opt,
			wavy,
			log_file);
	}
}

void test_core_dens() {
	using namespace std;
	Thakkar T_Rb(37);
	string base = "def2-ECP";
	Gaussian_Atom G_Rb(37, base);
	double PI = 3.14159265358979323846;
	double TWO_PI = 2 * PI;
	double FOUR_PI = 4 * PI;
	ofstream out("core_dens.dat", ios::out);
	WFN ECP_way(9);
	ECP_way.read_gbw("Atomic_densities\\atom_atom37.gbw", out, true, true);


	WFN multi_ecp_wavy(9);
	multi_ecp_wavy.read_gbw("Rb4.gbw", out, true, true);	

	//exit(0);
	WFN wavy(9);
	wavy.read_gbw("Rb.gbw", cout, false);
	wavy.delete_unoccupied_MOs();
	WFN wavy2(9);
	wavy2.read_gbw("Rb_+9.gbw", cout, false);
	wavy2.delete_unoccupied_MOs();

	vector<vec> d;
	d.resize(16);
	for (int i = 0; i < 16; i++)
		d[i].resize(wavy.get_ncen(), 0.0);

	vec phi(wavy.get_nmo(), 0.0);

	for (int i = 1; i < 1300000; i++) {
		double r = i * 0.001;
		double sr = r * 0.01;

		double tsr = (sr);
		cout << fixed << r << " " << T_Rb.get_core_form_factor(r, 28) << " " << G_Rb.get_core_form_factor(r, 28);
		cout << " " << T_Rb.get_core_density(sr, 28) << " " << G_Rb.get_radial_density(tsr);
		cout << " " << T_Rb.get_radial_density(sr) << " " << wavy.compute_dens(sr, 0, 0, d, phi, false) << " " << wavy2.compute_dens(sr, 0, 0, d, phi, false);
		cout << " " << sr << " " << ECP_way.compute_dens(sr, 0, 0, d, phi, false);
		cout << "\n";
	}
	cout << flush;

}

