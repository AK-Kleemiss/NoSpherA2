/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * main.cc
 * Copyright (C) 2016 Florian Kleemiss <florian.kleemiss@uni-bremen.de>
 * 
 * wfn-cpp is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * wfn-cpp is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fstream>
#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#include <sys/wait.h>
#include <termios.h>
#endif
#include <limits>
#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <iomanip>
#include <vector>
#include <omp.h>
#include <limits>
#include <cstddef>
#include <thread>

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
int main(int argc, char **argv){
	vector <string> known_keywords;
	//known_keywords.push_back("-wfn");
	//known_keywords.push_back("-fchk");
	//known_keywords.push_back("-b");
	//known_keywords.push_back("-d");
	//known_keywords.push_back("-v");
	//known_keywords.push_back("-hkl");
	//known_keywords.push_back("-cif");
	//known_keywords.push_back("-acc");
	//known_keywords.push_back("-mult");
	//known_keywords.push_back("-method");
	//known_keywords.push_back("-symm");
	//known_keywords.push_back("-asym_cif");
	//known_keywords.push_back("-asym-cif");
	//known_keywords.push_back("-wfn-cif");
	//known_keywords.push_back("-wfn_cif");
	if(debug_main) 
		cout << "argc:"<< argc << endl;
	vector<WFN> wavy;
	string wavename;
	string gaussian_path;
	string turbomole_path;
	string basis_set_path;
	int ncpus=0;
	float mem=0.0;
	int mult = 0;
	/*int config=program_confi(gaussian_path, turbomole_path, basis_set_path, ncpus, mem, debug_main, expert);
	if(config==-1){
		log_file << "No .cuQCT.conf found, do you want to continue without the employment of external programs?" << endl;
		if(!yesno()){
			log_file << "Then please start again and write a config file or use the template choice." << endl;
			log_file.flush();
			log_file.close();
			return -1;
		}
	}
	else
		if(config==0) 
			log_file << "The programs file was newly written, please check if everything is correct!" << endl;
	*/
	string wfn("");
	string fchk("");
	string basis_set("");
	string hkl("");
	string cif("");
	string symm("");
	string asym_cif("");
	string method("rhf");
	string xyz_file("");
	string temp;
	string fract_name("");
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
	int hirsh_number = 0;
	bool scnd = true, thrd = true, frth = true;
	double MinMax[6];
	double NbSteps[3];
	vector <int> MOs;
	vector <int> groups;
	vector < vector <double> > twin_law;
	bool all_mos = false;
	for (int i=0; i<argc; i++){
		temp = argv[i];
		if (temp.find("-") > 0) continue;
		if(temp.find("-wfn") < 1)
			wfn = argv[i+1];
		else if(temp.find("-fchk") < 1)
			fchk = argv[i+1];
		else if(temp.find("-b") < 1)
			basis_set = argv[i+1];
		else if (temp.find("-d") < 1)
			basis_set_path = argv[i + 1];
		else if (temp.find("-hkl") < 1)
			hkl = argv[i+1];
		else if (temp.find("-cif") < 1)
			cif = argv[i+1];
		else if (temp.find("-acc") < 1)
			accuracy = stoi(argv[i+1]);
		else if (temp.find("-mult") < 1)
			mult = stoi(argv[i + 1]);
		else if (temp.find("-method") < 1)
			method = argv[i + 1];
		else if (temp.find("-symm") < 1)
			symm = argv[i + 1];
		else if (temp.find("-asym_cif") < 1)
			asym_cif = argv[i + 1];
		else if (temp.find("-asym-cif") < 1)
			asym_cif = argv[i + 1];
		else if (temp.find("-wfn-cif") < 1)
			asym_cif = argv[i + 1];
		else if (temp.find("-wfn_cif") < 1)
			asym_cif = argv[i + 1];
		else if (temp.find("-cpus") < 1)
			threads = stoi(argv[i + 1]);
		else if (temp.find("-pbc") < 1)
			pbc = stoi(argv[i + 1]);
		else if (temp.find("-ED") < 1)
			electron_diffraction = true;
		else if (temp.find("-Olex2_1_3") < 1)
			Olex2_1_3_switch = true;
		else if (temp.find("-v2") < 1) {
			cout << "Turning on verbose mode 2!" << endl;
			debug_all = true;
			debug_main = true;
		}
		else if (temp.find("-v") < 1) {
			cout << "Turning on verbose mode!" << endl;
			debug_main = true;
		}
		else if (temp.find("-eli") < 1) {
			calc = true;
			eli = true;
		}
		else if (temp.find("-elf") < 1) {
			calc = true;
			elf = true;
		}
		else if (temp.find("-lap") < 1) {
			calc = true;
			lap = true;
		}
		else if (temp.find("-esp") < 1) {
			calc = true;
			esp = true;
		}
		else if (temp.find("-rdg") < 1) {
			calc = true;
			rdg = true;
		}
		else if (temp.find("-hirsh") < 1) {
			calc = true;
			hirsh = true;
			hirsh_number = stoi(argv[i + 1]);
		}
		else if (temp.find("-resolution") < 1) {
			resolution = stod(argv[i + 1]);
		}
		else if (temp.find("-radius") < 1) {
			radius = stod(argv[i + 1]);
		}
		else if (temp.find("-MO") < 1) {
			if (string(argv[i + 1]) != "all") 
				MOs.push_back(stoi(argv[i + 1]));
			else
				all_mos = true;
			calc = true;
		}
		else if (temp.find("-fractal") < 1) {
			fract = true;
			fract_name = argv[i + 1];
		}
		else if (temp.find("-HDEF") < 1) {
			hdef = true;
			calc = true;
		}
		else if (temp.find("-def") < 1) {
			def = true;
			calc = true;
		}
		else if (temp.find("-skpts") < 1) {
			save_k_pts = true;
		}
		else if (temp.find("-rkpts") < 1) {
			read_k_pts = true;
		}
		else if (temp.find("-group") < 1) {
			while (string(argv[i + 1]).find("-") != string::npos) {
				groups.push_back(stoi(argv[i + 1]));
			}
		}
		else if (temp.find("-twin") < 1) {
			twin_law.resize(twin_law.size() + 1);
			twin_law[twin_law.size() - 1].resize(9);
			for (int twl = 0; twl < 9; twl++) {
				twin_law[twin_law.size() - 1][twl] = stod(argv[i + 1 + twl]);
			}
			if (debug_main) {
				cout << "twin_law: ";
				for (int twl = 0; twl < 9; twl++)
					cout << setw(7) << setprecision(2) << twin_law[twin_law.size() - 1][twl];
				cout << endl;
			}
			i += 9;
		}
		else if (temp.find("-IAM") != string::npos) {
			iam_switch = true;
		}
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
		}
	}
	if (threads != -1) {
		omp_set_num_threads(threads);
		omp_set_dynamic(0);
	}
	string help_message = "----------------------------------------------------------------------------\n";
	help_message.append("		These commands and arguments are known by wfn2fchk:\n");
	help_message.append("----------------------------------------------------------------------------\n\n");
	help_message.append("	-wfn  <FILENAME>.wfn/ffn/wfx	Read the following file\n");
	help_message.append("	-fchk <FILENAME>.fchk			Write a wavefunction to the given filename\n");
	help_message.append("	-b    <FILENAME>				Read this basis set\n");
	help_message.append("	-d    <PATH>					Path to basis_sets directory with basis_sets in tonto style\n");
	help_message.append("	--help 							print this help\n");
	help_message.append("	-v 		    					Turn on Verbose (debug) Mode (Slow and a LOT of output!)\n");
	help_message.append("   -v2                             Even more stuff\n");
	help_message.append("	-mult	<NUMBER>				Input multiplicity of wavefunction\n");
	help_message.append("   -method <METHOD NAME>           Can be RKS or RHF to distinguish between DFT and HF\n");
	help_message.append("   -cif      <FILENAME>.cif        CIF to use for calculation of scatteriung factors\n");
	help_message.append("   -asym-cif <FILENAME>.cif        CIF to use to identify the asymetric unit atoms from\n");
	help_message.append("   -IAM                            Make scattering factors based on Thakkar functions for elements in CIF\n");
	help_message.append("   -xyz                            Read atom positions from this xyz file for IAM\n");
	help_message.append("   -hkl      <FILENAME>.hkl        hkl file (ideally merged) to use for calculation of form factors\n");
	help_message.append("   -twin     3x3 floating-matrix in the form -1 0 0 0 -1 0 0 0 -1 which contains the twin matrix to use\n");
	help_message.append("              If there is more than a single twin law to be used use the twin command multiple times.\n");
	help_message.append("   -merge    <List of .tsc files>  Names/Paths to .tsc files to be merged. They need to have identical hkl values.\n");
	//"    -e                          Turn on expert mode (Disable a lot of assumptions and enter paths manually)" << endl
	help_message.append("*****************************************************************************************************\n");
	help_message.append("	Explanation: A .ffn file is an extended wfn file containing also unoccupied orbitals\n");
	help_message.append("*****************************************************************************************************\n");
	string NoSpherA2_message = "    _   __     _____       __              ___   ___\n";
	NoSpherA2_message.append("   / | / /___ / ___/____  / /_  ___  _____/   | |__ \\\n");
	NoSpherA2_message.append("  /  |/ / __ \\\\__ \\/ __ \\/ __ \\/ _ \\/ ___/ /| | __/ /\n");
	NoSpherA2_message.append(" / /|  / /_/ /__/ / /_/ / / / /  __/ /  / ___ |/ __/\n");
	NoSpherA2_message.append("/_/ |_/\\____/____/ .___/_/ /_/\\___/_/  /_/  |_/____/\n");
	NoSpherA2_message.append("                /_/\n");
	NoSpherA2_message.append("This software is part of the cuQCT software suite developed by Florian Kleemiss.\n");
	NoSpherA2_message.append("Please give credit and cite corresponding pieces!\n");
	NoSpherA2_message.append("NoSpherA2 was published at : Kleemiss et al. Chem.Sci., 2021, 12, 1675 - 1692\n");
	if (argc > 1) {
		string keyword = argv[1];
		if (keyword.find("--help") != string::npos)
			cout << help_message << endl;
	}
	if (fract) {
		cube residual(fract_name,true,debug_all);
		residual.fractal_dimension(0.01);
	}
	if (!calc && wfn != "" && hkl == "" && fchk == "") {
		ofstream log_file("NoSpherA2.log", ios::out);
		if (debug_main)
			for (int i = 0; i < argc; i++)
				log_file << argv[i] << endl;

		log_file << NoSpherA2_message;
		log_file.flush();
		if (wfn.find(".wfn") != string::npos) wavy.push_back(WFN(2));
		else if (wfn.find(".ffn") != string::npos) wavy.push_back(WFN(4));
		else if (wfn.find(".wfx") != string::npos) wavy.push_back(WFN(6));
		else if (wfn.find(".fch") != string::npos) wavy.push_back(WFN(4));
		else {
			log_file << "Argument to -wfn does not look like wfn, wfx or ffn file!" << endl;
			return -1;
		}
		log_file.flush();
		if (wfn.find(".wfx") != string::npos)
			wavy[0].read_wfx(wfn, debug_all, log_file);
		else if (wfn.find(".fchk") != string::npos) {
			if (debug_all) log_file << "Reading FCHK" << endl;
			wavy[0].read_fchk(wfn, log_file, debug_all);
			if (debug_all) wavy[0].write_wfn(wfn+"_test.wfn", false, false);
		}
		else
			wavy[0].read_wfn(wfn, debug_all, log_file);
	}

	if (iam_switch) {
		ofstream log_file("NoSpherA2.log", ios::out);
		if (debug_main)
			for (int i = 0; i < argc; i++)
				log_file << argv[i] << endl;

		log_file << NoSpherA2_message;
		log_file.flush();
		if (argc > 1) {
			string keyword = argv[1];
			if (keyword.find("--help") != string::npos) {
				log_file << help_message << endl;
				log_file.close();
				return 0;
			}
		}

		log_file.flush();
		if (debug_main)
			log_file << "status: " << cif << "&" << hkl << "&" << asym_cif << "&" << xyz_file << endl;
		if (xyz_file == "") {
			log_file << "No .xyz file specified! Please provide an xyz file for calculation of IAM .tscs" << endl;
			return -1;
		}
		else {
			if (exists(xyz_file)) {
				wavy.push_back(WFN(7));
				wavy[0].read_xyz(xyz_file, log_file, debug_main);
			}
			else {
				log_file << ".xyz file does not exist!" << endl;
				return -1;
			}
		}

		if (electron_diffraction)
			log_file << "Making Electron diffraction scattering factors, be carefull what you are doing!" << endl;
		if (cif != "" || hkl != "") {
			if (debug_main)
				log_file << "Entering Structure Factor Calculation!" << endl;
			if (!thakkar_sfac(hkl, cif, asym_cif, debug_main, log_file, groups, twin_law, wavy[0], threads, electron_diffraction))
				log_file << "Error during SF Calculation!" << endl;
		}
	}
	else if (hkl != "" || basis_set != "" || fchk != "") {
		ofstream log_file("NoSpherA2.log", ios::out);
		if (debug_main)
			for (int i = 0; i < argc; i++)
				log_file << argv[i] << endl;

		log_file << NoSpherA2_message;
		log_file.flush();
		if (argc > 1) {
			string keyword = argv[1];
			if (keyword.find("--help") != string::npos) {
				log_file << help_message << endl;
				log_file.close();
				return 0;
			}
		}
		if (argc < 4) {
			cout << "Not enough arguments given, at least provide -wfn <FILENAME>.wfn/.wfx -b <basis_set>" << endl;
			return -1;
		}

		log_file.flush();
		if (debug_main)
			log_file << "status:" << wfn << "&" << fchk << "&" << basis_set << "&" << basis_set_path << "&" << cif << "&" << hkl << "&" << asym_cif << endl;
		if (wfn == "") {
			log_file << "Error, no wfn file specified!";
			return -1;
		}
		else {
			if (wfn.find(".wfn") != string::npos) wavy.push_back(WFN(2));
			else if (wfn.find(".ffn") != string::npos) wavy.push_back(WFN(4));
			else if (wfn.find(".wfx") != string::npos) wavy.push_back(WFN(6));
			else if (wfn.find(".fch") != string::npos) wavy.push_back(WFN(4));
			else {
				log_file << "Argument to -wfn does not look like wfn, wfx or ffn file!" << endl;
				return -1;
			}
		}
		log_file << "Reading: " << wfn;
		log_file.flush();
		if (wfn.find(".wfx") != string::npos)
			wavy[0].read_wfx(wfn, debug_all, log_file);
		else if (wfn.find(".fchk") != string::npos) {
			wavy[0].read_fchk(wfn, log_file, debug_all);
			if(debug_all) wavy[0].write_wfn("test.wfn", false, true);
		}
		else
			wavy[0].read_wfn(wfn, debug_all, log_file);
		wavy[0].set_method(method);
		log_file << "                   ...done!" << endl;
		if (electron_diffraction)
			log_file << "Making Electron diffraction scattering factors, be carefull what you are doing!" << endl;

		if (basis_set != "" || fchk != "") {
			join_path(basis_set_path, basis_set);
			if (!exists(basis_set_path)) {
				log_file << "Basis set file does not exist!" << endl;
				return -1;
			}
			wavy[0].set_basis_set_name(basis_set_path);

			vector<string> endings;
			endings.push_back(".fchk");
			endings.push_back(".Fchk");
			endings.push_back(".FChk");
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
				if (where >= outputname.length() && where != string::npos) {
					log_file << "Cannot make output file name!";
					return -1;
				}
				else
					outputname.erase(where, 3);
			}

			wavy[0].assign_charge(wavy[0].calculate_charge());
			if (mult == 0) {
				wavy[0].guess_multiplicity(log_file, expert);
			}
			else
				wavy[0].assign_multi(mult);
			free_fchk(log_file, outputname, "", wavy[0], debug_main, true);
		}
		if (cif != "" || hkl != "") {
			if (debug_main)
				log_file << "Entering Structure Factor Calculation!" << endl;
			if (!calculate_structure_factors_HF(
					hkl, 
					cif, 
					asym_cif, 
					symm, 
					wavy[0], 
					debug_main, 
					accuracy, 
					log_file, 
					groups, 
					twin_law, 
					threads, 
					electron_diffraction, 
					pbc, 
					Olex2_1_3_switch,
					save_k_pts,
					read_k_pts)
				)
				log_file << "Error during SF Calculation!" << endl;
		}
	}
	if (calc) {
		ofstream log2("NoSpherA2_cube.log", ios::out);
		if (debug_main)
			for (int i = 0; i < argc; i++)
				log2 << argv[i] << endl;
		log2 << NoSpherA2_message;
		log2.flush();
#ifdef _WIN32
		if(debug_main){
			log2 << "float min: " << dec << FLT_MIN << endl;
			log2 << "float max: " << dec << FLT_MAX << endl;
			log2 << "float Min_exp: " << dec << FLT_MIN_EXP << endl;
			log2 << "float Max_exp: " << dec << FLT_MAX_EXP << endl;

			log2 << "double min: " << dec << DBL_MIN << endl;
			log2 << "double max: " << dec << DBL_MAX << endl;
			log2 << "double Min_exp: " << dec << DBL_MIN_EXP << endl;
			log2 << "double Max_exp: " << dec << DBL_MAX_EXP << endl;

			log2 << "long double min: " << dec << LDBL_MIN << endl;
			log2 << "long double max: " << dec << LDBL_MAX << endl;
			log2 << "long double Min_exp: " << dec << LDBL_MIN_EXP << endl;
			log2 << "long double Max_exp: " << dec << LDBL_MAX_EXP << endl;
		}
#endif

		if (wfn == "") {
			log2 << "Error, no wfn file specified!";
			return -1;
		}
		else {
			if (wfn.find(".wfn") != string::npos) wavy.push_back(WFN(2));
			else if (wfn.find(".ffn") != string::npos) wavy.push_back(WFN(4));
			else if (wfn.find(".wfx") != string::npos) wavy.push_back(WFN(6));
			else {
				log2 << "Argument to -wfn does not look like wfn, wfx or ffn file!" << endl;
				return -1;
			}
		}
		log2.flush();
		if (wfn.find(".wfx") == string::npos)
			wavy[0].read_wfn(wfn, debug_all, log2);
		else
			wavy[0].read_wfx(wfn, debug_all, log2);
		if (debug_main)
			log2 << "Starting calcualtion of properties" << endl;
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
			cube temp(NbSteps[0], NbSteps[1], NbSteps[2], wavy[0].get_ncen(), hdef||hirsh);
			for (int i = 0; i < 3; i++) {
				temp.set_origin(i, MinMax[i]);
				for (int j = 0; j < 3; j++)
					temp.set_vector(i, j, cell_matrix[i][j]);
			}
			if (hdef || hirsh) {
				log2 << "Calcualting spheircal Rho...";
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

			if(hdef)
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
			Calc_Prop(Rho, RDG, Elf, Eli, Lap, ESP, wavy[0], ncpus, radius, log2);

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
	}
	return 0;
}

