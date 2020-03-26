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

#include "convenience.h"
#include "fchk.h"
#include "cube.h"
#include "basis_set.h"
#include "structure_factors.h"
#include "properties.h"

using namespace std;
bool debug_main = false;
bool debug_all = false;
bool expert = false;
int main(int argc, char **argv){
	vector <string> known_keywords;
	known_keywords.push_back("-wfn");
	known_keywords.push_back("-fchk");
	known_keywords.push_back("-b");
	known_keywords.push_back("-d");
	known_keywords.push_back("-v");
	known_keywords.push_back("-hkl");
	known_keywords.push_back("-cif");
	known_keywords.push_back("-acc");
	known_keywords.push_back("-mult");
	known_keywords.push_back("-method");
	known_keywords.push_back("-symm");
	known_keywords.push_back("-asym_cif");
	known_keywords.push_back("-asym-cif");
	known_keywords.push_back("-wfn-cif");
	known_keywords.push_back("-wfn_cif");
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
	string temp;
	int accuracy = 2;
	int threads = -1;
	bool becke = false;
	bool electron_diffraction = false;
	int pbc = 0;
	double resolution = 0.1;
	bool calc = false;
	bool eli = false;
	bool esp = false;
	bool elf = false;
	bool lap = false;
	bool rdg = false;
	double MinMax[6];
	double NbSteps[3];
	for (int i=0; i<argc; i++){
		temp = argv[i];
		if(temp.find(known_keywords[0]) != string::npos)
			wfn = argv[i+1];
		if(temp.find(known_keywords[1]) != string::npos)
			fchk = argv[i+1];
		if(temp.find(known_keywords[2]) != string::npos)
			basis_set = argv[i+1];
		if (temp.find(known_keywords[3]) != string::npos)
			basis_set_path = argv[i + 1];
		if (temp.find(known_keywords[5]) != string::npos)
			hkl = argv[i+1];
		if (temp.find(known_keywords[6]) != string::npos)
			cif = argv[i+1];
		if (temp.find(known_keywords[7]) != string::npos)
			accuracy = stoi(argv[i+1]);
		if (temp.find(known_keywords[8]) != string::npos)
			mult = stoi(argv[i + 1]);
		if (temp.find(known_keywords[9]) != string::npos)
			method = argv[i + 1];
		if (temp.find(known_keywords[10]) != string::npos)
			symm = argv[i + 1];
		if (temp.find(known_keywords[11]) != string::npos)
			asym_cif = argv[i + 1];
		if (temp.find(known_keywords[12]) != string::npos)
			asym_cif = argv[i + 1];
		if (temp.find(known_keywords[13]) != string::npos)
			asym_cif = argv[i + 1];
		if (temp.find(known_keywords[14]) != string::npos)
			asym_cif = argv[i + 1];
		if (temp.find("-cpus") != string::npos)
			threads = stoi(argv[i + 1]);
		if (temp.find("-pbc") != string::npos)
			pbc = stoi(argv[i + 1]);
		if (temp.find("-ED") != string::npos)
			electron_diffraction = true;
		if (temp.find("-v") != string::npos) {
			cout << "Turning on verbose mode!" << endl;
			debug_main = true;
		}
		if (temp.find("-v2") != string::npos) {
			cout << "Turning on verbose mode 2!" << endl;
			debug_all = true;
		}
		if (temp.find("-eli") != string::npos) {
			calc = true;
			eli = true;
		}
		if (temp.find("-elf") != string::npos) {
			calc = true;
			elf = true;
		}
		if (temp.find("-lap") != string::npos) {
			calc = true;
			lap = true;
		}
		if (temp.find("-esp") != string::npos) {
			calc = true;
			esp = true;
		}
		if (temp.find("-rdg") != string::npos) {
			calc = true;
			rdg = true;
		}
		if (temp.find("-resolution") != string::npos) {
			resolution = stod(argv[i + 1]);
		}
	}
	if (threads != -1) {
		omp_set_num_threads(threads);
		omp_set_dynamic(0);
	}
	if (hkl != "" || basis_set != "" || fchk != "") {
		ofstream log_file("NoSpherA2.log", ios::out);
		log_file << "    _   __     _____       __              ___   ___\n";
		log_file << "   / | / /___ / ___/____  / /_  ___  _____/   | |__ \\\n";
		log_file << "  /  |/ / __ \\\\__ \\/ __ \\/ __ \\/ _ \\/ ___/ /| | __/ /\n";
		log_file << " / /|  / /_/ /__/ / /_/ / / / /  __/ /  / ___ |/ __/\n";
		log_file << "/_/ |_/\\____/____/ .___/_/ /_/\\___/_/  /_/  |_/____/\n";
		log_file << "                /_/\n" << flush;
		log_file << "This software is part of the cuQCT software suite developed by Florian Kleemiss.\nPlease give credit and cite corresponding pieces!\n";
		log_file.flush();
		if (argc > 1) {
			string keyword = argv[1];
			if (keyword.find("--help") != string::npos) {
				log_file << "----------------------------------------------------------------------------" << endl
					<< "		These commands and arguments are known by wfn2fchk:" << endl
					<< "----------------------------------------------------------------------------" << endl << endl
					<< "	-wfn  <FILENAME>.wfn/ffn	Read the following file" << endl
					<< "	-fchk <FILENAME>.fchk		Write a wavefunction to the given filename" << endl
					<< "	-b <FILENAME>				Read this basis set" << endl
					<< "	-d <PATH>                   Path to basis_sets directory with basis_sets in tonto style" << endl
					<< "	--help 					    print this help" << endl
					<< "	-v 		    			    Turn on Verbose (debug) Mode (Slow and a LOT of output!)" << endl
					<< "	-mult                       Input multiplicity of wavefunction" << endl
					//				<< "    -e                          Turn on expert mode (Disable a lot of assumptions and enter paths manually)" << endl
					<< "*****************************************************************************************************" << endl
					<< "	Explanation: A .ffn file is an extended wfn file containing also unoccupied orbitals" << endl
					<< "*****************************************************************************************************" << endl;
				log_file.close();
				return 0;
			}
		}
		if (argc < 4) {
			cout << "Not enough arguments given, at least provide -wfn <FILENAME>.wfn -b <basis_set>" << endl;
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
			else wavy.push_back(WFN(4));
		}
		if (wfn.find(".wfn") == string::npos && wfn.find(".ffn") == string::npos) {
			log_file << "Argument to -wfn does not look like wfn or ffn file!" << endl;
			return -1;
		}
		log_file.flush();
		readwfn(wavy[0], wfn, debug_all, log_file);
		wavy[0].set_method(method);
		if (electron_diffraction)
			log_file << "Making Electron diffraction scattering factors, be carefull what you are doing!" << endl;

		if (basis_set != "" || fchk != "") {
			if (basis_set == "") {
				string filename;
				vector <string> temp;
				temp.resize(1);
				temp[0] = "All filetypes | *";
				/*    if(!open_file_dialog(filename,debug_main,temp))
					  log_file << "Error encountered!" << endl;
					else*/
				basis_set_path = filename;
			}
			else
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
				else where = outputname.find("ffn");
				if (where >= outputname.length() && where != string::npos) {
					/*if(!expert)
						save_file_dialog(outputname,debug_main,endings);
					else {*/
					log_file << "Enter filename: ";
					cin >> outputname;
					while (exists(outputname)) {
						log_file << outputname << " exists, do you want to overwrite it? ";
						if (!yesno()) {
							log_file << "Then try again: ";
							cin >> outputname;
						}
					}
					//}
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
			//debug_main = true;
			if (debug_main)
				log_file << "Entering Structure Factor Calculation!" << endl;
			if (!calculate_structure_factors(hkl, cif, asym_cif, symm, wavy[0], debug_main, accuracy, log_file, threads, electron_diffraction, pbc))
				log_file << "Error during SF Calculation!" << endl;
		}
	}
	if (calc) {
		ofstream log2("NoSpherA2_cube.log", ios::out);
		log2 << "    _   __     _____       __              ___   ___\n";
		log2 << "   / | / /___ / ___/____  / /_  ___  _____/   | |__ \\\n";
		log2 << "  /  |/ / __ \\\\__ \\/ __ \\/ __ \\/ _ \\/ ___/ /| | __/ /\n";
		log2 << " / /|  / /_/ /__/ / /_/ / / / /  __/ /  / ___ |/ __/\n";
		log2 << "/_/ |_/\\____/____/ .___/_/ /_/\\___/_/  /_/  |_/____/\n";
		log2 << "                /_/\n" << flush;
		log2 << "This software is part of the cuQCT software suite developed by Florian Kleemiss.\nPlease give credit and cite corresponding pieces!\n";

		if (wfn == "") {
			log2 << "Error, no wfn file specified!";
			return -1;
		}
		else {
			if (wfn.find(".wfn") != string::npos) wavy.push_back(WFN(2));
			else wavy.push_back(WFN(4));
		}
		if (wfn.find(".wfn") == string::npos && wfn.find(".ffn") == string::npos) {
			log2 << "Argument to -wfn does not look like wfn or ffn file!" << endl;
			return -1;
		}
		log2.flush();
		readwfn(wavy[0], wfn, debug_all, log2);
		if (debug_main)
			log2 << "Starting calcualtion of properties" << endl;

		vector < vector < double > > cell_matrix;
		cell_matrix.resize(3);
		for (int i = 0; i < 3; i++)
			cell_matrix[i].resize(3);

		readxyzMinMax_fromCIF(cif, MinMax, NbSteps, cell_matrix, resolution, debug_main);
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
		if (debug_main)
			log2 << "Rho is there" << endl;
		cube RDG(NbSteps[0], NbSteps[1], NbSteps[2], wavy[0].get_ncen(), rdg);
		if (debug_main)
			log2 << "RDG is there" << endl;
		cube Elf(NbSteps[0], NbSteps[1], NbSteps[2], wavy[0].get_ncen(), elf);
		if (debug_main)
			log2 << "Elf is there" << endl;
		cube Eli(NbSteps[0], NbSteps[1], NbSteps[2], wavy[0].get_ncen(), eli);
		if (debug_main)
			log2 << "Eli is there" << endl;
		cube Lap(NbSteps[0], NbSteps[1], NbSteps[2], wavy[0].get_ncen(), lap);
		if (debug_main)
			log2 << "Lap is there" << endl;
		cube ESP(NbSteps[0], NbSteps[1], NbSteps[2], wavy[0].get_ncen(), esp);
		if (debug_main)
			log2 << "ESP is there" << endl;
		for (int i = 0; i < 3; i++) {
			Rho.set_origin(i, MinMax[i]);
			RDG.set_origin(i, MinMax[i]);
			Elf.set_origin(i, MinMax[i]);
			Eli.set_origin(i, MinMax[i]);
			Lap.set_origin(i, MinMax[i]);
			ESP.set_origin(i, MinMax[i]);
			for (int j = 0; j < 3; j++) {
				Rho.set_vector(i, j, cell_matrix[i][j]);
				RDG.set_vector(i, j, cell_matrix[i][j]);
				Elf.set_vector(i, j, cell_matrix[i][j]);
				Eli.set_vector(i, j, cell_matrix[i][j]);
				Lap.set_vector(i, j, cell_matrix[i][j]);
				ESP.set_vector(i, j, cell_matrix[i][j]);
			}
		}
		Rho.set_comment1("Calculated density using NoSpherA2");
		RDG.set_comment1("Calculated reduced density gradient using NoSpherA2");
		Elf.set_comment1("Calculated electron localization function using NoSpherA2");
		Eli.set_comment1("Calculated same-spin electron localizability indicator using NoSpherA2");
		Lap.set_comment1("Calculated laplacian of electron density using NoSpherA2");
		ESP.set_comment1("Calculated electrostatic potential using NoSpherA2");
		Rho.set_comment2("from " + wavy[0].get_path());
		RDG.set_comment2("from " + wavy[0].get_path());
		Elf.set_comment2("from " + wavy[0].get_path());
		Eli.set_comment2("from " + wavy[0].get_path());
		Lap.set_comment2("from " + wavy[0].get_path());
		ESP.set_comment2("from " + wavy[0].get_path());
		Rho.path = get_basename_without_ending(wavy[0].get_path()) + "_rho.cube";
		RDG.path = get_basename_without_ending(wavy[0].get_path()) + "_rdg.cube";
		Elf.path = get_basename_without_ending(wavy[0].get_path()) + "_elf.cube";
		Eli.path = get_basename_without_ending(wavy[0].get_path()) + "_eli.cube";
		Lap.path = get_basename_without_ending(wavy[0].get_path()) + "_lap.cube";
		ESP.path = get_basename_without_ending(wavy[0].get_path()) + "_esp.cube";

		if (debug_main)
			log2 << "Everything is set up; starting calculation..." << endl;

		log2 << "Calculating for " << fixed << setprecision(0) << NbSteps[0] * NbSteps[1] * NbSteps[2] << " Gridpoints." << endl;

		if (lap || eli || elf || rdg)
			Calc_Prop(Rho, RDG, Elf, Eli, Lap, wavy[0], ncpus, log2);

		log2 << "Writing cubes to Disk..." << flush;
		if (rdg) {
			Rho.path = get_basename_without_ending(wavy[0].get_path()) + "_signed_rho.cube";
			Rho.write_file(wavy[0], true);
			Rho.path = get_basename_without_ending(wavy[0].get_path()) + "_rho.cube";
			Rho.write_file(wavy[0], true, true);
		}
		else if (lap || eli || elf) Rho.write_file(wavy[0], true);
		if (rdg) RDG.write_file(wavy[0], true);
		if (lap) Lap.write_file(wavy[0], true);
		if (elf) Elf.write_file(wavy[0], true);
		if (eli) Eli.write_file(wavy[0], true);

		log2 << " done!" << endl;

		if (esp) {
			Calc_ESP(ESP, wavy[0], ncpus, log2);
			log2 << "Writing cube to Disk..." << flush;
			ESP.write_file(wavy[0], true);
			log2 << "  done!" << endl;
		}
	}
	return 0;
}

