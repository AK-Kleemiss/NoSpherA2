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
#include "basis_set.h"
#include "structure_factors.h"

using namespace std;
bool debug_main=false;
bool expert=false;
int main(int argc, char **argv)
{
	ofstream log_file("wfn_2_fchk.log", ios::out);
	log_file << "             ____   _____ _______ " << std::endl;
	log_file << "            / __ \\ / ____|__   __|" << std::endl;
	log_file << "  ___ _   _| |  | | |       | |   " << std::endl;
	log_file << " / __| | | | |  | | |       | |   " << std::endl;
	log_file << "| (__| |_| | |__| | |____   | |   " << std::endl;
	log_file << " \\___|\\__,_|\\___\\_\\\\_____|  |_|   " << std::endl;
	log_file << "                                  " << std::endl;
	log_file << "This software is part of the cuQCT software suite developed by Florian Kleemiss.\nPlease give credit and cite corresponding pieces!\n";
	log_file.flush();
	if(argc>1){
		string keyword=argv[1];
		if(keyword.find("--help")!=string::npos){
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
  if(argc<4){
    log_file << "Not enough arguments given, at least provide -wfn <FILENAME>.wfn -b <basis_set>" << endl;
	log_file.flush();
	log_file.close();
    return -1;
  }
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
		log_file << "argc:"<< argc << endl;
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
	if (temp.find("-mode-2") != string::npos)
	  becke = true;
	if (temp.find("-ED") != string::npos)
	  electron_diffraction = true;
    if (temp.find("-v") != string::npos) {
      log_file << "Turning on verbose mode!" << endl;
      debug_main = true;
    }
  }
  if (threads != -1) {
	  omp_set_num_threads(threads);
	  omp_set_dynamic(0);
  }

  log_file.flush();
  if(debug_main)
  	log_file << "status:" << wfn << "&" << fchk << "&" << basis_set << "&" << basis_set_path << "&" << cif << "&" << hkl << endl;
  if(wfn == ""){
  	string filename;
  	vector <string> temp;
  	temp.resize(2);
  	temp[0]="Wavefunction files (wfn,ffn) | *.wfn *.ffn";
  	temp[1]="All filetypes | *";
  	if(!open_file_dialog(filename,debug_main,temp))
  		log_file << "Error encountered!" << endl;
  	else{
  		int origin=filetype_identifier(filename, debug_main);
  		if (debug_main) log_file << "Origin: " << origin << endl;
  		wavy.push_back(WFN(origin));
    }
  }
  else{
  	if(wfn.find(".wfn")!=string::npos) wavy.push_back(WFN(2));
  	else wavy.push_back(WFN(4));
  }
  if(wfn.find(".wfn")==string::npos&&wfn.find(".ffn")==string::npos){
      log_file << "Argument to -wfn does not look like wfn or ffn file!" << endl;
      return -1;
  }
  log_file.flush();
  //readwfn(wavy[0], wfn, debug_main, log_file);
  bool temp_debug = false;
  readwfn(wavy[0], wfn, debug_main, log_file);
  wavy[0].set_method(method);
  if (electron_diffraction)
	  log_file << "Making Electron diffraction scattering factors, be carefull what you are doing!" << endl;

  if(basis_set!=""||fchk!=""){
    if(basis_set == ""){
      string filename;
      vector <string> temp;
      temp.resize(1);
      temp[0]="All filetypes | *";
      if(!open_file_dialog(filename,debug_main,temp))
        log_file << "Error encountered!" << endl;
      else
        basis_set_path = filename;
    }
		else
			join_path(basis_set_path,basis_set);
		if(!exists(basis_set_path)){
			log_file << "Basis set file does not exist!" << endl;
			return -1;
		}
		wavy[0].set_basis_set_name(basis_set_path);

		vector<string> endings;
		endings.push_back(".fchk");
		endings.push_back(".Fchk");
		endings.push_back(".FChk");
		string outputname;
		if(fchk != "")
			outputname = fchk;
		else{
		    outputname = wavy[0].get_path();
		    if(debug_main) log_file << "Loaded path..." << endl;
			size_t where;
			if(wavy[0].get_origin()==2) where = outputname.find("wfn");
			else where = outputname.find("ffn");
			if(where>=outputname.length()&&where!=string::npos){
				if(!expert)
					save_file_dialog(outputname,debug_main,endings);
				else {
					log_file << "Enter filename: ";
					cin >> outputname;
					while(exists(outputname)){
						log_file << outputname << " exists, do you want to overwrite it? ";
						if(!yesno()){
							log_file << "Then try again: ";
							cin >> outputname;
						}
					}
   				}
			}
			else
				outputname.erase(where,3);
		}

		wavy[0].assign_charge(wavy[0].calculate_charge());
		if (mult == 0) {
			wavy[0].guess_multiplicity(log_file,expert);
		}
		else
			wavy[0].assign_multi(mult);
		free_fchk(log_file,outputname,"", wavy[0], debug_main,true);
	}
	if(cif!=""||hkl!=""){
		//debug_main = true;
		if(debug_main)
			log_file << "Entering Structure Factor Calculation!" << endl;
		if(!calculate_structure_factors(hkl,cif,asym_cif,symm,wavy[0],debug_main,accuracy, log_file, becke, threads, electron_diffraction, pbc))
			log_file << "Error during SF Calculation!" << endl;
	}
	return 0;
}

