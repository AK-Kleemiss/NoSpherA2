#ifndef __c_structure_factors__
#define __c_structure_factors__

#include <string>
#include <fstream>


class WFN;

bool merge_tscs(
	const std::string &mode,
	const std::vector<std::string> &files,
	const bool debug = false
	);

bool thakkar_sfac(
	const int atom_number,
	const double radial_acc = 1E-25,
	const int angular_level = 3);

bool calculate_structure_factors_HF(
	std::string &hkl_filename, 
	std::string &cif_file, 
	std::string &asym_cif, 
	std::string &symm, 
	WFN &wave, 
	bool debug, 
	int accuracy, 
	std::ofstream &file, 
	std::vector <int> &input_groups,
	std::vector <std::vector <double> > &twin_law,
	int cpus = -1,
	bool electron_diffraction = false,
	int pbc = 0,
	bool Olex2_1_3_switch = false);

bool calculate_structure_factors_RI(
	std::string& hkl_filename,
	std::string& cif_file,
	std::string& asym_cif,
	std::string& symm,
	WFN& wave,
	bool debug,
	int accuracy,
	std::ofstream& file,
	int cpus = -1,
	bool electron_diffraction = false,
	int pbc = 0);

#include "wfn_class.h"

#endif
