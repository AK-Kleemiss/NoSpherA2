#ifndef __c_structure_factors__
#define __c_structure_factors__

#include <string>
#include <fstream>


class WFN;
class tsc_block;

bool merge_tscs(
	const std::string &mode,
	const std::vector<std::string> &files,
	const bool debug = false
	);

bool merge_tscs_without_checks(
	const std::string& mode,
	const std::vector<std::string>& files,
	const bool debug = false
);

bool thakkar_sfac(
	std::string& hkl_filename,
	std::string& cif,
	bool debug,
	std::ofstream& file,
	std::vector <int>& input_groups,
	std::vector < std::vector <double> >& twin_law,
	WFN &wave,
	int cpus = -1,
	bool electron_diffraction = false,
	bool save_k_pts = false,
	bool read_k_pts = false);

tsc_block MTC_thakkar_sfac(
	std::string& hkl_filename,
	std::string& cif,
	bool debug,
	std::ofstream& file,
	std::vector <int>& input_groups,
	std::vector < std::vector <double> >& twin_law,
	WFN& wave,
	int cpus = -1,
	bool electron_diffraction = false,
	bool save_k_pts = false,
	bool read_k_pts = false);

bool calculate_structure_factors_HF(
	std::string &hkl_filename, 
	std::string &cif_file, 
	WFN &wave, 
	bool debug, 
	int accuracy, 
	std::ofstream &file, 
	std::vector <int> &input_groups,
	std::vector <std::vector <double> > &twin_law,
	int cpus = -1,
	bool electron_diffraction = false,
	int pbc = 0,
	bool Olex2_1_3_switch = false,
	bool save_k_pts = false,
	bool read_k_pts = false);

tsc_block calculate_structure_factors_MTC(
	std::string& hkl_filename,
	std::string& cif_file,
	WFN& wave,
	bool debug,
	int accuracy,
	std::ofstream& file,
	std::vector <int>& input_groups,
	std::vector <std::vector <double> >& twin_law,
	std::vector <std::string> known_atoms,
	int cpus = -1,
	bool electron_diffraction = false,
	bool save_k_pts = false,
	bool read_k_pts = false);

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
#include "tsc_block.h"

#endif
