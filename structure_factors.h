#pragma once

#include <string>
#include <fstream>

class WFN;
class tsc_block;

bool merge_tscs(
  const std::string& mode,
  const std::vector<std::string>& files,
  const bool debug = false
);

bool merge_tscs_without_checks(
  const std::string& mode,
  const std::vector<std::string>& files,
  const bool debug = false
);

bool thakkar_sfac(
  const std::string& hkl_filename,
  const std::string& cif,
  bool debug,
  std::ofstream& file,
  const std::vector <int>& input_groups,
  const std::vector < std::vector <double> >& twin_law,
  WFN& wave,
  int cpus = -1,
  bool electron_diffraction = false,
  bool save_k_pts = false,
  bool read_k_pts = false);

tsc_block MTC_thakkar_sfac(
  const std::string& hkl_filename,
  const std::string& cif,
  bool debug,
  std::ofstream& file,
  const std::vector <int>& input_groups,
  const std::vector < std::vector <double> >& twin_law,
  std::vector < std::string >& known_atoms,
  WFN& wave,
  int cpus = -1,
  bool electron_diffraction = false,
  bool save_k_pts = false,
  bool read_k_pts = false);

bool calculate_structure_factors_HF(
  const std::string& hkl_filename,
  const std::string& cif_file,
  WFN& wave,
  bool debug,
  int accuracy,
  std::ofstream& file,
  const std::vector <int>& input_groups,
  const std::vector <std::vector <double> >& twin_law,

  const int cpus = -1,
  const bool electron_diffraction = false,
  const int pbc = 0,
  const bool Olex2_1_3_switch = false,
  const bool save_k_pts = false,
  const bool read_k_pts = false,
  const int& ECP_mode = 0,
  const bool no_date = false);

tsc_block calculate_structure_factors_MTC(
  const std::string& hkl_filename,
  const std::string& cif_file,
  WFN& wave,
  bool debug,
  int accuracy,
  std::ofstream& file,
  const std::vector <int>& input_groups,
  const std::vector <std::vector <double> >& twin_law,
  std::vector <std::string>& known_atoms,

  const int cpus = -1,
  const bool electron_diffraction = false,
  const int pbc = 0,
  const bool save_k_pts = false,
  const bool read_k_pts = false,
  const int& ECP_mode = 0,
  const bool no_date = false);
/*
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
  */
#include "wfn_class.h"
#include "tsc_block.h"
