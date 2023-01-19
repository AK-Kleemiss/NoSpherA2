#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <complex>

class WFN;
class tsc_block;
class cell;
struct options;

bool thakkar_sfac(
  const options& opt,
  std::ofstream& file,
  WFN& wave);

tsc_block MTC_thakkar_sfac(
  options& opt,
  std::ofstream& file,
  std::vector < std::string >& known_atoms,
  std::vector<WFN>& wave,
  const int& nr);

bool calculate_structure_factors_HF(
  const options& opt,
  WFN& wave,
  std::ofstream& file);

tsc_block calculate_structure_factors_MTC(
  options& opt,
  std::vector<WFN>& wave,
  std::ofstream& file,
  std::vector <std::string>& known_atoms,
  const int& nr,
  std::vector<std::vector<double>>* kpts=NULL);

std::complex<double> convert_to_ED_single(const int& charge,
  std::complex<double>& sf,
  const double& k_vector);

std::complex<double> convert_to_ED_single(const int& neutralcharge,
  const int& charge,
  std::complex<double>& sf,
  const double& k_vector);

std::complex<double> convert_to_ED_single(const int& neutralcharge,
  const double& charge,
  std::complex<double>& sf,
  const double& k_vector);

void read_atoms_from_CIF(std::ifstream& cif_input,
  const std::vector <int>& input_groups,
  const cell& unit_cell,
  WFN& wave,
  const std::vector <std::string>& known_atoms,
  std::vector <int>& atom_type_list,
  std::vector <int>& asym_atom_to_type_list,
  std::vector <int>& asym_atom_list,
  std::vector <bool>& needs_grid,
  std::ostream& file,
  const bool debug = false);

int make_hirshfeld_grids(const int& pbc,
  const int& accuracy,
  const std::vector <int>& input_groups,
  cell& unit_cell,
  const WFN& wave,
  const std::vector <int>& atom_type_list,
  const std::vector <int>& asym_atom_to_type_list,
  const std::vector <int>& asym_atom_list,
  std::vector <bool>& needs_grid,
  std::vector<std::vector<double>>& d1,
  std::vector<std::vector<double>>& d2,
  std::vector<std::vector<double>>& d3,
  std::vector<std::vector<double>>& dens,
  std::ostream& file,
#ifdef _WIN64
  time_t& start,
  time_t& end_becke,
  time_t& end_prototypes,
  time_t& end_spherical,
  time_t& end_prune,
  time_t& end_aspherical,
#else
  timeval& t1,
  timeval& t2,
#endif
  bool debug = false,
  bool no_date = false);

void calc_SF(const int& points,
  std::vector<std::vector<double>>& k_pt,
  std::vector<std::vector<double>>& d1,
  std::vector<std::vector<double>>& d2,
  std::vector<std::vector<double>>& d3,
  std::vector<std::vector<double>>& dens,
  std::vector<std::vector<std::complex<double>>>& sf,
  std::ostream& file,
#ifdef _WIN64
  time_t& start,
  time_t& end1,
#else
  timeval& t1,
  timeval& t2,
#endif
  bool debug = false);

#include "wfn_class.h"
#include "tsc_block.h"
