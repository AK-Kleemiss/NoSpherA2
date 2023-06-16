#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <complex>
#include "convenience.h"

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

bool calculate_structure_factors_RI(
  const options& opt,
  WFN& wave,
  std::ofstream& file,
  const int exp_coefs);

bool calculate_structure_factors_RI_No_H(
  const options& opt,
  WFN& wave,
  std::ofstream& file,
  const int exp_coefs);

tsc_block calculate_structure_factors_MTC(
  options& opt,
  std::vector<WFN>& wave,
  std::ofstream& file,
  std::vector <std::string>& known_atoms,
  const int& nr,
  std::vector<vec>* kpts=NULL);

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
  cell& unit_cell,
  const WFN& wave,
  const std::vector <int>& atom_type_list,
  const std::vector <int>& asym_atom_list,
  std::vector <bool>& needs_grid,
  std::vector<vec>& d1,
  std::vector<vec>& d2,
  std::vector<vec>& d3,
  std::vector<vec>& dens,
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
  std::vector<vec>& k_pt,
  std::vector<vec>& d1,
  std::vector<vec>& d2,
  std::vector<vec>& d3,
  std::vector<vec>& dens,
  std::vector<cvec>& sf,
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
