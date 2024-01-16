#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <complex>
#include "convenience.h"

class WFN;
template class tsc_block<int,cdouble>;
typedef tsc_block<int, cdouble> itsc_block;
class cell;
struct options;

bool thakkar_sfac(
  const options& opt,
  std::ostream& file,
  WFN& wave);

itsc_block MTC_thakkar_sfac(
  options& opt,
  std::ostream& file,
  std::vector < std::string >& known_atoms,
  std::vector<WFN>& wave,
  const int& nr);

bool calculate_structure_factors_HF(
  const options& opt,
  WFN& wave,
  std::ostream& file);

bool calculate_structure_factors_RI(
  const options& opt,
  WFN& wave,
  std::ostream& file,
  const int exp_coefs);

bool calculate_structure_factors_RI_No_H(
  const options& opt,
  WFN& wave,
  std::ostream& file,
  const int exp_coefs);

itsc_block calculate_structure_factors_MTC(
  options& opt,
  std::vector<WFN>& wave,
  std::ostream& file,
  std::vector <std::string>& known_atoms,
  const int& nr,
  std::vector<vec>* kpts=NULL);

void generate_hkl(const double& dmin,
  hkl_list& hkl,
  const std::vector<vec>& twin_law,
  cell& unit_cell,
  std::ostream& file,
  bool debug = false);

void generate_fractional_hkl(const double& dmin,
  hkl_list_d& hkl,
  const std::vector<vec>& twin_law,
  cell& unit_cell,
  std::ostream& file,
  double stepsize,
  bool debug);

template <typename NumType>
std::complex<double> convert_to_ED_single(const int& neutralcharge,
	std::complex<double>& sf,
	const double& k_vector,
	const NumType& charge = 0) {
	const double h2 = pow(k_vector, 2);
	std::complex<double> neutral(constants::ED_fact * (neutralcharge - sf.real()) / h2, -constants::ED_fact * sf.imag() / h2);
	if (charge == 0) return neutral;
	return neutral + constants::ED_fact * charge / h2;
}

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

void sfac_diffuse(options& opt, std::ofstream& log_file);

#include "wfn_class.h"
#include "tsc_block.h"
