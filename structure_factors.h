#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <complex>

class WFN;
class tsc_block;
struct options;

bool thakkar_sfac(
  const options& opt,
  std::ofstream& file,
  WFN& wave);

tsc_block MTC_thakkar_sfac(
  const options& opt,
  std::ofstream& file,
  std::vector < std::string >& known_atoms,
  std::vector<std::vector<int>>& known_indices,
  std::vector<WFN>& wave,
  const int& nr);

bool calculate_structure_factors_HF(
  const options& opt,
  WFN& wave,
  std::ofstream& file);

tsc_block calculate_structure_factors_MTC(
  const options& opt,
  std::vector<WFN>& wave,
  std::ofstream& file,
  std::vector <std::string>& known_atoms,
  std::vector<std::vector<int>>& known_indices,
  const int& nr);

inline std::complex<double> convert_to_ED_single(const int& charge,
  std::complex<double>& sf,
  const double& k_vector)
{
  const double fact = 0.023934;
  const double h2 = pow(k_vector, 2);
  return std::complex<double>(fact * (charge - sf.real()) / h2, -fact * sf.imag() / h2);
}

#include "wfn_class.h"
#include "tsc_block.h"
