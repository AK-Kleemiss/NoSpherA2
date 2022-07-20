#pragma once

#include <string>
#include <vector>
#include <fstream>

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

#include "wfn_class.h"
#include "tsc_block.h"
