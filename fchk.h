#pragma once

#include <fstream>

class WFN;

//std::string prepare_gaussian(const std::string& basis_set_path, const std::string& fchkname,
//  WFN& wave, const int& ncpus, const float& mem, bool debug);
bool modify_fchk(const std::string& fchk_name, const std::string& basis_set_path, WFN& wave, bool& debug, const bool& read);
bool free_fchk(std::ofstream& file, const std::string& fchk_name, const std::string& basis_set_path, WFN& wave, bool& debug, bool force_overwrite = false);

#include "wfn_class.h"