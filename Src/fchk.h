#pragma once

#include "convenience.h"

class WFN;

//std::string prepare_gaussian(const std::string& basis_set_path, const std::string& fchkname,
//  WFN& wave, const int& ncpus, const float& mem, bool debug);
bool modify_fchk(const std::string& fchk_name, const std::filesystem::path& basis_set_path, WFN& wave, bool& debug, const bool& read);
bool free_fchk(std::ostream& file, const std::filesystem::path& fchk_name, const std::filesystem::path& basis_set_path, WFN& wave, bool& debug, bool force_overwrite = false);
int read_fchk_integer(const std::string& in);
double read_fchk_double(const std::string& in);
bool read_fchk_integer_block(std::ifstream& in, const char* heading, ivec& result, bool rewind = true);
bool read_fchk_double_block(std::ifstream& in, const char* heading, vec& result, bool rewind = true);
int read_fchk_integer(std::ifstream& in, const char* search, bool rewind = true);
double read_fchk_double(std::ifstream& in, const char* search, bool rewind = true);

#include "wfn_class.h"