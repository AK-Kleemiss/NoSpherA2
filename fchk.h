
#ifndef FCHK_H__
#define FCHK_H__

#include <fstream>

class WFN;

bool chk2fchk(const std::string &outputname, const std::string &gaussian_path);

bool gaussian(const std::string &programPath, const bool &debug);

std::string prepare_gaussian(const std::string &basis_set_path, const std::string &fchkname, 
                     WFN &wave, const int &ncpus, const float &mem, bool debug);

bool new_fchk(const std::string &basis_set_path, const std::string &fchkname, 
              const std::string &gaussian_path, WFN &wave, const int &ncpus, const float &mem, bool debug);

bool modify_fchk(const std::string &fchk_name, const std::string &basis_set_path, WFN &wave, bool &debug, const bool &read);
bool free_fchk(const std::string &fchk_name, const std::string &basis_set_path, WFN &wave, bool &debug, bool force_overwrite = false);
bool free_fchk(std::ofstream& file, const std::string& fchk_name, const std::string& basis_set_path, WFN& wave, bool& debug, bool force_overwrite = false);

#include "wfn_class.h"

#endif
