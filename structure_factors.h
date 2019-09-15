#ifndef __c_structure_factors__
#define __c_structure_factors__

#include <string>
#include <fstream>


class WFN;

bool calculate_structure_factors(std::string &hkl_filename, std::string &cif_file, WFN &wave, bool debug, int accuracy, std::ofstream &file);

#include "wfn_class.h"

#endif
