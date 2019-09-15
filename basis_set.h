#ifndef __BASIS_SET_H__
#define __BASIS_SET_H__

class WFN;
#include <string>

//bool read_basis_set(const std::string &basis_set_path, WFN &wave, bool debug);  THIS FUNCTION IS BUGGY!

bool read_basis_set_vanilla(const std::string &basis_set_path, WFN &wave, const bool &debug, bool manual);

bool read_basis_set_missing(const std::string &basis_set_path, WFN &wave, bool debug);

bool delete_basis_set_vanilla(const std::string &basis_set_path, WFN &wave, bool debug);

#include "wfn_class.h"

#endif
