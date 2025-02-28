#pragma once

class WFN;
#include <filesystem>

// bool read_basis_set(const std::string &basis_set_path, WFN &wave, bool debug);  THIS FUNCTION IS BUGGY!

bool read_basis_set_vanilla(const std::filesystem::path &basis_set_path, WFN &wave, const bool &debug);

bool read_basis_set_missing(const std::filesystem::path &basis_set_path, WFN &wave, bool debug);

#include "wfn_class.h"