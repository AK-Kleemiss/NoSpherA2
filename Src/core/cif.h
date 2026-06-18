#include "pch.h"

class WFN;

void write_wfn_CIF(WFN& wavy, const std::filesystem::path& fileName);
void write_wfn_CIF(WFN& wavy, const std::filesystem::path& fileName, tsc_block<int, cdouble>& tsc, options& opt);
void write_wfn_CIF(std::vector<WFN>& wavy, const std::filesystem::path& fileName, tsc_block<int, cdouble>& tsc, options& opt);