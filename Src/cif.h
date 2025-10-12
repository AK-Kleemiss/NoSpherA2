#include "pch.h"

class WFN;

void write_wfn_CIF(WFN& wavy, const std::filesystem::path& fileName, std::string additional_info = "");
void write_wfn_CIF(std::vector<WFN>& wavy, const std::filesystem::path& fileName, std::string additional_info = "");