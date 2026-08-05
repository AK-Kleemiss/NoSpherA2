#pragma once

#include <filesystem>
#include <iosfwd>

#include "tsc_block.h"

tsc_block<int, cdouble> read_tsc_table(
    const std::filesystem::path& table_file);

bool convert_tsc_ids_to_labels(
    const std::filesystem::path& table_file,
    const std::filesystem::path& cif_file,
    const std::filesystem::path& output_file,
    std::ostream& log);

bool write_label_tsc_from_block(
    const tsc_block<int, cdouble>& input,
    const std::filesystem::path& cif_file,
    const std::filesystem::path& output_file,
    std::ostream& log);
