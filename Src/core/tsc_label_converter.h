#pragma once

#include <filesystem>
#include <iosfwd>

bool convert_tsc_ids_to_labels(
    const std::filesystem::path& table_file,
    const std::filesystem::path& cif_file,
    const std::filesystem::path& output_file,
    std::ostream& log);
