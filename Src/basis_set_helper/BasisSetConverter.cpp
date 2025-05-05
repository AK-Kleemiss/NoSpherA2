// BasisSetConverter.cpp : This file contains the 'main' function. Program execution begins and ends there.

#include <iostream>
#include <fstream>
#include <filesystem>
#include "basis_set_helper.h"


const std::vector<std::filesystem::path> get_all_basis_set_paths(std::filesystem::path basis_path) {
    std::vector<std::filesystem::path> files;
    for (const auto& entry : std::filesystem::directory_iterator(basis_path))
    {
        if (entry.is_regular_file() && entry.path().extension() == ".csv")
        {
            files.push_back(entry.path());

        }
    }
    return files;
};

bool needs_rewrite(const std::filesystem::path& basis_path, const std::vector<std::filesystem::path>& files, std::ostream& log_file) {
    const auto aux_file = basis_path / "auxiliary_basis.cpp";
    const auto checkpoint_file_path = basis_path / "checkpoint.txt";

    // Check if the checkpoint file and auxiliary file exist
    if (!std::filesystem::exists(checkpoint_file_path) || !std::filesystem::exists(aux_file))
        return true;

    // Read checkpoint
    std::ifstream checkpoint_file(checkpoint_file_path);
    if (!checkpoint_file.is_open()) {
        log_file << "Failed to open checkpoint.txt, assuming rewrite necessary.\n";
        return true;
    }

    // Read the first line to get the number of files
    std::string line;
    if (!std::getline(checkpoint_file, line)) {
        log_file << "Empty checkpoint.txt, assuming rewrite necessary.\n";
        return true;
    }

    int nr_files = std::stoi(line.substr(line.find(':') + 1));
    if (nr_files != static_cast<int>(files.size())) {
        return true;
    }

    for (int i = 0; i < nr_files; ++i) {
        if (!std::getline(checkpoint_file, line)) {
            log_file << "Unexpected end of checkpoint.txt.\n";
            return true;
        }

        const auto sep = line.find(':');
        const std::string file_name = line.substr(0, sep);
        int expected_size = std::stoi(line.substr(sep + 1));

        const auto file_path = basis_path / file_name;

        //Check if file_path is in files
        if (std::find(files.begin(), files.end(), file_path) == files.end()) {
            log_file << "File from checkpoint not found in basis set directory: " << file_name << "\n";
            return true;
        }

        if (expected_size != static_cast<int>(std::filesystem::file_size(file_path))) {
            log_file << "File size mismatch: " << file_name << "            REWRITING!" << "\n";
            return true;
        }
    }

    return false;
}


void write_checkpoint_file(std::filesystem::path basis_path, std::vector<std::filesystem::path> files) {
    std::ofstream checkpoint_file(basis_path / "checkpoint.txt");
    checkpoint_file << "Nr Files:"<< files.size() << "\n";
    for (const auto& file : files)
    {
        checkpoint_file << file.filename().generic_string() << ":" << std::filesystem::file_size(file) << "\n";
    }
    checkpoint_file.close();
}



int main(int argc, char** argv)
{
    //--------------------Extract Directory form input----------------
    std::ofstream log_file("test.log");
    log_file << "Starting BasisSetConverter..." << std::endl; 
    std::filesystem::path basis_path(argv[1]);
    basis_path = basis_path.parent_path().parent_path();
    basis_path /= std::filesystem::path("Src/basis_set_helper/basis_sets/");
    basis_path.make_preferred();
    log_file << "Basis path: " << basis_path << std::endl;

    //--------------------Find all basis set files----------------
    const std::vector<std::filesystem::path> files = get_all_basis_set_paths(basis_path);
    log_file << "Number of files found: " << files.size() << std::endl;


    std::filesystem::path src_path = basis_path.parent_path().parent_path().parent_path();
    if (!needs_rewrite(src_path, files, log_file)) {
        log_file << "No need to rewrite auxiliary_basis.cpp, exiting..." << std::endl;
        return 0;
    }

    write_checkpoint_file(src_path, files);

    std::unordered_map<std::string, std::array<std::vector<primitive>, 118>> basis_sets;
    //Read all files and convert them to the new format
    for (const auto& file : files)
    {
        std::string basis_name = file.stem().string();
        std::transform(basis_name.begin(), basis_name.end(), basis_name.begin(), ::tolower);
        std::array<std::vector<primitive>, 118> basis_set = read_basis_set(file);
        basis_sets[basis_name] = basis_set;
    }


    //Write the basis sets to a file
    std::ofstream aux_file(src_path / "auxiliary_basis.cpp");
    aux_file << "#include \"JKFit.h\" \n";
    std::vector<std::string> basis_names_internal;
    aux_file << "namespace {\n";

    aux_file << "constexpr std::array<int, 118> compute_prefix_sum(const std::array<int, 118>& counts) {\n";
    aux_file << "    std::array<int, 118> offsets{};\n";
    aux_file << "    int sum = 0;\n";
    aux_file << "    for (int i = 0; i < 118; ++i) {\n";
    aux_file << "        offsets[i] = sum;\n";
    aux_file << "        sum += counts[i];\n";
    aux_file << "    }\n";
    aux_file << "    return offsets;\n";
    aux_file << "}\n";

    for (const auto& [basis_name, basis_set] : basis_sets)
    {
        basis_names_internal.push_back(write_basis_set(aux_file, basis_name, basis_set));
    }

    aux_file << "}\n";
    aux_file << "constexpr BasisSetMetadata aux_basis_sets[] = {\n";
    for (const std::string& name : basis_names_internal) {
        aux_file << "\t" << name << "_metadata,\n";
    }
    aux_file << "};\n";
    aux_file << "constexpr std::size_t aux_basis_set_count = std::size(aux_basis_sets);\n";

    aux_file.close();
    log_file.close();
}