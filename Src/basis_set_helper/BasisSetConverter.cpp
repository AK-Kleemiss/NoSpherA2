// BasisSetConverter.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <filesystem>
#include "basis_set_helper.h"

int main(int argc, char** argv)
{
    std::ofstream log_file("test.log");
    log_file << "Starting BasisSetConverter..." << std::endl; 
    std::filesystem::path basis_path(argv[1]);
    basis_path = basis_path.parent_path().parent_path();
    basis_path /= std::filesystem::path("Src/basis_set_helper/basis_sets/");
    basis_path.make_preferred();
    log_file << "Basis path: " << basis_path << std::endl;

    //Find all files in the directory
    std::vector<std::filesystem::path> files;
    for (const auto& entry : std::filesystem::directory_iterator(basis_path))
    {
        if (entry.is_regular_file() && entry.path().extension() == ".csv")
        {
            files.push_back(entry.path());
            log_file << "Found file: " << entry.path() << std::endl;
        }
    }
    log_file << "Number of files found: " << files.size() << std::endl;


    std::unordered_map<std::string, std::array<std::vector<primitive>, 118>> basis_sets;
    //Print all files
    for (const auto& file : files)
    {
        std::string basis_name = file.stem().string();
        std::transform(basis_name.begin(), basis_name.end(), basis_name.begin(), ::tolower);
        std::array<std::vector<primitive>, 118> basis_set = read_basis_set(file);
        basis_sets[basis_name] = basis_set;
    }

    std::ofstream aux_file(basis_path / "auxiliary_basis.cpp");
    aux_file << "#include \"../../JKFit.h\"\n";
    std::vector<std::string> basis_names_internal;
    aux_file << "namespace {\n";
    for (const auto& [basis_name, basis_set] : basis_sets)
    {
        basis_names_internal.push_back(write_basis_set(aux_file, basis_name, basis_set));
    }

    aux_file << "}\n";
    aux_file << "constexpr BasisSetMetadata aux_basis_sets[] = {\n";
    for (const std::string& name : basis_names_internal) {
        aux_file << "\t" << name << "_metadata,\n";
    }
    aux_file << "};\n\n";
    aux_file << "constexpr std::size_t aux_basis_set_count = std::size(aux_basis_sets);\n";

    aux_file.close();
    log_file.close();
}