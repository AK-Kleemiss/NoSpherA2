#include "basis_set_helper.h"

const std::array<std::vector<primitive>, 118> read_basis_set(std::filesystem::path basis_path) {
    std::array<std::vector<primitive>, 118> basis_set; // creates an array composed of 118 vectors composed of primitives
    //open file
    std::ifstream file(basis_path);          // data
    std::string line, delimiter = ",";             // data delimiter

    std::getline(file, line); // first line read before beginning (H)

    for (int elem_idx = 0; elem_idx < 118; elem_idx++) {   // loop over all elements

        //std::cout << line << std::endl; // prints last line read before loop (visualization of elements) (shows last "omitted" line)
        int basis_idx = -1;
        while (true) {                   // loop over all primitives
            std::getline(file, line);   // read line
            if (line.empty()) break;    // delimiter between elements - empty lines

            std::vector<double> res = split_string<double>(line, delimiter);    //splits line into components using delimiter
            if (basis_idx == static_cast<int>(res[4])) {
                basis_set[elem_idx][basis_set[elem_idx].size() - 1].exp.emplace_back(res[2]);
                basis_set[elem_idx][basis_set[elem_idx].size() - 1].coefficient.emplace_back(res[3]);
            }
            else {
                basis_set[elem_idx].emplace_back(
                    static_cast<int>(res[0]),
                    static_cast<int>(res[1]),
                    res[2], 
                    res[3]);
            }
        
            basis_idx = static_cast<int>(res[4]);
        }

        std::getline(file, line);                               //after while-break, one line omitted (element line)

    }
    file.close();
    return basis_set;
}


//constexpr Primitive sto3g_primitives[] = {
//    // H
//    {0, 13.01, 0.019685, 0},
//    {0, 1.962, 0.137977, 0},
//    {0, 0.4446, 0.478148, 0},
//    // He
//    {0, 38.36, 0.023766, 0},
//    {0, 5.77, 0.154679, 0},
//    {0, 1.24, 0.469630, 0}
//};
//
//constexpr ElementRange sto3g_ranges[118] = {
//    {0, 3},  // H
//    {3, 3},  // He
//    // Remaining are {0, 0}
//};
//
//constexpr BasisSetMetadata sto3g_metadata = {
//    "STO-3G",
//    sto3g_primitives,
//    std::size(sto3g_primitives),
//    sto3g_ranges
//};

std::string write_basis_set(std::ofstream& file, const std::string basis_name, std::array<std::vector<primitive>, 118> basis_set) {
    //replace - with _ and delete ( and ) 
    std::string basis_name_copy = basis_name;
    std::replace(basis_name_copy.begin(), basis_name_copy.end(), '-', '_');
    std::replace(basis_name_copy.begin(), basis_name_copy.end(), '(', '_');
    std::replace(basis_name_copy.begin(), basis_name_copy.end(), ')', '_');
    std::replace(basis_name_copy.begin(), basis_name_copy.end(), '+', 'P');

    //If basis_name_copy starts with numeric character, add a prefix
    if (std::isdigit(basis_name_copy[0])) {
        basis_name_copy = "basis_" + basis_name_copy;
    }

    file << "constexpr SimplePrimitive " << basis_name_copy << "_primitives[] = {\n";
    for (int i = 0; i < 118; i++) {
        int func = 0;
        for (const auto& p : basis_set[i]) {
            for (int j = 0; j < p.exp.size(); j++) {
                file << "    {" << p.center << ", " << p.type << ", " << p.exp[j] << ", " << p.coefficient[j] << ", " << func << "},\n";
            }
            func++;
        }
    }
    file << "};\n\n";
    file << "constexpr std::array<int, 118> " << basis_name_copy << "_counts {{";
    int running_idx = 0;
    for (int i = 0; i < 118; i++) {
        if (basis_set[i].size() == 0) {
            file << "0,";
            continue;
        }
        unsigned int s = 0;
        for (int j = 0; j < basis_set[i].size(); j++) {
            s += basis_set[i][j].exp.size();
        }
        file << "\n    " << s << ",";
        running_idx += s;
    }
    file << "}};\n\n";

    file << "constexpr std::array<int, 118> " << basis_name_copy << "_offsets = " << "compute_prefix_sum( " <<  basis_name_copy << "_counts );\n\n";

    file << "constexpr BasisSetMetadata " << basis_name_copy << "_metadata = {\n";
    file << "    \"" << basis_name << "\",\n";
    file << "    " << basis_name_copy << "_primitives,\n";
    file << "    std::size(" << basis_name_copy << "_primitives),\n";
    file << "    " << basis_name_copy << "_counts,\n";
    file << "    " << basis_name_copy << "_offsets\n";
    file << "};\n\n";
    return basis_name_copy;
}

