#include "basis_set_helper.h"

const std::array<std::vector<primitive>, 118> read_basis_set(std::filesystem::path basis_path) {
	std::array<std::vector<primitive>, 118> basis_set; // creates an array composed of 118 vectors composed of primitives
	//open file
	std::ifstream file(basis_path);          // data
	std::string line, delimiter = ",";             // data delimiter

	std::getline(file, line); // first line read before beginning (H)

	for (int elem_idx = 0; elem_idx < 118; elem_idx++) {   // loop over all elements

		//std::cout << line << std::endl; // prints last line read before loop (visualization of elements) (shows last "omitted" line)

		while (true) {                   // loop over all primitives
			std::getline(file, line);   // read line
			if (line.empty()) break;    // delimiter between elements - empty lines

			std::vector<double> res = split_string<double>(line, delimiter);    //splits line into components using delimiter

			//basis_set[elem_idx].push_back(tmp);                 //introduce primitive inside of the according vectors of the array
			basis_set[elem_idx].push_back({ static_cast<int>(res[0]),
											static_cast<int>(res[1]),
											 res[2],res[3] });
			//std::cout << line << std::endl;                     //print every data line as it is processes (no lines lost??)
		}

		std::getline(file, line);                               //after while-break, one line omitted (element line)

	}
	file.close();
	return basis_set;
}


//const std::unordered_map<std::string, BasisSet> built_in_basis_sets = {
//	{ "sto-3g", /* ... actual BasisSet value */ },
//	{ "6-31g", /* ... */ }
//};

void write_basis_sets(const std::filesystem::path basis_dir, const std::unordered_map<std::string, std::array<std::vector<primitive>, 118>> basis_sets) {
    std::ofstream file(basis_dir / "auxiliary_basis.hpp");
	//write a hpp file that contains the basis sets as unordered map
    file << "const std::unordered_map<std::string, BasisSet> built_in_basis_sets = {\n";
    for (const auto& [name, basis_set] : basis_sets) {
        file << "    { \"" << name << "\", ";
        file << "{ ";
        for (int i = 0; i < 118; ++i) {
            file << "{ ";
            for (const auto& primitive : basis_set[i]) {
                file << "{" << primitive.get_center() << ", " << primitive.get_type() << ", " << primitive.get_exp() << ", " << primitive.get_coef() << "}, ";
            }
            file << "},\n ";
        }
        file << "} },\n";
    }
    file << "};\n";

    file.close();

}

