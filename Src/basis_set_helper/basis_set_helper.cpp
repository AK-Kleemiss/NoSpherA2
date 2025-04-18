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


//const std::unordered_map<std::string, std::array<std::vector<primitive>, 118>> built_in_basis_sets = {
//	{ "basis name",
//		{{
//			{ {0, 0, 1.0, 0.5}, {0, 1, 2.0, 0.4} }, For lines containing data
//			{ },									For lines containing no data
//			{ {0, 0, 1.0, 0.5}, {0, 1, 2.0, 0.4} },
//			...
//		}}
//	}
//};

void write_basis_sets(const std::filesystem::path basis_dir, const std::unordered_map<std::string, std::array<std::vector<primitive>, 118>> basis_sets) {
    std::ofstream file(basis_dir / "auxiliary_basis.hpp");
	//write a hpp file that contains the basis sets as unordered map
	file << "#pragma once\n";
	file << "#define P(c, t, e, coef) primitive(c, t, e, coef)\n";
    file << "const std::unordered_map<std::string, std::array<std::vector<primitive>, 118>> built_in_basis_sets = {\n";
	for (const auto& [basis_name, basis_set] : basis_sets) {
		file << "\t{ \"" << basis_name << "\",\n";
		file << "\t\tstd::array<std::vector<primitive>, 118>{ {\n";
		for (int i = 0; i < 118; ++i) {
			file << "\t\t\t{ ";
			for (const auto& prim : basis_set[i]) {
				file << "P" << "(" << prim.get_center() << ", " << prim.get_type() << ", " << prim.get_exp() << ", " << prim.get_coef() << "), ";
			}
			file << "},\n";
		}
		file << "\t\t} }\n";
		file << "\t},\n";
	}
	file << "};\n";

    file.close();

}

