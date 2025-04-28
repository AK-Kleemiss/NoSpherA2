#include "pch.h"
#include "JKFit.h"


const std::array<std::vector<primitive>, 118>& BasisSet::get_data() {
	if (_converted_data[0].size() == 0) {
		for (int i = 0; i < 118; i++) {
			_converted_data[i].resize(_elementRanges[i].count);
			for (int j = 0; j < _elementRanges[i].count; j++) {
				_converted_data[i][j] = primitive(_primitives[_elementRanges[i].start + j]);
			}
		}
	}
	return _converted_data;
};

const std::span<const SimplePrimitive> BasisSet::operator[](int element) const {
	err_checkf(_elementRanges[element].count != 0, "Basis set for element " + std::to_string(element) + " is not defined!", std::cout);
	return { _primitives + _elementRanges[element].start, _elementRanges[element].count };
}


BasisSet BasisSetLibrary::get_basis_set(std::string basis_name) {
	//Check if the supplied basis name is contained in part of a given basis set name
	std::string found_basis = "";
	if (basis_name == "") {
		std::cout << "No Basis Name Supplied! Aborting!!!" << std::endl;
		exit(1);
	}
	std::replace(basis_name.begin(), basis_name.end(), '_', '-');
	//Cast basis_name to lowercase
	std::transform(basis_name.begin(), basis_name.end(), basis_name.begin(), ::tolower);

    int selected_idx = 0;
	for (; selected_idx < aux_basis_set_count; selected_idx++) {
		std::string_view name = aux_basis_sets[selected_idx].name;
        if (aux_basis_sets[selected_idx].name.find(basis_name) != std::string::npos) {
            found_basis = aux_basis_sets[selected_idx].name;
            break;
        }
    }
	err_checkf(found_basis != "", "Basis set " + basis_name + " not defined in BasisSetLibrary!", std::cout);
	std::cout << "Using basis set: " << found_basis << std::endl;
    return aux_basis_sets[selected_idx];
}


bool BasisSetLibrary::check_basis_set_exists(std::string basis_name) {
	//Check if the supplied basis name is in the basisSets map
    bool found_basis = false;
    std::replace(basis_name.begin(), basis_name.end(), '_', '-');
	for (int basis_set_idx = 0; basis_set_idx < aux_basis_set_count; basis_set_idx++) {
		if (aux_basis_sets[basis_set_idx].name.find(basis_name) != std::string::npos) {
            found_basis = true;
			break;
		}
	}
    return found_basis;
}


////Implementation following https://pubs.acs.org/doi/10.1021/acs.jctc.4c01555
//// This function generates even-tempered auxiliary basis sets for the given regular basis set
//// Returns a new BasisSet object with the auxiliary basis sets
//BasisSet BasisSetLibrary::gen_aux(const WFN& orbital_wfn, double& n) {
//	const int MAXLAUX = 6;
//	if (n >= MAXLAUX-1) {
//        std::cout << "Invalid auxiliary basis set choice: n = " << n << " must be smaller than " << MAXLAUX - 1 << ". Setting n to: " << MAXLAUX - 2 << std::endl;
//        n = MAXLAUX - 2;
//	}
//
//	BasisSet aux_basis;
//	aux_basis.set_name("aux-basis");
//
//	//atomic_numbers of all seen elements
//	ivec seen_elements;
//	for (const atom& atm : orbital_wfn.get_atoms()){
//        //Check if the element is already in the seen_elements vector
//		if (std::find(seen_elements.begin(), seen_elements.end(), atm.get_charge()) != seen_elements.end()) continue;
//        seen_elements.push_back(atm.get_charge());
//
//		std::pair<double, double> overall_min_max_exponents = { std::numeric_limits<double>::max() , std::numeric_limits<double>::min() };
//        std::pair<int, int> min_max_angular_momentum = { std::numeric_limits<int>::max() , std::numeric_limits<int>::min() };
//		//Find max and min angular momentum(-1 everywhere because orbital basis sets are stored as 1 = s, 2 = p, 3 = d ...)
//		std::unordered_map<int, std::pair<double, double>> angular_momentum_exponents;
//
//		for (const basis_set_entry& p : atm.get_basis_set()) {
//			int l = p.get_type() - 1;
//			double exp = p.get_exponent();
//
//			overall_min_max_exponents.first = std::min(overall_min_max_exponents.first, exp);
//            overall_min_max_exponents.second = std::max(overall_min_max_exponents.second, exp);
//
//			// Update min/max angular momentum
//			min_max_angular_momentum.first = std::min(min_max_angular_momentum.first, l);
//			min_max_angular_momentum.second = std::max(min_max_angular_momentum.second, l);
//
//			// Update exponent range
//			auto& [min_exp, max_exp] = angular_momentum_exponents[l];
//			if (min_exp == 0.0 && max_exp == 0.0) {
//				min_exp = max_exp = exp;
//			}
//			else {
//				min_exp = std::min(min_exp, exp);
//				max_exp = std::max(max_exp, exp);
//			}
//		}
//		int l_max = std::min(2 * min_max_angular_momentum.second, MAXLAUX), beta = MAXLAUX -n, n_a = 1;
//		
//        //Find upper and lower bounds for the auxiliary function
//        double a_0 = 2 * overall_min_max_exponents.first; //start auxiliary function exponent
//		double a_up = a_0, a_down = a_0;
//		while (a_up < 2 * overall_min_max_exponents.second) {
//            a_up *= beta;
//            n_a++;
//		}
//        while (a_down > overall_min_max_exponents.first) {
//            a_down /= beta;
//            n_a++;
//        }
//
//        //Calculate the auxiliary function exponents
//		//To do so, we initialize the auxiliary function set angularity to zero, Lset = 0 i.e. to a s-type auxiliary function set.
//		//The next loop, tests if the current processed auxiliary function exponent lies in a range of OBS product exponents, a_min and a_max, with higher angular momentum index.
//		//If this is true, the auxiliary function set angular momentum index is increased to min(2l, 6) where 2l is the angular momentum index associated with the OBS product.
//		//After the angular momentum index loop, we add the set with exponent α and angularity Lset to the ABS and decrease the α value by dividing it by beta.
//		double a = a_up;
//        std::vector<std::pair<int, double>> data;
//		for (int set = 0; set < n_a; set++) {
//			int l_set = 0;
//			for (int l = min_max_angular_momentum.first; l <= min_max_angular_momentum.second; l++) {
//				int L = std::min(2 * l, MAXLAUX);
//				double a_min = (2 * angular_momentum_exponents[l].first) / beta;
//				double a_max = (2 * angular_momentum_exponents[l].second) * beta;
//				if (a_min < a && a < a_max) {
//					l_set = std::max(L, l_set);
//				}
//			}
//			data.push_back({ l_set, a });
//            a /= beta;
//		}
//        // Sort the data by angular momentum and exponent
//		std::sort(data.begin(), data.end(), [](const auto& a, const auto& b) {
//			if (a.first != b.first)
//				return a.first < b.first;
//			return a.second > b.second;
//			});
//
//        for (int i = 0; i < n_a; i++) aux_basis.add_owned_primitive({ 0, data[i].first, data[i].second, 1.0, i });
//
//		aux_basis.set_element_range_for_element(atm.get_charge()-1, aux_basis.get_owned_primitive_count() - n_a, n_a);
//	}
//
//	return aux_basis;
//}

BasisSet BasisSetLibrary::gen_aux(const WFN& orbital_wfn, double& n) {
	ivec seen_elements;
	for (const atom& atm : orbital_wfn.get_atoms()) {
        const int nuc_charge = atm.get_charge();
		//Check if the element is already in the seen_elements vector
		if (std::find(seen_elements.begin(), seen_elements.end(), nuc_charge) != seen_elements.end()) continue;
		seen_elements.push_back(nuc_charge);

		std::array<int, 4> configuration = constants::GROUND_STATE_CONFIGURATION[nuc_charge];
        int max_shells = 4 - std::count(configuration.begin(), configuration.end(), 0);
        std::unordered_map<int, std::pair<double, double>> angular_momentum_exponents;  //Min, Max exponents per angular momentum
		int l_max = 0;

		for (const basis_set_entry& p : atm.get_basis_set()) {
			int l = p.get_type() - 1;
            l_max = std::max(l_max, l);
			if (l >= max_shells + 1) continue;

            double exp = p.get_exponent();
            if (std::abs(p.get_coefficient()) < 1e-3) continue;

			auto& [min_exp, max_exp] = angular_momentum_exponents[l];
			if (min_exp == 0.0 && max_exp == 0.0) min_exp = max_exp = exp;
            else {
                min_exp = std::min(min_exp, exp);
                max_exp = std::max(max_exp, exp);
            }
		}
        const int l_max2_m1 = 2 * l_max - 1;


        std::unordered_map<int, std::pair<double, double>> angular_exponent_ranges;  //ranges to produce the auxiliary function from
		for (int ll = 0; ll < 4 * l_max - 3; ll++) {
            double min_exp = std::numeric_limits<double>::max(), max_exp = std::numeric_limits<double>::min();
			for (int i = 0; i < 2 * l_max - 1; i++) {
				const int j = ll - i;
				if ((0 <= i < n) && (0 <= j < n) && (i <= j)) {
                    min_exp = std::min(min_exp, std::sqrt(angular_momentum_exponents[i].first * angular_momentum_exponents[j].first));
                    max_exp = std::max(max_exp, std::sqrt(angular_momentum_exponents[i].second * angular_momentum_exponents[j].second));
				}
			}
			std::cout << "DOOF" << std::endl;
		}


	}

    return BasisSet();
}
