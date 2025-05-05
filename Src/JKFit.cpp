#include "pch.h"
#include "JKFit.h"
#include "integration_params.h"

std::shared_ptr<std::array<std::vector<primitive>, 118>> BasisSet::get_data() {
    if (_convertedData[0].size() == 0) {
        int prim_ctr = 0;
        for (int i = 0; i < 118; i++) {
            _convertedData[i].resize(_elementCounts[i]);
            for (int j = 0; j < _elementCounts[i]; j++, prim_ctr++) {
                _convertedData[i][j] = primitive(_primitives[prim_ctr]);
            }
        }
    }
    return std::make_shared<std::array<std::vector<primitive>, 118>>(_convertedData);
};

const std::span<const SimplePrimitive> BasisSet::operator[](const int& element) const {
	err_checkf(_elementCounts[element] != 0, "Basis set for element " + std::to_string(element) + " is not defined!", std::cout);
    return { _primitives + _elementOffsets[element], static_cast<std::size_t>(_elementCounts[element])};
}

double primitve_normalization(const double exp, const int l) {
    double numerator = std::pow(2 * exp / constants::PI, 3.0 / 2.0) * std::pow(4 * exp, l);
    double denominator = ((l == 0) ? 1 : constants::double_ft[2 * l - 1]);
    return std::sqrt(numerator / denominator);
}

void BasisSet::gen_aux(const WFN& orbital_wfn) {
    if (_ownedPrimitives.size() != 0) {
        std::cout << "Basis set already generated! Using this one!" << std::endl;
        return;
    }
    _name = "aux-basis";

    ivec seen_elements;
    for (const atom& atm : orbital_wfn.get_atoms()) {
        const int nuc_charge = atm.get_charge();
        //Check if the element is already in the seen_elements vector
        if (std::find(seen_elements.begin(), seen_elements.end(), nuc_charge) != seen_elements.end()) continue;
        err_chekf(atm.get_basis_set().size() != 0, "Can not generate auto-aux! Orbital Basis for Element: " + std::to_string(nuc_charge) + " is not defined!", std::cout);
        seen_elements.push_back(nuc_charge);

        std::array<int, 4> configuration = constants::GROUND_STATE_CONFIGURATION[nuc_charge];
        int max_shells = 4 - std::count(configuration.begin(), configuration.end(), 0);
        std::unordered_map<int, std::pair<double, double>> angular_momentum_exponents;  //Min, Max exponents per angular momentum

        int l_max = 0;

        for (const basis_set_entry& p : atm.get_basis_set()) {
            int l = p.get_type() - 1;


            double exp = p.get_exponent();
            double coef = p.get_coefficient() / primitve_normalization(exp, l);

            // Skip if the angular momentum is not in the range of 0 to max_shells	
            if (l >= max_shells + 1) continue;
            if (std::abs(coef) < 1e-3) continue;

            l_max = std::max(l_max, l);
            auto& [min_exp, max_exp] = angular_momentum_exponents[l];
            if (min_exp == 0.0 && max_exp == 0.0) min_exp = max_exp = exp;
            else {
                min_exp = std::min(min_exp, exp);
                max_exp = std::max(max_exp, exp);
            }
        }

        // Traverse upper triangle of matrix of exponents in a sideways fashion to determine the lowest and highest exponent per angular momentum
        // Example:
        // l=0: 1 2 3 4
        // l=1:   5 6 7
        // l=2:     8 9
        // l=3:       10
        // Then the algorithm will traverse the matrix like this:
        // 1, 2, 35, 46, 78, 9, 10
        ivec n_functions(l_max * 2 + 1);
        std::unordered_map<int, std::pair<double, double>> angular_exponent_ranges;  //ranges to produce the auxiliary function from
        for (int ll = 0; ll <= l_max * 2; ll++) {
            int i = (ll <= l_max) ? 0 : ll - l_max;
            int j = (ll <= l_max) ? ll : l_max;
            auto& [min_exp, max_exp] = angular_exponent_ranges[ll];
            min_exp = std::numeric_limits<double>::max(), max_exp = std::numeric_limits<double>::min();
            while (i <= j) {
                min_exp = std::min(min_exp, std::sqrt(angular_momentum_exponents[i].first * angular_momentum_exponents[j].first));
                max_exp = std::max(max_exp, std::sqrt(angular_momentum_exponents[i].second * angular_momentum_exponents[j].second));
                i++, j--;
            }
            min_exp *= 2, max_exp *= 2;

            //Calculate number of neccecary exponents in aux basis
            n_functions[ll] = std::ceil(std::log((min_exp + max_exp) / min_exp) / std::log(_beta));
        }

        //Set index of new basis set
        set_count_for_element(nuc_charge - 1, std::accumulate(n_functions.begin(), n_functions.end(), 0));

        //Generate exponents following:
        //e = e_min * \beta^{i-1} for i = 1 .. n
        int shell = 0;
        for (int lam = 0; lam < n_functions.size(); lam++) {
            for (int i = n_functions[lam] - 1; i >= 0; i--, shell++) {
                add_owned_primitive({ 0, lam, angular_exponent_ranges[lam].first * std::pow(_beta, i), 1.0, shell });
            }
        }
    }
}

std::shared_ptr<BasisSet> BasisSetLibrary::get_basis_set(std::string basis_name) {
    //Check if the supplied basis name is contained in part of a given basis set name
    std::string found_basis = "";
    if (basis_name == "") {
        std::cout << "No Basis Name Supplied! Aborting!!!" << std::endl;
        exit(1);
    }
    std::replace(basis_name.begin(), basis_name.end(), '_', '-');
    //Cast basis_name to lowercase
    std::transform(basis_name.begin(), basis_name.end(), basis_name.begin(), ::tolower);

    //Check if the supplied basis name is one of the precompiled basis sets
    int selected_idx = 0;
    for (; selected_idx < aux_basis_set_count; selected_idx++) {
        std::string_view name = aux_basis_sets[selected_idx].name;
        if (aux_basis_sets[selected_idx].name.find(basis_name) != std::string::npos) {
            found_basis = aux_basis_sets[selected_idx].name;
            break;
        }
    }
    err_checkf(found_basis != "", "Basis set " + basis_name + " not defined in BasisSetLibrary!", std::cout);
    return std::make_shared<BasisSet>(aux_basis_sets[selected_idx]);
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



