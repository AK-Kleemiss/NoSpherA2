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

long double primitve_normalization(const double exp, const int l) {
    long double numerator = std::pow(2 * (long double) exp / constants::PI, 3.0 / 2.0) * std::pow(4 * (long double) exp, l);
    long double denominator = ((l == 0) ? 1 : constants::double_ft[2 * l - 1]);
    return std::sqrt(numerator / denominator);
}


//https://pubs.acs.org/doi/10.1021/acs.jctc.6b01041
//From https://pyscf.org/_modules/pyscf/df/addons.html
void BasisSet::gen_etb_aux(const WFN& orbital_wfn) {
    if (_ownedPrimitives.size() != 0) {
        std::cout << "Basis set already generated! Using this one!" << std::endl;
        return;
    }
    _name = "aux-etb-basis";

    ivec seen_elements;
    for (const atom& atm : orbital_wfn.get_atoms()) {
        const int nuc_charge = atm.get_charge();
        //Check if the element is already in the seen_elements vector
        if (std::find(seen_elements.begin(), seen_elements.end(), nuc_charge) != seen_elements.end()) continue;
        err_chekf(atm.get_basis_set().size() != 0, "Can not generate auto-aux! Orbital Basis for Element: " + std::to_string(nuc_charge) + " is not defined!", std::cout);
        seen_elements.emplace_back(nuc_charge);

        std::unordered_map<int, std::pair<double, double>> angular_momentum_exponents;  //Min, Max exponents per angular momentum

        int l_max = 0;

        for (const basis_set_entry& p : atm.get_basis_set()) {
            int l = p.get_type() - 1;
            l_max = std::max(l_max, l);

            double exp = p.get_exponent();
            double coef = p.get_coefficient() / primitve_normalization(exp, l);

            if (std::abs(coef) < 1e-3) continue;

            auto& [min_exp, max_exp] = angular_momentum_exponents[l];
            if (min_exp == 0.0 && max_exp == 0.0) min_exp = max_exp = exp;
            else {
                min_exp = std::min(min_exp, exp);
                max_exp = std::max(max_exp, exp);
            }
        }

        std::array<int, 4> configuration = constants::GROUND_STATE_CONFIGURATION[nuc_charge];
        int max_shells = 4 - std::count(configuration.begin(), configuration.end(), 0);
        int l_max_aux = std::min(l_max, max_shells) * 2;

        // Traverse upper triangle of matrix of exponents in a sideways fashion to determine the lowest and highest exponent per angular momentum
        // Example:
        // l=0: 1 2 3 4
        // l=1:   5 6 7
        // l=2:     8 9
        // l=3:       10
        // Then the algorithm will traverse the matrix like this:
        // 1, 2, 35, 46, 78, 9, 10
        ivec n_functions(l_max_aux + 1);
        std::unordered_map<int, std::pair<double, double>> angular_exponent_ranges;  //ranges to produce the auxiliary function from
        for (int ll = 0; ll <= l_max_aux; ll++) {
            int i = (ll <= l_max) ? 0 : ll - l_max;
            int j = (ll <= l_max) ? ll : l_max;

            auto& [min_exp, max_exp] = angular_exponent_ranges[ll];
            min_exp = std::numeric_limits<double>::max(), max_exp = std::numeric_limits<double>::min();
            while (i <= j) {
                min_exp = std::min(min_exp, (angular_momentum_exponents[i].first + angular_momentum_exponents[j].first));
                max_exp = std::max(max_exp, (angular_momentum_exponents[i].second + angular_momentum_exponents[j].second));
                i++, j--;
            }

            //Calculate number of neccecary exponents in aux basis
            n_functions[ll] = std::ceil(std::log((min_exp + max_exp) / min_exp) / std::log(_beta));
        }

        //Set index of new basis set
        set_count_for_element(nuc_charge - 1, std::accumulate(n_functions.begin(), n_functions.end(), 0));

        std::cout << "-Element: " << nuc_charge << " Max_Lam: " << l_max << " Max_aux_lam: " << l_max_aux << " Max_shells: " << max_shells << std::endl;
        //Generate exponents following:
        //e = e_min * \beta^{i-1} for i = 1 .. n
        int shell = 0;
        for (int lam = 0; lam < n_functions.size(); lam++) {
            for (int i = n_functions[lam] - 1; i >= 0; i--, shell++) {
                add_owned_primitive({ 0, lam, angular_exponent_ranges[lam].first * std::pow(_beta, i), 1.0, shell });
            }
            std::cout << "Element: " << atm.get_charge() << " Shell: " << shell << " Angular Momentum: " << lam << " Exponents: " << n_functions[lam] << std::endl;
        }
    }
}


//https://pyscf.org/_modules/pyscf/df/autoaux.html
//The AutoAux algorithm by ORCA
//JCTC, 13 (2016), 554
namespace auto_aux_constants {
    inline constexpr double BETA_SMALL = 2.0;
    // Index by total L (0..), fallback used if L exceeds size
    inline const std::vector<double> BETA_BIG = { 1.8, 2.0, 2.2, 2.2, 2.2, 2.3, 3.0, 3.0 };
    // Cap factors for compact L (0..2*l_val)
    inline const std::vector<double> F_LAUX = {20, 7.0, 4.0, 4.0, 3.5, 2.5, 2.0, 2.0 };
}

void BasisSet::gen_auto_aux(const WFN& orbital_wfn) {
    if (_ownedPrimitives.size() != 0) {
        std::cout << "Basis set already generated! Using this one!" << std::endl;
        return;
    }
    _name = "auto-aux-basis";

    const double COEF_CUTOFF = 1e-3;

    auto get_or = [](const std::vector<double>& v, size_t i, double fallback) -> double {
        return (i < v.size() ? v[i] : fallback);
        };

    ivec seen_elements;
    for (const atom& atm : orbital_wfn.get_atoms()) {
        const int Z = atm.get_charge();
        if (std::find(seen_elements.begin(), seen_elements.end(), Z) != seen_elements.end()) continue;
        err_chekf(atm.get_basis_set().size() != 0,
            "Can not generate auto-aux! Orbital Basis for Element: " + std::to_string(Z) + " is not defined!",
            std::cout);
        seen_elements.emplace_back(Z);

        // ---- Gather primitive-level per-l stats: min, max, and "effective" max ----
        int l_max = 0;
        // We’ll size after we learn l_max; first pass to detect it:
        for (const basis_set_entry& p : atm.get_basis_set()) {
            l_max = std::max(l_max, (int)(p.get_type() - 1));
        }
        const int l_max1 = l_max + 1;

        std::vector<double> a_min_by_l(l_max1, std::numeric_limits<double>::infinity());
        std::vector<double> a_max_by_l(l_max1, 0.0);
        std::vector<double> a_eff_by_l(l_max1, 0.0); // can be tuned; we use max of "effective" primitives

        for (const basis_set_entry& p : atm.get_basis_set()) {
            const int l = p.get_type() - 1;
            const double exp = p.get_exponent();
            const double coef = p.get_coefficient() / primitve_normalization(exp, l);

            // Keep primitive-level stats; treat |coef|<cutoff as ineffective for mins/maxes
            if (std::abs(coef) < COEF_CUTOFF) continue;

            a_min_by_l[l] = std::min(a_min_by_l[l], exp);
            a_max_by_l[l] = std::max(a_max_by_l[l], exp);
            a_eff_by_l[l] = std::max(a_eff_by_l[l], exp); // "effective" max
        }
        // If any l had no effective primitives, fall back to any primitive present
        for (int l = 0; l < l_max1; ++l) {
            if (!std::isfinite(a_min_by_l[l]) || a_max_by_l[l] == 0.0) {
                // second pass: accept all primitives for this l as fallback
                for (const basis_set_entry& p : atm.get_basis_set()) {
                    if (p.get_type() - 1 != l) continue;
                    const double exp = p.get_exponent();
                    a_min_by_l[l] = std::min(a_min_by_l[l], exp);
                    a_max_by_l[l] = std::max(a_max_by_l[l], exp);
                    a_eff_by_l[l] = std::max(a_eff_by_l[l], exp);
                }
            }
            // Harden against pathological input
            if (!std::isfinite(a_min_by_l[l]) || a_min_by_l[l] <= 0.0) a_min_by_l[l] = 1.0;
            if (a_max_by_l[l] <= 0.0) a_max_by_l[l] = a_min_by_l[l];
            if (a_eff_by_l[l] <= 0.0) a_eff_by_l[l] = a_max_by_l[l];
        }

        // ---- Determine L_max for aux following the Python heuristic ----
        int l_val = 0;
        if (Z <= 2) l_val = 0;
        else if (Z <= 20) l_val = 1;
        else if (Z <= 56) l_val = 2;
        else              l_val = 3;

        int l_inc = (Z > 18 ? 2 : 1);
        int l_max_aux = std::min(std::max(l_val * 2, l_max + l_inc), l_max * 2);

        // ---- Build pairwise ranges using the triangle coupling rule |i-j| <= L <= i+j ----
        std::vector<double> a_min_by_L(l_max_aux + 1, std::numeric_limits<double>::infinity());
        std::vector<double> a_max_by_L(l_max_aux + 1, 0.0);
        std::vector<double> a_aux_by_L(l_max_aux + 1, 0.0);

        auto consider_pair = [&](int i, int j) {
            const double emin_ij = a_min_by_l[i] + a_min_by_l[j];
            const double emax_ij = a_max_by_l[i] + a_max_by_l[j];
            const double eaux_ij = a_eff_by_l[i] + a_eff_by_l[j];
            const int Lmin = std::abs(i - j);
            const int Lmax = i + j;
            for (int L = Lmin; L <= Lmax && L <= l_max_aux; ++L) {
                a_min_by_L[L] = std::min(a_min_by_L[L], emin_ij);
                a_max_by_L[L] = std::max(a_max_by_L[L], emax_ij);
                a_aux_by_L[L] = std::max(a_aux_by_L[L], eaux_ij);
            }
            };

        for (int i = 0; i <= l_max; ++i)
            for (int j = 0; j <= l_max; ++j)
                consider_pair(i, j);

        // Fill any untouched L (shouldn’t happen often) with conservative fallbacks
        for (int L = 0; L <= l_max_aux; ++L) {
            if (!std::isfinite(a_min_by_L[L])) a_min_by_L[L] = *std::min_element(a_min_by_l.begin(), a_min_by_l.end());
            if (a_max_by_L[L] <= 0.0)          a_max_by_L[L] = *std::max_element(a_max_by_l.begin(), a_max_by_l.end()) * 2.0;
            if (a_aux_by_L[L] <= 0.0)          a_aux_by_L[L] = a_max_by_L[L];
        }

        // ---- Cap compact regions with F_LAUX, per _auto_aux_element ----
        std::vector<double> a_max_adjust(l_max_aux + 1, 0.0);
        const int L_compact_cap = std::min(l_val * 2, l_max_aux);
        for (int L = 0; L <= L_compact_cap; ++L) {
            const double cap = get_or(auto_aux_constants::F_LAUX, (size_t)L, 1.0);
            a_max_adjust[L] = std::min(cap * a_aux_by_L[L], a_max_by_L[L]);
        }
        for (int L = L_compact_cap + 1; L <= l_max_aux; ++L) {
            a_max_adjust[L] = a_aux_by_L[L];
        }

        // ---- Decide n(L) and betas (BETA_SMALL for low L, BETA_BIG[L] for high L) ----
        std::vector<int> n_functions(l_max_aux + 1, 0);
        std::vector<double> beta_by_L(l_max_aux + 1, auto_aux_constants::BETA_SMALL);

        for (int L = 0; L <= l_max_aux; ++L) {
            const bool use_small = (L <= L_compact_cap);
            const double beta = use_small ? auto_aux_constants::BETA_SMALL
                : get_or(auto_aux_constants::BETA_BIG, (size_t)L, auto_aux_constants::BETA_SMALL);
            beta_by_L[L] = beta;

            const double emin = std::max(1e-16, a_min_by_L[L]);
            const double emax = std::max(emin * (1.0 + 1e-12), a_max_adjust[L]);

            // +1 guard so the largest generated exponent exceeds the target max
            const double ratio = emax / emin;
            int n = (ratio <= 1.0) ? 1 : static_cast<int>(std::ceil(std::log(ratio) / std::log(beta))) + 1;
            n_functions[L] = std::max(1, n);
        }

        // ---- Allocate and emit primitives ----
        const int n_total = std::accumulate(n_functions.begin(), n_functions.end(), 0);
        set_count_for_element(Z - 1, n_total);

        std::cout << "-Element: " << Z
            << " Max_l: " << l_max
            << " Max_aux_l: " << l_max_aux
            << " (val l=" << l_val << ", inc=" << l_inc << ")\n";

        int shell = 0;
        for (int L = 0; L <= l_max_aux; ++L) {
            const int nL = n_functions[L];
            const double emin = a_min_by_L[L];
            const double beta = beta_by_L[L];

            for (int i = nL - 1; i >= 0; --i, ++shell) {
                const double exp = emin * std::pow(beta, static_cast<double>(i));
                add_owned_primitive({ 0, L, exp, 1.0, shell });
            }
            std::cout << "Element: " << Z
                << " Shell: " << shell
                << " L: " << L
                << " Exponents: " << nL
                << " beta: " << beta << std::endl;
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



