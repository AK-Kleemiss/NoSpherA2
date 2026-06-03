#include "pch.h"
#include "basis_set.h"
#include "convenience.h"
#include "libCintMain.h"

std::shared_ptr<std::array<std::vector<primitive>, 118>> BasisSet::get_data() {
    if (_convertedData[0].size() == 0) {
        for (int i = 0; i < 118; i++) {
            _convertedData[i].resize(_elementCounts[i]);
            int prim_offset = _elementOffsets[i];
            for (int j = 0; j < _elementCounts[i]; j++) {
                _convertedData[i][j] = primitive(_primitives[j + prim_offset]);
            }
        }
    }
    return std::make_shared<std::array<std::vector<primitive>, 118>>(_convertedData);
};

const std::span<const SimplePrimitive> BasisSet::operator[](const int& element) const {
    err_checkf(_elementCounts[element] != 0, "Basis set for element " + std::to_string(element) + " is not defined!", std::cout);
    return { _primitives + _elementOffsets[element], static_cast<std::size_t>(_elementCounts[element]) };
}

constexpr std::array<int, 118> compute_prefix_sum(const std::array<int, 118>& counts) {
    std::array<int, 118> offsets{};
    int sum = 0;
    for (int i = 0; i < 118; ++i) {
        offsets[i] = (counts[i] == 0) ? -1 : sum;
        sum += counts[i];
    }
    return offsets;
}
//Look through all primitives in the main basisset and see if they are defined in the other basis set, if so, add them to this basis set
void BasisSet::operator+=(const BasisSet& other) {
    _name += "_plus_" + other._name;
    std::vector<SimplePrimitive> new_primitives;
    for (int elem = 0; elem < 118; elem++) {
        //If this basis set has primitives for this element, add them
        int count_this = _elementCounts[elem];
        if (count_this != 0) {
            new_primitives.insert(new_primitives.end(), _primitives + _elementOffsets[elem], _primitives + _elementOffsets[elem] + count_this);
            continue;
        }

        //If not, look if the other basis set has primitives for this element, if so, add them
        int count_other = other._elementCounts[elem];
        if (count_other != 0) {
            new_primitives.insert(new_primitives.end(), other._primitives + other._elementOffsets[elem], other._primitives + other._elementOffsets[elem] + count_other);
            _elementCounts[elem] = count_other;
            continue;
        }
    }
    _ownedPrimitives = new_primitives;
    _primitiveCount = static_cast<uint32_t>(_ownedPrimitives.size());
    _elementOffsets = compute_prefix_sum(_elementCounts);
    _primitives = _ownedPrimitives.data();
}


occ::qm::AOBasis BasisSet::to_AOBasis(const std::vector<occ::core::Atom>& atoms) const {
    err_checkf(_primitiveCount != 0, "No basis data available in BasisSet Object, please load something first!", std::cout);
    std::vector<occ::gto::Shell> shells;
    std::vector<occ::gto::Shell> ecp_shells;

    std::vector<int> ecp_electrons(atoms.size(), 0);
    //int nsh_ecp = 0;

    for (size_t a = 0; a < atoms.size(); ++a) {
        std::array<double, 3> origin = { atoms[a].x, atoms[a].y, atoms[a].z };
        const std::size_t Z = atoms[a].atomic_number - 1;
        const std::span<const SimplePrimitive> primitives = this->operator[](Z);
        int primitive_idx = 0;
        for (int shell_nr = 0; shell_nr <= primitives.back().shell; shell_nr++) {
            int l = primitives[primitive_idx].type;
            std::vector<double> exponents;
            std::vector<double> coefficients;
            while (shell_nr == primitives[primitive_idx].shell) {
                exponents.push_back(primitives[primitive_idx].exp);
                coefficients.push_back(primitives[primitive_idx].coefficient);
                primitive_idx++;
            }
            occ::gto::Shell shell(l, exponents, { coefficients }, origin);
            shell.kind = occ::gto::Shell::Kind::Spherical;
            shell.incorporate_shell_norm();
            shells.push_back(shell);
        }
        // TODO: Handle ECPs if needed, currently not supported
        //if (element_basis.ecp_electrons > 0) {
        //    occ::log::debug("Setting ECPs on atom {}", a);
        //    if (element_basis.ecp_shells.size() < 1) {
        //        std::string errmsg =
        //            fmt::format("Element (z={}) in basis '{}' has ECP electrons but "
        //                "no defined ECP shells",
        //                Z, json_filepath);
        //        throw std::runtime_error(errmsg);
        //    }
        //    for (const auto& s : element_basis.ecp_shells) {
        //        // handle general contractions by splitting
        //        for (int i = 0; i < s.angular_momentum.size(); i++) {
        //            ecp_shells.push_back(Shell(s.angular_momentum[i], s.exponents,
        //                { s.coefficients[i] }, origin));
        //            const auto& n = s.r_exponents;
        //            auto& shell = ecp_shells[nsh_ecp];
        //            shell.ecp_r_exponents = Eigen::Map<const IVec>(n.data(), n.size());
        //            nsh_ecp++;
        //        }
        //    }
        //    ecp_electrons[a] = element_basis.ecp_electrons;
        //}
    }
    occ::qm::AOBasis result(atoms, shells, _name, ecp_shells);
    //result.set_ecp_electrons(ecp_electrons);
    return result;
}

//https://pyscf.org/_modules/pyscf/df/autoaux.html
//The AutoAux algorithm by ORCA
//JCTC, 13 (2016), 554
namespace auto_aux_constants {
    inline constexpr double BETA_SMALL = 1.8;
    // Index by total L (0..), fallback used if L exceeds size
    inline const std::vector<double> BETA_BIG = { 1.8, 2.0, 2.2, 2.2, 2.2, 2.3, 3.0, 3.0, 3.0, 3.0 };
    // Cap factors for compact L (0..2*l_val)
    inline const std::vector<double> F_LAUX = { 20, 7.0, 4.0, 4.0, 3.5, 2.5, 2.0, 2.0, 2.0, 2.0 };

    inline double gaussian_int(int n, double exp) {
        double n1 = (n + 1) * 0.5;
        return std::tgamma(n1) / (2.0 * pow(exp, n1));
    }
}

void BasisSet::gen_auto_aux(const WFN& orbital_wfn) {
    if (_ownedPrimitives.size() != 0) {
        std::cout << "Basis set already generated! Using this one!" << std::endl;
        return;
    }
    _name = "auto-aux-basis";

    ivec seen_elements;
    for (const atom& atm : orbital_wfn.get_atoms()) {
        const int Z = atm.get_charge();
        if (std::find(seen_elements.begin(), seen_elements.end(), Z) != seen_elements.end()) continue;
        err_chekf(atm.get_basis_set().size() != 0,
            "Can not generate auto-aux! Orbital Basis for Element: " + std::to_string(Z) + " is not defined!",
            std::cout);
        seen_elements.emplace_back(Z);

        gen_auto_aux_for_element(atm);

    }
    //SVD_prune_aux_basis(orbital_wfn);
}

/// Deduplicate candidates within tolerance on same atom
std::vector<double> dedup_exponents(std::vector<double> exps, double tol = 0.1) {
    std::sort(exps.begin(), exps.end(), std::greater<double>());
    std::vector<double> out;
    for (double e : exps) {
        bool dup = false;
        for (double x : out) {
            if (std::abs(e - x) / std::max(e, x) < tol) {
                dup = true;
                break;
            }
        }
        if (!dup) out.push_back(e);
    }
    return out;
}

/// Pivoted Cholesky decomposition for selecting linearly independent functions
std::vector<int> pivoted_cholesky(dMatrix2& A, double threshold) {
    const int n = A.extent(0);
    std::vector<int> pivot_idx;
    vec diag(n);
    for (int i = 0; i < n; ++i) {
        diag[i] = A(i, i);
    }
    dMatrix2 L(n, n);
    std::vector<int> perm(n);
    std::iota(perm.begin(), perm.end(), 0);

    for (int k = 0; k < n; ++k) {
        int max_idx = k;
        double max_val = diag[k];
        for (int i = k + 1; i < n; ++i) {
            if (diag[i] > max_val) {
                max_val = diag[i];
                max_idx = i;
            }
        }

        if (max_val < threshold)
            break;

        if (max_idx != k) {
            std::swap(perm[k], perm[max_idx]);
            std::swap(diag[k], diag[max_idx]);
            for (int j = 0; j < k; ++j) {
                std::swap(L(k, j), L(max_idx, j));
            }
        }

        pivot_idx.push_back(perm[k]);
        L(k, k) = std::sqrt(max_val);

        for (int i = k + 1; i < n; ++i) {
            int orig_i = perm[i];
            int orig_k = perm[k];
            double sum = 0.0;
            for (int j = 0; j < k; ++j) {
                sum += L(i, j) * L(k, j);
            }
            L(i, k) = (A(orig_i, orig_k) - sum) / L(k, k);

            double diag_sum = 0.0;
            for (int j = 0; j <= k; ++j) {
                diag_sum += L(i, j) * L(i, j);
            }
            diag[i] = std::max(0.0, A(orig_i, orig_i) - diag_sum);
        }
    }

    return pivot_idx;
}

std::vector<double> prune_element_candidates_for_L(
    int L,
    const std::vector<double>& exponents,
    double threshold)
{
    if (exponents.empty()) return {};

    std::vector<double> unique_exps = dedup_exponents(exponents, 0.1);
    const int n = static_cast<int>(unique_exps.size());
    if (n <= 1) return unique_exps;

    atom tmp_atom("H", "H", 1, 0.0, 0.0, 0.0, 1);
    for (int i = 0; i < unique_exps.size(); i++) {
        tmp_atom.push_back_basis_set(unique_exps[i], 1.0, L, i);
    }
    WFN tmp_wfn(e_origin::NOT_YET_DEFINED);
    tmp_wfn.push_back_atom(tmp_atom);

    Int_Params tmp_params(tmp_wfn);
    vec res;
    compute2C<Coulomb2C_SPH>(tmp_params, res);
    dMatrixRef2 V(res.data(), tmp_params.get_nao(), tmp_params.get_nao());

    const int funcs_per_shell = 2 * L + 1;
    dMatrix2 V_shell(n, n);

    for (int i = 0; i < n; ++i) {
        int i0 = i * funcs_per_shell;
        for (int j = 0; j < n; ++j) {
            int j0 = j * funcs_per_shell;
            double acc = 0.0;
            for (int m = 0; m < funcs_per_shell; ++m) {
                acc += V(i0 + m, j0 + m);
            }
            V_shell(i, j) = acc / funcs_per_shell;
        }
    }

    vec d(n);
    for (int i = 0; i < n; i++) {
        d[i] = V_shell(i, i);
    }

    for (int i = 0; i < n; ++i) {
        d[i] = std::sqrt(std::max(d[i], 1e-16));
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            V_shell(i, j) /= d[i] * d[j];
        }
    }

    auto keep = pivoted_cholesky(V_shell, threshold);

    std::vector<double> pruned;
    pruned.reserve(keep.size());
    for (int idx : keep) {
        pruned.push_back(unique_exps[idx]);
    }

    std::sort(pruned.begin(), pruned.end(), std::greater<double>());
    return pruned;
}

void BasisSet::gen_auto_aux_for_element(const atom& atm) {
    const int Z = atm.get_charge();
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

    int prim_idx = 0;
    for (unsigned int shell = 0; shell < atm.get_shellcount_size(); shell++)
    {
        const int shelltype = atm.get_basis_set_type(prim_idx) - 1;
        const int shellsize = atm.get_shellcount(shell);
        vec coefs(shellsize), exps(shellsize);
        for (int i = 0; i < shellsize; i++) {
            coefs[i] = atm.get_basis_set_coefficient(i + prim_idx);
            exps[i] = atm.get_basis_set_exponent(i + prim_idx);
            a_min_by_l[shelltype] = std::min(a_min_by_l[shelltype], exps[i]);
            a_max_by_l[shelltype] = std::max(a_max_by_l[shelltype], exps[i]);
        }
        prim_idx += shellsize;
        coefs = Int_Params::normalize_gto(coefs, exps, shelltype);  //TODO: What norm is better, we have to still figure out.... 
        //for (int i = 0; i < coefs.size(); i++)
        //{
        //    coefs[i] *= std::sqrt(constants::PI * 4 / constants::double_ft[2 * shelltype + 1]); // Conversion factor from GBW to libcint  ... something something, spherical harmonics...
        //}

        vec r_exp(shellsize, 0.0);
        for (int i = 0; i < shellsize; i++) {
            for (int j = 0; j < shellsize; j++) {
                r_exp[i] += coefs[i] * auto_aux_constants::gaussian_int(shelltype * 2 + 2, exps[i] + exps[j]) * coefs[j];
            }
        }
        const double k = (std::pow(2.0, (2 * shelltype + 1)) * constants::ft[shelltype + 1] * constants::ft[shelltype + 1]) /
            static_cast<double>(constants::ft[2 * shelltype + 2]);
        const double kk2 = 2.0 * k * k;
        for (int i = 0; i < shellsize; i++) {
            const double e_eff = (kk2) / (constants::PI * r_exp[i] * r_exp[i]);
            a_eff_by_l[shelltype] = std::max(e_eff, a_eff_by_l[shelltype]);
        }
    }
    vec2 a_min_prim(l_max1, vec(l_max1));
    vec2 a_max_prim(l_max1, vec(l_max1));
    vec2 a_max_aux(l_max1, vec(l_max1));

    for (int li = 0; li < l_max1; li++) {
        for (int lj = 0; lj < l_max1; lj++) {
            a_min_prim[li][lj] = a_min_by_l[li] + a_min_by_l[lj];
            a_max_prim[li][lj] = a_max_by_l[li] + a_max_by_l[lj];
            a_max_aux[li][lj] = a_eff_by_l[li] + a_eff_by_l[lj];
        }
    }
    int l_occ_max = 3;
    if (Z <= 2) {
        l_occ_max = 0;
    }
    else if (Z <= 18) {
        l_occ_max = 1;
    }
    else if (Z <= 54) {
        l_occ_max = 2;
    }
    int l_inc = 1;
    if (Z > 18) {
        l_inc = 2;
    }

    int l_keep_max = std::max(2 * l_occ_max, l_occ_max + l_max + l_inc);
    int l_max_aux = std::min(l_keep_max, 2 * l_max);

    vec a_min_by_l_aux(l_max_aux + 1);
    vec a_max_by_l_aux(l_max_aux + 1);
    vec a_aux_by_l_aux(l_max_aux + 1);
    for (int ll = 0; ll <= l_max_aux; ll++) {
        double min_val = std::numeric_limits<double>::infinity();
        double max_val = 0.0;
        double aux_val = 0.0;
        for (int li = 0; li < l_max1; li++) {
            for (int lj = 0; lj < l_max1; lj++) {
                if (std::abs(li - lj) <= ll && ll <= (li + lj)) {
                    min_val = std::min(min_val, a_min_prim[li][lj]);
                    max_val = std::max(max_val, a_max_prim[li][lj]);
                    aux_val = std::max(aux_val, a_max_aux[li][lj]);
                }
            }
        }
        a_min_by_l_aux[ll] = min_val;
        a_max_by_l_aux[ll] = max_val;
        a_aux_by_l_aux[ll] = aux_val;
    }
    vec a_max_adjusted(l_max_aux + 1);
    for (int l = 0; l <= l_occ_max * 2; ++l) {
        a_max_adjusted[l] = std::min(
            auto_aux_constants::F_LAUX[l] * a_aux_by_l_aux[l],
            a_max_by_l_aux[l]
        );
    }
    for (int l = l_occ_max * 2 + 1; l <= l_max_aux; ++l) {
        if (a_aux_by_l_aux[l] > 1e7) { //Small guard against rediculous exponents...
            a_max_adjusted[l] = a_max_by_l_aux[l];
        }
        else {
            a_max_adjusted[l] = a_aux_by_l_aux[l];
        }
    }

    int added_functions = 0;


    int shell = 0, n_funcs;
    double beta, ns;
    for (int l = 0; l <= l_max_aux; l++) {
        beta = (l < l_occ_max * 2 + 1) ? auto_aux_constants::BETA_SMALL : auto_aux_constants::BETA_BIG[l];
        ns = std::log(a_max_adjusted[l] / a_min_by_l_aux[l]) / std::log(beta);
        n_funcs = static_cast<int>(std::ceil(ns)) + 1;

        //vec candidate_exps(n_funcs);
        //for (int i = n_funcs - 1; i >= 0; --i) {
        //    candidate_exps[i] = a_min_by_l_aux[l] * std::pow(beta, i);
        //}
        //auto kept = prune_element_candidates_for_L(l, candidate_exps, 5e-5);
        //std::cout << "Element " << Z << ", L = " << l << " with beta = " << beta << ": Generated " << n_funcs << " candidates, kept " << kept.size() << " after pruning." << std::endl;
        //n_funcs = static_cast<int>(kept.size());
        //if (n_funcs <= 0) continue;
        //for (double exp : kept) {
        //    add_owned_primitive({ 0, l, exp, 1.0, shell });
        //    ++shell;
        //    ++added_functions;
        //}

        std::cout << "Beta for Element " << Z << " and l " << l << " : " << beta << " with " << n_funcs << " Functions." << std::endl;
        if (n_funcs <= 0) continue;
        for (int i = n_funcs - 1; i >= 0; --i, ++added_functions, shell++) {
            double exp = a_min_by_l_aux[l] * std::pow(beta, i);
            add_owned_primitive({ 0, l,
                exp,
                1.0, shell
                });
        }
    }
    set_count_for_element(Z - 1, added_functions);
    _elementOffsets[Z - 1] -= added_functions; //Small hack to set the correct offset
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
    for (; selected_idx < basis_set_count; selected_idx++) {
        std::string_view name = basis_sets[selected_idx].name;
        if (basis_sets[selected_idx].name.find(basis_name) != std::string::npos) {
            found_basis = basis_sets[selected_idx].name;
            break;
        }
    }
    err_checkf(found_basis != "", "Basis set " + basis_name + " not defined in BasisSetLibrary!", std::cout);
    return std::make_shared<BasisSet>(basis_sets[selected_idx]);
}


bool BasisSetLibrary::check_basis_set_exists(std::string basis_name) {
    //Check if the supplied basis name is in the basisSets map
    bool found_basis = false;
    std::replace(basis_name.begin(), basis_name.end(), '_', '-');
    std::transform(basis_name.begin(), basis_name.end(), basis_name.begin(), ::tolower);
    for (int basis_set_idx = 0; basis_set_idx < basis_set_count; basis_set_idx++) {
        if (basis_sets[basis_set_idx].name.find(basis_name) != std::string::npos) {
            found_basis = true;
            break;
        }
    }
    return found_basis;
}



int load_basis_into_WFN(WFN& wavy, std::shared_ptr<BasisSet> b, bool decontract)
{
    wavy.set_basis_set_ptr((*b).get_data());
    int nr_coefs = 0;
    for (int i = 0; i < wavy.get_ncen(); i++)
    {
        int current_charge = wavy.get_atom_charge(i) - 1;
        const std::span<const SimplePrimitive> basis = (*b)[current_charge];
        int size = (int)basis.size();
        if (decontract) {
            for (int e = 0; e < size; e++)
            {
                wavy.push_back_atom_basis_set(i, basis[e].exp, 1.0, basis[e].type, e);  //We decontract the basis, this makes it easier to handle in the long run.
                nr_coefs += 2 * basis[e].type + 1; //HIER WEITER MACHEN!!
            }
        }
        else {
            int shell = -1;
            for (int e = 0; e < size; e++)
            {
                if (basis[e].shell != shell)
                {
                    shell = basis[e].shell;
                    nr_coefs += 2 * basis[e].type + 1;
                }
                wavy.push_back_atom_basis_set(i, basis[e].exp, basis[e].coefficient, basis[e].type, basis[e].shell);
            }
        }
    }
    return nr_coefs;
}

// Supposed to complete setup of WFN from a basis by adding types, exponents and nex
void complete_WFN_basis(WFN& wavy)
{
    int nex_ = 0;
    vec exponents_;
    ivec types_;
    for (int a = 0; a < wavy.get_ncen(); a++) {
		atom atom_ = wavy.get_atom(a);
        for (int b = 0; b < atom_.get_basis_set_size(); b++) {
            basis_set_entry bf_ = wavy.get_atom_basis_set_entry(a, b);
            int temp_type = bf_.get_type();
            double temp_exp = bf_.get_exponent();
            int effective_type = 0;
            int end = 2 * temp_type + 1;
            switch (temp_type) {
            case(0):
                effective_type = 1;
                break;
            case(1):
                effective_type = 2;
                break;
            case(2):
                effective_type = 5;
                break;
            case(3):
                effective_type = 11;
                break;
            case(4):
                effective_type = 21;
                break;
            case(5):
                effective_type = 36;
                break;
            }
			for (int idx = 0; idx < end; idx++, nex_++, effective_type++) {
				exponents_.push_back(temp_exp);
                types_.push_back(effective_type);
			}
        }
    }
	wavy.set_exponents(exponents_);
    wavy.set_nex(nex_);
    wavy.set_types(types_);
}

void gen_missing_basis_auto_aux(const WFN& orbital_wfn, std::shared_ptr<BasisSet>& current_basis) {
    std::shared_ptr<BasisSet> missing_basis = std::make_shared<BasisSet>();
    for (const atom& atm : orbital_wfn.get_atoms()) {
        const int Z = atm.get_charge();
        if (!current_basis->has_element(Z) && !missing_basis->has_element(Z)) {
            missing_basis->gen_auto_aux_for_element(atm);
            std::cout << "Element " << Z << " was missing in the auxiliary basis set! Generated using auto-aux!" << std::endl;
        }
    }
    (*current_basis) += (*missing_basis);
}

bool basis_set_complete(const WFN& orbital_wfn, std::shared_ptr<BasisSet> aux_basis) {
    for (const atom& atm : orbital_wfn.get_atoms()) {
        if (!aux_basis->has_element(atm.get_charge())) return false;
    }
    return true;
}

WFN generate_aux_wfn(const WFN& orbital_wfn, std::vector<std::shared_ptr<BasisSet>>& aux_basis) {
    err_checkf(aux_basis.size() > 0, "Aux-Basis Vecor == 0, something went wrong... try calling -ri_fit and see what happens", std::cout);
    //Check every basis set in the vector, if it is empty, generate it using auto_aux
    for (int basis_nr = 0; basis_nr < aux_basis.size(); basis_nr++) { if ((*aux_basis[basis_nr]).get_primitive_count() == 0)(*aux_basis[basis_nr]).gen_auto_aux(orbital_wfn); }

    //If two or more basis sets are supplied, combine them into one
    std::shared_ptr<BasisSet> combined_aux_basis = aux_basis[0];
    for (int basis_nr = 1; basis_nr < aux_basis.size(); basis_nr++) { (*combined_aux_basis) += (*aux_basis[basis_nr]); }

    if (!basis_set_complete(orbital_wfn, combined_aux_basis)) {
        gen_missing_basis_auto_aux(orbital_wfn, combined_aux_basis);
    }

    WFN wavy_aux(e_origin::NOT_YET_DEFINED);
    wavy_aux.set_atoms(orbital_wfn.get_atoms());
    wavy_aux.set_ncen(orbital_wfn.get_ncen());
    wavy_aux.delete_basis_set();
    load_basis_into_WFN(wavy_aux, combined_aux_basis);

    return std::move(wavy_aux);
}

/**
 * Reads a basis set from a file and sets it in the given WFN object.
 * Does not perform basis set changes or checks.
 *
 * @param basis_set_path The path to the directory containing the basis set files.
 * @param wave The WFN object to set the basis set in.
 * @param debug A boolean flag indicating whether to enable debug mode.
 * @param manual A boolean flag indicating whether to manually input the basis set name.
 * @return Returns true if the basis set is successfully read and set, false otherwise.
 */
bool BasisSetLibrary::read_basis_set_vanilla(const std::filesystem::path& basis_set_path, WFN& wave, const bool& debug)
{
    using namespace std;
    string basis_set_name;
    std::filesystem::path temp_name;
    bool end = false;
    while (!end)
    {
        basis_set_name = wave.get_basis_set_name();
        // assemble basis set name and look if file exists
        temp_name = basis_set_path / basis_set_name;
        if (exists(temp_name))
        {
            if (debug)
               std::cout << "basis set is valid, continuing..." << endl;
            end = true;
        }
        else
        {
           std::cout << "sorry, could not find this basis set in the basis set directory specified in the programs.config file!" << endl;
            return false;
            //std::cout << "What is the name of the basis set in the directory: ";
            // cin >> basis_set_name;
        }
    }
    wave.set_basis_set_name(basis_set_name);
    if (debug)
       std::cout << "File of basis set to load: " << temp_name << endl;
    ifstream ifile(temp_name, ios::in);
    //  Looking for all the types of atoms we need to find
    svec elements_list;
    bool found = false;
    for (int i = 0; i < wave.get_ncen(); i++)
    {
        if (debug)
           std::cout << "i: " << i << endl;
        for (int j = 0; j < elements_list.size(); j++)
        {
            if (elements_list[j].compare(wave.get_atom_label(i)) == 0)
                found = true;
            if (debug)
               std::cout << "   j: " << j << " Atom label: " << wave.get_atom_label(i) << endl;
        }
        if (!found)
        {
            elements_list.emplace_back(wave.get_atom_label(i));
            if (debug)
               std::cout << "Added an atom which was not there yet! " << wave.get_atom_label(i) << "!" << endl;
        }
        found = false;
    }
    if (debug)
    {
       std::cout << "Number of elements in elements_list: " << elements_list.size() << endl;
       std::cout << "This is the elements list:" << endl;
        for (int l = 0; l < elements_list.size(); l++)
           std::cout << l << ": " << elements_list[l] << "," << endl;
    }
    // int found_counter = 0;
    for (int i = 0; i < elements_list.size(); i++)
    {
        if (debug)
           std::cout << "before: " << elements_list[i] << " " << i << endl;
        while (elements_list[i].find(" ") != -1)
        {
            elements_list[i].erase(elements_list[i].find(" "), 1);
        } // Cut out spaces
        elements_list[i].append(":");
        if (debug)
        {
           std::cout << "after: " << elements_list[i] << " " << i << endl;
        }
        // scan the tonto style basis set file for the entries we are looking or:
        string line;
        getline(ifile, line);
        int file_type = 0;
        // check if we support that type of basis set
        while (line.find("keys=") == -1 && !ifile.eof())
            getline(ifile, line);
        if (debug)
        {
           std::cout << "Line after looking for keys=: " << line << endl;
        }
        if (line.find("keys=") < line.size() && debug)
           std::cout << "Found keys=!" << endl;
        if (line.find("turbomole") < line.size())
        {
            file_type = 1;
            if (debug)
               std::cout << "This file is written in turbomole type!" << endl;
        }
        else if (line.find("gamess-us") < line.size())
        {
            file_type = 2;
            if (debug)
               std::cout << "This file is written in gamess-us type!" << endl;
        }
        else if (line.find("gaussian") < line.size())
        {
            file_type = 3;
            if (debug)
               std::cout << "This file is written in gaussian type!" << endl;
        }
        else if (line.find("CRYSTAL") < line.size())
        {
            file_type = 1;
            wave.set_d_f_switch(true);
            if (debug)
               std::cout << "This file is written in CRYSTAL type!" << endl;
        }
        else
        {
           std::cout << "This type of file is not supported, please provide another basis set!" << endl;
            return false;
        }
        if (ifile.eof())
        {
           std::cout << "Please provide a basis set in the turbomole, gaussian or gamess-us format compatible with tonto."
                << "Look at the example files \"examble.basis\" and \"examble2.basis\" in the wfn_cpp folder if you want to see how it has to look like" << endl;
            return false;
        }
        while (!(line.find(elements_list[i]) < line.size()) && !ifile.eof())
            getline(ifile, line);
        if (debug)
           std::cout << "line while search for " << elements_list[i] << " :" << line << endl;
        if (debug && line.find(elements_list[i]) != -1)
        {
           std::cout << "I found an entry i know from the element list!" << endl;
           std::cout << "The line is: " << line << endl;
        }
        if (ifile.eof())
        {
           std::cout << "Could not find the atom you were looking for in the basis set file... " << endl;
            return false;
        }
        unsigned int shell = 0;
        if (line.find("{") == -1)
        {
            getline(ifile, line);
            if (debug)
            {
               std::cout << "I read an additional line!" << endl;
            }
        }
        while (line.find("}") == string::npos && !ifile.eof())
        {
            getline(ifile, line);
            stringstream stream;
            stream << line;
            if (line.find("}") != string::npos)
                break;
            int count = 0;
            //int nr_exp = 0;
            char c_temp = '?';
            double temp_vals[2]{ 0, 0 };
            int dum = 0;
            if (file_type == 1)
            {
                stream >> count >> c_temp;
                if (debug)
                   std::cout << "count: " << count << " type: " << c_temp << endl;
            }
            else if (file_type == 2 || file_type == 3)
            {
                stream >> c_temp >> count;
                if (debug)
                   std::cout << "count: " << count << " type: " << c_temp << endl;
            }
            for (int j = 0; j < count; j++)
            {
                getline(ifile, line);
                if (debug)
                {
                   std::cout << "read the " << j << ". line: " << line << endl;
                }
                stringstream stream2;
                stream2 << line;
                if (file_type == 1)
                    stream2 >> temp_vals[0] >> temp_vals[1];
                else if (file_type == 2)
                    stream2 >> dum >> temp_vals[0] >> temp_vals[1];
                else if (file_type == 3)
                    stream2 >> temp_vals[0] >> temp_vals[1];
                // this is where i started copying
                for (int h = 0; h < wave.get_ncen(); h++)
                {
                    string temp_label;
                    temp_label = wave.get_atom_label(h);
                    while (temp_label.find(" ") != -1)
                        temp_label.erase(temp_label.find(" "), 1);
                    temp_label.append(":");
                    if (elements_list[i].find(temp_label) != -1)
                    {
                        int type = 0;
                        switch (c_temp)
                        {
                        case 's':
                        case 'S':
                            type = 1;
                            break;
                        case 'p':
                        case 'P':
                            type = 2;
                            break;
                        case 'd':
                        case 'D':
                            type = 3;
                            break;
                        case 'f':
                        case 'F':
                            type = 4;
                            break;
                        default:
                           std::cout << "Sorry, orbital types higher than f-type are not yet supported!" << endl;
                            return false;
                        } // end switch of types
                        if (!wave.push_back_atom_basis_set(h, temp_vals[0], temp_vals[1], type, shell))
                        {
                           std::cout << "ERROR while pushing back atoms basis set" << endl;
                        }
                        if (debug)
                           std::cout << "Pushing back on atom: " << h + 1 << " with coef: " << temp_vals[1]
                            << " and exp: " << temp_vals[0] << " and type " << type << endl;
                    } // end if(find atom_label + :
                }     // end for h = ncen
                //nr_exp++;
                if (debug)
                   std::cout << "recapitulation[" << j << "]... type: " << c_temp << " coef: " << temp_vals[0] << " exp: " << temp_vals[1] << endl;
                if (dum > count)
                {
                   std::cout << "this should not happen, lets stop before i do something silly!" << endl;
                    return false;
                }
            }
            shell++;
        } // end while line != }
        if (debug)
           std::cout << "I found }: " << line << endl;
        ifile.seekg(0);
    } // end for element_list.size()
    if (debug)
    {
       std::cout << "FINISHED WITH READING BASIS SET!" << endl;
    }
    ifile.close();
    return true;
};

/**
 * Reads a basis set from a specified file path and updates the basis set information in the given WFN object.
 * Only applies it to atoms wihtout basis set info.
 *
 * @param basis_set_path The file path of the basis set.
 * @param wave The WFN object to update with the basis set information.
 * @param debug Flag indicating whether to enable debug mode.
 * @return Returns true if the basis set was successfully read and updated, false otherwise.
 */
bool BasisSetLibrary::read_basis_set_missing(const std::filesystem::path &basis_set_path, WFN &wave, bool debug)
{
    using namespace std;
    filesystem::path temp_p;
    bool end = false;
    while (!end)
    {
        // assemble basis set name and look if file exists
        temp_p = basis_set_path / wave.get_basis_set_name();
        if (exists(temp_p))
        {
            if (debug)
               std::cout << "basis set is valid, continueing..." << endl;
            end = true;
        }
        else
        {
           std::cout << "sorry, could not find this basis set in the basis set directory specified in the programs.config file!" << endl;
            return false;
            // manual = true;
            //std::cout << "What is the name of the basis set in the directory: ";
            // cin >> basis_set_name;
        }
    }
    if (debug)
       std::cout << "File of basis set to load: " << temp_p << endl;
    ifstream ifile(temp_p, ios::in);
    //  Looking for all the types of atoms we need to find
    svec elements_list;
    bool found = false;
    for (int i = 0; i < wave.get_ncen(); i++)
    {
        if (debug)
           std::cout << "i: " << i << endl;
        for (int j = 0; j < elements_list.size(); j++)
        {
            if (elements_list[j].find(wave.get_atom_label(i)) != -1)
                found = true;
            if (debug)
               std::cout << "   j: " << j << endl;
        }
        if (!found)
        {
            elements_list.emplace_back(wave.get_atom_label(i));
            if (debug)
               std::cout << "Added an atom which was not there yet! " << wave.get_atom_label(i) << endl;
        }
        found = false;
    }
    if (debug)
    {
       std::cout << "Number of elements in elements_list: " << elements_list.size() << endl;
       std::cout << "This is the elements list:" << endl;
        for (int l = 0; l < elements_list.size(); l++)
           std::cout << l << ": " << elements_list[l] << endl;
    }
    // int found_counter = 0;
    for (int i = 0; i < elements_list.size(); i++)
    {
        if (debug)
           std::cout << "before: " << elements_list[i] << " " << i << endl;
        if (elements_list[i].find(" "))
            elements_list[i].erase(elements_list[i].find(" "), 1);
        elements_list[i].append(":");
        if (debug)
        {
           std::cout << "after: " << elements_list[i] << " " << i << endl;
        }
        // scan the tonto style basis set file for the entries we are looking or:
        string line;
        getline(ifile, line);
        int file_type = 0;
        // check if we support that type of basis set
        while (line.find("keys=") == -1 && !ifile.eof())
        {
            if (debug)
            {
               std::cout << "line.size of first line: " << line.size() << "line.find(\"keys=\"): " << line.find("keys=") << endl;
            }
            getline(ifile, line);
        }
        if (debug)
        {
           std::cout << "Line after looking for keys=: " << line << endl;
        }
        if (line.find("keys=") < line.size() && debug)
           std::cout << "Found keys=!" << endl;
        if (line.find("turbomole") < line.size())
        {
            file_type = 1;
            if (debug)
               std::cout << "This file is written in turbomole type!" << endl;
        }
        else if (line.find("gamess-us") < line.size())
        {
            file_type = 2;
            if (debug)
               std::cout << "This file is written in gamess-us type!" << endl;
        }
        else if (line.find("gaussian") < line.size())
        {
            file_type = 3;
            if (debug)
               std::cout << "This file is written in gaussian type!" << endl;
        }
        else
        {
           std::cout << "This type of file is not supported, please provide another basis set!" << endl;
            return false;
        }
        if (ifile.eof())
        {
           std::cout << "Please provide a basis set in the turbomole, gaussian or gamess-us format compatible with tonto."
                 << "Look at the example files \"examble.basis\" and \"examble2.basis\" in the wfn_cpp folder if you want to see how it has to look like" << endl;
            return false;
        }
        while (!(line.find(elements_list[i]) < line.size()) && !ifile.eof())
        {
            getline(ifile, line);
            if (debug)
               std::cout << "line while search for " << elements_list[i] << " :" << line << endl;
        }
        if (debug && line.find(elements_list[i]) != -1)
        {
           std::cout << "I found an entry i know from the element list!" << endl;
           std::cout << "The line is: " << line << endl;
        }
        if (ifile.eof())
        {
           std::cout << "Could not find the atom you were looking for in the basis set file... " << endl;
            return false;
        }
        unsigned int shell = 0;
        if (line.find("{") == -1)
        {
            getline(ifile, line);
            if (debug)
               std::cout << "I read an additional line!" << endl;
        }
        while (line.find("}") == -1 && !ifile.eof())
        {
            getline(ifile, line);
            stringstream stream;
            stream << line;
            int count = 0;
            //int nr_exp = 0;
            char c_temp = '?';
            double temp_num[2]{0, 0};
            int dum = 0;
            if (file_type == 1)
            {
                stream >> count >> c_temp;
                if (debug)
                   std::cout << "count: " << count << " type: " << c_temp << endl;
            }
            else if (file_type == 2)
            {
                stream >> c_temp >> count;
                if (debug)
                   std::cout << "count: " << count << " type: " << c_temp << endl;
            }
            else if (file_type == 3)
            {
                stream >> c_temp >> count;
                if (debug)
                   std::cout << "count: " << count << " type: " << c_temp << endl;
            }
            for (int j = 0; j < count; j++)
            {
                getline(ifile, line);
                if (debug)
                {
                   std::cout << "read the " << j << ". line: " << line << endl;
                }
                stringstream stream2;
                stream2 << line;
                if (file_type == 1)
                    stream2 >> temp_num[0] >> temp_num[1];
                else if (file_type == 2)
                    stream2 >> dum >> temp_num[0] >> temp_num[1];
                else if (file_type == 3)
                    stream2 >> temp_num[0] >> temp_num[1];
                // this is where i started copying
                for (int h = 0; h < wave.get_ncen(); h++)
                {
                    // skip atoms taht already have a basis set!
                    if (wave.get_atom_basis_set_loaded(h))
                        continue;
                    string temp_label;
                    temp_label = wave.get_atom_label(h);
                    if (temp_label.find(" "))
                        temp_label.erase(temp_label.find(" "), 1);
                    temp_label.append(":");
                    if (elements_list[i].find(temp_label) != -1)
                    {
                        if (debug)
                        {
                           std::cout << "It's a match!" << endl;
                           std::cout << "element_label: " << elements_list[i] << " temp_label: " << temp_label << endl;
                        }
                        switch (c_temp)
                        {
                        case 's':
                        case 'S':
                            if (!wave.push_back_atom_basis_set(h, temp_num[0], temp_num[1], 1, shell))
                            {
                               std::cout << "ERROR while pushing back atoms basis set" << endl;
                            }
                            if (debug)
                               std::cout << "Pushing back on atom: " << h + 1 << " with coef: " << temp_num[1]
                                     << " and exp: " << temp_num[0] << " and type S" << endl;
                            break;
                        case 'p':
                        case 'P':
                            if (!wave.push_back_atom_basis_set(h, temp_num[0], temp_num[1], 2, shell))
                            {
                               std::cout << "ERROR while pushing back atoms basis set" << endl;
                            }
                            if (debug)
                               std::cout << "Pushing back on atom: " << h + 1 << " with coef: " << temp_num[1]
                                     << " and exp: " << temp_num[0] << " and type P" << endl;
                            break;
                        case 'd':
                        case 'D':
                            if (!wave.push_back_atom_basis_set(h, temp_num[0], temp_num[1], 3, shell))
                            {
                               std::cout << "ERROR while pushing back atoms basis set" << endl;
                            }
                            if (debug)
                               std::cout << "Pushing back on atom: " << h + 1 << " with coef: " << temp_num[1]
                                     << " and exp: " << temp_num[0] << " and type D" << endl;
                            break;
                        case 'f':
                        case 'F':
                            if (!wave.push_back_atom_basis_set(h, temp_num[0], temp_num[1], 4, shell))
                            {
                               std::cout << "ERROR while pushing back atoms basis set" << endl;
                            }
                            if (debug)
                               std::cout << "Pushing back on atom: " << h + 1 << " with coef: " << temp_num[1]
                                     << " and exp: " << temp_num[0] << " and type F" << endl;
                            break;
                        default:
                           std::cout << "Sorry, orbital types higher than f-type are not yet supported!" << endl;
                            return false;
                        } // end switch of types
                    }     // end if(find atom_label + :
                }         // end for h = ncen
                //nr_exp++;
                if (debug)
                   std::cout << "recapitulation[" << j << "]... type: " << c_temp << " coef: " << temp_num[0] << " exp: " << temp_num[1] << endl;
                if (dum > count)
                {
                   std::cout << "this should not happen, lets stop before i do something silly!" << endl;
                    return false;
                }
            }
            shell++;
        } // end while line != }
        if (debug)
           std::cout << "I found }!" << endl;
        ifile.seekg(0);
    } // end for element_list.size()
    if (debug)
    {
       std::cout << "FINISHED WITH READING BASIS SET!" << endl;
    }
    ifile.close();
    return true;
};
