#include "pch.h"
#include "JKFit.h"
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
    return { _primitives + _elementOffsets[element], static_cast<std::size_t>(_elementCounts[element])};
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


//https://pyscf.org/_modules/pyscf/df/autoaux.html
//The AutoAux algorithm by ORCA
//JCTC, 13 (2016), 554
namespace auto_aux_constants {
    inline constexpr double BETA_SMALL = 1.8;
    // Index by total L (0..), fallback used if L exceeds size
    inline const std::vector<double> BETA_BIG = { 1.8, 2.0, 2.2, 2.2, 2.2, 2.3, 3.0, 3.0, 3.0, 3.0 };
    // Cap factors for compact L (0..2*l_val)
    inline const std::vector<double> F_LAUX = {20, 7.0, 4.0, 4.0, 3.5, 2.5, 2.0, 2.0, 2.0, 2.0 };

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

    atom tmp_atom(std::to_string('H'), std::to_string('H'), 1, 0.0, 0.0, 0.0, 1);
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
        coefs = Int_Params::normalize_gto(coefs, exps, shelltype);

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
        a_max_adjusted[l] = a_aux_by_l_aux[l];
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



int load_basis_into_WFN(WFN& wavy, std::shared_ptr<BasisSet> b)
{
    wavy.set_basis_set_ptr((*b).get_data());
    int nr_coefs = 0;
    for (int i = 0; i < wavy.get_ncen(); i++)
    {
        int current_charge = wavy.get_atom_charge(i) - 1;
        const std::span<const SimplePrimitive> basis = (*b)[current_charge];
        int size = (int)basis.size();
        for (int e = 0; e < size; e++)
        {
            wavy.push_back_atom_basis_set(i, basis[e].exp, 1.0, basis[e].type, e);  //We decontract the basis, this makes it easier to handle in the long run.
            nr_coefs += 2 * basis[e].type + 1; //HIER WEITER MACHEN!!
        }

        //Different loop to keep the original contraction coefficients
        //int shell = -1;
        //for (int e = 0; e < size; e++)
        //{
        //    if (basis[e].shell != shell)
        //    {
        //        shell = basis[e].shell;
        //        nr_coefs += 2 * basis[e].type + 1;
        //    }
        //    wavy.push_back_atom_basis_set(i, basis[e].exp, basis[e].coefficient, basis[e].type, basis[e].shell); 
        //}
    }
    return nr_coefs;
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
    for (int basis_nr = 0; basis_nr < aux_basis.size(); basis_nr++) { if ((*aux_basis[basis_nr]).get_primitive_count() == 0)(*aux_basis[basis_nr]).gen_auto_aux(orbital_wfn);}

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
}//////////////////////