#include "pch.h"
#include "cell.h"
#include "convenience.h"
#include "nos_math.h"
#include <bit>
#include <cctype>
#include <Eigen/Dense>

namespace {
	/**
	 * @brief Evaluates the numerical factor of a single term of a symmetry operation.
	 *
	 * Understands an empty string (the implicit 1 of "x"), a plain number ("2", "0.5")
	 * and a fraction ("1/2").
	 *
	 * @param number the factor as written in the CIF, without its sign
	 * @param value [out] the evaluated factor
	 * @return whether the string could be evaluated
	 */
	bool eval_symop_factor(const std::string& number, double& value) {
		if (number.empty()) {
			value = 1.0;
			return true;
		}
		auto to_double = [](const std::string& s, double& out) {
			if (s.empty())
				return false;
			try {
				size_t used = 0;
				out = std::stod(s, &used);
				return used == s.length();
			}
			catch (const std::exception&) {
				return false;
			}
			};
		const size_t slash = number.find('/');
		if (slash == std::string::npos)
			return to_double(number, value);
		double numerator = 0.0, denominator = 0.0;
		if (!to_double(number.substr(0, slash), numerator))
			return false;
		if (!to_double(number.substr(slash + 1), denominator))
			return false;
		if (denominator == 0.0)
			return false;
		value = numerator / denominator;
		return true;
	}
}

void cell::parse_symop(const std::string& operation,
	const std::filesystem::path& filename,
	int rot[3][3],
	double translation[3],
	std::ostream& file) {
	const std::string where = " of symmetry operation \"" + operation + "\" in " + filename.string() + "!";
	// Split into the three comma separated components, dropping any whitespace
	svec components(3);
	int column = 0;
	for (const char c : operation) {
		if (c == ',') {
			column++;
			err_checkf(column < 3, "Found more than 3 comma separated components" + where, file);
		}
		else if (!std::isspace(static_cast<unsigned char>(c)))
			components[column].push_back(c);
	}
	err_checkf(column == 2, "Expected 3 comma separated components" + where, file);

	for (int comp = 0; comp < 3; comp++) {
		const std::string& s = components[comp];
		err_checkf(!s.empty(), "Component " + std::to_string(comp + 1) + " is empty" + where, file);
		rot[comp][0] = rot[comp][1] = rot[comp][2] = 0;
		translation[comp] = 0.0;
		// Walk the signed terms, each of which is either an axis (optionally scaled) or a translation
		size_t pos = 0;
		while (pos < s.length()) {
			double sign = 1.0;
			if (s[pos] == '+')
				pos++;
			else if (s[pos] == '-')
				sign = -1.0, pos++;
			size_t end = pos;
			while (end < s.length() && s[end] != '+' && s[end] != '-')
				end++;
			err_checkf(end > pos, "Found an empty term" + where, file);
			std::string term = s.substr(pos, end - pos);
			pos = end;
			// Pull out the axis name, if this term has one; what remains is its factor
			int axis = -1;
			for (size_t k = 0; k < term.length() && axis == -1; k++) {
				switch (term[k]) {
				case 'x': case 'X': axis = 0; break;
				case 'y': case 'Y': axis = 1; break;
				case 'z': case 'Z': axis = 2; break;
				default: continue;
				}
				term.erase(k, 1);
			}
			// "2*x" carries the same information as "2x"
			term.erase(std::remove(term.begin(), term.end(), '*'), term.end());
			double factor = 0.0;
			err_checkf(eval_symop_factor(term, factor), "Could not interpret the factor \"" + term + "\"" + where, file);
			if (axis == -1) {
				translation[comp] += sign * factor;
				continue;
			}
			const double coefficient = sign * factor;
			const int rounded = static_cast<int>(std::lround(coefficient));
			err_checkf(std::abs(coefficient - rounded) < 1e-6,
				"The rotation coefficient " + std::to_string(coefficient) + " is not an integer" + where, file);
			rot[comp][axis] += rounded;
		}
	}
}

vec cell::apply_symmetry(const vec& pos, const int sym_op) {
	const vec trans_temp = { trans[0][sym_op], trans[1][sym_op], trans[2][sym_op] };
	const vec2 rot_temp = { { (double)sym[0][0][sym_op], (double)sym[0][1][sym_op], (double)sym[0][2][sym_op] },
							 { (double)sym[1][0][sym_op], (double)sym[1][1][sym_op], (double)sym[1][2][sym_op] },
							 { (double)sym[2][0][sym_op], (double)sym[2][1][sym_op], (double)sym[2][2][sym_op] } };
	vec temp_pos = self_dot(rot_temp, pos, true);
	return { temp_pos[0] + trans_temp[0], temp_pos[1] + trans_temp[1], temp_pos[2] + trans_temp[2] };
	// closing function
}

std::string cell::symop_string(const int sym_op) const {
	std::string s;
	for (int i = 0; i < 3; i++) {
		std::string c;
		//sym holds R transposed, see apply_symmetry
		for (int j = 0; j < 3; j++) {
			const int r = sym[j][i][sym_op];
			if (r == 0) continue;
			c += r < 0 ? "-" : "+";
			if (std::abs(r) != 1) c += std::to_string(std::abs(r));
			c += "xyz"[j];
		}
		const double t = trans[i][sym_op];
		if (std::abs(t) > 1e-8) {
			int d = 1;
			while (d < 12 && std::abs(t * d - std::round(t * d)) > 1e-6) d++;
			const long n = std::lround(t * d);
			c += (n < 0 ? "-" : "+") + std::to_string(std::abs(n)) + (d > 1 ? "/" + std::to_string(d) : "");
		}
		if (c.empty()) c = "0";
		if (c[0] == '+') c.erase(0, 1);
		s += (i ? "," : "") + c;
	}
	return s;
}

void cell::convert_to_fracs(std::vector<asym_atom>& atoms, const std::string input_unit) {
	for (asym_atom& atom : atoms) {
		occ::Vec temp_cart_pos(3);
		temp_cart_pos << atom.pos[0], atom.pos[1], atom.pos[2];
		occ::Mat3N inv_cell(3, 3);
		if (input_unit == "bohr") {
			inv_cell << cm[0][0], cm[1][0], cm[2][0], cm[0][1], cm[1][1], cm[2][1], cm[0][2], cm[1][2], cm[2][2];
		}
		else if (input_unit == "angstrom") {
			inv_cell << constants::bohr2ang(cm[0][0]), constants::bohr2ang(cm[1][0]), constants::bohr2ang(cm[2][0]), constants::bohr2ang(cm[0][1]), constants::bohr2ang(cm[1][1]), constants::bohr2ang(cm[2][1]), constants::bohr2ang(cm[0][2]), constants::bohr2ang(cm[1][2]), constants::bohr2ang(cm[2][2]);
		}
		else {
			std::cerr << "Unknown input unit. Choose 'bohr' or 'angstrom'.\n";
		}
		inv_cell = inv_cell.inverse();
		occ::Vec frac_vec = inv_cell * temp_cart_pos;
		atom.frac_pos = { frac_vec(0), frac_vec(1), frac_vec(2) };
	}
}

void cell::grow_asym_atoms(std::vector<asym_atom>& asym_atoms, std::vector<asym_atom>& xyz_atoms) {
	// Add fractional coordinates to the atoms from the xyz file since these are way easier to compare
	convert_to_fracs(xyz_atoms, "bohr");
	/* Find out which atoms are already present in the asym_atoms vector and which are not
	This only handles the growing of the asym_atoms vector, symmetry equivalency is checked later on to link the asymetric unit atoms with the grown ones
	Important: Label is a dummy and should not be used before being properly set (hopefully I will remember to set it later on ^^)*/
	for (asym_atom& xyz_atom : xyz_atoms) {
		bool found = false;
		const vec xyz_pos = { xyz_atom.frac_pos[0], xyz_atom.frac_pos[1], xyz_atom.frac_pos[2] };
		for (const asym_atom& asym_atom : asym_atoms) {
			const vec asym_unit_pos = { asym_atom.frac_pos[0], asym_atom.frac_pos[1], asym_atom.frac_pos[2] };
			if (check_special(xyz_pos, asym_unit_pos, 1e-2)) {
				found = true;
				break;
			}
		}
		if (!found) {
			xyz_atom.grown = true;
			asym_atoms.push_back(xyz_atom);
		}
	}
}

void cell::resolve_internal_symmetry(std::vector<asym_atom>& asym_atoms) {
	// Loop over the grown atoms
	for (asym_atom& grown_atom : asym_atoms) {
		if (!grown_atom.grown) {
			continue;
		}
		const vec grown_pos = { grown_atom.frac_pos[0], grown_atom.frac_pos[1], grown_atom.frac_pos[2] };
		// Now we have a grown atom, let's find out which asymmetric unit atom it was generated from and by which symmetry
		bool found = false;
		for (int sym_op = 0; sym_op < sym[0][0].size(); sym_op++) {
			if (check_identity(sym_op)) {
				continue;
			}
			for (asym_atom& asymm_atom : asym_atoms) {
				if (asymm_atom.grown) {
					continue;
				}
				const vec asymm_pos = { asymm_atom.frac_pos[0], asymm_atom.frac_pos[1], asymm_atom.frac_pos[2] };
				const vec symmetry_generated_pos = apply_symmetry(asymm_pos, sym_op);
				if (check_special(grown_pos, symmetry_generated_pos, 1e-2)) {
					grown_atom.grown_from = &asymm_atom - &asym_atoms[0];
					if (grown_atom.grown_by[0] == -1) {
						grown_atom.grown_by[0] = sym_op;
					}
					else {
						grown_atom.grown_by.push_back(sym_op);
					}
					grown_atom.cryst_mult = asymm_atom.cryst_mult;
				}
				if (check_special(asymm_pos, apply_symmetry(grown_pos, sym_op), 1e-2)) {
					grown_atom.maps_to_parent.push_back(sym_op);
				}
			}
		}
	}
}

void cell::eval_symm(std::vector<asym_atom>& asym_atoms, const int& asymmetric_atoms, ivec3& linking_list) {
	const int total_atoms = asym_atoms.size();
	auto frac_pos = [&](int i) -> vec {
		return { asym_atoms[i].frac_pos[0], asym_atoms[i].frac_pos[1], asym_atoms[i].frac_pos[2] };
		};
	linking_list.resize(asymmetric_atoms);
	const int num_sym_ops = sym[0][0].size();
	int idx1 = 0;
	for (int idx1 = 0; idx1 < asymmetric_atoms; idx1++) {
		linking_list[idx1].resize(total_atoms);
		for (int sym_op = 0; sym_op < num_sym_ops; sym_op++) {
			const vec pos1 = apply_symmetry(frac_pos(idx1), sym_op);
			for (int idx2 = 0; idx2 < total_atoms; idx2++) {
				const vec pos2 = frac_pos(idx2);
				if (check_special(pos1, pos2, 1e-4)) {
					linking_list[idx1][idx2].push_back(sym_op);
					break;
				}
			}
		}
	}
	// Closing function
}

bool cell::check_special(const vec& pos1, const vec& pos2, const double& tolerance) {
	for (int i = 0; i < 3; ++i)
	{
		double diff = pos2[i] - pos1[i];
		if (std::abs(diff - std::round(diff)) > tolerance)
			return false;
	}
	return true;
	//closing function
}

bool cell::check_identity(const int& sym_op) {
	bool is_identity = false;
	vec trans_identity = { 0 ,0, 0 };
	vec2 rot_identity = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
	vec actual_trans = { trans[0][sym_op], trans[1][sym_op], trans[2][sym_op] };
	vec2 actual_rot = { { (double)sym[0][0][sym_op], (double)sym[0][1][sym_op], (double)sym[0][2][sym_op] },
						{ (double)sym[1][0][sym_op], (double)sym[1][1][sym_op], (double)sym[1][2][sym_op] },
						{ (double)sym[2][0][sym_op], (double)sym[2][1][sym_op], (double)sym[2][2][sym_op] } };
	if (trans_identity == actual_trans && rot_identity == actual_rot) {
		is_identity = true;
	}
	return is_identity;
}

// This function does not check for improperly applied symmetry operations anymore, hoping that the new symmetry resolution can handle such cases
ivec cell::confirm_applied_symmetry(const std::vector<asym_atom>& asym_atoms) {
	const int asymmetric_atoms = asym_atoms.size();
	const int num_sym_ops = sym[0][0].size();
	ivec applied_symmetry;

	// Count how often a symmetry operation is applied (check if it is applied to all asymmetric atoms)
	for (int sym_op = 0; sym_op < num_sym_ops; sym_op++) {
		if (check_identity(sym_op)) {
			continue;
		}
		int count = 0;
		// This loop counts the asymmetric atoms whose image under the symmetry operation is in the model, either the atom itself
		// (a special position) or a grown atom, whichever operation that grown atom was credited to in resolve_internal_symmetry
		for (const asym_atom& asymm_atom : asym_atoms) {
			if (asymm_atom.grown) {
				continue;
			}
			const vec original_pos = { asymm_atom.frac_pos[0], asymm_atom.frac_pos[1], asymm_atom.frac_pos[2] };
			const vec symmetry_generated_pos = apply_symmetry(original_pos, sym_op);
			bool present = check_special(original_pos, symmetry_generated_pos, 1e-4);
			for (const asym_atom& grown_atom : asym_atoms) {
				if (present) {
					break;
				}
				if (grown_atom.grown) {
					present = check_special({ grown_atom.frac_pos[0], grown_atom.frac_pos[1], grown_atom.frac_pos[2] }, symmetry_generated_pos, 1e-2);
				}
			}
			if (present) {
				count++;
			}
		}
		// Now we have counted the number of atoms in the model, for which the symmetry operation was applied
		// If this number is equal to the number of asymmetric atoms, then the symmetry operation is fully applied
		int original_ncen = 0;
		for (const asym_atom& asym_atom : asym_atoms) {
			if (!asym_atom.grown) {
				original_ncen++;
			}
		}
		if (count == original_ncen) {
			applied_symmetry.push_back(sym_op);
		}
	}
	return applied_symmetry;
}

void cell::asym_unit_cryst_mult(std::vector<asym_atom>& asym_atoms) {
	for (asym_atom& a : asym_atoms) {
		int counter = 1;
		vec originial_pos = { a.frac_pos[0], a.frac_pos[1], a.frac_pos[2] };
		for (int sym_op = 0; sym_op < sym[0][0].size(); sym_op++) {
			if (check_identity(sym_op)) {
				continue;
			}
			vec symmetry_pos = apply_symmetry(originial_pos, sym_op);
			if (check_special(originial_pos, symmetry_pos, 1e-4)) {
				counter++;
			}
		}
		a.cryst_mult = counter;
	}
}

namespace {
	// The subgroup the generators generate: every operation reached from the identity by applying them again and again
	ivec generated_subgroup(const ivec2& product, const int identity, const ivec& generators) {
		bvec member(product.size(), false);
		member[identity] = true;
		ivec subgroup = { identity };
		for (int next = 0; next < subgroup.size(); next++)
			for (int generator = 0; generator < generators.size(); generator++) {
				const int op = product[subgroup[next]][generators[generator]];
				if (!member[op]) {
					member[op] = true;
					subgroup.push_back(op);
				}
			}
		std::sort(subgroup.begin(), subgroup.end());
		return subgroup;
	}

	// Every subgroup: the cyclic ones, joined with each other until no join gives a new subgroup
	ivec2 all_subgroups(const ivec2& product, const int identity) {
		std::map<ivec, ivec> generators_of;
		for (int op = 0; op < product.size(); op++)
			generators_of.emplace(generated_subgroup(product, identity, { op }), ivec{ op });
		const std::vector<std::pair<ivec, ivec>> cyclic(generators_of.begin(), generators_of.end());
		std::vector<std::pair<ivec, ivec>> frontier = cyclic;
		while (!frontier.empty()) {
			std::vector<std::pair<ivec, ivec>> next_frontier;
			for (int known = 0; known < frontier.size(); known++)
				for (int added = 0; added < cyclic.size(); added++) {
					const ivec& members = frontier[known].first;
					if (std::includes(members.begin(), members.end(), cyclic[added].first.begin(), cyclic[added].first.end()))
						continue;
					ivec generators = frontier[known].second;
					generators.insert(generators.end(), cyclic[added].second.begin(), cyclic[added].second.end());
					const ivec joined = generated_subgroup(product, identity, generators);
					if (generators_of.emplace(joined, generators).second)
						next_frontier.push_back({ joined, generators });
				}
			frontier = next_frontier;
		}
		ivec2 subgroups;
		for (const auto& [members, generators] : generators_of)
			subgroups.push_back(members);
		return subgroups;
	}

	// Each family needs a model atom that some operation of every right coset L.h of H puts on a model atom
	bool families_fit(const ivec2& product, const ivec& coset_subgroup, const ivec& kept_subgroup, const std::vector<std::vector<bvec>>& lands_in_model) {
		ivec2 right_cosets;
		bvec assigned(product.size(), false);
		for (int outer = 0; outer < coset_subgroup.size(); outer++) {
			if (assigned[coset_subgroup[outer]]) continue;
			ivec coset;
			for (int inner = 0; inner < kept_subgroup.size(); inner++) {
				coset.push_back(product[kept_subgroup[inner]][coset_subgroup[outer]]);
				assigned[coset.back()] = true;
			}
			right_cosets.push_back(coset);
		}
		for (int family = 0; family < lands_in_model.size(); family++) {
			bool fits = false;
			for (int member = 0; member < lands_in_model[family].size() && !fits; member++) {
				fits = true;
				for (int coset = 0; coset < right_cosets.size() && fits; coset++) {
					bool reaches_model = false;
					for (int op = 0; op < right_cosets[coset].size(); op++)
						reaches_model = reaches_model || lands_in_model[family][member][right_cosets[coset][op]];
					fits = reaches_model;
				}
			}
			if (!fits) return false;
		}
		return true;
	}

	// Branch and bound: every coset takes one footprint (the classes of its d.L), the union of all as small as possible
	void cover_classes(const std::vector<std::map<unsigned, int>>& options_per_coset, const unsigned covered, unsigned& best_covered, int& best_count) {
		if (std::popcount(covered) >= best_count) return;
		int open_coset = -1;
		for (int coset = 0; coset < options_per_coset.size(); coset++) {
			bool done = false;
			for (const auto& [footprint, op] : options_per_coset[coset])
				done = done || (footprint & ~covered) == 0;
			if (!done && (open_coset < 0 || options_per_coset[coset].size() < options_per_coset[open_coset].size()))
				open_coset = coset;
		}
		if (open_coset < 0) {
			best_covered = covered;
			best_count = std::popcount(covered);
			return;
		}
		std::vector<unsigned> footprints;
		for (const auto& [footprint, op] : options_per_coset[open_coset])
			footprints.push_back(footprint);
		std::sort(footprints.begin(), footprints.end(), [covered](const unsigned first, const unsigned second) {
			return std::popcount(first & ~covered) < std::popcount(second & ~covered); });
		for (int option = 0; option < footprints.size(); option++)
			cover_classes(options_per_coset, covered | footprints[option], best_covered, best_count);
	}

	// K = D.L with one d per left coset d.H, no d.L holding a forbidden operation, as few classes as possible; empty if none exists
	ivec cheapest_kept_set(const ivec2& product, const ivec& rotation_class, const int identity, const ivec& coset_subgroup, const ivec& kept_subgroup, const bvec& forbidden) {
		std::vector<std::map<unsigned, int>> options_per_coset;
		bvec assigned(product.size(), false);
		for (int start = 0; start < product.size(); start++) {
			if (assigned[start]) continue;
			ivec coset;
			for (int outer = 0; outer < coset_subgroup.size(); outer++) {
				coset.push_back(product[start][coset_subgroup[outer]]);
				assigned[coset.back()] = true;
			}
			const bool holds_identity = std::find(coset.begin(), coset.end(), identity) != coset.end();
			std::map<unsigned, int> options;
			for (int candidate = 0; candidate < coset.size(); candidate++) {
				if (holds_identity && coset[candidate] != identity) continue;
				unsigned footprint = 0;
				bool allowed = true;
				for (int inner = 0; inner < kept_subgroup.size(); inner++) {
					const int op = product[coset[candidate]][kept_subgroup[inner]];
					allowed = allowed && !forbidden[op];
					footprint |= 1u << rotation_class[op];
				}
				if (allowed) options.emplace(footprint, coset[candidate]);
			}
			if (options.empty()) return {};
			options_per_coset.push_back(options);
		}
		unsigned best_covered = 0;
		int best_count = 33;
		cover_classes(options_per_coset, 0u, best_covered, best_count);
		ivec kept;
		for (int coset = 0; coset < options_per_coset.size(); coset++)
			for (const auto& [footprint, representative] : options_per_coset[coset])
				if ((footprint & ~best_covered) == 0) {
					for (int inner = 0; inner < kept_subgroup.size(); inner++)
						kept.push_back(product[representative][kept_subgroup[inner]]);
					break;
				}
		std::sort(kept.begin(), kept.end());
		kept.erase(std::unique(kept.begin(), kept.end()), kept.end());
		return kept;
	}
}

// K = D.L for subgroups L inside H, D one operation per left coset of H: a model atom per right coset L.h of H gives a set C with K.C = D.H = G
ivec cell::choose_kept_operations(const std::vector<asym_atom>& asym_atoms, const ivec& applied_symmetry) {
	const int num_sym_ops = static_cast<int>(sym[0][0].size());
	int identity = 0;
	while (!check_identity(identity)) identity++;
	ivec2 product(num_sym_ops, ivec(num_sym_ops));
	for (int first = 0; first < num_sym_ops; first++) {
		for (int second = 0; second < num_sym_ops; second++) {
			product[first][second] = equal_to_concatenation(first, second);
		}
	}
	// Operations with equal or opposite rotations share one k-point after inversion screening: one class, the cost is the number of classes
	ivec rotation_class(num_sym_ops, -1), rotation_sign(num_sym_ops, 1);
	int num_classes = 0;
	for (int op = 0; op < num_sym_ops; op++) {
		for (int earlier = 0; earlier < op && rotation_class[op] < 0; earlier++) {
			for (int sign = -1; sign <= 1 && rotation_class[op] < 0; sign += 2) {
				bool same = true;
				for (int i = 0; i < 3 && same; i++) {
					for (int j = 0; j < 3 && same; j++) {
						same = sym[i][j][op] == sign * sym[i][j][earlier];
					}
				}
				if (same) {
					rotation_class[op] = rotation_class[earlier];
					rotation_sign[op] = sign * rotation_sign[earlier];
				}
			}
		}
		if (rotation_class[op] < 0) rotation_class[op] = num_classes++;
	}
	err_checkf(num_classes <= 32, "More than 32 rotation classes in the symmetry operations", std::cout);
	auto class_count = [&](const ivec& ops) {
		unsigned classes = 0;
		for (int op = 0; op < ops.size(); op++) classes |= 1u << rotation_class[ops[op]];
		return std::popcount(classes);
		};
	// Kept operations whose opposite rotation is kept too: the pairs inversion screening integrates once
	auto inversion_partners = [&](const ivec& ops) {
		int partners = 0;
		for (int op = 0; op < ops.size(); op++) {
			bool has_partner = false;
			for (int other = 0; other < ops.size() && !has_partner; other++)
				has_partner = rotation_class[ops[other]] == rotation_class[ops[op]] && rotation_sign[ops[other]] != rotation_sign[ops[op]];
			partners += has_partner;
		}
		return partners;
		};
	// Per family (an asymmetric atom and the atoms grown from it), which operations put each member on a member
	bvec3 lands_in_model;
	for (int parent = 0; parent < asym_atoms.size(); parent++) {
		if (asym_atoms[parent].grown) continue;
		ivec family = { parent };
		for (int atom = 0; atom < asym_atoms.size(); atom++)
			if (asym_atoms[atom].grown && asym_atoms[atom].grown_from == parent) family.push_back(atom);
		bvec2 member_lands(family.size(), bvec(num_sym_ops, false));
		for (int member = 0; member < family.size(); member++)
			for (int op = 0; op < num_sym_ops; op++) {
				const vec image = apply_symmetry({ asym_atoms[family[member]].frac_pos[0], asym_atoms[family[member]].frac_pos[1], asym_atoms[family[member]].frac_pos[2] }, op);
				for (int target = 0; target < family.size() && !member_lands[member][op]; target++)
					member_lands[member][op] = check_special({ asym_atoms[family[target]].frac_pos[0], asym_atoms[family[target]].frac_pos[1], asym_atoms[family[target]].frac_pos[2] }, image, 1e-2);
			}
		lands_in_model.push_back(member_lands);
	}
	// Pairs L inside H with the fewest classes any D.L can have: classes(G) * classes(L) / classes(H)
	const ivec2 subgroups = all_subgroups(product, identity);
	std::vector<std::tuple<double, int, int>> pairs;
	for (int coset_subgroup = 0; coset_subgroup < subgroups.size(); coset_subgroup++)
		for (int kept_subgroup = 0; kept_subgroup < subgroups.size(); kept_subgroup++)
			if (std::includes(subgroups[coset_subgroup].begin(), subgroups[coset_subgroup].end(), subgroups[kept_subgroup].begin(), subgroups[kept_subgroup].end()))
				pairs.push_back({ static_cast<double>(num_classes) * class_count(subgroups[kept_subgroup]) / class_count(subgroups[coset_subgroup]), coset_subgroup, kept_subgroup });
	std::sort(pairs.begin(), pairs.end());
	// Fewest classes, then most inversion partners, then fewest operations; first every applied operation stays out, then only those whose powers are all applied, then none
	bvec forbidden(num_sym_ops, false);
	for (int op = 0; op < applied_symmetry.size(); op++) forbidden[applied_symmetry[op]] = true;
	for (int attempt = 0; attempt < 3; attempt++) {
		if (attempt == 1) {
			const bvec applied = forbidden;
			for (int op = 0; op < num_sym_ops; op++)
				for (int power = op; power != identity && forbidden[op]; power = product[power][op])
					forbidden[op] = applied[power];
		}
		if (attempt == 2) forbidden.assign(num_sym_ops, false);
		ivec best_kept;
		int best_cost = num_classes + 1, best_partners = -1;
		for (int pair = 0; pair < pairs.size(); pair++) {
			const auto& [lowest_cost, coset_subgroup, kept_subgroup] = pairs[pair];
			if (lowest_cost > best_cost + 1e-9) break;
			bool deletes_none = true;
			for (int op = 0; op < subgroups[kept_subgroup].size(); op++)
				deletes_none = deletes_none && !forbidden[subgroups[kept_subgroup][op]];
			if (!deletes_none || !families_fit(product, subgroups[coset_subgroup], subgroups[kept_subgroup], lands_in_model)) continue;
			const ivec kept = cheapest_kept_set(product, rotation_class, identity, subgroups[coset_subgroup], subgroups[kept_subgroup], forbidden);
			if (kept.empty()) continue;
			const int cost = class_count(kept), partners = inversion_partners(kept);
			if (cost < best_cost || (cost == best_cost && (partners > best_partners || (partners == best_partners && kept.size() < best_kept.size())))) {
				best_kept = kept;
				best_cost = cost;
				best_partners = partners;
			}
		}
		if (best_kept.empty()) continue;
		if (attempt > 0)
			std::cout << "No kept symmetry operations count every position once with all applied operations deleted, some of them stay in the sum\n";
		return best_kept;
	}
	err_checkf(false, "No kept symmetry operations count every position once", std::cout);
	return {};
}

// Deletes every operation not in kept and enlarges hkl by the kept rotations; applied_symmetry returns the deleted operations
void cell::keep_operations(const ivec& kept, ivec& applied_symmetry, hkl_list& hkl_enlarged, const hkl_list& hkl) {
	const int num_sym_ops = static_cast<int>(sym[0][0].size());
	applied_symmetry.clear();
	for (int sym_op = 0; sym_op < num_sym_ops; sym_op++)
		if (std::find(kept.begin(), kept.end(), sym_op) == kept.end()) {
			applied_symmetry.push_back(sym_op);
		}
	for (const int sym_op : std::views::reverse(applied_symmetry)) {
		for (ivec2& middle : sym) {
			for (ivec& inner : middle) {
				inner.erase(inner.begin() + sym_op);
			}
		}
		for (vec& inner : trans) {
			inner.erase(inner.begin() + sym_op);
		}
	}

	vec3 rotations;
	for (int sym_op = 0; sym_op < sym[0][0].size(); ++sym_op) {
		rotations.push_back({ {
				static_cast<double>(sym[0][0][sym_op]),
				static_cast<double>(sym[0][1][sym_op]),
				static_cast<double>(sym[0][2][sym_op])},{
				static_cast<double>(sym[1][0][sym_op]),
				static_cast<double>(sym[1][1][sym_op]),
				static_cast<double>(sym[1][2][sym_op])},{
				static_cast<double>(sym[2][0][sym_op]),
				static_cast<double>(sym[2][1][sym_op]),
				static_cast<double>(sym[2][2][sym_op])} });
	}

	const int nr = hkl.size();
	std::vector<i3> hkl_vec(hkl.begin(), hkl.end());
	hkl_list new_enlarged;
	for (const auto& h : hkl) {
		vec hkl_temp = { static_cast<double>(h[0]), static_cast<double>(h[1]), static_cast<double>(h[2]) };
		for (const vec2& rot : rotations) {
			vec new_hkl = self_dot(rot, hkl_temp, false);
			i3 new_hkl_int = { (int)std::round(new_hkl[0]), (int)std::round(new_hkl[1]), (int)std::round(new_hkl[2]) };
			new_enlarged.insert(new_hkl_int);
		}
	}
	hkl_enlarged = new_enlarged;
	//closing function
}

void cell::grow_ADPs(std::vector<asym_atom>& asym_atoms, vec3& ADPs) {
	vec3 new_ADPs(asym_atoms.size());
	for (int i = 0; i < ADPs.size(); i++) {
		new_ADPs[i] = ADPs[i];
	}
	// The grown atoms follow the asymmetric ones; each takes its parent's values, the ADPs rotated by the linking operation
	int idx = 0;
	for (asym_atom& grown_atom : asym_atoms) {
		if (!grown_atom.grown) {
			idx++;
			continue;
		}
		const int parent_idx = grown_atom.grown_from;
		const int sym_op = grown_atom.grown_by[0];
		new_ADPs[idx] = ADPs[parent_idx];
		rotate_grown_ADPs(new_ADPs[idx], sym_op);
		grown_atom.U_iso = asym_atoms[parent_idx].U_iso;
		grown_atom.anom = asym_atoms[parent_idx].anom;
		idx++;
	}
	ADPs = new_ADPs;
};

// U*, C and D are contravariant tensors in the fractional basis, so an image's are T' = R T R^T with R
// the rotation of sym_op (x' = R x + t); sym holds R transposed, which is what transform_ADPs takes.
// U comes as the CIF gives it, U*_ij / (a*_i a*_j), so it is rotated as U* and scaled back.
void cell::rotate_grown_ADPs(vec2& ADPs, const int sym_op) const {
	vec2 M(3, vec(3));
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			M[i][j] = sym[i][j][sym_op];
	transform_ADPs(ADPs, M);
}

// Position of a sorted index triple/quadruple in the Voigt storage of C (10) and D (15)
static int get_voigt_index(const ivec& indices) {
	static const ivec2 map3{ { 0, 0, 0 }, { 0, 0, 1 }, { 0, 0, 2 }, { 0, 1, 1 }, {0, 1, 2}, {0, 2, 2}, {1, 1, 1}, { 1, 1, 2 }, { 1, 2, 2 }, { 2, 2, 2 } };
	static const ivec2 map4{ { 0, 0, 0, 0 }, { 0, 0, 0, 1 }, { 0, 0, 0, 2 }, { 0, 0, 1, 1 }, { 0, 0, 1, 2 }, { 0, 0, 2, 2 }, { 0, 1, 1, 1 }, { 0, 1, 1, 2 }, { 0, 1, 2, 2 }, { 0, 2, 2, 2 }, { 1, 1, 1, 1 }, { 1, 1, 1, 2 }, { 1, 1, 2, 2 }, { 1, 2, 2, 2 }, { 2, 2, 2, 2 } };
	const auto& map = indices.size() == 3 ? map3 : map4;
	return std::find(map.begin(), map.end(), indices) - map.begin();
}

//Maybe an alternative to transform_ADPs but I don't like it yet since it has many auxiliary functions and is really hard to read.
//template <std::size_t N>
//int flat_index(const std::array<int, N>& t) {
//	int f = 0;
//	for (int i : t) f = 3 * f + i;
//	return f;
//}
//
//template <std::size_t N, std::size_t K>
//vec transform_symmetric(const vec& packed, const vec2& M, const std::array<std::array<int, N>, K>& order) {
//	assert(packed.size() == K);
//	constexpr int full = ipow3(N);
//
//	// Unpack into the full tensor: every permutation of a stored tuple gets its value.
//	std::array<double, full> T{};
//	for (std::size_t c = 0; c < K; ++c) {
//		auto t = order[c];  // tuples are stored sorted, so next_permutation hits all of them
//		do T[flat_index(t)] = packed[c];
//		while (std::next_permutation(t.begin(), t.end()));
//	}
//
//	// Contract M into one axis at a time: N * 3^(N+1) flops instead of 3^(2N).
//	for (std::size_t axis = 0; axis < N; ++axis) {
//		const int stride = ipow3(N - 1 - axis);
//		std::array<double, full> R{};
//		for (int f = 0; f < full; ++f) {
//			const int i = (f / stride) % 3;
//			const int base = f - i * stride;
//			for (int p = 0; p < 3; ++p)
//				R[f] += M[p][i] * T[base + p * stride];
//		}
//		T = R;
//	}
//
//	// Pack back into the independent components.
//	vec out(K);
//	for (std::size_t c = 0; c < K; ++c)
//		out[c] = T[flat_index(order[c])];
//	return out;
//}
//
//
//constexpr int ipow3(std::size_t n) { return n == 0 ? 1 : 3 * ipow3(n - 1); }
//void cell::transform_ADPs(vec2& ADPs, const vec2& M) {
//	constexpr std::array<std::array<int, 2>, 6> voigt2{ { {0,0}, {1,1}, {2,2}, {0,1}, {0,2}, {1,2} } };
//	constexpr std::array<std::array<int, 3>, 10> voigt3{ {
//		{0,0,0}, {0,0,1}, {0,0,2}, {0,1,1}, {0,1,2}, {0,2,2}, {1,1,1}, {1,1,2}, {1,2,2}, {2,2,2} } };
//	constexpr std::array<std::array<int, 4>, 15> voigt4{ {{0,0,0,0}, {0,0,0,1}, {0,0,0,2}, {0,0,1,1}, {0,0,1,2}, {0,0,2,2}, {0,1,1,1}, {0,1,1,2},
//		{0,1,2,2}, {0,2,2,2}, {1,1,1,1}, {1,1,1,2}, {1,1,2,2}, {1,2,2,2}, {2,2,2,2} } };
//	if (ADPs.size() > 0 && !ADPs[0].empty()) ADPs[0] = transform_symmetric(ADPs[0], M, voigt2);
//	if (ADPs.size() > 1 && !ADPs[1].empty()) ADPs[1] = transform_symmetric(ADPs[1], M, voigt3);
//	if (ADPs.size() > 2 && !ADPs[2].empty()) ADPs[2] = transform_symmetric(ADPs[2], M, voigt4);
//}



// T'_{ij..} = sum M_pi M_qj .. T_pq.. for the U (rank 2), C (rank 3) and D (rank 4) tensors in their
// Voigt storage: structure_factors::U_star2U_cart hands in the cell matrix, rotate_grown_ADPs the transposed symmetry operation
void cell::transform_ADPs(vec2& ADPs, const vec2& M) {
	if (ADPs.size() > 0 && ADPs[0].size() > 0) {
		vec2 U(3, vec(3));
		U[0][0] = ADPs[0][0];
		U[0][1] = ADPs[0][3];
		U[0][2] = ADPs[0][4];
		U[1][0] = ADPs[0][3];
		U[1][1] = ADPs[0][1];
		U[1][2] = ADPs[0][5];
		U[2][0] = ADPs[0][4];
		U[2][1] = ADPs[0][5];
		U[2][2] = ADPs[0][2];
		U = self_dot(self_dot(M, U, true, false), M, false, false);
		ADPs[0][0] = U[0][0];
		ADPs[0][1] = U[1][1];
		ADPs[0][2] = U[2][2];
		ADPs[0][3] = U[0][1];
		ADPs[0][4] = U[0][2];
		ADPs[0][5] = U[1][2];
	}
	if (ADPs.size() > 1 && ADPs[1].size() > 0) {
		int running_idx = 0;
		vec C_out(10);
		for (int i = 0, idx = 0; i < 3; i++) {
			for (int j = i; j < 3; j++) {
				for (int k = j; k < 3; k++) {
					double sum = 0;
					for (int p = 0; p < 3; p++) {
						for (int q = 0; q < 3; q++) {
							for (int r = 0; r < 3; r++) {
								ivec sorted_idx = { p, q, r };
								std::sort(sorted_idx.begin(), sorted_idx.end());
								int ADP_idx;
								ADP_idx = get_voigt_index(sorted_idx);
								sum += M[p][i] * M[q][j] * M[r][k] * ADPs[1][ADP_idx];
							}
						}
					}
					C_out[idx++] = sum;
				}
			}
		}
		ADPs[1] = C_out;
	}
	if (ADPs.size() > 2 && ADPs[2].size() > 0) {
		int running_idx = 0;
		vec D_out(15);
		for (int i = 0; i < 3; i++) {
			for (int j = i; j < 3; j++) {
				for (int k = j; k < 3; k++) {
					for (int l = k; l < 3; l++) {
						double sum = 0;
						for (int p = 0; p < 3; p++) {
							for (int q = 0; q < 3; q++) {
								for (int r = 0; r < 3; r++) {
									for (int s = 0; s < 3; s++) {
										ivec sorted_idx = { p, q, r, s };
										std::sort(sorted_idx.begin(), sorted_idx.end());
										int ADP_idx;
										ADP_idx = get_voigt_index(sorted_idx);
										sum += M[p][i] * M[q][j] * M[r][k] * M[s][l] * ADPs[2][ADP_idx];
									}
								}
							}
						}
						D_out[running_idx] = sum;
						running_idx++;
					}
				}
			}
		}
		ADPs[2] = D_out;
	}
}

bool cell::check_inversion(const i3& h1, const i3& h2) {
	return h1[0] == -h2[0] && h1[1] == -h2[1] && h1[2] == -h2[2];
}

// -h is looked up among all entries, not only among images under an inversion operation: two images of
// different operations can be each other's inverse, and a grown structure need not keep the inversion at all
ivec2 cell::remove_inv_vectors(const hkl_list& hkl_enlarged, hkl_list& inversion_cleaned_hkl) {
	const std::vector<i3> h(hkl_enlarged.begin(), hkl_enlarged.end());
	const int n = static_cast<int>(h.size());
	bvec paired(n, false);
	ivec2 link;
	inversion_cleaned_hkl.clear();
	for (int i = 0; i < n; i++) {
		if (paired[i]) continue;
		link.push_back({ i });
		inversion_cleaned_hkl.insert(inversion_cleaned_hkl.end(), h[i]);
		// The entries are sorted and unique, so a vector has one inverse at most and a binary search finds it
		const i3 minus_h = { -h[i][0], -h[i][1], -h[i][2] };
		const int j = static_cast<int>(std::lower_bound(h.begin(), h.end(), minus_h) - h.begin());
		if (j > i && j < n && check_inversion(h[i], h[j])) {
			link.back().push_back(j);
			paired[j] = true;
		}
	}
	return link;
}

ivec cell::generate_asym_lookup(const int r, const hkl_list& hkl, const hkl_list& hkl_enlarged) const {
	ivec asym_list;
	auto it = hkl.begin();
	std::advance(it, r);
	i3 tempv;
	const i3& hkl_temp = *it;
	for (int s = 0; s < sym[0][0].size(); s++) {
		tempv = { 0, 0, 0 };
		for (int h = 0; h < 3; h++) {
			for (int j = 0; j < 3; j++) {
				tempv[j] += hkl_temp[h] * sym[j][h][s];
			}
		}
		//a miss would silently read reflection 0's scattering factors
		auto idx = hkl_enlarged.find(tempv);
		err_checkf(idx != hkl_enlarged.end(), "Reflection " + std::to_string(hkl_temp[0]) + " " + std::to_string(hkl_temp[1]) + " " +
			std::to_string(hkl_temp[2]) + " rotated by symmetry operation " + std::to_string(s) + " is not in the enlarged reflection list", std::cout);
		asym_list.push_back(static_cast<int>(std::distance(hkl_enlarged.begin(), idx)));
	}
	return asym_list;
}

ivec cell::apply_grown(const hkl_list& hkl, hkl_list& hkl_enlarged, std::vector<asym_atom>& asym_atoms) {
	// Check which symmetry operations are fully applied to the asymmetric unit, ONLY those can be deleted
	ivec applied_symmetry = confirm_applied_symmetry(asym_atoms);
	const int operations_before_deletion = static_cast<int>(sym[0][0].size());
	keep_operations(choose_kept_operations(asym_atoms, applied_symmetry), applied_symmetry, hkl_enlarged, hkl);
	set_sym_fact_from_position_counts(asym_atoms, operations_before_deletion);
	return applied_symmetry;
	// closing function
}

// The structure factor is right if the sym_facts of all atoms a kept operation puts on a position add up to 1 there
void cell::set_sym_fact_from_position_counts(std::vector<asym_atom>& asym_atoms, const int operations_before_deletion) {
	const int kept_operations = static_cast<int>(sym[0][0].size());
	for (int parent_idx = 0; parent_idx < asym_atoms.size(); parent_idx++) {
		if (asym_atoms[parent_idx].grown) {
			continue;
		}
		// An asymmetric atom and the atoms grown from it are the only atoms that can land on each other's positions
		ivec family = { parent_idx };
		for (int atom_idx = 0; atom_idx < asym_atoms.size(); atom_idx++) {
			if (asym_atoms[atom_idx].grown && asym_atoms[atom_idx].grown_from == parent_idx) {
				family.push_back(atom_idx);
			}
		}
		// Every distinct position a family member is put on, and how often each member lands there
		vec2 positions;
		ivec2 landings_per_position;
		for (int member = 0; member < family.size(); member++) {
			const asym_atom& atom = asym_atoms[family[member]];
			for (int sym_op = 0; sym_op < kept_operations; sym_op++) {
				const vec image = apply_symmetry({ atom.frac_pos[0], atom.frac_pos[1], atom.frac_pos[2] }, sym_op);
				int position_idx = 0;
				while (position_idx < positions.size() && !check_special(positions[position_idx], image, 1e-2)) {
					position_idx++;
				}
				if (position_idx == positions.size()) {
					positions.push_back(image);
					landings_per_position.push_back(ivec(family.size(), 0));
				}
				landings_per_position[position_idx][member]++;
			}
		}
		// All operations together put the asymmetric atom on (all operations / operations leaving it in place) positions
		err_checkf(positions.size() * asym_atoms[parent_idx].cryst_mult == operations_before_deletion,
			"The kept symmetry operations do not reach every position of atom " + asym_atoms[parent_idx].label, std::cout);
		// One condition per position: the weights of all members landing there add up to 1
		Eigen::MatrixXd landing_counts(positions.size(), family.size());
		for (int position_idx = 0; position_idx < positions.size(); position_idx++) {
			for (int member = 0; member < family.size(); member++) {
				landing_counts(position_idx, member) = landings_per_position[position_idx][member];
			}
		}
		const Eigen::VectorXd one_per_position = Eigen::VectorXd::Ones(positions.size());
		// When several weightings fulfil the conditions, the minimum norm solution shares the weight evenly between members
		const Eigen::VectorXd weights = landing_counts.completeOrthogonalDecomposition().solve(one_per_position);
		const double largest_deviation = (landing_counts * weights - one_per_position).cwiseAbs().maxCoeff();
		err_checkf(largest_deviation < 1e-8, "No sym_facts count every position of atom " + asym_atoms[parent_idx].label + " exactly once under the kept symmetry operations", std::cout);
		for (int member = 0; member < family.size(); member++) {
			asym_atoms[family[member]].sym_fact = weights[member];
		}
	}
}

int cell::equal_to_concatenation(const int op_a, const int op_b) {
	int rot[3][3];
	double t[3];
	for (int i = 0; i < 3; i++) {
		rot[i][0] = rot[i][1] = rot[i][2] = 0;
		t[i] = trans[i][op_a];
		for (int k = 0; k < 3; k++) {
			t[i] += sym[k][i][op_a] * trans[k][op_b];
			for (int j = 0; j < 3; j++)
				rot[i][j] += sym[k][i][op_a] * sym[j][k][op_b];
		}
	}
	const int num_sym_ops = static_cast<int>(sym[0][0].size());
	for (int c = 0; c < num_sym_ops; c++) {
		bool same = true;
		for (int i = 0; i < 3 && same; i++) {
			const double dt = trans[i][c] - t[i];
			same = std::abs(dt - std::round(dt)) < 1e-6;
			for (int j = 0; j < 3 && same; j++)
				same = sym[j][i][c] == rot[i][j];
		}
		if (same) return c;
	}
	return -1;
}