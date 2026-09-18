#include "pch.h"
#include "cell.h"
#include "convenience.h"
#include "nos_math.h"
#include <cctype>

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

void cell::convert_to_fracs(std::vector<asym_atom>& atoms, const std::string input_unit) {
	for (asym_atom& atom : atoms) {
		occ::Vec temp_cart_pos(3);
		temp_cart_pos << atom.pos[0], atom.pos[1], atom.pos[2];
		occ::Mat3N inv_cell(3, 3);
		if (input_unit == "bohr") {
			inv_cell << cm[0][0], cm[1][0], cm[2][0], cm[0][1], cm[1][1], cm[2][1], cm[0][2], cm[1][2], cm[2][2];
		}
		else if (input_unit == "angstrom") {
			inv_cell << constants::bohr2ang(cm[0][0]), constants::bohr2ang(cm[0][1]), constants::bohr2ang(cm[0][2]), constants::bohr2ang(cm[1][0]), constants::bohr2ang(cm[1][1]), constants::bohr2ang(cm[1][2]), constants::bohr2ang(cm[2][0]), constants::bohr2ang(cm[2][1]), constants::bohr2ang(cm[2][2]);
		}
		else {
			std::cerr << "Unknown input unit. Choose 'bohr' or 'angstrom'." << std::endl;
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
		const vec pos1 = { xyz_atom.frac_pos[0], xyz_atom.frac_pos[1], xyz_atom.frac_pos[2] };
		for (const asym_atom& asym_atom : asym_atoms) {
			const vec pos2 = { asym_atom.frac_pos[0], asym_atom.frac_pos[1], asym_atom.frac_pos[2] };
			// This is not a real solution, only a quick fix. Should use the same procedure as in the cif reader (scattering_factors.cpp)

			if (check_special(pos1, pos2, 1e-2)) {
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

// Handles the processing of grown structures
ivec cell::apply_grown(ivec3& linking_list) {
	return confirm_applied_symmetry(linking_list);
	// closing function
}

// Number of explicit atoms in asym_atoms that are symmetry images of the same
// asymmetric atom as atom_links belongs to, itself included. Using the full,
// undeleted symmetry operation set to rebuild each atom's orbit always visits
// each of its physical sites exactly self_links.size() times,
// uniformly, whether or not another explicit atom already covers part
// of that orbit; dividing by orbit_copies as well corrects for that overlap.
int cell::orbit_copies(const ivec2& atom_links) {
	int count = 0;
	for (const ivec& link : atom_links) {
		if (!link.empty()) count++;
	}
	return count;
}

void cell::set_symmetry_factors(std::vector<asym_atom>& asym_atoms, const ivec3& linking_list) {
	int idx1 = 0;
	for (asym_atom& a : asym_atoms) {
		if (!a.grown) {
			a.asym_fact = 1.0 / (linking_list[idx1][idx1].size() * orbit_copies(linking_list[idx1]));
			idx1++;
			continue;
		}
		for (int idx2 = 0; idx2 < linking_list.size(); idx2++) {
			if (linking_list[idx2][idx1].size() != 0) {
				a.asym_fact = 1.0 / (linking_list[idx2][idx2].size() * orbit_copies(linking_list[idx2]));
			}
		}
		idx1++;
	}
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

ivec cell::confirm_applied_symmetry(ivec3& linking_list) {
	ivec applied_symmetry;
	const int asymmetric_atoms = linking_list.size();
	const int num_sym_ops = sym[0][0].size();
	for (int sym_op = 0; sym_op < num_sym_ops; sym_op++) {
		if (check_identity(sym_op)) {
			continue;
		}
		int counter = 0;
		for (int idx1 = 0; idx1 < linking_list.size(); idx1++) {
			for (int idx2 = 0; idx2 < linking_list[idx1].size(); idx2++) {
				for (int idx3 = 0; idx3 < linking_list[idx1][idx2].size(); idx3++) {
					if (linking_list[idx1][idx2][idx3] == sym_op) {
						counter++;
					}
				}
			}
		}
		if (counter == asymmetric_atoms) {
			applied_symmetry.push_back(sym_op);
		}
		else if (counter != 0) {
			int counter2 = 0;
			for (int idx1 = 0; idx1 < linking_list.size(); idx1++) {
				for (int idx2 = 0; idx2 < linking_list[idx1][idx1].size(); idx2++) {
					if (linking_list[idx1][idx1][idx2] == sym_op) {
						counter2++;
					}
				}
			}
			if (counter != counter2) {
				std::cerr << "Warning: Symmetry operation not fully matched. Structure seems to be grown improperly!" << std::endl;
			}
		}
	}
	return applied_symmetry;
}

