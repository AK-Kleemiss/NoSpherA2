#include "pch.h"
#include "cell.h"
#include "convenience.h"
#include "nos_math.h"

vec cell::apply_symmetry(const vec& pos, const int sym_op) {
	const vec trans_temp = { trans[0][sym_op], trans[1][sym_op], trans[2][sym_op] };
	const vec2 rot_temp = { { (double)sym[0][0][sym_op], (double)sym[0][1][sym_op], (double)sym[0][2][sym_op] },
							 { (double)sym[1][0][sym_op], (double)sym[1][1][sym_op], (double)sym[1][2][sym_op] },
							 { (double)sym[2][0][sym_op], (double)sym[2][1][sym_op], (double)sym[2][2][sym_op] } };
	vec temp_pos = self_dot(rot_temp, pos, false);
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
			if (check_special(pos1, pos2, 1e-4)) {
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

void cell::link_symmetry_atoms(std::vector<asym_atom>& asym_atoms, ivec3& linking_list, const int& asymmetric_atoms) {
	auto frac_pos = [&](int i) -> vec {
		return { asym_atoms[i].frac_pos[0], asym_atoms[i].frac_pos[1], asym_atoms[i].frac_pos[2] };
		};
	linking_list.resize(asymmetric_atoms);
	const int num_sym_ops = sym[0][0].size();
	int idx1 = 0;
	for (asym_atom a : asym_atoms) {
		if (a.grown) {
			idx1++;
			continue;
		}
		for (int sym_op = 0; sym_op < num_sym_ops; sym_op++) {
			const vec pos1 = apply_symmetry(frac_pos(idx1), sym_op);
			int idx2 = 0;
			for (asym_atom b : asym_atoms) {
				if (!b.grown) {
					idx2++;
					continue;
				}
				const vec pos2 = frac_pos(idx2);
				if (check_special(pos1, pos2, 1e-4)) {
					linking_list[idx1].push_back({ idx2, sym_op });
					asym_atoms[idx2].label = asym_atoms[idx1].label + "_g";
					break;
				}
				idx2++;
			}
		}
		idx1++;
	}
	// Closing function
}

void cell::eval_symm(std::vector<asym_atom>& asym_atoms, const int& asymmetric_atoms, ivec3& linking_list, const bool& grown) {
	const int ncen = asym_atoms.size();
	vec pos(3);
	vec new_pos(3);
	const int num_sym = trans[0].size();
	int idx = 0;
	for (asym_atom a : asym_atoms) {
		if (a.grown) {
			continue;
		}
		int count = 0;
		pos[0] = a.frac_pos[0];
		pos[1] = a.frac_pos[1];
		pos[2] = a.frac_pos[2];
		for (int t = 0; t < num_sym; t++) {
			const vec trans_temp = { trans[0][t], trans[1][t], trans[2][t] };
			const vec2 rot_temp = { { (double)sym[0][0][t], (double)sym[0][1][t], (double)sym[0][2][t] },
									 { (double)sym[1][0][t], (double)sym[1][1][t], (double)sym[1][2][t] },
									 { (double)sym[2][0][t], (double)sym[2][1][t], (double)sym[2][2][t] } };
			vec temp_pos = self_dot(rot_temp, pos, false);
			new_pos[0] = temp_pos[0] + trans_temp[0];
			new_pos[1] = temp_pos[1] + trans_temp[1];
			new_pos[2] = temp_pos[2] + trans_temp[2];
			if (check_special(pos, new_pos)) {
				count++;
			}
		}
		asym_atoms[idx].asym_fact = 1.0 / count;
		idx++;
	}
	if (grown) {
		link_symmetry_atoms(asym_atoms, linking_list, asymmetric_atoms);
		for (int i = 0; i < linking_list.size(); i++) {
			for (int j = 0; j < linking_list[i].size(); j++) {
				asym_atoms[linking_list[i][j][0]].asym_fact = asym_atoms[i].asym_fact;
			}
		}
	}

	//closing function
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
void cell::apply_grown(const hkl_list& hkl, hkl_list& hkl_enlarged, std::vector<asym_atom>& asym_atoms, const ivec3& linking_list) {
	const ivec applied_symmetry = confirm_applied_symmetry(asym_atoms, linking_list);
	delete_symmetry(applied_symmetry, hkl_enlarged, hkl);
	// closing function
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

ivec cell::confirm_applied_symmetry(std::vector<asym_atom>& asym_atoms, const ivec3& linking_list) {
	ivec applied_symmetry;
	const int num_sym_ops = sym[0][0].size();
	const int ncen = asym_atoms.size();
	auto frac_pos = [&](int i) -> vec {
		return { asym_atoms[i].frac_pos[0], asym_atoms[i].frac_pos[1], asym_atoms[i].frac_pos[2] };
		};

	// Generate the seperate lists of asymmetric unit atoms and grown atoms
	std::vector<asym_atom> unique_atoms;
	for (int i = 0; i < ncen; i++) {
		if (!asym_atoms[i].grown) {
			unique_atoms.push_back(asym_atoms[i]);
		}
	}
	std::vector<asym_atom> grown_atoms;
	for (int i = 0; i < ncen; i++) {
		if (asym_atoms[i].grown) {
			grown_atoms.push_back(asym_atoms[i]);
		}
	}

	for (int sym_op = 0; sym_op < num_sym_ops; sym_op++) {
		// Exclude the identity operation
		if (check_identity(sym_op)) {
			continue;
		}
		int match_count = 0;
		int idx1 = 0;
		for (asym_atom a : unique_atoms) {
			const vec pos1 = apply_symmetry(frac_pos(idx1), sym_op);
			for (asym_atom b : grown_atoms) {
				const vec pos2 = { b.frac_pos[0], b.frac_pos[1], b.frac_pos[2] };
				if (check_special(pos1, pos2, 1e-4)) {
					match_count++;
					break;
				}
			}
			// This is an additional check for atoms on special positions since these will not result in a grown atom and would break the confirmation
			const vec pos2 = { a.frac_pos[0], a.frac_pos[1], a.frac_pos[2] };
			if (check_special(pos1, pos2, 1e-4)) {
				match_count++;
			}
			idx1++;
		}
		if (match_count == unique_atoms.size()) {
			applied_symmetry.push_back(sym_op);
		}
		else if (match_count > 0) {
			std::cerr << "Warning: Symmetry operation not fully matched. Structure seems to be grown improperly!" << std::endl;
		}
	}
	return applied_symmetry;
}

void cell::delete_symmetry(const ivec& applied_symmetry, hkl_list& hkl_enlarged, const hkl_list& hkl) {
	const int nr = hkl.size();
	std::vector<i3> hkl_vec(hkl.begin(), hkl.end());
	for (int r = 0; r < nr; r++) {
		vec hkl_temp = { (double)hkl_vec[r][0], (double)hkl_vec[r][1], (double)hkl_vec[r][2] };
		for (int sym_op : applied_symmetry) {
			const vec2 rot_temp = { { (double)sym[0][0][sym_op], (double)sym[0][1][sym_op], (double)sym[0][2][sym_op] },
									{ (double)sym[1][0][sym_op], (double)sym[1][1][sym_op], (double)sym[1][2][sym_op] },
									{ (double)sym[2][0][sym_op], (double)sym[2][1][sym_op], (double)sym[2][2][sym_op] } };
			vec new_hkl = self_dot(rot_temp, hkl_temp, false);
			i3 new_hkl_int = { (int)std::round(new_hkl[0]), (int)std::round(new_hkl[1]), (int)std::round(new_hkl[2]) };
			hkl_enlarged.erase(new_hkl_int);
		}
	}
	for (int sym_op : applied_symmetry) {
		for (ivec2& middle : sym) {
			for (ivec& inner : middle) {
				inner.erase(inner.begin() + sym_op);
			}
		}
		for (vec& inner : trans) {
			inner.erase(inner.begin() + sym_op);
		}
	}
	// closing function
}