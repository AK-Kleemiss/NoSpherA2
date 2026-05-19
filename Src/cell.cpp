#include "cell.h"
#include "convenience.h"
#include "wfn_class.h"
#include "nos_math.h"

std::vector<cell::asym_atom> cell::get_asym_atoms(WFN &wave, svec& labels, ivec& atom_type_list, ivec& asym_atom_to_type_list, ivec& asym_atom_list) {
	std::vector<cell::asym_atom> asym_atoms(labels.size());
	for (int i = 0; i < labels.size(); i++) {
		std::string label = labels[i];
		int type = atom_type_list[asym_atom_to_type_list[i]];
		int idx = asym_atom_list[i];
		d3 pos = { wave.get_atom(idx).get_coordinate(0), wave.get_atom(idx).get_coordinate(1), wave.get_atom(idx).get_coordinate(2) };
		d3 frac_pos = { wave.get_atom(idx).get_frac_coordinate(0), wave.get_atom(idx).get_frac_coordinate(1), wave.get_atom(idx).get_frac_coordinate(2) };
		cell::asym_atom asym_atom;
		asym_atom.label.push_back(label);
		asym_atom.type.push_back(type);
		asym_atom.pos = pos;
		asym_atom.frac_pos = frac_pos;
		asym_atom.asym_fact = 0.0;
		asym_atom.anom = 0;
		asym_atoms[i] = asym_atom;
	}
	return asym_atoms;
}

void cell::eval_symm(std::vector<cell::asym_atom>& asym_atoms) {
	const int ncen = asym_atoms.size();
	vec pos(3);
	vec new_pos(3);
	auto sym_ = sym;
	auto trans_ = trans;
	int num_sym = trans[0].size();
	int idx = 0;
	for (asym_atom a : asym_atoms) {
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
			if(check_special(pos, new_pos)) {
				count++;
			}
		}
		asym_atoms[idx].asym_fact = 1.0 / count;
		idx++;
	}
	//closing function
}

bool cell::check_special(const vec& pos1, const vec& pos2) {
	bool special = true;
	double dist1 = (pos2[0] - pos1[0]);
	double dist2 = (pos2[1] - pos1[1]);
	double dist3 = (pos2[2] - pos1[2]);
	double mod1 = std::fmod(dist1, 1);
	double mod2 = std::fmod(dist2, 1);
	double mod3 = std::fmod(dist3, 1);
	if (std::abs(mod1) > 1e-10 || std::abs(mod2) > 1e-10 || std::abs(mod3) > 1e-10) {
			special = false;
	}
	return special;
	//closing function
}