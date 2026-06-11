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

void cell::eval_symm(std::vector<asym_atom> &asym_atoms) {
    const int ncen = asym_atoms.size();
    vec pos(3);
    vec new_pos(3);
    const int num_sym = trans[0].size();
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
            if (check_special(pos, new_pos)) {
                count++;
            }
        }
        asym_atoms[idx].asym_fact = 1.0 / count;
        idx++;
    }
    //closing function
}

bool cell::check_special(const vec &pos1, const vec &pos2) {
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

// Handles the processing of grown structures
void cell::apply_grown(const hkl_list& hkl, hkl_list& hkl_enlarged, std::vector<asym_atom>& asym_atoms) {
	const ivec applied_symmetry = find_applied_symmetry(asym_atoms);
	delete_symmetry(applied_symmetry, hkl_enlarged, hkl);
	// closing function
}

ivec cell::find_applied_symmetry(std::vector<asym_atom>& asym_atoms) {
	ivec applied_symmetry;
	const int num_sym_ops = sym[0][0].size();
	const int ncen = asym_atoms.size();
	for (int a = 0; a < ncen; a++) {
        const vec pos_temp = { asym_atoms[a].frac_pos[0], asym_atoms[a].frac_pos[1], asym_atoms[a].frac_pos[2] };
        for (int sym_op = 0; sym_op < num_sym_ops; sym_op++) {
			const vec new_pos = apply_symmetry(pos_temp, sym_op);
			for (int b = a + 1; b < ncen; b++) {
				const vec pos_temp2 = { asym_atoms[b].frac_pos[0], asym_atoms[b].frac_pos[1], asym_atoms[b].frac_pos[2] };
				if (check_special(new_pos, pos_temp2)) {
					applied_symmetry.push_back(sym_op);
					break;
				}
			}
        }
	}
	
	return applied_symmetry;
	//closing function
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
            sym.erase(sym.begin() + sym_op);
        }
    }
    // closing function
}