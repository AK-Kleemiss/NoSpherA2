#include "pch.h"
#include "XCW.h"
#include "convenience.h"
#include "scattering_factors.h"
#include "nos_math.h"
#include "basis_set.h"

void XCW::construct(const options& opt_in) {
	opt = &opt_in;
	std::filesystem::path hkl_filename = opt_in.hkl;
	std::filesystem::path cif = opt_in.cif;
	std::ifstream cif_input(cif.c_str(), std::ios::in);
	unit_cell = cell(cif, std::cout, opt_in.debug, opt_in.do_XCW);
	hkl_enlarged = read_hkl_full(hkl_filename, hkl, opt_in.twin_law, unit_cell, std::cout, obs, opt_in.debug);
	std::ofstream log3("log3.txt", std::ios::out);
	bvec needs_grid;
	read_atoms_from_CIF(cif_input, unit_cell, cryst.ncen, needs_grid, asym_atoms, opt_in.debug);

	// Generate WFN object from asym_atoms
	dummy_wave.assign_charge(opt->charge);
	dummy_wave.assign_multi(opt->mult);
	for (int at = 0; at < cryst.ncen; at++) {
		asym_atom_list.push_back(at);
		atom temp_atom;
		temp_atom.set_coordinate(0, asym_atoms[at].pos[0]);
		temp_atom.set_coordinate(1, asym_atoms[at].pos[1]);
		temp_atom.set_coordinate(2, asym_atoms[at].pos[2]);
		temp_atom.set_charge(asym_atoms[at].type);
		dummy_wave.push_back_atom(temp_atom);
	}
	std::shared_ptr<BasisSet> basis = BasisSetLibrary::get_basis_set(settings.basis_set_name);
	load_basis_into_WFN(dummy_wave, basis, false, true);

	cryst.U_iso = read_U_iso_from_CIF(cif, dummy_wave, unit_cell, log3, opt_in.debug);
	unit_cell.eval_symm(asym_atoms);
	if (settings.grown) {
		unit_cell.apply_grown(hkl, hkl_enlarged, asym_atoms);
	}
	cryst.nr = hkl_enlarged.size();
	cryst.nr_small = hkl.size();
	make_k_pts(cryst.nr != 0 && hkl.size() == 0, opt_in.save_k_pts, unit_cell, hkl_enlarged, k_pt, std::cout, opt_in.debug);
	read_fracs_ADPs_from_CIF(cif, dummy_wave, log3, opt_in.debug);
	XCW_log.open("XCW.log");
	cryst.inv_scale = 1.0 / (cryst.nr_small - settings.n_params);
	F_calc.resize(2);
	for (int i = 0; i < 2; i++) {
		F_calc[i].resize(cryst.nr_small, 0);
	}
}

XCW::SCF_settings XCW::loadSettings() {
	SCF_settings settings;
	settings.quant_diff = 1e-6;
	settings.diis_stop_damping = 0.01;
	settings.max_diis_error = 1e-5;
	settings.gradient = 1e-5;
	settings.MaxP_diff = 1e-7;
	settings.RMSP_diff = 5e-9;
	settings.diis_stop_shift = 0.01;
	settings.basis_set_name = "def2-svp";
	settings.grown = false;
	settings.n_params = 161;
	// 1: Refine against F values
	// 2: Refine against F^2 values
	settings.refine_against = 1;
	settings.hf_type = occ::qm::SpinorbitalKind::Restricted;
	settings.alpha = 0.5;
	settings.level_shift = 0.0;
	settings.num_xcw_steps = 10;
	settings.xcw_step_size = 0.01;
	settings.max_scf_iterations = 100;
	settings.charge = 0;
	settings.multiplicity = 1;
	return settings;
}

void XCW::U_cif2U_star() {
	vec norm(3);
	vec2 rec_matrix(3, vec(3));
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			rec_matrix[i][j] = unit_cell.get_rcm(i, j);
		}
	}
	const double scale = constants::ang2bohr(1) / constants::TWO_PI;
	std::transform(rec_matrix.begin(), rec_matrix.end(), rec_matrix.begin(), [scale](std::vector<double>& vec) {
		std::transform(vec.begin(), vec.end(), vec.begin(), [scale](double x) { return x * scale; });
		return vec; });
	norm[0] = std::sqrt(rec_matrix[0][0] * rec_matrix[0][0] + rec_matrix[1][0] * rec_matrix[1][0] + rec_matrix[2][0] * rec_matrix[2][0]);
	norm[1] = std::sqrt(rec_matrix[0][1] * rec_matrix[0][1] + rec_matrix[1][1] * rec_matrix[1][1] + rec_matrix[2][1] * rec_matrix[2][1]);
	norm[2] = std::sqrt(rec_matrix[0][2] * rec_matrix[0][2] + rec_matrix[1][2] * rec_matrix[1][2] + rec_matrix[2][2] * rec_matrix[2][2]);
	vec transform(6);

	transform[0] = norm[0] * norm[0];
	transform[1] = norm[1] * norm[1];
	transform[2] = norm[2] * norm[2];
	transform[3] = norm[0] * norm[1];
	transform[4] = norm[0] * norm[2];
	transform[5] = norm[1] * norm[2];

	for (int a = 0; a < cryst.ncen; a++) {
		if (dummy_wave.get_atom(a).get_ADPs()[0].size() > 0) {
			vec2 ADPs = dummy_wave.get_atom(a).get_ADPs();
			for (int i = 0; i < 6; i++) {
				ADPs[0][i] *= transform[i];
			}
			dummy_wave.set_atom_ADPs(a, ADPs);
		}
	}
}

void XCW::U_star2U_cart() {
	const double scale = constants::bohr2ang(1);
	vec2 cart_matrix(3, vec(3));
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			cart_matrix[i][j] = unit_cell.get_cm(i, j);
		}
	}
	std::transform(cart_matrix.begin(), cart_matrix.end(), cart_matrix.begin(), [scale](std::vector<double>& vec) {
		std::transform(vec.begin(), vec.end(), vec.begin(), [scale](double x) { return x * scale; });
		return vec; });
	for (int a = 0; a < cryst.ncen; a++) {
		vec2 ADPs = dummy_wave.get_atom(a).get_ADPs();
		if (dummy_wave.get_atom(a).get_ADPs()[0].size() > 0) {
			vec2 U_star(3, vec(3));
			U_star[0][0] = ADPs[0][0];
			U_star[0][1] = ADPs[0][3];
			U_star[0][2] = ADPs[0][4];
			U_star[1][0] = ADPs[0][3];
			U_star[1][1] = ADPs[0][1];
			U_star[1][2] = ADPs[0][5];
			U_star[2][0] = ADPs[0][4];
			U_star[2][1] = ADPs[0][5];
			U_star[2][2] = ADPs[0][2];
			U_star = self_dot(self_dot(cart_matrix, U_star, true, false), cart_matrix, false, false);
			ADPs[0][0] = U_star[0][0];
			ADPs[0][1] = U_star[1][1];
			ADPs[0][2] = U_star[2][2];
			ADPs[0][3] = U_star[0][1];
			ADPs[0][4] = U_star[0][2];
			ADPs[0][5] = U_star[1][2];
			dummy_wave.set_atom_ADPs(a, ADPs);
		}
		if (dummy_wave.get_atom(a).get_ADPs()[1].size() > 0) {
			int running_idx = 0;
			vec C_cart(10);
			for (int i = 0; i < 3; i++) {
				for (int j = i; j < 3; j++) {
					for (int k = j; k < 3; k++) {
						double sum = 0;
						for (int p = 0; p < 3; p++) {
							for (int q = 0; q < 3; q++) {
								for (int r = 0; r < 3; r++) {
									ivec sorted_idx = { p, q, r };
									std::sort(sorted_idx.begin(), sorted_idx.end());
									int ADP_idx;
									get_voigt_index(sorted_idx, ADP_idx);
									sum += cart_matrix[p][i] * cart_matrix[q][j] * cart_matrix[r][k] * ADPs[1][ADP_idx];
								}
							}
						}
						C_cart[running_idx] = sum;
						running_idx++;
					}
				}
			}
			for (int i = 0; i < 10; i++) {
				ADPs[1][i] = C_cart[i];
			}
			dummy_wave.set_atom_ADPs(a, ADPs);
		}
		if (dummy_wave.get_atom(a).get_ADPs()[2].size() > 0) {
			int running_idx = 0;
			vec D_cart(15);
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
											get_voigt_index(sorted_idx, ADP_idx);
											sum += cart_matrix[p][i] * cart_matrix[q][j] * cart_matrix[r][k] * cart_matrix[s][l] * ADPs[2][ADP_idx];
										}
									}
								}
							}
							D_cart[running_idx] = sum;
							running_idx++;
						}
					}
				}
			}
			for (int i = 0; i < 15; i++) {
				ADPs[2][i] = D_cart[i];
			}
			dummy_wave.set_atom_ADPs(a, ADPs);
		}
	}
}

void XCW::get_voigt_index(const ivec& indices, int& ADP_idx) {
	ivec2 map3, map4;
	ivec mult3, mult4;
	map3 = { { 0, 0, 0 }, { 0, 0, 1 }, { 0, 0, 2 }, { 0, 1, 1 }, {0, 1, 2}, {0, 2, 2}, {1, 1, 1}, { 1, 1, 2 }, { 1, 2, 2 }, { 2, 2, 2 } };
	map4 = { { 0, 0, 0, 0 }, { 0, 0, 0, 1 }, { 0, 0, 0, 2 }, { 0, 0, 1, 1 }, { 0, 0, 1, 2 }, { 0, 0, 2, 2 }, { 0, 1, 1, 1 }, { 0, 1, 1, 2 }, { 0, 1, 2, 2 }, { 0, 2, 2, 2 }, { 1, 1, 1, 1 }, { 1, 1, 1, 2 }, { 1, 1, 2, 2 }, { 1, 2, 2, 2 }, { 2, 2, 2, 2 } };
	if (indices.size() == 3) {
		int idx = 0;
		while (indices != map3[idx]) {
			idx++;
		}
		ADP_idx = idx;
	}
	if (indices.size() == 4) {
		int idx = 0;
		while (indices != map4[idx]) {
			idx++;
		}
		ADP_idx = idx;
	}
}

void XCW::eval_DW(cvec2& DW_fact) {
	DW_fact.resize(cryst.ncen, cvec(cryst.nr, 0));
	//Converts angstrom to bohr OR MORE IMPORTANTLY reciprocal bohr to reciprocal angstrom
	const double angstrom2bohr = constants::ang2bohr(1);
	const double scale = angstrom2bohr / constants::TWO_PI;
	std::vector<int> level;
	level.reserve(cryst.ncen);
	//Figure out which level of anisotropic displacements parameters are avaialable
	for (int a = 0; a < cryst.ncen; a++) {
		vec2 ADPs = dummy_wave.get_atom(a).get_ADPs();
		if (ADPs.size() != 3) {
			ADPs.resize(3);
			dummy_wave.set_atom_ADPs(a, ADPs);
			level.emplace_back(0);
		}
		else if (ADPs[2].size() != 0) {
			level.emplace_back(3);
		}
		else if (ADPs[1].size() != 0) {
			level.emplace_back(2);
		}
		else if (ADPs[0].size() != 0) {
			level.emplace_back(1);
		}
		else {
			level.emplace_back(0);
		}
	}
	// Convert ADPs from cif format to Cartesian coordinates
	U_cif2U_star();
	U_star2U_cart();
	vec2 q(cryst.nr, vec(3));
	for (int h = 0; h < cryst.nr; h++) {
		q[h][0] = k_pt[0][h];
		q[h][1] = k_pt[1][h];
		q[h][2] = k_pt[2][h];
	}
	std::transform(q.begin(), q.end(), q.begin(), [scale](std::vector<double>& vec) {
		std::transform(vec.begin(), vec.end(), vec.begin(), [scale](double x) { return x * scale; });
		return vec; });
	for (int a = 0; a < cryst.ncen; a++) {
		vec2 ADPs = dummy_wave.get_atom(a).get_ADPs();
		vec2 Uij;
		if (level[a] > 0) {
			Uij = { { ADPs[0][0], ADPs[0][3], ADPs[0][4] },
						 { ADPs[0][3], ADPs[0][1], ADPs[0][5] },
						 { ADPs[0][4], ADPs[0][5], ADPs[0][2] } };
		}
		switch (level[a]) {
		case 0: {
			// Isotropic
			double U = cryst.U_iso[a], temp;
			for (int r = 0; r < cryst.nr; r++) {
				temp = -0.5 * (constants::TWO_PI * constants::TWO_PI) * U * (q[r][0] * q[r][0] + q[r][1] * q[r][1] + q[r][2] * q[r][2]);
				DW_fact[a][r] = std::exp(temp);
			}
			break;
		}
		case 1: {
			// Anisotropic U_ij
			double temp1;
			for (int h = 0; h < cryst.nr; h++) {
				vec q_ = { q[h][0], q[h][1], q[h][2] };
				temp1 = -0.5 * (constants::TWO_PI * constants::TWO_PI) * dot_BLAS(dot(Uij, q_, true), q_, false);
				DW_fact[a][h] = std::exp(temp1);
			}
			break;
		}
		case 2: {
			// Anisotropic C_ijk
			double temp1, temp2;
			for (int h = 0; h < cryst.nr; h++) {
				vec q_ = { q[h][0], q[h][1], q[h][2] };
				temp1 = -0.5 * (constants::TWO_PI * constants::TWO_PI) * dot_BLAS(dot(Uij, q_, true), q_, false);
				temp2 = -1.0 / 6.0 * (constants::TWO_PI * constants::TWO_PI * constants::TWO_PI) * (ADPs[1][0] * q_[0] * q_[0] * q_[0] + ADPs[1][6] * q_[1] * q_[1] * q_[1] + ADPs[1][9] * q_[2] * q_[2] * q_[2]
					+ 3 * ADPs[1][1] * q_[0] * q_[0] * q_[1] + 3 * ADPs[1][2] * q_[0] * q_[0] * q_[2] + 3 * ADPs[1][3] * q_[0] * q_[1] * q_[1] + 3 * ADPs[1][5] * q_[0] * q_[2] * q_[2] + 3 * ADPs[1][7] * q_[1] * q_[1] * q_[2] + 3 * ADPs[1][8] * q_[1] * q_[2] * q_[2]
					+ 6 * ADPs[1][4] * q_[0] * q_[1] * q_[2]);
				DW_fact[a][h] = std::exp(temp1) * cdouble(1, temp2);
			}
			break;
		}
		case 3: {
			// Anisotropic D_ijkl
			double temp1, temp2, temp3;
			for (int h = 0; h < cryst.nr; h++) {
				vec q_ = { q[h][0], q[h][1], q[h][2] };
				temp1 = -0.5 * (constants::TWO_PI * constants::TWO_PI) * dot_BLAS(dot(Uij, q_, true), q_, false);
				temp2 = -1.0 / 6.0 * (constants::TWO_PI * constants::TWO_PI * constants::TWO_PI) * (ADPs[1][0] * q_[0] * q_[0] * q_[0] + ADPs[1][6] * q_[1] * q_[1] * q_[1] + ADPs[1][9] * q_[2] * q_[2] * q_[2]
					+ 3 * ADPs[1][1] * q_[0] * q_[0] * q_[1] + 3 * ADPs[1][2] * q_[0] * q_[0] * q_[2] + 3 * ADPs[1][3] * q_[0] * q_[1] * q_[1] + 3 * ADPs[1][5] * q_[0] * q_[2] * q_[2] + 3 * ADPs[1][7] * q_[1] * q_[1] * q_[2] + 3 * ADPs[1][8] * q_[1] * q_[2] * q_[2]
					+ 6 * ADPs[1][4] * q_[0] * q_[1] * q_[2]);
				temp3 = (1.0 / 24.0) * (constants::TWO_PI * constants::TWO_PI * constants::TWO_PI * constants::TWO_PI) * (ADPs[2][0] * q_[0] * q_[0] * q_[0] * q_[0] + 4.0 * ADPs[2][1] * q_[0] * q_[0] * q_[0] * q_[1] + 4.0 * ADPs[2][2] * q_[0] * q_[0] * q_[0] * q_[2]
					+ 6.0 * ADPs[2][3] * q_[0] * q_[0] * q_[1] * q_[1] + 12.0 * ADPs[2][4] * q_[0] * q_[0] * q_[1] * q_[2] + 6.0 * ADPs[2][5] * q_[0] * q_[0] * q_[2] * q_[2] + 4.0 * ADPs[2][6] * q_[0] * q_[1] * q_[1] * q_[1] + 12.0 * ADPs[2][7] * q_[0] * q_[1] * q_[1] * q_[2]
					+ 12.0 * ADPs[2][8] * q_[0] * q_[1] * q_[2] * q_[2] + 4.0 * ADPs[2][9] * q_[0] * q_[2] * q_[2] * q_[2] + ADPs[2][10] * q_[1] * q_[1] * q_[1] * q_[1] + 4.0 * ADPs[2][11] * q_[1] * q_[1] * q_[1] * q_[2] + 6.0 * ADPs[2][12] * q_[1] * q_[1] * q_[2] * q_[2]
					+ 4.0 * ADPs[2][13] * q_[1] * q_[2] * q_[2] * q_[2] + ADPs[2][14] * q_[2] * q_[2] * q_[2] * q_[2]);
				DW_fact[a][h] = std::exp(temp1) * cdouble(1 + temp3, temp2);
			}
			break;
		}
		}
	}
	// closing function
}

void XCW::eval_phase(cvec2& phase_fact) {
	phase_fact.resize(cryst.ncen, cvec(cryst.nr, 0));
	const double bohr2angstrom = constants::bohr2ang(1);
	const double angstrom2bohr = constants::ang2bohr(1);
	const double scale = angstrom2bohr / constants::TWO_PI;
	cdouble exponent;
	vec2 cm = { { unit_cell.get_cm(0,0), unit_cell.get_cm(0,1), unit_cell.get_cm(0,2)},
							{ unit_cell.get_cm(1,0), unit_cell.get_cm(1,1), unit_cell.get_cm(1,2)},
							{ unit_cell.get_cm(2,0), unit_cell.get_cm(2,1), unit_cell.get_cm(2,2)} };
	std::transform(cm.begin(), cm.end(), cm.begin(), [bohr2angstrom](std::vector<double>& vec) {
		std::transform(vec.begin(), vec.end(), vec.begin(), [bohr2angstrom](double x) { return x * bohr2angstrom; });
		return vec; });
	for (int at = 0; at < cryst.ncen; at++) {
		vec pos_frac = { asym_atoms[at].frac_pos[0], asym_atoms[at].frac_pos[1], asym_atoms[at].frac_pos[2] };
		vec new_pos_cart = dot(cm, pos_frac, true);
		for (int r = 0; r < cryst.nr; r++) {
			vec q = { k_pt[0][r], k_pt[1][r], k_pt[2][r] };
			std::transform(q.begin(), q.end(), q.begin(), [scale](double x) { return x * scale; });
			exponent = cdouble(0, constants::TWO_PI * dot_BLAS(q, new_pos_cart, false));
			phase_fact[at][r] = std::exp(exponent);
		}
	}
}

void XCW::eval_translation_phase(cvec2& translation_phase) {
	translation_phase.resize(cryst.nr_small, cvec(unit_cell.get_trans()[0].size(), 0));
	const double angstrom2bohr = constants::ang2bohr(1);
	const double bohr2angstrom = constants::bohr2ang(1);
	const double scale = angstrom2bohr / constants::TWO_PI;
	vec2 trans = unit_cell.get_trans();
	vec2 cm = { { unit_cell.get_cm(0,0), unit_cell.get_cm(0,1), unit_cell.get_cm(0,2)},
								  { unit_cell.get_cm(1,0), unit_cell.get_cm(1,1), unit_cell.get_cm(1,2)},
								  { unit_cell.get_cm(2,0), unit_cell.get_cm(2,1), unit_cell.get_cm(2,2)} };
	std::transform(cm.begin(), cm.end(), cm.begin(), [bohr2angstrom](std::vector<double>& vec) {
		std::transform(vec.begin(), vec.end(), vec.begin(), [bohr2angstrom](double x) { return x * bohr2angstrom; });
		return vec; });
	for (int r = 0; r < cryst.nr_small; r++) {
		ivec asym_list = generate_asym_lookup(r);
		vec q_temp = { k_pt[0][asym_list[0]], k_pt[1][asym_list[0]], k_pt[2][asym_list[0]] };
		std::transform(q_temp.begin(), q_temp.end(), q_temp.begin(), [scale](double x) { return x * scale; });
		for (int t = 0; t < trans[0].size(); t++) {
			vec trans_temp = { trans[0][t], trans[1][t], trans[2][t] };
			trans_temp = dot(cm, trans_temp, true);
			cdouble exponent(0, constants::TWO_PI * dot_BLAS(q_temp, trans_temp, false));
			translation_phase[r][t] = std::exp(exponent);
		}
	}
	// closing function
}

void XCW::parse_anom_atoms(std::vector<anom_atom> &anom_atoms) {
    std::ifstream file(opt->anom_disp_path);
    if (!file) {
        std::cout << "Could not open anomalous dispersion file. Continuing without anomalous dispersions." << std::endl;
    }
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty())
            continue;
        std::istringstream iss(line);
        std::string symbol;
        double real_part, imag_part;
        if (iss >> symbol >> real_part >> imag_part) {
            if (!symbol.empty() && symbol[0] != '_' && symbol != "loop_") {
                anom_atoms.push_back({ symbol, cdouble(real_part, imag_part) });
            }
        }
    }
}

void XCW::eval_anom_disp(cvec2& DW_fact, cvec2& phase_fact, cvec2& translation_phase) {
	std::vector<anom_atom> anom_atoms;
	parse_anom_atoms(anom_atoms);
	int r, at, r_asym;
	for (int at = 0; at < cryst.ncen; at++) {
		const char* symbol = constants::atnr2letter(asym_atoms[at].type);
		for (const anom_atom& anom_atom : anom_atoms) {
			if (symbol == anom_atom.identifier) {
				asym_atoms[at].anom = anom_atom.dispersion;
				break;
			}
		}
	}
	std::vector<ivec> asym_lookup(cryst.nr_small);
	for (r = 0; r < cryst.nr_small; r++) {
		asym_lookup[r] = generate_asym_lookup(r);
	}
	for (r = 0; r < cryst.nr_small; r++) {
		const ivec& lookup = asym_lookup[r];
		for (at = 0; at < cryst.ncen; at++) {
			cdouble temp1 = 0;
			for (r_asym = 0; r_asym < lookup.size(); r_asym++) {
				temp1 += phase_fact[at][lookup[r_asym]] * DW_fact[at][lookup[r_asym]] * translation_phase[r][r_asym];
			}
			F_calc[1][r] += temp1 * asym_atoms[at].asym_fact * asym_atoms[at].anom;
		}
	}
}

void XCW::eval_scale() {
    double numerator = 0.0;
    double denominator = 0.0;

#pragma omp parallel for reduction(+:numerator, denominator)
	for (int i = 0; i < cryst.nr_small; ++i) {
		const double calc = std::abs(F_calc[0][i]);
		numerator += calc * obs[i].F_obs;
		denominator += calc * calc;
	}

	cryst.F_scale = (denominator != 0.0) ? numerator / denominator : 1.0;
}

void XCW::calc_criteria() {
	double prefactor = 1.0 / static_cast<double>(cryst.nr_small - settings.n_params);
	double sum_chi = 0;
	double sum_goof = 0;
#pragma omp parallel for reduction(+:sum_chi, sum_goof)
	for (int i = 0; i < cryst.nr_small; i++) {
		const double scaled_F_calc = cryst.F_scale * std::abs(F_calc[0][i]);
		const double diff = scaled_F_calc - obs[i].abs_F_obs;
		const double diff2 = scaled_F_calc * scaled_F_calc - obs[i].F_obs2;
		const double weighted_chi = diff * diff / obs[i].sigma_obs;
		const double weighted_goof = diff2 * diff2 / (obs[i].sigma_obs2 * obs[i].sigma_obs2);
		sum_chi += weighted_chi;
		sum_goof += weighted_goof;
	}
	cryst.chi2 = prefactor * sum_chi;
	cryst.GooF = std::sqrt(prefactor * sum_goof);
}

void XCW::create_prims(std::vector<ao_data>& ao_data_shells, occ::qm::AOBasis& occ_basis_set) {
	for (int atm = 0; atm < cryst.ncen; atm++) {
		d3 pos = { occ_basis_set.atoms()[atm].x, occ_basis_set.atoms()[atm].y, occ_basis_set.atoms()[atm].z };
		const int first_shell = *occ_basis_set.atom_to_shell()[atm].begin();
		const int last_shell = occ_basis_set.atom_to_shell()[atm].back();
		for (int shell = first_shell; shell <= last_shell; shell++) {
			occ::gto::Shell current_shell = occ_basis_set.shells()[shell];
			const int shell_type = current_shell.l;
			std::vector<primitive> tmp_prims;
			for (int prim_idx = 0; prim_idx < current_shell.exponents.size(); prim_idx++) {
				const double alpha = current_shell.exponents[prim_idx];
				const double coeff = current_shell.contraction_coefficients(prim_idx);
				tmp_prims.emplace_back(0, shell_type, alpha, coeff);
			}
			for (int m = -shell_type; m <= shell_type; m++) {
				ao_data_shells.push_back({ tmp_prims, pos, m });
			}
		}
	}
	cryst.nmo = ao_data_shells.size();
}

ivec XCW::generate_asym_lookup(const int r) {
    ivec asym_list;
    auto it = hkl.begin();
    std::advance(it, r);
    ivec3 rots = unit_cell.get_sym();
    i3 tempv;
    const i3 &hkl_temp = *it;
    for (int s = 0; s < rots[0][0].size(); s++) {
        tempv = { 0, 0, 0 };
        for (int h = 0; h < 3; h++) {
            for (int j = 0; j < 3; j++) {
                tempv[j] += hkl_temp[h] * rots[j][h][s];
            }
        }
        int idx_ = 0;
        auto idx = hkl_enlarged.find(tempv);
        if (idx != hkl_enlarged.end()) {
            idx_ = std::distance(hkl_enlarged.begin(), idx);
        }
        asym_list.push_back(idx_);
    }
    return asym_list;
    // closing function
}

cvec2 XCW::calculateXCWintegral(GridManager &grid_manager, const ao_data &mu_data, const ao_data &nu_data, const ivec &asym_atom_list, vec2 &k_pt, bool &equal, vec2 &mu_vals) {
    GridData &GD = grid_manager.getGridData();
    vec2 *grids = grid_manager.getNeedsHelper() ? GD.helper_grids.data() : GD.atomic_grids.data();
    const int n_grids = grid_manager.getNeedsHelper() ? GD.helper_grids.size() : GD.atomic_grids.size();
    const int *points = grid_manager.getNeedsHelper() ? GD.helper_num_points_per_atom.data() : GD.num_points_per_atom.data();
    const double mp0 = mu_data.pos[0];
    const double mp1 = mu_data.pos[1];
    const double mp2 = mu_data.pos[2];
    const double np0 = nu_data.pos[0];
    const double np1 = nu_data.pos[1];
    const double np2 = nu_data.pos[2];
    const std::vector<primitive> mu_prims = mu_data.prims;
    const std::vector<primitive> nu_prims = nu_data.prims;
    if (equal) {
        mu_vals.resize(n_grids);
        for (int g = 0; g < n_grids; g++) {
            const int num_points = points[g];
            mu_vals[g].resize(num_points);
            vec2 &atom_grid = grids[g];
            const double *x_ptr = atom_grid[GridData::GridIndex::X].data();
            const double *y_ptr = atom_grid[GridData::GridIndex::Y].data();
            const double *z_ptr = atom_grid[GridData::GridIndex::Z].data();
            double *overlap_ptr = atom_grid[GridData::GridIndex::WFN_DENSITY].data();
#pragma omp parallel for schedule(dynamic, 4)
			for (int p = 0; p < num_points; p++) {
				std::array<double, 4> d_mu;
				double* d_mu_ptr = d_mu.data();
				*d_mu_ptr = x_ptr[p] - mp0;
				*(d_mu_ptr + 1) = y_ptr[p] - mp1;
				*(d_mu_ptr + 2) = z_ptr[p] - mp2;
				*(d_mu_ptr + 3) = std::hypot(*(d_mu_ptr), *(d_mu_ptr + 1), *(d_mu_ptr + 2));
				mu_vals[g][p] = dummy_wave.eval_ao(d_mu, mu_prims, mu_data.m);
				overlap_ptr[p] = std::pow(mu_vals[g][p], 2);
			}
		}
	}
	else {
		for (int g = 0; g < n_grids; g++) {
			const int num_points = points[g];
			vec2& atom_grid = grids[g];
			const double* x_ptr = atom_grid[GridData::GridIndex::X].data();
			const double* y_ptr = atom_grid[GridData::GridIndex::Y].data();
			const double* z_ptr = atom_grid[GridData::GridIndex::Z].data();
			double* overlap_ptr = atom_grid[GridData::GridIndex::WFN_DENSITY].data();
#pragma omp parallel for schedule(dynamic, 4)
			for (int p = 0; p < num_points; p++) {
				std::array<double, 4> d_mu, d_nu;
				double* d_mu_ptr = d_mu.data();
				double* d_nu_ptr = d_nu.data();
				*d_mu_ptr = x_ptr[p] - mp0;
				*(d_mu_ptr + 1) = y_ptr[p] - mp1;
				*(d_mu_ptr + 2) = z_ptr[p] - mp2;
				*(d_mu_ptr + 3) = std::hypot(*(d_mu_ptr), *(d_mu_ptr + 1), *(d_mu_ptr + 2));
				*d_nu_ptr = x_ptr[p] - np0;
				*(d_nu_ptr + 1) = y_ptr[p] - np1;
				*(d_nu_ptr + 2) = z_ptr[p] - np2;
				*(d_nu_ptr + 3) = std::hypot(*(d_nu_ptr), *(d_nu_ptr + 1), *(d_nu_ptr + 2));
				overlap_ptr[p] = (*d_mu_ptr + 3) > 12 || (*d_nu_ptr + 3) > 12 ? 0 : mu_vals[g][p] * dummy_wave.eval_ao(d_nu, nu_prims, nu_data.m);
			}
		}
	}
	vec2 d1, d2, d3, dens;
	std::vector<_time_point> time_points({ get_time() });
	_time_point end;
	grid_manager.getDensityVectors(dummy_wave, asym_atom_list, d1, d2, d3, dens);
	const int points_ = grid_manager.getTotalGridPoints();
	cvec2 sf;
	calc_SF(points_, k_pt, d1, d2, d3, dens, sf, std::cout, time_points.front(), end, opt->debug, true, true);
	return sf;
}

size_t XCW::tri_index(int mu, int nu) const noexcept {
	return mu * cryst.nmo - (mu * (mu - 1)) / 2 + (nu - mu);
}

size_t XCW::flattened_idx(int r, int mu, int nu) const noexcept {
	return r * cryst.nmo * (cryst.nmo + 1) / 2 + tri_index(mu, nu);
}

void XCW::eval_I_anom_disp(std::vector<ao_data>& ao_data_shells, bool read) {
	cvec2 DW_fact, phase_fact, translation_phase;
	eval_phase(phase_fact);
	eval_DW(DW_fact);
	eval_translation_phase(translation_phase);
	if (read) {
		std::ifstream in("I_tensor", std::ios::binary);
		if (!in)
			throw std::runtime_error("Cannot open file for reading");
		int nr_safe;
		int nmo_safe;
		int num_elements_safe;
		int total_size_safe;
		in.read(reinterpret_cast<char*>(&nr_safe), sizeof(nr_safe));
		in.read(reinterpret_cast<char*>(&nmo_safe), sizeof(nmo_safe));
		in.read(reinterpret_cast<char*>(&num_elements_safe), sizeof(num_elements_safe));
		in.read(reinterpret_cast<char*>(&total_size_safe), sizeof(total_size_safe));
		I.resize(total_size_safe);
		in.read(reinterpret_cast<char*>(I.data()),
			total_size_safe * sizeof(cdouble));
	}
	else {
		eval_I(ao_data_shells, DW_fact, phase_fact, translation_phase);
	}
	eval_anom_disp(DW_fact, phase_fact, translation_phase);
	// closing function
}

void XCW::eval_I(std::vector<ao_data>& ao_data_shells, cvec2& DW_fact, cvec2& phase_fact, cvec2& translation_phase) {

	I.resize(cryst.nr_small * (cryst.nmo * (cryst.nmo + 1)) / 2);
	int at = 0, mu = 0, nu = 0, r = 0, s = 0, r_asym = 0;

    cvec2 XCW_integral;

	GridConfiguration config;
	config.accuracy = opt->accuracy;
	config.partition_type = opt->partition_type;
	config.no_density_eval = true;
	config.pbc = opt->pbc;
	config.debug = opt->debug;
	config.all_charges = opt->all_charges;
	GridManager grid_manager(config);
	dummy_wave.delete_unoccupied_MOs();
	bvec needs_grid(cryst.ncen, true);
	grid_manager.setup3DGridsForMolecule(dummy_wave, asym_atom_list, needs_grid, unit_cell);
	int cofs;
	bool equal = false;
	double cutoff = 0;
	std::vector<ivec> asym_lookup(cryst.nr_small);
	for (r = 0; r < cryst.nr_small; r++) {
		asym_lookup[r] = generate_asym_lookup(r);
	}
	for (mu = 0; mu < cryst.nmo; mu++) {
		vec2 mu_vals;
		const ao_data& mu_prims = ao_data_shells[mu];
		const std::vector<primitive>& mu_primitives = mu_prims.prims;
		const double& mp0 = mu_prims.pos[0];
		const double& mp1 = mu_prims.pos[1];
		const double& mp2 = mu_prims.pos[2];
		for (nu = mu; nu < cryst.nmo; nu++) {
			if (nu == mu) {
				equal = true;
			}
			else {
				equal = false;
			}
			const ao_data& nu_prims = ao_data_shells[nu];
			const std::vector<primitive>& nu_primitives = nu_prims.prims;
			const double& np0 = nu_prims.pos[0];
			const double& np1 = nu_prims.pos[1];
			const double& np2 = nu_prims.pos[2];
			// Screening
			const double dist0 = mp0 - np0;
			const double dist1 = mp1 - np1;
			const double dist2 = mp2 - np2;
			const double dist = std::sqrt(dist0 * dist0 + dist1 * dist1 + dist2 * dist2);
			if (!equal) {
				double c_tot = 0;
				for (cofs = 0; cofs < mu_primitives.size(); cofs++) {
					double alpha = mu_primitives[cofs].get_exp();
					double beta = nu_primitives[cofs].get_exp();
					double alpha_k_l = alpha * beta / (alpha + beta);
					c_tot += std::abs(mu_primitives[cofs].get_coef()) * std::abs(nu_primitives[cofs].get_coef()) * std::pow(constants::PI / alpha_k_l, 1.5);
				}
				c_tot *= std::sqrt((2 * mu_primitives[0].get_type() + 1) * (2 * nu_primitives[0].get_type() + 1) / (constants::TWO_PI * constants::TWO_PI));
				double alpha_min = mu_primitives[mu_primitives.size() - 1].get_exp();
				double beta_min = nu_primitives[nu_primitives.size() - 1].get_exp();
				double a_min = (alpha_min + beta_min) / (alpha_min * beta_min);
				const double screen = 0.001;
				cutoff = std::sqrt(a_min * std::log(c_tot / screen));
			}
			if (dist > cutoff && !equal) {
				continue;
			}

			cdouble temp1;
			double asym_fact;
			int idx;
			XCW_integral = calculateXCWintegral(grid_manager, mu_prims, nu_prims, asym_atom_list, k_pt, equal, mu_vals);
			for (r = 0; r < cryst.nr_small; r++) {
				const cvec& translation_phase_r = translation_phase[r];
				const ivec& asym_list = asym_lookup[r];
				for (at = 0; at < cryst.ncen; at++) {
					asym_fact = asym_atoms[at].asym_fact;
					const cvec& XCW_at = XCW_integral[at];
					const cvec& DW_at = DW_fact[at];
					const cvec& phase_at = phase_fact[at];
					temp1 = 0;
					for (r_asym = 0; r_asym < asym_list.size(); r_asym++) {
						idx = asym_list[r_asym];
						temp1 += XCW_at[idx] * DW_at[idx] * phase_at[idx] * translation_phase_r[r_asym];
					}
					I[flattened_idx(r, mu, nu)] += temp1 * asym_fact;
				}
			}
		}
		// Print out that contributions from mu are done
		std::cout << "Finished contributions from mu = " << mu << std::endl;
	}
}

void XCW::calc_F_calc(const dMatrix2& D) {
	// Density matrix from occ is half of what I need, so times 2 and times (2x2)=4
#pragma omp parallel for schedule(static)
	for (int r = 0; r < cryst.nr_small; ++r) {
		cdouble sum = F_calc[1][r];
		size_t base = r * (cryst.nmo * (cryst.nmo + 1)) / 2;
		for (int mu = 0; mu < cryst.nmo; mu++) {
			sum += 2.0 * I[base] * D(mu, mu);
			base++;
			for (int nu = mu + 1; nu < cryst.nmo; nu++, base++) {
				sum += 4.0 * I[base] * D(mu, nu);
			}
		}
		F_calc[0][r] = sum;
	}
}

void XCW::calc_perturb(occ::Mat& perturb, const occ::qm::SCF<occ::qm::HartreeFock>& scf) {
	switch (settings.refine_against) {
	case 1: {
		perturb.setZero(cryst.nmo, cryst.nmo);
		const double prefactor = 2.0 * cryst.F_scale / (cryst.nr_small - settings.n_params);
		const int packed_size = cryst.nmo * (cryst.nmo + 1) / 2;
#pragma omp parallel
		{
			occ::Mat local = occ::Mat::Zero(cryst.nmo, cryst.nmo);
#pragma omp for nowait
			for (int r = 0; r < cryst.nr_small; r++) {
				const double F_calc_abs = std::abs(F_calc[0][r]);
				// There should be another sigma but somehow that does not seem to work
				const cdouble precompute = std::conj(F_calc[0][r]) * (cryst.F_scale * F_calc_abs - obs[r].abs_F_obs) / (obs[r].sigma_obs * F_calc_abs);
				const size_t base = r * packed_size;
				for (int mu = 0; mu < cryst.nmo; mu++) {
					size_t offset = base + tri_index(mu, mu);
					for (int nu = mu; nu < cryst.nmo; ++nu) {
						const double temp = std::real(precompute * I[offset]);
						local(mu, nu) += temp;
						++offset;
					}
				}
			}
#pragma omp critical
			{
				perturb += local;
			}
		}
		perturb *= prefactor;
		for (int mu = 0; mu < cryst.nmo; mu++) {
			for (int nu = mu + 1; nu < cryst.nmo; nu++) {
				perturb(nu, mu) = perturb(mu, nu);
			}
		}
		if (scf.ctx.mo.kind == occ::qm::SpinorbitalKind::Unrestricted) {
			perturb.conservativeResize(2 * cryst.nmo, Eigen::NoChange);
			perturb.bottomRows(cryst.nmo) = perturb.topRows(cryst.nmo);
		}
		break;
	}
	case 2: {
		perturb.setZero(cryst.nmo, cryst.nmo);
		const double scale_sq = cryst.F_scale * cryst.F_scale;
		const double prefactor = 4.0 * scale_sq / (cryst.nr_small - settings.n_params);
		const int packed_size = cryst.nmo * (cryst.nmo + 1) / 2;
#pragma omp parallel
		{
			occ::Mat local = occ::Mat::Zero(cryst.nmo, cryst.nmo);
#pragma omp for nowait
			for (int r = 0; r < cryst.nr_small; r++) {
				const double F_calc_abs = std::abs(F_calc[0][r]);
				// There should be another sigma but somehow that does not seem to work
				const cdouble precompute = std::conj(F_calc[0][r]) * (scale_sq * F_calc_abs * F_calc_abs - obs[r].F_obs2) / (obs[r].sigma_obs2);
				const size_t base = r * packed_size;
				for (int mu = 0; mu < cryst.nmo; mu++) {
					size_t offset = base + tri_index(mu, mu);
					for (int nu = mu; nu < cryst.nmo; ++nu) {
						const double temp = std::real(precompute * I[offset]);
						local(mu, nu) += temp;
						++offset;
					}
				}
			}
#pragma omp critical
			{
				perturb += local;
			}
		}
		perturb *= prefactor;
		for (int mu = 0; mu < cryst.nmo; mu++) {
			for (int nu = mu + 1; nu < cryst.nmo; nu++) {
				perturb(nu, mu) = perturb(mu, nu);
			}
		}
		if (scf.ctx.mo.kind == occ::qm::SpinorbitalKind::Unrestricted) {
			perturb.conservativeResize(2 * cryst.nmo, Eigen::NoChange);
			perturb.bottomRows(cryst.nmo) = perturb.topRows(cryst.nmo);
		}
		break;
	}
	default:
		XCW_log << "Invalid refinement option" << std::endl;
	}
}

void XCW::setup_SCF_mol(occ::core::Molecule &mol) {
    double bohr2angstrom = constants::bohr2ang(1);
    std::ostringstream init_stream;

	init_stream << cryst.ncen << "\n\n";

	for (int i = 0; i < cryst.ncen; ++i) {
		init_stream
			<< constants::atnr2letter(asym_atoms[i].type) << " "
			<< asym_atoms[i].pos[0] * bohr2angstrom << " "
			<< asym_atoms[i].pos[1] * bohr2angstrom << " "
			<< asym_atoms[i].pos[2] * bohr2angstrom;

		if (i != cryst.ncen - 1)
			init_stream << "\n";
	}

    std::string init = init_stream.str();
    mol = occ::io::molecule_from_xyz_string(init);
    mol.set_charge(opt->charge);
    mol.set_multiplicity(opt->mult);
}

void XCW::setup_basis(occ::core::Molecule& mol, std::string& basis_set_name, occ::qm::AOBasis& occ_basis_set) {
	std::shared_ptr<BasisSet> basis_set = BasisSetLibrary::get_basis_set(basis_set_name);
	occ_basis_set = basis_set->to_AOBasis(mol.atoms());
}

double XCW::dynamic_damping(const occ::qm::SCF<occ::qm::HartreeFock>& scf, const double& current_alpha, const double& e_diff, double& e_diff_mem) {
	double new_alpha = current_alpha;
	if (e_diff < e_diff_mem / 10) {
		new_alpha *= 0.75;
		e_diff_mem = e_diff;
		if (e_diff < 10 * scf.convergence_settings.energy_threshold) {
			print_centered_message("***Turned off damping***", 76, XCW_log);
			new_alpha = 0;
		}
		else {
			std::stringstream print_;
			print_ << "***Decreased damping to " << std::fixed << std::setprecision(3) << new_alpha << "***";
			print_centered_message(print_.str(), 76, XCW_log);
		}
	}
	return new_alpha;
	// closing function
}

void XCW::apply_level_shift(const occ::Mat& C_old, const occ::qm::SCF<occ::qm::HartreeFock>& scf, occ::Mat& F_diis) {
	const double temp_shift = scf.convergence_settings.effective_level_shift(scf.diis_error);
	const int nocc = scf.ctx.mo.Cocc.cols();
	if (scf.ctx.mo.kind == occ::qm::SpinorbitalKind::Restricted) {
		const occ::Mat SC_virt = scf.ctx.S * C_old.rightCols(cryst.nmo - nocc);
		F_diis.noalias() += temp_shift * SC_virt * SC_virt.transpose();
	}
	else {
		const int nao = C_old.rows() / 2;
		const auto S_ao = scf.ctx.S.topRows(nao);
		const occ::Mat SC_virt_a = S_ao * C_old.topRows(nao).rightCols(cryst.nmo - nocc);
		const occ::Mat SC_virt_b = S_ao * C_old.bottomRows(nao).rightCols(cryst.nmo - nocc);
		F_diis.topRows(nao).noalias() += temp_shift * SC_virt_a * SC_virt_a.transpose();
		F_diis.bottomRows(nao).noalias() += temp_shift * SC_virt_b * SC_virt_b.transpose();
	}
}

void XCW::build_effective_dm(const occ::qm::SCF<occ::qm::HartreeFock>& scf, dMatrix2& dm_ref, const occ::Mat& dm_old) {
	if (scf.ctx.mo.kind == occ::qm::SpinorbitalKind::Unrestricted) {
		for (int i = 0; i < cryst.nmo; i++) {
			dm_ref(i, i) = dm_old(i, i);
			dm_ref(i, i) += dm_old(i + cryst.nmo, i);
			for (int j = i + 1; j < cryst.nmo; j++) {
				dm_ref(i, j) = dm_old(i, j);
				dm_ref(i, j) += dm_old(i + cryst.nmo, j);
				dm_ref(j, i) = dm_ref(i, j);
			}
		}
	}
	else {
		for (int i = 0; i < cryst.nmo; i++) {
			dm_ref(i, i) = dm_old(i, i);
			for (int j = i + 1; j < cryst.nmo; j++) {
				dm_ref(i, j) = dm_old(i, j);
				dm_ref(j, i) = dm_ref(i, j);
			}
		}
	}
}

void XCW::do_SCF(const double& lambda, double& alpha, occ::qm::SCF<occ::qm::HartreeFock>& scf, occ::qm::Wavefunction& last_wfn, bool& has_guess) {

	settings.clear();

	XCW_log << "Starting XCW SCF solver with lambda = " << std::fixed << std::setprecision(5) << lambda << "\n";
	XCW_log << "____________________________________________________________________________\n";
	XCW_log << " Iteration	Chi^2	GooF	Total Energy	Perturbation	Target quantity \n";
	XCW_log << "								(Eh)		   (a. u.)			(a. u.)\n";
	XCW_log << "____________________________________________________________________________\n";

	// Compute first guess and update the energy according to this guess
	if (has_guess) {
		scf.set_initial_guess_from_wfn(last_wfn);
	}
	else {
		scf.compute_initial_guess();
		has_guess = true;
	}
	scf.ctx.K = scf.m_procedure.compute_schwarz_ints();
	scf.update_scf_energy(false);

	scf.ctx.H = scf.ctx.T + scf.ctx.V;
	bool converged;
	double quant;
	double last_quant = 0;
	double e_diff_mem = 0;
	occ::Mat dm_last = scf.ctx.mo.D;

    do {

		converged = SCF_iteration(scf, lambda, alpha, e_diff_mem, quant, last_quant, dm_last);

    } while (!converged && scf.iter < scf.maxiter);

    if (converged) {
        XCW_log << "____________________________________________________________________________\n";
        std::stringstream print_;
        print_ << "***SCF converged in " << scf.iter + 1 << " iterations***";
        print_centered_message(print_.str(), 76, XCW_log);

		std::cout << "\t" << std::fixed << std::setprecision(3) << lambda << "\t" << std::fixed << std::setprecision(3) << cryst.chi2 << "\t" << cryst.GooF << "\t" << std::fixed << std::setprecision(9) << scf.ctx.energy["total"] << "\t\t" << std::fixed << std::setprecision(3) << lambda * cryst.chi2 << "\t\t" << std::fixed << std::setprecision(9) << quant << std::endl;

		create_tscb(scf, lambda);
	}
	else {
		XCW_log << "____________________________________________________________________________\n";
		print_centered_message("***SCF did not converge***", 76, XCW_log);
	}
	// closing function
}

double XCW::compute_orbital_gradient(const occ::qm::SCF<occ::qm::HartreeFock>& scf) {
	occ::Mat C = scf.molecular_orbitals().C;
	occ::Mat Cocc = scf.molecular_orbitals().Cocc;
	occ::Mat Cvir = C.rightCols(C.cols() - Cocc.cols());
	occ::Mat G = 2.0 * Cvir.transpose() * scf.ctx.F * Cocc;
	return G.norm();
	// closing funciton
}

void XCW::get_density_criteria(double& RMSP_diff, double& maxP_diff, const occ::Mat& dm, const occ::Mat& dm_last) {
	occ::Mat difference = dm - dm_last;
	RMSP_diff = std::sqrt(difference.squaredNorm() / difference.size());
	maxP_diff = difference.cwiseAbs().maxCoeff();
	// closing function
}

bool XCW::SCF_iteration(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double& lambda, double& alpha, double& e_diff_mem, double& quant, double& last_quant, occ::Mat& dm_last) {
	// Set up energy values & crystallographic information
	scf.iter++;
	const double ehf_last = scf.ctx.energy["electronic"];
	const occ::Mat dm_old = scf.ctx.mo.D;
	dMatrix2 dm_eff(cryst.nmo, cryst.nmo);
	build_effective_dm(scf, dm_eff, dm_old);
	calc_F_calc(dm_eff);
	eval_scale();
	calc_criteria();

	// Generates the perturbation matrix
	occ::Mat perturbation;
	calc_perturb(perturbation, scf);

	// Build perturbed Fock matrix
	// Maybe necessary to update the Hamiltoian if a potential changes depending on the density, but that does not happen in normal HF
	//scf.ctx.H = scf.ctx.T + scf.ctx.V + scf.ctx.Vecp + scf.ctx.V_ext;
	scf.m_procedure.update_core_hamiltonian(scf.ctx.mo, scf.ctx.H);
	scf.ctx.F = scf.ctx.H;
	scf.ctx.F += scf.m_procedure.compute_fock(scf.ctx.mo, scf.ctx.K);
	scf.update_scf_energy(false);
	const double ehf = scf.ctx.energy["electronic"];
	const double e_diff = std::abs(ehf - ehf_last);
	quant = scf.ctx.energy["total"] + lambda * cryst.chi2;
	scf.ctx.F += perturbation * lambda;

	// Prints output line for iteration
	XCW_log << "\t" << scf.iter << "\t\t" << std::fixed << std::setprecision(3) << cryst.chi2 << "\t" << cryst.GooF << "\t" << std::fixed << std::setprecision(9) << scf.ctx.energy["total"] << "\t\t" << std::fixed << std::setprecision(3) << lambda * cryst.chi2 << "\t\t" << std::fixed << std::setprecision(9) << quant << std::endl;

	// DIIS extrapolation
	occ::Mat F_diis = scf.convergence_accelerator.update(scf.ctx.mo.kind, scf.ctx.S, scf.ctx.mo.D, scf.ctx.F, scf.ctx.energy["electronic"]);
	scf.diis_error = scf.convergence_accelerator.max_error();
	settings.update(scf.diis_error, XCW_log, alpha);

	// Convergence check
	const double gradient = compute_orbital_gradient(scf);
	const double quant_diff = std::abs(quant - last_quant);
	if (SCF_convergence_check(quant_diff, gradient, scf, dm_last)) {
		return true;
	}
	last_quant = quant;

	// Apply level shift
	if (settings.apply_shift) {
		const occ::Mat& C_old = scf.ctx.mo.C;
		apply_level_shift(C_old, scf, F_diis);
	}

	// Solves central eigenvalue problem
	scf.ctx.orthogonalizer.orthogonalize_molecular_orbitals(scf.ctx.mo, F_diis);

	// Apply damping
	if (settings.apply_damping) {
		if (scf.iter == 2) {
			e_diff_mem = e_diff;
		}
		if (scf.iter > 2) {
			alpha = dynamic_damping(scf, alpha, e_diff, e_diff_mem);
		}
	}

	scf.ctx.mo.D *= (1 - alpha);
	scf.ctx.mo.D += alpha * dm_old;
	dm_last = scf.ctx.mo.D;
	return false;

    // Apply damping
    scf.ctx.mo.D *= (1 - alpha);
    scf.ctx.mo.D += alpha * dm_old;

    return false;

    //closing function
}

bool XCW::SCF_convergence_check(const double& quant_diff, const double& gradient, occ::qm::SCF<occ::qm::HartreeFock>& scf, occ::Mat& dm_last) {
	double RMSP_diff, maxP_diff;
	get_density_criteria(RMSP_diff, maxP_diff, scf.ctx.mo.D, dm_last);
	if (quant_diff < settings.quant_diff) {
		settings.conv_quant_diff = true;
	}
	if (scf.diis_error < settings.max_diis_error) {
		settings.conv_max_diis_error = true;
	}
	if (gradient < settings.gradient) {
		settings.conv_gradient = true;
	}
	if (RMSP_diff < settings.RMSP_diff) {
		settings.conv_RMSP_diff = true;
	}
	if (maxP_diff < settings.MaxP_diff) {
		settings.conv_MaxP_diff = true;
	}
	return settings.convergence_check();
	// closing function
}

void XCW::create_tscb(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double& lambda) {
	XCW_log << "Creating .tscb file from converged SCF calculation..." << std::endl;
	std::vector<WFN> sf_wave_vec(1, { scf.wavefunction(), false });
	svec known_atoms_;
	tsc_block<int, cdouble> result;
	vec2 known_kpts_;
	options* opt_ = const_cast<options*>(opt);
	result.append(calculate_scattering_factors<itsc_block, std::vector<WFN>&>(
		*opt_,
		sf_wave_vec,
		XCW_log,
		known_atoms_,
		0,
		&k_pt),
		XCW_log);
	int value = static_cast<int>(std::round(lambda * 100));
	std::ostringstream oss;
	oss << "NA2_" << std::setw(3) << std::setfill('0') << value << ".tscb";
	result.write_tscb_file("test.cif", oss.str());
	std::ostringstream oss2;
	oss2 << "NA2_" << std::setw(3) << std::setfill('0') << value << ".wfn";
	sf_wave_vec[0].write_wfn(oss2.str(), false, true);
	//Roby_information Roby(sf_wave_vec[0]);
}

occ::qm::HartreeFock XCW::setup_XCW_procedure(bool read, bool safe) {
	std::vector<ao_data> ao_data_shells;
	occ::core::Molecule mol;
	setup_SCF_mol(mol);
	occ::qm::AOBasis occ_basis_set;
	setup_basis(mol, settings.basis_set_name, occ_basis_set);
	occ::qm::HartreeFock hf(occ_basis_set);
	create_prims(ao_data_shells, occ_basis_set);
	eval_I_anom_disp(ao_data_shells, read);
	if (safe) {
		std::ofstream out("I_tensor", std::ios::binary);
		if (!out)
			throw std::runtime_error("Cannot open file for writing");
		int nr_safe = cryst.nr_small;
		int nmo_safe = cryst.nmo;
		int num_elements_safe = (cryst.nmo * (cryst.nmo + 1)) / 2;
		int total_size_safe = static_cast<int>(I.size());
		out.write(reinterpret_cast<const char*>(&nr_safe), sizeof(nr_safe));
		out.write(reinterpret_cast<const char*>(&nmo_safe), sizeof(nmo_safe));
		out.write(reinterpret_cast<const char*>(&num_elements_safe), sizeof(num_elements_safe));
		out.write(reinterpret_cast<const char*>(&total_size_safe), sizeof(total_size_safe));
		out.write(reinterpret_cast<const char*>(I.data()),
			total_size_safe * sizeof(cdouble));
	}
	return hf;
	// closing function
}

// Needs rework
//void XCW::calc_F_calc_fast() {
//	eval_phase();
//	//eval_DW();
//	eval_translation_phase();
//	eval_anom_disp();
//
//	cvec2 atomic_scattering_factors(ncen, cvec(nr, 0));
//	F_calc.resize(nr_small, 0);
//
//	// Calculate atomic scattering factors for each atom and symmetry generated reflexes
//	GridConfiguration config;
//	config.accuracy = opt->accuracy;
//	config.partition_type = opt->partition_type;
//	config.pbc = opt->pbc;
//	config.no_density_eval = false;
//	config.debug = opt->debug;
//	config.all_charges = opt->all_charges;
//	GridManager grid_manager(config);
//	WFN temp = wave;
//	vec2 d1, d2, d3, dens;
//	std::vector<_time_point> time_points({ get_time() });
//	_time_point end;
//	temp.delete_unoccupied_MOs();
//	grid_manager.setup3DGridsForMolecule(temp, asym_atom_list, needs_grid, unit_cell);
//	grid_manager.calculateNonSphericalDensities(temp, unit_cell);
//	svec time_descriptions;
//	grid_manager.addTimingInfoToVecs(time_points, time_descriptions);
//	PartitionResults results = grid_manager.calculatePartitionedCharges(temp, unit_cell);
//	grid_manager.printChargeTable(labels, temp, asym_atom_list, std::cout, results);
//	time_points.push_back(get_time());
//	time_descriptions.push_back("calculate charges");
//	grid_manager.getDensityVectors(temp, asym_atom_list, d1, d2, d3, dens);
//	const int points = grid_manager.getTotalGridPoints();
//	calc_SF(points, k_pt, d1, d2, d3, dens, atomic_scattering_factors, std::cout, time_points.front(), end, opt->debug, true, true);
//	// Calculate F_calc
//#pragma omp parallel for
//	for (int r = 0; r < nr_small; r++) {
//		const ivec& lookup = generate_asym_lookup(r);
//		for (int at = 0; at < ncen; at++) {
//			for (int r_asym = 0; r_asym < lookup.size(); r_asym++) {
//				F_calc[r] += atomic_scattering_factors[at][lookup[r_asym]] * DW_fact[at][lookup[r_asym]] * phase_fact[at][lookup[r_asym]] * translation_phase[r][r_asym] * asym_atoms[at].asym_fact;
//			}
//		}
//		// Add anomalous dispersion correction
//		F_calc[r] += anom_corr[r];
//		//std::cout << std::fixed << std::setprecision(5) << std::pow(std::abs(F_calc[r]), 2) << std::endl;
//	}
//	//dump F_calc values as binary file called F_calc
//
//	std::ofstream fout("F_calc.bin", std::ios::out | std::ios::binary);
//	//First byte is the number of bytes per double, the next one is the size of a compelx double, to understand how to read the data.
//	//After that an int64 (8 byte) of the number of F.calc values to be expected after that.
//	//Finally, the dump of all F_calc values as cdouble (A,B)
//	char size = sizeof(double);
//	fout.write(reinterpret_cast<const char*>(&size), sizeof(size));
//	size = sizeof(cdouble);
//	fout.write(reinterpret_cast<const char*>(&size), sizeof(size));
//	size_t vec_size = F_calc.size();
//	fout.write(reinterpret_cast<const char*>(&vec_size), sizeof(size_t));
//	fout.write(reinterpret_cast<const char*>(F_calc.data()), vec_size * sizeof(cdouble));
//	fout.close();
//
//}

void XCW::run_XCW_fitting() {

	bool read = true;
	bool safe = false;
	occ::qm::HartreeFock hf = setup_XCW_procedure(read, safe);
	occ::qm::SCF scf(hf, settings.hf_type);
	occ::qm::Wavefunction last_wfn;
	bool has_guess = false;

	std::cout << "More detailed output in XCW.log file..." << std::endl;
	std::cout << "____________________________________________________________________________\n";
	std::cout << " Lambda\t\tChi^2\tGooF\tTotal Energy\tPerturbation\tTarget quantity \n";
	std::cout << "								(Eh)		   (a. u.)			(a. u.)\n";
	std::cout << "____________________________________________________________________________\n";

	// Runs the lambda steps for XCW fitting
	for (int step = 0; step < settings.num_xcw_steps; step++) {

		occ::qm::SCF scf(hf, settings.hf_type);
		double alpha = settings.alpha;
		const double lambda = step * settings.xcw_step_size;
		scf.set_charge_multiplicity(settings.charge, settings.multiplicity);
		scf.maxiter = settings.max_scf_iterations;
		scf.convergence_settings.level_shift = settings.level_shift;
		scf.convergence_settings.level_shift_threshold = 0;
		scf.update_occupied_orbital_count();
		scf.convergence_accelerator.set_strategy(scf.convergence_settings.diis_strategy);
		scf.convergence_accelerator.set_switch_threshold(scf.convergence_settings.diis_switch_threshold);
 		do_SCF(lambda, alpha, scf, last_wfn, has_guess);
		last_wfn = scf.wavefunction();
	}

    std::cout << "Finished XCW fitting procedure." << std::endl;
    exit(0);

}
