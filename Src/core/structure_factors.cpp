#include "pch.h"
#include "structure_factors.h"
#if defined(NOSPHERA2_USE_GPU) || defined(NOSPHERA2_USE_METAL)
#include "itensor_gpu.h"
#endif
#include "wfn_class.h"
#include "scattering_factors.h"
#include "GridManager.h"
#include "nos_math.h"
#include "basis_set.h"
#include <mutex>
#include <limits>


structure_factors::structure_factors(options& opt_in) {
	opt = &opt_in;

	// Load the settings from the options object
	model_data.n_params = opt_in.xcw_settings.n_params;
	extinction_settings.extinction_model = opt_in.xcw_settings.extinction_model;
	I_tens.read_tensor = opt_in.xcw_settings.read_tensor;
	I_tens.i_tensor_file_path = opt_in.xcw_settings.i_tensor_file_path;
	I_tens.i_tensor_save_path = opt_in.xcw_settings.i_tensor_save_path;
	I_tens.tensor_single = opt_in.xcw_settings.i_tensor_single;
	I_tens.tensor_double = opt_in.xcw_settings.i_tensor_double;
	I_tens.basis_set_name = opt_in.xcw_settings.basis_set_name;
	extinction_settings.aniso = opt_in.xcw_settings.extinction_aniso;
	extinction_settings.refine = opt_in.xcw_settings.extinction_refine;
	extinction_settings.start = opt_in.xcw_settings.extinction_start;
	quality_criteria.refine_against = opt_in.xcw_settings.refine_against;
	quality_criteria.goof_type = opt_in.xcw_settings.XWR_type;
	I_tens.i_tensor_max_mb = opt_in.xcw_settings.i_tensor_max_mb;


	// Setup file paths
	std::filesystem::path hkl_filename = opt_in.hkl;
	std::filesystem::path cif = opt_in.cif;
	std::ifstream cif_input(cif.c_str(), std::ios::in);

	// Read the hkl and cif file and generate the unit cell
	unit_cell = cell(cif, std::cout, opt_in.debug, opt_in.do_XCW);
	scatter_data.hkl_enlarged = read_hkl_full(hkl_filename, scatter_data.hkl, opt_in.twin_law, unit_cell, std::cout, scatter_data, opt_in.debug);
	double wavelength_ = read_CIF(cif_input, unit_cell, model_data.ncen, asym_atoms, ADPs, opt_in.debug);
	// Directly convert into reciprocal space so rotation for grown structures can be done directly
	U_cif2U_star();
	model_data.wavelength = opt_in.xcw_settings.wavelength > 0.0 ? opt_in.xcw_settings.wavelength : wavelength_;
	err_checkf(model_data.ncen > 0, "No atoms were read from " + cif.string() + "! Is there an _atom_site loop with labels, type symbols and fractional coordinates?", std::cout);

	{
		ivec3 symmetry_linking_list;
		ivec applied_symmetry;
		// Handle grown structures
		std::vector<asym_atom> xyz_atoms;
		if (opt_in.xcw_settings.grown) {
			err_checkf(!opt_in.xyz_file.empty(), "Grown structures require an xyz file with the grown structure. Please provide one with the `-xyz` option.", std::cout);
			// Read xyz file and grow the asymmetric unit
			const std::filesystem::path xyz_path = opt_in.xyz_file;
			WFN dummy_wave(xyz_path, opt_in.debug);
			dummy_wave.read_xyz(xyz_path, std::cout, opt_in.debug);
			xyz_atoms = dummy_wave.extract_xyz("bohr");
			unit_cell.grow_asym_atoms(asym_atoms, xyz_atoms);
		}
		/*Generate the symmetry linking list
		The linking list is ordered like this: Asymmetric atom, list with all atoms, then index of symmetry operation that generated it
		"diagonal elements" have to have size equivalent to multiplicity, otherwise something broke */
		unit_cell.eval_symm(asym_atoms, model_data.ncen, symmetry_linking_list);
		model_data.ncen = asym_atoms.size();
		if (opt_in.xcw_settings.grown) {
			// Copy U_iso, the dispersion and the ADPs of each grown atom's parent, the ADPs rotated onto the image
			unit_cell.grow_ADPs(asym_atoms, symmetry_linking_list, ADPs);
			// Project grown structure into its symmetry subgroup and update everything accordingly
			applied_symmetry = unit_cell.apply_grown(scatter_data.hkl, scatter_data.hkl_enlarged, asym_atoms, symmetry_linking_list);
		}
		// Set the symmetry factors for each atom
		unit_cell.set_symmetry_factors(asym_atoms, symmetry_linking_list, applied_symmetry);
		if (std::getenv("NOSPHERA2_DEBUG_ASYMFACT")) { // Flawfinder: ignore
			std::cerr << "applied_symmetry (deleted):";
			for (int s : applied_symmetry) std::cerr << " " << s;
			std::cerr << std::endl << "surviving sym ops: " << unit_cell.get_trans()[0].size() << std::endl;
			for (size_t i = 0; i < asym_atoms.size(); i++)
				std::cerr << i << " grown=" << asym_atoms[i].grown << " sym_op=" << asym_atoms[i].sym_op
				<< " asym_fact=" << asym_atoms[i].asym_fact << std::endl;
			std::cerr << "hkl_enlarged size: " << scatter_data.hkl_enlarged.size() << std::endl;
		}
	}

	// Convert the ADPs from reciprocal space to Cartesian coordinates
	U_star2U_cart();

	// Generate k_pts and set the number of reflections
	make_k_pts(model_data.nr_enlarged != 0 && scatter_data.hkl.size() == 0, opt_in.save_k_pts, unit_cell, scatter_data.hkl_enlarged, k_pt, std::cout, opt_in.debug);
	model_data.nr_enlarged = scatter_data.hkl_enlarged.size();
	model_data.nr = scatter_data.hkl.size();

	// Find all reflections that are above the I/sigma cutoff
	const double i_sigma_cutoff = opt_in.xcw_settings.i_sigma_cutoff;
	scatter_data.hkl_mask.resize(model_data.nr, 0);
	model_data.nr_fit = 0;
	for (int r = 0; r < model_data.nr; r++) {
		const double I_over_sigma = (scatter_data.F_obs[r] < 0 ? -scatter_data.F_obs2[r] : scatter_data.F_obs2[r]) / scatter_data.sigma_obs2[r];
		scatter_data.hkl_mask[r] = scatter_data.sigma_obs2[r] > 0 && I_over_sigma >= i_sigma_cutoff;
		model_data.nr_fit += scatter_data.hkl_mask[r];
	}

	// Extinction correction
	setup_extinction(cif);
	err_checkf(model_data.nr_fit > n_params(), "Fewer reflections above the I/sigma cutoff than parameters", std::cout);
	std::cout << "XCW: I/sigma(I) >= " << i_sigma_cutoff << " (F/sigma(F) >= " << 2 * i_sigma_cutoff << "): " << model_data.nr_fit << " of " << model_data.nr << " reflections in the fit; R1 and Criterion are over these, R1(all) and Crit(all) over all" << std::endl;

	// Set F_calc sizes
	scatter_data.F_calc.resize(model_data.nr, 0);
	scatter_data.anom_correction.resize(model_data.nr, 0);

	// Initialize DW factors and phase factors, so that F_calc can be calculated without DW factors or phase factors if they are not requested
	DW_facts.resize(model_data.ncen, cvec(model_data.nr_enlarged, 1));
	phase_facts.resize(model_data.ncen, cvec(model_data.nr_enlarged, 1));
	translation_phase_facts.resize(model_data.nr, cvec(unit_cell.get_trans()[0].size(), 1));

}

void structure_factors::setup_extinction(const std::filesystem::path& cif) {
	if (extinction_settings.extinction_model == extinction::model::none) return;
	const double lambda = model_data.wavelength;
	err_checkf(lambda > 0.0, "Extinction needs a wavelength: put `wavelength <lambda>` in the XCW "
		"settings file, or _diffrn_radiation_wavelength in " + cif.string(), std::cout);
	ensure_hkl_ordered();
	const size_t np = extinction_settings.aniso ? 6 : 1;
	//the anisotropic tensor starts isotropic, where x(h) is the start value for every h
	ext_p_.assign(np, 0.0);
	for (size_t p = 0; p < (extinction_settings.aniso ? 3u : 1u); p++) ext_p_[p] = extinction_settings.start;
	ext_c_.resize(model_data.nr);
	ext_cos2t_.resize(model_data.nr);
	if (extinction_settings.aniso) ext_a_.resize(static_cast<size_t>(model_data.nr) * np);
	for (int r = 0; r < model_data.nr; r++) {
		const double stl = unit_cell.get_stl_of_hkl(hkl_ordered_[r]);
		ext_c_[r] = extinction::geometry_constant(lambda, stl);
		ext_cos2t_[r] = extinction::cos_2theta(lambda, stl);
		if (!extinction_settings.aniso) continue;
		//the scattering vector in Cartesian, |h| = 1/d: rcm's rows are the Cartesian
		//components, its columns the reciprocal basis vectors
		std::array<double, 3> h_unit{ 0.0, 0.0, 0.0 };
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++) h_unit[i] += unit_cell.get_rcm_angs(i, j) * hkl_ordered_[r][j];
		const double norm = std::sqrt(h_unit[0] * h_unit[0] + h_unit[1] * h_unit[1] + h_unit[2] * h_unit[2]);
		if (norm > 0.0) for (double& v : h_unit) v /= norm;
		std::array<double, 6> a{};
		extinction::aniso_coefficients(h_unit, a);
		for (size_t p = 0; p < np; p++) ext_a_[static_cast<size_t>(r) * np + p] = a[p];
	}
	ext_y_.assign(model_data.nr, 1.0);
	ext_sqrt_y_.assign(model_data.nr, 1.0);
	ext_g_.assign(model_data.nr, 1.0);
	ext_m_.assign(model_data.nr, 1.0);
	ext_dyc_.assign(model_data.nr, 0.0);
	std::ostringstream banner;
	banner << "XCW extinction: " << extinction::name(extinction_settings.extinction_model)
		<< (extinction_settings.aniso ? ", anisotropic (azimuth-averaged, 6 parameters)" : ", isotropic (1 parameter)")
		<< (extinction_settings.refine ? ", refined with the scale" : ", held fixed")
		<< ", lambda = " << lambda << " A, start value " << extinction_settings.start;
	std::cout << banner.str() << std::endl;
}

void structure_factors::ensure_hkl_ordered() {
	if (!hkl_ordered_.empty() || scatter_data.hkl.empty()) {
		return;
	}
	hkl_ordered_.reserve(scatter_data.hkl.size());
	for (const i3& h : scatter_data.hkl) {
		hkl_ordered_.push_back(h);
	}
}

void structure_factors::U_cif2U_star() {
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

	for (int a = 0; a < model_data.ncen; a++) {
		vec2 ADPs_ = ADPs[a];
		if (ADPs_.size() > 0 && ADPs_[0].size() > 0) {
			for (int i = 0; i < 6; i++) {
				ADPs_[0][i] *= transform[i];
			}
			ADPs[a] = ADPs_;
		}
	}
}

void structure_factors::U_star2U_cart() {
	const double scale = constants::bohr2ang(1);
	vec2 cart_matrix(3, vec(3));
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			cart_matrix[i][j] = unit_cell.get_cm(i, j) * scale;
		}
	}
	for (int a = 0; a < model_data.ncen; a++) {
		vec2 ADPs_ = ADPs[a];
		cell::transform_ADPs(ADPs_, cart_matrix);
		ADPs[a] = ADPs_;
	}
}

template <int N>
static const std::vector<structure_factors::gc_term>& gc_terms() {
	static const std::vector<structure_factors::gc_term> terms = [] {
		const int fact[5] = { 1, 1, 2, 6, 24 };
		int total = 1;
		for (int i = 0; i < N; i++) total *= 3;
		std::vector<structure_factors::gc_term> t;
		for (int n = 0; n < total; n++) {
			int idx[N], m = n;
			for (int i = N - 1; i >= 0; i--) { idx[i] = m % 3; m /= 3; }
			structure_factors::gc_term g = { { 0, 0, 0 }, 0.0 };
			bool sorted = true;
			for (int i = 0; i < N; i++) {
				g.e[idx[i]]++;
				if (i > 0 && idx[i] < idx[i - 1]) sorted = false;
			}
			if (!sorted) continue;
			g.mult = static_cast<double>(fact[N]) / (fact[g.e[0]] * fact[g.e[1]] * fact[g.e[2]]);
			t.push_back(g);
		}
		return t;
	}();
	return terms;
}

//Sum of mult * c * x^e0 y^e1 z^e2 over the packed tensor c, with p[d][axis] the d-th power of the point
template <int N>
static double gc_sum(const vec& c, const double (*p)[3]) {
	const std::vector<structure_factors::gc_term>& t = gc_terms<N>();
	double sum = 0.0;
	for (int n = 0; n < static_cast<int>(t.size()); n++)
		sum += t[n].mult * c[n] * p[t[n].e[0]][0] * p[t[n].e[1]][1] * p[t[n].e[2]][2];
	return sum;
}

//0 isotropic, 1 U, 2 U and C, 3 U, C and D; an atom without ADPs gets empty ones
static int adp_level(vec2& adp) {
	if (adp.size() != 3) {
		adp.resize(3);
		return 0;
	}
	if (!adp[2].empty()) return 3;
	if (!adp[1].empty()) return 2;
	if (!adp[0].empty()) return 1;
	return 0;
}

void structure_factors::set_DW() {
	//Converts angstrom to bohr OR MORE IMPORTANTLY reciprocal bohr to reciprocal angstrom
	const double angstrom2bohr = constants::ang2bohr(1);
	ivec level;
	level.reserve(model_data.ncen);
	for (int a = 0; a < model_data.ncen; a++)
		level.emplace_back(adp_level(ADPs[a]));
	vec2 q(model_data.nr_enlarged, vec(3));
	for (int h = 0; h < model_data.nr_enlarged; h++) {
		q[h][0] = k_pt[0][h];
		q[h][1] = k_pt[1][h];
		q[h][2] = k_pt[2][h];
	}
	std::transform(q.begin(), q.end(), q.begin(), [angstrom2bohr](std::vector<double>& vec) {
		std::transform(vec.begin(), vec.end(), vec.begin(), [angstrom2bohr](double x) { return x * angstrom2bohr; });
		return vec; });
	for (int a = 0; a < model_data.ncen; a++) {
		vec2 ADPs_ = ADPs[a];
		if (level[a] == 0) {
			// Isotropic
			double U = asym_atoms[a].U_iso, temp;
			for (int r = 0; r < model_data.nr_enlarged; r++) {
				temp = -0.5 * U * (q[r][0] * q[r][0] + q[r][1] * q[r][1] + q[r][2] * q[r][2]);
				DW_facts[a][r] = std::exp(temp);
			}
			continue;
		}
		// Anisotropic U_ij, with the C_ijk (level 2) and D_ijkl (level 3) terms of the Gram-Charlier expansion
		vec2 Uij = { { ADPs_[0][0], ADPs_[0][3], ADPs_[0][4] },
					 { ADPs_[0][3], ADPs_[0][1], ADPs_[0][5] },
					 { ADPs_[0][4], ADPs_[0][5], ADPs_[0][2] } };
		for (int h = 0; h < model_data.nr_enlarged; h++) {
			vec q_ = { q[h][0], q[h][1], q[h][2] };
			const double temp1 = -0.5 * dot_BLAS(dot(Uij, q_, true), q_, false);
			// p stores the powers of q up to 4
			double p[5][3] = { { 1.0, 1.0, 1.0 } };
			for (int d = 1; d < 5; d++)
				for (int i = 0; i < 3; i++) {
					p[d][i] = p[d - 1][i] * q_[i];
				}
			const double c3 = level[a] >= 2 ? -1.0 / 6.0 * gc_sum<3>(ADPs_[1], p) : 0.0;
			const double d4 = level[a] >= 3 ? 1.0 / 24.0 * gc_sum<4>(ADPs_[2], p) : 0.0;
			DW_facts[a][h] = std::exp(temp1) * cdouble(1 + d4, c3);
		}
	}
	DW_set_ = true;
}

void structure_factors::set_phases() {
	cdouble exponent;
	for (int at = 0; at < model_data.ncen; at++) {
		vec pos_cart = { asym_atoms[at].pos[0], asym_atoms[at].pos[1], asym_atoms[at].pos[2] };
		for (int r = 0; r < model_data.nr_enlarged; r++) {
			vec q = { k_pt[0][r], k_pt[1][r], k_pt[2][r] };
			exponent = cdouble(0, dot_BLAS(q, pos_cart, false));
			phase_facts[at][r] = std::exp(exponent);
		}
	}
	const double angstrom2bohr = constants::ang2bohr(1);
	const double bohr2angstrom = constants::bohr2ang(1);
	vec2 trans = unit_cell.get_trans();
	vec2 cm = { { unit_cell.get_cm(0,0), unit_cell.get_cm(0,1), unit_cell.get_cm(0,2)},
								  { unit_cell.get_cm(1,0), unit_cell.get_cm(1,1), unit_cell.get_cm(1,2)},
								  { unit_cell.get_cm(2,0), unit_cell.get_cm(2,1), unit_cell.get_cm(2,2)} };
	std::transform(cm.begin(), cm.end(), cm.begin(), [bohr2angstrom](std::vector<double>& vec) {
		std::transform(vec.begin(), vec.end(), vec.begin(), [bohr2angstrom](double x) { return x * bohr2angstrom; });
		return vec; });
	for (int r = 0; r < model_data.nr; r++) {
		ivec asym_list = generate_asym_lookup(r);
		vec q_temp = { k_pt[0][asym_list[0]], k_pt[1][asym_list[0]], k_pt[2][asym_list[0]] };
		std::transform(q_temp.begin(), q_temp.end(), q_temp.begin(), [angstrom2bohr](double x) { return x * angstrom2bohr; });
		for (int t = 0; t < trans[0].size(); t++) {
			vec trans_temp = { trans[0][t], trans[1][t], trans[2][t] };
			trans_temp = dot(cm, trans_temp, true);
			cdouble exponent(0, dot_BLAS(q_temp, trans_temp, false));
			translation_phase_facts[r][t] = std::exp(exponent);
		}
	}
	phases_set_ = true;
	// closing function
}

void structure_factors::set_anom() {
	int r, at, r_asym;
	ivec2 asym_lookup(model_data.nr);
	for (r = 0; r < model_data.nr; r++) {
		asym_lookup[r] = generate_asym_lookup(r);
	}
	for (r = 0; r < model_data.nr; r++) {
		const ivec& lookup = asym_lookup[r];
		cdouble sum = 0;
		for (at = 0; at < model_data.ncen; at++) {
			cdouble temp1 = 0;
			for (r_asym = 0; r_asym < lookup.size(); r_asym++) {
				temp1 += phase_facts[at][lookup[r_asym]] * DW_facts[at][lookup[r_asym]] * translation_phase_facts[r][r_asym];
			}
			sum += temp1 * asym_atoms[at].asym_fact * asym_atoms[at].anom;
		}
		scatter_data.anom_correction[r] = sum;
	}
	anom_set_ = true;
}

//y_r and the chain factors, from the current F_calc and the current coefficients
void structure_factors::update_extinction() {
	if (ext_p_.empty()) return;
	const extinction::model m = extinction_settings.extinction_model;
#pragma omp parallel for
	for (int r = 0; r < model_data.nr; r++) {
		const double u = std::norm(scatter_data.F_calc[r]);
		const double t = ext_c_[r] * ext_x(r) * u;
		double dydt = 0.0;
		const double y = extinction::correction(m, ext_cos2t_[r], t, &dydt);
		ext_y_[r] = y;
		ext_sqrt_y_[r] = std::sqrt(y);
		//I = y u, so dI/du = y + t dy/dt, and the amplitude's slope is that over sqrt(y)
		ext_g_[r] = y + t * dydt;
		ext_m_[r] = ext_sqrt_y_[r] > 0.0 ? ext_g_[r] / ext_sqrt_y_[r] : 1.0;
		ext_dyc_[r] = dydt * ext_c_[r] * u;
	}
}

//One Gauss-Newton step on the coefficients against the criterion's own weighted residual,
//with the scale held where solve_scale put it. Accepted only if it lowers that residual and
//leaves x_r >= 0 everywhere, since a negative coefficient is not extinction.
bool structure_factors::refine_extinction_step() {
	const Eigen::Index np = static_cast<Eigen::Index>(ext_p_.size());
	const bool against_F2 = quality_criteria.refine_against == 2, weighted = quality_criteria.goof_type == 2;
	const double k = scatter_data.scale, s = k * k;
	auto residual_sum = [&]() {
		update_extinction();
		double sum = 0.0;
		for (int r = 0; r < model_data.nr; r++) {
			if (!scatter_data.hkl_mask[r]) continue;
			if (ext_x(r) < 0.0) return std::numeric_limits<double>::infinity();
			const double w = weighted ? inv_H2_[r] : 1.0;
			const double d = against_F2
				? (s * ext_y_[r] * std::norm(scatter_data.F_calc[r]) - scatter_data.F_obs2[r]) / scatter_data.sigma_obs2[r]
				: (k * ext_sqrt_y_[r] * std::abs(scatter_data.F_calc[r]) - scatter_data.abs_F_obs[r]) / scatter_data.sigma_obs[r];
			sum += w * d * d;
		}
		return sum;
		};

	occ::Mat JtJ = occ::Mat::Zero(np, np);
	occ::Vec JtR = occ::Vec::Zero(np), J(np);
	double chi2 = 0.0;
	for (int r = 0; r < model_data.nr; r++) {
		if (!scatter_data.hkl_mask[r]) continue;
		const double w = weighted ? inv_H2_[r] : 1.0, Fm = std::abs(scatter_data.F_calc[r]);
		double resid = 0.0, dmodel = 0.0;   //d(model)/dP_p = dmodel * a_{r,p}
		if (against_F2) {
			resid = (s * ext_y_[r] * Fm * Fm - scatter_data.F_obs2[r]) / scatter_data.sigma_obs2[r];
			dmodel = s * Fm * Fm * ext_dyc_[r] / scatter_data.sigma_obs2[r];
		}
		else {
			if (ext_sqrt_y_[r] <= 0.0) continue;
			resid = (k * ext_sqrt_y_[r] * Fm - scatter_data.abs_F_obs[r]) / scatter_data.sigma_obs[r];
			dmodel = k * Fm * ext_dyc_[r] / (2.0 * ext_sqrt_y_[r] * scatter_data.sigma_obs[r]);
		}
		chi2 += w * resid * resid;
		if (dmodel == 0.0) continue;
		for (Eigen::Index p = 0; p < np; p++) J(p) = dmodel * ext_a(r, static_cast<size_t>(p));
		JtJ += w * J * J.transpose();
		JtR += w * resid * J;
	}
	//Levenberg damping, so a direction the data barely sees does not throw the step
	for (Eigen::Index p = 0; p < np; p++) JtJ(p, p) *= 1.001;
	const occ::Vec step = JtJ.ldlt().solve(-JtR);
	if (!step.allFinite() || step.norm() == 0.0) return false;
	const vec start = ext_p_;
	for (int half = 0; half < 8; half++) {
		const double f = std::pow(0.5, half);
		for (Eigen::Index p = 0; p < np; p++) ext_p_[p] = start[p] + f * step(p);
		if (residual_sum() < chi2) return true;
	}
	ext_p_ = start;
	update_extinction();
	return false;
}

std::string structure_factors::extinction_report() const {
	if (ext_p_.empty()) return "";
	std::ostringstream out;
	out << "extinction(" << extinction::name(extinction_settings.extinction_model)
		<< (extinction_settings.aniso ? ", aniso)" : ")") << std::scientific << std::setprecision(4);
	for (const double p : ext_p_) out << " " << p;
	return out.str();
}

//The scale that minimises the criterion the SCF descends: k over the fit set from
//Sum w (k|Fc| - |Fo|)^2 / sigma^2, or k^2 from Sum w (k^2|Fc|^2 - Fo^2)^2 / sigma(I)^2 against
//F^2, with w the 1/|H|^2 weights of XWR_type 2. perturbation takes the scale as given, so the
//scale had better be stationary for the criterion, or the SCF descends a different
//functional than the one it prints: an unweighted fit of k put Fe(phen)2(SCN)2 (chi^2 12.96
//vs 12.30 at the weighted k, dchi^2/dk 2e3) 4 mEh above the previous step's orbitals at
//every converged lambda >= 0.03.
void structure_factors::eval_scale() {
	ensure_inv_H2_weights();
	update_extinction();
	solve_scale();
	if (ext_p_.empty() || !extinction_settings.refine) return;
	//the scale and the extinction coefficients are coupled through the same residual, so
	//alternate: a Gauss-Newton step on the coefficients, then the closed-form scale again
	for (int it = 0; it < 5; it++) {
		if (!refine_extinction_step()) break;
		update_extinction();
		solve_scale();
	}
}

void structure_factors::solve_scale() {
	const bool against_F2 = quality_criteria.refine_against == 2, weighted = quality_criteria.goof_type == 2;
	const int chunk = 128, nchunk = (model_data.nr + chunk - 1) / chunk;
	vec numerators(nchunk), denominators(nchunk);
#pragma omp parallel for schedule(static)
	for (int c = 0; c < nchunk; c++) {
		const int first = c * chunk, last = std::min(first + chunk, model_data.nr);
		for (int i = first; i < last; i++) {
			if (!scatter_data.hkl_mask[i]) continue;
			const double w = weighted ? inv_H2_[i] : 1.0;
			const double calc = ext_sqrt_y(i) * std::abs(scatter_data.F_calc[i]);
			if (against_F2) {
				const double calc2 = calc * calc, wi = w / (scatter_data.sigma_obs2[i] * scatter_data.sigma_obs2[i]);
				numerators[c] += wi * calc2 * scatter_data.F_obs2[i];
				denominators[c] += wi * calc2 * calc2;
			}
			else {
				const double wi = w / (scatter_data.sigma_obs[i] * scatter_data.sigma_obs[i]);
				numerators[c] += wi * calc * scatter_data.F_obs[i];
				denominators[c] += wi * calc * calc;
			}
		}
	}
	double numerator = 0.0, denominator = 0.0;
	for (int c = 0; c < nchunk; c++) {
		numerator += numerators[c];
		denominator += denominators[c];
	}
	const double ratio = (denominator != 0.0) ? numerator / denominator : 1.0;
	scatter_data.scale = against_F2 ? std::sqrt(std::max(ratio, 0.0)) : ratio;
}

void structure_factors::calc_criteria() {
	ensure_inv_H2_weights();
	//index 0: the fit set, 1: all reflections
	const double prefactor[2] = { 1.0 / static_cast<double>(model_data.nr_fit - n_params()), 1.0 / static_cast<double>(model_data.nr - n_params()) };
	const int chunk = 128, nchunk = (model_data.nr + chunk - 1) / chunk;
	vec2 goof1(2, vec(nchunk)), goof2(2, vec(nchunk)), wgoof1(2, vec(nchunk)), wgoof2(2, vec(nchunk)), r1_num(2, vec(nchunk)), r1_den(2, vec(nchunk));
	const double scale = scatter_data.scale;
	const cdouble* F_calc_0 = scatter_data.F_calc.data();
#pragma omp parallel for schedule(static)
	for (int c = 0; c < nchunk; c++) {
		const int first = c * chunk, last = std::min(first + chunk, model_data.nr);
		for (int i = first; i < last; i++) {
			const double scaled_F_calc = scale * ext_sqrt_y(i) * std::abs(F_calc_0[i]);
			const double scaled_difference = scaled_F_calc - scatter_data.F_obs[i];
			const double diff2 = (scaled_F_calc * scaled_F_calc) - scatter_data.F_obs2[i];
			const double weighted_diff1 = scaled_difference / scatter_data.sigma_obs[i];
			const double weighted_diff2 = diff2 / scatter_data.sigma_obs2[i];
			const double weighted_diff1_sq = weighted_diff1 * weighted_diff1;
			const double weighted_diff2_sq = weighted_diff2 * weighted_diff2;
			const double w = quality_criteria.goof_type == 2 ? inv_H2_[i] : 0.0;
			for (int set = scatter_data.hkl_mask[i] ? 0 : 1; set < 2; set++) {
				goof1[set][c] += weighted_diff1_sq;
				goof2[set][c] += weighted_diff2_sq;
				wgoof1[set][c] += weighted_diff1_sq * w;
				wgoof2[set][c] += weighted_diff2_sq * w;
				r1_num[set][c] += std::abs(scaled_F_calc - scatter_data.abs_F_obs[i]);
				r1_den[set][c] += scatter_data.abs_F_obs[i];
			}
		}
	}
	double sum[6][2] = {};
	for (int set = 0; set < 2; set++) {
		for (int c = 0; c < nchunk; c++) {
			sum[0][set] += goof1[set][c];
			sum[1][set] += goof2[set][c];
			sum[2][set] += wgoof1[set][c];
			sum[3][set] += wgoof2[set][c];
			sum[4][set] += r1_num[set][c];
			sum[5][set] += r1_den[set][c];
		}
	}
	quality_criteria.R1 = sum[5][0] > 0.0 ? sum[4][0] / sum[5][0] : 0.0;
	quality_criteria.R1_all = sum[5][1] > 0.0 ? sum[4][1] / sum[5][1] : 0.0;
	quality_criteria.GooF1 = std::sqrt(prefactor[0] * sum[0][0]);
	quality_criteria.GooF2 = std::sqrt(prefactor[0] * sum[1][0]);
	quality_criteria.weighted_GooF1 = std::sqrt(prefactor[0] * sum[2][0]);
	quality_criteria.weighted_GooF2 = std::sqrt(prefactor[0] * sum[3][0]);
	quality_criteria.GooF1_all = std::sqrt(prefactor[1] * sum[0][1]);
	quality_criteria.GooF2_all = std::sqrt(prefactor[1] * sum[1][1]);
	quality_criteria.weighted_GooF1_all = std::sqrt(prefactor[1] * sum[2][1]);
	quality_criteria.weighted_GooF2_all = std::sqrt(prefactor[1] * sum[3][1]);
}

double structure_factors::criterion(const bool all) const {
	const bool weighted = quality_criteria.goof_type == 2, against_F2 = quality_criteria.refine_against == 2;
	if (all) return weighted ? (against_F2 ? quality_criteria.weighted_GooF2_all : quality_criteria.weighted_GooF1_all) : (against_F2 ? quality_criteria.GooF2_all : quality_criteria.GooF1_all);
	return weighted ? (against_F2 ? quality_criteria.weighted_GooF2 : quality_criteria.weighted_GooF1) : (against_F2 ? quality_criteria.GooF2 : quality_criteria.GooF1);
}

//Per-reflection 1/|H|^2 weights for the residual self-energy criterion,
//U_res ~ Sum_h |dF_h|^2/|H_h|^2, with |H| = 1/d = 2*sin(theta)/lambda.
//(0,0,0) is already excluded from hkl at read time; depends only on geometry.
void structure_factors::ensure_inv_H2_weights() {
	if (quality_criteria.goof_type == 1 || !inv_H2_.empty()) {
		return;
	}
	ensure_hkl_ordered();
	inv_H2_.resize(model_data.nr);
#pragma omp parallel for
	for (int r = 0; r < model_data.nr; r++) {
		const double stl = unit_cell.get_stl_of_hkl(hkl_ordered_[r]);
		const double H2 = 4.0 * stl * stl;
		inv_H2_[r] = (H2 > 0.0) ? 1.0 / H2 : 0.0;
	}
}

void structure_factors::create_prims(std::vector<ao_data>& ao_data_shells, const occ::qm::AOBasis& occ_basis_set) {
	for (int atm = 0; atm < model_data.ncen; atm++) {
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
	model_data.nmo = ao_data_shells.size();
}

ivec structure_factors::generate_asym_lookup(const int r) {
	ivec asym_list;
	auto it = scatter_data.hkl.begin();
	std::advance(it, r);
	ivec3 rots = unit_cell.get_sym();
	i3 tempv;
	const i3& hkl_temp = *it;
	for (int s = 0; s < rots[0][0].size(); s++) {
		tempv = { 0, 0, 0 };
		for (int h = 0; h < 3; h++) {
			for (int j = 0; j < 3; j++) {
				tempv[j] += hkl_temp[h] * rots[j][h][s];
			}
		}
		int idx_ = 0;
		auto idx = scatter_data.hkl_enlarged.find(tempv);
		if (idx != scatter_data.hkl_enlarged.end()) {
			idx_ = std::distance(scatter_data.hkl_enlarged.begin(), idx);
		}
		asym_list.push_back(idx_);
	}
	return asym_list;
	// closing function
}

size_t structure_factors::tri_index(int mu, int nu) const noexcept {
	return mu * model_data.nmo - (mu * (mu - 1)) / 2 + (nu - mu);
}

structure_factors::I_tensor& structure_factors::eval_I(const occ::gto::AOBasis& aobasis) {
	if (!DW_set_) std::cout << "Debye-Waller factors are not computed and set to 1. Continuing." << std::endl;
	if (!phases_set_) std::cout << "Phase factors are not computed and set to 1. Continuing." << std::endl;
	if (!anom_set_) std::cout << "Anomalous dispersion corrections are not computed and set to 0. Continuing." << std::endl;

	// Creates the primitive data for each shell, needed to compute the I tensor
	std::vector<structure_factors::ao_data> ao_data_shells;
	create_prims(ao_data_shells, aobasis);

	bool single_on_disk = false;
	// Path for reading the I tensor from disk
	if (I_tens.read_tensor && !I_tens.i_tensor_file_path.empty()
		&& i_tensor_file::matches(i_tensor_path(), model_data.nr, model_data.nmo, I_tens.i_compact_, single_on_disk)) {
		if (!I_tens.i_tensor_save_path.empty())
			throw std::runtime_error("The I tensor is read from " + i_tensor_path().string() + " and should not be saved again to "
				+ I_tens.i_tensor_save_path.string() + ": a second copy costs as much memory as the first. Remove `save`/`safe` from the XCW settings.");
		//open() checks the header and throws if the shape does not match, which is what
		//stops a tensor from a different structure being used by accident.
		const char* source = "";
		bool automatic = false;
		const size_t w = items_within_budget(static_cast<size_t>(model_data.nr),
			i_tensor_file::block_bytes(I_tens.i_compact_, single_on_disk), i_budget(source, automatic));
		I_tens.i_streamed_ = w != 0;
		//Streamed, the window the budget allows, as decide_i_storage gives a built tensor; held, the
		//read below only passes through the window when it narrows
		I_tens.i_window_ = I_tens.i_streamed_ ? static_cast<int>(std::min(w, static_cast<size_t>(model_data.nr))) : std::max(1, std::min(model_data.nr, 64));
		open_i_stream_for_reading(i_tensor_path());
		// Figure out which pair of mu,nu each compact index corresponds to
		I_tens.i_pair_mu_ = I_tens.i_file_.pair_mu();
		I_tens.i_pair_nu_ = I_tens.i_file_.pair_nu();
		//The file's element type is kept as it is: a single-precision tensor cannot regain
		//anything by widening, and a double one is narrowed only on request
		const char* f = std::getenv("NOSPHERA2_XCW_I_FLOAT"); // Flawfinder: ignore
		I_tens.i_float_ = single_on_disk || (!I_tens.i_streamed_ && (I_tens.tensor_single || (f && std::atoi(f) != 0)));
		std::cout << "I tensor read from " << i_tensor_path().string()
			<< " (" << (i_tensor_file::total_bytes(model_data.nr, I_tens.i_compact_, single_on_disk) / 1048576.0)
			<< " MB" << (single_on_disk ? ", single precision" : "") << "), not recomputed"
			<< (I_tens.i_streamed_ ? ", read a window at a time" : ", held in memory");
		if (I_tens.i_streamed_) std::cout << " (" << I_tens.i_window_ << " of " << model_data.nr << " reflections resident)";
		std::cout << std::endl;
		if (I_tens.i_float_ && !single_on_disk)
			std::cout << "NOTE: the tensor on disk is double precision; it is narrowed to single as i_float asks" << std::endl;
		if (single_on_disk && I_tens.tensor_double)
			std::cout << "NOTE: the tensor on disk is single precision; i_double cannot widen it, it is used as stored" << std::endl;
		if (!I_tens.i_streamed_) {
			//Read straight into the tensor in the file's element type; only narrowing a double file goes through the window
			const size_t n = I_tens.i_compact_;
			if (I_tens.i_float_) I_tens.I32.resize(static_cast<size_t>(model_data.nr) * n);
			else I_tens.I.resize(static_cast<size_t>(model_data.nr) * n);
			for (int r0 = 0; r0 < model_data.nr; r0 += I_tens.i_window_) {
				const int r1 = std::min(model_data.nr, r0 + I_tens.i_window_);
				const size_t o = static_cast<size_t>(r0) * n;
				if (!I_tens.i_float_) I_tens.i_file_.read(r0, r1, I_tens.I.data() + o);
				else if (single_on_disk) I_tens.i_file_.read(r0, r1, I_tens.I32.data() + o);
				else {
					I_tens.i_file_.load(r0, r1);
					std::copy_n(I_tens.i_file_.block(r0), static_cast<size_t>(r1 - r0) * n, I_tens.I32.data() + o);
				}
			}
			I_tens.i_file_.close();
		}
	}
	// Path for computing the I tensor (through build_I)
	else {
		double time_taken;
		long long screen_counter = 0;
		long long skipped_grids = 0;
		build_I(ao_data_shells, time_taken, screen_counter, skipped_grids);
		if (!(opt->no_date)) {
			std::cout << std::fixed << std::setprecision(2) << "Time taken for XCW integrals: " << time_taken << " seconds. \n";
		}
		start_i_save();
		std::cout << std::fixed << std::setprecision(2) << "Screened out " << screen_counter << " unique pairs of mu, nu (" << static_cast<size_t>(screen_counter) / (static_cast<double>(model_data.nmo * (model_data.nmo + 1)) / 2) * 100.00 << "%) \n";
		std::cout << std::fixed << std::setprecision(2) << "Skipped evaluation of " << skipped_grids << " grids (" << static_cast<double>(skipped_grids) / ((static_cast<double>(model_data.nmo * (model_data.nmo + 1)) / 2) * model_data.nr_enlarged * model_data.ncen) * 100.00 << "%) \n";

	}
	return I_tens;
	// closing function
}

//Whether the I tensor is held or streamed, and the largest window that fits the budget
//The tensor is written on a thread while the refinement runs. Nothing in the SCF modifies
//it - both walks only read - so a reader alongside them needs no lock, and the run pays only
//the disk bandwidth, which it is not competing for while it works out of memory.
void structure_factors::start_i_save()
{
	//A streamed tensor was built straight into the save file, see decide_i_storage
	if (I_tens.i_tensor_save_path.empty() || I_tens.i_streamed_) return;
	const size_t packed = I_tens.i_compact_;
	const int nr = model_data.nr;
	const std::filesystem::path path = I_tens.i_tensor_save_path;
	const bool from_float = I_tens.i_float_;
	std::cout << "Writing the I tensor to " << path.string()
		<< " in the background; a later run can `read " << path.string()
		<< "` instead of building it" << std::endl;
	i_writer_ = std::thread([this, path, nr, packed, from_float]() {
		try {
			i_tensor_file out;
			out.create(path, nr, model_data.nmo, I_tens.i_pair_mu_, I_tens.i_pair_nu_, from_float);
			for (int r = 0; r < nr; r++) {
				if (from_float) out.write_block(r, I_tens.I32.data() + static_cast<size_t>(r) * packed);
				else out.write_block(r, I_tens.I.data() + static_cast<size_t>(r) * packed);
			}
			out.finish_write();
		}
		catch (const std::exception& e) { i_writer_error_ = e.what(); }
		catch (...) { i_writer_error_ = "unknown error"; }
		});
}

//Called before the run ends, and before anything that could invalidate the tensor. A thread
//left running past main is the bug 939268f was about; this one also holds a file handle.
void structure_factors::finish_i_save()
{
	if (!i_writer_.joinable()) return;
	i_writer_.join();
	if (!i_writer_error_.empty())
		std::cout << "Could not write the I tensor to "
			<< I_tens.i_tensor_save_path.string() << ": " << i_writer_error_
			<< " (the refinement itself is unaffected)" << std::endl;
	else
		std::cout << "I tensor written to " << I_tens.i_tensor_save_path.string() << std::endl;
}

size_t structure_factors::i_budget(const char*& source, bool& automatic) const {
	size_t budget = I_tens.i_tensor_max_mb * 1024ULL * 1024ULL;
	source = "i_tensor_mb";
	automatic = false;
	if (budget == 0 && opt->mem_given && opt->mem > 0.0) {
		budget = static_cast<size_t>(opt->mem * 1024.0 * 1024.0);
		source = "-mem";
	}
	if (budget == 0) {
		const size_t avail = available_memory_bytes();
		if (avail > 0) {
			//Four fifths: the SCF matrices, the grids and OCC's own allocations live in the
			//rest, and a tensor that only just fits would page rather than run. What the
			//process can have is a platform question - a cgroup here, a job object on
			//Windows, page classes on a Mac - and lives in convenience.cpp.
			budget = avail / 5 * 4;
			source = "four fifths of the memory this job can have";
			automatic = true;
		}
	}
	return budget;
}

void structure_factors::decide_i_storage() {
	const size_t per_block = i_tensor_file::block_bytes(I_tens.i_compact_, I_tens.i_float_);
	const size_t total = i_tensor_file::total_bytes(model_data.nr, I_tens.i_compact_, I_tens.i_float_);

	//The settings file budget wins, then -mem; with neither, what the process can actually
	//have. Left to a keyword this is the single most expensive decision in an XCW run and
	//the wrong answer is silent: both SCF walks re-read the whole tensor every iteration, so
	//streaming a tensor that would have fit is many times slower on the stage a
	//lambda scan spends its life in.
	//Nobody should have to know that to get it right.
	const char* source = "";
	bool automatic = false;
	const size_t budget = i_budget(source, automatic);
	//items_within_budget returns 0 for "hold everything": no file, no re-read twice per SCF iteration
	const size_t w = items_within_budget(static_cast<size_t>(model_data.nr), per_block, budget);
	I_tens.i_streamed_ = (w != 0);
	if (!I_tens.i_streamed_) {
		//Announced only when someone asked about memory: saying it unconditionally
		//shifts every reference output by a line, and the automatic budget would say it on
		//every run - including the reference tests, which is why it stays quiet there.
		if ((budget > 0 && !automatic) || ProgressBar::report_counts) {
			std::cout << std::fixed << std::setprecision(2)
				<< "I tensor held in memory: " << (total / 1048576.0) << " MB"
				<< (I_tens.i_float_ ? " (single precision)" : "");
			if (budget > 0)
				std::cout << " (fits the " << (budget / 1048576.0) << " MB " << source << " budget)";
			std::cout << std::endl;
		}
		return;
	}
	I_tens.i_window_ = static_cast<int>(std::min(w, static_cast<size_t>(model_data.nr)));
	//Streamed straight into the `save` file when there is one, so saving it costs nothing
	const std::filesystem::path path = I_tens.i_tensor_save_path.empty() ? i_tensor_path() : I_tens.i_tensor_save_path;
	I_tens.i_file_.create(path, model_data.nr, model_data.nmo, I_tens.i_pair_mu_, I_tens.i_pair_nu_, I_tens.i_float_);
	std::cout << std::fixed << std::setprecision(2)
		<< "I tensor streamed to disk: " << (total / 1048576.0) << " MB" << (I_tens.i_float_ ? " (single precision)" : "") << " total, "
		<< I_tens.i_window_ << " of " << model_data.nr << " reflections resident ("
		<< (I_tens.i_window_ * per_block / 1048576.0) << " MB) to fit " << source
		<< " (" << (budget / 1048576.0) << " MB)" << std::endl;
	if (I_tens.i_window_ == 1 && per_block > budget)
		std::cout << "  NOTE: one reflection alone is " << (per_block / 1048576.0)
		<< " MB, over the budget. Running one at a time." << std::endl;
	if (!I_tens.i_tensor_save_path.empty())
		std::cout << "The streamed I tensor is saved in " << path.string() << "; a later run can `read "
		<< path.string() << "` instead of building it" << std::endl;
}

std::filesystem::path structure_factors::i_tensor_path() const {
	return I_tens.i_tensor_file_path.empty()
		? std::filesystem::path(i_tensor_default)
		: I_tens.i_tensor_file_path;
}

void structure_factors::open_i_stream_for_reading(const std::filesystem::path& p) {
	I_tens.i_file_.open(p, static_cast<size_t>(I_tens.i_window_));
}

//One tile of the CPU I tensor, C = A * B^T row-major with k the block's points
static void tile_gemm(const int m, const int n, const int k, const double* a, const double* b, double* c)
{
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, n, k, 1.0, a, k, b, k, 0.0, c, n);
}

static void tile_gemm(const int m, const int n, const int k, const float* a, const float* b, float* c)
{
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, n, k, 1.0f, a, k, b, k, 0.0f, c, n);
}

void structure_factors::build_I(const std::vector<ao_data>& ao_data_shells, double& time_taken, long long& screen_counter, long long& skipped_grids_) {
	long long skipped_grids = 0;
	const int packed_size = (model_data.nmo * (model_data.nmo + 1)) / 2;
	int at = 0, mu = 0, nu = 0, r = 0, s = 0, r_asym = 0;

	cvec XCW_integrals;
	cvec2 XCW_integral_old;

	// Grid setup
	GridConfiguration config;
	config.accuracy = opt->accuracy;
	config.partition_type = opt->partition_type;
	config.no_density_eval = true;
	config.pbc = opt->pbc;
	config.debug = opt->debug;
	config.all_charges = opt->all_charges;
	GridManager grid_manager(config);
	WFN dummy_wave;
	for (int at = 0; at < model_data.ncen; at++) {
		atom temp_atom;
		temp_atom.set_coordinate(0, asym_atoms[at].pos[0]);
		temp_atom.set_coordinate(1, asym_atoms[at].pos[1]);
		temp_atom.set_coordinate(2, asym_atoms[at].pos[2]);
		temp_atom.set_charge(asym_atoms[at].type);
		dummy_wave.push_back_atom(temp_atom);
	}
	std::shared_ptr<BasisSet> basis = BasisSetLibrary::get_basis_set(I_tens.basis_set_name);
	load_basis_into_WFN(dummy_wave, basis, false, true);
	dummy_wave.delete_unoccupied_MOs();
	bvec needs_grid(model_data.ncen, true);
	ivec asym_atom_list(model_data.ncen);
	for (int at = 0; at < model_data.ncen; at++) {
		asym_atom_list[at] = at;
	}
	grid_manager.setup3DGridsForMolecule(dummy_wave, asym_atom_list, needs_grid, unit_cell);

	bool equal = false;

	GridData& GD = grid_manager.getGridData();
	vec2* grids = grid_manager.getNeedsHelper() ? GD.helper_grids.data() : GD.atomic_grids.data();
	vec2 d1, d2, d3, weights;
	const int n_grids = grid_manager.getNeedsHelper() ? GD.helper_grids.size() : GD.atomic_grids.size();
	for (int g = 0; g < n_grids; g++) {
		std::fill(grids[g][GridData::GridIndex::WFN_DENSITY].begin(), grids[g][GridData::GridIndex::WFN_DENSITY].end(), 1.0);
	}
	grid_manager.getDensityVectors(dummy_wave, asym_atom_list, d1, d2, d3, weights);
	const int* points = grid_manager.getNeedsHelper() ? GD.helper_num_points_per_atom.data() : GD.num_points_per_atom.data();
	const int total_points = grid_manager.getTotalGridPoints();
	std::cout << "Total number of grid points after pruning: " << total_points << std::endl;

	ivec2 asym_lookup(model_data.nr);
	for (r = 0; r < model_data.nr; r++) {
		asym_lookup[r] = generate_asym_lookup(r);
	}
	const unsigned int num_syms = asym_lookup[0].size();

	if (std::getenv("NOSPHERA2_DEBUG_LOOKUP")) { // Flawfinder: ignore
		long long misses = 0, total = 0;
		for (r = 0; r < model_data.nr; r++) {
			for (int s = 0; s < static_cast<int>(num_syms); s++) {
				total++;
				if (asym_lookup[r][s] == 0) {
					// index 0 is ambiguous: a genuine hit on hkl_enlarged's first entry, or the
					// silent "not found" fallback in generate_asym_lookup. Recompute by hand to
					// tell them apart.
					auto it = scatter_data.hkl.begin();
					std::advance(it, r);
					ivec3 rots = unit_cell.get_sym();
					i3 tempv{ 0,0,0 };
					const i3& hkl_temp = *it;
					for (int h = 0; h < 3; h++)
						for (int j = 0; j < 3; j++)
							tempv[j] += hkl_temp[h] * rots[j][h][s];
					if (scatter_data.hkl_enlarged.find(tempv) == scatter_data.hkl_enlarged.end()) misses++;
				}
			}
		}
		std::cerr << "asym_lookup misses: " << misses << " of " << total
			<< "  (hkl_enlarged size " << scatter_data.hkl_enlarged.size() << ")" << std::endl;
		std::cerr.flush();
		std::exit(0);
	}

	vec2 grid_positions(model_data.ncen);
	for (int at = 0; at < model_data.ncen; at++) {
		grid_positions[at] = { dummy_wave.get_atom_pos(at)[0], dummy_wave.get_atom_pos(at)[1], dummy_wave.get_atom_pos(at)[2] };
	}

	// Precompute screening
	std::chrono::high_resolution_clock::time_point screening_start = std::chrono::high_resolution_clock::now();
	ivec2 skip(model_data.nmo, ivec(model_data.nmo, 0));
	{
		const double e_tol = 0.0005;
		double estimate = 0.0;
		const double sqrt_inv_four_pi = std::sqrt(constants::INV_FOUR_PI);
		for (mu = 0; mu < model_data.nmo; mu++) {
			const std::vector<primitive>& mu_primitives = ao_data_shells[mu].prims;
			const int mu_type = mu_primitives[0].get_type();
			const double mu_type_half = static_cast<double>(mu_type) * 0.5;
			const double min_mu_type_half_exp = std::exp(-mu_type_half);
			const double sph_harmonic_max_mu = constants::spherical_harmonic_max(mu_type, ao_data_shells[mu].m);
			for (nu = mu; nu < model_data.nmo; nu++) {
				const double dist2 = (ao_data_shells[mu].pos[0] - ao_data_shells[nu].pos[0]) * (ao_data_shells[mu].pos[0] - ao_data_shells[nu].pos[0]) +
					(ao_data_shells[mu].pos[1] - ao_data_shells[nu].pos[1]) * (ao_data_shells[mu].pos[1] - ao_data_shells[nu].pos[1]) +
					(ao_data_shells[mu].pos[2] - ao_data_shells[nu].pos[2]) * (ao_data_shells[mu].pos[2] - ao_data_shells[nu].pos[2]);
				if (dist2 < 1e-5) {
					continue;
				}
				const std::vector<primitive>& nu_primitives = ao_data_shells[nu].prims;
				const double nu_type = static_cast<double>(nu_primitives[0].get_type());
				const double nu_type_half = static_cast<double>(nu_type) * 0.5;
				const double min_nu_type_half_exp = std::exp(-nu_type_half);
				const double sph_harmonic_max_nu = constants::spherical_harmonic_max(nu_type, ao_data_shells[nu].m);
				estimate = 0.0;
				for (int k = 0; k < mu_primitives.size(); k++) {
					const double k_coef = mu_primitives[k].get_coef();
					const double k_exp = mu_primitives[k].get_exp();
					const double N_k = mu_type == 0 ? sqrt_inv_four_pi : sph_harmonic_max_mu * min_mu_type_half_exp * std::pow(mu_type / k_exp, mu_type_half);
					for (int l = 0; l < nu_primitives.size(); l++) {
						const double l_coef = nu_primitives[l].get_coef();
						const double l_exp = nu_primitives[l].get_exp();
						const double inv_exp_sum = 1.0 / (k_exp + l_exp);
						const double N_l = nu_type == 0 ? sqrt_inv_four_pi : sph_harmonic_max_nu * min_nu_type_half_exp * std::pow(nu_type / l_exp, nu_type_half);
						const double combined_coeffs = std::abs(k_coef * l_coef);
						const double N_kl = N_k * N_l * std::pow(constants::TWO_PI * inv_exp_sum, 1.5);
						const double gamma = 0.5 * k_exp * l_exp * inv_exp_sum;
						const double gaussian = std::exp(-gamma*dist2);
						estimate += combined_coeffs * N_kl * gaussian;
					}
				}
				if (estimate < e_tol) {
					skip[mu][nu] = 1;
				}
			}
		}
	}

	std::chrono::high_resolution_clock::time_point screening_end = std::chrono::high_resolution_clock::now();
	if (!opt->no_date) {
		std::cout << std::fixed << std::setprecision(5) << "Time taken for AO screening: " << std::chrono::duration_cast<std::chrono::microseconds>(screening_end - screening_start).count() << " microseconds." << "\n";
	}

	// Grid screening
	constexpr double maximum_ao_grid_cutoff = 12;
	constexpr double minimum_ao_grid_cutoff = 11;
	double minimum_primitive_exponent = std::numeric_limits<double>::max();
	vec ao_grid_cutoff_squared(model_data.nmo);
	{
		vec ao_minimum_exponent(model_data.nmo);
		for (int ao = 0; ao < model_data.nmo; ao++) {
			double min_exp = std::numeric_limits<double>::max();
			for (const primitive& prim : ao_data_shells[ao].prims) {
				min_exp = std::min(min_exp, prim.get_exp());
			}
			ao_minimum_exponent[ao] = min_exp;
			minimum_primitive_exponent = std::min(minimum_primitive_exponent, min_exp);
		}
		for (int ao = 0; ao < model_data.nmo; ao++) {
			const double adaptive_cutoff = maximum_ao_grid_cutoff * std::sqrt(minimum_primitive_exponent / ao_minimum_exponent[ao]);
			const double cutoff = std::clamp(adaptive_cutoff, minimum_ao_grid_cutoff, maximum_ao_grid_cutoff);
			ao_grid_cutoff_squared[ao] = cutoff * cutoff;
		}
	}

	const int n_atom_grids = std::min(n_grids, model_data.ncen);
	ivec2 ao_prefix_end(model_data.nmo, ivec(n_atom_grids));
	bvec2 ao_within_cutoff(model_data.nmo, bvec(n_atom_grids));

	// Compute radial distance for every grid point
	vec2 grid_radial_distances(n_atom_grids);
	for (int g = 0; g < n_atom_grids; g++) {
		vec& radial_distances = grid_radial_distances[g];
		radial_distances.resize(points[g]);
		const double* x_ptr = grids[g][GridData::GridIndex::X].data();
		const double* y_ptr = grids[g][GridData::GridIndex::Y].data();
		const double* z_ptr = grids[g][GridData::GridIndex::Z].data();
		for (int p = 0; p < points[g]; p++) {
			const double dx = x_ptr[p] - grid_positions[g][0];
			const double dy = y_ptr[p] - grid_positions[g][1];
			const double dz = z_ptr[p] - grid_positions[g][2];
			radial_distances[p] = std::sqrt(dx * dx + dy * dy + dz * dz);
		}
	}

#pragma omp parallel for schedule(static)
	for (int ao = 0; ao < model_data.nmo; ao++) {
		const ao_data& ao_shell = ao_data_shells[ao];
		const double cutoff2 = ao_grid_cutoff_squared[ao];
		const double cutoff = std::sqrt(cutoff2);
		bvec& within_row = ao_within_cutoff[ao];
		ivec& prefix_row = ao_prefix_end[ao];
		// Compute distance between grid center and AO center
		for (int g = 0; g < n_atom_grids; g++) {
			const double dx = grid_positions[g][0] - ao_shell.pos[0];
			const double dy = grid_positions[g][1] - ao_shell.pos[1];
			const double dz = grid_positions[g][2] - ao_shell.pos[2];
			const double d2 = dx * dx + dy * dy + dz * dz;
			within_row[g] = d2 <= cutoff2;
			prefix_row[g] = (d2 < 1e-12)
				? static_cast<int>(std::upper_bound(grid_radial_distances[g].begin(), grid_radial_distances[g].end(), cutoff) - grid_radial_distances[g].begin())
				: points[g];
		}
	}

	// Precompute AO values
	// add a timer for computation time of AO values
	std::chrono::high_resolution_clock::time_point start_AOs = std::chrono::high_resolution_clock::now();
	vec3 mu_vals(model_data.nmo, vec2(n_grids));
#pragma omp parallel for schedule(dynamic)
	for (mu = 0; mu < model_data.nmo; mu++) {
		const ao_data& mu_prims = ao_data_shells[mu];
		const std::vector<primitive>& mu_primitives = mu_prims.prims;
		const double mp0 = mu_prims.pos[0];
		const double mp1 = mu_prims.pos[1];
		const double mp2 = mu_prims.pos[2];
		for (int g = 0; g < n_grids; g++) {
			vec2& atom_grid = grids[g];
			const double* x_ptr = atom_grid[GridData::GridIndex::X].data();
			const double* y_ptr = atom_grid[GridData::GridIndex::Y].data();
			const double* z_ptr = atom_grid[GridData::GridIndex::Z].data();
			mu_vals[mu][g].resize(points[g]);
			double* local_mu_vals_ptr = mu_vals[mu][g].data();
			const int prefix_end = g < n_atom_grids ? ao_prefix_end[mu][g] : points[g];
			for (int p = 0; p < prefix_end; p++) {
				d4 d_mu{ x_ptr[p] - mp0, y_ptr[p] - mp1 , z_ptr[p] - mp2 , 0 };
				d_mu[3] = d_mu[0] * d_mu[0] + d_mu[1] * d_mu[1] + d_mu[2] * d_mu[2];
				const double root_d3 = std::sqrt(d_mu[3]);
				if (d_mu[3] > ao_grid_cutoff_squared[mu]) {
					local_mu_vals_ptr[p] = 0.0;
				}
				else {
					local_mu_vals_ptr[p] = dummy_wave.eval_ao(d_mu, mu_primitives, mu_prims.m, root_d3);
				}
			}
		}
	}
	std::chrono::high_resolution_clock::time_point end_AOs = std::chrono::high_resolution_clock::now();
	std::cout << "AO values calculated for all grids." << std::endl;
	if (!(opt->no_date))
		std::cout << "Time taken for AO values computation: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_AOs - start_AOs).count() << " milliseconds." << std::endl;
	//Morton-order every atom grid's points, so that a run of consecutive points is a
	//compact ball rather than a spherical shell. This is what OCC does
	//(occ/qm/spatial_grid_hierarchy.h) and what grid-based codes do generally, and the
	//permutation is the point of it: the grid arrives sorted by radius, so consecutive
	//points span a whole sphere and every AO reaching any part of it stays active. Measured
	//on the twisted ethylene at def2-TZVP, cutting the radial bands into chunks without
	//reordering moved the work by 2.7% and left n_active at 852; the innermost eighth of a
	//grid, which is compact because its radius is small, needs 370 AOs against 803.
	//
	//Work is sum over blocks of n_active^2 * points, so this is quadratic in what it saves.
	//Sums over points are order independent, so the reordering changes no result.
	//Compaction and the AO threshold are worth nothing apart and a great deal together:
	//reordering alone is slower (the blocks shrink and n_active does not), the threshold
	//alone barely moves the work, and together they are faster with the GooF, energies and
	//convergence lines identical.
	//So only reorder when the threshold can actually prune: at -acc 4 cutoff() is 1e-30 and
	//nothing would be dropped, and paying the compaction cost for that would make asking for
	//more accuracy slower for no reason.
	const double ao_block_threshold = [&] {
		const char* e = std::getenv("NOSPHERA2_ITENSOR_AO_TOL"); // Flawfinder: ignore
		if (e) { const double v = std::atof(e); return v >= 0.0 ? v : 0.0; }
		return cutoff(opt->accuracy);
	}();
	const bool morton_applied = (std::getenv("NOSPHERA2_ITENSOR_NO_MORTON") == nullptr) // Flawfinder: ignore
		&& ao_block_threshold >= 1e-20;
	if (morton_applied) {
#pragma omp parallel for schedule(dynamic)
		for (int g = 0; g < n_atom_grids; g++) {
			const int npts = points[g];
			if (npts < 2) continue;
			vec2& grid = grids[g];
			double* xs = grid[GridData::GridIndex::X].data();
			double* ys = grid[GridData::GridIndex::Y].data();
			double* zs = grid[GridData::GridIndex::Z].data();
			std::array<double, 3> lo{ 1e30, 1e30, 1e30 }, hi{ -1e30, -1e30, -1e30 };
			for (int p = 0; p < npts; p++) {
				lo[0] = std::min(lo[0], xs[p]); hi[0] = std::max(hi[0], xs[p]);
				lo[1] = std::min(lo[1], ys[p]); hi[1] = std::max(hi[1], ys[p]);
				lo[2] = std::min(lo[2], zs[p]); hi[2] = std::max(hi[2], zs[p]);
			}
			auto spread = [](const unsigned int v) {
				unsigned long long x = v & 0x1fffffu;   //21 bits, three of them interleave to 63
				x = (x | (x << 32)) & 0x1f00000000ffffull;
				x = (x | (x << 16)) & 0x1f0000ff0000ffull;
				x = (x | (x << 8))  & 0x100f00f00f00f00full;
				x = (x | (x << 4))  & 0x10c30c30c30c30c3ull;
				x = (x | (x << 2))  & 0x1249249249249249ull;
				return x;
			};
			std::vector<std::pair<unsigned long long, int>> keyed(npts);
			for (int p = 0; p < npts; p++) {
				unsigned int c[3];
				const double v[3] = { xs[p], ys[p], zs[p] };
				for (int d = 0; d < 3; d++) {
					const double span = hi[d] - lo[d];
					const double t = span > 1e-12 ? (v[d] - lo[d]) / span : 0.0;
					c[d] = static_cast<unsigned int>(std::min(2097151.0, std::max(0.0, t * 2097151.0)));
				}
				keyed[p] = { spread(c[0]) | (spread(c[1]) << 1) | (spread(c[2]) << 2), p };
			}
			std::sort(keyed.begin(), keyed.end());
			ivec perm(npts);
			for (int p = 0; p < npts; p++) perm[p] = keyed[p].second;
			auto apply = [&](double* a) {
				vec tmp(npts);
				for (int p = 0; p < npts; p++) tmp[p] = a[perm[p]];
				std::copy(tmp.begin(), tmp.end(), a);
			};
			apply(xs); apply(ys); apply(zs);
			apply(grid[GridData::GridIndex::WEIGHT].data());
			//The coordinates and weights the phase factor and the GEMM actually use are
			//these, taken from getDensityVectors above and not the grid arrays: reordering
			//the AO values without them pairs each value with another point's coordinate,
			//which is wrong everywhere rather than only where a screening decision was made.
			if (static_cast<int>(d1[g].size()) >= npts) apply(d1[g].data());
			if (static_cast<int>(d2[g].size()) >= npts) apply(d2[g].data());
			if (static_cast<int>(d3[g].size()) >= npts) apply(d3[g].data());
			if (static_cast<int>(weights[g].size()) >= npts) apply(weights[g].data());
			for (int mu = 0; mu < model_data.nmo; mu++) {
				vec& v = mu_vals[mu][g];
				if (static_cast<int>(v.size()) == npts) apply(v.data());
				else if (!v.empty()) {
					//Values were only filled to the radial prefix; the tail is zero and the
					//permutation mixes the two, so grow it before reordering.
					v.resize(npts, 0.0);
					apply(v.data());
				}
			}
			//Radial distance follows its point, and the band bounds below are recomputed
			//from it - they are no longer monotone, which is what the chunking wants.
			vec& rd = grid_radial_distances[g];
			if (static_cast<int>(rd.size()) == npts) apply(rd.data());
		}
	}



	//NOSPHERA2_ITENSOR_AOSTATS=1: how much of each block's AO set is actually carrying
	//anything. The active set comes from a cutoff clamped into an 11-12 bohr band
	//(std::clamp above), so it is set by the distance between two atom centres and barely
	//by the block - a 266-point block keeps 756 of 852 AOs and a 7968-point one keeps 803.
	//Work is sum over blocks of na^2 * points, so what an OCC-style per-batch bounding
	//sphere would save is quadratic in whatever this measures. The AO values are already
	//computed here, so the honest number is a max over the points they hold, not an estimate.
	if (std::getenv("NOSPHERA2_ITENSOR_AOSTATS")) { // Flawfinder: ignore
		for (int g = 0; g < n_atom_grids; g++) {
			const int npts = points[g];
			if (npts <= 0) continue;
			//max |chi| per AO over this grid, and over the first eighth of it as a stand-in
			//for a compact spatial batch
			const int batch = std::max(1, npts / 8);
			long long active_full = 0, active_batch = 0, kept = 0;
			for (int ao = 0; ao < model_data.nmo; ao++) kept += ao_within_cutoff[ao][g] ? 1 : 0;
			for (int ao = 0; ao < model_data.nmo; ao++) {
				const vec& v = mu_vals[ao][g];
				if (v.empty()) continue;
				double mx_full = 0.0, mx_batch = 0.0;
				const int end = std::min<int>(static_cast<int>(v.size()), npts);
				for (int p = 0; p < end; p++) {
					const double a = std::abs(v[p]);
					mx_full = std::max(mx_full, a);
					if (p < batch) mx_batch = std::max(mx_batch, a);
				}
				if (mx_full > 1e-10) active_full++;
				if (mx_batch > 1e-10) active_batch++;
			}
			std::fprintf(stderr, "aostats grid %d: points %d  nmo %d  kept by cutoff %d"
				"  carrying |chi|>1e-10: whole grid %lld  first eighth %lld\n",
				g, npts, model_data.nmo, static_cast<int>(kept),
				active_full, active_batch);
		}
	}

	ivec2 active_grids(packed_size);
	ivec skipped_grids_per_pair(packed_size, 0);
	for (mu = 0; mu < model_data.nmo; mu++) {
		const std::vector<bool>& mu_within = ao_within_cutoff[mu];
		for (nu = mu; nu < model_data.nmo; nu++) {
			if (skip[mu][nu]) {
				continue;
			}
			const bvec& nu_within = ao_within_cutoff[nu];
			const size_t pair_idx = tri_index(mu, nu);
			ivec& pair_grids = active_grids[pair_idx];
			pair_grids.reserve(n_atom_grids);
			for (int g = 0; g < n_atom_grids; g++) {
				if (mu_within[g] && nu_within[g]) {
					pair_grids.push_back(g);
				}
			}
			skipped_grids_per_pair[pair_idx] = n_atom_grids - static_cast<int>(pair_grids.size());
		}
	}
	//Only the pairs that survive the screening are stored, in (mu, nu) order; tri_compact
	//takes a packed slot to its stored index, -1 when screened out
	ivec tri_compact(packed_size, -1);
	I_tens.i_pair_mu_.clear();
	I_tens.i_pair_nu_.clear();
	for (mu = 0; mu < model_data.nmo; mu++)
		for (nu = mu; nu < model_data.nmo; nu++)
			if (!skip[mu][nu]) {
				tri_compact[tri_index(mu, nu)] = static_cast<int>(I_tens.i_pair_mu_.size());
				I_tens.i_pair_mu_.push_back(mu);
				I_tens.i_pair_nu_.push_back(nu);
			}
	I_tens.i_compact_ = I_tens.i_pair_mu_.size();
	ivec2 grid_active_aos(n_atom_grids);
	vec2 grid_ao_values(n_atom_grids);
#pragma omp parallel for schedule(dynamic)
	for (int g = 0; g < n_atom_grids; g++) {
		ivec& active_aos = grid_active_aos[g];
		vec& values = grid_ao_values[g];
		for (int mu = 0; mu < model_data.nmo; mu++) {
			if (!ao_within_cutoff[mu][g]) {
				continue;
			}
			active_aos.push_back(mu);
			const vec& values_for_ao = mu_vals[mu][g];
			values.insert(values.end(), values_for_ao.begin(), values_for_ao.end());
		}
	}

	// Tile the grid points for each atom into blocks of size 64
	struct MatrixTile {
		int row_start;
		int row_count;
		int col_start;
		int col_count;
		size_t result_offset;
	};
	struct GridBlock {
		int point_start;
		int point_count;
		ivec active_aos;
		vec ao_values;
		std::vector<float> ao_values_f;
		std::vector<MatrixTile> matrix_tiles;
		int tile_result_size = 0;
	};
	//128 rows a tile: at 64 the GEMM calls are overhead-bound in both precisions, and above
	//it a double tile pair leaves the core's L2 while single precision stays flat to 256
	constexpr int screened_tile_size = 128;
	std::vector<std::vector<GridBlock>> grid_blocks(n_atom_grids);
	auto make_matrix_tiles = [&](GridBlock& block) {
		const int n_active = static_cast<int>(block.active_aos.size());
		size_t result_offset = 0;
		for (int row_start = 0; row_start < n_active; row_start += screened_tile_size) {
			const int row_count = std::min(screened_tile_size, n_active - row_start);
			for (int col_start = row_start; col_start < n_active; col_start += screened_tile_size) {
				const int col_count = std::min(screened_tile_size, n_active - col_start);
				bool needed = false;
				for (int row = row_start; row < row_start + row_count && !needed; ++row) {
					const int first_col = (row_start == col_start) ? row : col_start;
					for (int col = first_col; col < col_start + col_count; ++col) {
						if (!skip[block.active_aos[row]][block.active_aos[col]]) {
							needed = true;
							break;
						}
					}
				}
				if (needed) {
					block.matrix_tiles.push_back({ row_start, row_count, col_start, col_count, result_offset });
					result_offset += static_cast<size_t>(row_count) * col_count;
				}
			}
		}
		block.tile_result_size = static_cast<int>(result_offset);
		};
	//NOSPHERA2_ITENSOR_SKIPSTATS=1: what a spatial reordering of the active AOs would buy.
	//Half the mu,nu pairs are screened out and none of the 64x64 tiles are, because AO index
	//order is atom order as the CIF lists them and dead pairs land scattered. Sorting the
	//active AOs of a block along a Morton curve over their centres puts distant atoms in
	//distant tiles, which is the only way a tile becomes wholly dead. Measured here, per
	//block, before anyone writes a kernel that depends on it.
	auto skipstats = [&](const int g, const ivec& active, const int npoints) {
		if (!std::getenv("NOSPHERA2_ITENSOR_SKIPSTATS")) return; // Flawfinder: ignore
		const int na = static_cast<int>(active.size());
		if (na < 2) return;
		auto tiles_alive = [&](const ivec& order, const int T) {
			int alive = 0, total = 0;
			for (int r0 = 0; r0 < na; r0 += T)
				for (int c0 = r0; c0 < na; c0 += T) {
					total++;
					bool needed = false;
					for (int r = r0; r < std::min(r0 + T, na) && !needed; r++) {
						const int first = (r0 == c0) ? r : c0;
						for (int c = first; c < std::min(c0 + T, na); c++)
							if (!skip[order[r]][order[c]]) { needed = true; break; }
					}
					if (needed) alive++;
				}
			return std::pair<int, int>{ alive, total };
		};
		//Morton key over the AO centres, 10 bits per axis on the block's own bounding box
		std::array<double, 3> lo{ 1e30, 1e30, 1e30 }, hi{ -1e30, -1e30, -1e30 };
		for (const int ao : active)
			for (int d = 0; d < 3; d++) {
				lo[d] = std::min(lo[d], ao_data_shells[ao].pos[d]);
				hi[d] = std::max(hi[d], ao_data_shells[ao].pos[d]);
			}
		auto spread = [](unsigned int v) {
			unsigned long long x = v & 0x3ffu;
			x = (x | (x << 16)) & 0x30000ffull; x = (x | (x << 8)) & 0x300f00full;
			x = (x | (x << 4)) & 0x30c30c3ull;  x = (x | (x << 2)) & 0x9249249ull;
			return x;
		};
		std::vector<std::pair<unsigned long long, int>> keyed;
		keyed.reserve(na);
		for (const int ao : active) {
			unsigned int c[3];
			for (int d = 0; d < 3; d++) {
				const double span = hi[d] - lo[d];
				const double t = span > 1e-12 ? (ao_data_shells[ao].pos[d] - lo[d]) / span : 0.0;
				c[d] = static_cast<unsigned int>(std::min(1023.0, std::max(0.0, t * 1023.0)));
			}
			keyed.emplace_back(spread(c[0]) | (spread(c[1]) << 1) | (spread(c[2]) << 2), ao);
		}
		std::sort(keyed.begin(), keyed.end());
		ivec sorted_order(na), plain_order(na);
		for (int i = 0; i < na; i++) { sorted_order[i] = keyed[i].second; plain_order[i] = active[i]; }
		long long pairs = 0, dead = 0;
		for (int i = 0; i < na; i++)
			for (int j = i; j < na; j++) { pairs++; dead += skip[active[i]][active[j]] ? 1 : 0; }
		std::fprintf(stderr, "skipstats grid %d: na %d points %d  pairs dead %.1f%%", g, na, npoints,
			100.0 * (double)dead / (double)pairs);
		for (const int T : { 32, 64, 128 }) {
			const auto [a0, t0] = tiles_alive(plain_order, T);
			const auto [a1, t1] = tiles_alive(sorted_order, T);
			std::fprintf(stderr, "  | T=%d tiles pruned: as-is %.1f%% morton %.1f%%", T,
				100.0 * (1.0 - (double)a0 / t0), 100.0 * (1.0 - (double)a1 / t1));
		}
		std::fprintf(stderr, "\n");
	};

	//Counted the way the mu,nu screening is, so a run says what this cost it as well as
	//what it saved: how many AO-block entries were dropped, and what that did to the work
	//the GEMMs actually do.
	long long ao_slots_carrying = 0, ao_slots_kept = 0;
	//ao_block_threshold is defined above, with the reordering it enables. What counts as
	//nothing is the run's -acc setting, not a number invented here: cutoff() is the same
	//ladder the scattering-factor code screens on, 1e-10 up to -acc 2 and 1e-14 at 3.
#pragma omp parallel for schedule(dynamic) reduction(+:ao_slots_carrying, ao_slots_kept)
	for (int g = 0; g < n_atom_grids; g++) {
		const ivec& active_aos = grid_active_aos[g];
		const vec& full_ao_values = grid_ao_values[g];
		const vec& radial_distances = grid_radial_distances[g];
		const int inner_end = static_cast<int>(std::upper_bound(radial_distances.begin(), radial_distances.end(), minimum_ao_grid_cutoff) - radial_distances.begin());
		const int middle_end = static_cast<int>(std::upper_bound(radial_distances.begin(), radial_distances.end(), maximum_ao_grid_cutoff) - radial_distances.begin());
		const std::array<int, 4> block_bounds{ 0, inner_end, middle_end, points[g] };
		//A block keeps every AO that is non-zero anywhere in it, and the work is
		//sum over blocks of n_active^2 * points, so the block's spatial extent is what sets
		//the cost. Three radial bands make the first one nearly the whole grid: measured on
		//the twisted ethylene at def2-TZVP, a 7934-point band keeps 803 of 852 AOs while its
		//innermost eighth needs 370. Cutting the bands into chunks is what OCC does with its
		//Morton leaves (occ/qm/spatial_grid_hierarchy.h, 128 points a leaf) and what every
		//grid-based code does for the same reason. The points come radially sorted, so
		//consecutive chunks are already spatially compact and nothing has to be permuted.
		//
		//NOSPHERA2_ITENSOR_CHUNK sets the target; 0 restores the three whole bands.
		const int chunk = [] {
			const char* e = std::getenv("NOSPHERA2_ITENSOR_CHUNK"); // Flawfinder: ignore
			return e ? std::atoi(e) : 1024;
		}();
		//Even chunks rather than a short tail: a 40-point remainder is a GEMM that costs a
		//launch and returns almost nothing.
		auto cut = [&](const int from, const int to, std::vector<std::pair<int, int>>& out) {
			const int n = to - from;
			if (n <= 0) return;
			if (chunk <= 0) { out.emplace_back(from, to); return; }
			const int pieces = std::max(1, (n + chunk - 1) / chunk);
			const int per = (n + pieces - 1) / pieces;
			for (int p0 = from; p0 < to; p0 += per) out.emplace_back(p0, std::min(p0 + per, to));
		};
		std::vector<std::pair<int, int>> spans;
		if (morton_applied) {
			//The three radial bands are what the point order was for, and after Morton
			//ordering it is gone: block_bounds comes from upper_bound over the radial
			//distances, which needs a sorted range and no longer has one. Left in, the
			//bounds come back arbitrary, a band with end below start is skipped, and its
			//points drop out of the integration entirely - the structure factors then move
			//far more than any screening would explain (GooF 3.82 -> 26.34 on the twisted
			//ethylene, which is how this was found). Cut the grid itself instead: the bands
			//existed to group points by cutoff regime and a compact chunk does that better.
			cut(0, points[g], spans);
		}
		else {
			for (int block_index = 0; block_index < 3; block_index++)
				cut(block_bounds[block_index], block_bounds[block_index + 1], spans);
		}
		for (const auto& [point_start, point_end] : spans) {
			const int point_count = point_end - point_start;
			if (point_count == 0) {
				continue;
			}
			GridBlock block{ point_start, point_count };
			for (int local_ao = 0; local_ao < static_cast<int>(active_aos.size()); local_ao++) {
				const double* full_row = full_ao_values.data() + static_cast<size_t>(local_ao) * points[g];
				//Not "is it exactly zero" but "does it carry anything here". The values were
				//only zeroed where the 11-12 bohr cutoff cut them off, so a function whose
				//value on this block is 1e-40 was counted as active and multiplied at full
				//cost: n_active stayed at 852 of 852 where the AOs actually carrying more
				//than 1e-10 numbered 370. Work is n_active^2 * points, so this is quadratic.
				//What every grid-based code does, and the threshold is the same kind of
				//number as the 5e-4 the pair screening above already accepts.
				double largest = 0.0;
				for (int p = point_start; p < point_end; p++)
					largest = std::max(largest, std::abs(full_row[p]));
				if (largest > 0.0) ao_slots_carrying++;
				if (largest > ao_block_threshold) {
					block.active_aos.push_back(active_aos[local_ao]);
					block.ao_values.insert(block.ao_values.end(), full_row + point_start, full_row + point_end);
				}
			}
			if (!block.active_aos.empty()) {
				ao_slots_kept += static_cast<long long>(block.active_aos.size());
				make_matrix_tiles(block);
				if (opt->cpu_itensor_fp32) block.ao_values_f.assign(block.ao_values.begin(), block.ao_values.end());
				skipstats(g, block.active_aos, block.point_count);
				grid_blocks[g].push_back(std::move(block));
			}
		}
	}
	//Said next to "Screened out ... unique pairs of mu, nu", because it is the same kind of
	//saving measured on the other axis: that one drops pairs whose product cannot reach the
	//grid, this one drops an AO from a block where it carries nothing. Gated on no_date like
	//the timing lines, so the reference outputs keep their shape.
	if (!(opt->no_date) && ao_slots_carrying > 0) {
		const long long dropped = ao_slots_carrying - ao_slots_kept;
		std::cout << std::fixed << std::setprecision(2)
			<< "Screened out " << dropped << " of " << ao_slots_carrying
			<< " AO-block entries (" << 100.0 * static_cast<double>(dropped)
			/ static_cast<double>(ao_slots_carrying) << "%) below "
			<< std::scientific << std::setprecision(0) << ao_block_threshold
			<< std::fixed << std::setprecision(2) << " on their block\n";

		//The screenings in one number. Each of the lines above counts what it removed on its
		//own axis - pairs, AO-block entries, whole grids - and none of them says what the
		//run will actually cost. This does: the I tensor's work is the sum over blocks of
		//n_active^2 times points, and the same sum with every AO on every point is what it
		//would be with no screening at all. The ratio is what the GEMMs were spared.
		double work_done = 0.0, work_unscreened = 0.0;
		for (int g = 0; g < n_atom_grids; g++)
			for (const GridBlock& b : grid_blocks[g]) {
				const double na = static_cast<double>(b.active_aos.size());
				work_done += na * na * b.point_count;
				work_unscreened += static_cast<double>(model_data.nmo) * model_data.nmo * b.point_count;
			}
		//Per reflection and symmetry operation the sum is a small number and says nothing;
		//what the run costs is that times both, so scale it before printing or the figure
		//reads as a thousandth of the truth.
		const double runs = static_cast<double>(model_data.nr) * static_cast<double>(num_syms);
		if (work_done > 0.0)
			std::cout << std::fixed << std::setprecision(2)
				<< "I tensor work after all screening: " << (work_done * runs / 1e12)
				<< " of " << (work_unscreened * runs / 1e12) << " Tflop-equivalents ("
				<< std::setprecision(2) << 100.0 * work_done / work_unscreened << "%, "
				<< (work_unscreened / work_done) << "x less than unscreened)\n";
	}

	//The whole cost of the device path in one number: sum over blocks of n_active^2 times
	//points, which is what the GEMMs do per reflection and symmetry operation. Printed under
	//-gflops so a chunk size can be judged without running a reflection.
	if (throughput::enabled()) {
		double work = 0.0;
		long long nblocks = 0, na_min = 1LL << 60, na_max = 0, pts_min = 1LL << 60, pts_max = 0;
		for (int g = 0; g < n_atom_grids; g++)
			for (const GridBlock& b : grid_blocks[g]) {
				const long long na = static_cast<long long>(b.active_aos.size());
				work += static_cast<double>(na) * na * b.point_count;
				nblocks++;
				na_min = std::min(na_min, na); na_max = std::max(na_max, na);
				pts_min = std::min<long long>(pts_min, b.point_count);
				pts_max = std::max<long long>(pts_max, b.point_count);
			}
		std::fprintf(stderr, "I tensor blocks: %lld, n_active %lld-%lld, points %lld-%lld, "
			"sum n_active^2 * points = %.3e (lower is less work per reflection)\n",
			nblocks, na_min, na_max, pts_min, pts_max, work);
	}
	ivec2().swap(grid_active_aos);
	vec2().swap(grid_ao_values);
	vec2().swap(grid_radial_distances);
	vec3().swap(mu_vals);

	std::optional<ProgressBar> pb;
	auto start = std::chrono::high_resolution_clock::now();

	// Bookkeeping skipped pairs and grids
	for (int mu = 0; mu < model_data.nmo; mu++) {
		for (int nu = mu; nu < model_data.nmo; nu++) {
			screen_counter += skip[mu][nu];
			if (!skip[mu][nu]) {
				skipped_grids += static_cast<long long>(num_syms) * skipped_grids_per_pair[tri_index(mu, nu)];
			}
		}
	}
	bool itensor_on_gpu = false;
	double itensor_gpu_dense_flops = 0.0;
	//The device and the CPU threads draw reflections from one counter, the CPU stopping
	//once the device would finish what is left before a thread finished one more.
	std::atomic<int> next_refl{0};
	std::atomic<long long> gpu_ns_per_refl{0};
	std::mutex i_write_mutex;
	std::thread gpu_thread;
#if defined(NOSPHERA2_USE_GPU) || defined(NOSPHERA2_USE_METAL)
	//Read with use_gpu rather than on its own, so -no_gpu means what it says. Checking the
	//pair here rather than clearing the flag at parse time keeps it order-independent.
	if (opt->gpu_itensor && opt->use_gpu) {
		//Flatten what the device needs: the AO values never change with the reflection,
		//so they are uploaded once and every reflection reuses them.
		ivec bg, bps, bpc, bna, goff(n_atom_grids + 1, 0);
		std::vector<long long> bao, baos;
		vec ao_all;
		ivec aos_all;
		for (int gg = 0; gg < n_atom_grids; gg++) goff[gg + 1] = goff[gg] + points[gg];
		for (int gg = 0; gg < n_atom_grids; gg++) {
			for (const GridBlock& blkk : grid_blocks[gg]) {
				bg.push_back(gg);
				bps.push_back(blkk.point_start);
				bpc.push_back(blkk.point_count);
				bna.push_back(static_cast<int>(blkk.active_aos.size()));
				bao.push_back(static_cast<long long>(ao_all.size()));
				baos.push_back(static_cast<long long>(aos_all.size()));
				ao_all.insert(ao_all.end(), blkk.ao_values.begin(), blkk.ao_values.end());
				aos_all.insert(aos_all.end(), blkk.active_aos.begin(), blkk.active_aos.end());
			}
		}
		ivec compact_flat(static_cast<size_t>(model_data.nmo) * model_data.nmo, -1);
		for (int m = 0; m < model_data.nmo; m++)
			for (int n = m; n < model_data.nmo; n++)
				compact_flat[static_cast<size_t>(m) * model_data.nmo + n] = tri_compact[tri_index(m, n)];
		vec fd1, fd2, fd3, fw;
		for (int gg = 0; gg < n_atom_grids; gg++) {
			fd1.insert(fd1.end(), d1[gg].begin(), d1[gg].begin() + points[gg]);
			fd2.insert(fd2.end(), d2[gg].begin(), d2[gg].begin() + points[gg]);
			fd3.insert(fd3.end(), d3[gg].begin(), d3[gg].begin() + points[gg]);
			fw.insert(fw.end(), weights[gg].begin(), weights[gg].begin() + points[gg]);
		}
		itensor_gpu_layout L;
		L.nmo = model_data.nmo; L.packed = static_cast<int>(I_tens.i_compact_); L.n_grids = n_atom_grids;
		L.n_blocks = static_cast<int>(bg.size());
		L.blk_grid = bg.data(); L.blk_point_start = bps.data(); L.blk_point_count = bpc.data();
		L.blk_n_active = bna.data(); L.blk_ao_off = bao.data(); L.blk_aos_off = baos.data();
		L.ao_all = ao_all.data(); L.ao_all_len = static_cast<long long>(ao_all.size());
		L.aos_all = aos_all.data(); L.aos_all_len = static_cast<long long>(aos_all.size());
		L.compact = compact_flat.data(); L.grid_point_off = goff.data();
		L.d1 = fd1.data(); L.d2 = fd2.data(); L.d3 = fd3.data(); L.weights = fw.data();
		L.n_points = static_cast<long long>(fd1.size());
		//-gpu_fp64 raises the whole device path to double. It is worth asking for on a card
		//with real double-precision units and expensive on one without, which is why it is
		//asked for rather than detected.
		const sf_precision iprec = opt->gpu_fp64 ? sf_precision::FP64 : sf_precision::FP32;
		const auto init_start = std::chrono::high_resolution_clock::now();
		itensor_on_gpu = itensor_gpu_init(L, iprec, opt->gpu_itensor_tensor);
		if (itensor_on_gpu && throughput::enabled())
			std::fprintf(stderr, "I tensor GPU: %.3f s upload and plan\n",
				std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - init_start).count());
		//What the device path actually issues, the blocks padded to their batch shapes,
		//real and imaginary halves together. Counted the way the path runs, or the GFLOP/s
		//row is fiction.
		if (itensor_on_gpu)
			itensor_gpu_dense_flops = itensor_gpu_issued_flops()
				* static_cast<double>(model_data.nr) * static_cast<double>(num_syms);
		//Say which processor produced the numbers, and which GEMM: the three do not agree
		//in the last digits, so a log that does not name one cannot be compared with
		//another. Gated like the other timing lines so the golden-file tests, which run
		//with no_date, keep their reference output.
		//
		//stderr, not cout, for the reason the shape diagnostic in itensor_gpu.cu gives:
		//cout is redirected into the log and moved again later in the run, so anything
		//written here never reached either the terminal or the file.
		if (!(opt->no_date)) {
			std::cerr << "GPU in use: XCW I tensor on ";
			if (itensor_on_gpu)
				std::cerr << "the device (" << (opt->gpu_fp64 ? "double" : "single")
				<< "-precision " << itensor_gpu_gemm_name() << " GEMM)"
				<< (opt->itensor_hybrid ? " with the CPU threads taking reflections alongside" : "");
			else
				std::cerr << "the CPU - device unavailable or problem too large";
			std::cerr << std::endl;
		}
	}
#endif
	//Held in the precision it is built in unless the settings say otherwise: single when
	//any path that contributes runs single. nr_small * i_compact_ deliberately in size_t,
	//the product passes 2^31 at nmo = 500 with 20k reflections.
	{
		const char* f = std::getenv("NOSPHERA2_XCW_I_FLOAT"); // Flawfinder: ignore
		const bool single_build = (itensor_on_gpu && !opt->gpu_fp64)
			|| ((!itensor_on_gpu || opt->itensor_hybrid) && opt->cpu_itensor_fp32);
		I_tens.i_float_ = I_tens.tensor_single || (f && std::atoi(f) != 0) || (single_build && !I_tens.tensor_double);
		if (single_build && !I_tens.i_float_)
			std::cout << "NOTE: the I tensor is built in single precision and held in double as i_double asks" << std::endl;
		if (!single_build && I_tens.i_float_)
			std::cout << "NOTE: the I tensor is built in double precision and narrowed to single as i_float asks" << std::endl;
	}
	decide_i_storage();
	if (!I_tens.i_streamed_) {
		if (I_tens.i_float_)
			I_tens.I32.assign(static_cast<size_t>(model_data.nr) * I_tens.i_compact_, std::complex<float>{});
		else
			I_tens.I.assign(static_cast<size_t>(model_data.nr) * I_tens.i_compact_, cdouble{});
	}
	//After the storage line: the bar owns the console from here until the last reflection
	if (!(opt->no_date)) {
		pb.emplace((unsigned long long)model_data.nr, 60, "=", "|", "Calculating XCW integrals...", std::cout);
	}
#if defined(NOSPHERA2_USE_GPU) || defined(NOSPHERA2_USE_METAL)
	if (itensor_on_gpu) gpu_thread = std::thread([&]() {
		const auto gpu_start = std::chrono::high_resolution_clock::now();
		//The device takes reflections in batches; the CPU threads keep taking them one at
		//a time from the same counter
		const int gpu_batch = itensor_gpu_batch(static_cast<int>(num_syms));
		cvec blk_gpu;
		if (I_tens.i_streamed_ || I_tens.i_float_) blk_gpu.assign(static_cast<size_t>(gpu_batch) * I_tens.i_compact_, cdouble{});
		vec kxs(static_cast<size_t>(gpu_batch) * num_syms), kys(kxs.size()), kzs(kxs.size());
		cvec facs(kxs.size() * n_atom_grids);
		int done = 0;
		auto collect_gpu = [&](const int rr0, const int n, const int slot) {
			if (I_tens.i_streamed_ || I_tens.i_float_) std::fill(blk_gpu.begin(), blk_gpu.end(), cdouble{});
			cdouble* const I_rr = (I_tens.i_streamed_ || I_tens.i_float_) ? blk_gpu.data()
									  : I_tens.I.data() + static_cast<size_t>(rr0) * I_tens.i_compact_;
			if (!itensor_gpu_collect(slot, I_rr, static_cast<long long>(I_tens.i_compact_)))
				err_checkf(false, "I tensor GPU read-back failed", std::cout);
			for (int r = 0; r < n && (I_tens.i_streamed_ || I_tens.i_float_); r++) {
				const cdouble* const row = blk_gpu.data() + static_cast<size_t>(r) * I_tens.i_compact_;
				if (I_tens.i_streamed_) {
					std::lock_guard<std::mutex> lock(i_write_mutex);
					I_tens.i_file_.write_block(rr0 + r, row);
				}
				else if (I_tens.i_float_) {
					std::complex<float>* const dst = I_tens.I32.data() + static_cast<size_t>(rr0 + r) * I_tens.i_compact_;
					for (size_t i = 0; i < I_tens.i_compact_; i++)
						dst[i] = std::complex<float>(static_cast<float>(row[i].real()),
													 static_cast<float>(row[i].imag()));
				}
			}
			done += n;
			gpu_ns_per_refl = std::chrono::duration_cast<std::chrono::nanoseconds>(
				std::chrono::high_resolution_clock::now() - gpu_start).count() / done;
			if (!(opt->no_date) && pb) for (int r = 0; r < n; r++) pb->update();
		};
		//Two result slots: a batch is collected after the next one has been submitted, so
		//the read-back overlaps that calculation
		int prev = -1, prev_n = 0, slot = 0;
		for (;;) {
			const int rr0 = next_refl.fetch_add(gpu_batch);
			if (rr0 >= model_data.nr) break;
			const int n = std::min(gpu_batch, model_data.nr - rr0);
			for (int r = 0; r < n; r++)
				for (int sy = 0; sy < static_cast<int>(num_syms); sy++) {
					const int rr = rr0 + r, c = r * static_cast<int>(num_syms) + sy;
					kxs[c] = k_pt[0][asym_lookup[rr][sy]];
					kys[c] = k_pt[1][asym_lookup[rr][sy]];
					kzs[c] = k_pt[2][asym_lookup[rr][sy]];
					for (int gg = 0; gg < n_atom_grids; gg++)
						facs[static_cast<size_t>(c) * n_atom_grids + gg] =
						asym_atoms[gg].asym_fact * DW_facts[gg][asym_lookup[rr][sy]]
						* phase_facts[gg][asym_lookup[rr][sy]] * translation_phase_facts[rr][sy];
				}
			if (!itensor_gpu_submit(slot, n, static_cast<int>(num_syms), kxs.data(), kys.data(), kzs.data(), facs.data()))
				err_checkf(false, "I tensor GPU evaluation failed", std::cout);
			if (prev >= 0) collect_gpu(prev, prev_n, slot ^ 1);
			prev = rr0; prev_n = n;
			slot ^= 1;
		}
		if (prev >= 0) collect_gpu(prev, prev_n, slot ^ 1);
		if (throughput::enabled())
			std::fprintf(stderr, "I tensor GPU: %.3f s wall time, %d of %d reflections\n",
				std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - gpu_start).count(),
				done, model_data.nr);
		itensor_gpu_free();
		//No bookkeeping here: the loop above runs for both paths and eval_I multiplies the
		//total by nr_small on the way out, so anything added here counts twice.
	});
#endif
	//Counted serially from the block structure both paths walk, so the CPU and GPU rows are
	//the same work measured two ways and no counter is touched by two threads.
	double itensor_flops = 0.0;
	double itensor_unscreened_tile_flops = 0.0;
	for (int g = 0; g < n_atom_grids; g++)
		for (const GridBlock& block : grid_blocks[g]) {
			for (const MatrixTile& tile : block.matrix_tiles)
				//real and imaginary passes, hence the factor of two
				itensor_flops += 2.0 * throughput::flops_gemm(tile.row_count, tile.col_count,
					block.point_count);
			const int n_active = static_cast<int>(block.active_aos.size());
			for (int row = 0; row < n_active; row += screened_tile_size) {
				const int row_count = std::min(screened_tile_size, n_active - row);
				for (int col = row; col < n_active; col += screened_tile_size)
					itensor_unscreened_tile_flops += 2.0 * throughput::flops_gemm(row_count,
						std::min(screened_tile_size, n_active - col), block.point_count);
			}
		}
	itensor_flops *= static_cast<double>(model_data.nr) * static_cast<double>(num_syms);
	itensor_unscreened_tile_flops *= static_cast<double>(model_data.nr) * static_cast<double>(num_syms);
	if (itensor_on_gpu && throughput::enabled()) {
		const double screened = itensor_unscreened_tile_flops > 0.0
			? 100.0 * (1.0 - itensor_flops / itensor_unscreened_tile_flops) : 0.0;
		std::fprintf(stderr, "I tensor GPU: %.3f dense GEMM GFLOP; CPU tiles %.3f unscreened, %.3f after overlap screening (%.1f%% pruned)\n",
			itensor_gpu_dense_flops / 1.0e9, itensor_unscreened_tile_flops / 1.0e9,
			itensor_flops / 1.0e9, screened);
	}

	if (!itensor_on_gpu || opt->itensor_hybrid)
	{
		//One core stays free for the thread feeding the device, which must not queue
		//behind a tile GEMM to submit the next reflection
		const int cpu_threads = itensor_on_gpu ? std::max(1, omp_get_max_threads() - 1) : omp_get_max_threads();
#pragma omp parallel num_threads(cpu_threads) reduction(+:skipped_grids)
		{
			vec2 single_k_pts(num_syms, vec(3));
			vec phase_angles;
			vec phase_sines;
			vec phase_cosines;
			vec w, c;
			std::vector<float> wf, cf;
			//One reflection's block while streaming or holding the tensor in single. No ordering is needed on
			//the way out: the file is reflection-major and the writer seeks to r's offset
			cvec blk;
			if (I_tens.i_streamed_ || I_tens.i_float_) blk.assign(I_tens.i_compact_, cdouble{});

#if !defined(__APPLE__)
			mkl_set_num_threads_local(1);
#endif
#if defined(__SSE2__) || defined(_M_X64)
			//Single precision underflows into subnormals on this data - AO tails of 1e-40 are
			//ordinary here - and an x86 core handles those a hundred times slower. Flush
			//them, as the device does; the double path never gets near 1e-308.
			const unsigned int csr_before = _mm_getcsr();
			if (opt->cpu_itensor_fp32) _mm_setcsr(csr_before | 0x8040);
#endif

			size_t max_points = 0;
			for (int g = 0; g < n_grids; g++) {
				max_points = std::max(max_points, static_cast<size_t>(points[g]));
			}
			cvec phase_buffer(max_points);
			phase_angles.resize(max_points);
			phase_sines.resize(max_points);
			phase_cosines.resize(max_points);
			size_t max_ao_block_size = 0, max_tile_result_size = 0;
			for (int g = 0; g < n_grids; g++) {
				for (const GridBlock& block : grid_blocks[g]) {
					max_ao_block_size = std::max(max_ao_block_size, block.ao_values.size());
					for (const MatrixTile& tile : block.matrix_tiles) {
						max_tile_result_size = std::max(max_tile_result_size, tile.result_offset + static_cast<size_t>(tile.row_count) * tile.col_count);
					}
				}
			}
			if (opt->cpu_itensor_fp32) {
				wf.resize(2 * max_ao_block_size);
				cf.resize(2 * max_tile_result_size);
			}
			else {
				w.resize(2 * max_ao_block_size);
				c.resize(2 * max_tile_result_size);
			}

			long long my_ns = 0;
			int my_done = 0;
			for (;;) {
				int r = next_refl.load();
				if (r >= model_data.nr) break;
				if (itensor_on_gpu && my_done > 0 && gpu_ns_per_refl.load() > 0 &&
					static_cast<long long>(model_data.nr - r) * gpu_ns_per_refl.load() < my_ns / my_done) break;
				if (!next_refl.compare_exchange_strong(r, r + 1)) continue;
				const auto r_start = std::chrono::high_resolution_clock::now();
				if (I_tens.i_streamed_ || I_tens.i_float_) std::fill(blk.begin(), blk.end(), cdouble{});
				cdouble* const I_r = (I_tens.i_streamed_ || I_tens.i_float_) ? blk.data()
					: I_tens.I.data() + static_cast<size_t>(r) * I_tens.i_compact_;
				const int* asym_lookup_r = asym_lookup[r].data();
				// Precompute weighted phase factors for integration
				for (int syms = 0; syms < num_syms; syms++) {
					single_k_pts[syms] = { k_pt[0][asym_lookup_r[syms]], k_pt[1][asym_lookup_r[syms]], k_pt[2][asym_lookup_r[syms]] };
					const int idx = asym_lookup_r[syms];
					for (int g = 0; g < n_grids; g++) {
						const int np_g = points[g];
						double* const angles = phase_angles.data();
						double* const sines = phase_sines.data();
						double* const cosines = phase_cosines.data();
						for (int p = 0; p < points[g]; p++) {
							angles[p] = single_k_pts[syms][0] * d1[g][p] + single_k_pts[syms][1] * d2[g][p] + single_k_pts[syms][2] * d3[g][p];
						}
#if defined(__APPLE__)
						for (int p = 0; p < points[g]; p++) {
							__sincos(angles[p], &sines[p], &cosines[p]);
						}
#else
						vdSinCos(np_g, angles, sines, cosines);
#endif
						cdouble* const phase_g = phase_buffer.data();
						const double* w_g = weights[g].data();
						for (int p = 0; p < np_g; p++) {
							phase_g[p] = cdouble(w_g[p] * cosines[p], w_g[p] * sines[p]);
						}
						const double asym_fact = asym_atoms[g].asym_fact;
						const cdouble* DW_fact_g = DW_facts[g].data();
						const double DW_im = DW_fact_g[idx].imag();
						const cdouble* phase_fact_g = phase_facts[g].data();
						const double phase_im = phase_fact_g[idx].imag();
						const double DW_re = DW_fact_g[idx].real();
						const double phase_re = phase_fact_g[idx].real();
						// Precompute basis function independent factors
						const cdouble grid_factor = cdouble(asym_fact * (DW_re * phase_re - DW_im * phase_im),
							asym_fact * (DW_re * phase_im + DW_im * phase_re));
						const cdouble factor = grid_factor * translation_phase_facts[r][syms];
						// This is where the magic happens
						for (const GridBlock& block : grid_blocks[g]) {
							const ivec& active_aos = block.active_aos;
							const int n_active = static_cast<int>(active_aos.size());
							const int np = block.point_count;
							const cdouble* phase_values = phase_g + block.point_start;
							//Weighted rows interleaved, 2j real and 2j + 1 imaginary of AO j, so a
							//tile's B rows are contiguous and one GEMM of twice the width returns
							//both parts; the result then alternates real, imaginary per column.
							auto tiles = [&](const auto* ao, auto* wv, auto* cv) {
								for (int local_mu = 0; local_mu < n_active; local_mu++) {
									const auto* ao_row = ao + static_cast<size_t>(local_mu) * np;
									auto* rw = wv + static_cast<size_t>(2 * local_mu) * np;
									auto* iw = rw + np;
									for (int p = 0; p < np; p++) {
										rw[p] = ao_row[p] * phase_values[p].real();
										iw[p] = ao_row[p] * phase_values[p].imag();
									}
								}
								for (const MatrixTile& tile : block.matrix_tiles)
									tile_gemm(tile.row_count, 2 * tile.col_count, np, ao + static_cast<size_t>(tile.row_start) * np,
										wv + static_cast<size_t>(2 * tile.col_start) * np, cv + 2 * tile.result_offset);
								for (const MatrixTile& tile : block.matrix_tiles) {
									for (int tile_row = 0; tile_row < tile.row_count; tile_row++) {
										const int mu = active_aos[tile.row_start + tile_row];
										const int first_tile_col = (tile.row_start == tile.col_start) ? tile_row : 0;
										const auto* crow = cv + 2 * tile.result_offset + static_cast<size_t>(tile_row) * 2 * tile.col_count;
										for (int tile_col = first_tile_col; tile_col < tile.col_count; tile_col++) {
											const int nu = active_aos[tile.col_start + tile_col];
											const int t = tri_compact[tri_index(mu, nu)];
											if (t >= 0)
												I_r[t] += cdouble(crow[2 * tile_col], crow[2 * tile_col + 1]) * factor;
										}
									}
								}
							};
							if (opt->cpu_itensor_fp32) tiles(block.ao_values_f.data(), wf.data(), cf.data());
							else tiles(block.ao_values.data(), w.data(), c.data());
						}
					}
				}
				if (I_tens.i_streamed_) {
					//A write is packed_size * 16 bytes against a whole reflection's worth
					//of integration, so the lock is not on the hot path
					std::lock_guard<std::mutex> lock(i_write_mutex);
					I_tens.i_file_.write_block(r, blk.data());
				}
				else if (I_tens.i_float_) {
					std::complex<float>* const dst = I_tens.I32.data() + static_cast<size_t>(r) * I_tens.i_compact_;
					for (size_t i = 0; i < I_tens.i_compact_; i++)
						dst[i] = std::complex<float>(static_cast<float>(blk[i].real()), static_cast<float>(blk[i].imag()));
				}
				my_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - r_start).count();
				my_done++;
				if (!(opt->no_date) && pb) {
					pb->update();
				}
			}
#if defined(__SSE2__) || defined(_M_X64)
			_mm_setcsr(csr_before);
#endif

		}
	}
	if (gpu_thread.joinable()) gpu_thread.join();
	if (I_tens.i_streamed_) {
		I_tens.i_file_.finish_write();
		open_i_stream_for_reading(I_tens.i_file_.path());
	}
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = end - start;
	skipped_grids_ = skipped_grids * model_data.nr;
	time_taken = std::chrono::duration<double>(duration).count();
	throughput::record("XCW I tensor", itensor_on_gpu, itensor_flops,
		1.0e3 * std::chrono::duration<double>(duration).count());
}

//F_r = Sum_at asym_fact Sum_s f_at(h R_s) exp(i h R_s x_at) T_at(h R_s) exp(2 pi i h t_s) plus the
//anomalous dispersion correction, with the atomic scattering factors laid out like DW_facts
void structure_factors::calc_F_calc(const cvec2& scattering_factors) {
#pragma omp parallel for
	for (int r = 0; r < model_data.nr; r++) {
		const ivec lookup = generate_asym_lookup(r);
		cdouble sum = scatter_data.anom_correction[r];
		for (int at = 0; at < model_data.ncen; at++) {
			cdouble temp1 = 0;
			for (int r_asym = 0; r_asym < lookup.size(); r_asym++) {
				temp1 += scattering_factors[at][lookup[r_asym]] * phase_facts[at][lookup[r_asym]] * DW_facts[at][lookup[r_asym]] * translation_phase_facts[r][r_asym];
			}
			sum += temp1 * asym_atoms[at].asym_fact;
		}
		scatter_data.F_calc[r] = sum;
	}
}

void structure_factors::calc_F_calc(I_tensor& I, const dMatrix2& D) {
	// Density matrix from occ is half of what I need, so times 2 and times (2x2)=4
	//Streamed or resident the walk is the same; the outer loop is one window when
	//the tensor is resident
	const int step = std::max(1, I.i_streamed_ ? I.i_window_ : model_data.nr);
	//The parallel region wraps the window loop: entering one per window would pay
	//team startup and a barrier for a few reflections of work. omp single does the
	//read, and its implicit barrier stops a thread entering a window not yet loaded.
	//load() can throw and an exception must not leave an OpenMP structured block,
	//so it is recorded and rethrown after the region.
	std::string io_error;
#if defined(NOSPHERA2_USE_GPU) || defined(NOSPHERA2_USE_METAL)
	if (I.i_on_device_) {
		vec w(I.i_compact_);
		for (size_t k = 0; k < I.i_compact_; k++)
			w[k] = (I.i_pair_mu_[k] == I.i_pair_nu_[k] ? 2.0 : 4.0) * D(I.i_pair_mu_[k], I.i_pair_nu_[k]);
		err_checkf(itensor_gpu_rows(w.data(), scatter_data.anom_correction.data(), scatter_data.F_calc.data()), "I tensor walk on the device failed", std::cout);
		return;
	}
#endif
#pragma omp parallel
	{
		for (int r0 = 0; r0 < model_data.nr; r0 += step) {
			const int r1 = std::min(r0 + step, model_data.nr);
			if (I.i_streamed_) {
#pragma omp single
				{
					try { I.i_file_.load(r0, r1); }
					catch (const std::exception& e) { io_error = e.what(); }
				}
			}
#pragma omp for schedule(static)
			for (int r = r0; r < r1; ++r) {
				if (!io_error.empty()) continue;
				//One walk, either element type: the accumulation stays in double whatever
				//the tensor is stored as, so float storage costs precision in the stored
				//value and nothing in the sum.
				auto accumulate = [&](const auto* I_r) {
					cdouble sum = scatter_data.anom_correction[r];
					const int* pmu = I.i_pair_mu_.data();
					const int* pnu = I.i_pair_nu_.data();
					for (size_t k = 0; k < I.i_compact_; k++)
						sum += (pmu[k] == pnu[k] ? 2.0 : 4.0) * cdouble(I_r[k]) * D(pmu[k], pnu[k]);
					scatter_data.F_calc[r] = sum;
				};
				if (I.i_float_) accumulate(I.i_block32(r)); else accumulate(I.i_block(r));
			}
		}
	}
	if (!io_error.empty()) throw std::runtime_error(io_error);
}
