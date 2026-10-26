#include "pch.h"
#include "tuning.h"
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
	orbital_basis = opt_in.xcw_settings.basis_overrides.empty() ? BasisSetLibrary::get_basis_set(opt_in.xcw_settings.basis_set_name) :
		BasisSetLibrary::get_basis_set_with_overrides(opt_in.xcw_settings.basis_set_name, opt_in.xcw_settings.basis_overrides);
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
	model_data.nr = scatter_data.hkl.size();
	// Set F_calc sizes
	scatter_data.F_calc.resize(model_data.nr, 0);
	scatter_data.anom_correction.resize(model_data.nr, 0);
	// Use a mask that is directly saved in the anomalous dispersion correction, that way is most computationally efficient
	if(opt_in.xcw_settings.mask) {
		err_checkf(!opt_in.mask_file.empty(), "Please provide an fab file with the `-mask` option or remove the mask keyword from the XCW settings file.", std::cout);
		// Read out the mask corrections
		std::filesystem::path mask_filename = opt_in.mask_file;
		std::cout << "Reading mask from " << mask_filename << std::endl;
		vec2 mask_info;
		read_mask(mask_filename, scatter_data.hkl, scatter_data.hkl_enlarged, unit_cell, std::cout, mask_info, model_data.nr, opt_in.debug);
		for (int r = 0; r < model_data.nr; r++) {
			scatter_data.anom_correction[r] = cdouble(mask_info[r][0], mask_info[r][1]);
		}
	}
	// Directly convert into reciprocal space so rotation for grown structures can be done directly
	U_cif2U_star();
	model_data.wavelength = opt_in.xcw_settings.wavelength > 0.0 ? opt_in.xcw_settings.wavelength : wavelength_;
	err_checkf(model_data.ncen > 0, "No atoms were read from " + cif.string() + "! Is there an _atom_site loop with labels, type symbols and fractional coordinates?", std::cout);

	// First step of the symmetry resolution: Calculate the crystallographic multiplicities for the asymmetric unit
	unit_cell.asym_unit_cryst_mult(asym_atoms);

	{
		//ivec3 symmetry_linking_list;
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
			// Second step of the symmetry resolution: Grow the asymmetric unit (this will now grow from the complete xyz file since we hope that we can deal with improperly grown structures)
			// No check for symmetry equivalency is done here anymore
			unit_cell.grow_asym_atoms(asym_atoms, xyz_atoms);
			// Third step of the symmetry resolution: Modify the symmetry factors for the grown atoms, leave ungrown atoms untouched, set the symmetry operation that generated the grown atom from its asymmetric parent (grown_by and grown_from)
			unit_cell.resolve_internal_symmetry(asym_atoms);
		}
		model_data.ncen = asym_atoms.size();
		if (opt_in.xcw_settings.grown) {
			// Fourth step of the symmetry resolution: Copy U_iso, the dispersion and the ADPs of each grown atom's parent, the ADPs rotated onto the symmetry generated atoms
			unit_cell.grow_ADPs(asym_atoms, ADPs);
			// Fifth step of the symmetry resolution: Delete the fully applied operations and solve the symmetry factors for the operations that are kept
			applied_symmetry = unit_cell.apply_grown(scatter_data.hkl, scatter_data.hkl_enlarged, asym_atoms);
		}
		else {
			// Symmetry factors of an ungrown model, also for an atom on a special position
			unit_cell.set_sym_fact_from_mult(asym_atoms);
		}
		if (tuning("NOSPHERA2_DEBUG_ASYMFACT")) {
			std::cerr << "applied_symmetry (deleted):";
			for (int s : applied_symmetry) std::cerr << " " << s;
			std::cerr << std::endl << "surviving sym ops: " << unit_cell.get_trans()[0].size() << std::endl;
			for (size_t i = 0; i < asym_atoms.size(); i++)
				std::cerr << i << " grown=" << asym_atoms[i].grown << " sym_op=" << asym_atoms[i].grown_by[0]
				<< " sym_fact=" << asym_atoms[i].sym_fact << std::endl;
			std::cerr << "hkl_enlarged size: " << scatter_data.hkl_enlarged.size() << std::endl;
		}
	}
	// occ only finds this inside the first SCF, and its message names neither the atom count nor the charge
	int n_electrons = -opt_in.xcw_settings.charge;
	for (int at = 0; at < model_data.ncen; at++) n_electrons += asym_atoms[at].type;
	err_checkf(n_electrons > 0 && (n_electrons + opt_in.xcw_settings.multiplicity - 1) % 2 == 0,
		"The XCW model has " + std::to_string(model_data.ncen) + " atoms and " + std::to_string(n_electrons) +
		" electrons at charge " + std::to_string(opt_in.xcw_settings.charge) + ", which cannot have multiplicity " +
		std::to_string(opt_in.xcw_settings.multiplicity) + ". Check the charge and multiplicity, missing H atoms and disorder parts; a molecule on a special position needs the grown model (keyword grown and -xyz).",
		std::cout);

	// Convert the ADPs from reciprocal space to Cartesian coordinates
	U_star2U_cart();

	// Generate k_pts and set the number of reflections
	make_k_pts(model_data.nr_enlarged != 0 && scatter_data.hkl.size() == 0, opt_in.save_k_pts, unit_cell, scatter_data.hkl_enlarged, k_pt, std::cout, opt_in.debug);
	model_data.nr_enlarged = scatter_data.hkl_enlarged.size();

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
		ivec asym_list = unit_cell.generate_asym_lookup(r, scatter_data.hkl, scatter_data.hkl_enlarged);
		vec q_temp = { k_pt[0][asym_list[0]], k_pt[1][asym_list[0]], k_pt[2][asym_list[0]] };
		std::transform(q_temp.begin(), q_temp.end(), q_temp.begin(), [angstrom2bohr](double x) { return x * angstrom2bohr; });
		for (int sym_op = 0; sym_op < trans[0].size(); sym_op++) {
			vec trans_temp = { trans[0][sym_op], trans[1][sym_op], trans[2][sym_op] };
			trans_temp = dot(cm, trans_temp, true);
			cdouble exponent(0, dot_BLAS(q_temp, trans_temp, false));
			translation_phase_facts[r][sym_op] = std::exp(exponent);
		}
	}
	phases_set_ = true;
	// closing function
}

void structure_factors::set_anom() {
	if(anom_set_) return;
	int r, at, sym_op;
	ivec2 asym_lookup(model_data.nr);
	for (r = 0; r < model_data.nr; r++) {
		asym_lookup[r] = unit_cell.generate_asym_lookup(r, scatter_data.hkl, scatter_data.hkl_enlarged);
	}
	for (r = 0; r < model_data.nr; r++) {
		//lookup[sym_op]: the hkl_enlarged entry of reflection r under symmetry operation sym_op
		const ivec& lookup = asym_lookup[r];
		cdouble sum = 0;
		for (at = 0; at < model_data.ncen; at++) {
			cdouble temp1 = 0;
			for (sym_op = 0; sym_op < lookup.size(); sym_op++) {
				temp1 += phase_facts[at][lookup[sym_op]] * DW_facts[at][lookup[sym_op]] * translation_phase_facts[r][sym_op];
			}
			sum += temp1 * asym_atoms[at].sym_fact * asym_atoms[at].anom;
		}
		scatter_data.anom_correction[r] += sum;
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
//scale has to be stationary for the criterion, or the SCF descends a different functional
//than the one it prints; with weights, an unweighted k is not.
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
				+ I_tens.i_tensor_save_path.string() + ": a second copy costs as much memory as the first. Remove `safe` from the XCW settings.");
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
		const char* f = tuning("NOSPHERA2_XCW_I_FLOAT");
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
	if (I_tens.i_tensor_save_path.empty() || I_tens.i_streamed_ || XG.mode == 1) return;
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
	if (budget == 0) {
		const size_t avail = free_memory();
		if (avail > 0) {
			//Four fifths: the SCF matrices, the grids and OCC's own allocations live in the
			//rest, and a tensor that only just fits would page rather than run. What the
			//process can have is a platform question - a cgroup here, a job object on
			//Windows, page classes on a Mac - and lives in convenience.cpp.
			//At least 1, since 0 means "no budget, hold it all".
			budget = std::max<size_t>(avail / 5 * 4, 1);
			automatic = !(opt->mem_given && opt->mem > 0.0);
			source = automatic ? "four fifths of the memory this job can have" : "-mem";
		}
	}
	return budget;
}

//What the next large allocation may take: the platform's and scheduler's headroom and, with -mem,
//-mem minus what this process already holds. The density-fitting integrals, the I tensor and the
//stored ERIs are sized in that order from it, so together they stay under -mem (they used to take
//half the node, all of -mem and four fifths of the node, and a SLURM job went OOM). OCC fills
//the density-fitting store only at the first Fock build, after all three decisions, so its
//size is held back here by hand.
size_t structure_factors::free_memory() const {
	size_t avail = available_memory_bytes();
	if (opt->mem_given && opt->mem > 0.0) {
		const size_t cap = static_cast<size_t>(opt->mem * 1024.0 * 1024.0), held = resident_memory_bytes();
		const size_t left = cap > held ? cap - held : 1;
		if (avail == 0 || left < avail) avail = left;
	}
	if (avail != 0) avail = avail > df_reserved_ ? avail - df_reserved_ : 1;
	return avail;
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
	//Streamed straight into the `safe` file when there is one, so saving it costs nothing
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
	// Construct a dummy wavefunction from the asym atoms (needed for grids)
	WFN dummy_wave;
	for (int at = 0; at < model_data.ncen; at++) {
		atom temp_atom;
		temp_atom.set_coordinate(0, asym_atoms[at].pos[0]);
		temp_atom.set_coordinate(1, asym_atoms[at].pos[1]);
		temp_atom.set_coordinate(2, asym_atoms[at].pos[2]);
		temp_atom.set_charge(asym_atoms[at].type);
		dummy_wave.push_back_atom(temp_atom);
	}
	// Load the basis set into the dummy wavefunction
	load_basis_into_WFN(dummy_wave, orbital_basis, false, true);
	dummy_wave.delete_unoccupied_MOs();
	bvec needs_grid(model_data.ncen, true);
	ivec asym_atom_list(model_data.ncen);
	for (int at = 0; at < model_data.ncen; at++) {
		asym_atom_list[at] = at;
	}
	grid_manager.setup3DGridsForMolecule(dummy_wave, asym_atom_list, needs_grid, unit_cell);

	// Evaluation of density on the grids for grid pruning
	GridData& GD = grid_manager.getGridData();
	vec2* grids = grid_manager.getNeedsHelper() ? GD.helper_grids.data() : GD.atomic_grids.data();
	vec2 d1, d2, d3, weights;
	const int n_grids = grid_manager.getNeedsHelper() ? GD.helper_grids.size() : GD.atomic_grids.size();
	for (int g = 0; g < n_grids; g++) {
		std::fill(grids[g][GridData::GridIndex::WFN_DENSITY].begin(), grids[g][GridData::GridIndex::WFN_DENSITY].end(), 1.0);
	}
	//where each point sits in its prototype grid, for the grid route's Fourier-Bessel transform
	sf_shell_layout layout;
	grid_manager.getDensityVectors(dummy_wave, asym_atom_list, d1, d2, d3, weights, false, &layout);
	const int* points = grid_manager.getNeedsHelper() ? GD.helper_num_points_per_atom.data() : GD.num_points_per_atom.data();
	const int total_points = grid_manager.getTotalGridPoints();
	std::cout << "Total number of grid points after pruning: " << total_points << std::endl;

	// Precompute the lookup table linking each unique reflection to its symmetry equivalents
	ivec2 asym_lookup(model_data.nr);
	for (r = 0; r < model_data.nr; r++) {
		asym_lookup[r] = unit_cell.generate_asym_lookup(r, scatter_data.hkl, scatter_data.hkl_enlarged);
	}
	const unsigned int num_syms = unit_cell.get_num_syms();
	// The AO values and weights are real, so the grid integral sum_p w_p phi_mu phi_nu exp(i q d_p) of -q is the
	// conjugate of that of q: only the entries of inversion_cleaned_hkl are integrated
	const bool no_screening = tuning("NOSPHERA2_NO_INVERSION_SCREENING") != nullptr;
	ivec2 inversion_link;
	if (no_screening)
		for (int e = 0; e < model_data.nr_enlarged; e++) inversion_link.push_back({ e });
	else
		inversion_link = unit_cell.remove_inv_vectors(scatter_data.hkl_enlarged, scatter_data.inversion_cleaned_hkl);
	const bool inversion_screening = static_cast<int>(inversion_link.size()) != model_data.nr_enlarged;
	ivec clean_of(model_data.nr_enlarged);
	bvec is_inverse(model_data.nr_enlarged, false);
	for (int c = 0; c < inversion_link.size(); c++)
		for (int i = 0; i < inversion_link[c].size(); i++) {
			clean_of[inversion_link[c][i]] = c;
			is_inverse[inversion_link[c][i]] = i == 1;
		}
	// A reflection whose inverse is measured too has the conjugate of the inverse's row, the images of -h being minus those
	// of h with conjugate phases: only the lower of the two is integrated
	ivec row_partner(model_data.nr, -1);
	bvec row_copy(model_data.nr, false);
	if (!no_screening) {
		const std::vector<i3> hv(scatter_data.hkl.begin(), scatter_data.hkl.end());
		for (r = 0; r < model_data.nr; r++) {
			if (row_copy[r]) continue;
			const i3 minus_h = { -hv[r][0], -hv[r][1], -hv[r][2] };
			const int j = static_cast<int>(std::lower_bound(hv.begin(), hv.end(), minus_h) - hv.begin());
			if (j > r && j < model_data.nr && cell::check_inversion(hv[r], hv[j])) {
				row_partner[r] = j;
				row_copy[j] = true;
			}
		}
	}
	// The integrations of reflection r: with inversion screening one per cleaned entry its k-points land on, carrying
	// the summed translation phases of the operations landing on the entry and on its inverse; without, one per operation
	//c: the cleaned entry integrated; self / inv: whether an operation lands on it or on its inverse, t_*: their summed translation phases
	struct k_job { int c; bool self, inv; cdouble t_self, t_inv; };
	auto make_jobs = [&](const int r, std::vector<k_job>& jobs) {
		jobs.clear();
		if (row_copy[r]) return;
		for (int sym_op = 0; sym_op < static_cast<int>(num_syms); sym_op++) {
			const int e = asym_lookup[r][sym_op];
			int j = inversion_screening ? 0 : static_cast<int>(jobs.size());
			while (j < jobs.size() && jobs[j].c != clean_of[e]) j++;
			if (j == jobs.size()) jobs.push_back({ clean_of[e], false, false, 0.0, 0.0 });
			if (is_inverse[e]) {
				jobs[j].inv = true;
				jobs[j].t_inv += translation_phase_facts[r][sym_op];
			}
			else {
				jobs[j].self = true;
				jobs[j].t_self += translation_phase_facts[r][sym_op];
			}
		}
	};
	// The basis function independent factor of atom grid g at k-point e
	auto grid_factor = [&](const int g, const int e) {
		const double sym_fact = asym_atoms[g].sym_fact;
		const cdouble DW = DW_facts[g][e], ph = phase_facts[g][e];
		return cdouble(sym_fact * (DW.real() * ph.real() - DW.imag() * ph.imag()),
			sym_fact * (DW.real() * ph.imag() + DW.imag() * ph.real()));
	};
	long long n_integrations = 0;
	{
		std::vector<k_job> jobs;
		for (r = 0; r < model_data.nr; r++) {
			make_jobs(r, jobs);
			n_integrations += jobs.size();
		}
	}
	if (inversion_screening && !(opt->no_date)) {
		const size_t n_inv = scatter_data.hkl_enlarged.size() - scatter_data.inversion_cleaned_hkl.size();
		std::cout << std::fixed << std::setprecision(2) << "Inversion screening: Screened out " << n_inv << " of " << scatter_data.hkl_enlarged.size()
			<< " symmetry generated reflections (" << 100.0 * n_inv / scatter_data.hkl_enlarged.size() << "%)" << std::endl;
	}

	// Save the grid positions for each atom in the dummy wavefunction (used for screening)
	vec2 grid_positions(model_data.ncen);
	for (int at = 0; at < model_data.ncen; at++) {
		grid_positions[at] = { dummy_wave.get_atom_pos(at)[0], dummy_wave.get_atom_pos(at)[1], dummy_wave.get_atom_pos(at)[2] };
	}

	// AO pair screening: an AO pair mu, nu on different centres is skipped when the bound on its overlap, summed over
	// the primitives k of mu and l of nu, falls below e_tol
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
			vec N_k_cache(mu_primitives.size());
			for (int k = 0; k < mu_primitives.size(); k++) {
				const double k_exp = mu_primitives[k].get_exp();
				N_k_cache[k] = mu_type == 0 ? sqrt_inv_four_pi : sph_harmonic_max_mu * min_mu_type_half_exp * std::pow(mu_type / k_exp, mu_type_half);
			}
			for (nu = mu + 1; nu < model_data.nmo; nu++) {
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
				vec N_l_cache(nu_primitives.size());
				for (int l = 0; l < nu_primitives.size(); l++) {
					const double l_exp = nu_primitives[l].get_exp();
					N_l_cache[l] = nu_type == 0 ? sqrt_inv_four_pi : sph_harmonic_max_nu * min_nu_type_half_exp * std::pow(nu_type / l_exp, nu_type_half);
				}
				estimate = 0.0;
				for (int k = 0; k < mu_primitives.size(); k++) {
					const double k_coef = mu_primitives[k].get_coef();
					const double k_exp = mu_primitives[k].get_exp();
					const double N_k = N_k_cache[k];
					for (int l = 0; l < nu_primitives.size(); l++) {
						const double l_coef = nu_primitives[l].get_coef();
						const double l_exp = nu_primitives[l].get_exp();
						const double inv_exp_sum = 1.0 / (k_exp + l_exp);
						const double N_l = N_l_cache[l];
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
	//compact ball rather than a spherical shell.
	//Work is sum over blocks of n_active^2 * points, so this is quadratic in what it saves.
	//Sums over points are order independent, so the reordering changes no result.
	//Compaction and the AO threshold are worth nothing apart and a great deal together:
	//reordering alone is slower (the blocks shrink and n_active does not), the threshold
	//alone barely moves the work, and together they are faster with the GooF, energies and
	//convergence lines identical.
	//So only reorder when the threshold can actually prune
	const double ao_block_threshold = [&] {
		const char* e = tuning("NOSPHERA2_ITENSOR_AO_TOL");
		if (e) { const double v = std::atof(e); return v >= 0.0 ? v : 0.0; }
		return cutoff(opt->accuracy);
	}();
	const bool morton_applied = (tuning("NOSPHERA2_ITENSOR_NO_MORTON") == nullptr)
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
			if (g < static_cast<int>(layout.idx.size()) && static_cast<int>(layout.idx[g].size()) >= npts) {
				ivec t(npts);
				for (int p = 0; p < npts; p++) t[p] = layout.idx[g][perm[p]];
				layout.idx[g].swap(t);
			}
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
	if (tuning("NOSPHERA2_ITENSOR_AOSTATS")) {
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

	// Tile the grid points for each atom into blocks
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
		if (!tuning("NOSPHERA2_ITENSOR_SKIPSTATS")) return;
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
		//the cost.
		//NOSPHERA2_ITENSOR_CHUNK sets the target; 0 restores the three whole bands.
		const int chunk = [] {
			const char* e = tuning("NOSPHERA2_ITENSOR_CHUNK");
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
	if (!(opt->no_date) && ao_slots_carrying > 0) {
		const long long dropped = ao_slots_carrying - ao_slots_kept;
		std::cout << std::fixed << std::setprecision(2)
			<< "Screened out " << dropped << " of " << ao_slots_carrying
			<< " AO-block entries (" << 100.0 * static_cast<double>(dropped)
			/ static_cast<double>(ao_slots_carrying) << "%) below "
			<< std::scientific << std::setprecision(0) << ao_block_threshold
			<< std::fixed << std::setprecision(2) << " on their block\n";

		double work_done = 0.0, work_unscreened = 0.0;
		for (int g = 0; g < n_atom_grids; g++)
			for (const GridBlock& b : grid_blocks[g]) {
				const double na = static_cast<double>(b.active_aos.size());
				work_done += na * na * b.point_count;
				work_unscreened += static_cast<double>(model_data.nmo) * model_data.nmo * b.point_count;
			}
		//Unscreened is one integration for every entry of hkl_enlarged, screened only the integrations inversion screening leaves
		work_done *= static_cast<double>(n_integrations);
		work_unscreened *= static_cast<double>(model_data.nr_enlarged);
		if (work_done > 0.0)
			std::cout << std::fixed << std::setprecision(2)
				<< "I tensor work after all screening: " << (work_done / 1e12)
				<< " of " << (work_unscreened / 1e12) << " Tflop-equivalents ("
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
	//The grid route is the default; the tensor is built only for NOSPHERA2_XCW_GRID=0 or when it is to be saved
	const char* xg = tuning("NOSPHERA2_XCW_GRID");
	if (const int xg_mode = xg ? std::atoi(xg) : (I_tens.i_tensor_save_path.empty() ? 1 : 0); xg_mode != 0) {
		XG.mode = xg_mode;
		const bool keep = XG.mode == 2;
		auto take = [keep](auto& a) -> std::decay_t<decltype(a)> { if (keep) return a; return std::move(a); };
		XG.blocks = take(grid_blocks);
		XG.d1 = take(d1); XG.d2 = take(d2); XG.d3 = take(d3); XG.w = take(weights);
		XG.skip = take(skip); XG.lookup = take(asym_lookup); XG.inversion_link = take(inversion_link);
		XG.clean_of = take(clean_of); XG.is_inverse = take(is_inverse);
		plan_grid_fb(layout);
		std::cout << "XCW on the grid: the I tensor is " << (keep ? "built as well, for comparison" : "not built") << std::endl;
		if (!keep) {
			skipped_grids_ = skipped_grids * model_data.nr;
			time_taken = 0;
			return;
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
	//Row r of the tensor and, when the inverse of r is measured too, the conjugate row of that one, whichever way the tensor is held
	auto store_rows = [&](const int r, const cdouble* row) {
		const size_t n = I_tens.i_compact_;
		cvec flipped;
		for (int i = 0; i < 2; i++) {
			const int rr = i == 0 ? r : row_partner[r];
			if (rr < 0) break;
			const cdouble* src = row;
			if (i == 1) {
				flipped.resize(n);
				for (size_t t = 0; t < n; t++) flipped[t] = std::conj(row[t]);
				src = flipped.data();
			}
			const size_t at = static_cast<size_t>(rr) * n;
			if (I_tens.i_streamed_) {
				//A write is packed_size * 16 bytes against a whole reflection's worth
				//of integration, so the lock is not on the hot path
				std::lock_guard<std::mutex> lock(i_write_mutex);
				I_tens.i_file_.write_block(rr, src);
			}
			else if (I_tens.i_float_) {
				std::complex<float>* const dst = I_tens.I32.data() + at;
				for (size_t t = 0; t < n; t++) dst[t] = std::complex<float>(static_cast<float>(src[t].real()), static_cast<float>(src[t].imag()));
			}
			else if (src != I_tens.I.data() + at)
				std::copy(src, src + n, I_tens.I.data() + at);
		}
	};
#if defined(NOSPHERA2_USE_GPU) || defined(NOSPHERA2_USE_METAL)
	//GPU: sets the device up for the I tensor GEMMs; the CPU threads below integrate whatever it does not take.
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
			itensor_gpu_dense_flops = itensor_gpu_issued_flops() * static_cast<double>(n_integrations);
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
		const char* f = tuning("NOSPHERA2_XCW_I_FLOAT");
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
	//GPU: one host thread feeds the device batches of reflections next to the CPU workers
	if (itensor_on_gpu) gpu_thread = std::thread([&]() {
		const auto gpu_start = std::chrono::high_resolution_clock::now();
		//The device takes reflections in batches; the CPU threads keep taking them one at
		//a time from the same counter
		const int gpu_batch = itensor_gpu_batch(static_cast<int>(num_syms));
		cvec blk_gpu(static_cast<size_t>(gpu_batch) * I_tens.i_compact_, cdouble{});
		//A column per integration, at most num_syms per reflection
		vec kxs(static_cast<size_t>(gpu_batch) * num_syms), kys(kxs.size()), kzs(kxs.size());
		cvec facs(kxs.size() * n_atom_grids), facs_inv(facs.size());
		ivec col_off(gpu_batch + 1);
		std::vector<k_job> jobs_gpu;
		int done = 0;
		auto collect_gpu = [&](const int rr0, const int n, const int slot) {
			std::fill(blk_gpu.begin(), blk_gpu.end(), cdouble{});
			if (!itensor_gpu_collect(slot, blk_gpu.data(), static_cast<long long>(I_tens.i_compact_)))
				err_checkf(false, "I tensor GPU read-back failed", std::cout);
			//A row whose inverse is measured too had no columns and comes back empty; the row it is copied from writes it
			for (int r = 0; r < n; r++)
				if (!row_copy[rr0 + r]) store_rows(rr0 + r, blk_gpu.data() + static_cast<size_t>(r) * I_tens.i_compact_);
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
			//The integrations of the CPU path: at the cleaned entry, its inverse taking the conjugate
			int c = 0;
			for (int r = 0; r < n; r++) {
				col_off[r] = c;
				make_jobs(rr0 + r, jobs_gpu);
				for (int jb = 0; jb < static_cast<int>(jobs_gpu.size()); jb++, c++) {
					const k_job& job = jobs_gpu[jb];
					const int idx = inversion_link[job.c][0];
					const int inv_idx = job.inv ? inversion_link[job.c][1] : idx;
					kxs[c] = k_pt[0][idx];
					kys[c] = k_pt[1][idx];
					kzs[c] = k_pt[2][idx];
					for (int gg = 0; gg < n_atom_grids; gg++) {
						const size_t f = static_cast<size_t>(c) * n_atom_grids + gg;
						facs[f] = job.self ? grid_factor(gg, idx) * job.t_self : cdouble(0.0, 0.0);
						facs_inv[f] = job.inv ? grid_factor(gg, inv_idx) * job.t_inv : cdouble(0.0, 0.0);
					}
				}
			}
			col_off[n] = c;
			if (c == 0) {
				//Every row of the batch is the conjugate of one integrated elsewhere
				done += n;
				if (!(opt->no_date) && pb) for (int r = 0; r < n; r++) pb->update();
				continue;
			}
			if (!itensor_gpu_submit(slot, n, col_off.data(), kxs.data(), kys.data(), kzs.data(), facs.data(), facs_inv.data()))
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
	itensor_flops *= static_cast<double>(n_integrations);
	itensor_unscreened_tile_flops *= static_cast<double>(n_integrations);
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
#ifdef __ANDROID__
		//Once its thread count is set explicitly, OpenMP OpenBLAS threads calls made inside a
		//parallel region too, and they all queue for its one buffer slot spinning in sched_yield
		const int blas_threads = openblas_get_num_threads();
		openblas_set_num_threads(1);
#endif
#pragma omp parallel num_threads(cpu_threads) reduction(+:skipped_grids)
		{
			std::vector<k_job> jobs;
			vec phase_angles;
			vec phase_sines;
			vec phase_cosines;
			vec w, c;
			std::vector<float> wf, cf;
			//One reflection's block while streaming or holding the tensor in single. No ordering is needed on
			//the way out: the file is reflection-major and the writer seeks to r's offset
			cvec blk;
			if (I_tens.i_streamed_ || I_tens.i_float_) blk.assign(I_tens.i_compact_, cdouble{});

#if !defined(__APPLE__) && !defined(NSA2_OPENBLAS)
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
				if (row_copy[r]) {
					if (!(opt->no_date) && pb) pb->update();
					continue;
				}
				const auto r_start = std::chrono::high_resolution_clock::now();
				if (I_tens.i_streamed_ || I_tens.i_float_) std::fill(blk.begin(), blk.end(), cdouble{});
				cdouble* const I_r = (I_tens.i_streamed_ || I_tens.i_float_) ? blk.data()
					: I_tens.I.data() + static_cast<size_t>(r) * I_tens.i_compact_;
				make_jobs(r, jobs);
				// Precompute weighted phase factors for integration
				for (int jb = 0; jb < jobs.size(); jb++) {
					const k_job& job = jobs[jb];
					// Integrated at the cleaned entry, the operations landing on its inverse take the conjugate
					const int idx = inversion_link[job.c][0];
					const int inv_idx = job.inv ? inversion_link[job.c][1] : idx;
					const double kx = k_pt[0][idx], ky = k_pt[1][idx], kz = k_pt[2][idx];
					for (int g = 0; g < n_grids; g++) {
						const int np_g = points[g];
						double* const angles = phase_angles.data();
						double* const sines = phase_sines.data();
						double* const cosines = phase_cosines.data();
						for (int p = 0; p < points[g]; p++) {
							angles[p] = kx * d1[g][p] + ky * d2[g][p] + kz * d3[g][p];
						}
#if defined(__APPLE__)
						for (int p = 0; p < points[g]; p++) {
							__sincos(angles[p], &sines[p], &cosines[p]);
						}
#elif defined(NSA2_OPENBLAS)
						sincos_shared_n(np_g, angles, sines, cosines);
#else
						vdSinCos(np_g, angles, sines, cosines);
#endif
						cdouble* const phase_g = phase_buffer.data();
						const double* w_g = weights[g].data();
						for (int p = 0; p < np_g; p++) {
							phase_g[p] = cdouble(w_g[p] * cosines[p], w_g[p] * sines[p]);
						}
						// Precompute basis function independent factors
						const cdouble factor = grid_factor(g, idx) * job.t_self, factor_inv = grid_factor(g, inv_idx) * job.t_inv;
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
											if (t < 0) continue;
											const cdouble v(crow[2 * tile_col], crow[2 * tile_col + 1]);
											if (job.self) I_r[t] += v * factor;
											if (job.inv) I_r[t] += std::conj(v) * factor_inv;
										}
									}
								}
							};
							if (opt->cpu_itensor_fp32) tiles(block.ao_values_f.data(), wf.data(), cf.data());
							else tiles(block.ao_values.data(), w.data(), c.data());
						}
					}
				}
				store_rows(r, I_r);
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
#ifdef __ANDROID__
		openblas_set_num_threads(blas_threads);
#endif
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

//F_r = Sum_at sym_fact Sum_s f_at(h R_s) exp(i h R_s x_at) T_at(h R_s) exp(2 pi i h t_s) plus the
//anomalous dispersion correction, with the atomic scattering factors laid out like DW_facts
void structure_factors::calc_F_calc(const cvec2& scattering_factors) {
#pragma omp parallel for
	for (int r = 0; r < model_data.nr; r++) {
		const ivec lookup = unit_cell.generate_asym_lookup(r, scatter_data.hkl, scatter_data.hkl_enlarged);
		cdouble sum = scatter_data.anom_correction[r];
		for (int at = 0; at < model_data.ncen; at++) {
			cdouble temp1 = 0;
			for (int sym_op = 0; sym_op < lookup.size(); sym_op++) {
				temp1 += scattering_factors[at][lookup[sym_op]] * phase_facts[at][lookup[sym_op]] * DW_facts[at][lookup[sym_op]] * translation_phase_facts[r][sym_op];
			}
			sum += temp1 * asym_atoms[at].sym_fact;
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
	//GPU: the same walk over the device copy of the tensor, the pair weights folded into D
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

//F_calc of D through the I tensor or on the grid, see xcw_grid; mode 2 runs both and prints how far
//the grid lands from the tensor
void structure_factors::calc_F_calc(const dMatrix2& D) {
	if (XG.mode == 1) { grid_F(D); return; }
	calc_F_calc(I_tens, D);
	if (XG.mode != 2) return;
	const cvec tensor = scatter_data.F_calc;
	grid_F(D);
	double d = 0, m = 0;
	for (size_t r = 0; r < tensor.size(); r++) { d = std::max(d, std::abs(tensor[r] - scatter_data.F_calc[r])); m = std::max(m, std::abs(tensor[r])); }
	std::cout << "XCW grid check: max |F_grid - F_tensor| " << std::scientific << d << " of max |F| " << m << std::defaultfloat << std::endl;
	scatter_data.F_calc = tensor;
}

//F_calc from the scattering factors of a tsc table: rows matched by hkl, scatterers by label, by
//position where the labels differ
void structure_factors::calc_F_calc(const itsc_block& tsc) {
	std::map<i3, size_t> row_of;
	for (size_t r = 0; r < tsc.reflection_size(); r++)
		row_of[{ (int)tsc.get_index(0, r), (int)tsc.get_index(1, r), (int)tsc.get_index(2, r) }] = r;
	const svec names = tsc.get_scatterers_string();
	cvec2 sfs(model_data.ncen);
	for (int at = 0; at < model_data.ncen; at++) {
		size_t s = std::find(names.begin(), names.end(), asym_atoms[at].label) - names.begin();
		if (s == names.size()) {
			if (at == 0) std::cout << "F_calc from a tsc table: scatterer " << names[0] << " for " << asym_atoms[0].label << ", matched by position" << std::endl;
			s = at;
		}
		for (const i3& h : scatter_data.hkl_enlarged)
			sfs[at].push_back(tsc.get_sf_for_scatterer(s)[row_of.at(h)]);
	}
	calc_F_calc(sfs);
}

//The tsc table of the wavefunctions at every entry of hkl_enlarged, by the usual tsc route
itsc_block structure_factors::tsc_table(std::vector<WFN>& wavefunctions, GridManager* grids, std::ostream& log) {
	svec known_atoms;
	opt->m_hkl_list = scatter_data.hkl_enlarged;
	opt->grid_cache = grids;
	return calculate_scattering_factors<itsc_block, std::vector<WFN>&>(*opt, wavefunctions, log, known_atoms, 0, &k_pt);
}

void structure_factors::calc_F_calc(std::vector<WFN>& wavefunctions, GridManager* grids, std::ostream& log) {
	calc_F_calc(tsc_table(wavefunctions, grids, log));
}

void structure_factors::perturbation(occ::Mat& perturb, const bool unrestricted, std::ostream& log) {
	ensure_inv_H2_weights();

	//The four (XWR_type, refine_against) combinations differ only in the per-reflection
	//scalar and the prefactor, so one walk over I serves all of them
	const int xwr = static_cast<int>(opt->xcw_settings.XWR_type), ref = static_cast<int>(opt->xcw_settings.refine_against);
	const bool against_F2 = (ref == 2);
	const bool weighted = (xwr == 2);
	const bool valid = (xwr == 1 || xwr == 2) && (ref == 1 || ref == 2);
	if (!valid) log << "Invalid refinement option" << std::endl;
	const double scale_sq = scatter_data.scale * scatter_data.scale;
	const double prefactor = against_F2
		? 4.0 * scale_sq / (model_data.nr_fit - n_params())
		: 2.0 * scatter_data.scale / (model_data.nr_fit - n_params());

	cvec pre(model_data.nr);
#pragma omp parallel for
	for (int r = 0; r < model_data.nr; r++) {
		if (!valid || !scatter_data.hkl_mask[r]) continue;
		cdouble precompute;
		//with extinction the model is I = y |Fc|^2, so the residual carries sqrt(y)|Fc| and the
		//carrier d/dD picks up dI/d|Fc|^2 (F^2) or d(sqrt(y)|Fc|)/d|Fc| (F); both are 1 without
		//a model, and the expressions below are then exactly the ones this always used
		if (against_F2) {
			const double I = ext_y(r) * std::pow(std::abs(scatter_data.F_calc[r]), 2);
			precompute = ext_g(r) * std::conj(scatter_data.F_calc[r]) * (scale_sq * I - scatter_data.F_obs2[r]) / (scatter_data.sigma_obs2[r] * scatter_data.sigma_obs2[r]);
		}
		else {
			const double F_calc_abs = std::abs(scatter_data.F_calc[r]);
			precompute = ext_m(r) * std::conj(scatter_data.F_calc[r]) * (scatter_data.scale * ext_sqrt_y(r) * F_calc_abs - scatter_data.abs_F_obs[r]) / (scatter_data.sigma_obs[r] * scatter_data.sigma_obs[r] * F_calc_abs);
		}
		if (weighted) precompute *= inv_H2_[r];
		pre[r] = precompute;
	}
	contract_I(perturb, pre);
	perturb *= prefactor;
	if (unrestricted) {
		perturb.conservativeResize(2 * model_data.nmo, Eigen::NoChange);
		perturb.bottomRows(model_data.nmo) = perturb.topRows(model_data.nmo);
	}
	//closing function
}

occ::Mat structure_factors::perturbation_response(const dMatrix2& dD_eff, const double lambda) {
	//dF_r = Sum c I_r dD_eff: calc_F_calc's walk with the fixed part zeroed
	const cvec F0 = scatter_data.F_calc, fixed = scatter_data.anom_correction;
	std::fill(scatter_data.anom_correction.begin(), scatter_data.anom_correction.end(), cdouble(0.0, 0.0));
	calc_F_calc(dD_eff);
	const cvec dFr = scatter_data.F_calc;
	scatter_data.F_calc = F0;
	scatter_data.anom_correction = fixed;
	//the scale's response from the weighted least-squares scale of eval_scale
	const bool against_F2 = opt->xcw_settings.refine_against == 2, weighted = opt->xcw_settings.XWR_type == 2;
	const double k = scatter_data.scale, s = k * k;
	const int chunk = 128, nchunk = (model_data.nr + chunk - 1) / chunk;
	vec dnum(nchunk), dden(nchunk), den(nchunk);
#pragma omp parallel for schedule(static)
	for (int ch = 0; ch < nchunk; ch++) {
		for (int r = ch * chunk; r < std::min((ch + 1) * chunk, model_data.nr); r++) {
			if (!scatter_data.hkl_mask[r]) continue;
			//ponytail: the extinction shape is frozen over the step (dy/d|Fc| dropped), so
			//sqrt(y)|Fc| and m d|Fc| reproduce I and dI exactly but their own curvature is
			//neglected. Only the Hessian's step proposal degrades - TRAH accepts on the
			//exact energy rebuild_at returns, and perturbation's gradient stays exact.
			const double w = weighted ? inv_H2_[r] : 1.0, Fa = ext_sqrt_y(r) * std::abs(F0[r]);
			const double dFa = std::abs(F0[r]) > 0 ? ext_m(r) * std::real(std::conj(F0[r]) * dFr[r]) / std::abs(F0[r]) : 0.0;
			if (against_F2) {
				const double wi = w / (scatter_data.sigma_obs2[r] * scatter_data.sigma_obs2[r]);
				dnum[ch] += wi * 2.0 * Fa * dFa * scatter_data.F_obs2[r];
				dden[ch] += wi * 4.0 * Fa * Fa * Fa * dFa;
				den[ch] += wi * Fa * Fa * Fa * Fa;
			}
			else {
				const double wi = w / (scatter_data.sigma_obs[r] * scatter_data.sigma_obs[r]);
				dnum[ch] += wi * dFa * scatter_data.F_obs[r];
				dden[ch] += wi * 2.0 * Fa * dFa;
				den[ch] += wi * Fa * Fa;
			}
		}
	}
	double dnum_sum = 0, dden_sum = 0, den_sum = 0;
	for (int ch = 0; ch < nchunk; ch++) { dnum_sum += dnum[ch]; dden_sum += dden[ch]; den_sum += den[ch]; }
	//F: dk = (dN - k dD) / D; F^2: ds = (dN - s dD) / D for s = k^2
	const double dscale = den_sum != 0 ? (dnum_sum - (against_F2 ? s : k) * dden_sum) / den_sum : 0.0;
	//d(q_r) for q_r = prefactor's scale power times perturbation's scalar, both moving:
	//F: q = w conj(F) (k^2 - k Fo / |F|) / sigma^2; F^2: q = w conj(F) (s^2 |F|^2 - s Fo^2) / sigma_I^2
	cvec dq(model_data.nr);
#pragma omp parallel for
	for (int r = 0; r < model_data.nr; r++) {
		if (!scatter_data.hkl_mask[r]) continue;
		const double w = weighted ? inv_H2_[r] : 1.0, Fm = std::abs(F0[r]);
		if (Fm == 0) continue;
		//as above: the model amplitude is sqrt(y)|Fc| and its response m d|Fc|, with the
		//shape frozen. The carrier conj(F) keeps perturbation's chain factor.
		const double Fa = ext_sqrt_y(r) * Fm, dFa = ext_m(r) * std::real(std::conj(F0[r]) * dFr[r]) / Fm;
		const cdouble carrier = (against_F2 ? ext_g(r) : ext_m(r)) * std::conj(F0[r]);
		const cdouble dcarrier = (against_F2 ? ext_g(r) : ext_m(r)) * std::conj(dFr[r]);
		if (against_F2) {
			const double wi = w / (scatter_data.sigma_obs2[r] * scatter_data.sigma_obs2[r]), Fo2 = scatter_data.F_obs2[r];
			dq[r] = wi * (dcarrier * (s * s * Fa * Fa - s * Fo2)
				+ carrier * (2.0 * s * dscale * Fa * Fa + 2.0 * s * s * Fa * dFa - dscale * Fo2));
		}
		else {
			const double wi = w / (scatter_data.sigma_obs[r] * scatter_data.sigma_obs[r]), Fo = scatter_data.abs_F_obs[r];
			dq[r] = wi * (dcarrier * (k * k * ext_sqrt_y(r) - k * Fo / Fm)
				+ carrier * (dscale * (2.0 * k * ext_sqrt_y(r) - Fo / Fm) + k * Fo * std::real(std::conj(F0[r]) * dFr[r]) / (Fm * Fm * Fm)));
		}
	}
	occ::Mat dP;
	contract_I(dP, dq);
	dP *= (against_F2 ? 4.0 : 2.0) / (model_data.nr_fit - n_params()) * lambda;
	return dP;
}

void structure_factors::contract_I(occ::Mat& out, const cvec& pre) {
	if (XG.mode == 1) { grid_contract(out, pre); return; }
	if (XG.mode == 2) {
		occ::Mat grid;
		grid_contract(grid, pre);
		contract_tensor(out, pre);
		std::cout << "XCW grid check: max |P_grid - P_tensor| " << std::scientific << (grid - out).cwiseAbs().maxCoeff()
			<< " of max |P| " << out.cwiseAbs().maxCoeff() << std::defaultfloat << std::endl;
		return;
	}
	contract_tensor(out, pre);
}

void structure_factors::contract_tensor(occ::Mat& out, const cvec& pre) {
	out.setZero(model_data.nmo, model_data.nmo);
	const int step = std::max(1, I_tens.i_streamed_ ? I_tens.i_window_ : model_data.nr);
	bool on_device = false;
#if defined(NOSPHERA2_USE_GPU) || defined(NOSPHERA2_USE_METAL)
	//GPU: the same contraction over the device copy of the tensor; the CPU walk below is skipped
	if (I_tens.i_on_device_) {
		vec dev(I_tens.i_compact_);
		err_checkf(itensor_gpu_cols(pre.data(), dev.data()), "I tensor walk on the device failed", std::cout);
		for (size_t k = 0; k < I_tens.i_compact_; k++) out(I_tens.i_pair_mu_[k], I_tens.i_pair_nu_[k]) = dev[k];
		on_device = true;
	}
#endif
	// See calc_F_calc: an exception must not leave an OpenMP structured block.
	std::string io_error;
	//One region, one accumulator per thread, one reduction, however many windows the
	//budget implies: an nmo x nmo matrix per window is what made narrow windows dear
	std::vector<occ::Mat> parts(omp_get_max_threads());
#pragma omp parallel if (!on_device)
	{
		occ::Mat local = occ::Mat::Zero(model_data.nmo, model_data.nmo);
		double* local_ptr = local.data();
		for (int r0 = 0; !on_device && r0 < model_data.nr; r0 += step) {
			const int r1 = std::min(r0 + step, model_data.nr);
			if (I_tens.i_streamed_) {
#pragma omp single
				{
					try { I_tens.i_file_.load(r0, r1); }
					catch (const std::exception& e) { io_error = e.what(); }
				}
			}
			//No nowait: the next window's read must not start until every thread has read this one
#pragma omp for
			for (int r = r0; r < r1; r++) {
				if (!io_error.empty()) continue;
				const cdouble precompute = pre[r];
				if (precompute == cdouble(0.0, 0.0)) continue;
				//As in calc_F_calc: one walk over whichever element type is resident, with
				//the accumulation in double either way.
				auto accumulate = [&](const auto* I_r) {
					const int* pmu = I_tens.i_pair_mu_.data();
					const int* pnu = I_tens.i_pair_nu_.data();
					for (size_t k = 0; k < I_tens.i_compact_; k++) {
						const double vr = static_cast<double>(I_r[k].real());
						const double vi = static_cast<double>(I_r[k].imag());
						local_ptr[pnu[k] * model_data.nmo + pmu[k]] += precompute.real() * vr - precompute.imag() * vi;
					}
				};
				if (I_tens.i_float_) accumulate(I_tens.i_block32(r)); else accumulate(I_tens.i_block(r));
			}
		}
		parts[omp_get_thread_num()].swap(local);
	}
	if (!io_error.empty()) throw std::runtime_error(io_error);
	//The partials outweigh the matrix many times over, so their sum is parallel too
	if (!on_device) {
#pragma omp parallel for schedule(static)
		for (int mu = 0; mu < model_data.nmo; mu++) {
			for (int nu = mu; nu < model_data.nmo; nu++) {
				double sum = 0.0;
				for (int t = 0; t < static_cast<int>(parts.size()); t++)
					if (parts[t].size() != 0) sum += parts[t](mu, nu);
				out(mu, nu) = sum;
			}
		}
	}
	for (int mu = 0; mu < model_data.nmo; mu++) {
		for (int nu = mu + 1; nu < model_data.nmo; nu++) {
			out(nu, mu) = out(mu, nu);
		}
	}
}

//The tsc table of the converged density D_wfn for the tscb output. The grid route already holds rho's
//partition: its per-atom transform at D is the table, so after the first step, whose tsc-route table
//gives the layout, no second grid, no WFN pass and no tsc machinery run per lambda.
//NOSPHERA2_XCW_GRID_TSCB_OLD keeps the tsc route on every step
itsc_block structure_factors::result_table(const occ::Mat& D_wfn, std::vector<WFN>& wavefunctions, std::ostream& log) {
	itsc_block result;
	const bool grid_tscb = XG.mode == 1 && !opt->electron_diffraction && !tuning("NOSPHERA2_XCW_GRID_TSCB_OLD");
	cvec sfc;
	const int nc = static_cast<int>(XG.inversion_link.size());
	auto grid_value = [&](const int at, const int e) {
		const cdouble v = sfc[static_cast<size_t>(at) * nc + XG.clean_of[e]];
		return XG.is_inverse[e] ? std::conj(v) : v;
	};
	if (grid_tscb) {
		//the spin blocks summed, as the SCF evaluates the structure factors
		const int n = static_cast<int>(D_wfn.cols());
		const bool unrestricted = D_wfn.rows() == 2 * D_wfn.cols();
		dMatrix2 D(n, n);
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++) D(i, j) = D_wfn(i, j) + (unrestricted ? D_wfn(i + n, j) : 0.0);
		sfc = grid_sfc(D);
	}
	if (grid_tscb && !tscb_scat_.empty()) {
		cvec2 table(tscb_scat_.size(), cvec(tscb_entry_.size()));
#pragma omp parallel for
		for (int scatterer = 0; scatterer < static_cast<int>(table.size()); scatterer++)
			for (size_t r = 0; r < tscb_entry_.size(); r++) table[scatterer][r] = grid_value(tscb_atom_[scatterer], tscb_entry_[r]);
		result = itsc_block(table, tscb_scat_, tscb_hkl_);
	}
	else {
		result.append(tsc_table(wavefunctions, &tsc_grids, log), log);
		//first step on the grid route: learn the table's layout, each scatterer's atom by the closest grid column
		if (grid_tscb) {
			std::map<i3, int> e_of;
			for (const i3& h : scatter_data.hkl_enlarged) e_of.emplace(h, static_cast<int>(e_of.size()));
			const int nrow = static_cast<int>(result.reflection_size()), ns = static_cast<int>(result.get_scatterers().size());
			ivec entry(nrow), atom(ns);
			bool ok = true;
			for (int r = 0; r < nrow && ok; r++) {
				const auto it = e_of.find({ (int)result.get_index(0, r), (int)result.get_index(1, r), (int)result.get_index(2, r) });
				ok = it != e_of.end();
				if (ok) entry[r] = it->second;
			}
			double worst = 0.0;
			for (int scatterer = 0; scatterer < ns && ok; scatterer++) {
				const cvec& f = result.get_sf_for_scatterer(scatterer);
				double best = 1e300;
				for (int at = 0; at < model_data.ncen; at++) {
					double d = 0;
					for (int r = 0; r < std::min(nrow, 256); r++) d += std::norm(f[r] - grid_value(at, entry[r]));
					if (d < best) { best = d; atom[scatterer] = at; }
				}
				double dmax = 0, fmax = 0;
				for (int r = 0; r < nrow; r++) {
					dmax = std::max(dmax, std::abs(f[r] - grid_value(atom[scatterer], entry[r])));
					fmax = std::max(fmax, std::abs(f[r]));
				}
				worst = std::max(worst, dmax / fmax);
			}
			//1e-2 only catches a wrong layout; how well grid and tsc agree is what NOSPHERA2_XCW_GRID_FCHECK measures
			ok = ok && worst < 1e-2;
			std::cout << "XCW grid tscb: max |f_grid - f_tsc| / max |f| " << std::scientific << std::setprecision(2) << worst << std::defaultfloat
				<< (ok ? ", later steps write the tscb from the grid" : ", layout not matched, the tsc route stays") << std::endl;
			if (ok) {
				tscb_scat_ = result.get_scatterers();
				tscb_hkl_ = result.get_index_vector();
				tscb_entry_ = entry;
				tscb_atom_ = atom;
			}
		}
	}
	return result;
}

void structure_factors::save_F_calc() {
	saved_F_calc_ = scatter_data.F_calc;
	saved_scale_ = scatter_data.scale;
}

void structure_factors::restore_F_calc() {
	scatter_data.F_calc = saved_F_calc_;
	scatter_data.scale = saved_scale_;
}
