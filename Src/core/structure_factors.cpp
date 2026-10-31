#include "pch.h"
#include "tuning.h"
#include "structure_factors.h"
#include "wfn_class.h"
#include "scattering_factors.h"
#include "GridManager.h"
#include "nos_math.h"
#include "basis_set.h"
#include <limits>
#include <numeric>


structure_factors::structure_factors(options& opt_in) {
	opt = &opt_in;

	// Load the settings from the options object
	model_data.n_params = opt_in.xcw_settings.n_params;
	extinction_settings.extinction_model = opt_in.xcw_settings.extinction_model;
	orbital_basis = opt_in.xcw_settings.basis_overrides.empty() ? BasisSetLibrary::get_basis_set(opt_in.xcw_settings.basis_set_name) :
		BasisSetLibrary::get_basis_set_with_overrides(opt_in.xcw_settings.basis_set_name, opt_in.xcw_settings.basis_overrides);
	extinction_settings.aniso = opt_in.xcw_settings.extinction_aniso;
	extinction_settings.refine = opt_in.xcw_settings.extinction_refine;
	extinction_settings.start = opt_in.xcw_settings.extinction_start;
	quality_criteria.refine_against = opt_in.xcw_settings.refine_against;
	quality_criteria.goof_type = opt_in.xcw_settings.XWR_type;


	// Setup file paths
	std::filesystem::path hkl_filename = opt_in.hkl;
	std::filesystem::path cif = opt_in.cif;
	std::ifstream cif_input(cif.c_str(), std::ios::in);

	// Read the hkl and cif file and generate the unit cell
	unit_cell = cell(cif, std::cout, opt_in.debug, opt_in.do_XCW);
	scatter_data.hkl_enlarged = read_hkl_full(hkl_filename, scatter_data.hkl, opt_in.twin_law, unit_cell, std::cout, scatter_data, opt_in.debug);
	double wavelength_ = read_CIF(cif_input, unit_cell, model_data.ncen, asym_atoms, ADPs, opt_in.debug);
	// The ADPs as the CIF gives them, for the model summary; U_cif2U_star below rescales them in place
	std::ostringstream adp_summary;
	adp_summary << std::fixed << std::setprecision(6)
		<< "ADPs used (Angstrom^2, U11 U22 U33 U12 U13 U23 or Uiso); grown atoms take their parent's, rotated by its operation\n";
	for (int a = 0; a < asym_atoms.size(); a++) {
		adp_summary << "  " << std::left << std::setw(6) << asym_atoms[a].label << std::right;
		if (ADPs[a].empty() || ADPs[a][0].empty()) adp_summary << "  Uiso " << asym_atoms[a].U_iso;
		else for (double u : ADPs[a][0]) adp_summary << std::setw(11) << u;
		adp_summary << "\n";
	}
	for (int order = 1; order < 3; order++) {
		bool header = false;
		for (int a = 0; a < asym_atoms.size(); a++) {
			if (ADPs[a].size() <= order || ADPs[a][order].empty()) continue;
			if (!header) adp_summary << "Gram-Charlier " << (order == 1 ? "C (111 112 113 122 123 133 222 223 233 333)" : "D (1111 ... 3333, CIF order)") << " as read\n";
			header = true;
			adp_summary << "  " << std::left << std::setw(6) << asym_atoms[a].label << std::right << std::scientific << std::setprecision(4);
			for (double c : ADPs[a][order]) adp_summary << std::setw(12) << c;
			adp_summary << std::fixed << std::setprecision(6) << "\n";
		}
	}
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
		// Here, before apply_grown deletes operations, grown_by still indexes the CIF's operations
		std::ostringstream summary;
		summary << std::fixed << std::setprecision(5);
		int n_asym = 0;
		for (const asym_atom& at : asym_atoms) n_asym += !at.grown;
		summary << "XCW model: " << n_asym << " asymmetric unit atoms, " << asym_atoms.size() - n_asym << " grown, "
			<< asym_atoms.size() << " in the calculation; " << unit_cell.get_num_syms() << " symmetry operations\n"
			<< "Asymmetric unit:\n";
		for (int a = 0; a < asym_atoms.size(); a++) {
			const asym_atom& at = asym_atoms[a];
			if (at.grown) continue;
			summary << std::setw(5) << a << " " << std::left << std::setw(6) << at.label << std::setw(3) << constants::atnr2letter(at.type) << std::right
				<< std::setw(10) << at.frac_pos[0] << std::setw(10) << at.frac_pos[1] << std::setw(10) << at.frac_pos[2]
				<< "  site multiplicity 1/" << at.cryst_mult << "\n";
		}
		summary << adp_summary.str();
		if (n_asym < asym_atoms.size()) {
			summary << "Grown atoms: parent, operation, lattice translation, fractional coordinates\n";
			const vec2 t = unit_cell.get_trans();
			for (int a = 0; a < asym_atoms.size(); a++) {
				const asym_atom& at = asym_atoms[a];
				if (!at.grown) continue;
				summary << std::setw(5) << a << " " << std::setw(3) << constants::atnr2letter(at.type);
				if (at.grown_from < 0) {
					summary << "  NO PARENT: no operation maps an asymmetric unit atom onto it";
				}
				else {
					const int op = at.grown_by[0];
					const asym_atom& p = asym_atoms[at.grown_from];
					std::ostringstream shift;
					for (int i = 0; i < 3; i++) {
						double image = t[i][op];
						for (int j = 0; j < 3; j++) image += unit_cell.get_sym(j, i, op) * p.frac_pos[j];
						shift << (i ? "," : "") << std::lround(at.frac_pos[i] - image);
					}
					summary << "  from " << std::left << std::setw(6) << p.label << std::setw(22) << unit_cell.symop_string(op) << std::right
						<< std::setw(10) << "(" + shift.str() + ")";
					if (at.grown_by.size() > 1) summary << " +" << at.grown_by.size() - 1 << " more op(s)";
				}
				summary << std::setw(10) << at.frac_pos[0] << std::setw(10) << at.frac_pos[1] << std::setw(10) << at.frac_pos[2] << "\n";
			}
		}
		if (opt_in.xcw_settings.grown) {
			// Fourth step of the symmetry resolution: Copy U_iso, the dispersion and the ADPs of each grown atom's parent, the ADPs rotated onto the symmetry generated atoms
			unit_cell.grow_ADPs(asym_atoms, ADPs);
		}
		if (opt_in.xcw_settings.grown && tuning("NOSPHERA2_XCW_GROWN_COPYSUM")) {
			// Old route: delete the fully applied operations and sum every grown copy's own f over the kept ones. A site's copies
			// are then not averaged over the operations, which is not the F Olex2 refines from the tscb (R1 0.777 vs 0.748 % on PPh4UF6)
			applied_symmetry = unit_cell.apply_grown(scatter_data.hkl_table, scatter_data.hkl_enlarged, asym_atoms);
		}
		else {
			// Symmetry factors of an ungrown model, also for an atom on a special position. A grown cluster only shapes the
			// density and its partition: F is the asymmetric unit's f over all operations, as Olex2 computes it from the tscb
			unit_cell.set_sym_fact_from_mult(asym_atoms);
			for (asym_atom& at : asym_atoms) if (at.grown) at.sym_fact = 0.0;
		}
		summary << "Operations F_calc sums over (" << unit_cell.get_num_syms() << "):";
		for (int op = 0; op < unit_cell.get_num_syms(); op++) summary << "  " << unit_cell.symop_string(op);
		summary << "\nsym_fact:";
		for (int a = 0; a < asym_atoms.size(); a++)
			summary << (a % 8 ? "  " : "\n  ") << std::left << std::setw(6) << (asym_atoms[a].grown ? "#" + std::to_string(a) : asym_atoms[a].label) << std::right
			<< std::setprecision(4) << asym_atoms[a].sym_fact;
		summary << "\n";
		model_summary = summary.str();
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
	const ivec2 lookup = unit_cell.generate_asym_lookup(scatter_data.hkl, scatter_data.hkl_enlarged);
	for (int r = 0; r < model_data.nr; r++) {
		const ivec& asym_list = lookup[r];
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
	const ivec2 asym_lookup = unit_cell.generate_asym_lookup(scatter_data.hkl, scatter_data.hkl_enlarged);
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
				? (s * ext_y_[r] * std::norm(scatter_data.F_calc[r]) - scatter_data.F_obs2[r]) / sig_I_[r]
				: (k * ext_sqrt_y_[r] * std::abs(scatter_data.F_calc[r]) - scatter_data.abs_F_obs[r]) / sig_F_[r];
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
			resid = (s * ext_y_[r] * Fm * Fm - scatter_data.F_obs2[r]) / sig_I_[r];
			dmodel = s * Fm * Fm * ext_dyc_[r] / sig_I_[r];
		}
		else {
			if (ext_sqrt_y_[r] <= 0.0) continue;
			resid = (k * ext_sqrt_y_[r] * Fm - scatter_data.abs_F_obs[r]) / sig_F_[r];
			dmodel = k * Fm * ext_dyc_[r] / (2.0 * ext_sqrt_y_[r] * sig_F_[r]);
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
	ensure_cryst_weights();
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
				const double calc2 = calc * calc, wi = w / (sig_I_[i] * sig_I_[i]);
				numerators[c] += wi * calc2 * scatter_data.F_obs2[i];
				denominators[c] += wi * calc2 * calc2;
			}
			else {
				const double wi = w / (sig_F_[i] * sig_F_[i]);
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
	ensure_cryst_weights();
	//index 0: the fit set, 1: all reflections
	const double prefactor[2] = { 1.0 / static_cast<double>(model_data.nr_fit - n_params()), 1.0 / static_cast<double>(model_data.nr - n_params()) };
	const int chunk = 128, nchunk = (model_data.nr + chunk - 1) / chunk;
	vec2 goof1(2, vec(nchunk)), goof2(2, vec(nchunk)), wgoof1(2, vec(nchunk)), wgoof2(2, vec(nchunk)), r1_num(2, vec(nchunk)), r1_den(2, vec(nchunk)), wr2_num(2, vec(nchunk)), wr2_den(2, vec(nchunk));
	const double scale = scatter_data.scale;
	const cdouble* F_calc_0 = scatter_data.F_calc.data();
#pragma omp parallel for schedule(static)
	for (int c = 0; c < nchunk; c++) {
		const int first = c * chunk, last = std::min(first + chunk, model_data.nr);
		for (int i = first; i < last; i++) {
			const double scaled_F_calc = scale * ext_sqrt_y(i) * std::abs(F_calc_0[i]);
			const double scaled_difference = scaled_F_calc - scatter_data.F_obs[i];
			const double diff2 = (scaled_F_calc * scaled_F_calc) - scatter_data.F_obs2[i];
			const double weighted_diff1 = scaled_difference / sig_F_[i];
			const double weighted_diff2 = diff2 / sig_I_[i];
			const double weighted_diff1_sq = weighted_diff1 * weighted_diff1;
			const double weighted_diff2_sq = weighted_diff2 * weighted_diff2;
			const double w = quality_criteria.goof_type == 2 ? inv_H2_[i] : 0.0;
			//wR2 takes Fo^2 with its sign, as SHELXL does; F_obs2 holds |Fo^2|
			const double fo2 = scatter_data.F_obs[i] < 0 ? -scatter_data.F_obs2[i] : scatter_data.F_obs2[i];
			const double wr_d = (scaled_F_calc * scaled_F_calc - fo2) / sig_I_[i], wr_o = fo2 / sig_I_[i];
			for (int set = scatter_data.hkl_mask[i] ? 0 : 1; set < 2; set++) {
				wr2_num[set][c] += wr_d * wr_d;
				wr2_den[set][c] += wr_o * wr_o;
				goof1[set][c] += weighted_diff1_sq;
				goof2[set][c] += weighted_diff2_sq;
				wgoof1[set][c] += weighted_diff1_sq * w;
				wgoof2[set][c] += weighted_diff2_sq * w;
				r1_num[set][c] += std::abs(scaled_F_calc - scatter_data.abs_F_obs[i]);
				r1_den[set][c] += scatter_data.abs_F_obs[i];
			}
		}
	}
	double sum[8][2] = {};
	for (int set = 0; set < 2; set++) {
		for (int c = 0; c < nchunk; c++) {
			sum[6][set] += wr2_num[set][c];
			sum[7][set] += wr2_den[set][c];
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
	quality_criteria.wR2 = sum[7][0] > 0.0 ? std::sqrt(sum[6][0] / sum[7][0]) : 0.0;
	quality_criteria.wR2_all = sum[7][1] > 0.0 ? std::sqrt(sum[6][1] / sum[7][1]) : 0.0;
}

//cctbx smtbx/refinement/weighting_schemes.h, mainstream_shelx_weighting::compute
double shelx_weight(const double fo_sq, const double sigma, const double fc_sq, const double K, const double a, const double b) {
	const double p = (std::max(fo_sq, 0.) + 2 * K * fc_sq) / 3.;
	return 1. / (sigma * sigma + std::pow(a * p, 2) + b * K * p);
}

void structure_factors::ensure_cryst_weights() {
	if (!sig_I_.empty()) return;
	sig_F_ = scatter_data.sigma_obs;
	sig_I_ = scatter_data.sigma_obs2;
	using ws = options::SF_settings::weight_scheme;
	const ws scheme = opt->xcw_settings.weight;
	if (scheme == ws::sigma && opt->xcw_settings.robust_c <= 0) return;
	rescale_pending_ = opt->xcw_settings.weight_rescale && (scheme == ws::unit || scheme == ws::stl || scheme == ws::chebyshev);
	solve_scale();
	update_cryst_weights();
}

double chebyshev_weight(const double t, const vec& A) {
	double sum = 0.0, T0 = 1.0, T1 = t;
	for (size_t r = 0; r < A.size(); r++) {
		sum += A[r] * (r == 0 ? T0 : T1);
		if (r > 0) { const double T2 = 2.0 * t * T1 - T0; T0 = T1; T1 = T2; }
	}
	return sum > 0.0 ? 1.0 / sum : 0.0;
}

vec fit_chebyshev(const vec& t, const vec& d2, const int n) {
	Eigen::MatrixXd M(static_cast<Eigen::Index>(t.size()), n);
	Eigen::VectorXd b(static_cast<Eigen::Index>(t.size()));
	for (Eigen::Index i = 0; i < M.rows(); i++) {
		double T0 = 1.0, T1 = t[i];
		for (int r = 0; r < n; r++) {
			M(i, r) = r == 0 ? T0 : T1;
			if (r > 0) { const double T2 = 2.0 * t[i] * T1 - T0; T0 = T1; T1 = T2; }
		}
		b(i) = d2[i];
	}
	const Eigen::VectorXd A = M.colPivHouseholderQr().solve(b);
	return vec(A.data(), A.data() + A.size());
}

vec tukey_biweight(const vec& z, const ivec& in_scale, const double c) {
	vec a, f(z.size(), 1.0);
	for (size_t i = 0; i < z.size(); i++)
		if (in_scale[i] && std::isfinite(z[i])) a.push_back(std::abs(z[i]));
	if (a.empty()) return f;
	const size_t h = a.size() / 2;
	std::nth_element(a.begin(), a.begin() + h, a.end());
	const double med = a.size() % 2 ? a[h] : 0.5 * (a[h] + *std::max_element(a.begin(), a.begin() + h));
	const double cs = c * med / 0.6745;
	if (!(cs > 0.0)) return f;
	for (size_t i = 0; i < z.size(); i++) {
		if (!std::isfinite(z[i])) continue;
		const double u = z[i] / cs, v = 1.0 - u * u;
		f[i] = std::abs(u) < 1.0 ? v * v : 0.0;
	}
	return f;
}

//Delta of the refined metric, as the criterion has it: k sqrt(y)|Fc| - Fo, or k^2 y|Fc|^2 - |Fo^2|
double structure_factors::fit_residual(const int r) const {
	const double c = scatter_data.scale * ext_sqrt_y(r) * std::abs(scatter_data.F_calc[r]);
	return quality_criteria.refine_against == 2 ? c * c - scatter_data.F_obs2[r] : c - scatter_data.F_obs[r];
}

//Chebyshev's variable from Fc whatever is refined: t = 2 Fc/max(Fc) - 1, Fc = sqrt(y)|F_calc| over all reflections
vec structure_factors::chebyshev_t() const {
	vec t(model_data.nr);
	double mx = 0.0;
	for (int r = 0; r < model_data.nr; r++) mx = std::max(mx, t[r] = ext_sqrt_y(r) * std::abs(scatter_data.F_calc[r]));
	for (double& x : t) x = mx > 0.0 ? 2.0 * x / mx - 1.0 : -1.0;
	return t;
}

std::string structure_factors::weight_report() const {
	using ws = options::SF_settings::weight_scheme;
	static const char* names[] = { "sigma", "shelx", "shelx_opt", "unit", "stl", "chebyshev", "chebyshev_opt" };
	const auto& xs = opt->xcw_settings;
	std::ostringstream s;
	s << std::setprecision(6);
	if (xs.robust_c > 0) s << "robust(c " << xs.robust_c << ") ";
	s << names[static_cast<int>(xs.weight)];
	const bool cheb = xs.weight == ws::chebyshev || xs.weight == ws::chebyshev_opt;
	const vec& p = weight_ab_.empty() ? xs.weight_params : weight_ab_;
	if (!p.empty())
		s << (xs.weight == ws::shelx || xs.weight == ws::shelx_opt ? " WGHT" : xs.weight == ws::stl ? " a" : (cheb && weight_ab_.empty() && xs.weight == ws::chebyshev_opt) ? " n" : " A");
	for (const double v : p) s << " " << v;
	if (weight_k_ != 1.0) s << " k " << weight_k_;
	return s.str();
}

//The default keeps the read sigmas. SHELX: w on Fo^2 from cctbx's formula with Fc^2 = y|Fc|^2 and K = scale^2,
//both scales put Fo^2 ~ scale^2 y |Fc|^2. Against F the same w is propagated as the reader propagates sigma,
//sigma(F) = sigma(Fo^2)/(2|Fo|), i.e. w_F = 4 Fo^2 w: first-order error propagation, SHELXL's for L.S. on F
//and cctbx's f_sq_as_f away from Fo ~ 0 (there sigma(F) is infinite and the reflection weighs nothing, as now)
//unit, stl and Chebyshev: w in the refined metric (on F against F, on Fo^2 against F^2), the other sigma propagated
//the same way, sigma(Fo^2) = 2|Fo| sigma(F); w = 0 is sigma = inf, so the reflection drops out of every sum.
//robust: Tukey's biweight on z = Delta sqrt(w) of the refined metric, its scale over the fit set, w' = f w
void structure_factors::update_cryst_weights() {
	using ws = options::SF_settings::weight_scheme;
	const ws scheme = opt->xcw_settings.weight;
	const int nr = model_data.nr;
	const bool F2 = quality_criteria.refine_against == 2;
	const double inf = std::numeric_limits<double>::infinity();
	if (scheme == ws::chebyshev_opt && weight_ab_.empty()) {
		optimise_cryst_weights();
		return;
	}
	if (weight_ab_.empty()) weight_ab_ = opt->xcw_settings.weight_params;
	if (scheme == ws::sigma) {
		sig_F_ = scatter_data.sigma_obs;
		sig_I_ = scatter_data.sigma_obs2;
	}
	else if (scheme == ws::shelx || scheme == ws::shelx_opt) {
		const double K = scatter_data.scale * scatter_data.scale, a = weight_ab_[0], b = weight_ab_[1];
		sig_F_.resize(nr);
		sig_I_.resize(nr);
		for (int r = 0; r < nr; r++) {
			const double fo2 = scatter_data.F_obs[r] < 0 ? -scatter_data.F_obs2[r] : scatter_data.F_obs2[r];
			sig_I_[r] = 1.0 / std::sqrt(shelx_weight(fo2, scatter_data.sigma_obs2[r], ext_y(r) * std::norm(scatter_data.F_calc[r]), K, a, b));
			sig_F_[r] = sig_I_[r] * 0.5 / scatter_data.abs_F_obs[r];
		}
	}
	else {
		const vec t = scheme == ws::chebyshev || scheme == ws::chebyshev_opt ? chebyshev_t() : vec();
		if (scheme == ws::stl) ensure_hkl_ordered();
		sig_F_.resize(nr);
		sig_I_.resize(nr);
		for (int r = 0; r < nr; r++) {
			const double w = scheme == ws::unit ? 1.0 : scheme == ws::stl ? std::pow(unit_cell.get_stl_of_hkl(hkl_ordered_[r]), weight_ab_[0]) : chebyshev_weight(t[r], weight_ab_);
			const double s = w > 0.0 ? 1.0 / std::sqrt(w) : inf;
			if (F2) {
				sig_I_[r] = s;
				sig_F_[r] = s * 0.5 / scatter_data.abs_F_obs[r];
			}
			else {
				const double si = 2.0 * scatter_data.abs_F_obs[r] * s;
				sig_F_[r] = s;
				sig_I_[r] = si > 0.0 ? si : inf;
			}
		}
		//REWEIGHT: k = 1/S^2 of the base scheme at its own scale (which k leaves alone), S^2 as criterion() has it
		if (rescale_pending_) {
			rescale_pending_ = false;
			solve_scale();
			const bool weighted = quality_criteria.goof_type == 2;
			double chi2 = 0.0;
			for (int r = 0; r < nr; r++) {
				if (!scatter_data.hkl_mask[r]) continue;
				const double z = fit_residual(r) / (F2 ? sig_I_[r] : sig_F_[r]);
				chi2 += z * z * (weighted ? inv_H2_[r] : 1.0);
			}
			const double S2 = chi2 / static_cast<double>(model_data.nr_fit - n_params());
			if (S2 > 0.0 && std::isfinite(S2)) weight_k_ = 1.0 / S2;
			std::ostringstream s;
			s << "XCW weight rescale k=" << std::setprecision(6) << weight_k_ << " (S=1 at scan start)";
			std::cout << s.str() << std::endl;
			if (weight_log_) *weight_log_ << s.str() << std::endl;
		}
		if (weight_k_ != 1.0) {
			const double f = 1.0 / std::sqrt(weight_k_);
			for (int r = 0; r < nr; r++) {
				sig_F_[r] *= f;
				sig_I_[r] *= f;
			}
		}
	}
	if (opt->xcw_settings.robust_c <= 0) return;
	vec z(nr);
	for (int r = 0; r < nr; r++) z[r] = fit_residual(r) / (F2 ? sig_I_[r] : sig_F_[r]);
	const vec f = tukey_biweight(z, scatter_data.hkl_mask, opt->xcw_settings.robust_c);
	for (int r = 0; r < nr; r++) {
		if (f[r] == 1.0) continue;
		const double k = f[r] > 0.0 ? 1.0 / std::sqrt(f[r]) : inf;
		sig_F_[r] *= k;
		sig_I_[r] *= k;
	}
}

namespace {
	//The data of new_shelx_weighting.optimise_parameters after its sort: Fo^2 clipped at 0, sigma(Fo^2) and
	//Fc^2 in ascending Fc^2, the 10 bin limits, K and its n = N // (N - n_params)
	struct shelx_bins {
		vec fo, sig, fc;
		size_t lim[11];
		double K;
		double n;
		//calcres / the grid's variance: Sum over the bins of (GooF_bin - 1)^2, GooF_bin = sqrt(chi^2 n / count)
		double variance(const double a, const double b) const {
			double v = 0.0;
			for (int i = 0; i < 10; i++) {
				const size_t count = lim[i + 1] - lim[i];
				if (count == 0) continue;
				double chi = 0.0;
				for (size_t j = lim[i]; j < lim[i + 1]; j++) {
					const double d = fo[j] - K * fc[j];
					chi += shelx_weight(fo[j], sig[j], fc[j], K, a, b) * d * d;
				}
				const double g = std::sqrt(chi * n / static_cast<double>(count)) - 1.0;
				v += g * g;
			}
			return v;
		}
		void residuals(const double a, const double b, double e[10]) const {
			for (int i = 0; i < 10; i++) {
				const size_t count = lim[i + 1] - lim[i];
				double chi = 0.0;
				for (size_t j = lim[i]; j < lim[i + 1]; j++) {
					const double d = fo[j] - K * fc[j];
					chi += shelx_weight(fo[j], sig[j], fc[j], K, a, b) * d * d;
				}
				e[i] = count ? std::sqrt(chi * n / static_cast<double>(count)) - 1.0 : 0.0;
			}
		}
	};

	//optimise_parameters_old: SHELXL's 9 x 9 grid search, step for step (the last point of a tie wins)
	vec shelx_grid_search(const shelx_bins& B, double start_a, double start_b) {
		double a_step = 0.2 * start_a, b_step = 0.4 * start_b;
		while (a_step > 1e-4 && b_step > 5e-3) {
			start_a = std::max(start_a, 4 * a_step) - 4 * a_step;
			start_b = std::max(start_b, 4 * b_step) - 4 * b_step;
			double min_variance = 9e9;
			int j_min = 0, k_min = 0;
			double b = start_b;
			for (int j = 0; j < 9; j++) {
				double a = start_a;
				b += b_step;
				for (int k = 0; k < 9; k++) {
					a += a_step;
					const double variance = B.variance(a, b);
					min_variance = std::min(variance, min_variance);
					if (variance == min_variance) { j_min = j; k_min = k; }
				}
			}
			start_a += k_min * a_step;
			start_b += j_min * b_step;
			if (k_min == 8) { a_step *= 2; continue; }
			else if (k_min != 0) a_step /= 4;
			if (j_min == 8) { b_step *= 2; continue; }
			else if (j_min != 0) b_step /= 4;
			if (start_a <= 1e-4) a_step /= 4;
			if (start_b <= 1e-3) b_step /= 4;
		}
		if (start_a > 0.2) { start_a = 0.2; start_b = 0; }
		return { start_a, start_b };
	}
}

//cctbx smtbx/refinement/weighting_schemes.py, new_shelx_weighting.optimise_parameters: the start values, the
//sort, the 10 bins, n and calcres exactly. Where it calls scipy's SLSQP with a, b >= 0, this runs a bounded
//Levenberg-Marquardt on the 10 bin residuals GooF_bin - 1 (whose square sum is calcres): the same minimum of
//the same function from the same start; where SLSQP fails it falls back to optimise_parameters_old, so does this
vec optimise_shelx_weights(const vec& fo_sq, const vec& sigma, const vec& fc_sq, const double K, const int n_params) {
	const size_t N = fo_sq.size();
	if (N == 0 || K <= 0.0) return { 0.1, 0.0 };
	shelx_bins B;
	B.K = K;
	//start values, SHELXL's formulae on Fo^2/K, sigma/K and Fc^2
	double x = 0, y = 0, z = 0;
	for (size_t i = 0; i < N; i++) {
		const double fo2 = std::max(fo_sq[i], 0.0) / K, s = sigma[i] / K, s2 = s * s, fc2 = fc_sq[i];
		const double p = (fo2 + 2. * fc2) / 3., p2 = p * p;
		x += (std::pow(fo2 - fc2, 2) - s) * (p2 / s2);
		y += p2 * p2 / s2;
		z += p;
	}
	const double start_a = std::sqrt(std::max(0.0001, 0.64 * x / std::max(1e-8, y)));
	const double start_b = 0.5 * z * start_a * start_a / static_cast<double>(N);
	//sorted by Fc^2 / max(Fc^2)
	std::vector<size_t> perm(N);
	std::iota(perm.begin(), perm.end(), size_t(0));
	std::stable_sort(perm.begin(), perm.end(), [&](const size_t i, const size_t j) { return fc_sq[i] < fc_sq[j]; });
	B.fo.resize(N); B.sig.resize(N); B.fc.resize(N);
	for (size_t i = 0; i < N; i++) {
		B.fo[i] = std::max(fo_sq[perm[i]], 0.0);
		B.sig[i] = sigma[perm[i]];
		B.fc[i] = fc_sq[perm[i]];
	}
	B.lim[0] = 0;
	for (int i = 0; i < 10; i++)
		B.lim[i + 1] = static_cast<size_t>(std::ceil(static_cast<double>((i + 1) * N) / 10.0));
	const long long Nl = static_cast<long long>(N), denom = Nl - n_params;
	B.n = static_cast<double>(denom > 0 ? Nl / denom : 1);

	//bounded Levenberg-Marquardt, the Jacobian by differences, steps clipped onto a, b >= 0
	double p[2] = { start_a, start_b }, e[10], f = B.variance(p[0], p[1]);
	double mu = 1e-3;
	bool ok = std::isfinite(f);
	for (int it = 0; ok && it < 500; it++) {
		double J[10][2], ep[10], em[10];
		B.residuals(p[0], p[1], e);
		for (int c = 0; c < 2; c++) {
			const double h = 1e-7 * std::max(std::abs(p[c]), 1e-3);
			double q[2] = { p[0], p[1] }, lo = std::max(p[c] - h, 0.0);
			q[c] = p[c] + h; B.residuals(q[0], q[1], ep);
			q[c] = lo; B.residuals(q[0], q[1], em);
			for (int i = 0; i < 10; i++) J[i][c] = (ep[i] - em[i]) / (p[c] + h - lo);
		}
		double A[2][2] = {}, g[2] = {};
		for (int i = 0; i < 10; i++)
			for (int c = 0; c < 2; c++) {
				g[c] += J[i][c] * e[i];
				for (int d = 0; d < 2; d++) A[c][d] += J[i][c] * J[i][d];
			}
		const double floor = 1e-12 * (A[0][0] + A[1][1]) + 1e-300;
		bool moved = false;
		double np[2] = { p[0], p[1] }, nf = f;
		for (; mu < 1e12; mu *= 10) {
			const double d00 = A[0][0] + mu * std::max(A[0][0], floor), d11 = A[1][1] + mu * std::max(A[1][1], floor);
			const double det = d00 * d11 - A[0][1] * A[1][0];
			if (!(det > 0)) continue;
			np[0] = std::max(0.0, p[0] + (-g[0] * d11 + g[1] * A[0][1]) / det);
			np[1] = std::max(0.0, p[1] + (-g[1] * d00 + g[0] * A[1][0]) / det);
			nf = B.variance(np[0], np[1]);
			if (nf < f) { moved = true; break; }
		}
		if (!moved) break;
		const double dp = std::abs(np[0] - p[0]) + std::abs(np[1] - p[1]);
		p[0] = np[0]; p[1] = np[1];
		const double df = f - nf;
		f = nf;
		mu = std::max(mu / 10, 1e-12);
		if (dp <= 1e-12 * (1.0 + std::abs(p[0]) + std::abs(p[1])) || df <= 1e-16 * std::max(f, 1e-300)) break;
	}
	ok = ok && std::isfinite(f) && std::isfinite(p[0]) && std::isfinite(p[1]);
	if (ok) return { p[0], p[1] };
	std::cout << "Warning: optimisation of weighting parameters failed, falling back to old method." << std::endl;
	return shelx_grid_search(B, start_a, start_b);
}

void structure_factors::optimise_cryst_weights() {
	if (opt->xcw_settings.weight == options::SF_settings::weight_scheme::chebyshev_opt) {
		const vec t = chebyshev_t();
		vec tf, d2;
		for (int r = 0; r < model_data.nr; r++) {
			if (!scatter_data.hkl_mask[r]) continue;
			const double d = fit_residual(r);
			tf.push_back(t[r]);
			d2.push_back(d * d);
		}
		weight_ab_ = fit_chebyshev(tf, d2, static_cast<int>(opt->xcw_settings.weight_params[0]));
		update_cryst_weights();
		return;
	}
	std::vector<double> fo, sig, fc;
	for (int r = 0; r < model_data.nr; r++) {
		if (!scatter_data.hkl_mask[r]) continue;
		fo.push_back(scatter_data.F_obs[r] < 0 ? -scatter_data.F_obs2[r] : scatter_data.F_obs2[r]);
		sig.push_back(scatter_data.sigma_obs2[r]);
		fc.push_back(ext_y(r) * std::norm(scatter_data.F_calc[r]));
	}
	weight_ab_ = optimise_shelx_weights(fo, sig, fc, scatter_data.scale * scatter_data.scale, n_params());
	update_cryst_weights();
}

double structure_factors::criterion(const bool all) const {
	const bool weighted = quality_criteria.goof_type == 2, against_F2 = quality_criteria.refine_against == 2;
	const double goof = all ? (weighted ? (against_F2 ? quality_criteria.weighted_GooF2_all : quality_criteria.weighted_GooF1_all) : (against_F2 ? quality_criteria.GooF2_all : quality_criteria.GooF1_all))
		: (weighted ? (against_F2 ? quality_criteria.weighted_GooF2 : quality_criteria.weighted_GooF1) : (against_F2 ? quality_criteria.GooF2 : quality_criteria.GooF1));
	return goof * goof;
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

void structure_factors::setup_grid(const occ::gto::AOBasis& aobasis) {
	if (!DW_set_) std::cout << "Debye-Waller factors are not computed and set to 1. Continuing." << std::endl;
	if (!phases_set_) std::cout << "Phase factors are not computed and set to 1. Continuing." << std::endl;
	if (!anom_set_) std::cout << "Anomalous dispersion corrections are not computed and set to 0. Continuing." << std::endl;
	std::vector<ao_data> ao_data_shells;
	create_prims(ao_data_shells, aobasis);
	XG.fp32 = opt->cpu_itensor_fp32 && !opt->xcw_settings.grid_double;
	link_symmetry_equivalents();
	screen_ao_pairs(ao_data_shells);
	//TFVC, MBIS and EMBIS shares need rho: their grid is Hirshfeld's, with plain atom weights and pruned on the
	//promolecule, and setup_grid_weights takes it from there
	const bool share = opt->partition_type == PartitionType::TFVC || opt->partition_type == PartitionType::MBIS || opt->partition_type == PartitionType::EMBIS;
	GridConfiguration config;
	config.accuracy = opt->accuracy;
	config.partition_type = share ? PartitionType::Hirshfeld : opt->partition_type;
	config.no_density_eval = true;
	config.pbc = opt->pbc;
	config.debug = opt->debug;
	config.all_charges = opt->all_charges;
	GridManager grid_manager(config);
	WFN dummy_wave;
	sf_shell_layout layout;
	setup_atom_grids(grid_manager, dummy_wave, layout);
	build_grid_blocks(ao_data_shells, grid_manager, dummy_wave, layout);
	//the Fourier-Bessel transform of every atom grid, planned once on its prototype
	plan_grid_fb(layout);
	if (share) setup_grid_weights(aobasis);
}

void structure_factors::link_symmetry_equivalents() {
	XG.lookup = unit_cell.generate_asym_lookup(scatter_data.hkl, scatter_data.hkl_enlarged);
	//The AO values and weights are real, so the grid integral of -q is the conjugate of that of q: only the entries of
	//inversion_cleaned_hkl are integrated
	if (tuning("NOSPHERA2_NO_INVERSION_SCREENING")) {
		XG.inversion_link.clear();
		for (int e = 0; e < model_data.nr_enlarged; e++) XG.inversion_link.push_back({ e });
	}
	else
		XG.inversion_link = unit_cell.remove_inv_vectors(scatter_data.hkl_enlarged, scatter_data.inversion_cleaned_hkl);
	XG.clean_of.assign(model_data.nr_enlarged, 0);
	XG.is_inverse.assign(model_data.nr_enlarged, false);
	for (int c = 0; c < static_cast<int>(XG.inversion_link.size()); c++)
		for (int i = 0; i < static_cast<int>(XG.inversion_link[c].size()); i++) {
			XG.clean_of[XG.inversion_link[c][i]] = c;
			XG.is_inverse[XG.inversion_link[c][i]] = i == 1;
		}
	if (static_cast<int>(XG.inversion_link.size()) != model_data.nr_enlarged && !(opt->no_date)) {
		const size_t n_inv = scatter_data.hkl_enlarged.size() - scatter_data.inversion_cleaned_hkl.size();
		std::cout << std::fixed << std::setprecision(2) << "Inversion screening: Screened out " << n_inv << " of " << scatter_data.hkl_enlarged.size()
			<< " symmetry generated reflections (" << 100.0 * n_inv / scatter_data.hkl_enlarged.size() << "%)" << std::endl;
	}
}

//An AO pair mu, nu on different centres is skipped when the bound on its overlap, summed over the primitives k of mu and l
//of nu, falls below e_tol
void structure_factors::screen_ao_pairs(const std::vector<ao_data>& ao_data_shells) {
	const _time_point screening_start = get_time();
	const int nmo = model_data.nmo;
	XG.skip.assign(nmo, ivec(nmo, 0));
	const double e_tol = 0.0005;
	const double sqrt_inv_four_pi = std::sqrt(constants::INV_FOUR_PI);
	//the largest value of each primitive's radial part times the largest of its spherical harmonic
	auto peak_factors = [&](const ao_data& ao) {
		const int type = ao.prims[0].get_type();
		const double type_half = type * 0.5, sph_max = constants::spherical_harmonic_max(type, ao.m);
		vec N(ao.prims.size());
		for (int k = 0; k < static_cast<int>(ao.prims.size()); k++)
			N[k] = type == 0 ? sqrt_inv_four_pi : sph_max * std::exp(-type_half) * std::pow(type / ao.prims[k].get_exp(), type_half);
		return N;
	};
	long long screened = 0;
	for (int mu = 0; mu < nmo; mu++) {
		const std::vector<primitive>& mu_primitives = ao_data_shells[mu].prims;
		const vec N_k = peak_factors(ao_data_shells[mu]);
		for (int nu = mu + 1; nu < nmo; nu++) {
			const d3& a = ao_data_shells[mu].pos, & b = ao_data_shells[nu].pos;
			const double dist2 = (a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]);
			if (dist2 < 1e-5) continue;
			const std::vector<primitive>& nu_primitives = ao_data_shells[nu].prims;
			const vec N_l = peak_factors(ao_data_shells[nu]);
			double estimate = 0.0;
			for (int k = 0; k < static_cast<int>(mu_primitives.size()); k++)
				for (int l = 0; l < static_cast<int>(nu_primitives.size()); l++) {
					const double k_exp = mu_primitives[k].get_exp(), l_exp = nu_primitives[l].get_exp();
					const double inv_exp_sum = 1.0 / (k_exp + l_exp);
					const double N_kl = N_k[k] * N_l[l] * std::pow(constants::TWO_PI * inv_exp_sum, 1.5);
					estimate += std::abs(mu_primitives[k].get_coef() * nu_primitives[l].get_coef()) * N_kl * std::exp(-0.5 * k_exp * l_exp * inv_exp_sum * dist2);
				}
			if (estimate < e_tol) {
				XG.skip[mu][nu] = 1;
				screened++;
			}
		}
	}
	if (!opt->no_date)
		std::cout << "Time taken for AO screening: " << get_msec(screening_start, get_time()) << " ms" << std::endl;
	std::cout << std::fixed << std::setprecision(2) << "Screened out " << screened << " unique pairs of mu, nu ("
		<< screened / (static_cast<double>(nmo) * (nmo + 1) / 2) * 100.0 << "%) " << std::endl;
}

void structure_factors::setup_atom_grids(GridManager& grid_manager, WFN& dummy_wave, sf_shell_layout& layout) {
	//the atoms of the asymmetric unit with the orbital basis, which the grids are built for
	for (int at = 0; at < model_data.ncen; at++) {
		atom temp_atom;
		for (int x = 0; x < 3; x++) temp_atom.set_coordinate(x, asym_atoms[at].pos[x]);
		temp_atom.set_charge(asym_atoms[at].type);
		dummy_wave.push_back_atom(temp_atom);
	}
	load_basis_into_WFN(dummy_wave, orbital_basis, false, true);
	dummy_wave.delete_unoccupied_MOs();
	ivec asym_atom_list(model_data.ncen);
	for (int at = 0; at < model_data.ncen; at++) asym_atom_list[at] = at;
	grid_manager.setup3DGridsForMolecule(dummy_wave, asym_atom_list, bvec(model_data.ncen, true), unit_cell);
	//a flat density for the pruning, which is all the partition weights need here
	GridData& GD = grid_manager.getGridData();
	std::vector<vec2>& grids = grid_manager.getNeedsHelper() ? GD.helper_grids : GD.atomic_grids;
	for (vec2& grid : grids)
		std::fill(grid[GridData::GridIndex::WFN_DENSITY].begin(), grid[GridData::GridIndex::WFN_DENSITY].end(), 1.0);
	grid_manager.getDensityVectors(dummy_wave, asym_atom_list, XG.d1, XG.d2, XG.d3, XG.w, false, &layout);
	std::cout << "Total number of grid points after pruning: " << grid_manager.getTotalGridPoints() << std::endl;
}

void structure_factors::build_grid_blocks(const std::vector<ao_data>& ao_data_shells, GridManager& grid_manager, const WFN& dummy_wave, sf_shell_layout& layout) {
	const int nmo = model_data.nmo;
	GridData& GD = grid_manager.getGridData();
	vec2* grids = grid_manager.getNeedsHelper() ? GD.helper_grids.data() : GD.atomic_grids.data();
	const int* points = grid_manager.getNeedsHelper() ? GD.helper_num_points_per_atom.data() : GD.num_points_per_atom.data();
	const int n_grids = static_cast<int>(grid_manager.getNeedsHelper() ? GD.helper_grids.size() : GD.atomic_grids.size());
	const int n_atom_grids = std::min(n_grids, model_data.ncen);
	//Each AO counts out to a cutoff of 11 to 12 bohr, the wider the more diffuse its widest primitive
	constexpr double maximum_ao_grid_cutoff = 12, minimum_ao_grid_cutoff = 11;
	vec ao_min_exp(nmo, std::numeric_limits<double>::max()), ao_grid_cutoff_squared(nmo);
	for (int ao = 0; ao < nmo; ao++)
		for (const primitive& prim : ao_data_shells[ao].prims) ao_min_exp[ao] = std::min(ao_min_exp[ao], prim.get_exp());
	const double min_exp = *std::min_element(ao_min_exp.begin(), ao_min_exp.end());
	for (int ao = 0; ao < nmo; ao++) {
		const double cutoff = std::clamp(maximum_ao_grid_cutoff * std::sqrt(min_exp / ao_min_exp[ao]), minimum_ao_grid_cutoff, maximum_ao_grid_cutoff);
		ao_grid_cutoff_squared[ao] = cutoff * cutoff;
	}
	//the distance of every point from its atom, and per AO and atom grid whether the AO reaches the atom and up to which
	//point (points are radially ordered) an AO on that atom is not cut off
	vec2 grid_radial_distances(n_atom_grids);
	for (int g = 0; g < n_atom_grids; g++) {
		grid_radial_distances[g].resize(points[g]);
		for (int p = 0; p < points[g]; p++) {
			const double dx = grids[g][GridData::GridIndex::X][p] - asym_atoms[g].pos[0];
			const double dy = grids[g][GridData::GridIndex::Y][p] - asym_atoms[g].pos[1];
			const double dz = grids[g][GridData::GridIndex::Z][p] - asym_atoms[g].pos[2];
			grid_radial_distances[g][p] = std::sqrt(dx * dx + dy * dy + dz * dz);
		}
	}
	ivec2 ao_prefix_end(nmo, ivec(n_atom_grids));
	bvec2 ao_within_cutoff(nmo, bvec(n_atom_grids));
#pragma omp parallel for schedule(static)
	for (int ao = 0; ao < nmo; ao++) {
		const double cutoff2 = ao_grid_cutoff_squared[ao];
		for (int g = 0; g < n_atom_grids; g++) {
			const double dx = asym_atoms[g].pos[0] - ao_data_shells[ao].pos[0];
			const double dy = asym_atoms[g].pos[1] - ao_data_shells[ao].pos[1];
			const double dz = asym_atoms[g].pos[2] - ao_data_shells[ao].pos[2];
			const double d2 = dx * dx + dy * dy + dz * dz;
			ao_within_cutoff[ao][g] = d2 <= cutoff2;
			ao_prefix_end[ao][g] = d2 < 1e-12
				? static_cast<int>(std::upper_bound(grid_radial_distances[g].begin(), grid_radial_distances[g].end(), std::sqrt(cutoff2)) - grid_radial_distances[g].begin())
				: points[g];
		}
	}
	//AO values on every point within their cutoff
	const _time_point start_AOs = get_time();
	vec3 mu_vals(nmo, vec2(n_grids));
#pragma omp parallel for schedule(dynamic)
	for (int mu = 0; mu < nmo; mu++) {
		const ao_data& ao = ao_data_shells[mu];
		for (int g = 0; g < n_grids; g++) {
			const double* x_ptr = grids[g][GridData::GridIndex::X].data();
			const double* y_ptr = grids[g][GridData::GridIndex::Y].data();
			const double* z_ptr = grids[g][GridData::GridIndex::Z].data();
			mu_vals[mu][g].resize(points[g]);
			double* values = mu_vals[mu][g].data();
			const int prefix_end = g < n_atom_grids ? ao_prefix_end[mu][g] : points[g];
			for (int p = 0; p < prefix_end; p++) {
				d4 d_mu{ x_ptr[p] - ao.pos[0], y_ptr[p] - ao.pos[1], z_ptr[p] - ao.pos[2], 0 };
				d_mu[3] = d_mu[0] * d_mu[0] + d_mu[1] * d_mu[1] + d_mu[2] * d_mu[2];
				values[p] = d_mu[3] > ao_grid_cutoff_squared[mu] ? 0.0 : dummy_wave.eval_ao(d_mu, ao.prims, ao.m, std::sqrt(d_mu[3]));
			}
		}
	}
	std::cout << "AO values calculated for all grids." << std::endl;
	if (!(opt->no_date))
		std::cout << "Time taken for AO values computation: " << get_msec(start_AOs, get_time()) << " ms" << std::endl;
	//Morton order makes a run of consecutive points a compact ball rather than a spherical shell, so that a block keeps
	//few AOs; the work is sum over blocks of n_active^2 points. Sums over points do not depend on the order, and only
	//an AO threshold that can prune makes the reordering worth it
	const double ao_block_threshold = [&] {
		const char* e = tuning("NOSPHERA2_XCW_GRID_AO_TOL");
		return e ? std::max(0.0, std::atof(e)) : cutoff(opt->accuracy);
	}();
	const bool morton_applied = !tuning("NOSPHERA2_XCW_GRID_NO_MORTON") && ao_block_threshold >= 1e-20;
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
			//21 bits per axis, interleaved to a 63-bit key
			auto spread = [](const unsigned int v) {
				unsigned long long x = v & 0x1fffffu;
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
			//every array that is paired with a point moves with it
			if (static_cast<int>(XG.d1[g].size()) >= npts) apply(XG.d1[g].data());
			if (static_cast<int>(XG.d2[g].size()) >= npts) apply(XG.d2[g].data());
			if (static_cast<int>(XG.d3[g].size()) >= npts) apply(XG.d3[g].data());
			if (static_cast<int>(XG.w[g].size()) >= npts) apply(XG.w[g].data());
			if (g < static_cast<int>(layout.idx.size()) && static_cast<int>(layout.idx[g].size()) >= npts) {
				ivec t(npts);
				for (int p = 0; p < npts; p++) t[p] = layout.idx[g][perm[p]];
				layout.idx[g].swap(t);
			}
			for (int mu = 0; mu < nmo; mu++) {
				vec& v = mu_vals[mu][g];
				if (v.empty()) continue;
				//values filled only to the radial prefix: the zero tail is mixed in by the permutation
				v.resize(npts, 0.0);
				apply(v.data());
			}
			vec& rd = grid_radial_distances[g];
			if (static_cast<int>(rd.size()) == npts) apply(rd.data());
		}
	}
	//the values of the AOs that reach each atom grid, one row of points per AO
	ivec2 grid_active_aos(n_atom_grids);
	vec2 grid_ao_values(n_atom_grids);
#pragma omp parallel for schedule(dynamic)
	for (int g = 0; g < n_atom_grids; g++)
		for (int mu = 0; mu < nmo; mu++) {
			if (!ao_within_cutoff[mu][g]) continue;
			grid_active_aos[g].push_back(mu);
			grid_ao_values[g].insert(grid_ao_values[g].end(), mu_vals[mu][g].begin(), mu_vals[mu][g].end());
		}
	vec3().swap(mu_vals);
	//Blocks of about NOSPHERA2_XCW_GRID_CHUNK points (1024, 0: the three radial bands), each keeping the AOs above
	//ao_block_threshold somewhere on it; the threshold is the -acc ladder the scattering factors screen on
	const int chunk = [] {
		const char* e = tuning("NOSPHERA2_XCW_GRID_CHUNK");
		return e ? std::atoi(e) : 1024;
	}();
	XG.blocks.assign(n_atom_grids, {});
	long long ao_slots_carrying = 0, ao_slots_kept = 0;
#pragma omp parallel for schedule(dynamic) reduction(+:ao_slots_carrying, ao_slots_kept)
	for (int g = 0; g < n_atom_grids; g++) {
		const ivec& active_aos = grid_active_aos[g];
		const vec& radial_distances = grid_radial_distances[g];
		const int inner_end = static_cast<int>(std::upper_bound(radial_distances.begin(), radial_distances.end(), minimum_ao_grid_cutoff) - radial_distances.begin());
		const int middle_end = static_cast<int>(std::upper_bound(radial_distances.begin(), radial_distances.end(), maximum_ao_grid_cutoff) - radial_distances.begin());
		const std::array<int, 4> band_bounds{ 0, inner_end, middle_end, points[g] };
		//even chunks rather than a short tail, which would cost a GEMM call for almost nothing
		auto cut = [&](const int from, const int to, std::vector<std::pair<int, int>>& out) {
			const int n = to - from;
			if (n <= 0) return;
			if (chunk <= 0) { out.emplace_back(from, to); return; }
			const int pieces = std::max(1, (n + chunk - 1) / chunk);
			const int per = (n + pieces - 1) / pieces;
			for (int p0 = from; p0 < to; p0 += per) out.emplace_back(p0, std::min(p0 + per, to));
		};
		std::vector<std::pair<int, int>> spans;
		if (morton_applied)
			cut(0, points[g], spans);
		else
			for (int band = 0; band < 3; band++) cut(band_bounds[band], band_bounds[band + 1], spans);
		for (const auto& [point_start, point_end] : spans) {
			if (point_end == point_start) continue;
			GridBlock block{ point_start, point_end - point_start };
			for (int local_ao = 0; local_ao < static_cast<int>(active_aos.size()); local_ao++) {
				const double* full_row = grid_ao_values[g].data() + static_cast<size_t>(local_ao) * points[g];
				double largest = 0.0;
				for (int p = point_start; p < point_end; p++) largest = std::max(largest, std::abs(full_row[p]));
				if (largest > 0.0) ao_slots_carrying++;
				if (largest > ao_block_threshold) {
					block.active_aos.push_back(active_aos[local_ao]);
					block.ao_values.insert(block.ao_values.end(), full_row + point_start, full_row + point_end);
				}
			}
			if (block.active_aos.empty()) continue;
			ao_slots_kept += static_cast<long long>(block.active_aos.size());
			if (XG.fp32) block.ao_values_f.assign(block.ao_values.begin(), block.ao_values.end());
			XG.blocks[g].push_back(std::move(block));
		}
	}
	if (!(opt->no_date) && ao_slots_carrying > 0) {
		const long long dropped = ao_slots_carrying - ao_slots_kept;
		std::cout << std::fixed << std::setprecision(2) << "Screened out " << dropped << " of " << ao_slots_carrying << " AO-block entries ("
			<< 100.0 * static_cast<double>(dropped) / static_cast<double>(ao_slots_carrying) << "%) below "
			<< std::scientific << std::setprecision(0) << ao_block_threshold << std::fixed << std::setprecision(2) << " on their block\n";
	}
	//-gflops: the work of one evaluation, sum over blocks of n_active^2 points, so a chunk size can be judged
	if (throughput::enabled()) {
		double work = 0.0;
		long long nblocks = 0;
		for (int g = 0; g < n_atom_grids; g++)
			for (const GridBlock& b : XG.blocks[g]) {
				work += static_cast<double>(b.active_aos.size()) * b.active_aos.size() * b.point_count;
				nblocks++;
			}
		std::fprintf(stderr, "XCW grid blocks: %lld, sum n_active^2 * points = %.3e\n", nblocks, work);
	}
	//the integration weights alone, before any density-following share multiplies them in
	XG.w0.clear();
	const bool share = opt->partition_type == PartitionType::TFVC || opt->partition_type == PartitionType::MBIS || opt->partition_type == PartitionType::EMBIS;
	for (int g = 0; share && g < n_atom_grids; g++)
		XG.w0.emplace_back(grids[g][GridData::GridIndex::WEIGHT].begin(), grids[g][GridData::GridIndex::WEIGHT].begin() + points[g]);
}

//What the next large allocation may take: the platform's and scheduler's headroom and, with -mem,
//-mem minus what this process already holds. The density-fitting integrals and the stored integrals
//are sized in that order from it, so together they stay under -mem. OCC fills the density-fitting
//store only at its first Fock build, after both decisions, so its size is held back here by hand.
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

//F_r = Sum_at sym_fact Sum_s f_at(h R_s) exp(i h R_s x_at) T_at(h R_s) exp(2 pi i h t_s) plus the
//anomalous dispersion correction, with the atomic scattering factors laid out like DW_facts
void structure_factors::calc_F_calc(const cvec2& scattering_factors) {
	const ivec2 lookups = unit_cell.generate_asym_lookup(scatter_data.hkl, scatter_data.hkl_enlarged);
#pragma omp parallel for
	for (int r = 0; r < model_data.nr; r++) {
		const ivec& lookup = lookups[r];
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

//F_calc of D on XCW's grid, see xcw_grid
void structure_factors::calc_F_calc(const dMatrix2& D) {
	grid_F(D);
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
	ensure_cryst_weights();

	//The four (XWR_type, refine_against) combinations differ only in the per-reflection
	//scalar and the prefactor, so one contraction serves all of them
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
			precompute = ext_g(r) * std::conj(scatter_data.F_calc[r]) * (scale_sq * I - scatter_data.F_obs2[r]) / (sig_I_[r] * sig_I_[r]);
		}
		else {
			const double F_calc_abs = std::abs(scatter_data.F_calc[r]);
			precompute = ext_m(r) * std::conj(scatter_data.F_calc[r]) * (scatter_data.scale * ext_sqrt_y(r) * F_calc_abs - scatter_data.abs_F_obs[r]) / (sig_F_[r] * sig_F_[r] * F_calc_abs);
		}
		if (weighted) precompute *= inv_H2_[r];
		pre[r] = precompute;
	}
	grid_contract(perturb, pre);
	perturb *= prefactor;
	if (unrestricted) {
		perturb.conservativeResize(2 * model_data.nmo, Eigen::NoChange);
		perturb.bottomRows(model_data.nmo) = perturb.topRows(model_data.nmo);
	}
	//closing function
}

occ::Mat structure_factors::perturbation_response(const dMatrix2& dD_eff, const double lambda) {
	//dF_r = dF_r/dD dD_eff: calc_F_calc with the fixed part zeroed
	const cvec F0 = scatter_data.F_calc, fixed = scatter_data.anom_correction;
	std::fill(scatter_data.anom_correction.begin(), scatter_data.anom_correction.end(), cdouble(0.0, 0.0));
	XR.response = true;
	calc_F_calc(dD_eff);
	XR.response = false;
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
				const double wi = w / (sig_I_[r] * sig_I_[r]);
				dnum[ch] += wi * 2.0 * Fa * dFa * scatter_data.F_obs2[r];
				dden[ch] += wi * 4.0 * Fa * Fa * Fa * dFa;
				den[ch] += wi * Fa * Fa * Fa * Fa;
			}
			else {
				const double wi = w / (sig_F_[r] * sig_F_[r]);
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
			const double wi = w / (sig_I_[r] * sig_I_[r]), Fo2 = scatter_data.F_obs2[r];
			dq[r] = wi * (dcarrier * (s * s * Fa * Fa - s * Fo2)
				+ carrier * (2.0 * s * dscale * Fa * Fa + 2.0 * s * s * Fa * dFa - dscale * Fo2));
		}
		else {
			const double wi = w / (sig_F_[r] * sig_F_[r]), Fo = scatter_data.abs_F_obs[r];
			dq[r] = wi * (dcarrier * (k * k * ext_sqrt_y(r) - k * Fo / Fm)
				+ carrier * (dscale * (2.0 * k * ext_sqrt_y(r) - Fo / Fm) + k * Fo * std::real(std::conj(F0[r]) * dFr[r]) / (Fm * Fm * Fm)));
		}
	}
	occ::Mat dP;
	grid_contract(dP, dq);
	dP *= (against_F2 ? 4.0 : 2.0) / (model_data.nr_fit - n_params()) * lambda;
	return dP;
}

//The tsc table of the converged density D_wfn for the tscb output. The grid route already holds rho's
//partition: its per-atom transform at D is the table, so after the first step, whose tsc-route table
//gives the layout, no second grid, no WFN pass and no tsc machinery run per lambda.
//NOSPHERA2_XCW_GRID_TSCB_OLD keeps the tsc route on every step
itsc_block structure_factors::result_table(const occ::Mat& D_wfn, std::vector<WFN>& wavefunctions, std::ostream& log, const bool refit_shares) {
	itsc_block result;
	const bool grid_tscb = !opt->electron_diffraction && !tuning("NOSPHERA2_XCW_GRID_TSCB_OLD");
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
		XG.refresh = XG.ws && refit_shares;
		sfc = grid_sfc(D);
		XG.refresh = false;
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
			const svec names = result.get_scatterers_string();
			for (int scatterer = 0; scatterer < ns && ok; scatterer++) {
				const cvec& f = result.get_sf_for_scatterer(scatterer);
				double best = 1e300;
				//by the atom's id or label, as write_tshb keys them: two atoms of one element can be closer in f across a partition
				//change (RI against the tsc route's Hirshfeld) than to themselves; else the closest grid column
				atom[scatterer] = -1;
				for (int at = 0; at < model_data.ncen && at < wavefunctions[0].get_ncen() && !tuning("NOSPHERA2_XCW_TSCB_BY_F"); at++)
					if (wavefunctions[0].get_id_for_atom(at).to_hex_string() == names[scatterer] || asym_atoms[at].label == names[scatterer]) atom[scatterer] = at;
				const bool named = atom[scatterer] >= 0;
				for (int at = 0; at < model_data.ncen && !named; at++) {
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
			//the RI partition is not the tsc route's, nor are shares fit tighter than the stand-alone MBIS and EMBIS stop, so there
			//the closest atom is all the layout can be matched by
			ok = ok && (worst < 1e-2 || XR.on || XG.ws);
			std::cout << "XCW grid tscb: max |f_grid - f_tsc| / max |f| " << std::scientific << std::setprecision(2) << worst << std::defaultfloat
				<< (ok ? ", later steps write the tscb from the grid" : ", layout not matched, the tsc route stays") << std::endl;
			if (ok) {
				tscb_scat_ = result.get_scatterers();
				tscb_hkl_ = result.get_index_vector();
				tscb_entry_ = entry;
				tscb_atom_ = atom;
			}
			if (ok && (XR.on || XG.ws)) {
				cvec2 table(tscb_scat_.size(), cvec(tscb_entry_.size()));
				for (int scatterer = 0; scatterer < static_cast<int>(table.size()); scatterer++)
					for (size_t r = 0; r < tscb_entry_.size(); r++) table[scatterer][r] = grid_value(tscb_atom_[scatterer], tscb_entry_[r]);
				result = itsc_block(table, tscb_scat_, tscb_hkl_);
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

//F_obs, sigma, scale * |F_calc| and the phase of every reflection, for an R-factor against another route
//(an ORCA tsc through Olex2, or a Laue check between symmetry-equivalent rows)
void structure_factors::write_F_calc(const std::string& path) {
	ensure_hkl_ordered();
	std::ofstream fc(path);
	fc << "#    h    k    l          F_obs        sig(F)   scale*|F_calc|     phase(deg)   R1(gt) = " << std::setprecision(5) << quality_criteria.R1 << " R1(all) = " << quality_criteria.R1_all << " scale = " << std::setprecision(10) << scatter_data.scale << "\n";
	for (int r = 0; r < model_data.nr; r++) {
		const cdouble& f = scatter_data.F_calc[r];
		fc << std::setw(5) << hkl_ordered_[r][0] << std::setw(5) << hkl_ordered_[r][1] << std::setw(5) << hkl_ordered_[r][2]
			<< std::fixed << std::setprecision(4) << std::setw(15) << scatter_data.F_obs[r] << std::setw(14) << scatter_data.sigma_obs[r]
			<< std::setw(17) << scatter_data.scale * ext_sqrt_y(r) * std::abs(f) << std::setw(15) << std::arg(f) * 180.0 / constants::PI << "\n";
	}
}
