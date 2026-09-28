#include "pch.h"
#include "structure_factors.h"
#include "wfn_class.h"
#include "scattering_factors.h"
#include "basis_set.h"
#include "XCW.h"


structure_factors::structure_factors(const options& opt_in) {
	opt = &opt_in;

	// Read hkl and load cell
	std::filesystem::path hkl_filename = opt->hkl;
	std::filesystem::path cif = opt->cif;
	std::ifstream cif_input(cif.c_str(), std::ios::in);
	std::vector<asym_atom> xyz_atoms;
	if (!opt_in.xyz_file.empty()) {
		const std::filesystem::path xyz_path = opt_in.xyz_file;
		WFN dummy_wave(xyz_path, opt_in.debug);
		dummy_wave.read_xyz(xyz_path, std::cout, opt_in.debug);
		xyz_atoms = dummy_wave.extract_xyz("bohr");
	}
	unit_cell = cell(cif, std::cout, opt->debug, opt->do_XCW);
	std::vector<scattering_data> obs;
	scatter_data.hkl_enlarged = read_hkl_full(hkl_filename, scatter_data.hkl, opt_in.twin_law, unit_cell, std::cout, obs, opt_in.debug);
	std::ofstream log3("log3.txt", std::ios::out);
	for (int i = 0; i < obs.size(); i++) {
		scatter_data.F_obs.push_back(obs[i].F_obs);
		scatter_data.F_obs2.push_back(obs[i].F_obs2);
		scatter_data.sigma_obs.push_back(obs[i].sigma_obs);
		scatter_data.sigma_obs2.push_back(obs[i].sigma_obs2);
		scatter_data.abs_F_obs.push_back(obs[i].abs_F_obs);
	}
	bvec needs_grid;
	read_atoms_from_CIF(cif_input, unit_cell, model_data.ncen, needs_grid, asym_atoms, opt_in.debug);
	err_checkf(model_data.ncen > 0, "No atoms were read from " + cif.string() + "! Is there an _atom_site loop with labels, type symbols and fractional coordinates?", std::cout);

	// Adds symmetry generated atoms
	if (!opt_in.xyz_file.empty()) {
		unit_cell.grow_asym_atoms(asym_atoms, xyz_atoms);
	}

	// Evaluate symmetry and assign asymmetry factors to each atom (also update ncen)
	//The linking list is ordered like this: Asymmetric atom, list with all atoms, then index of symmetry operation that generated it
	// "diagonal elements" have to have size equivalent to multiplicity, otherwise something broke
	ivec3 symmetry_linking_list;
	unit_cell.eval_symm(asym_atoms, model_data.ncen, symmetry_linking_list);
	model_data.ncen = asym_atoms.size();

	// Warn if a grown structure's explicit atoms don't consistently cover the same
	// symmetry operations for every asymmetric atom
	ivec applied_symmetry;
	if (!opt_in.xyz_file.empty()) {
		// Below is working
		//unit_cell.apply_grown(symmetry_linking_list);
		applied_symmetry = unit_cell.apply_grown(scatter_data.hkl, scatter_data.hkl_enlarged, asym_atoms, symmetry_linking_list, original_rotations);
	}

	//unit_cell.set_symmetry_factors(asym_atoms, symmetry_linking_list);
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

	// Generate WFN object from asym_atoms
	for (int at = 0; at < model_data.ncen; at++) {
		asym_atom_list.push_back(at);
		atom temp_atom;
		temp_atom.set_coordinate(0, asym_atoms[at].pos[0]);
		temp_atom.set_coordinate(1, asym_atoms[at].pos[1]);
		temp_atom.set_coordinate(2, asym_atoms[at].pos[2]);
		temp_atom.set_charge(asym_atoms[at].type);
		dummy_wave.push_back_atom(temp_atom);
	}

	// Load basis set & generate basis for each atom
	std::shared_ptr<BasisSet> basis = BasisSetLibrary::get_basis_set(basis_set_name);
	load_basis_into_WFN(dummy_wave, basis, false, true);

	// Read isotropic displacement parameters
	vec U_iso = read_U_iso_from_CIF(cif, dummy_wave, unit_cell, log3, opt_in.debug);
	for (int i = 0; i < asym_atoms.size(); i++) {
		asym_atoms[i].U_iso = U_iso[i];
	}


	if (!opt_in.xyz_file.empty()) {
		unit_cell.grow_U_iso(asym_atoms, symmetry_linking_list);
	}

	// Generate k_pts and set the number of reflections
	model_data.nr_enlarged = scatter_data.hkl_enlarged.size();
	model_data.nr = scatter_data.hkl.size();
	make_k_pts(model_data.nr_enlarged != 0 && scatter_data.hkl.size() == 0, opt_in.save_k_pts, unit_cell, scatter_data.hkl_enlarged, k_pt, std::cout, opt_in.debug);

	// Read ADPs
	bool grown = !opt_in.xyz_file.empty();
	read_fracs_ADPs_from_CIF(cif, dummy_wave, log3, opt_in.debug, grown, symmetry_linking_list);

	// Prepare output files
	XCW_log.open("XCW.log");
	std::cout << "XCW orbital basis set: " << basis->get_name() << std::endl;
	XCW_log << "XCW orbital basis set: " << basis->get_name() << std::endl;

	// The fit set, see i_sigma_cutoff. F_obs2 is |I|, the sign lives in F_obs
	scatter_data.hkl_mask.resize(model_data.nr, 0);
	nr_fit = 0;
	for (int r = 0; r < model_data.nr; r++) {
		const double I_over_sigma = (scatter_data.F_obs[r] < 0 ? -scatter_data.F_obs2[r] : scatter_data.F_obs2[r]) / scatter_data.sigma_obs2[r];
		scatter_data.hkl_mask[r] = scatter_data.sigma_obs2[r] > 0 && I_over_sigma >= i_sigma_cutoff;
		nr_fit += scatter_data.hkl_mask[r];
	}
	setup_extinction(cif);
	err_checkf(nr_fit > n_params(), "Fewer reflections above the I/sigma cutoff than parameters", std::cout);
	std::cout << "XCW: I/sigma(I) >= " << i_sigma_cutoff << " (F/sigma(F) >= " << 2 * i_sigma_cutoff << "): " << nr_fit << " of " << model_data.nr << " reflections in the fit; R1 and Criterion are over these, R1(all) and Crit(all) over all" << std::endl;
	XCW_log << "XCW: I/sigma(I) >= " << i_sigma_cutoff << ": " << nr_fit << " of " << model_data.nr << " reflections in the fit" << std::endl;

	// Precompute GooF scaling factor
	inv_scale = 1.0 / (nr_fit - n_params());

	// Set F_calc sizes
	scatter_data.F_calc.resize(model_data.nr, 0);
	scatter_data.anom_correction.resize(model_data.nr, 0);
}

void structure_factors::setup_extinction(const std::filesystem::path& cif) {
	if (ext_model == extinction::model::none) return;
	const double lambda = wavelength > 0.0 ? wavelength : read_cif_wavelength(cif);
	err_checkf(lambda > 0.0, "Extinction needs a wavelength: put `wavelength <lambda>` in the XCW "
		"settings file, or _diffrn_radiation_wavelength in " + cif.string(), std::cout);
	ensure_hkl_ordered();
	const size_t np = extinction_aniso ? 6 : 1;
	//the anisotropic tensor starts isotropic, where x(h) is the start value for every h
	ext_p_.assign(np, 0.0);
	for (size_t p = 0; p < (extinction_aniso ? 3u : 1u); p++) ext_p_[p] = extinction_start;
	ext_c_.resize(model_data.nr);
	ext_cos2t_.resize(model_data.nr);
	if (extinction_aniso) ext_a_.resize(static_cast<size_t>(model_data.nr) * np);
	for (int r = 0; r < model_data.nr; r++) {
		const double stl = unit_cell.get_stl_of_hkl(hkl_ordered_[r]);
		ext_c_[r] = extinction::geometry_constant(lambda, stl);
		ext_cos2t_[r] = extinction::cos_2theta(lambda, stl);
		if (!extinction_aniso) continue;
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
	banner << "XCW extinction: " << extinction::name(ext_model)
		<< (extinction_aniso ? ", anisotropic (azimuth-averaged, 6 parameters)" : ", isotropic (1 parameter)")
		<< (extinction_refine ? ", refined with the scale" : ", held fixed")
		<< ", lambda = " << lambda << " A, start value " << extinction_start;
	std::cout << banner.str() << std::endl;
	XCW_log << banner.str() << std::endl;
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

double structure_factors::read_cif_wavelength(const std::filesystem::path& cif) {
	std::ifstream input(cif, std::ios::in);
	std::string line;
	while (input.good() && !input.eof()) {
		getline_universal(input, line);
		std::istringstream words(line);
		std::string tag, value;
		if (!(words >> tag) || tag != "_diffrn_radiation_wavelength") continue;
		if (!(words >> value)) continue;
		try { return std::stod(value.substr(0, value.find('('))); }
		catch (const std::exception&) { return 0.0; }
	}
	return 0.0;
}