#include "pch.h"
#include "XCW.h"
#if defined(NOSPHERA2_USE_GPU) || defined(NOSPHERA2_USE_METAL)
#include "itensor_gpu.h"
#endif
#include "convenience.h"
#include "scattering_factors.h"
#include "nos_math.h"
#include "basis_set.h"
#include "bondwise_analysis.h"
#include <mutex>
#include <random>

void XCW::construct(const options& opt_in) {
	opt = &opt_in;

	// Read hkl and load cell
	std::filesystem::path hkl_filename = opt->hkl;
	std::filesystem::path cif = opt->cif;
	std::ifstream cif_input(cif.c_str(), std::ios::in);
	std::optional<std::filesystem::path> xyz_path;
	std::optional<std::vector<asym_atom>> xyz_atoms;
	if (settings.grown) {
		if ((opt->xyz_file.empty())) {
			std::cerr << "I need an xyz file to grow the crystal, but none was provided. Exiting." << std::endl;
		}
		xyz_path = opt->xyz_file;
		WFN dummy_wave;
		dummy_wave.read_xyz(*xyz_path, std::cout, opt->debug);
		xyz_atoms = dummy_wave.extract_xyz("bohr");
	}
	unit_cell = cell(cif, std::cout, opt->debug, opt->do_XCW);
	hkl_enlarged = read_hkl_full(hkl_filename, hkl, opt->twin_law, unit_cell, std::cout, obs, opt->debug);
	std::ofstream log3("log3.txt", std::ios::out);
	bvec needs_grid;
	read_atoms_from_CIF(cif_input, unit_cell, cryst.ncen, needs_grid, asym_atoms, opt->debug);
	err_checkf(cryst.ncen > 0, "No atoms were read from " + cif.string() + "! Is there an _atom_site loop with labels, type symbols and fractional coordinates?", std::cout);

	// Adds symmetry generated atoms
	if (settings.grown) {
		unit_cell.grow_asym_atoms(asym_atoms, xyz_atoms.value());
	}

	// Evaluate symmetry and assign asymmetry factors to each atom (also update ncen)
	//The linking list is ordered like this: Asymmetric atom, list with all atoms, then index of symmetry operation that generated it
	// "diagonal elements" have to have size equivalent to multiplicity, otherwise something broke
	ivec3 symmetry_linking_list;
	unit_cell.eval_symm(asym_atoms, cryst.ncen, symmetry_linking_list);
	cryst.ncen = asym_atoms.size();

	// Warn if a grown structure's explicit atoms don't consistently cover the same
	// symmetry operations for every asymmetric atom
	ivec applied_symmetry;
	if (settings.grown) {
		// Below is working
		//unit_cell.apply_grown(symmetry_linking_list);
		applied_symmetry = unit_cell.apply_grown(hkl, hkl_enlarged, asym_atoms, symmetry_linking_list, original_rotations);
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
		std::cerr << "hkl_enlarged size: " << hkl_enlarged.size() << std::endl;
	}

	// Below is working
	// Structure factors sum over every operation, unless the grown cluster is a union of complete
	// orbits of a subgroup H: then one operation per coset of H covers the cell with |H| times fewer
	// terms and the cluster's own symmetry is not applied a second time
	//sym_ops_.resize(unit_cell.get_trans()[0].size());
	//std::iota(sym_ops_.begin(), sym_ops_.end(), 0);
	//if (settings.grown) {
	//	const ivec subgroup = unit_cell.grown_subgroup(symmetry_linking_list);
	//	if (subgroup.size() < 2)
	//		std::cout << "XCW: grown cluster is mapped onto itself by no symmetry operation, summing all " << sym_ops_.size() << " operations" << std::endl;
	//	else {
	//		sym_ops_ = unit_cell.coset_representatives(subgroup);
	//		unit_cell.set_subgroup_factors(asym_atoms, symmetry_linking_list, subgroup);
	//		std::cout << "XCW: grown cluster is mapped onto itself by a subgroup of order " << subgroup.size() << ", summing "
	//			<< sym_ops_.size() << " coset representatives instead of " << unit_cell.get_trans()[0].size() << " operations" << std::endl;
	//	}
	//}

	// Generate WFN object from asym_atoms
	dummy_wave.assign_charge(settings.charge);
	dummy_wave.assign_multi(settings.multiplicity);
	for (int at = 0; at < cryst.ncen; at++) {
		asym_atom_list.push_back(at);
		atom temp_atom;
		temp_atom.set_coordinate(0, asym_atoms[at].pos[0]);
		temp_atom.set_coordinate(1, asym_atoms[at].pos[1]);
		temp_atom.set_coordinate(2, asym_atoms[at].pos[2]);
		temp_atom.set_charge(asym_atoms[at].type);
		dummy_wave.push_back_atom(temp_atom);
	}

	// Load basis set & generate basis for each atom
	std::shared_ptr<BasisSet> basis = BasisSetLibrary::get_basis_set(settings.basis_set_name);
	load_basis_into_WFN(dummy_wave, basis, false, true);

	// Read isotropic displacement parameters
	cryst.U_iso = read_U_iso_from_CIF(cif, dummy_wave, unit_cell, log3, opt->debug);
	if (settings.grown) {
		cryst.grow_U_iso(asym_atoms, symmetry_linking_list);
	}

	// Generate k_pts and set the number of reflections
	cryst.nr = hkl_enlarged.size();
	cryst.nr_small = hkl.size();
	make_k_pts(cryst.nr != 0 && hkl.size() == 0, opt->save_k_pts, unit_cell, hkl_enlarged, k_pt, std::cout, opt->debug);

	// Read ADPs
	read_fracs_ADPs_from_CIF(cif, dummy_wave, log3, opt->debug, settings.grown, symmetry_linking_list);

	// Prepare output files
	XCW_log.open("XCW.log");
	std::cout << "XCW orbital basis set: " << basis->get_name() << std::endl;
	XCW_log << "XCW orbital basis set: " << basis->get_name() << std::endl;

	// The fit set, see i_sigma_cutoff. F_obs2 is |I|, the sign lives in F_obs
	fit_mask_.assign(cryst.nr_small, false);
	cryst.n_fit = 0;
	for (int r = 0; r < cryst.nr_small; r++) {
		const double I_over_sigma = (obs[r].F_obs < 0 ? -obs[r].F_obs2 : obs[r].F_obs2) / obs[r].sigma_obs2;
		fit_mask_[r] = obs[r].sigma_obs2 > 0 && I_over_sigma >= settings.i_sigma_cutoff;
		cryst.n_fit += fit_mask_[r];
	}
	err_checkf(cryst.n_fit > settings.n_params, "Fewer reflections above the I/sigma cutoff than parameters", std::cout);
	std::cout << "XCW: I/sigma(I) >= " << settings.i_sigma_cutoff << " (F/sigma(F) >= " << 2 * settings.i_sigma_cutoff << "): " << cryst.n_fit << " of " << cryst.nr_small << " reflections in the fit; R1 and Criterion are over these, R1(all) and Crit(all) over all" << std::endl;
	XCW_log << "XCW: I/sigma(I) >= " << settings.i_sigma_cutoff << ": " << cryst.n_fit << " of " << cryst.nr_small << " reflections in the fit" << std::endl;

	// Precompute GooF scaling factor
	cryst.inv_scale = 1.0 / (cryst.n_fit - settings.n_params);

	// Set F_calc sizes
	F_calc.resize(2);
	F_calc[0].resize(cryst.nr_small, 0);
	F_calc[1].resize(cryst.nr_small, 0);
}

XCW::SCF_settings XCW::loadSettings(const std::filesystem::path& settings_path) {
	auto lowercase = [](std::string s) {
		std::transform(s.begin(), s.end(), s.begin(),
			[](unsigned char c) { return std::tolower(c); });
		return s;
		};

	SCF_settings settings;

	settings.grown = false;
	std::string conv_preset = "normal";
	bool soscf = false, check_hessian = false;
	std::string speed_preset = "normal_conv";

	/* 1: Classical Jayatilaka XWR
	   2: Ewald sum XWR */
	settings.XWR_type = 1;
	double quant_diff = 32768, diis_stop_damping = 32768, diis_stop_shift = 32768, max_diis_error = 32768, gradient = 32768, MaxP_diff = 32768, RMSP_diff = 32768, alpha = 32768, level_shift = 32768, start = 32768, end = 32768, step_size = 32768;
	int max_scf_iterations = 32768, charge = 32768, multiplicity = 32768, n_params = 32768, refine_against = 32768;
	std::string basis_set_name = "Undefined";
	std::string df_basis_name, guess_basis_name;
	bool grown = false, read_tensor = false, read_first_guess = false, nbo_output = false;
	bool i_tensor_single = false, i_tensor_double = false;
	std::filesystem::path i_tensor_file_path;
	std::filesystem::path i_tensor_save_path;
	// 0 = hold the whole tensor, which is what every run did before this existed.
	size_t i_tensor_max_mb = 0;
	occ::qm::SpinorbitalKind hf_type = occ::qm::SpinorbitalKind::Restricted;
	if (!std::filesystem::exists(settings_path)) {
		throw std::runtime_error("Settings file not found! Aborting run!");
	}
	else {

		std::ifstream input(settings_path);
		using Handler = std::function<void(std::istream&)>;
		std::unordered_map<std::string, Handler> handlers;
		handlers["conv"] = [&](std::istream& is) {
			if (!(is >> quant_diff))
				throw std::runtime_error("Expected value after 'conv'");
			};

		handlers["diis_damping"] = [&](std::istream& is) {
			if (!(is >> diis_stop_damping))
				throw std::runtime_error("Expected value after 'diis_damping'");
			};

		handlers["diis_shift"] = [&](std::istream& is) {
			if (!(is >> diis_stop_shift))
				throw std::runtime_error("Expected value after 'diis_shift'");
			};

		handlers["conv_diis"] = [&](std::istream& is) {
			if (!(is >> max_diis_error))
				throw std::runtime_error("Expected value after 'conv_diis'");
			};

		handlers["gradient"] = [&](std::istream& is) {
			if (!(is >> gradient))
				throw std::runtime_error("Expected value after 'gradient'");
			};

		handlers["maxp_diff"] = [&](std::istream& is) {
			if (!(is >> MaxP_diff))
				throw std::runtime_error("Expected value after 'MaxP_diff'");
			};

		handlers["rmsp_diff"] = [&](std::istream& is) {
			if (!(is >> RMSP_diff))
				throw std::runtime_error("Expected value after 'RMSP_diff'");
			};

		handlers["params"] = [&](std::istream& is) {
			if (!(is >> n_params))
				throw std::runtime_error("Expected value after 'params'");
			};

		handlers["damp"] = [&](std::istream& is) {
			if (!(is >> alpha))
				throw std::runtime_error("Expected value after 'damp'");
			};

		handlers["shift"] = [&](std::istream& is) {
			if (!(is >> level_shift))
				throw std::runtime_error("Expected value after 'shift'");
			};

		handlers["max_iter"] = [&](std::istream& is) {
			if (!(is >> max_scf_iterations))
				throw std::runtime_error("Expected value after 'max_iter'");
			};

		handlers["charge"] = [&](std::istream& is) {
			if (!(is >> charge))
				throw std::runtime_error("Expected value after 'charge'");
			};

		handlers["mult"] = [&](std::istream& is) {
			if (!(is >> multiplicity))
				throw std::runtime_error("Expected value after 'mult'");
			};

		handlers["f"] = [&](std::istream&) {
			refine_against = 1;
			};

		handlers["f2"] = [&](std::istream&) {
			refine_against = 2;
			};

		handlers["i_sigma"] = [&](std::istream& is) {
			if (!(is >> settings.i_sigma_cutoff))
				throw std::runtime_error("Expected value after 'i_sigma'");
			};

		handlers["weighted"] = [&](std::istream&) {
			settings.XWR_type = 2;
			};

		handlers["basis_set"] = [&](std::istream& is) {
			if (!(is >> basis_set_name))
				throw std::runtime_error("Expected basis set name");
			};
		handlers["df_basis"] = [&](std::istream& is) {
			if (!(is >> df_basis_name))
				throw std::runtime_error("Expected a fitting basis name after 'df_basis'");
			};
		handlers["guess_basis"] = [&](std::istream& is) {
			if (!(is >> guess_basis_name))
				throw std::runtime_error("Expected a basis name after 'guess_basis'");
			};

		handlers["start"] = [&](std::istream& is) {
			if (!(is >> start))
				throw std::runtime_error("Expected value after 'start'");
			};

		handlers["end"] = [&](std::istream& is) {
			if (!(is >> end))
				throw std::runtime_error("Expected value after 'end'");
			};

		handlers["step_size"] = [&](std::istream& is) {
			if (!(is >> step_size))
				throw std::runtime_error("Expected value after 'step_size'");
			};

		handlers["grown"] = [&](std::istream&) {
			grown = true;
			};

		handlers["rhf"] = [&](std::istream&) {
			hf_type = occ::qm::SpinorbitalKind::Restricted;
			};

		handlers["uhf"] = [&](std::istream&) {
			hf_type = occ::qm::SpinorbitalKind::Unrestricted;
			};

		handlers["sloppy"] = [&](std::istream&) {
			conv_preset = "sloppy";
			};

		handlers["normal"] = [&](std::istream&) {
			conv_preset = "normal";
			};

		handlers["tight"] = [&](std::istream&) {
			conv_preset = "tight";
			};

		handlers["very_tight"] = [&](std::istream&) {
			conv_preset = "very_tight";
			};

		handlers["slow_conv"] = [&](std::istream&) {
			speed_preset = "slow_conv";
			};

		handlers["normal_conv"] = [&](std::istream&) {
			speed_preset = "normal_conv";
			};

		handlers["soscf"] = [&](std::istream&) {
			soscf = true;
			};

		handlers["check_hessian"] = [&](std::istream&) {
			check_hessian = true;
			};

		handlers["fast_conv"] = [&](std::istream&) {
			speed_preset = "fast_conv";
			};

		//`safe` is `save` to the default path, which is where `read` without a path looks.
		handlers["safe"] = [&](std::istream&) {
			i_tensor_save_path = i_tensor_default;
			};

		//`read <path>` reuses the streamed tensor at that path, which is the expensive thing
		//a run produces and which depends only on the geometry, the basis and the
		//reflections, not on any refinement setting. So trying another lambda range or
		//convergence preset need not rebuild it. It is held in memory when it fits.
		//
		//The path is optional, and the settings file is one whitespace-separated stream of
		//tokens, so a bare `read` followed by another keyword must not swallow it: take the
		//next token, put it back if it is a keyword.
		handlers["read"] = [&](std::istream& is) {
			read_tensor = true;
			i_tensor_file_path = i_tensor_default;
			const std::streampos before = is.tellg();
			std::string token;
			if (!(is >> token)) { is.clear(); return; }
			std::string lowered = token;
			std::transform(lowered.begin(), lowered.end(), lowered.begin(),
				[](unsigned char c) { return std::tolower(c); });
			if (handlers.find(lowered) != handlers.end()) {
				is.clear();
				is.seekg(before);
				return;
			}
			i_tensor_file_path = token;
			};
		handlers["load_wfn"] = [&](std::istream&) {
			read_first_guess = true;
			};
		handlers["nbo"] = [&](std::istream&) {
			nbo_output = true;
			};
		// The tensor is nr_small blocks of nmo(nmo+1)/2 complex doubles and grows
		// quadratically with the basis, so on anything past a minimal basis it is
		// the largest thing in the process. `stream` puts it on disk with a
		// default budget; `i_tensor_mb <n>` names the budget.
		handlers["stream"] = [&](std::istream&) {
			if (i_tensor_max_mb == 0) i_tensor_max_mb = 2048;
			};

		//The tensor comes off the device in single precision and was stored in double, so
		//half of every byte the SCF loop reads back was padding. Holding it as computed is
		//half the memory and, measured on the full twisted ethylene, 1.63x on the two walks
		//each SCF iteration makes over it. Opt-in because it is a change of stored
		//precision: the lambda scan agrees to 5e-13 and iteration for iteration, but that
		//is a measurement on one system rather than a proof.
		//Writing the tensor is worth ~40 minutes to a later run and costs this one nothing:
		//the refinement only reads the tensor, so a thread can push it to disk while the SCF
		//gets on with it. The path is what `read <path>` will want afterwards.
		handlers["save"] = [&](std::istream& is) {
			std::string path;
			if (!(is >> path))
				throw std::runtime_error("Expected a path after 'save'");
			i_tensor_save_path = path;
			};

		handlers["i_float"] = [&](std::istream&) {
			i_tensor_single = true;
			};

		handlers["i_double"] = [&](std::istream&) {
			i_tensor_double = true;
			};

		handlers["i_tensor_mb"] = [&](std::istream& in2) {
			long long mb = 0;
			in2 >> mb;
			i_tensor_max_mb = (mb > 0) ? static_cast<size_t>(mb) : 0;
			};

		std::string keyword;

		while (input >> keyword)
		{
			std::transform(keyword.begin(), keyword.end(), keyword.begin(),
				[](unsigned char c) { return std::tolower(c); });

			auto it = handlers.find(keyword);
			if (it == handlers.end()) {
				throw std::runtime_error("Unknown keyword '" + keyword + "'");
			}

			it->second(input);
		}
	}

	if (conv_preset == "sloppy") {
		settings.quant_diff = 3e-5;
		settings.max_diis_error = 1e-4;
		settings.gradient = 5e-4;
		settings.MaxP_diff = 1e-4;
		settings.RMSP_diff = 1e-5;
		settings.max_scf_iterations = 100;
	}
	else if (conv_preset == "normal") {
		settings.quant_diff = 1e-6;
		settings.max_diis_error = 1e-5;
		settings.gradient = 7e-5;
		settings.MaxP_diff = 1e-5;
		settings.RMSP_diff = 1e-6;
		settings.max_scf_iterations = 100;
	}
	else if (conv_preset == "tight") {
		settings.quant_diff = 5e-7;
		settings.max_diis_error = 5e-6;
		settings.gradient = 3e-5;
		settings.MaxP_diff = 1e-6;
		settings.RMSP_diff = 1e-7;
		settings.max_scf_iterations = 100;
	}
	else if (conv_preset == "very_tight") {
		settings.quant_diff = 1e-7;
		settings.max_diis_error = 1e-6;
		settings.gradient = 1e-5;
		settings.MaxP_diff = 1e-7;
		settings.RMSP_diff = 1e-8;
		settings.max_scf_iterations = 100;
	}

	if (speed_preset == "slow_conv") {
		settings.slow_conv = true;
		settings.alpha = 0.8;
		settings.level_shift = 1;
		settings.diis_stop_damping = 1e-5;
		settings.diis_stop_shift = 1e-5;
	}
	else if (speed_preset == "normal_conv") {
		settings.alpha = 0.5;
		settings.level_shift = 0.5;
		settings.diis_stop_damping = 1e-3;
		settings.diis_stop_shift = 1e-2;
	}
	else if (speed_preset == "fast_conv") {
		settings.method_apply_damping = false;
		settings.method_apply_shift = false;
	}

	if (basis_set_name == "Undefined") {
		throw std::runtime_error("Basis set name not specified in settings file! Aborting run!");
	}
	settings.basis_set_name = basis_set_name;
	if (quant_diff != 32768) settings.quant_diff = quant_diff;
	if (diis_stop_damping != 32768) settings.diis_stop_damping = diis_stop_damping;
	if (diis_stop_shift != 32768) settings.diis_stop_shift = diis_stop_shift;
	if (max_diis_error != 32768) settings.max_diis_error = max_diis_error;
	if (gradient != 32768) settings.gradient = gradient;
	if (MaxP_diff != 32768) settings.MaxP_diff = MaxP_diff;
	if (RMSP_diff != 32768) settings.RMSP_diff = RMSP_diff;
	if (alpha != 32768) settings.alpha = alpha;
	if (level_shift != 32768) settings.level_shift = level_shift;
	if (max_scf_iterations != 32768) settings.max_scf_iterations = max_scf_iterations;
	settings.charge = (charge != 32768) ? charge : settings.charge;
	settings.multiplicity = (multiplicity != 32768) ? multiplicity : settings.multiplicity;
	if (n_params != 32768) settings.n_params = n_params;
	if (refine_against != 32768) settings.refine_against = refine_against;
	settings.xcw_start_value = (start != 32768) ? start : 0;
	settings.xcw_step_size = (step_size != 32768) ? step_size : 0.01;
	if (end != 32768) {
		settings.num_xcw_steps = static_cast<int>(std::round((end - settings.xcw_start_value) / settings.xcw_step_size)) + 1;
	}
	else {
		settings.num_xcw_steps = static_cast<int>(std::round((1 - settings.xcw_start_value) / settings.xcw_step_size)) + 1;
	}
	settings.grown = grown;
	settings.hf_type = hf_type;
	settings.read_tensor = read_tensor;
	settings.read_first_guess = read_first_guess;
	settings.i_tensor_max_mb = i_tensor_max_mb;
	settings.i_tensor_single = i_tensor_single;
	settings.i_tensor_double = i_tensor_double;
	settings.i_tensor_file_path = i_tensor_file_path;
	settings.i_tensor_save_path = i_tensor_save_path;
	settings.nbo_output = nbo_output;
	settings.df_basis_name = df_basis_name;
	settings.guess_basis_name = guess_basis_name;
	settings.soscf = soscf;
	settings.check_hessian = check_hessian;

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

// Position of a sorted index triple/quadruple in the Voigt storage of C (10) and D (15)
static void get_voigt_index(const ivec& indices, int& ADP_idx) {
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

// T'_{ij..} = sum M_pi M_qj .. T_pq.. for the U (rank 2), C (rank 3) and D (rank 4) tensors in their
// Voigt storage: U_star2U_cart hands in the cell matrix, the grown-atom rotation the transposed symmetry operation
void transform_ADPs(vec2& ADPs, const vec2& M) {
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
								sum += M[p][i] * M[q][j] * M[r][k] * ADPs[1][ADP_idx];
							}
						}
					}
					C_out[running_idx] = sum;
					running_idx++;
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
										get_voigt_index(sorted_idx, ADP_idx);
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

void XCW::U_star2U_cart() {
	const double scale = constants::bohr2ang(1);
	vec2 cart_matrix(3, vec(3));
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			cart_matrix[i][j] = unit_cell.get_cm(i, j) * scale;
		}
	}
	for (int a = 0; a < cryst.ncen; a++) {
		vec2 ADPs = dummy_wave.get_atom(a).get_ADPs();
		transform_ADPs(ADPs, cart_matrix);
		dummy_wave.set_atom_ADPs(a, ADPs);
	}
}

// A grown atom carries a copy of its parent's ADPs (read_fracs_ADPs_from_CIF); U*, C and D are
// contravariant tensors in the fractional basis, so the image's are T' = R T R^T with R the rotation
// of the linking operation (x' = R x + t). cell stores R transposed, which is what transform_ADPs takes.
void XCW::rotate_grown_ADPs() {
	for (int a = 0; a < cryst.ncen; a++) {
		const int op = asym_atoms[a].sym_op;
		if (op < 0) continue;
		vec2 M(3, vec(3));
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				M[i][j] = original_rotations[i][j][op];
				//M[i][j] = unit_cell.get_sym(i, j, op);
			}
		}
		vec2 ADPs = dummy_wave.get_atom(a).get_ADPs();
		transform_ADPs(ADPs, M);
		dummy_wave.set_atom_ADPs(a, ADPs);
	}
}

void XCW::eval_DW(cvec2& DW_fact) {
	DW_fact.resize(cryst.ncen, cvec(cryst.nr, 0));
	//Converts angstrom to bohr OR MORE IMPORTANTLY reciprocal bohr to reciprocal angstrom
	const double angstrom2bohr = constants::ang2bohr(1);
	ivec level;
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
	rotate_grown_ADPs();
	U_star2U_cart();
	vec2 q(cryst.nr, vec(3));
	for (int h = 0; h < cryst.nr; h++) {
		q[h][0] = k_pt[0][h];
		q[h][1] = k_pt[1][h];
		q[h][2] = k_pt[2][h];
	}
	std::transform(q.begin(), q.end(), q.begin(), [angstrom2bohr](std::vector<double>& vec) {
		std::transform(vec.begin(), vec.end(), vec.begin(), [angstrom2bohr](double x) { return x * angstrom2bohr; });
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
				temp = -0.5 * U * (q[r][0] * q[r][0] + q[r][1] * q[r][1] + q[r][2] * q[r][2]);
				DW_fact[a][r] = std::exp(temp);
			}
			break;
		}
		case 1: {
			// Anisotropic U_ij
			double temp1;
			for (int h = 0; h < cryst.nr; h++) {
				vec q_ = { q[h][0], q[h][1], q[h][2] };
				temp1 = -0.5 * dot_BLAS(dot(Uij, q_, true), q_, false);
				DW_fact[a][h] = std::exp(temp1);
			}
			break;
		}
		case 2: {
			// Anisotropic C_ijk
			double temp1, temp2;
			for (int h = 0; h < cryst.nr; h++) {
				vec q_ = { q[h][0], q[h][1], q[h][2] };
				//temp2 = -1.0 / 6.0 * (constants::TWO_PI * constants::TWO_PI * constants::TWO_PI) * (ADPs[1][0] * q_[0] * q_[0] * q_[0] + ADPs[1][6] * q_[1] * q_[1] * q_[1] + ADPs[1][9] * q_[2] * q_[2] * q_[2]
				//	+ 3 * ADPs[1][1] * q_[0] * q_[0] * q_[1] + 3 * ADPs[1][2] * q_[0] * q_[0] * q_[2] + 3 * ADPs[1][3] * q_[0] * q_[1] * q_[1] + 3 * ADPs[1][5] * q_[0] * q_[2] * q_[2] + 3 * ADPs[1][7] * q_[1] * q_[1] * q_[2] + 3 * ADPs[1][8] * q_[1] * q_[2] * q_[2]
				//	+ 6 * ADPs[1][4] * q_[0] * q_[1] * q_[2]);
				temp1 = -0.5 * dot_BLAS(dot(Uij, q_, true), q_, false);
				temp2 = -1.0 / 6.0 * (ADPs[1][0] * q_[0] * q_[0] * q_[0] + ADPs[1][6] * q_[1] * q_[1] * q_[1] + ADPs[1][9] * q_[2] * q_[2] * q_[2]
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
				//temp2 = -1.0 / 6.0 * (constants::TWO_PI * constants::TWO_PI * constants::TWO_PI) * (ADPs[1][0] * q_[0] * q_[0] * q_[0] + ADPs[1][6] * q_[1] * q_[1] * q_[1] + ADPs[1][9] * q_[2] * q_[2] * q_[2]
				//	+ 3 * ADPs[1][1] * q_[0] * q_[0] * q_[1] + 3 * ADPs[1][2] * q_[0] * q_[0] * q_[2] + 3 * ADPs[1][3] * q_[0] * q_[1] * q_[1] + 3 * ADPs[1][5] * q_[0] * q_[2] * q_[2] + 3 * ADPs[1][7] * q_[1] * q_[1] * q_[2] + 3 * ADPs[1][8] * q_[1] * q_[2] * q_[2]
				//	+ 6 * ADPs[1][4] * q_[0] * q_[1] * q_[2]);
				//temp3 = (1.0 / 24.0) * (constants::TWO_PI * constants::TWO_PI * constants::TWO_PI * constants::TWO_PI) * (ADPs[2][0] * q_[0] * q_[0] * q_[0] * q_[0] + 4.0 * ADPs[2][1] * q_[0] * q_[0] * q_[0] * q_[1] + 4.0 * ADPs[2][2] * q_[0] * q_[0] * q_[0] * q_[2]
				//	+ 6.0 * ADPs[2][3] * q_[0] * q_[0] * q_[1] * q_[1] + 12.0 * ADPs[2][4] * q_[0] * q_[0] * q_[1] * q_[2] + 6.0 * ADPs[2][5] * q_[0] * q_[0] * q_[2] * q_[2] + 4.0 * ADPs[2][6] * q_[0] * q_[1] * q_[1] * q_[1] + 12.0 * ADPs[2][7] * q_[0] * q_[1] * q_[1] * q_[2]
				//	+ 12.0 * ADPs[2][8] * q_[0] * q_[1] * q_[2] * q_[2] + 4.0 * ADPs[2][9] * q_[0] * q_[2] * q_[2] * q_[2] + ADPs[2][10] * q_[1] * q_[1] * q_[1] * q_[1] + 4.0 * ADPs[2][11] * q_[1] * q_[1] * q_[1] * q_[2] + 6.0 * ADPs[2][12] * q_[1] * q_[1] * q_[2] * q_[2]
				//	+ 4.0 * ADPs[2][13] * q_[1] * q_[2] * q_[2] * q_[2] + ADPs[2][14] * q_[2] * q_[2] * q_[2] * q_[2]);
				temp1 = -0.5 * dot_BLAS(dot(Uij, q_, true), q_, false);
				temp2 = -1.0 / 6.0 * (ADPs[1][0] * q_[0] * q_[0] * q_[0] + ADPs[1][6] * q_[1] * q_[1] * q_[1] + ADPs[1][9] * q_[2] * q_[2] * q_[2]
					+ 3 * ADPs[1][1] * q_[0] * q_[0] * q_[1] + 3 * ADPs[1][2] * q_[0] * q_[0] * q_[2] + 3 * ADPs[1][3] * q_[0] * q_[1] * q_[1] + 3 * ADPs[1][5] * q_[0] * q_[2] * q_[2] + 3 * ADPs[1][7] * q_[1] * q_[1] * q_[2] + 3 * ADPs[1][8] * q_[1] * q_[2] * q_[2]
					+ 6 * ADPs[1][4] * q_[0] * q_[1] * q_[2]);
				temp3 = (1.0 / 24.0) * (ADPs[2][0] * q_[0] * q_[0] * q_[0] * q_[0] + 4.0 * ADPs[2][1] * q_[0] * q_[0] * q_[0] * q_[1] + 4.0 * ADPs[2][2] * q_[0] * q_[0] * q_[0] * q_[2]
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
	cdouble exponent;
	for (int at = 0; at < cryst.ncen; at++) {
		vec pos_cart = { asym_atoms[at].pos[0], asym_atoms[at].pos[1], asym_atoms[at].pos[2] };
		for (int r = 0; r < cryst.nr; r++) {
			vec q = { k_pt[0][r], k_pt[1][r], k_pt[2][r] };
			exponent = cdouble(0, dot_BLAS(q, pos_cart, false));
			phase_fact[at][r] = std::exp(exponent);
		}
	}
}

// Below is working
//void XCW::eval_translation_phase(cvec2& translation_phase) {
//	translation_phase.resize(cryst.nr_small, cvec(sym_ops_.size(), 0));
//	const double angstrom2bohr = constants::ang2bohr(1);
//	const double bohr2angstrom = constants::bohr2ang(1);
//	vec2 trans = unit_cell.get_trans();
//	vec2 cm = { { unit_cell.get_cm(0,0), unit_cell.get_cm(0,1), unit_cell.get_cm(0,2)},
//								  { unit_cell.get_cm(1,0), unit_cell.get_cm(1,1), unit_cell.get_cm(1,2)},
//								  { unit_cell.get_cm(2,0), unit_cell.get_cm(2,1), unit_cell.get_cm(2,2)} };
//	std::transform(cm.begin(), cm.end(), cm.begin(), [bohr2angstrom](std::vector<double>& vec) {
//		std::transform(vec.begin(), vec.end(), vec.begin(), [bohr2angstrom](double x) { return x * bohr2angstrom; });
//		return vec; });
//	for (int r = 0; r < cryst.nr_small; r++) {
//		ivec asym_list = generate_asym_lookup(r);
//		vec q_temp = { k_pt[0][asym_list[0]], k_pt[1][asym_list[0]], k_pt[2][asym_list[0]] };
//		std::transform(q_temp.begin(), q_temp.end(), q_temp.begin(), [angstrom2bohr](double x) { return x * angstrom2bohr; });
//		for (int t = 0; t < sym_ops_.size(); t++) {
//			const int op = sym_ops_[t];
//			vec trans_temp = { trans[0][op], trans[1][op], trans[2][op] };
//			trans_temp = dot(cm, trans_temp, true);
//			cdouble exponent(0, dot_BLAS(q_temp, trans_temp, false));
//			translation_phase[r][t] = std::exp(exponent);
//		}
//	}
//	// closing function
//}

// sym_ops_ is the opp
void XCW::eval_translation_phase(cvec2& translation_phase) {
	translation_phase.resize(cryst.nr_small, cvec(unit_cell.get_trans()[0].size(), 0));
	const double angstrom2bohr = constants::ang2bohr(1);
	const double bohr2angstrom = constants::bohr2ang(1);
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
		std::transform(q_temp.begin(), q_temp.end(), q_temp.begin(), [angstrom2bohr](double x) { return x * angstrom2bohr; });
		for (int t = 0; t < trans[0].size(); t++) {
			vec trans_temp = { trans[0][t], trans[1][t], trans[2][t] };
			trans_temp = dot(cm, trans_temp, true);
			cdouble exponent(0, dot_BLAS(q_temp, trans_temp, false));
			translation_phase[r][t] = std::exp(exponent);
		}
	}
	// closing function
}

void XCW::parse_anom_atoms(std::vector<anom_atom>& anom_atoms) {
	std::ifstream file(opt->anom_disp_path);
	if (!file) {
		std::cout << "Could not open anomalous dispersion file. Continuing without anomalous dispersions." << std::endl;
	}
	std::string line;
	while (getline_universal(file, line)) {
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
	ivec2 asym_lookup(cryst.nr_small);
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

//The scale that minimises the criterion the SCF descends: k over the fit set from
//Sum w (k|Fc| - |Fo|)^2 / sigma^2, or k^2 from Sum w (k^2|Fc|^2 - Fo^2)^2 / sigma(I)^2 against
//F^2, with w the 1/|H|^2 weights of XWR_type 2. calc_perturb takes the scale as given, so the
//scale had better be stationary for the criterion, or the SCF descends a different
//functional than the one it prints: an unweighted fit of k put Fe(phen)2(SCN)2 (chi^2 12.96
//vs 12.30 at the weighted k, dchi^2/dk 2e3) 4 mEh above the previous step's orbitals at
//every converged lambda >= 0.03.
void XCW::eval_scale() {
	ensure_inv_H2_weights();
	const bool against_F2 = settings.refine_against == 2, weighted = settings.XWR_type == 2;
	const int chunk = 128, nchunk = (cryst.nr_small + chunk - 1) / chunk;
	vec numerators(nchunk), denominators(nchunk);
#pragma omp parallel for schedule(static)
	for (int c = 0; c < nchunk; c++) {
		const int first = c * chunk, last = std::min(first + chunk, cryst.nr_small);
		for (int i = first; i < last; i++) {
			if (!fit_mask_[i]) continue;
			const double w = weighted ? inv_H2_[i] : 1.0;
			const double calc = std::abs(F_calc[0][i]);
			if (against_F2) {
				const double calc2 = calc * calc, wi = w / (obs[i].sigma_obs2 * obs[i].sigma_obs2);
				numerators[c] += wi * calc2 * obs[i].F_obs2;
				denominators[c] += wi * calc2 * calc2;
			}
			else {
				const double wi = w / (obs[i].sigma_obs * obs[i].sigma_obs);
				numerators[c] += wi * calc * obs[i].F_obs;
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
	cryst.F_scale = against_F2 ? std::sqrt(std::max(ratio, 0.0)) : ratio;
}

void XCW::calc_criteria() {
	ensure_inv_H2_weights();
	//index 0: the fit set, 1: all reflections
	const double prefactor[2] = { 1.0 / static_cast<double>(cryst.n_fit - settings.n_params), 1.0 / static_cast<double>(cryst.nr_small - settings.n_params) };
	const int chunk = 128, nchunk = (cryst.nr_small + chunk - 1) / chunk;
	vec2 goof1(2, vec(nchunk)), goof2(2, vec(nchunk)), wgoof1(2, vec(nchunk)), wgoof2(2, vec(nchunk)), r1_num(2, vec(nchunk)), r1_den(2, vec(nchunk));
	const double scale = cryst.F_scale;
	const cdouble* F_calc_0 = F_calc[0].data();
#pragma omp parallel for schedule(static)
	for (int c = 0; c < nchunk; c++) {
		const int first = c * chunk, last = std::min(first + chunk, cryst.nr_small);
		for (int i = first; i < last; i++) {
			const scattering_data& obs_ptr = obs[i];
			const double scaled_F_calc = scale * std::abs(F_calc_0[i]);
			const double scaled_difference = scaled_F_calc - obs_ptr.F_obs;
			const double diff2 = (scaled_F_calc * scaled_F_calc) - obs_ptr.F_obs2;
			const double weighted_diff1 = scaled_difference / obs_ptr.sigma_obs;
			const double weighted_diff2 = diff2 / obs_ptr.sigma_obs2;
			const double weighted_diff1_sq = weighted_diff1 * weighted_diff1;
			const double weighted_diff2_sq = weighted_diff2 * weighted_diff2;
			const double w = settings.XWR_type == 2 ? inv_H2_[i] : 0.0;
			for (int set = fit_mask_[i] ? 0 : 1; set < 2; set++) {
				goof1[set][c] += weighted_diff1_sq;
				goof2[set][c] += weighted_diff2_sq;
				wgoof1[set][c] += weighted_diff1_sq * w;
				wgoof2[set][c] += weighted_diff2_sq * w;
				r1_num[set][c] += std::abs(scaled_F_calc - obs_ptr.abs_F_obs);
				r1_den[set][c] += obs_ptr.abs_F_obs;
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
	cryst.R1 = sum[5][0] > 0.0 ? sum[4][0] / sum[5][0] : 0.0;
	cryst.R1_all = sum[5][1] > 0.0 ? sum[4][1] / sum[5][1] : 0.0;
	cryst.GooF1 = std::sqrt(prefactor[0] * sum[0][0]);
	cryst.GooF2 = std::sqrt(prefactor[0] * sum[1][0]);
	cryst.weighted_GooF1 = std::sqrt(prefactor[0] * sum[2][0]);
	cryst.weighted_GooF2 = std::sqrt(prefactor[0] * sum[3][0]);
	cryst.GooF1_all = std::sqrt(prefactor[1] * sum[0][1]);
	cryst.GooF2_all = std::sqrt(prefactor[1] * sum[1][1]);
	cryst.weighted_GooF1_all = std::sqrt(prefactor[1] * sum[2][1]);
	cryst.weighted_GooF2_all = std::sqrt(prefactor[1] * sum[3][1]);
}

double XCW::criterion(const bool all) const {
	const bool weighted = settings.XWR_type == 2, against_F2 = settings.refine_against == 2;
	if (all) return weighted ? (against_F2 ? cryst.weighted_GooF2_all : cryst.weighted_GooF1_all) : (against_F2 ? cryst.GooF2_all : cryst.GooF1_all);
	return weighted ? (against_F2 ? cryst.weighted_GooF2 : cryst.weighted_GooF1) : (against_F2 ? cryst.GooF2 : cryst.GooF1);
}

void XCW::ensure_hkl_ordered() {
	if (!hkl_ordered_.empty() || hkl.empty()) {
		return;
	}
	hkl_ordered_.reserve(hkl.size());
	for (const i3& h : hkl) {
		hkl_ordered_.push_back(h);
	}
}

//Per-reflection 1/|H|^2 weights for the residual self-energy criterion,
//U_res ~ Sum_h |dF_h|^2/|H_h|^2, with |H| = 1/d = 2*sin(theta)/lambda.
//(0,0,0) is already excluded from hkl at read time; depends only on geometry.
void XCW::ensure_inv_H2_weights() {
	if (settings.XWR_type == 1 || !inv_H2_.empty()) {
		return;
	}
	ensure_hkl_ordered();
	inv_H2_.resize(cryst.nr_small);
#pragma omp parallel for
	for (int r = 0; r < cryst.nr_small; r++) {
		const double stl = unit_cell.get_stl_of_hkl(hkl_ordered_[r]);
		const double H2 = 4.0 * stl * stl;
		inv_H2_[r] = (H2 > 0.0) ? 1.0 / H2 : 0.0;
	}
}

//z_h = (|F_obs,h| - |F_calc,h|) / sigma_h for the converged F_calc/F_scale over
//strong reflections (|F_obs|/sigma >= opt->xcw_strong_cutoff), tested against
//N(0,1) after a global shape/scale rescale. Full reflection set only; the
//free/working-set cross-validation variant is not implemented.
//Background: tests/P1_test/XCW_plan.md, Src/core/xcw_halting.h.
void XCW::evaluate_gaussian_halting(const double lambda) {
	ensure_hkl_ordered();

	const int n_total = cryst.nr_small;
	vec z_raw;
	vec resolution;
	vec abs_F;
	z_raw.reserve(n_total);
	resolution.reserve(n_total);
	abs_F.reserve(n_total);

	for (int i = 0; i < n_total; i++) {
		if (obs[i].sigma_obs <= 0.0) {
			continue;
		}
		const double f_over_sigma = obs[i].abs_F_obs / obs[i].sigma_obs;
		if (f_over_sigma < opt->xcw_strong_cutoff) {
			continue;
		}
		const double scaled_F_calc = cryst.F_scale * std::abs(F_calc[0][i]);
		const double diff = scaled_F_calc - obs[i].abs_F_obs;
		z_raw.push_back(diff / obs[i].sigma_obs);
		resolution.push_back(unit_cell.get_stl_of_hkl(hkl_ordered_[i]));
		abs_F.push_back(obs[i].abs_F_obs);
	}

	GaussianHaltEntry entry;
	entry.lambda = lambda;
	entry.n_total = n_total;
	entry.n_used = static_cast<int>(z_raw.size());

	if (entry.n_used < 8) {
		XCW_log << "Gaussian halting criterion: only " << entry.n_used
			<< " strong reflections (|F|/sigma >= " << opt->xcw_strong_cutoff
			<< ") at lambda=" << lambda << ", skipping (need >= 8)." << std::endl;
		gaussian_halt_history_.push_back(entry);
		return;
	}

	//Decouple shape from scale: one global factor rescaling z so <z^2> ~ 1, not the
	//resolution-uniform weighting of XCW_plan.md 4.1
	double mean_z2 = 0.0;
	for (const double v : z_raw) {
		mean_z2 += v * v;
	}
	mean_z2 /= entry.n_used;
	entry.sigma_scale = (mean_z2 > 0.0) ? std::sqrt(mean_z2) : 1.0;

	vec z(entry.n_used);
	for (int i = 0; i < entry.n_used; i++) {
		z[i] = z_raw[i] / entry.sigma_scale;
	}

	entry.A2 = anderson_darling_statistic(z);
	entry.ad_reject_5pct = entry.A2 > ANDERSON_DARLING_CRITICAL_5PCT;

	const ProbabilityPlotFit pp = normal_probability_plot_fit(z);
	entry.pp_slope = pp.slope;
	entry.pp_intercept = pp.intercept;

	entry.skewness = sample_skewness(z);
	entry.excess_kurtosis = sample_excess_kurtosis(z);
	entry.jarque_bera = jarque_bera_statistic(z);

	const int n_bins = std::clamp(entry.n_used / 20, 2, 10);
	const BinnedTrend res_trend = binned_z_squared_trend(z, resolution, n_bins);
	entry.resolution_trend_slope = res_trend.slope;
	entry.resolution_trend_r = res_trend.spearman_r;
	entry.resolution_trend_flagged = res_trend.flagged;

	const BinnedTrend int_trend = binned_z_squared_trend(z, abs_F, n_bins);
	entry.intensity_trend_slope = int_trend.slope;
	entry.intensity_trend_r = int_trend.spearman_r;
	entry.intensity_trend_flagged = int_trend.flagged;

	XCW_log << "Gaussian halting criterion at lambda=" << std::fixed << std::setprecision(5) << lambda << ":\n"
		<< "  n_used=" << entry.n_used << "/" << entry.n_total << " (|F|/sigma >= " << opt->xcw_strong_cutoff
		<< "), sigma_scale=" << entry.sigma_scale << "\n"
		<< "  A^2=" << entry.A2 << (entry.ad_reject_5pct ? " (rejects N(0,1) at 5%)" : " (consistent with N(0,1) at 5%)") << "\n"
		<< "  probability-plot slope=" << entry.pp_slope << " intercept=" << entry.pp_intercept << "\n"
		<< "  skewness=" << entry.skewness << " excess_kurtosis=" << entry.excess_kurtosis
		<< " Jarque-Bera=" << entry.jarque_bera << "\n"
		<< "  resolution-binned <z^2> trend: slope=" << entry.resolution_trend_slope
		<< " spearman_r=" << entry.resolution_trend_r << (entry.resolution_trend_flagged ? " [FLAGGED]" : "") << "\n"
		<< "  |F|-binned <z^2> trend: slope=" << entry.intensity_trend_slope
		<< " spearman_r=" << entry.intensity_trend_r << (entry.intensity_trend_flagged ? " [FLAGGED]" : "") << std::endl;

	gaussian_halt_history_.push_back(entry);
}

//Full per-lambda table to XCW_log, then the final recommendation
void XCW::report_gaussian_halting_summary() {
	if (gaussian_halt_history_.empty()) {
		return;
	}

	XCW_log << "\n____________________________________________________________________________\n"
		<< "Gaussian halting criterion summary (tests/P1_test/XCW_plan.md)\n"
		<< " Lambda\t\tA^2\treject5%\tpp_slope\tpp_intercept\tskew\tkurt\tres_trend_r\tint_trend_r\tn_used\n";
	for (const GaussianHaltEntry& e : gaussian_halt_history_) {
		XCW_log << "\t" << std::fixed << std::setprecision(5) << e.lambda
			<< "\t" << std::setprecision(4) << e.A2
			<< "\t" << (e.ad_reject_5pct ? "yes" : "no")
			<< "\t\t" << e.pp_slope << "\t\t" << e.pp_intercept
			<< "\t" << e.skewness << "\t" << e.excess_kurtosis
			<< "\t" << e.resolution_trend_r << "\t\t" << e.intensity_trend_r
			<< "\t" << e.n_used << "\n";
	}

	report_halting_progress_estimate(true);
}

//lambda* is the argmin of A^2 over steps with enough strong reflections. A flagged
//binned-trend test means spatially correlated residuals that a marginal normality
//test would miss, so it is warned about. See XCW.h for the full behaviour.
void XCW::report_halting_progress_estimate(bool is_final) {
	const GaussianHaltEntry* best = nullptr;
	double max_valid_lambda = 0.0;
	vec fit_lambda, fit_A2;
	for (const GaussianHaltEntry& e : gaussian_halt_history_) {
		if (e.n_used < 8) {
			continue;
		}
		if (!best || e.A2 < best->A2) {
			best = &e;
		}
		max_valid_lambda = std::max(max_valid_lambda, e.lambda);
		fit_lambda.push_back(e.lambda);
		fit_A2.push_back(e.A2);
	}
	if (!best) {
		return;
	}

	//A minimum at the last evaluated lambda is a scan-boundary artifact, not a
	//reached lambda*, so it is flagged rather than reported as the answer
	const bool at_boundary = (best->lambda >= max_valid_lambda - 1e-12) && (fit_lambda.size() > 1);
	const bool trend_ok = !best->resolution_trend_flagged && !best->intensity_trend_flagged;

	//AIC picks between candidate forms of the A^2(lambda) trend; quartic only enters
	//once the scan has enough steps, via fit_polynomial's degree + 3 minimum
	std::vector<PolynomialFit> candidates;
	const PolynomialFit fit = choose_best_polynomial_fit(fit_lambda, fit_A2, { 2, 4 }, &candidates);

	std::ostream* streams[2] = { &XCW_log, &std::cout };
	for (std::ostream* s : streams) {
		*s << "____________________________________________________________________________\n";
		if (!is_final) {
			*s << "Gaussian halting criterion: progress update after " << fit_lambda.size() << " lambda steps\n";
		}
		*s << "Recommended halting lambda* = " << std::fixed << std::setprecision(5) << best->lambda
			<< " (A^2=" << std::setprecision(4) << best->A2 << ")";
		if (!trend_ok) {
			*s << " -- WARNING: binned <z^2> trend test flagged at this lambda; "
				<< "residuals may be spatially correlated, inspect before trusting lambda*.";
		}
		*s << "\n";

		for (const PolynomialFit& c : candidates) {
			*s << "  candidate fit: degree=" << c.degree;
			if (!c.valid) {
				*s << " -- not enough points yet (need >= degree+3 evaluated lambda steps)\n";
				continue;
			}
			*s << " RSS=" << std::setprecision(4) << c.rss << " R^2=" << c.r_squared
				<< " AIC=" << c.aic << ((fit.valid && c.degree == fit.degree) ? " [chosen]" : "")
				<< "\n";
		}

		if (at_boundary) {
			*s << "WARNING: A^2 is still falling at the last evaluated lambda (" << std::setprecision(5)
				<< max_valid_lambda << "); this is a scan-boundary value, not a confirmed interior minimum. "
				<< "Extend the scan (larger max_value in -do_XCW stepsize max_value) to find the true optimum.";
			if (fit.valid && fit.has_minimum && fit.vertex_x > max_valid_lambda) {
				*s << " Degree-" << fit.degree << " polynomial extrapolation (best fit by AIC, R^2="
					<< std::setprecision(4) << fit.r_squared << ") of the A^2(lambda) trend so far estimates "
					<< "the minimum near lambda ~= " << std::setprecision(5) << fit.vertex_x
					<< " (extrapolated, not yet observed -- treat as a rough guide to where to extend the scan, not a final answer).";
			}
			else if (!fit.valid) {
				*s << " Not enough points yet to extrapolate an estimate (need >= 5 evaluated lambda steps).";
			}
			else if (!fit.has_minimum) {
				*s << " The trend so far is not curving upward yet within the search window; "
					<< "cannot extrapolate a stopping estimate, extend the scan further.";
			}
			else {
				//fit.vertex_x <= max_valid_lambda: the fitted curve turns before the last step
				//while the observed A^2 is still falling there - a flat, noisy trend, so no
				//estimate is offered
				*s << " The degree-" << fit.degree << " fit places its minimum at lambda ~= "
					<< std::setprecision(5) << fit.vertex_x << ", inside the range already scanned, "
					<< "which contradicts A^2 still falling at the last step; the trend is too flat "
					<< "here to extrapolate from, only further steps can decide it.";
			}
			*s << "\n";
		}
	}
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
	const i3& hkl_temp = *it;
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


size_t XCW::tri_index(int mu, int nu) const noexcept {
	return mu * cryst.nmo - (mu * (mu - 1)) / 2 + (nu - mu);
}

void XCW::eval_I_anom_disp(std::vector<ao_data>& ao_data_shells, bool read) {
	cvec2 DW_fact, phase_fact, translation_phase;
	eval_phase(phase_fact);
	eval_DW(DW_fact);
	eval_translation_phase(translation_phase);
	size_t kept_on_disk = 0;
	bool single_on_disk = false;
	if (settings.read_tensor && !settings.i_tensor_file_path.empty()
		&& i_tensor_file::matches(i_tensor_path(), cryst.nr_small, cryst.nmo, kept_on_disk, single_on_disk)) {
		//A streamed tensor already there and big enough for this problem. It depends on the
		//geometry, the basis and the reflections and on none of the refinement settings, so
		//a second run that changes those can read it rather than spend the build again.
		//open() checks the header and throws if the shape does not match, which is what
		//stops a tensor from a different structure being used by accident.
		i_compact_ = kept_on_disk;
		const size_t packed = i_compact_;
		const char* source = "";
		bool automatic = false;
		i_streamed_ = items_within_budget(static_cast<size_t>(cryst.nr_small),
			i_tensor_file::block_bytes(i_compact_, single_on_disk), i_budget(source, automatic)) != 0;
		i_window_ = std::max(1, std::min(cryst.nr_small, 64));
		open_i_stream_for_reading();
		i_pair_mu_ = i_file_.pair_mu();
		i_pair_nu_ = i_file_.pair_nu();
		//The file's element type is kept as it is: a single-precision tensor cannot regain
		//anything by widening, and a double one is narrowed only on request
		const char* f = std::getenv("NOSPHERA2_XCW_I_FLOAT"); // Flawfinder: ignore
		i_float_ = single_on_disk || (!i_streamed_ && (settings.i_tensor_single || (f && std::atoi(f) != 0)));
		std::cout << "I tensor read from " << i_tensor_path().string()
			<< " (" << (i_tensor_file::total_bytes(cryst.nr_small, i_compact_, single_on_disk) / 1048576.0)
			<< " MB" << (single_on_disk ? ", single precision" : "") << "), not recomputed"
			<< (i_streamed_ ? ", read a window at a time" : ", held in memory") << std::endl;
		if (i_float_ && !single_on_disk)
			std::cout << "NOTE: the tensor on disk is double precision; it is narrowed to single as i_float asks" << std::endl;
		if (single_on_disk && settings.i_tensor_double)
			std::cout << "NOTE: the tensor on disk is single precision; i_double cannot widen it, it is used as stored" << std::endl;
		if (!i_streamed_) {
			if (i_float_)
				I32.assign(static_cast<size_t>(cryst.nr_small) * packed, std::complex<float>{});
			else
				I.resize(static_cast<size_t>(cryst.nr_small) * packed);
			for (int r0 = 0; r0 < cryst.nr_small; r0 += i_window_) {
				const int r1 = std::min(cryst.nr_small, r0 + i_window_);
				i_file_.load(r0, r1);
				for (int r = r0; r < r1; r++) {
					if (single_on_disk)
						std::copy(i_file_.block32(r), i_file_.block32(r) + packed, I32.data() + static_cast<size_t>(r) * packed);
					else if (i_float_)
						for (size_t i = 0; i < packed; i++)
							I32[static_cast<size_t>(r) * packed + i] = std::complex<float>(static_cast<float>(i_file_.block(r)[i].real()), static_cast<float>(i_file_.block(r)[i].imag()));
					else
						std::copy(i_file_.block(r), i_file_.block(r) + packed, I.data() + static_cast<size_t>(r) * packed);
				}
			}
			i_file_.close();
		}
	}
	else {
		double time_taken;
		long long screen_counter = 0;
		long long skipped_grids = 0;
		eval_I(ao_data_shells, DW_fact, phase_fact, translation_phase, time_taken, screen_counter, skipped_grids);
		if (!(opt->no_date)) {
			std::cout << std::fixed << std::setprecision(2) << "Time taken for XCW integrals: " << time_taken << " seconds. \n";
		}
		start_i_save();
		std::cout << std::fixed << std::setprecision(2) << "Screened out " << screen_counter << " unique pairs of mu, nu (" << static_cast<size_t>(screen_counter) / (static_cast<double>(cryst.nmo * (cryst.nmo + 1)) / 2) * 100.00 << "%) \n";
		std::cout << std::fixed << std::setprecision(2) << "Skipped evaluation of " << skipped_grids << " grids (" << static_cast<double>(skipped_grids) / ((static_cast<double>(cryst.nmo * (cryst.nmo + 1)) / 2) * cryst.nr * cryst.ncen) * 100.00 << "%) \n";

	}
	eval_anom_disp(DW_fact, phase_fact, translation_phase);
	// closing function
}

//Whether the I tensor is held or streamed, and the largest window that fits the budget
//The tensor is written on a thread while the refinement runs. Nothing in the SCF modifies
//it - both walks only read - so a reader alongside them needs no lock, and the run pays only
//the disk bandwidth, which it is not competing for while it works out of memory.
void XCW::start_i_save()
{
	if (settings.i_tensor_save_path.empty() || i_streamed_) return;
	const size_t packed = i_compact_;
	const int nr = cryst.nr_small;
	const std::filesystem::path path = settings.i_tensor_save_path;
	const bool from_float = i_float_;
	std::cout << "Writing the I tensor to " << path.string()
		<< " in the background; a later run can `read " << path.string()
		<< "` instead of building it" << std::endl;
	i_writer_ = std::thread([this, path, nr, packed, from_float]() {
		try {
			i_tensor_file out;
			out.create(path, nr, cryst.nmo, i_pair_mu_, i_pair_nu_, from_float);
			for (int r = 0; r < nr; r++) {
				if (from_float) out.write_block(r, I32.data() + static_cast<size_t>(r) * packed);
				else out.write_block(r, I.data() + static_cast<size_t>(r) * packed);
			}
			out.finish_write();
		}
		catch (const std::exception& e) { i_writer_error_ = e.what(); }
		catch (...) { i_writer_error_ = "unknown error"; }
		});
}

//Called before the run ends, and before anything that could invalidate the tensor. A thread
//left running past main is the bug 939268f was about; this one also holds a file handle.
void XCW::finish_i_save()
{
	if (!i_writer_.joinable()) return;
	i_writer_.join();
	if (!i_writer_error_.empty())
		std::cout << "Could not write the I tensor to "
			<< settings.i_tensor_save_path.string() << ": " << i_writer_error_
			<< " (the refinement itself is unaffected)" << std::endl;
	else
		std::cout << "I tensor written to " << settings.i_tensor_save_path.string() << std::endl;
}

size_t XCW::i_budget(const char*& source, bool& automatic) const {
	size_t budget = settings.i_tensor_max_mb * 1024ULL * 1024ULL;
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

void XCW::decide_i_storage() {
	const size_t per_block = i_tensor_file::block_bytes(i_compact_, i_float_);
	const size_t total = i_tensor_file::total_bytes(cryst.nr_small, i_compact_, i_float_);

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
	const size_t w = items_within_budget(static_cast<size_t>(cryst.nr_small), per_block, budget);
	i_streamed_ = (w != 0);
	if (!i_streamed_) {
		//Announced only when someone asked about memory: saying it unconditionally
		//shifts every reference output by a line, and the automatic budget would say it on
		//every run - including the reference tests, which is why it stays quiet there.
		if ((budget > 0 && !automatic) || ProgressBar::report_counts) {
			std::cout << std::fixed << std::setprecision(2)
				<< "I tensor held in memory: " << (total / 1048576.0) << " MB"
				<< (i_float_ ? " (single precision)" : "");
			if (budget > 0)
				std::cout << " (fits the " << (budget / 1048576.0) << " MB " << source << " budget)";
			std::cout << std::endl;
		}
		return;
	}
	i_window_ = static_cast<int>(std::min(w, static_cast<size_t>(cryst.nr_small)));
	i_file_.create(i_tensor_path(), cryst.nr_small, cryst.nmo, i_pair_mu_, i_pair_nu_, i_float_);
	std::cout << std::fixed << std::setprecision(2)
		<< "I tensor streamed to disk: " << (total / 1048576.0) << " MB" << (i_float_ ? " (single precision)" : "") << " total, "
		<< i_window_ << " of " << cryst.nr_small << " reflections resident ("
		<< (i_window_ * per_block / 1048576.0) << " MB) to fit " << source
		<< " (" << (budget / 1048576.0) << " MB)" << std::endl;
	if (i_window_ == 1 && per_block > budget)
		std::cout << "  NOTE: one reflection alone is " << (per_block / 1048576.0)
		<< " MB, over the budget. Running one at a time." << std::endl;
}

std::filesystem::path XCW::i_tensor_path() const {
	return settings.i_tensor_file_path.empty()
		? std::filesystem::path(i_tensor_default)
		: settings.i_tensor_file_path;
}

void XCW::open_i_stream_for_reading() {
	i_file_.open(i_tensor_path(), static_cast<size_t>(i_window_));
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

//C = op(A) op(B) for column-major occ matrices through MKL, Src/core's Eigen being serial
static occ::Mat gemm(const occ::Mat& A, const occ::Mat& B, const bool ta = false, const bool tb = false)
{
	const int m = static_cast<int>(ta ? A.cols() : A.rows()), k = static_cast<int>(ta ? A.rows() : A.cols()), n = static_cast<int>(tb ? B.rows() : B.cols());
	occ::Mat C(m, n);
	cblas_dgemm(CblasColMajor, ta ? CblasTrans : CblasNoTrans, tb ? CblasTrans : CblasNoTrans, m, n, k, 1.0, A.data(), static_cast<int>(A.rows()), B.data(), static_cast<int>(B.rows()), 0.0, C.data(), m);
	return C;
}

void XCW::eval_I(std::vector<ao_data>& ao_data_shells, cvec2& DW_fact, cvec2& phase_fact, cvec2& translation_phase, double& time_taken, long long& screen_counter, long long& skipped_grids_) {
	long long skipped_grids = 0;
	const int packed_size = (cryst.nmo * (cryst.nmo + 1)) / 2;
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
	dummy_wave.delete_unoccupied_MOs();
	bvec needs_grid(cryst.ncen, true);
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

	ivec2 asym_lookup(cryst.nr_small);
	for (r = 0; r < cryst.nr_small; r++) {
		asym_lookup[r] = generate_asym_lookup(r);
	}
	const unsigned int num_syms = asym_lookup[0].size();

	if (std::getenv("NOSPHERA2_DEBUG_LOOKUP")) { // Flawfinder: ignore
		long long misses = 0, total = 0;
		for (r = 0; r < cryst.nr_small; r++) {
			for (int s = 0; s < static_cast<int>(num_syms); s++) {
				total++;
				if (asym_lookup[r][s] == 0) {
					// index 0 is ambiguous: a genuine hit on hkl_enlarged's first entry, or the
					// silent "not found" fallback in generate_asym_lookup. Recompute by hand to
					// tell them apart.
					auto it = hkl.begin();
					std::advance(it, r);
					ivec3 rots = unit_cell.get_sym();
					i3 tempv{ 0,0,0 };
					const i3& hkl_temp = *it;
					for (int h = 0; h < 3; h++)
						for (int j = 0; j < 3; j++)
							tempv[j] += hkl_temp[h] * rots[j][h][s];
					if (hkl_enlarged.find(tempv) == hkl_enlarged.end()) misses++;
				}
			}
		}
		std::cerr << "asym_lookup misses: " << misses << " of " << total
			<< "  (hkl_enlarged size " << hkl_enlarged.size() << ")" << std::endl;
		std::cerr.flush();
		std::exit(0);
	}

	vec2 grid_positions(cryst.ncen);
	for (int at = 0; at < cryst.ncen; at++) {
		grid_positions[at] = { dummy_wave.get_atom_pos(at)[0], dummy_wave.get_atom_pos(at)[1], dummy_wave.get_atom_pos(at)[2] };
	}

	// Precompute screening
	ivec2 skip(cryst.nmo, ivec(cryst.nmo, 0));
	{
		double e_tol = 0.0005;
		const double root_inv_four_pi = std::sqrt(constants::INV_FOUR_PI);
		for (mu = 0; mu < cryst.nmo; mu++) {
			const ao_data& mu_prims = ao_data_shells[mu];
			const std::vector<primitive>& mu_primitives = mu_prims.prims;
			const double& mp0 = mu_prims.pos[0];
			const double& mp1 = mu_prims.pos[1];
			const double& mp2 = mu_prims.pos[2];
			for (nu = mu + 1; nu < cryst.nmo; nu++) {
				const ao_data& nu_prims = ao_data_shells[nu];
				const std::vector<primitive>& nu_primitives = nu_prims.prims;
				const double& np0 = nu_prims.pos[0];
				const double& np1 = nu_prims.pos[1];
				const double& np2 = nu_prims.pos[2];
				const double dist0 = mp0 - np0;
				const double dist1 = mp1 - np1;
				const double dist2 = mp2 - np2;
				const double dist = dist0 * dist0 + dist1 * dist1 + dist2 * dist2;
				if (dist < 1e-5) {
					continue;
				}

				double c = 0;
				double mu_min = std::numeric_limits<double>::max();
				double nu_min = std::numeric_limits<double>::max();
				const int mu_l = mu_prims.prims[0].get_type();
				const double mu_l_half = 0.5 * mu_l;
				const double temp_mu = std::sqrt((2 * mu_l + 1) * constants::INV_FOUR_PI) * std::exp(-mu_l_half);
				const int nu_l = nu_prims.prims[0].get_type();
				const double nu_l_half = 0.5 * nu_l;
				const double temp_nu = std::sqrt((2 * nu_l + 1) * constants::INV_FOUR_PI) * std::exp(-nu_l_half);
				std::vector<std::pair<double, double>> pairs;
				pairs.reserve(mu_primitives.size() * nu_primitives.size());
				for (int k = 0; k < mu_primitives.size(); k++) {
					mu_min = std::min(mu_min, mu_primitives[k].get_exp());
					const double c_k = std::abs(mu_primitives[k].get_coef());
					const double alpha_k = mu_primitives[k].get_exp();
					const double N_k = mu_l == 0 ? root_inv_four_pi : temp_mu * std::pow(mu_l / alpha_k, mu_l_half);
					for (int l = 0; l < nu_primitives.size(); l++) {
						nu_min = std::min(nu_min, nu_primitives[l].get_exp());
						const double c_l = std::abs(nu_primitives[l].get_coef());
						const double alpha_l = nu_primitives[l].get_exp();
						const double N_l = nu_l == 0 ? root_inv_four_pi : temp_nu * std::pow(nu_l / alpha_l, nu_l_half);
						const double N_kl = N_k * N_l * std::pow(constants::TWO_PI / (alpha_k + alpha_l), 1.5);
						const double temp1 = c_k * c_l * N_kl;
						c += temp1;
						pairs.emplace_back(temp1, alpha_k * alpha_l / (2.0 * (alpha_k + alpha_l)));
					}
				}
				const double gamma = 2 * (mu_min + nu_min) / (mu_min * nu_min);
				const double cutoff = std::log(c / e_tol) * gamma;
				//Newton method for finding correct cutoff
				double newton_cutoff;
				if (cutoff <= 0.0) {
					newton_cutoff = 0.0;
				}
				else {
					double lo = 0.0, hi = cutoff;
					newton_cutoff = 0.5 * (lo + hi);
					for (int iter = 0; iter < 50; iter++) {
						double upper_bound = 0.0, bound_derivative = 0.0;
						for (const auto& [weight, gamma_kl] : pairs) {
							const double upper_bound_temp = weight * std::exp(-gamma_kl * newton_cutoff);
							upper_bound += upper_bound_temp;
							bound_derivative -= gamma_kl * upper_bound_temp;
						}
						const double delta = upper_bound - e_tol;
						if (delta >= 0.0) lo = newton_cutoff; else hi = newton_cutoff;
						double next = newton_cutoff - delta / bound_derivative;
						if (!(next > lo) || !(next < hi)) next = 0.5 * (lo + hi);
						const double step = std::abs(next - newton_cutoff);
						newton_cutoff = next;
						if (step < 1e-12 * hi) break;
					}
				}
				if (dist > newton_cutoff) {
					skip[mu][nu] = 1;
				}
			}
		}
	}

	// Grid screening
	constexpr double maximum_ao_grid_cutoff = 12;
	constexpr double minimum_ao_grid_cutoff = 11;
	double minimum_primitive_exponent = std::numeric_limits<double>::max();
	vec ao_grid_cutoff_squared(cryst.nmo);
	{
		vec ao_minimum_exponent(cryst.nmo);
		for (int ao = 0; ao < cryst.nmo; ao++) {
			double min_exp = std::numeric_limits<double>::max();
			for (const primitive& prim : ao_data_shells[ao].prims) {
				min_exp = std::min(min_exp, prim.get_exp());
			}
			ao_minimum_exponent[ao] = min_exp;
			minimum_primitive_exponent = std::min(minimum_primitive_exponent, min_exp);
		}
		for (int ao = 0; ao < cryst.nmo; ao++) {
			const double adaptive_cutoff = maximum_ao_grid_cutoff * std::sqrt(minimum_primitive_exponent / ao_minimum_exponent[ao]);
			const double cutoff = std::clamp(adaptive_cutoff, minimum_ao_grid_cutoff, maximum_ao_grid_cutoff);
			ao_grid_cutoff_squared[ao] = cutoff * cutoff;
		}
	}

	const int n_atom_grids = std::min(n_grids, cryst.ncen);
	ivec2 ao_prefix_end(cryst.nmo, ivec(n_atom_grids));
	bvec2 ao_within_cutoff(cryst.nmo, bvec(n_atom_grids));

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
	for (int ao = 0; ao < cryst.nmo; ao++) {
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
	vec3 mu_vals(cryst.nmo, vec2(n_grids));
#pragma omp parallel for schedule(dynamic)
	for (mu = 0; mu < cryst.nmo; mu++) {
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
				d_mu[3] = std::hypot(d_mu[0], d_mu[1], d_mu[2]);
				if (d_mu[3] * d_mu[3] > ao_grid_cutoff_squared[mu]) {
					local_mu_vals_ptr[p] = 0.0;
				}
				else {
					local_mu_vals_ptr[p] = dummy_wave.eval_ao(d_mu, mu_primitives, mu_prims.m);
				}
			}
		}
	}
	std::cout << "AO values calculated for all grids." << std::endl;
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
			for (int mu = 0; mu < cryst.nmo; mu++) {
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
			for (int ao = 0; ao < cryst.nmo; ao++) kept += ao_within_cutoff[ao][g] ? 1 : 0;
			for (int ao = 0; ao < cryst.nmo; ao++) {
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
				g, npts, cryst.nmo, static_cast<int>(kept),
				active_full, active_batch);
		}
	}

	ivec2 active_grids(packed_size);
	ivec skipped_grids_per_pair(packed_size, 0);
	for (mu = 0; mu < cryst.nmo; mu++) {
		const std::vector<bool>& mu_within = ao_within_cutoff[mu];
		for (nu = mu; nu < cryst.nmo; nu++) {
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
	i_pair_mu_.clear();
	i_pair_nu_.clear();
	for (mu = 0; mu < cryst.nmo; mu++)
		for (nu = mu; nu < cryst.nmo; nu++)
			if (!skip[mu][nu]) {
				tri_compact[tri_index(mu, nu)] = static_cast<int>(i_pair_mu_.size());
				i_pair_mu_.push_back(mu);
				i_pair_nu_.push_back(nu);
			}
	i_compact_ = i_pair_mu_.size();
	ivec2 grid_active_aos(n_atom_grids);
	vec2 grid_ao_values(n_atom_grids);
#pragma omp parallel for schedule(dynamic)
	for (int g = 0; g < n_atom_grids; g++) {
		ivec& active_aos = grid_active_aos[g];
		vec& values = grid_ao_values[g];
		for (int mu = 0; mu < cryst.nmo; mu++) {
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
				work_unscreened += static_cast<double>(cryst.nmo) * cryst.nmo * b.point_count;
			}
		//Per reflection and symmetry operation the sum is a small number and says nothing;
		//what the run costs is that times both, so scale it before printing or the figure
		//reads as a thousandth of the truth.
		const double runs = static_cast<double>(cryst.nr_small) * static_cast<double>(num_syms);
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
	for (int mu = 0; mu < cryst.nmo; mu++) {
		for (int nu = mu; nu < cryst.nmo; nu++) {
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
		ivec compact_flat(static_cast<size_t>(cryst.nmo) * cryst.nmo, -1);
		for (int m = 0; m < cryst.nmo; m++)
			for (int n = m; n < cryst.nmo; n++)
				compact_flat[static_cast<size_t>(m) * cryst.nmo + n] = tri_compact[tri_index(m, n)];
		vec fd1, fd2, fd3, fw;
		for (int gg = 0; gg < n_atom_grids; gg++) {
			fd1.insert(fd1.end(), d1[gg].begin(), d1[gg].begin() + points[gg]);
			fd2.insert(fd2.end(), d2[gg].begin(), d2[gg].begin() + points[gg]);
			fd3.insert(fd3.end(), d3[gg].begin(), d3[gg].begin() + points[gg]);
			fw.insert(fw.end(), weights[gg].begin(), weights[gg].begin() + points[gg]);
		}
		itensor_gpu_layout L;
		L.nmo = cryst.nmo; L.packed = static_cast<int>(i_compact_); L.n_grids = n_atom_grids;
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
				* static_cast<double>(cryst.nr_small) * static_cast<double>(num_syms);
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
		i_float_ = settings.i_tensor_single || (f && std::atoi(f) != 0) || (single_build && !settings.i_tensor_double);
		if (single_build && !i_float_)
			std::cout << "NOTE: the I tensor is built in single precision and held in double as i_double asks" << std::endl;
		if (!single_build && i_float_)
			std::cout << "NOTE: the I tensor is built in double precision and narrowed to single as i_float asks" << std::endl;
	}
	decide_i_storage();
	if (!i_streamed_) {
		if (i_float_)
			I32.assign(static_cast<size_t>(cryst.nr_small) * i_compact_, std::complex<float>{});
		else
			I.assign(static_cast<size_t>(cryst.nr_small) * i_compact_, cdouble{});
	}
	//After the storage line: the bar owns the console from here until the last reflection
	if (!(opt->no_date)) {
		pb.emplace((unsigned long long)cryst.nr_small, 60, "=", "|", "Calculating XCW integrals...", std::cout);
	}
#if defined(NOSPHERA2_USE_GPU) || defined(NOSPHERA2_USE_METAL)
	if (itensor_on_gpu) gpu_thread = std::thread([&]() {
		const auto gpu_start = std::chrono::high_resolution_clock::now();
		//The device takes reflections in batches; the CPU threads keep taking them one at
		//a time from the same counter
		const int gpu_batch = itensor_gpu_batch(static_cast<int>(num_syms));
		cvec blk_gpu;
		if (i_streamed_ || i_float_) blk_gpu.assign(static_cast<size_t>(gpu_batch) * i_compact_, cdouble{});
		vec kxs(static_cast<size_t>(gpu_batch) * num_syms), kys(kxs.size()), kzs(kxs.size());
		cvec facs(kxs.size() * n_atom_grids);
		int done = 0;
		auto collect_gpu = [&](const int rr0, const int n, const int slot) {
			if (i_streamed_ || i_float_) std::fill(blk_gpu.begin(), blk_gpu.end(), cdouble{});
			cdouble* const I_rr = (i_streamed_ || i_float_) ? blk_gpu.data()
									  : I.data() + static_cast<size_t>(rr0) * i_compact_;
			if (!itensor_gpu_collect(slot, I_rr, static_cast<long long>(i_compact_)))
				err_checkf(false, "I tensor GPU read-back failed", std::cout);
			for (int r = 0; r < n && (i_streamed_ || i_float_); r++) {
				const cdouble* const row = blk_gpu.data() + static_cast<size_t>(r) * i_compact_;
				if (i_streamed_) {
					std::lock_guard<std::mutex> lock(i_write_mutex);
					i_file_.write_block(rr0 + r, row);
				}
				else if (i_float_) {
					std::complex<float>* const dst = I32.data() + static_cast<size_t>(rr0 + r) * i_compact_;
					for (size_t i = 0; i < i_compact_; i++)
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
			if (rr0 >= cryst.nr_small) break;
			const int n = std::min(gpu_batch, cryst.nr_small - rr0);
			for (int r = 0; r < n; r++)
				for (int sy = 0; sy < static_cast<int>(num_syms); sy++) {
					const int rr = rr0 + r, c = r * static_cast<int>(num_syms) + sy;
					kxs[c] = k_pt[0][asym_lookup[rr][sy]];
					kys[c] = k_pt[1][asym_lookup[rr][sy]];
					kzs[c] = k_pt[2][asym_lookup[rr][sy]];
					for (int gg = 0; gg < n_atom_grids; gg++)
						facs[static_cast<size_t>(c) * n_atom_grids + gg] =
						asym_atoms[gg].asym_fact * DW_fact[gg][asym_lookup[rr][sy]]
						* phase_fact[gg][asym_lookup[rr][sy]] * translation_phase[rr][sy];
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
				done, cryst.nr_small);
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
	itensor_flops *= static_cast<double>(cryst.nr_small) * static_cast<double>(num_syms);
	itensor_unscreened_tile_flops *= static_cast<double>(cryst.nr_small) * static_cast<double>(num_syms);
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
			if (i_streamed_ || i_float_) blk.assign(i_compact_, cdouble{});

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
				if (r >= cryst.nr_small) break;
				if (itensor_on_gpu && my_done > 0 && gpu_ns_per_refl.load() > 0 &&
					static_cast<long long>(cryst.nr_small - r) * gpu_ns_per_refl.load() < my_ns / my_done) break;
				if (!next_refl.compare_exchange_strong(r, r + 1)) continue;
				const auto r_start = std::chrono::high_resolution_clock::now();
				if (i_streamed_ || i_float_) std::fill(blk.begin(), blk.end(), cdouble{});
				cdouble* const I_r = (i_streamed_ || i_float_) ? blk.data()
					: I.data() + static_cast<size_t>(r) * i_compact_;
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
						const cdouble* DW_fact_g = DW_fact[g].data();
						const double DW_im = DW_fact_g[idx].imag();
						const cdouble* phase_fact_g = phase_fact[g].data();
						const double phase_im = phase_fact_g[idx].imag();
						const double DW_re = DW_fact_g[idx].real();
						const double phase_re = phase_fact_g[idx].real();
						// Precompute basis function independent factors
						const cdouble grid_factor = cdouble(asym_fact * (DW_re * phase_re - DW_im * phase_im),
							asym_fact * (DW_re * phase_im + DW_im * phase_re));
						const cdouble factor = grid_factor * translation_phase[r][syms];
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
				if (i_streamed_) {
					//A write is packed_size * 16 bytes against a whole reflection's worth
					//of integration, so the lock is not on the hot path
					std::lock_guard<std::mutex> lock(i_write_mutex);
					i_file_.write_block(r, blk.data());
				}
				else if (i_float_) {
					std::complex<float>* const dst = I32.data() + static_cast<size_t>(r) * i_compact_;
					for (size_t i = 0; i < i_compact_; i++)
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
	if (i_streamed_) {
		i_file_.finish_write();
		open_i_stream_for_reading();
	}
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = end - start;
	skipped_grids_ = skipped_grids * cryst.nr_small;
	time_taken = std::chrono::duration<double>(duration).count();
	throughput::record("XCW I tensor", itensor_on_gpu, itensor_flops,
		1.0e3 * std::chrono::duration<double>(duration).count());
	if (!(opt->no_date) && pb) {
		XCW_log << "Time taken for XCW integrals: " << std::fixed << std::setprecision(2) << std::chrono::duration<double>(duration).count() << " seconds." << std::endl;
	}
}

void XCW::calc_F_calc(const dMatrix2& D) {
	// Density matrix from occ is half of what I need, so times 2 and times (2x2)=4
	//Streamed or resident the walk is the same; the outer loop is one window when
	//the tensor is resident
	const int step = std::max(1, i_streamed_ ? i_window_ : cryst.nr_small);
	//The parallel region wraps the window loop: entering one per window would pay
	//team startup and a barrier for a few reflections of work. omp single does the
	//read, and its implicit barrier stops a thread entering a window not yet loaded.
	//load() can throw and an exception must not leave an OpenMP structured block,
	//so it is recorded and rethrown after the region.
	std::string io_error;
#if defined(NOSPHERA2_USE_GPU) || defined(NOSPHERA2_USE_METAL)
	if (i_on_device_) {
		vec w(i_compact_);
		for (size_t k = 0; k < i_compact_; k++)
			w[k] = (i_pair_mu_[k] == i_pair_nu_[k] ? 2.0 : 4.0) * D(i_pair_mu_[k], i_pair_nu_[k]);
		err_checkf(itensor_gpu_rows(w.data(), F_calc[1].data(), F_calc[0].data()), "I tensor walk on the device failed", std::cout);
		return;
	}
#endif
#pragma omp parallel
	{
		for (int r0 = 0; r0 < cryst.nr_small; r0 += step) {
			const int r1 = std::min(r0 + step, cryst.nr_small);
			if (i_streamed_) {
#pragma omp single
				{
					try { i_file_.load(r0, r1); }
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
					cdouble sum = F_calc[1][r];
					const int* pmu = i_pair_mu_.data();
					const int* pnu = i_pair_nu_.data();
					for (size_t k = 0; k < i_compact_; k++)
						sum += (pmu[k] == pnu[k] ? 2.0 : 4.0) * cdouble(I_r[k]) * D(pmu[k], pnu[k]);
					F_calc[0][r] = sum;
				};
				if (i_float_) accumulate(i_block32(r)); else accumulate(i_block(r));
			}
		}
	}
	if (!io_error.empty()) throw std::runtime_error(io_error);
}

void XCW::calc_perturb(occ::Mat& perturb, const occ::qm::SCF<occ::qm::HartreeFock>& scf) {
	ensure_inv_H2_weights();

	//The four (XWR_type, refine_against) combinations differ only in the per-reflection
	//scalar and the prefactor, so one walk over I serves all of them
	const int xwr = static_cast<int>(settings.XWR_type), ref = static_cast<int>(settings.refine_against);
	const bool against_F2 = (ref == 2);
	const bool weighted = (xwr == 2);
	const bool valid = (xwr == 1 || xwr == 2) && (ref == 1 || ref == 2);
	if (!valid) XCW_log << "Invalid refinement option" << std::endl;
	const double scale_sq = cryst.F_scale * cryst.F_scale;
	const double prefactor = against_F2
		? 4.0 * scale_sq / (cryst.n_fit - settings.n_params)
		: 2.0 * cryst.F_scale / (cryst.n_fit - settings.n_params);

	cvec pre(cryst.nr_small);
#pragma omp parallel for
	for (int r = 0; r < cryst.nr_small; r++) {
		if (!valid || !fit_mask_[r]) continue;
		cdouble precompute;
		if (against_F2) {
			const double F_calc_abs_sq = std::pow(std::abs(F_calc[0][r]), 2);
			precompute = std::conj(F_calc[0][r]) * (scale_sq * F_calc_abs_sq - obs[r].F_obs2) / (obs[r].sigma_obs2 * obs[r].sigma_obs2);
		}
		else {
			const double F_calc_abs = std::abs(F_calc[0][r]);
			precompute = std::conj(F_calc[0][r]) * (cryst.F_scale * F_calc_abs - obs[r].abs_F_obs) / (obs[r].sigma_obs * obs[r].sigma_obs * F_calc_abs);
		}
		if (weighted) precompute *= inv_H2_[r];
		pre[r] = precompute;
	}
	contract_I(perturb, pre);
	perturb *= prefactor;
	if (scf.ctx.mo.kind == occ::qm::SpinorbitalKind::Unrestricted) {
		perturb.conservativeResize(2 * cryst.nmo, Eigen::NoChange);
		perturb.bottomRows(cryst.nmo) = perturb.topRows(cryst.nmo);
	}
	//closing function
}

void XCW::contract_I(occ::Mat& out, const cvec& pre) {
	out.setZero(cryst.nmo, cryst.nmo);
	const int step = std::max(1, i_streamed_ ? i_window_ : cryst.nr_small);
	bool on_device = false;
#if defined(NOSPHERA2_USE_GPU) || defined(NOSPHERA2_USE_METAL)
	if (i_on_device_) {
		vec dev(i_compact_);
		err_checkf(itensor_gpu_cols(pre.data(), dev.data()), "I tensor walk on the device failed", std::cout);
		for (size_t k = 0; k < i_compact_; k++) out(i_pair_mu_[k], i_pair_nu_[k]) = dev[k];
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
		occ::Mat local = occ::Mat::Zero(cryst.nmo, cryst.nmo);
		double* local_ptr = local.data();
		for (int r0 = 0; !on_device && r0 < cryst.nr_small; r0 += step) {
			const int r1 = std::min(r0 + step, cryst.nr_small);
			if (i_streamed_) {
#pragma omp single
				{
					try { i_file_.load(r0, r1); }
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
					const int* pmu = i_pair_mu_.data();
					const int* pnu = i_pair_nu_.data();
					for (size_t k = 0; k < i_compact_; k++) {
						const double vr = static_cast<double>(I_r[k].real());
						const double vi = static_cast<double>(I_r[k].imag());
						local_ptr[pnu[k] * cryst.nmo + pmu[k]] += precompute.real() * vr - precompute.imag() * vi;
					}
				};
				if (i_float_) accumulate(i_block32(r)); else accumulate(i_block(r));
			}
		}
		parts[omp_get_thread_num()].swap(local);
	}
	if (!io_error.empty()) throw std::runtime_error(io_error);
	//The partials outweigh the matrix many times over, so their sum is parallel too
	if (!on_device) {
#pragma omp parallel for schedule(static)
		for (int mu = 0; mu < cryst.nmo; mu++) {
			for (int nu = mu; nu < cryst.nmo; nu++) {
				double sum = 0.0;
				for (int t = 0; t < static_cast<int>(parts.size()); t++)
					if (parts[t].size() != 0) sum += parts[t](mu, nu);
				out(mu, nu) = sum;
			}
		}
	}
	for (int mu = 0; mu < cryst.nmo; mu++) {
		for (int nu = mu + 1; nu < cryst.nmo; nu++) {
			out(nu, mu) = out(mu, nu);
		}
	}
}

void XCW::setup_SCF_mol(occ::core::Molecule& mol) {
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
	mol.set_charge(settings.charge);
	mol.set_multiplicity(settings.multiplicity);
}

//The same walk as OCC's SOAD guess with a converged Hartree-Fock density in place of the
//tabulated atoms: F = H + G(D_small) through the mixed-basis build, the orbitals from its
//diagonalisation. The spin blocks are summed, the first iteration polarises them again.
void XCW::small_basis_guess(occ::qm::SCF<occ::qm::HartreeFock>& scf) {
	occ::core::Molecule mol;
	setup_SCF_mol(mol);
	occ::qm::AOBasis small_bs;
	setup_basis(mol, settings.guess_basis_name, small_bs);
	occ::qm::HartreeFock hf_small(small_bs);
	occ::qm::SCF scf_small(hf_small, settings.hf_type);
	scf_small.set_charge_multiplicity(settings.charge, settings.multiplicity);
	scf_small.maxiter = settings.max_scf_iterations;
	const double e_small = scf_small.compute_scf_energy();
	XCW_log << "XCW: initial guess from a " << settings.guess_basis_name << " Hartree-Fock (" << small_bs.nbf()
		<< " functions), E = " << std::fixed << std::setprecision(8) << e_small << " Eh" << std::endl;

	occ::qm::MolecularOrbitals guess;
	guess.kind = scf.ctx.mo.kind;
	guess.n_ao = static_cast<int>(small_bs.nbf());
	guess.n_alpha = scf.n_alpha();
	guess.n_beta = scf.n_beta();
	guess.D = settings.hf_type == occ::qm::SpinorbitalKind::Unrestricted
		? occ::Mat(occ::qm::block::a(scf_small.ctx.mo.D) + occ::qm::block::b(scf_small.ctx.mo.D))
		: scf_small.ctx.mo.D;

	scf.update_occupied_orbital_count();
	scf.set_core_matrices();
	scf.ctx.F = scf.ctx.H;
	scf.set_conditioning_orthogonalizer();
	scf.ctx.F += scf.m_procedure.compute_fock_mixed_basis(guess, small_bs, false);
	scf.ctx.orthogonalizer.orthogonalize_molecular_orbitals(scf.ctx.mo, scf.ctx.F);
	scf.m_have_initial_guess = true;
}

void XCW::setup_basis(occ::core::Molecule& mol, std::string& basis_set_name, occ::qm::AOBasis& occ_basis_set) {
	std::shared_ptr<BasisSet> basis_set = BasisSetLibrary::get_basis_set(basis_set_name);
	occ_basis_set = basis_set->to_AOBasis(mol.atoms());
}

double XCW::dynamic_damping(const occ::qm::SCF<occ::qm::HartreeFock>& scf, const double& current_alpha, const double& quant_diff, double& quant_diff_mem) {
	double new_alpha = current_alpha;
	if (quant_diff < quant_diff_mem / 10) {
		new_alpha *= 0.75;
		quant_diff_mem = quant_diff;
		if (quant_diff < 10 * scf.convergence_settings.energy_threshold) {
			print_centered_message("***Turned off damping***", 84, XCW_log);
			new_alpha = 0;
			settings.apply_damping = false;
		}
		else {
			std::stringstream print_;
			print_ << "***Decreased damping to " << std::fixed << std::setprecision(3) << new_alpha << "***";
			print_centered_message(print_.str(), 84, XCW_log);
		}
	}
	return new_alpha;
	// closing function
}

void XCW::apply_level_shift(const occ::Mat& C_old, const occ::qm::SCF<occ::qm::HartreeFock>& scf, occ::Mat& F_diis) {
	const int nocc = scf.ctx.mo.Cocc.cols();
	if (scf.ctx.mo.kind == occ::qm::SpinorbitalKind::Restricted) {
		const occ::Mat SC_virt = scf.ctx.S * C_old.rightCols(cryst.nmo - nocc);
		F_diis.noalias() += settings.level_shift * SC_virt * SC_virt.transpose();
	}
	else {
		const int nao = C_old.rows() / 2;
		const auto S_ao = scf.ctx.S.topRows(nao);
		//Cocc has max(n_alpha, n_beta) columns, the spin counts differ for open shells
		const occ::Mat SC_virt_a = S_ao * C_old.topRows(nao).rightCols(cryst.nmo - scf.ctx.mo.n_alpha);
		const occ::Mat SC_virt_b = S_ao * C_old.bottomRows(nao).rightCols(cryst.nmo - scf.ctx.mo.n_beta);
		F_diis.topRows(nao).noalias() += settings.level_shift * SC_virt_a * SC_virt_a.transpose();
		F_diis.bottomRows(nao).noalias() += settings.level_shift * SC_virt_b * SC_virt_b.transpose();
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

bool XCW::do_SCF(const double& lambda, double& alpha, occ::qm::SCF<occ::qm::HartreeFock>& scf, occ::qm::Wavefunction& last_wfn, bool& has_guess, const bool write_result) {

	settings.clear();
	diis_F_.clear();
	diis_E_.clear();
	adiis_.reset();
	ediis_.reset();
	best_quant_ = std::numeric_limits<double>::infinity();
	rescues_ = 0;
	soscf_reset();

	XCW_log << "Starting XCW SCF solver with lambda = " << std::fixed << std::setprecision(5) << lambda << "\n";
	XCW_log << "____________________________________________________________________________________\n";
	XCW_log << " Iteration\t\tCriterion\tGooF(F^2)\tR1(gt)\t\tTotal Energy\t\tPerturbation\tTarget quantity\n";
	XCW_log << "\t\t\t\t\t\t\t\t\t(Eh)\t\t\t(a. u.)\t\t(a. u.)\n";
	XCW_log << "____________________________________________________________________________________\n";

	// Compute first guess and update the energy according to this guess
	const _time_point guess_t0 = get_time();
	if (has_guess) {
		scf.set_initial_guess_from_wfn(last_wfn);
	}
	else {
		if (settings.guess_basis_name.empty())
			scf.compute_initial_guess();
		else
			small_basis_guess(scf);
		has_guess = true;
	}
	scf.ctx.K = scf.m_procedure.compute_schwarz_ints();
	scf.update_scf_energy(false);
	G_last_.resize(0, 0);
	last_full_build_ = 0;
	next_full_build_error_ = 0.0;

	scf.ctx.H = scf.ctx.T + scf.ctx.V;
	throughput::record_time("XCW initial guess", false, get_msec(guess_t0, get_time()));
	bool converged;
	double quant;
	double last_quant = 0;
	double quant_diff_mem = 0;
	occ::Mat dm_last = scf.ctx.mo.D;

	do {
		const _time_point scf_t0 = get_time();
		converged = SCF_iteration(scf, lambda, alpha, quant_diff_mem, quant, last_quant, dm_last);
		throughput::record_time("XCW SCF iteration", false, get_msec(scf_t0, get_time()));
	} while (!converged && scf.iter < scf.maxiter);

	if (converged) {
		XCW_log << "____________________________________________________________________________________\n";
		std::stringstream print_;
		print_ << "***SCF converged in " << scf.iter << " iterations***";
		print_centered_message(print_.str(), 84, XCW_log);

		//Before the summary line below so its A^2 can be appended as an extra column
		if (opt->xcw_gaussian_halt) {
			evaluate_gaussian_halting(lambda);
		}


		const double current_criterion = criterion(false);
		std::cout << std::fixed << std::setprecision(5) << lambda << "\t\t" << std::fixed << std::setprecision(4) << current_criterion << "\t\t" << cryst.GooF2 << "\t\t" << std::setprecision(5) << cryst.R1 << "\t\t" << std::fixed << std::setprecision(9) << scf.ctx.energy["total"] << "\t\t" << std::fixed << std::setprecision(3) << lambda * current_criterion << "\t\t" << std::fixed << std::setprecision(9) << quant
			<< "\t\t" << std::setprecision(4) << criterion(true) << "\t\t" << std::setprecision(5) << cryst.R1_all;
		if (opt->xcw_gaussian_halt && !gaussian_halt_history_.empty()) {
			std::cout << "\t\t" << std::setprecision(4) << gaussian_halt_history_.back().A2;
		}
		std::cout << std::endl;

		if (write_result) {
			const _time_point tscb_t0 = get_time();
			create_tscb(scf, lambda);
			throughput::record_time("XCW tscb", false, get_msec(tscb_t0, get_time()));
		}
	}
	else {
		XCW_log << "____________________________________________________________________________________\n";
		print_centered_message("***SCF did not converge***", 84, XCW_log);
		std::ostringstream perturbed_energy;
		perturbed_energy << " for perturbed energy: " << std::scientific << settings.quant_diff << " (current: " << settings.current_quant_diff << ") \n";
		std::ostringstream diis_error;
		diis_error << " for DIIS error: " << std::scientific << settings.max_diis_error << " (current: " << settings.current_max_diis_error << ") \n";
		std::ostringstream orbital_gradient;
		orbital_gradient << " for orbital gradient: " << std::scientific << settings.gradient << " (current: " << settings.current_gradient << ") \n";
		std::ostringstream max_density_diff;
		max_density_diff << " for maximum difference in density matrix: " << std::scientific << settings.MaxP_diff << " (current: " << settings.current_MaxP_diff << ") \n";
		std::ostringstream rmsd_density;
		rmsd_density << " for RMSD of density matrix: " << std::scientific << settings.RMSP_diff << " (current: " << settings.current_RMSP_diff << ") \n";
		if (settings.conv_quant_diff) {
			XCW_log << "CONVERGED";
		}
		else {
			XCW_log << "NOT CONVERGED";
		}
		XCW_log << perturbed_energy.str();
		if (settings.conv_max_diis_error) {
			XCW_log << "CONVERGED";
		}
		else {
			XCW_log << "NOT CONVERGED";
		}
		XCW_log << diis_error.str();
		if (settings.conv_gradient) {
			XCW_log << "CONVERGED";
		}
		else {
			XCW_log << "NOT CONVERGED";
		}
		XCW_log << orbital_gradient.str();
		if (settings.conv_MaxP_diff) {
			XCW_log << "CONVERGED";
		}
		else {
			XCW_log << "NOT CONVERGED";
		}
		XCW_log << max_density_diff.str();
		if (settings.conv_RMSP_diff) {
			XCW_log << "CONVERGED";
		}
		else {
			XCW_log << "NOT CONVERGED";
		}
		XCW_log << rmsd_density.str();
	}
	return converged;
}

double XCW::compute_orbital_gradient(const occ::qm::SCF<occ::qm::HartreeFock>& scf) {
	if (settings.hf_type == occ::qm::SpinorbitalKind::Restricted) {
		const occ::Mat& C = scf.molecular_orbitals().C;
		const occ::Mat& Cocc = scf.molecular_orbitals().Cocc;
		const occ::Mat Cvir = C.rightCols(C.cols() - Cocc.cols());
		return 2.0 * gemm(Cvir, gemm(scf.ctx.F, Cocc), true).norm();
	}
	else if (settings.hf_type == occ::qm::SpinorbitalKind::Unrestricted) {
		occ::Mat C_alpha = scf.molecular_orbitals().C.topRows(cryst.nmo);
		occ::Mat C_beta = scf.molecular_orbitals().C.bottomRows(cryst.nmo);
		occ::Mat Cocc_alpha = C_alpha.leftCols(scf.ctx.mo.n_alpha);
		occ::Mat Cocc_beta = C_beta.leftCols(scf.ctx.mo.n_beta);
		occ::Mat Cvir_alpha = C_alpha.rightCols(C_alpha.cols() - Cocc_alpha.cols());
		occ::Mat Cvir_beta = C_beta.rightCols(C_beta.cols() - Cocc_beta.cols());
		occ::Mat G_alpha = Cvir_alpha.transpose() * scf.ctx.F.topRows(cryst.nmo) * Cocc_alpha;
		occ::Mat G_beta = Cvir_beta.transpose() * scf.ctx.F.bottomRows(cryst.nmo) * Cocc_beta;
		return (std::hypot(G_alpha.norm(), G_beta.norm()));
	}
	err_not_impl_f("Orbital gradient for a general spinorbital kind", std::cout);
	return 0.0;
}

//F C = e S C through the orthogonaliser's X, one spin block at a time: X^T F X diagonalised
//by MKL, C = X C'. occ's MolecularOrbitals::update does the same with Eigen, which occ
//compiles serial; the occupation, smearing and density steps stay occ's.
void XCW::solve_orbitals(occ::qm::SCF<occ::qm::HartreeFock>& scf, const occ::Mat& F) const {
	occ::qm::MolecularOrbitals& mo = scf.ctx.mo;
	const occ::Mat& X = scf.ctx.orthogonalizer.transformation_matrix();
	const int n = static_cast<int>(X.rows()), m = static_cast<int>(X.cols()), nb = mo.kind == occ::qm::SpinorbitalKind::Unrestricted ? 2 : 1, ldf = static_cast<int>(F.rows());
	occ::Mat FX(n, m), Fp(m, m);
	mo.C.resize(static_cast<Eigen::Index>(nb) * n, m);
	mo.energies.resize(static_cast<Eigen::Index>(nb) * m);
	for (int b = 0; b < nb; b++) {
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, m, n, 1.0, F.data() + b * n, ldf, X.data(), n, 0.0, FX.data(), n);
		cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, m, m, n, 1.0, X.data(), n, FX.data(), n, 0.0, Fp.data(), m);
#if defined(__APPLE__)
		int mm = m, lwork = -1, liwork = -1, iwq = 0, info = 0;
		double wq = 0.0;
		dsyevd_((char*)"V", (char*)"L", &mm, Fp.data(), &mm, mo.energies.data() + b * m, &wq, &lwork, &iwq, &liwork, &info);
		lwork = static_cast<int>(wq), liwork = iwq;
		vec work(lwork);
		ivec iwork(liwork);
		dsyevd_((char*)"V", (char*)"L", &mm, Fp.data(), &mm, mo.energies.data() + b * m, work.data(), &lwork, iwork.data(), &liwork, &info);
		err_checkf(info == 0, "Fock diagonalisation failed", std::cout);
#else
		err_checkf(LAPACKE_dsyevd(LAPACK_COL_MAJOR, 'V', 'L', m, Fp.data(), m, mo.energies.data() + b * m) == 0, "Fock diagonalisation failed", std::cout);
#endif
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, m, m, 1.0, X.data(), n, Fp.data(), m, 0.0, mo.C.data() + b * n, nb * n);
	}
	mo.update_occupied_orbitals();
	mo.smearing.smear_orbitals(mo);
	mo.update_density_matrix();
}

//occ's ConvergenceAccelerator::update with the CDIIS commutator S D F - F D S from MKL; with
//the three matrices symmetric it is T - T^T for T = S D F. CDIIS is cdiis_extrapolate below,
//ADIIS and EDIIS are occ's.
occ::Mat XCW::diis_update(occ::qm::SCF<occ::qm::HartreeFock>& scf) {
	const occ::Mat& S = scf.ctx.S;
	const occ::Mat& D = scf.ctx.mo.D;
	const occ::Mat& F = scf.ctx.F;
	const int n = static_cast<int>(D.cols()), nb = static_cast<int>(D.rows()) / n;
	occ::Mat comm(nb * n, n), T(n, n), SD(n, n);
	for (int b = 0; b < nb; b++) {
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, S.data() + b * n, nb * n, D.data() + b * n, nb * n, 0.0, SD.data(), n);
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, SD.data(), n, F.data() + b * n, nb * n, 0.0, T.data(), n);
		comm.middleRows(static_cast<Eigen::Index>(b) * n, n) = T - T.transpose();
	}
	scf.diis_error = comm.array().abs().maxCoeff();
	occ::Mat F_cdiis = cdiis_extrapolate(F, comm);
	const occ::qm::DiisStrategy strategy = scf.convergence_settings.diis_strategy;
	if (strategy == occ::qm::DiisStrategy::CDIIS || scf.diis_error <= scf.convergence_settings.diis_switch_threshold) return F_cdiis;
	if (strategy == occ::qm::DiisStrategy::ADIIS_CDIIS) return adiis_.update(scf.ctx.mo.kind, D, F);
	return ediis_.update(scf.ctx.mo.kind, D, F, scf.ctx.energy["electronic"]);
}

//Pulay's CDIIS on the Fock matrix, extrapolating from the second vector on: the perturbed
//Roothaan map is stiff and the plain steps occ's DIIS takes first already leave the linear
//regime. B c = 1 over B_ij = <E_i|E_j>, scaled by its largest diagonal element, is solved by
//the pseudo-inverse over the eigenvalues above eigenvalue_cutoff: near convergence the errors
//are almost collinear and a QR solve returns coefficients in the hundreds that amplify the
//noise of the stored Fock matrices. Should one still exceed max_coefficient the oldest vector
//is dropped and the system solved again.
occ::Mat XCW::cdiis_extrapolate(const occ::Mat& F, const occ::Mat& E) {
	diis_F_.push_back(F);
	diis_E_.push_back(E);
	if (diis_F_.size() > diis_subspace_) {
		diis_F_.pop_front();
		diis_E_.pop_front();
	}
	constexpr double eigenvalue_cutoff = 1e-14, max_coefficient = 100.0;
	while (diis_F_.size() > 1) {
		const int n = static_cast<int>(diis_F_.size());
		occ::Mat B(n, n);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j <= i; j++) {
				B(i, j) = B(j, i) = diis_E_[i].cwiseProduct(diis_E_[j]).sum();
			}
		}
		const double scale = B.diagonal().maxCoeff();
		if (!(scale > 0.0)) return F;
		B /= scale;
		const Eigen::SelfAdjointEigenSolver<occ::Mat> es(B);
		const occ::Vec& w = es.eigenvalues();
		const occ::Mat& V = es.eigenvectors();
		//c = B^+ 1 / (1^T B^+ 1) minimises c^T B c under sum(c) = 1
		occ::Vec c = occ::Vec::Zero(n);
		for (int k = 0; k < n; k++) {
			if (w(k) <= eigenvalue_cutoff * w(n - 1)) continue;
			c += (V.col(k).sum() / w(k)) * V.col(k);
		}
		const double norm = c.sum();
		if (std::abs(norm) > 0.0 && c.cwiseAbs().maxCoeff() <= max_coefficient * std::abs(norm)) {
			c /= norm;
			occ::Mat F_out = c(0) * diis_F_[0];
			for (int i = 1; i < n; i++) F_out += c(i) * diis_F_[i];
			return F_out;
		}
		diis_F_.pop_front();
		diis_E_.pop_front();
	}
	return F;
}

//The change from the density of the previous iteration, dm_last, to the one this iteration
//built its Fock matrix from
void XCW::get_density_criteria(double& RMSP_diff, double& maxP_diff, const occ::Mat& dm, const occ::Mat& dm_last) {
	occ::Mat difference = dm - dm_last;
	RMSP_diff = std::sqrt(difference.squaredNorm() / difference.size());
	maxP_diff = difference.cwiseAbs().maxCoeff();
	// closing function
}

bool XCW::SCF_iteration(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double& lambda, double& alpha, double& quant_diff_mem, double& quant, double& last_quant, occ::Mat& dm_last) {
	// Set up energy values & crystallographic information
	scf.iter++;
	const occ::Mat dm_old = scf.ctx.mo.D;
	dMatrix2 dm_eff(cryst.nmo, cryst.nmo);
	//This block is NoSpherA2 code, the Fock build below is OCC. Without the split the whole
	//remainder looks equally ours.
	const _time_point it_t0 = get_time();
	build_effective_dm(scf, dm_eff, dm_old);
	calc_F_calc(dm_eff);
	eval_scale();
	calc_criteria();

	// Generates the perturbation matrix
	occ::Mat perturbation;
	calc_perturb(perturbation, scf);
	const _time_point it_t1 = get_time();
	throughput::record_time("XCW structure factors + perturbation", i_on_device_, get_msec(it_t0, it_t1));

	// Build perturbed Fock matrix
	// Maybe necessary to update the Hamiltoian if a potential changes depending on the density, but that does not happen in normal HF
	//scf.ctx.H = scf.ctx.T + scf.ctx.V + scf.ctx.Vecp + scf.ctx.V_ext;
	scf.m_procedure.update_core_hamiltonian(scf.ctx.mo, scf.ctx.H);
	scf.ctx.F = scf.ctx.H;
	const _time_point fock_t0 = get_time();
	//G(D) = G(D_last) + G(D - D_last): the direct kernel screens shell quartets on the density's
	//shell-block norms and the stored one skips integral segments on the difference's elements,
	//and the difference shrinks as the SCF converges. Rebuilt in full every 8 iterations or once
	//the DIIS error has fallen tenfold since the last full build, as OCC's own loop does, so the
	//screening error does not accumulate. A device that holds the integrals contracts the whole
	//density each time, nothing to skip.
	//Off once TRAH is steering: a step is accepted or rejected on a rise of E + lambda chi^2
	//against the value at the orbitals it left, and trah_.noise puts that threshold at 1e-8 Eh,
	//far below the screening error of a difference build. Comparing one against the other made
	//good steps read as rises, and a rejection can only shrink the trust radius, so a single
	//artefact capped it for the rest of the lambda step. A rotation moves the density by a whole
	//step anyway, so there was little left to skip. OCC's own loop does the same (it forces a
	//full rebuild while the second-order step is active).
	const bool incremental = opt->xcw_incremental && !soscf_ && !eri_on_device_ && scf.m_procedure.fock_build_properties().density_screened && G_last_.size() > 0
		&& scf.iter - last_full_build_ < 8 && scf.diis_error > next_full_build_error_;
	if (incremental) {
		occ::Mat D_diff = scf.ctx.mo.D - D_last_build_;
		std::swap(scf.ctx.mo.D, D_diff);
		G_last_ += eri_ ? eri_fock(scf.ctx.mo, true) : scf.m_procedure.compute_fock(scf.ctx.mo, scf.ctx.K);
		std::swap(scf.ctx.mo.D, D_diff);
	}
	else {
		G_last_ = eri_ ? eri_fock(scf.ctx.mo, false) : scf.m_procedure.compute_fock(scf.ctx.mo, scf.ctx.K);
		last_full_build_ = scf.iter;
		next_full_build_error_ = scf.diis_error / 10.0;
	}
	D_last_build_ = scf.ctx.mo.D;
	scf.ctx.F += G_last_;
	throughput::record_time("XCW Fock build (OCC)", eri_on_device_, get_msec(fock_t0, get_time()));
	scf.update_scf_energy(false);

	const double current_criterion = criterion(false);
	const double temp_penalty = current_criterion * lambda;
	quant = scf.ctx.energy["total"] + temp_penalty;

	scf.ctx.F += perturbation * lambda;

	// Prints output line for iteration
	XCW_log << "\t" << scf.iter << "\t\t" << std::fixed << std::setprecision(4) << current_criterion << "\t\t" << cryst.GooF2 << "\t\t" << std::setprecision(5) << cryst.R1 << "\t\t" << std::fixed << std::setprecision(9) << scf.ctx.energy["total"] << "\t\t" << std::fixed << std::setprecision(3) << temp_penalty << "\t\t" << std::fixed << std::setprecision(9) << quant << std::endl;

	//calc_perturb is the gradient of lambda * criterion^2 (the chi^2 of Jayatilaka's functional),
	//so that, not the printed lambda * criterion, is what the SCF descends and what the rescue ranks by
	const double phi = scf.ctx.energy["total"] + lambda * current_criterion * current_criterion;
	if (!soscf_ && rescue_scf(scf, phi, alpha)) return false;

	// DIIS extrapolation
	occ::Mat F_diis = diis_update(scf);
	settings.current_max_diis_error = scf.diis_error;
	settings.update(XCW_log, alpha);

	// Convergence check
	settings.current_gradient = compute_orbital_gradient(scf);
	settings.current_quant_diff = std::abs(quant - last_quant);
	if (SCF_convergence_check(scf, dm_last)) {
		return true;
	}
	last_quant = quant;

	if (!soscf_) {
		if (scf.iter == 1 || settings.current_gradient < 0.5 * soscf_patience_grad_) {
			soscf_patience_grad_ = settings.current_gradient;
			soscf_patience_iter_ = scf.iter;
		}
		//with soscf requested the DIIS stage only has to reach the quadratic region; Fe_phen HS lambda 0.08
		//oscillated for 30 iterations above the 1e-2 gate while TRAH converged in 6 from that very point
		const int patience = settings.soscf ? trah_.patience_requested : trah_.patience;
		const bool stuck = scf.iter - soscf_patience_iter_ >= patience;
		if (stuck || (settings.soscf && scf.diis_error < trah_.start_threshold)) {
			soscf_ = true;
			std::ostringstream what;
			what << "***" << (stuck ? "Orbital gradient not halved in " + std::to_string(patience) + " iterations" : "DIIS error below 1e-2")
				<< ": second-order steps on the orbital rotations from here***";
			print_centered_message(what.str(), 84, XCW_log);
		}
	}
	if (soscf_) {
		soscf_step(scf, lambda, phi);
		dm_last = dm_old;
		return false;
	}

	// Apply level shift
	if (settings.apply_shift) {
		const occ::Mat& C_old = scf.ctx.mo.C;
		apply_level_shift(C_old, scf, F_diis);
	}

	// Solves central eigenvalue problem
	solve_orbitals(scf, F_diis);

	// Apply damping
	if (settings.apply_damping) {
		if (scf.iter == 2) {
			quant_diff_mem = settings.current_quant_diff;
		}
		if (scf.iter > 2) {
			alpha = dynamic_damping(scf, alpha, settings.current_quant_diff, quant_diff_mem);
		}
	}

	scf.ctx.mo.D *= (1 - alpha);
	scf.ctx.mo.D += alpha * dm_old;
	dm_last = dm_old;
	return false;

	// Apply damping
	scf.ctx.mo.D *= (1 - alpha);
	scf.ctx.mo.D += alpha * dm_old;

	return false;

	//closing function
}

void XCW::soscf_reset() {
	soscf_ = false;
	soscf_patience_iter_ = 0;
	soscf_patience_grad_ = 0;
	trah_B_.clear();
	trah_HB_.clear();
	soscf_kappa_.resize(0);
	soscf_grad_.resize(0);
	soscf_hdiag_.resize(0);
	soscf_C_.resize(0, 0);
	soscf_phi_ = std::numeric_limits<double>::infinity();
	soscf_pred_ = 0;
	soscf_boundary_ = false;
	soscf_floored_ = 0;
	//The radius the last lambda step earned carries into this one - consecutive lambda steps are
	//nearly the same problem, which is the point of ramping lambda at all, and re-learning the
	//radius from 0.5 costs a rejected macro step, hence up to micro_max Fock builds, per halving.
	//Never above the default, so an easy lambda cannot set a hard one up for a fall.
	soscf_trust_ = std::clamp(soscf_trust_, trah_.trust_min, trah_.trust_first);
	trah_micro_total_ = 0;
}

//The gradient of E + lambda chi^2 with respect to the rotation kappa_ai of occupied orbital i
//into virtual a, C -> C exp(kappa), one spin block after the other: 4 F_ai for the restricted
//density 2 C_occ C_occ^T, 2 F_ai per spin block otherwise, F the perturbed Fock matrix in the
//MO basis. diagonal_hessian gets the usual approximation of its diagonal, 4 (e_a - e_i) or
//2 (e_a - e_i) over the diagonal of that F, floored so a near-degenerate pair cannot blow up
//the step.
occ::Vec XCW::orbital_rotation_gradient(const occ::qm::SCF<occ::qm::HartreeFock>& scf, occ::Vec& diagonal_hessian) const {
	const occ::qm::MolecularOrbitals& mo = scf.ctx.mo;
	const int n = cryst.nmo, nb = mo.kind == occ::qm::SpinorbitalKind::Unrestricted ? 2 : 1;
	const double fac = nb == 1 ? 4.0 : 2.0, gap_floor = 0.02;
	Eigen::Index size = 0;
	for (int b = 0; b < nb; b++) {
		const Eigen::Index nocc = static_cast<Eigen::Index>(b == 0 ? mo.n_alpha : mo.n_beta);
		size += nocc * (n - nocc);
	}
	occ::Vec g(size);
	diagonal_hessian.resize(size);
	Eigen::Index at = 0;
	for (int b = 0; b < nb; b++) {
		const Eigen::Index nocc = static_cast<Eigen::Index>(b == 0 ? mo.n_alpha : mo.n_beta), nvir = n - nocc;
		const occ::Mat Cb = mo.C.middleRows(static_cast<Eigen::Index>(b) * n, n);
		const occ::Mat F_mo = Cb.transpose() * scf.ctx.F.middleRows(static_cast<Eigen::Index>(b) * n, n) * Cb;
		const occ::Vec e = F_mo.diagonal();
		for (Eigen::Index i = 0; i < nocc; i++) {
			for (Eigen::Index a = 0; a < nvir; a++, at++) {
				g(at) = fac * F_mo(nocc + a, i);
				diagonal_hessian(at) = fac * std::max(e(nocc + a) - e(i), gap_floor);
			}
		}
	}
	return g;
}

//C_from exp(kappa) by the Cayley transform (I - K/2)^-1 (I + K/2), exactly orthogonal for the
//antisymmetric K that carries kappa in its occupied-virtual blocks, so the orbitals stay
//S-orthonormal. The first nocc columns stay the occupied ones; the density follows.
void XCW::rotate_orbitals(occ::qm::SCF<occ::qm::HartreeFock>& scf, const occ::Mat& C_from, const occ::Vec& kappa) const {
	occ::qm::MolecularOrbitals& mo = scf.ctx.mo;
	const int n = cryst.nmo, nb = mo.kind == occ::qm::SpinorbitalKind::Unrestricted ? 2 : 1;
	Eigen::Index at = 0;
	for (int b = 0; b < nb; b++) {
		const Eigen::Index nocc = static_cast<Eigen::Index>(b == 0 ? mo.n_alpha : mo.n_beta), nvir = n - nocc;
		occ::Mat K = occ::Mat::Zero(n, n);
		for (Eigen::Index i = 0; i < nocc; i++) {
			for (Eigen::Index a = 0; a < nvir; a++, at++) {
				K(nocc + a, i) = 0.5 * kappa(at);
				K(i, nocc + a) = -0.5 * kappa(at);
			}
		}
		const occ::Mat I = occ::Mat::Identity(n, n);
		const occ::Mat U = (I - K).partialPivLu().solve(I + K);
		mo.C.middleRows(static_cast<Eigen::Index>(b) * n, n) = C_from.middleRows(static_cast<Eigen::Index>(b) * n, n) * U;
	}
	mo.update_occupied_orbitals();
	mo.update_density_matrix();
}

//One second-order macro iteration, see soscf_ in the header. phi is E + lambda chi^2 at the
//current orbitals, which the last step produced. Above the functional it left by more than
//the noise of a Fock build, the step is rejected: the step is re-solved at the smaller radius
//in the subspace the micro-iterations already built, from the orbitals it left. The radius
//itself follows occ::qm::trust_radius_update either way - the cube root of the model error the
//step revealed, so the region tracks how far the quadratic model is actually worth trusting,
//and it only moves for a step the region actually stopped or one the model got wrong.
//Two rejections at the smallest radius end the second-order phase and give DIIS the orbitals
//back. Otherwise a fresh gradient starts the next macro step.
void XCW::soscf_step(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double lambda, const double phi) {
	const bool stepped = soscf_kappa_.size() > 0;
	if (stepped) {
		const double actual = phi - soscf_phi_;
		XCW_log << "\t\tTRAH: predicted " << std::scientific << std::setprecision(3) << soscf_pred_
			<< ", actual " << actual << " Eh, rho " << std::fixed << std::setprecision(2)
			<< (soscf_pred_ < 0 ? actual / soscf_pred_ : 0.0) << ", model error "
			<< std::scientific << std::setprecision(1) << std::abs(actual - soscf_pred_) << std::endl;
		soscf_trust_ = occ::qm::trust_radius_update(soscf_trust_, soscf_kappa_.norm(), soscf_pred_, actual, soscf_boundary_, trah_);
	}
	if (stepped && phi > soscf_phi_ + trah_.noise && soscf_kappa_.cwiseAbs().maxCoeff() > 1e-6) {
		//A radius this small is the model saying it is worthless at these orbitals, not that the
		//step should be shorter again. Two in a row and DIIS gets the orbitals back, with the
		//rescue path available again, rather than the lambda step grinding out max_iter on steps
		//too short to move anything.
		soscf_floored_ = soscf_trust_ <= trah_.trust_min * 1.000001 ? soscf_floored_ + 1 : 0;
		if (soscf_floored_ >= 2) {
			std::ostringstream give_up;
			give_up << "***E + lambda chi^2 still rising at the smallest trust radius: back to DIIS***";
			print_centered_message(give_up.str(), 84, XCW_log);
			//hand DIIS the orbitals the last accepted step left, not the rejected ones
			rotate_orbitals(scf, soscf_C_, occ::Vec::Zero(soscf_kappa_.size()));
			soscf_ = false;
			soscf_kappa_.resize(0);
			return;
		}
		std::ostringstream what;
		what << "***E + lambda chi^2 " << std::scientific << std::setprecision(1) << phi - soscf_phi_
			<< " Eh above the orbitals the step left: trust radius " << std::scientific << std::setprecision(2) << soscf_trust_ << ", step re-solved***";
		print_centered_message(what.str(), 84, XCW_log);
		trah_solve(scf, lambda, false);
		rotate_orbitals(scf, soscf_C_, soscf_kappa_);
		return;
	}
	soscf_floored_ = 0;
	//The Roothaan iterations damp the density, so the Fock matrix of the iteration that hands
	//over belongs to a mix of densities, not to the orbitals: rebuilt from them once, so the
	//first gradient, Hessian and reference functional are consistent
	soscf_phi_ = stepped ? phi : rebuild_at(scf, lambda, scf.ctx.mo.C, occ::Vec());
	soscf_grad_ = orbital_rotation_gradient(scf, soscf_hdiag_);
	trah_B_.clear();
	trah_HB_.clear();
	soscf_C_ = scf.ctx.mo.C;
	if (settings.check_hessian && !stepped) check_hessian(scf, lambda);
	trah_solve(scf, lambda, true);
	rotate_orbitals(scf, soscf_C_, soscf_kappa_);
}

//The step of the current macro iteration: the lowest eigenvector of the augmented Hessian
//[[0, a g^T], [a g, H]] in the Davidson subspace trah_B_ (H trah_B_ in trah_HB_), scaled to
//kappa = x / (a x0), with a >= 1 raised by bisection until |kappa| fits the trust radius
//(a = 1 is the plain augmented-Hessian step). With extend, the subspace grows by the
//preconditioned residual of the level-shifted Newton equation (H - theta) kappa = -g, one
//exact Hessian-vector product per micro-iteration, until the residual is below a fraction of
//the gradient that shrinks with it, or trah_.micro_max is reached; without, the retained
//subspace is re-solved at the current radius (after a rejected step the Fock matrix belongs
//to the rejected orbitals, so no product could be added). Sets soscf_kappa_, the predicted
//decrease g.kappa + kappa.H kappa / 2 and whether the step lies on the boundary.
void XCW::trah_solve(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double lambda, const bool extend) {
	const occ::Vec& g = soscf_grad_;
	const double gnorm = g.norm();
	const double tol = gnorm * std::min(0.1, std::max(1e-3, std::sqrt(gnorm)));
	auto add = [&](occ::Vec b) {
		for (int pass = 0; pass < 2; pass++)
			for (const occ::Vec& Bk : trah_B_) b -= Bk.dot(b) * Bk;
		const double nb = b.norm();
		if (nb < 1e-10) return false;
		b /= nb;
		trah_HB_.push_back(hessian_vector(scf, lambda, b));
		trah_B_.push_back(b);
		return true;
	};
	if (extend && trah_B_.empty()) add(-g.cwiseQuotient(soscf_hdiag_));
	//the small problem at shift parameter a: theta, x0, xs and |kappa|
	occ::Mat Hs;
	occ::Vec gs;
	double theta = 0, x0 = 1, alpha = 1;
	occ::Vec xs;
	auto reduced = [&](const double a, double& th, double& v0, occ::Vec& v) {
		const Eigen::Index m = gs.size();
		occ::Mat A = occ::Mat::Zero(m + 1, m + 1);
		A(0, 0) = 0;
		A.block(0, 1, 1, m) = a * gs.transpose();
		A.block(1, 0, m, 1) = a * gs;
		A.block(1, 1, m, m) = Hs;
		Eigen::SelfAdjointEigenSolver<occ::Mat> es(A);
		th = es.eigenvalues()(0);
		const occ::Vec ev = es.eigenvectors().col(0);
		v0 = ev(0);
		v = ev.tail(m);
		return std::abs(v0) > 1e-14 ? v.norm() / (a * std::abs(v0)) : std::numeric_limits<double>::infinity();
	};
	int micro = 0;
	double knorm = 0, rnorm = 0;
	occ::Vec kappa, Hkappa;
	for (;;) {
		const Eigen::Index m = static_cast<Eigen::Index>(trah_B_.size());
		Hs.resize(m, m);
		gs.resize(m);
		for (Eigen::Index i = 0; i < m; i++) {
			gs(i) = trah_B_[i].dot(g);
			for (Eigen::Index j = 0; j <= i; j++) Hs(i, j) = Hs(j, i) = 0.5 * (trah_B_[i].dot(trah_HB_[j]) + trah_B_[j].dot(trah_HB_[i]));
		}
		alpha = 1;
		knorm = reduced(alpha, theta, x0, xs);
		soscf_boundary_ = knorm > soscf_trust_;
		if (soscf_boundary_) {
			//|kappa| falls monotonically with a: bracket, then bisect on log a
			double lo = 1, hi = 1;
			while (hi < 1e8 && reduced(hi, theta, x0, xs) > soscf_trust_) { lo = hi; hi *= 2; }
			for (int it = 0; it < 60 && hi / lo > 1 + 1e-6; it++) {
				const double mid = std::sqrt(lo * hi);
				if (reduced(mid, theta, x0, xs) > soscf_trust_) lo = mid; else hi = mid;
			}
			alpha = hi;
			knorm = reduced(alpha, theta, x0, xs);
		}
		kappa = occ::Vec::Zero(g.size());
		Hkappa = occ::Vec::Zero(g.size());
		for (Eigen::Index i = 0; i < m; i++) {
			kappa += xs(i) * trah_B_[i];
			Hkappa += xs(i) * trah_HB_[i];
		}
		if (std::isfinite(knorm)) {
			kappa /= alpha * x0;
			Hkappa /= alpha * x0;
			//the bisection tolerance may leave the step a hair outside: scale, do not resolve
			if (knorm > soscf_trust_) {
				kappa *= soscf_trust_ / knorm;
				Hkappa *= soscf_trust_ / knorm;
				knorm = soscf_trust_;
			}
		}
		else {
			//no finite step from the subspace: the preconditioned gradient, scaled into the radius.
			//Its curvature may come from a Fock build only while the Fock matrix still belongs to
			//these orbitals - after a rejected step (!extend) the diagonal is all there is
			kappa = -g.cwiseQuotient(soscf_hdiag_);
			kappa *= soscf_trust_ / kappa.norm();
			Hkappa = extend ? hessian_vector(scf, lambda, kappa) : soscf_hdiag_.cwiseProduct(kappa);
			knorm = soscf_trust_;
			soscf_boundary_ = true;
		}
		const occ::Vec r = g + Hkappa - theta * kappa;
		rnorm = r.norm();
		if (!extend || rnorm < tol || micro >= trah_.micro_max) break;
		if (!add(-r.cwiseQuotient((soscf_hdiag_.array() - theta).max(1e-2).matrix()))) break;
		micro++;
	}
	trah_micro_total_ += micro;
	soscf_kappa_ = kappa;
	soscf_pred_ = g.dot(kappa) + 0.5 * kappa.dot(Hkappa);
	XCW_log << "\t\tTRAH: " << micro << " micro-iterations (" << trah_micro_total_ << " in this lambda step), |g| " << std::scientific << std::setprecision(1) << gnorm
		<< ", residual " << rnorm << ", |kappa| " << knorm << (soscf_boundary_ ? " on" : " within") << " the trust radius " << std::fixed << std::setprecision(3) << soscf_trust_
		<< ", shift " << std::scientific << std::setprecision(1) << theta << ", predicted " << soscf_pred_ << " Eh" << std::endl;
}

//H v for the rotation direction v, the derivative of orbital_rotation_gradient along it:
//dC_occ = C_vir V, dC_vir = -C_occ V^T from C exp(kappa); dD from the density's own
//convention (C_occ C_occ^T, halved per spin block for UHF); dF = G(dD) by the same Fock
//build the SCF uses, plus lambda times the response of the perturbation, which is the
//derivative of calc_perturb's per-reflection scalar with the scale moving along (the scale
//minimises chi^2, so the gradient does not see it but the Hessian does); and
//H v = fac (dC_vir^T F C_occ + C_vir^T dF C_occ + C_vir^T F dC_occ) with F the perturbed
//Fock matrix of the current orbitals. The occupied-occupied and virtual-virtual parts of
//the second-order orbital change leave the functional alone, so this is the exact Hessian.
occ::Vec XCW::hessian_vector(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double lambda, const occ::Vec& v) {
	occ::qm::MolecularOrbitals& mo = scf.ctx.mo;
	const int n = cryst.nmo, nb = mo.kind == occ::qm::SpinorbitalKind::Unrestricted ? 2 : 1;
	const double fac = nb == 1 ? 4.0 : 2.0, dfac = nb == 1 ? 1.0 : 0.5;
	occ::Mat dC = occ::Mat::Zero(mo.C.rows(), n), dD = occ::Mat::Zero(mo.D.rows(), n);
	Eigen::Index at = 0;
	for (int b = 0; b < nb; b++) {
		const Eigen::Index nocc = static_cast<Eigen::Index>(b == 0 ? mo.n_alpha : mo.n_beta), nvir = n - nocc, row = static_cast<Eigen::Index>(b) * n;
		occ::Mat V(nvir, nocc);
		for (Eigen::Index i = 0; i < nocc; i++)
			for (Eigen::Index a = 0; a < nvir; a++, at++) V(a, i) = v(at);
		const occ::Mat Cocc = mo.C.block(row, 0, n, nocc), Cvir = mo.C.block(row, nocc, n, nvir);
		dC.block(row, 0, n, nocc) = Cvir * V;
		dC.block(row, nocc, n, nvir) = -Cocc * V.transpose();
		const occ::Mat dCocc = dC.block(row, 0, n, nocc);
		dD.middleRows(row, n) = dfac * (dCocc * Cocc.transpose() + Cocc * dCocc.transpose());
	}
	//G(dD): the Fock build reads mo.D, and a screen on the density would drop the small dD
	std::swap(mo.D, dD);
	occ::Mat dF = eri_ ? eri_fock(mo, false) : scf.m_procedure.compute_fock(mo, scf.ctx.K);
	std::swap(mo.D, dD);
	if (lambda != 0) {
		//dF_r = Sum c I_r dD_eff: calc_F_calc's walk with the fixed part zeroed
		dMatrix2 dD_eff(n, n);
		build_effective_dm(scf, dD_eff, dD);
		const cvec F0 = F_calc[0], fixed = F_calc[1];
		std::fill(F_calc[1].begin(), F_calc[1].end(), cdouble(0.0, 0.0));
		calc_F_calc(dD_eff);
		const cvec dFr = F_calc[0];
		F_calc[0] = F0;
		F_calc[1] = fixed;
		//the scale's response from the weighted least-squares scale of eval_scale
		const bool against_F2 = settings.refine_against == 2, weighted = settings.XWR_type == 2;
		const double k = cryst.F_scale, s = k * k;
		const int chunk = 128, nchunk = (cryst.nr_small + chunk - 1) / chunk;
		vec dnum(nchunk), dden(nchunk), den(nchunk);
#pragma omp parallel for schedule(static)
		for (int ch = 0; ch < nchunk; ch++) {
			for (int r = ch * chunk; r < std::min((ch + 1) * chunk, cryst.nr_small); r++) {
				if (!fit_mask_[r]) continue;
				const double w = weighted ? inv_H2_[r] : 1.0, Fa = std::abs(F0[r]);
				const double dFa = Fa > 0 ? std::real(std::conj(F0[r]) * dFr[r]) / Fa : 0.0;
				if (against_F2) {
					const double wi = w / (obs[r].sigma_obs2 * obs[r].sigma_obs2);
					dnum[ch] += wi * 2.0 * Fa * dFa * obs[r].F_obs2;
					dden[ch] += wi * 4.0 * Fa * Fa * Fa * dFa;
					den[ch] += wi * Fa * Fa * Fa * Fa;
				}
				else {
					const double wi = w / (obs[r].sigma_obs * obs[r].sigma_obs);
					dnum[ch] += wi * dFa * obs[r].F_obs;
					dden[ch] += wi * 2.0 * Fa * dFa;
					den[ch] += wi * Fa * Fa;
				}
			}
		}
		double dnum_sum = 0, dden_sum = 0, den_sum = 0;
		for (int ch = 0; ch < nchunk; ch++) { dnum_sum += dnum[ch]; dden_sum += dden[ch]; den_sum += den[ch]; }
		//F: dk = (dN - k dD) / D; F^2: ds = (dN - s dD) / D for s = k^2
		const double dscale = den_sum != 0 ? (dnum_sum - (against_F2 ? s : k) * dden_sum) / den_sum : 0.0;
		//d(q_r) for q_r = prefactor's scale power times calc_perturb's scalar, both moving:
		//F: q = w conj(F) (k^2 - k Fo / |F|) / sigma^2; F^2: q = w conj(F) (s^2 |F|^2 - s Fo^2) / sigma_I^2
		cvec dq(cryst.nr_small);
#pragma omp parallel for
		for (int r = 0; r < cryst.nr_small; r++) {
			if (!fit_mask_[r]) continue;
			const double w = weighted ? inv_H2_[r] : 1.0, Fa = std::abs(F0[r]);
			if (Fa == 0) continue;
			const double dFa = std::real(std::conj(F0[r]) * dFr[r]) / Fa;
			if (against_F2) {
				const double wi = w / (obs[r].sigma_obs2 * obs[r].sigma_obs2), Fo2 = obs[r].F_obs2;
				dq[r] = wi * (std::conj(dFr[r]) * (s * s * Fa * Fa - s * Fo2)
					+ std::conj(F0[r]) * (2.0 * s * dscale * Fa * Fa + 2.0 * s * s * Fa * dFa - dscale * Fo2));
			}
			else {
				const double wi = w / (obs[r].sigma_obs * obs[r].sigma_obs), Fo = obs[r].abs_F_obs;
				dq[r] = wi * (std::conj(dFr[r]) * (k * k - k * Fo / Fa)
					+ std::conj(F0[r]) * (dscale * (2.0 * k - Fo / Fa) + k * Fo * dFa / (Fa * Fa)));
			}
		}
		occ::Mat dP;
		contract_I(dP, dq);
		dP *= (against_F2 ? 4.0 : 2.0) / (cryst.n_fit - settings.n_params) * lambda;
		for (int b = 0; b < nb; b++) dF.middleRows(static_cast<Eigen::Index>(b) * n, n) += dP;
	}
	occ::Vec Hv(v.size());
	at = 0;
	for (int b = 0; b < nb; b++) {
		const Eigen::Index nocc = static_cast<Eigen::Index>(b == 0 ? mo.n_alpha : mo.n_beta), nvir = n - nocc, row = static_cast<Eigen::Index>(b) * n;
		const occ::Mat Cocc = mo.C.block(row, 0, n, nocc), Cvir = mo.C.block(row, nocc, n, nvir);
		const occ::Mat dCocc = dC.block(row, 0, n, nocc), dCvir = dC.block(row, nocc, n, nvir);
		const occ::Mat F = scf.ctx.F.middleRows(row, n), dFb = dF.middleRows(row, n);
		const occ::Mat M = fac * (dCvir.transpose() * F * Cocc + Cvir.transpose() * dFb * Cocc + Cvir.transpose() * F * dCocc);
		for (Eigen::Index i = 0; i < nocc; i++)
			for (Eigen::Index a = 0; a < nvir; a++, at++) Hv(at) = M(a, i);
	}
	return Hv;
}

//The orbitals C_from exp(kappa) (the current ones, density rebuilt from them, for an empty
//kappa) with everything that depends on them rebuilt: structure factors, scale, criteria,
//perturbation, Fock matrix and energy, the way SCF_iteration does it. Returns E + lambda chi^2
//there.
double XCW::rebuild_at(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double lambda, const occ::Mat& C_from, const occ::Vec& kappa) {
	if (kappa.size() > 0) rotate_orbitals(scf, C_from, kappa);
	else {
		scf.ctx.mo.update_occupied_orbitals();
		scf.ctx.mo.update_density_matrix();
	}
	dMatrix2 dm_eff(cryst.nmo, cryst.nmo);
	build_effective_dm(scf, dm_eff, scf.ctx.mo.D);
	calc_F_calc(dm_eff);
	eval_scale();
	calc_criteria();
	occ::Mat perturbation;
	calc_perturb(perturbation, scf);
	scf.ctx.F = scf.ctx.H + (eri_ ? eri_fock(scf.ctx.mo, false) : scf.m_procedure.compute_fock(scf.ctx.mo, scf.ctx.K));
	scf.update_scf_energy(false);
	const double crit = criterion(false);
	scf.ctx.F += lambda * perturbation;
	return scf.ctx.energy["total"] + lambda * crit * crit;
}

//orbital_rotation_gradient at C_from exp(kappa), rebuilt there and put back afterwards. Only
//the Hessian check pays for it.
occ::Vec XCW::gradient_at(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double lambda, const occ::Mat& C_from, const occ::Vec& kappa) {
	const occ::qm::MolecularOrbitals mo_saved = scf.ctx.mo;
	const occ::Mat F_saved = scf.ctx.F;
	const auto energy_saved = scf.ctx.energy;
	const cvec F0_saved = F_calc[0];
	const double scale_saved = cryst.F_scale;
	rebuild_at(scf, lambda, C_from, kappa);
	occ::Vec h;
	const occ::Vec g = orbital_rotation_gradient(scf, h);
	scf.ctx.mo = mo_saved;
	scf.ctx.F = F_saved;
	scf.ctx.energy = energy_saved;
	F_calc[0] = F0_saved;
	cryst.F_scale = scale_saved;
	return g;
}

//Central finite difference of the gradient along a fixed pseudo-random direction against
//hessian_vector, reported to XCW.log; `check_hessian` runs it on the first second-order step
void XCW::check_hessian(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double lambda) {
	std::mt19937 rng(7);
	std::normal_distribution<double> gauss;
	occ::Vec v(soscf_grad_.size());
	for (Eigen::Index i = 0; i < v.size(); i++) v(i) = gauss(rng);
	v /= v.norm();
	const occ::Vec Hv = hessian_vector(scf, lambda, v);
	const double eps = 1e-4;
	const occ::Mat C0 = scf.ctx.mo.C;
	const occ::Vec fd = (gradient_at(scf, lambda, C0, eps * v) - gradient_at(scf, lambda, C0, -eps * v)) / (2.0 * eps);
	XCW_log << "\t\tHessian check: |Hv - FD| / |FD| = " << std::scientific << std::setprecision(2) << (Hv - fd).norm() / fd.norm()
		<< " (|Hv| " << Hv.norm() << ", |FD| " << fd.norm() << ", v.Hv " << v.dot(Hv) << ", v.FD " << v.dot(fd) << ")" << std::endl;
}

//The rescue of a lost SCF step, see best_mo_ in the header. True when it restarted the step:
//the orbitals are the best ones again, the iteration ends here and the next one rebuilds
//their Fock matrix in full. The stabiliser thresholds stay tightened for the later steps,
//which are at least as hard.
bool XCW::rescue_scf(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double quant, double& alpha) {
	if (quant < best_quant_) {
		best_quant_ = quant;
		best_mo_ = scf.ctx.mo;
		return false;
	}
	if (quant - best_quant_ < rescue_rise_ || rescues_ >= 3) return false;
	rescues_++;
	settings.diis_stop_shift /= 10;
	settings.diis_stop_damping /= 10;
	std::ostringstream what;
	what << "***E + lambda chi^2 " << std::fixed << std::setprecision(3) << quant - best_quant_
		<< " Eh above its best: back to the best orbitals, level shift and damping on until DIIS error "
		<< std::scientific << std::setprecision(0) << settings.diis_stop_shift << " / " << settings.diis_stop_damping
		<< " (rescue " << rescues_ << "/3)***";
	print_centered_message(what.str(), 84, XCW_log);
	scf.ctx.mo = best_mo_;
	G_last_.resize(0, 0);
	diis_F_.clear();
	diis_E_.clear();
	adiis_.reset();
	ediis_.reset();
	settings.apply_shift = settings.method_apply_shift;
	settings.apply_damping = settings.method_apply_damping;
	alpha = settings.alpha;
	return true;
}

bool XCW::SCF_convergence_check(occ::qm::SCF<occ::qm::HartreeFock>& scf, occ::Mat& dm_last) {
	get_density_criteria(settings.current_RMSP_diff, settings.current_MaxP_diff, scf.ctx.mo.D, dm_last);
	settings.conv_quant_diff = settings.current_quant_diff < settings.quant_diff;
	settings.conv_max_diis_error = settings.current_max_diis_error < settings.max_diis_error;
	settings.conv_gradient = settings.current_gradient < settings.gradient;
	settings.conv_RMSP_diff = settings.current_RMSP_diff < settings.RMSP_diff;
	settings.conv_MaxP_diff = settings.current_MaxP_diff < settings.MaxP_diff;
	return settings.convergence_check();
	// closing function
}

//The sign of f(+-3), g(+-3) and g(+-4) in every orbital, in OCC's own m = -l..l order,
//and the density matrix rebuilt from them. Applied once on the way out and once on the
//way back in, it is its own inverse.
void XCW::flip_high_m_phases(occ::qm::Wavefunction& w) {
	int row = 0;
	const int spins = w.mo.kind == occ::qm::SpinorbitalKind::Unrestricted ? 2 : 1;
	for (const auto& shell : w.basis.shells()) {
		const int nsph = 2 * shell.l + 1;
		if (shell.l >= 3)
			for (int spin = 0; spin < spins; spin++) {
				const int base = spin * w.nbf + row;
				w.mo.C.row(base).array() *= -1.0;
				w.mo.C.row(base + nsph - 1).array() *= -1.0;
				if (shell.l >= 4) {
					w.mo.C.row(base + 1).array() *= -1.0;
					w.mo.C.row(base + nsph - 2).array() *= -1.0;
				}
			}
		row += nsph;
	}
	w.mo.update_occupied_orbitals();
	w.mo.update_density_matrix();
}

void XCW::create_tscb(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double& lambda) {
	XCW_log << "Creating .tscb file from converged SCF calculation..." << std::endl;
	std::vector<WFN> sf_wave_vec(1, { scf.wavefunction(), false });
	//The constructor marks anything taken from OCC as OCC-origin. What this refinement holds
	//is an OCC result over a basis this program loaded, so say that: Int_Params then reads the
	//shells with the convention they actually have.
	sf_wave_vec[0].set_origin(e_origin::XCW_fit);
	svec known_atoms_;
	tsc_block<int, cdouble> result;
	vec2 known_kpts_;
	options* opt_ = const_cast<options*>(opt);
	opt_->m_hkl_list = hkl_enlarged;
	opt_->grid_cache = &tsc_grids;
	result.append(calculate_scattering_factors<itsc_block, std::vector<WFN>&>(
		*opt_,
		sf_wave_vec,
		XCW_log,
		known_atoms_,
		0,
		&k_pt),
		XCW_log);
	std::string value = std::to_string(lambda);
	value.erase(std::remove(value.begin(), value.end(), '.'), value.end());
	while (value.length() < 7) {
		value += '0';
	}
	std::ostringstream oss;
	oss << "NA2_" << value << ".tscb";
	result.write_tscb_file("test.cif", oss.str());
	std::ostringstream oss2;
	oss2 << "NA2_" << value << ".wfn";
	sf_wave_vec[0].write_wfn(oss2.str(), false, true);
	//the structure factors this step was scored on, for an R-factor against another route
	//(an ORCA tsc through Olex2, or a Laue check between symmetry-equivalent rows)
	{
		ensure_hkl_ordered();
		std::ofstream fc("NA2_" + value + "_Fcalc.txt");
		fc << "#    h    k    l          F_obs        sig(F)   scale*|F_calc|     phase(deg)   R1(gt) = " << std::setprecision(5) << cryst.R1 << " R1(all) = " << cryst.R1_all << " scale = " << std::setprecision(10) << cryst.F_scale << "\n";
		for (int r = 0; r < cryst.nr_small; r++) {
			const cdouble& f = F_calc[0][r];
			fc << std::setw(5) << hkl_ordered_[r][0] << std::setw(5) << hkl_ordered_[r][1] << std::setw(5) << hkl_ordered_[r][2]
				<< std::fixed << std::setprecision(4) << std::setw(15) << obs[r].F_obs << std::setw(14) << obs[r].sigma_obs
				<< std::setw(17) << cryst.F_scale * std::abs(f) << std::setw(15) << std::arg(f) * 180.0 / constants::PI << "\n";
		}
	}
	std::ostringstream oss3;
	oss3 << "NA2_" << value << ".fchk";
	//OCC's fchk writer reorders to Gaussian's basis functions but keeps libcint's phases,
	//and Gaussian's f(+-3), g(+-3), g(+-4) are the opposite sign; the file has to carry
	//Gaussian's so that anything reading an fchk gets the density right
	{
		occ::qm::Wavefunction w = scf.wavefunction();
		flip_high_m_phases(w);
		w.save(oss3.str());
	}
	if (settings.nbo_output) {
		std::ostringstream oss4;
		oss4 << "NA2_" << value << ".47";
		sf_wave_vec[0].write_nbo(oss4.str(), opt->debug, &XCW_log);
	}
	//Neither file written above can carry this analysis - a .wfn has bare primitives and
	//the fchk reader keeps no shells - so -rgbi runs it here, on the refined wavefunction,
	//and its report goes to a file of its own per lambda
	if (opt->rgbi) {
		std::ostringstream oss5;
		oss5 << "NA2_" << value << "_RGBI.txt";
		std::ofstream rgbi_out(oss5.str());
		std::streambuf* const cout_buf = std::cout.rdbuf(rgbi_out.rdbuf());
		Roby_information Roby(sf_wave_vec[0], opt->rgbi_group_sets, !opt->rgbi_no_sym,
			opt->rgbi_orbital_basis == RGBIOrbitalBasis::ANO, opt->rgbi_EVs, opt->rgbi_theta);
		std::cout.rdbuf(cout_buf);
		XCW_log << "RGBI analysis written to " << oss5.str() << std::endl;
	}
}

//The stored integrals' two-electron Fock part, contracted on the device when it holds them
occ::Mat XCW::eri_fock(const occ::qm::MolecularOrbitals& mo, bool screen) const {
	stored_eri::jk_fn jk;
#if defined(NOSPHERA2_USE_GPU) || defined(NOSPHERA2_USE_METAL)
	if (eri_on_device_) {
		jk = [](const occ::Mat& D, occ::Mat& J, occ::Mat& K) {
			J.resize(D.rows(), D.cols());
			K.resize(D.rows(), D.cols());
			err_checkf(eri_gpu_JK(D.data(), J.data(), K.data()), "Fock build on the device failed", std::cout);
		};
	}
#endif
	return eri_.fock(mo, screen, jk);
}

occ::qm::HartreeFock XCW::setup_XCW_procedure(bool read_tensor) {
	std::vector<ao_data> ao_data_shells;
	occ::core::Molecule mol;
	setup_SCF_mol(mol);
	occ::qm::AOBasis occ_basis_set;
	setup_basis(mol, settings.basis_set_name, occ_basis_set);
	occ::qm::HartreeFock hf(occ_basis_set);
	if (!settings.df_basis_name.empty()) {
		//OCC loads a fitting basis by name from a data directory this build does not ship,
		//but reads a .json path as given: the library's set is written out once, for the
		//elements present, and handed over that way
		std::shared_ptr<BasisSet> aux = BasisSetLibrary::get_basis_set(settings.df_basis_name);
		//a Coulomb-only set fits J and leaves the exchange to a basis never meant for it
		if (aux->get_name().find("jkfit") == std::string::npos)
			std::cout << "WARNING: " << aux->get_name() << " is not a JK-fitting basis; Hartree-Fock exchange is fitted with it all the same. "
				"def2-universal-jkfit serves the def2 family, cc-pvXz-jkfit the cc-pVXZ family." << std::endl;
		ivec elements;
		for (int i = 0; i < static_cast<int>(mol.atoms().size()); i++)
			if (std::find(elements.begin(), elements.end(), mol.atoms()[i].atomic_number) == elements.end())
				elements.push_back(mol.atoms()[i].atomic_number);
		const std::string file = aux->get_name() + "_df.json";
		aux->write_occ_json(file, elements);
		hf.set_density_fitting_basis(file);
		//OCC keeps the three-index integrals only under a 512 MB limit and otherwise recomputes
		//them every iteration, which costs twice a direct build here. Held whenever they fit in
		//half of what the process can have; they are computed once for the whole lambda scan.
		const size_t naux = aux->to_AOBasis(mol.atoms()).nbf();
		const size_t nbf = occ_basis_set.nbf();
		const size_t store = naux * nbf * (nbf + 1) / 2 * sizeof(double);
		const size_t avail = available_memory_bytes();
		const bool stored = avail == 0 || store < avail / 2;
		hf.set_density_fitting_policy(stored ? occ::qm::IntegralEngineDF::Policy::Stored : occ::qm::IntegralEngineDF::Policy::Direct);
		std::cout << "XCW density fitting with " << aux->get_name() << ": " << naux << " functions, "
			<< (store / 1048576.0) << " MB of three-index integrals " << (stored ? "held in memory" : "recomputed every iteration") << std::endl;
	}
	if (opt->xcw_int_precision > 0.0) hf.set_precision(opt->xcw_int_precision);
	create_prims(ao_data_shells, occ_basis_set);
	eval_I_anom_disp(ao_data_shells, read_tensor);
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
	//OCC parallelises through TBB, which does not read OMP_NUM_THREADS. Not a speedup - it
	//already used every core - but it makes -cpus bind the 82% of a run that OCC owns.
	occ::parallel::set_num_threads(opt->threads > 0 ? opt->threads : omp_get_max_threads());
	occ::qm::HartreeFock hf = setup_XCW_procedure(settings.read_tensor);
#if defined(NOSPHERA2_USE_GPU) || defined(NOSPHERA2_USE_METAL)
	//The two walks of an iteration read the whole tensor and the host is bound by its memory
	//bandwidth doing so; the device reads it several times faster. The host copy stays for
	//the background writer.
	if (opt->gpu_itensor && opt->use_gpu && !i_streamed_) {
		i_on_device_ = i_float_ ? itensor_gpu_hold(I32.data(), cryst.nr_small, static_cast<int>(i_compact_))
			: itensor_gpu_hold(I.data(), cryst.nr_small, static_cast<int>(i_compact_));
		if (!(opt->no_date))
			std::cerr << "GPU in use: XCW structure factors and perturbation on "
			<< (i_on_device_ ? "the device" : "the CPU - device unavailable or the tensor too large") << std::endl;
	}
#endif
	//Four fifths of the free memory, the I tensor's budget, for the stored integrals
	eri_.clear();
	if (settings.hf_type != occ::qm::SpinorbitalKind::General) {
		const _time_point eri_t0 = get_time();
		const size_t avail = available_memory_bytes();
		if (eri_.build(hf, avail ? avail / 5 * 4 : 0, XCW_log))
			throughput::record_time("XCW two-electron integrals", false, get_msec(eri_t0, get_time()));
	}
#if defined(NOSPHERA2_USE_GPU) || defined(NOSPHERA2_USE_METAL)
	if (eri_ && opt->gpu_itensor && opt->use_gpu) {
		const auto up_t0 = get_time();
		eri_on_device_ = eri_gpu_hold(eri_.data(), eri_.nbf(), eri_.npairs(), eri_.pair_a().data(), eri_.pair_b().data(),
			eri_.first_pair().data(), eri_.pair_index().data());
		if (eri_on_device_) throughput::record_time("XCW two-electron integrals upload", true, get_msec(up_t0, get_time()));
		if (!(opt->no_date))
			std::cerr << "GPU in use: XCW Fock build from the stored integrals on "
			<< (eri_on_device_ ? "the device" : "the CPU - device unavailable or the integrals too large") << std::endl;
	}
#endif
	occ::qm::SCF scf(hf, settings.hf_type);
	bool has_guess = false;
	occ::qm::Wavefunction last_wfn, prev_wfn;
	const occ::Mat S_ao = hf.compute_overlap_matrix();
	if (settings.read_first_guess) {
		std::ostringstream oss2;
		std::string start_value_str = std::to_string(settings.xcw_start_value);
		start_value_str.erase(std::remove(start_value_str.begin(), start_value_str.end(), '.'), start_value_str.end());
		if (start_value_str.length() > 7) {
			start_value_str = start_value_str.substr(0, 7);
		}
		else if (start_value_str.length() < 7) {
			start_value_str.append(7 - start_value_str.length(), '0');
		}
		oss2 << "NA2_" << start_value_str << ".fchk";
		last_wfn = occ::qm::Wavefunction::load(oss2.str());
		//Written in Gaussian's phases by create_tscb; OCC's loader does not undo that
		flip_high_m_phases(last_wfn);
		has_guess = true;
	}

	std::cout << "More detailed output in XCW.log file..." << std::endl;
	if (settings.XWR_type == 2) {
		std::cout << "XCW: fitting against the 1/|H|^2-weighted residual self-energy criterion "
			<< "Criterion below is this weighted quantity, not the classical GoF." << std::endl;
		XCW_log << "XCW: fitting against the 1/|H|^2-weighted residual self-energy criterion "
			<< "Criterion below are this weighted quantity, not the classical GoF." << std::endl;
	}
	std::cout << "____________________________________________________________________________________\n";
	std::cout << " Lambda\t\tCriterion\tGooF(F2)\tR1(gt)\t\tTotal Energy\t\tPerturbation\tTarget quantity\t\tCrit(all)\tR1(all)";
	if (opt->xcw_gaussian_halt) {
		std::cout << "\t\tA^2 (halt)";
	}
	std::cout << "\n";
	std::cout << "\t\t\t\t\t\t\t\t(Eh)\t\t\t(a. u.)\t\t(a. u.)\n";
	std::cout << "____________________________________________________________________________________\n";

	// Runs the lambda steps for XCW fitting
	auto run_lambda = [&](const double lambda, occ::qm::Wavefunction guess, const bool use_guess, const bool write_result) {
		occ::qm::SCF scf(hf, settings.hf_type);
		double alpha = settings.alpha;
		bool has_local_guess = use_guess;
		scf.set_charge_multiplicity(settings.charge, settings.multiplicity);
		scf.maxiter = settings.max_scf_iterations;
		scf.convergence_settings.level_shift = settings.level_shift;
		scf.convergence_settings.level_shift_threshold = 0;
		scf.update_occupied_orbital_count();
		const bool converged = do_SCF(lambda, alpha, scf, guess, has_local_guess, write_result);
		return std::make_pair(converged, scf.wavefunction());
	};
	const double min_lambda_step = settings.xcw_step_size / 128.0;
	double last_lambda = settings.xcw_start_value;
	//slow_conv's damping (alpha 0.8, shift 1) is for the perturbed steps, which start from
	//converged orbitals; from the core guess it runs out of iterations (Fe(phen)2(SCN)2 UHF:
	//200 iterations, RMSD still 4e-5). The first step without a guess takes the normal
	//schedule, every later one the slow settings from the file
	const bool slow_start = settings.slow_conv && !has_guess;
	const double slow_alpha = settings.alpha, slow_shift = settings.level_shift, slow_stop_damping = settings.diis_stop_damping, slow_stop_shift = settings.diis_stop_shift;
	if (slow_start) {
		settings.alpha = 0.5; settings.level_shift = 0.5; settings.diis_stop_damping = 1e-3; settings.diis_stop_shift = 1e-2;
		std::cout << "XCW: slow_conv - the unperturbed first step runs the normal schedule, slow damping from the second step on" << std::endl;
		XCW_log << "XCW: slow_conv - the unperturbed first step runs the normal schedule, slow damping from the second step on" << std::endl;
	}
	//The scan used to stop at max_value even when A^2 was still falling there, i.e. on the
	//scan boundary instead of on lambda*, and only printed "extend the scan" afterwards.
	//The perturbed steps start from converged orbitals and TRAH takes the ones that DIIS
	//cannot, so running past the requested range is stable enough to just do it. Capped at
	//the requested number of steps, so a trend that never curves upward at most doubles the
	//scan rather than running away
	int planned_steps = settings.num_xcw_steps;
	const int max_extra_steps = settings.num_xcw_steps;
	int extra_steps = 0;
	for (int step = 0; step < planned_steps; step++) {
		const double lambda = step * settings.xcw_step_size + settings.xcw_start_value;
		if (slow_start && step == 1) {
			settings.alpha = slow_alpha; settings.level_shift = slow_shift; settings.diis_stop_damping = slow_stop_damping; settings.diis_stop_shift = slow_stop_shift;
		}
		const occ::qm::Wavefunction previous_wfn = last_wfn;
		occ::qm::Wavefunction guess = last_wfn;
		if (opt->xcw_extrapolate && step >= 2 && settings.hf_type == occ::qm::SpinorbitalKind::Restricted) {
			//The density extrapolated through the two previous steps, pulled back to
			//idempotency by two McWeeny steps D <- 3DSD - 2DSDSD; the orbitals stay those of
			//the last step, they only seed the level shift and the gradient
			occ::Mat D = 2.0 * last_wfn.mo.D - prev_wfn.mo.D;
			for (int k = 0; k < 2; k++) {
				const occ::Mat DS = D * S_ao;
				D = 3.0 * DS * D - 2.0 * DS * DS * D;
			}
			guess.mo.D = D;
		}
		auto result = run_lambda(lambda, guess, has_guess, true);
		if (!result.first && step > 0) {
			double lambda_step = lambda - last_lambda;
			while (!result.first && lambda_step > min_lambda_step) {
				lambda_step *= 0.5;
				const double trial_lambda = std::min(lambda, last_lambda + lambda_step);
				XCW_log << "XCW: retrying lambda " << std::fixed << std::setprecision(8) << trial_lambda
					<< " from converged lambda " << last_lambda << " with step " << lambda_step << std::endl;
				std::cout << "XCW: retrying lambda " << std::fixed << std::setprecision(8) << trial_lambda
					<< " with step " << lambda_step << std::endl;
				auto trial = run_lambda(trial_lambda, last_wfn, true, trial_lambda == lambda);
				if (trial.first) {
					if (trial_lambda == lambda) {
						result = std::move(trial);
						break;
					}
					last_lambda = trial_lambda;
					last_wfn = trial.second;
					lambda_step = std::min(2.0 * lambda_step, lambda - last_lambda);
					result = run_lambda(lambda, last_wfn, true, true);
				}
			}
		}
		if (!result.first) {
			//step 0 has no converged neighbour to continue from, so the halving loop never ran
			std::ostringstream why;
			why << "XCW: unable to converge lambda " << std::fixed << std::setprecision(8) << lambda;
			if (step == 0)
				why << " in " << settings.max_scf_iterations << " SCF iterations (raise max_iter or loosen the criteria); stopping scan.";
			else
				why << " with a continuation step above " << min_lambda_step << "; stopping scan.";
			XCW_log << why.str() << std::endl;
			std::cout << why.str() << std::endl;
			break;
		}
		prev_wfn = previous_wfn;
		last_wfn = result.second;
		last_lambda = lambda;
		has_guess = true;

		//At the end of the planned scan, keep going while the minimum is still ahead
		if (opt->xcw_gaussian_halt && step + 1 == planned_steps) {
			double estimated_minimum = 0.0;
			if (halting_minimum_beyond_scan(gaussian_halt_history_, estimated_minimum)) {
				std::ostringstream msg;
				msg << std::fixed << std::setprecision(5)
					<< "XCW: A^2 is still falling at lambda " << lambda << ", so the minimum is outside the requested range - ";
				if (extra_steps < max_extra_steps) {
					planned_steps++;
					extra_steps++;
					msg << "extending the scan to lambda " << (lambda + settings.xcw_step_size)
						<< " (extra step " << extra_steps << " of at most " << max_extra_steps << ")";
					if (estimated_minimum > 0.0) {
						msg << ", extrapolated minimum near lambda ~= " << estimated_minimum;
					}
				}
				else {
					msg << "stopping anyway, the scan has already been extended by its limit of "
						<< max_extra_steps << " steps; raise max_value in -do_XCW to continue";
				}
				XCW_log << msg.str() << std::endl;
				std::cout << msg.str() << std::endl;
			}
		}

		//Progress estimate every 5 lambda steps; the last step is skipped because the
		//summary below always prints a final one
		if (opt->xcw_gaussian_halt && (step + 1) % 5 == 0 && step + 1 < planned_steps) {
			report_halting_progress_estimate(false);
		}
	}

	if (opt->xcw_gaussian_halt) {
		report_gaussian_halting_summary();
	}

	//Before the run ends: the writer holds a file handle and reads the resident tensor, and
	//by now it has usually been finished for a long while - the refinement takes far longer
	//than the write. Joining is what keeps it from outliving the process.
	finish_i_save();
#if defined(NOSPHERA2_USE_GPU) || defined(NOSPHERA2_USE_METAL)
	itensor_gpu_release();
	eri_gpu_release();
	i_on_device_ = eri_on_device_ = false;
#endif

	std::cout << "Finished XCW fitting procedure." << std::endl;
}
