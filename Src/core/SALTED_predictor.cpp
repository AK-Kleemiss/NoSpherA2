#include "pch.h"
#include "SALTED_predictor.h"
#ifdef NOSPHERA2_USE_GPU
#include "salted_gpu.h"
#include "SALTED_equicomb.h"
#endif
#include "SALTED_utilities.h"
#include <occ/core/eeq.h>
#include "spherical_density.h"
#include "SALTED_equicomb.h"
#include "nos_math.h"
#include "constants.h"
#include "wfn_class.h"
#include "basis_set.h"
#include <filesystem>
#include <future>


SALTEDPredictor::SALTEDPredictor(WFN wavy_in, options& opt_in)
{
	std::filesystem::path _path = opt_in.salted_model_dir;
	SALTED_DIR = opt_in.salted_model_dir;
	debug = opt_in.debug;
	force_charge_constraint = opt_in.salted_charge_constraint;

	if (opt_in.salted_model_dir.empty() && opt_in.coef_file != "") {
		std::cout << "Using density coefficients found in: " << opt_in.coef_file << std::endl;
		wavy = generate_aux_wfn(wavy_in, opt_in.aux_basis);
		bbasis_set_loaded = true;
		config.dfbasis = opt_in.aux_basis[0]->get_name();
		config.salted_filename = "coefficient file";
		coef_file = opt_in.coef_file;
		return;
	}

	if (opt_in.salted_model_dirs.size() > 1) {
		build_merged(wavy_in, opt_in);
		return;
	}

	// A model may be named by its own file: two models often live in the same
	// folder, and then the directory alone cannot say which one is meant.
	if (std::filesystem::is_regular_file(_path)) {
		SALTED_DIR = _path.parent_path();
		config.salted_filename = _path.filename();
		_path = SALTED_DIR;
	}
	else {
		config.salted_filename = find_first_salted_file(opt_in.salted_model_dir);
		if (config.salted_filename.empty()) {
			std::cout << "No SALTED binary file found in directory: " << opt_in.salted_model_dir << std::endl;
			exit(1);
		}
	}

	wavy = wavy_in;
	if (opt_in.debug) std::cout << "Using SALTED Binary file: " << config.salted_filename << std::endl;
	_path = _path / config.salted_filename;
	SALTED_BINARY_FILE file = SALTED_BINARY_FILE(_path);
	file.populate_config(config);

	const std::vector<char> use_thakkar = SALTED_Utils::filter_input(wavy, opt_in, config);
	if (!use_thakkar.empty()) {
		spherical_fill_used = true;
		estimate_fill_charges(wavy_in, use_thakkar, opt_in);
	}

	//wavy.write_xyz("temp_rascaline.xyz"); //Also this
	//config.predict_filename = "temp_rascaline.xyz";

	natoms = wavy.get_ncen();

	if (wavy.get_nmo() != 0)
		wavy.clear_MOs(); // Delete unneccesarry MOs, since we are predicting anyway.

	wavy.delete_basis_set();
	if (file.basis_set_defined()) {
		model_basis = file.read_basis_set();
		// SALTED predicts one coefficient per DECONTRACTED auxiliary function, so the
		// basis has to be loaded that way whatever the default of the day is.
		load_basis_into_WFN(wavy, model_basis, true);
		bbasis_set_loaded = true;
	}
}

void SALTEDPredictor::estimate_fill_charges(const WFN& wavy_in, const std::vector<char>& use_thakkar, options& opt_in)
{
	{
		// The filled atoms get a NEUTRAL Thakkar density, which fixes how many
		// electrons they carry. Estimate what they should really carry, so the
		// size of that assumption can be reported rather than hidden. EEQ gives
		// smooth non-integer charges from geometry and honours the net charge,
		// which suits coordination chemistry far better than assigning a formal
		// oxidation state.
		const int ncen_in = wavy_in.get_ncen();
		n_filled = static_cast<int>(std::count(use_thakkar.begin(), use_thakkar.end(), (char)1));
		try
		{
			occ::IVec nums(ncen_in);
			occ::Mat3N pos(3, ncen_in);
			const bool bohr = wavy_in.get_isBohr();
			for (int a = 0; a < ncen_in; a++)
			{
				nums(a) = wavy_in.get_atom_charge(a);
				for (int ax = 0; ax < 3; ax++)
				{
					const double c = wavy_in.get_atom_coordinate(a, ax);
					pos(ax, a) = bohr ? constants::bohr2ang(c) : c;   // EEQ wants Angstrom
				}
			}
			const occ::Vec q = occ::core::charges::eeq_partial_charges(
				nums, pos, static_cast<double>(wavy_in.get_charge()));
			filled_eeq_charge = 0.0;
			applied_fill_charge = 0.0;
			opt_in.spherical_fill_charges.clear();
			for (int a = 0; a < ncen_in; a++)
			{
				if (!use_thakkar[a]) continue;
				filled_eeq_charge += q(a);
				// Only charge the fill can actually carry may be moved out of
				// the predicted region. If no ion is tabulated for this element
				// the fill stays neutral, so the target must stay neutral too -
				// otherwise the two disagree and the system total is wrong.
				const int Zf = wavy_in.get_atom_charge(a);
				const bool ion_ok = (q(a) > 0.0) ? Thakkar_Cation::available(Zf)
												 : Thakkar_Anion::available(Zf);
				if (ion_ok) applied_fill_charge += q(a);
				// Position-keyed, in the wavefunction's own units: the fill
				// rebuilds its wavefunction from the original file, so indices
				// there are not ours to assume.
				opt_in.spherical_fill_charges.push_back({
					wavy_in.get_atom_coordinate(a, 0),
					wavy_in.get_atom_coordinate(a, 1),
					wavy_in.get_atom_coordinate(a, 2),
					ion_ok ? q(a) : 0.0});
			}
		}
		catch (const std::exception &e)
		{
			std::cout << "Could not estimate the filled-region charge (" << e.what()
					  << "); reporting it as unknown." << std::endl;
			filled_eeq_charge = std::numeric_limits<double>::quiet_NaN();
		}
	}
}

const std::string SALTEDPredictor::get_dfbasis_name() const
{
	return config.dfbasis;
}

std::shared_ptr<BasisSet> SALTEDPredictor::get_model_basis() const
{
	return model_basis ? model_basis : BasisSetLibrary::get_basis_set(config.dfbasis);
}

// Match a predictor's (possibly filtered) structure back onto the full one. The
// filter only erases atoms and leaves the rest in order and untouched, so one
// forward walk with an exact comparison is enough.
static ivec map_atoms(const WFN& full, const WFN& part)
{
	ivec idx(full.get_ncen(), -1);
	int p = 0;
	for (int a = 0; a < full.get_ncen() && p < part.get_ncen(); a++)
	{
		if (full.get_atom_charge(a) != part.get_atom_charge(p)) continue;
		bool same = true;
		for (int ax = 0; ax < 3 && same; ax++)
			same = (full.get_atom_coordinate(a, ax) == part.get_atom_coordinate(p, ax));
		if (same) idx[a] = p++;
	}
	return idx;
}

// Where each atom's coefficients start in a wavefunction's auxiliary basis - the
// same table the density evaluation indexes with. One entry past the end.
static ivec coef_offsets(const WFN& w)
{
	const aux_density_table t(*w.get_atoms_ptr());
	ivec off(w.get_ncen() + 1, t.n_coef);
	for (int a = 0; a < w.get_ncen(); a++) off[a] = t.coef_off[t.sh_start[a]];
	return off;
}

void SALTEDPredictor::build_merged(const WFN& wavy_in, options& opt_in)
{
	// One prediction per model, each of them seeing the whole structure: an atom
	// another model owns is still a neighbour of the ones this model predicts.
	// The blocks are per atom and each carries its own model's auxiliary basis,
	// so putting them side by side is a concatenation, not a mixture.
	std::cout << "Combining " << opt_in.salted_model_dirs.size()
			  << " SALTED models. Each of them first reports what IT alone cannot predict;"
			  << " the assignment that counts is printed below." << std::endl;
	for (const auto& model : opt_in.salted_model_dirs)
	{
		options sub_opt = opt_in;               // its own spherical-fill bookkeeping
		sub_opt.salted_model_dirs.clear();
		sub_opt.salted_model_dir = model;
		sub_opt.needs_Thakkar_fill = false;
		sub_opt.spherical_fill_charges.clear();
		auto sub = std::make_unique<SALTEDPredictor>(wavy_in, sub_opt);
		sub->skip_charge_constraint = true;     // the stitched density is scaled once, at the end
		if (!sub->basis_set_loaded())
			load_basis_into_WFN(sub->wavy, sub->get_model_basis(), true);
		sub_models.push_back(std::move(sub));
	}

	// An element goes to the first model that predicts it, and all its atoms with
	// it: the auxiliary basis lives per element, so two models cannot share one.
	std::vector<ivec> sub_index;
	for (const auto& sub : sub_models) sub_index.push_back(map_atoms(wavy_in, sub->wavy));
	ivec element_model(118, -1);
	for (int a = 0; a < wavy_in.get_ncen(); a++)
	{
		const int Z = wavy_in.get_atom_charge(a) - 1;
		if (element_model[Z] >= 0) continue;
		for (int m = 0; m < (int)sub_models.size(); m++)
			if (sub_index[m][a] >= 0) { element_model[Z] = m; break; }
	}

	wavy = wavy_in;
	atom_model.assign(wavy_in.get_ncen(), -1);
	atom_in_model.assign(wavy_in.get_ncen(), -1);
	for (int a = 0; a < wavy_in.get_ncen(); a++)
	{
		const int m = element_model[wavy_in.get_atom_charge(a) - 1];
		// A model that dropped this atom for having no environment cannot predict
		// it either, however well it knows the element.
		if (m < 0 || sub_index[m][a] < 0) continue;
		atom_model[a] = m;
		atom_in_model[a] = sub_index[m][a];
	}

	for (int m = 0; m < (int)sub_models.size(); m++)
	{
		svec elements;
		int n = 0;
		for (int Z = 0; Z < 118; Z++)
			if (element_model[Z] == m) elements.push_back(constants::atnr2letter(Z + 1));
		for (int a = 0; a < wavy_in.get_ncen(); a++) if (atom_model[a] == m) n++;
		std::cout << "Model " << sub_models[m]->get_salted_filename() << " predicts " << n << " atom(s): ";
		for (const auto& e : elements) std::cout << e << " ";
		std::cout << std::endl;
	}

	// Only an atom that no model handles goes to the spherical fill.
	std::vector<char> use_thakkar(wavy_in.get_ncen(), 0);
	int n_uncovered = 0;
	for (int a = 0; a < wavy_in.get_ncen(); a++)
		if (atom_model[a] < 0) { use_thakkar[a] = 1; ++n_uncovered; }
	if (n_uncovered > 0)
	{
		std::cout << "No model predicts " << n_uncovered
			<< " atom(s); they are filled with spherical Thakkar densities." << std::endl;
		spherical_fill_used = true;
		opt_in.needs_Thakkar_fill = true;
		opt_in.spherical_fill_charges.clear();
		estimate_fill_charges(wavy_in, use_thakkar, opt_in);
		for (int a = wavy_in.get_ncen() - 1; a >= 0; a--)
			if (use_thakkar[a])
			{
				wavy.erase_atom(a);
				atom_model.erase(atom_model.begin() + a);
				atom_in_model.erase(atom_in_model.begin() + a);
			}
	}

	natoms = wavy.get_ncen();
	if (wavy.get_nmo() != 0) wavy.clear_MOs();
	wavy.delete_basis_set();

	// Each element keeps the basis of the model that predicts it; anything else
	// would read that model's coefficients with the wrong radial functions.
	auto merged = std::make_shared<BasisSet>();
	std::string name;
	for (int Z = 0; Z < 118; Z++)
	{
		const int m = element_model[Z];
		if (m < 0) continue;
		const std::span<const SimplePrimitive> prims = (*sub_models[m]->get_model_basis())[Z];
		if (prims.empty()) continue;
		merged->set_count_for_element(Z, (int)prims.size());
		for (const auto& p : prims) merged->add_owned_primitive(p);
	}
	for (const auto& sub : sub_models)
		name += (name.empty() ? "" : "_plus_") + sub->get_dfbasis_name();
	merged->set_name(name);
	load_basis_into_WFN(wavy, merged, true);
	model_basis = merged;
	bbasis_set_loaded = true;
	config.dfbasis = name;
	for (const auto& sub : sub_models)
		config.salted_filename += (config.salted_filename.empty() ? "" : "+")
			+ sub->get_salted_filename().string();
}

vec SALTEDPredictor::merge_predictions()
{
	std::vector<vec> parts;
	std::vector<ivec> offsets;
	for (auto& sub : sub_models)
	{
		parts.push_back(sub->gen_SALTED_densities());
		offsets.push_back(coef_offsets(sub->wavy));
		err_checkf(parts.back().size() == (size_t)offsets.back().back(),
			"Model " + sub_models[parts.size() - 1]->get_salted_filename().string() + " predicted "
			+ std::to_string(parts.back().size()) + " coefficients for a basis of "
			+ std::to_string(offsets.back().back()) + " - model and basis do not belong together", std::cout);
	}

	vec coefs;
	for (int a = 0; a < wavy.get_ncen(); a++)
	{
		const int m = atom_model[a], p = atom_in_model[a];
		err_checkf(offsets[m][p + 1] <= (int)parts[m].size(),
			"Predicted coefficients are shorter than the basis of " + sub_models[m]->get_salted_filename().string(), std::cout);
		coefs.insert(coefs.end(), parts[m].begin() + offsets[m][p], parts[m].begin() + offsets[m][p + 1]);
	}
	const ivec merged_offsets = coef_offsets(wavy);
	err_checkf(coefs.size() == (size_t)merged_offsets.back(),
		"Stitched coefficients do not fit the combined basis", std::cout);

	// Each model only constrained its own share, which is meaningless on its own;
	// the electron count belongs to the whole density.
	bool wanted = force_charge_constraint;
	for (const auto& sub : sub_models) wanted = wanted || sub->wants_charge_constraint();
	if (wanted)
		apply_charge_constraint(wavy.get_atoms(), coefs, wavy.get_charge(),
								spherical_fill_used, n_filled,
								filled_eeq_charge, applied_fill_charge, std::cout);
	return coefs;
}

void calculateConjugate(SALTEDDescriptors& v2)
{
#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(v2.values().size()); ++i)
	{
		v2.values()[i] = std::conj(v2.values()[i]);
	}
}

void SALTEDPredictor::setup_atomic_environment()
{
	const std::shared_ptr<std::array<std::vector<primitive>, 118>> bs = wavy.get_basis_set_ptr();
	SALTED_Utils::set_lmax_nmax(lmax, nmax, *bs, config.species);

	atomic_symbols.reserve(wavy.get_ncen());
	for (int i = 0; i < wavy.get_ncen(); i++)
	{
		std::string label = wavy.get_atom_label(i);
		// Deuterium is hydrogen for the electron density; without this a joint
		// X-ray/neutron structure is refused with "Excluded species: D"
		if (label == "D" || label == "d")
		{
			label = "H";
		}
		atomic_symbols.emplace_back(label);
	}

	// Print all Atomic symbols
	if (debug)
	{
		std::cout << "Atomic symbols: ";
		for (const auto &symbol : atomic_symbols)
		{
			std::cout << symbol << " ";
		}
		std::cout << std::endl;
	}

	natoms = static_cast<int>(atomic_symbols.size());
	for (int i = 0; i < atomic_symbols.size(); i++)
	{
		atom_idx[atomic_symbols[i]].push_back(i);
		natom_dict[atomic_symbols[i]] += 1;
	}


	SALTED_Utils::FeatomicHyperParameters hp{
		.cutoff_radius = config.rcut1,
		.max_radial = config.nrad1,
		.max_angular = config.nang1,
		.atomic_gaussian_width = config.sig1,
		.center_atom_weight = 1.0,
		.species = config.species,
		.neighspe = config.neighspe1,
		.radial_basis = {
			.type = "Gto",
			.spline_accuracy = 1e-6
		},
		.cutoff_function = {
			.type = "ShiftedCosine",
			.width = 0.1
		}
	};

	featomic::SimpleSystem featomic_system = SALTED_Utils::gen_featomic_system(wavy);
	// RASCALINE (Generate descriptors)
	const auto _t_desc = std::chrono::steady_clock::now();
	v1 = SALTED_Utils::calculate_SALTED_descriptors(featomic_system, hp);

	if ((config.nrad2 != config.nrad1) || (config.nang2 != config.nang1) || (config.sig2 != config.sig1) || (config.rcut2 != config.rcut1) || (config.neighspe2 != config.neighspe1))
	{
		hp.max_radial = config.nrad2;
		hp.max_angular = config.nang2;
		hp.atomic_gaussian_width = config.sig2;
		hp.cutoff_radius = config.rcut2;
		hp.neighspe = config.neighspe2;
		v2 = SALTED_Utils::calculate_SALTED_descriptors(featomic_system, hp);
	}
	else
	{
		// Same hyperparameters: v2 would duplicate v1 exactly, so record that instead.
		// equicomb then reads conj(v1); conjugation only flips the sign of the imaginary
		// part, exact in IEEE, so the result is bit-identical
		v2_is_conj_of_v1 = true;
		v2.clear();
	}

	// Conjugate v2 once here rather than per use in equicomb; skipped when v2 is conj(v1)
	if (ProgressBar::report_counts)
		std::cout << "[stages] featomic descriptors "
				  << std::chrono::duration<double>(std::chrono::steady_clock::now() - _t_desc).count()
				  << " s" << std::endl;
	if (!v2_is_conj_of_v1)
		calculateConjugate(v2);
	std::cout << "Descriptor sets " << (v2_is_conj_of_v1 ? "identical: sharing one copy"
														: "differ: two copies held") << std::endl;
	// END RASCALINE
}


void SALTEDPredictor::read_model_data() {
	const auto _t_model = std::chrono::steady_clock::now();
	const std::filesystem::path _SALTEDpath = SALTED_DIR / config.salted_filename;
	// Kept open for the whole prediction; the matrices are fetched lambda by lambda
	model_file = std::make_unique<SALTED_BINARY_FILE>(_SALTEDpath);
	SALTED_BINARY_FILE &file = *model_file;
	if (config.field) {
		err_not_impl_f("Calculations using 'Field = True' are not yet supported", std::cout);
	}
	weights = file.read_weights();
	wigner3j = file.read_wigners();

	if (config.average) av_coefs = file.read_averages();
	if (config.sparsify) vfps = file.read_fps();


	// Only present species need their model data; absent ones are needed for their
	// shape alone, because `weights` is one flat vector laid out over every species
	// the model knows and the absent widths in front shift a present species' offset
	std::unordered_set<std::string> present;
	for (const std::string &spe : config.species)
		if (atom_idx.find(spe) != atom_idx.end()) present.insert(spe);

	// Indexing only, no payload: offset and shape per (species, lambda). The matrices
	// are read in load_model_lambda() and dropped again, each block used once per run
	feat_index = file.index_lambda_based_data("FEATS");
	proj_index = file.index_lambda_based_data("PROJ");
	model_species = present;

	// From the shapes alone: Mspe, the number of sparse environments of a present
	// species, and the projector width of every species the model knows
	for (const auto &[k, ref] : proj_index)
		proj_dims[k] = { ref.rows, ref.cols };
	for (const std::string &spe : present)
	{
		const auto it = feat_index.find(spe + "0");
		if (it != feat_index.end()) Mspe[spe] = static_cast<int>(it->second.rows);
	}

	if (ProgressBar::report_counts)
	{
		auto mb = [](const size_t doubles) { return doubles * sizeof(double) / 1048576.0; };
		size_t pes = 0, vm = 0, wg = 0;
		std::map<int, size_t> per_lam;
		for (const std::string &spe : present)
			for (int lam = 0; lam < lmax[spe] + 1; lam++)
			{
				const std::string k = spe + std::to_string(lam);
				const auto f = feat_index.find(k), pr = proj_index.find(k);
				const size_t p = (f == feat_index.end()) ? 0 : f->second.rows * f->second.cols;
				const size_t v = (pr == proj_index.end()) ? 0 : pr->second.rows * pr->second.cols;
				pes += p; vm += v; per_lam[lam] += p + v;
			}
		for (const auto &[lam, w] : wigner3j) wg += w.size();
		size_t worst = 0;
		for (const auto &[lam, sz] : per_lam) worst = std::max(worst, sz);
		std::cout << "[model] features " << mb(pes) << " MB + projectors " << mb(vm)
				  << " MB + wigner " << mb(wg) << " MB + weights " << mb(weights.size())
				  << " MB = " << mb(pes + vm + wg + weights.size()) << " MB if held whole;"
				  << " lazily, the largest lambda is " << mb(worst) << " MB" << std::endl;
		std::cout << "[model] indexed in "
				  << std::chrono::duration<double>(std::chrono::steady_clock::now() - _t_model).count()
				  << " s" << std::endl;
		std::cout << "[model] per lambda:";
		for (const auto &[lam, sz] : per_lam) std::cout << " l" << lam << "=" << mb(sz);
		std::cout << " MB" << std::endl;
	}
}


// Fetch the model matrices for one lambda, use them, drop them. Each lambda is
// visited once and nothing above it reads them again: the weight accounting
// downstream works off psi_nm and the projector shapes, which outlive the matrices
void SALTEDPredictor::load_model_lambda(const int lam)
{
	if (!model_file) return;
	for (const std::string &spe : model_species)
	{
		if (lam > lmax[spe]) continue;
		const std::string key = spe + std::to_string(lam);
		if (power_env_sparse.find(key) != power_env_sparse.end()) continue;
		const auto pr = proj_index.find(key);
		const auto ft = feat_index.find(key);
		if (pr == proj_index.end() || ft == feat_index.end()) continue;
		Vmat[key] = model_file->load_block(pr->second);
		dMatrix2 feats = model_file->load_block(ft->second);
		if (config.zeta == 1.0)
			power_env_sparse[key] = dot(Vmat[key], feats, true, false);
		else
			power_env_sparse[key] = std::move(feats);
	}
}

void SALTEDPredictor::free_model_lambda(const int lam)
{
	if (!model_file) return;
	for (const std::string &spe : model_species)
	{
		const std::string key = spe + std::to_string(lam);
		power_env_sparse.erase(key);
		Vmat.erase(key);
	}
}

vec SALTEDPredictor::predict()
{
	using namespace std;
#ifdef NOSPHERA2_USE_GPU
	struct salted_gpu_cache_scope {
		salted_gpu_cache_scope() { salted_gpu_clear_cache(); }
		~salted_gpu_cache_scope() { salted_gpu_clear_cache(); }
	} gpu_cache_scope;
#endif
	const auto _t_predict_start = std::chrono::steady_clock::now();
	auto _elapsed = [](const std::chrono::steady_clock::time_point &from)
	{ return std::chrono::duration<double>(std::chrono::steady_clock::now() - from).count(); };
	double _t_equicomb = 0.0, _t_kernels = 0.0, _t_model_wait = 0.0, _t_model_work = 0.0;
	bool overlap_model_loading = false;
#ifdef NOSPHERA2_USE_GPU
	overlap_model_loading = equicomb_gpu_enabled() && salted_gpu_available();
#endif
	// Compute equivariant descriptors for each lambda value entering the SPH expansion of the electron density
	// How many lambda blocks are alive at once. A block is natoms * (2*lam+1) *
	// featsize doubles and holding all of them sums to (nang+1)^2 times a single
	// block; fewer bounds that, at the cost of revisiting each species per group
	const int lmax_max = SALTED_Utils::get_lmax_max(lmax);
	ivec featsize(lmax_max + 1);
	std::vector<std::vector<dMatrix2>> psi_nm(config.species.size());
	for (int spe_idx = 0; spe_idx < (int)config.species.size(); spe_idx++)
		psi_nm[spe_idx].resize(lmax[config.species[spe_idx]] + 1);
	// The only quantity that crosses lambda: set at lam = 0, reused above it when zeta != 1
	std::vector<dMatrix2> kernell0(config.species.size());
	for (int lam = 0; lam <= lmax_max; lam++)
	{
		vec p;
		std::future<double> model_loader;
		if (overlap_model_loading)
			model_loader = std::async(std::launch::async, [this, lam]() {
				const auto start = std::chrono::steady_clock::now();
				load_model_lambda(lam);
				return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
			});
		const auto _t_eq = std::chrono::steady_clock::now();
		int llmax = 0;
		unordered_map<int, ivec> lvalues{};
		for (int l1 = 0; l1 < config.nang1 + 1; l1++)
		{
			for (int l2 = 0; l2 < config.nang2 + 1; l2++)
			{
				// keep only even combination to enforce inversion symmetry
				if ((lam + l1 + l2) % 2 == 0)
				{
					if (abs(l2 - lam) <= l1 && l1 <= (l2 + lam))
					{
						lvalues[llmax] = { l1, l2 };
						llmax += 1;
					}
				}
			}
		}
		// Fill dense array from dictionary
		ivec2 llvec(llmax, ivec(2));
		for (int i = 0; i < llmax; i++)
		{
			llvec[i] = lvalues[i];
		}

		cvec2 c2r = SALTED_Utils::complex_to_real_transformation({ 2 * lam + 1 })[0];

		featsize[lam] = config.nspe1 * config.nspe2 * config.nrad1 * config.nrad2 * llmax;
		ivec2 llvec_t = transpose<int>(llvec);
		if (config.sparsify)
		{
			int nfps = static_cast<int>(vfps[lam].size());
			p.assign((size_t)natoms * ((size_t)2 * lam + 1) * nfps, 0.0);
			equicomb(natoms, (config.nspe1 * config.nrad1), (config.nspe2 * config.nrad2), v1, v2, wigner3j[lam], llvec_t, lam, c2r, featsize[lam], nfps, vfps[lam], p, v2_is_conj_of_v1);
			featsize[lam] = nfps;
		}
		else
		{
			p.assign((size_t)natoms * ((size_t)2 * lam + 1) * featsize[lam], 0.0);
			equicomb(natoms, (config.nspe1 * config.nrad1), (config.nspe2 * config.nrad2), v1, v2, wigner3j[lam], llmax, llvec_t, lam, c2r, featsize[lam], p, v2_is_conj_of_v1);
		}
		_t_equicomb += _elapsed(_t_eq);
		const auto _t_wait = std::chrono::steady_clock::now();
		if (overlap_model_loading)
			_t_model_work += model_loader.get();
		else
			load_model_lambda(lam);
		const double model_wait = _elapsed(_t_wait);
		_t_model_wait += model_wait;
		if (!overlap_model_loading) _t_model_work += model_wait;
		const auto _t_kn = std::chrono::steady_clock::now();
		// Species-outer within the group, so a species keeps its sparse matrices hot
		for (int spe_idx = 0; spe_idx < (int)config.species.size(); spe_idx++)
		{
			const string spe = config.species[spe_idx];
			if (atom_idx.find(spe) == atom_idx.end()) continue;
			if (lam > lmax[spe]) continue;

			int lam2_1 = 2 * lam + 1;
			int row_size = featsize[lam] * lam2_1; // Size of a block of rows

			dMatrix2 pvec_lam(atom_idx[spe].size() * lam2_1, featsize[lam]);
			dMatrixRef2 _pvec(p.data(), natoms, featsize[lam] * lam2_1);
			double* pvec_ptr = pvec_lam.data();
			for (const int idx : atom_idx[spe])
			{
				auto _temp = Kokkos::submdspan(_pvec, idx, Kokkos::full_extent);
				std::copy(_temp.data_handle(), _temp.data_handle() + row_size, pvec_ptr);
				pvec_ptr += row_size;
			}
			//The regression GEMM stays on the CPU: the device barely wins, because fp64 runs
			//at a sixty-fourth rate on a consumer part and each call ships its own operands.
			//It would pay on a datacentre part.
			dMatrix2 kernel_nm = dot(pvec_lam, power_env_sparse[spe + to_string(lam)], false, true);

			if (config.zeta == 1)
			{
				psi_nm[spe_idx][lam] = kernel_nm;
			}
			else {

				if (lam == 0)
				{
					kernell0[spe_idx] = kernel_nm;
					kernel_nm = elementWiseExponentiation(kernel_nm, config.zeta);
				}
				else
				{
					for (size_t i1 = 0; i1 < natom_dict[spe]; ++i1)
					{
						for (size_t i2 = 0; i2 < Mspe[spe]; ++i2)
						{
							double scale_factor = pow(kernell0[spe_idx](i1, i2), config.zeta - 1);
							size_t base_i = i1 * lam2_1;
							size_t base_j = i2 * lam2_1;
							for (size_t i = 0; i < lam2_1; ++i)
							{
								for (size_t j = 0; j < lam2_1; ++j)
								{
									kernel_nm(base_i + i, base_j + j) *= scale_factor;
								}
							}
						}
					}
				}
				psi_nm[spe_idx][lam] = dot(kernel_nm, Vmat[spe + to_string(lam)], false, false);
			}
		}
		_t_kernels += _elapsed(_t_kn);
		free_model_lambda(lam);
	}

	unordered_map<string, dMatrix1> C{};
	unordered_map<string, int> ispe{};
	int isize = 0;
	for (int spe_idx = 0; spe_idx < config.species.size(); spe_idx++)
	{
		const string spe = config.species[spe_idx];
		if (atom_idx.find(spe) == atom_idx.end())
		{
			for (int l = 0; l < lmax[spe] + 1; ++l)
			{
				// Never loaded for an absent species; its shape was, and that is all this needs
				const auto dim_it = proj_dims.find(spe + to_string(l));
				if (dim_it == proj_dims.end() || dim_it->second[1] == 0)
				{
				   std::cout << "The projector for species " << spe << " and l = " << l << " does not exist. This is a problem with the model, not NoSpherA2.\n";
				   std::cout << "Continuing with the next species..., make sure there is no: " << spe << " in the structure you are trying to predict!!!!\n";
					break;
				}

				// for (int n = 0; n < nmax[spe + to_string(l)]; ++n)
				//{
				//     isize += static_cast<int>(Vmat[spe + to_string(l)][0].size());
				// }
				isize += static_cast<int>(dim_it->second[1]) * nmax[spe + to_string(l)];
			}
			continue;
		}
		ispe[spe] = 0;
		for (int l = 0; l < lmax[spe] + 1; ++l)
		{
			for (int n = 0; n < nmax[spe + to_string(l)]; ++n)
			{
				// int Mcut = static_cast<int>(psi_nm[spe + to_string(l)][0].size());
				int Mcut = static_cast<int>(psi_nm[spe_idx][l].extent(1));
				// Check if isize + Mcut > weights.size()
				err_checkf(isize + Mcut <= weights.size(), "isize + Mcut > weights.size()", std::cout);

				dMatrix1 weights_subset(Mcut);
				std::copy(weights.data() + isize, weights.data() + isize + Mcut, weights_subset.data());

				C[spe + to_string(l) + to_string(n)] = dot(psi_nm[spe_idx][l], weights_subset, false);

				isize += Mcut;
			}
		}
	}
	psi_nm.clear();
	psi_nm.shrink_to_fit();


	int Tsize = 0;
	for (int iat = 0; iat < natoms; iat++)
	{
		string spe = atomic_symbols[iat];
		for (int l = 0; l < lmax[spe] + 1; l++)
		{
			for (int n = 0; n < nmax[spe + to_string(l)]; n++)
			{
				Tsize += 2 * l + 1;
			}
		}
	}

	// A model carrying a BASIS block that does not describe what it predicts is
	// broken; name the species that disagrees instead of failing later on a
	// total count that says nothing about where it went wrong.
	if (bbasis_set_loaded)
	{
		const aux_density_table t(*wavy.get_atoms_ptr());
		if (Tsize != t.n_coef)
		{
			std::cout << "The model predicts " << Tsize << " coefficients, the basis in the model file holds "
					  << t.n_coef << ":" << std::endl;
			std::unordered_set<std::string> reported;
			for (int iat = 0; iat < natoms; iat++)
			{
				const std::string& spe = atomic_symbols[iat];
				if (!reported.insert(spe).second) continue;
				int model_size = 0;
				for (int l = 0; l < lmax[spe] + 1; l++)
					model_size += nmax[spe + to_string(l)] * (2 * l + 1);
				const int basis_size = (iat + 1 < natoms ? t.coef_off[t.sh_start[iat + 1]] : t.n_coef)
									 - t.coef_off[t.sh_start[iat]];
				std::cout << "   " << spe << ": predicted " << model_size << ", basis " << basis_size
						  << (model_size == basis_size ? "" : "   <-- mismatch") << std::endl;
				if (model_size == basis_size) continue;
				// Say which angular momentum is short, which is what has to be fixed
				// in the model file - a total count alone does not point anywhere.
				ivec basis_l;
				const int sh1 = (iat + 1 < natoms ? t.sh_start[iat + 1] : t.n_sh);
				for (int s = t.sh_start[iat]; s < sh1; s++)
				{
					if ((int)basis_l.size() <= t.sh_l[s]) basis_l.resize(t.sh_l[s] + 1, 0);
					basis_l[t.sh_l[s]]++;
				}
				for (int l = 0; l <= std::max(lmax[spe], (int)basis_l.size() - 1); l++)
					std::cout << "      l = " << l << ": model "
							  << (l <= lmax[spe] ? nmax[spe + to_string(l)] : 0) << " shells, basis "
							  << (l < (int)basis_l.size() ? basis_l[l] : 0) << " shells" << std::endl;
			}
			err_checkf(false, "SALTED model and the basis set in its file do not belong together", std::cout);
		}
	}

	vec Av_coeffs(Tsize, 0.0);

	// fill vector of predictions
	int i = 0;
	vec pred_coefs(Tsize, 0.0);
	for (int iat = 0; iat < natoms; ++iat)
	{
		string spe = atomic_symbols[iat];
		for (int l = 0; l < lmax[spe] + 1; ++l)
		{
			for (int n = 0; n < nmax[spe + to_string(l)]; ++n)
			{
				// for (int ind = 0; ind < 2 * l + 1; ++ind)
				//{
				//     pred_coefs[i + ind] = C[spe + to_string(l) + to_string(n)][ispe[spe] * (2 * l + 1) + ind];
				// }
				std::copy_n(C[spe + to_string(l) + to_string(n)].data() + ispe[spe] * (2 * l + 1), 2 * l + 1, pred_coefs.begin() + i);

				if (config.average && l == 0)
				{
					Av_coeffs[i] = av_coefs[spe][n];
				}
				i += 2 * l + 1;
			}
		}
		ispe[spe] += 1;
	}

	if (config.average)
	{
		for (i = 0; i < Tsize; i++)
		{
			pred_coefs[i] += Av_coeffs[i];
		}
	}

	//std::cout << "          ... done!\nNumber of predicted coefficients: " << pred_coefs.size() << endl;
	// npy::npy_data<double> coeffs;
	// coeffs.data = pred_coefs;
	// coeffs.fortran_order = false;
	// coeffs.shape = { unsigned long(pred_coefs.size()) };
	// npy::write_npy("folder_model.npy", coeffs);
	if (ProgressBar::report_counts)
	{
		const double total = std::chrono::duration<double>(
			std::chrono::steady_clock::now() - _t_predict_start).count();
		std::cout << "[stages] predict " << total << " s = equicomb " << _t_equicomb
				  << " s (" << (100.0 * _t_equicomb / total) << "%) + kernels "
				  << _t_kernels << " s (" << (100.0 * _t_kernels / total)
				  << "%) + model wait " << _t_model_wait << " s ("
				  << (100.0 * _t_model_wait / total) << "%) + rest "
				  << (total - _t_equicomb - _t_kernels - _t_model_wait) << " s; model preparation "
				  << _t_model_work << " s" << std::endl;
	}
	return pred_coefs;
}

bool SALTEDPredictor::wants_charge_constraint() const
{
	if (force_charge_constraint) return true;
	if (!model_file || !model_file->charge_constraint_defined()) return false;
	const auto entries = model_file->read_charge_constraint();
	const auto mode_it = entries.find("MODE");
	const int mode = (mode_it != entries.end() && !mode_it->second.empty())
						 ? static_cast<int>(std::lround(mode_it->second[0]))
						 : 0;
	if (mode == 1) return true;
	if (mode != 0)
		std::cout << "Unknown charge-constraint mode " << mode
				  << " in the model file; leaving the density alone." << std::endl;
	return false;
}

vec SALTEDPredictor::gen_SALTED_densities()
{
	using namespace std;
	if (!sub_models.empty())
		return merge_predictions();
	if (coef_file != "")
	{
		vec coefs{};
		std::cout << "Reading coefficients from file: " << coef_file << endl;
		read_npy<double>(coef_file, coefs);
		vec double_coefs(coefs.size());
		for (int i = 0; i < coefs.size(); i++)
		{
			double_coefs[i] = static_cast<double>(coefs[i]);
		}
		return double_coefs;
	}

	// Run generation of tsc file
	_time_point start;
	if (debug)
		start = get_time();

	setup_atomic_environment();


	read_model_data();


	vec coefs = predict();

	// File VERSION 3 models carry an optional NORMC block asking for the
	// electron count to be constrained. Applied here rather than at each call
	// site so the tsc, the charge table and the cubes all see the same density.
	// V2 models have no such block, so they are untouched. A sub model of a
	// merged prediction only holds a share of the electrons, so it is skipped
	// and the stitched density is constrained once instead.
	if (!skip_charge_constraint && wants_charge_constraint())
		apply_charge_constraint(wavy.get_atoms(), coefs, wavy.get_charge(),
								spherical_fill_used, n_filled,
								filled_eeq_charge, applied_fill_charge, std::cout);

	shrink_intermediate_vectors();
	return coefs;
}

void SALTEDPredictor::shrink_intermediate_vectors()
{
	v1.clear();
	v2.clear();
	weights.clear();
	Vmat.clear();
	natom_dict.clear();
	lmax.clear();
	nmax.clear();
	Mspe.clear();
	vfps.clear();
	wigner3j.clear();
	av_coefs.clear();
	power_env_sparse.clear();
	featsize.clear();
	v1.shrink_to_fit();
	v2.shrink_to_fit();
	weights.shrink_to_fit();
	std::unordered_map<std::string, dMatrix2> umap;
	std::unordered_map<std::string, int> umap2;
	std::unordered_map<std::string, dMatrix1> umap3;
	Vmat.swap(umap);
	natom_dict.swap(umap2);
	lmax.swap(umap2);
	nmax.swap(umap2);
	power_env_sparse.swap(umap);
};
