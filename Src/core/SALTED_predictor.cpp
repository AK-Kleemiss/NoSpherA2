#include "pch.h"
#include "SALTED_predictor.h"
#ifdef NOSPHERA2_USE_GPU
#include "salted_gpu.h"
#include "sf_gpu.h"
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
#include "citations.h"
#include <cstring>
#include <filesystem>
#include <future>


SALTEDPredictor::SALTEDPredictor(WFN wavy_in, options& opt_in)
{
	std::filesystem::path _path = opt_in.salted_model_dir;
	SALTED_DIR = opt_in.salted_model_dir;
	debug = opt_in.debug;

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

	// A model may be named by its own file: two models in one folder are otherwise indistinguishable.
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

	featomic_system = SALTED_Utils::gen_featomic_system(wavy);

	natoms = 0;
	for (auto a : wavy.get_atoms()) {
		const std::string atom_symbol = a.get_label();
		if (std::find(config.species.begin(), config.species.end(), atom_symbol) != config.species.end())
			natoms++;
	}//natoms is the number of atoms, that featomic creates a descriptor as a central atom for.

	const auto _t_filter = std::chrono::steady_clock::now();
	const std::vector<char> use_thakkar = SALTED_Utils::filter_input(wavy, opt_in, config);
	if (ProgressBar::report_counts)
		std::cout << "[stages] filter input "
				  << std::chrono::duration<double>(std::chrono::steady_clock::now() - _t_filter).count()
				  << " s" << std::endl;
	if (!use_thakkar.empty()) {
		spherical_fill_used = true;
		estimate_fill_charges(wavy_in, use_thakkar, opt_in);
	}

	//wavy.write_xyz("temp_rascaline.xyz"); //Also this
	//config.predict_filename = "temp_rascaline.xyz";


	if (wavy.get_nmo() != 0)
		wavy.clear_MOs(); // Delete unneccesarry MOs, since we are predicting anyway.

	wavy.delete_basis_set();
	if (file.basis_set_defined()) {
		model_basis = file.read_basis_set();
		load_basis_into_WFN(wavy, model_basis);
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

// The predictor's filter only erases atoms and keeps the rest in order and untouched, so one forward
// walk with an exact comparison maps its structure back onto the full one.
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

// Each atom's first coefficient in the auxiliary basis, as the density evaluation indexes it; one entry past the end.
static ivec coef_offsets(const WFN& w)
{
	const aux_density_table t(*w.get_atoms_ptr());
	ivec off(w.get_ncen() + 1, t.n_coef);
	for (int a = 0; a < w.get_ncen(); a++) off[a] = t.coef_off[t.sh_start[a]];
	return off;
}

void SALTEDPredictor::build_merged(const WFN& wavy_in, options& opt_in)
{
	// One prediction per model on the whole structure: an atom another model owns is still a neighbour.
	// Each per-atom block carries its own model's auxiliary basis, so stitching them is a concatenation.
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
		if (!sub->basis_set_loaded())
			load_basis_into_WFN(sub->wavy, sub->get_model_basis());
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
		// A model that dropped this atom for having no environment cannot predict it.
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
	load_basis_into_WFN(wavy, merged);
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

	return coefs;
}

// predict() only wants C_n = psi_nm w_n = (K V) w_n; PROJW stores V W^T per (species, l), so the
// kernel product makes nmax columns instead of Mcut (V6 carbon l = 5: 1936) and C_n is column n.
// PROJW has the PROJ layout: datatype word, species count, then per species its tag, lambda count
// and one 2D dataset per lambda, Mspe (2l + 1) x nmax instead of x Mcut. V6: 782 -> 525 MB, PROJ
// 258 MB -> PROJW 0.5 MB. The weights are walked as read_model_data() walks them
void fold_salted_file(const std::filesystem::path& in, const std::filesystem::path& out)
{
	// One BLAS thread: a threaded GEMM sums in an order set by the thread count, and the file must
	// not depend on the machine that wrote it (Pi 4 OpenBLAS: 1.6e-14 rel. apart). The run ends here
	MKL_Set_Num_Threads(1);
	SALTED_BINARY_FILE file(in);
	err_checkf(!file.has_block("PROJW"), in.string() + " is folded already", std::cout);
	SALTEDConfig config;
	file.populate_config(config);
	const std::shared_ptr<BasisSet> basis = file.basis_set_defined() ? file.read_basis_set() : BasisSetLibrary::get_basis_set(config.dfbasis);
	std::unordered_map<std::string, int> lmax, nmax;
	SALTED_Utils::set_lmax_nmax(lmax, nmax, *basis, config.species);
	const vec weights = file.read_weights();
	const auto proj = file.index_lambda_based_data("PROJ");

	std::string block;
	auto put = [&block](const auto v) { block.append(reinterpret_cast<const char*>(&v), sizeof(v)); };
	put(int32_t(3));   // the datatype word SALTED writes for these blocks; nothing reads it
	put(int32_t(0));   // species count, set below
	int32_t n_species = 0;
	size_t isize = 0;
	for (const std::string& spe : config.species)
	{
		int nlam = 0;
		while (nlam < lmax.at(spe) + 1)
		{
			const auto it = proj.find(spe + std::to_string(nlam));
			if (it == proj.end() || it->second.cols == 0) break;
			nlam++;
		}
		if (nlam == 0) continue;
		n_species++;
		std::string tag = spe;
		tag.resize(5, '\0');
		block += tag;
		put(int32_t(nlam));
		for (int l = 0; l < nlam; l++)
		{
			const std::string key = spe + std::to_string(l);
			const dMatrix2 V = file.load_block(proj.at(key));
			err_checkf(isize + nmax.at(key) * V.extent(1) <= weights.size(), "The weights end before the projectors of " + key, std::cout);
			// nmax rows of V's width, n-major as the flat weights lie
			dMatrix2 Wt(nmax.at(key), V.extent(1));
			std::copy_n(weights.data() + isize, Wt.extent(0) * Wt.extent(1), Wt.data());
			const dMatrix2 VW = dot(V, Wt, false, true);
			isize += Wt.extent(0) * Wt.extent(1);
			put(int32_t(2));
			put(uint32_t(VW.extent(0)));
			put(uint32_t(VW.extent(1)));
			block.append(reinterpret_cast<const char*>(VW.data()), VW.extent(0) * VW.extent(1) * sizeof(double));
		}
	}
	std::memcpy(block.data() + sizeof(int32_t), &n_species, sizeof(int32_t));
	file.write_with_block_replaced(out, 4, "PROJ", "PROJW", block);
	std::cout << "Wrote " << out.string() << ": " << n_species << " species, projectors with the weights folded in" << std::endl;
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
	SALTED_Utils::set_lmax_nmax(lmax, nmax, *get_model_basis(), config.species);

	atomic_symbols.reserve(natoms);
	for (int i = 0; i < featomic_system.size(); i++)
	{
		std::string label = constants::Labels[featomic_system.types()[i]];
		if (std::find(config.species.begin(), config.species.end(), label) == config.species.end()) continue;
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

	for (int i = 0; i < natoms; i++)
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
	projector_folded = file.has_block("PROJW");
	if (!projector_folded) weights = file.read_weights();
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
	// are read in read_model_lambda() and dropped again, each block used once per run
	feat_index = file.index_lambda_based_data("FEATS");
	proj_index = file.index_lambda_based_data(projector_folded ? "PROJW" : "PROJ");
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
	// Where each present (species, l) starts in the flat weights, which run species, l, n over
	// every species the model knows; predict() contracts psi_nm with that block
	if (!projector_folded)
	{
		size_t isize = 0;
		for (const std::string &spe : config.species)
			for (int l = 0; l < lmax[spe] + 1; ++l)
			{
				const std::string k = spe + std::to_string(l);
				const auto dim_it = proj_dims.find(k);
				if (dim_it == proj_dims.end() || dim_it->second[1] == 0)
				{
					if (!present.count(spe))
					{
						std::cout << "The projector for species " << spe << " and l = " << l << " does not exist. This is a problem with the model, not NoSpherA2.\n";
						std::cout << "Continuing with the next species..., make sure there is no: " << spe << " in the structure you are trying to predict!!!!\n";
					}
					break;
				}
				if (present.count(spe)) weight_offset[k] = isize;
				isize += dim_it->second[1] * nmax[k];
			}
		err_checkf(isize <= weights.size(), "isize + Mcut > weights.size()", std::cout);
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


// Model matrices are loaded per lambda and dropped after use: nothing rereads them, the weight
// offsets come from the projector shapes
// File reads only, touching no member the kernels use, so it can run on another thread
SALTEDPredictor::lambda_blocks SALTEDPredictor::read_model_lambda(const int lam)
{
	lambda_blocks out;
	if (!model_file) return out;
	std::vector<std::string> keys;
	std::vector<SALTED_BINARY_FILE::block_ref> refs;   // projector, features per key
	for (const std::string &spe : model_species)
	{
		if (lam > lmax.at(spe)) continue;
		const std::string key = spe + std::to_string(lam);
		const auto pr = proj_index.find(key);
		const auto ft = feat_index.find(key);
		if (pr == proj_index.end() || ft == feat_index.end()) continue;
		keys.push_back(key);
		refs.push_back(pr->second);
		refs.push_back(ft->second);
	}
	std::vector<dMatrix2> blocks = model_file->load_blocks(refs);
	for (std::size_t k = 0; k < keys.size(); k++)
		out.emplace_back(keys[k], std::move(blocks[2 * k]), std::move(blocks[2 * k + 1]));
	return out;
}

void SALTEDPredictor::install_model_lambda(lambda_blocks blocks)
{
	for (auto &[key, V, feats] : blocks)
	{
		if (power_env_sparse.find(key) != power_env_sparse.end()) continue;
		Vmat[key] = std::move(V);
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
	bool gpu_equicomb = false;
#ifdef NOSPHERA2_USE_GPU
	gpu_equicomb = equicomb_gpu_enabled() && salted_gpu_available();
	//Context creation started with the run; the first kernel would otherwise pay for it
	if (gpu_equicomb)
		sf_gpu_warmup_wait();
#endif
	const int lmax_max = SALTED_Utils::get_lmax_max(lmax);
	// Reading a lambda's blocks takes about as long as its kernels (V6 sucrose: ~180 ms from page
	// cache on x64, 65-140 MB/s from Android flash), so read lambda + 1 while lambda's kernels run.
	// Plain file reads, so no OpenMP team of its own (see below). One thread on top of -cpus;
	// holds two lambdas at once, ~100 MB more at the peak for the V6 model. Output is identical
	std::future<lambda_blocks> next_lambda = std::async(std::launch::async, &SALTEDPredictor::read_model_lambda, this, 0);
	// The (l1, l2) pairs and complex-to-real matrix of every lambda, which the norm wants up front
	std::vector<ivec2> llvec_all(lmax_max + 1);
	std::vector<cvec2> c2r_all(lmax_max + 1);
	for (int lam = 0; lam <= lmax_max; lam++)
	{
		ivec2 llvec;
		for (int l1 = 0; l1 < config.nang1 + 1; l1++)
			for (int l2 = 0; l2 < config.nang2 + 1; l2++)
				// keep only even combination to enforce inversion symmetry
				if ((lam + l1 + l2) % 2 == 0 && abs(l2 - lam) <= l1 && l1 <= (l2 + lam))
					llvec.push_back({ l1, l2 });
		llvec_all[lam] = transpose<int>(llvec);
		c2r_all[lam] = SALTED_Utils::complex_to_real_transformation({ 2 * lam + 1 })[0];
	}
	// Every lambda's norm in one pass over the atoms, with no density matrices kept for all of
	// them; the device computes its own
	vec2 norms;
	if (config.sparsify && !gpu_equicomb)
	{
		const auto _t_norm = std::chrono::steady_clock::now();
		std::vector<const vec *> w3j_all(lmax_max + 1);
		for (int lam = 0; lam <= lmax_max; lam++)
			w3j_all[lam] = &wigner3j[lam];
		norms = equicomb_norms(natoms, config.nspe1 * config.nrad1, config.nspe2 * config.nrad2, v1, v2, w3j_all, llvec_all, c2r_all, v2_is_conj_of_v1);
		const double norm_s = _elapsed(_t_norm);
		_t_equicomb += norm_s;
		throughput::record_time("SALTED equicomb", false, 1000.0 * norm_s);
		if (ProgressBar::report_counts)
			std::cout << "[equicomb] norms of lambda 0.." << lmax_max << ": " << 1000.0 * norm_s << " ms" << std::endl;
	}
	// Compute equivariant descriptors for each lambda value entering the SPH expansion of the electron density
	ivec featsize(lmax_max + 1);
	std::vector<std::vector<dMatrix2>> psi_nm(config.species.size());
	for (int spe_idx = 0; spe_idx < (int)config.species.size(); spe_idx++)
		psi_nm[spe_idx].resize(lmax[config.species[spe_idx]] + 1);
	// The only quantity that crosses lambda: set at lam = 0, reused above it when zeta != 1
	std::vector<dMatrix2> kernell0(config.species.size());
	for (int lam = 0; lam <= lmax_max; lam++)
	{
		vec p;
		const auto _t_eq = std::chrono::steady_clock::now();
		ivec2 &llvec_t = llvec_all[lam];
		cvec2 &c2r = c2r_all[lam];
		const int llmax = static_cast<int>(llvec_t[0].size());
		featsize[lam] = config.nspe1 * config.nspe2 * config.nrad1 * config.nrad2 * llmax;
		if (config.sparsify)
		{
			int nfps = static_cast<int>(vfps[lam].size());
			p.assign((size_t)natoms * ((size_t)2 * lam + 1) * nfps, 0.0);
			equicomb(natoms, (config.nspe1 * config.nrad1), (config.nspe2 * config.nrad2), v1, v2, wigner3j[lam], llvec_t, lam, c2r, featsize[lam], nfps, vfps[lam], p, v2_is_conj_of_v1,
				norms.empty() ? nullptr : norms[lam].data());
			featsize[lam] = nfps;
		}
		else
		{
			p.assign((size_t)natoms * ((size_t)2 * lam + 1) * featsize[lam], 0.0);
			equicomb(natoms, (config.nspe1 * config.nrad1), (config.nspe2 * config.nrad2), v1, v2, wigner3j[lam], llmax, llvec_t, lam, c2r, featsize[lam], p, v2_is_conj_of_v1);
		}
		_t_equicomb += _elapsed(_t_eq);
		// The file reads of lambda ran during lambda - 1's kernels; installing the blocks (the
		// projector products) stays in line, also with the descriptors on the device: doing that
		// on a second thread gave that thread its own OpenMP/MKL team, whose spin-wait then cost
		// the kernels below more than the overlap could ever hide
		const auto _t_wait = std::chrono::steady_clock::now();
		install_model_lambda(next_lambda.get());
		if (lam < lmax_max)
			next_lambda = std::async(std::launch::async, &SALTEDPredictor::read_model_lambda, this, lam + 1);
		const double model_wait = _elapsed(_t_wait);
		_t_model_wait += model_wait;
		_t_model_work += model_wait;
		const auto _t_kn = std::chrono::steady_clock::now();
#if defined(NSA2_OPENBLAS) && !defined(NSA2_ARMPL) && !defined(__ANDROID__)
		// pch.h keeps a pthreads OpenBLAS serial, as it cannot see an enclosing OpenMP region
		// (Android's OpenMP OpenBLAS keeps the count pch.h caps to cpu0's cluster, so not there).
		// This loop is outside one (its parallel for calls no BLAS), so its GEMMs may take every core,
		// but one while lambda + 1 is read: OpenBLAS splits a GEMM evenly, so a thread the reader
		// pushes off its core holds back all of them (Pi 4, -cpus 4: kernels 1.9 -> 3.0 s)
		const int blas_n = lam < lmax_max ? std::max(1, std::min(omp_get_max_threads(), omp_get_num_procs() - 1))
		                                  : omp_get_max_threads();
		struct blas_threads_scope {
			explicit blas_threads_scope(const int n) { openblas_set_num_threads(n); }
			~blas_threads_scope() { openblas_set_num_threads(1); }
		} blas_threads(blas_n);
#endif
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
				psi_nm[spe_idx][lam] = std::move(kernel_nm);
			}
			else {

				if (lam == 0)
				{
					// Every lambda above scales by k0^(zeta-1), so that power is taken once, here
					kernell0[spe_idx] = elementWiseExponentiation(kernel_nm, config.zeta - 1);
					kernel_nm = elementWiseExponentiation(kernel_nm, config.zeta);
				}
				else
				{
					const dMatrix2 &k0 = kernell0[spe_idx];
					const int n_at = (int)natom_dict[spe];
					const size_t n_env = Mspe[spe];
#pragma omp parallel for
					for (int i1 = 0; i1 < n_at; ++i1)
					{
						for (size_t i2 = 0; i2 < n_env; ++i2)
						{
							const double scale_factor = k0(i1, i2);
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
	for (int spe_idx = 0; spe_idx < config.species.size(); spe_idx++)
	{
		const string spe = config.species[spe_idx];
		if (atom_idx.find(spe) == atom_idx.end()) continue;
		ispe[spe] = 0;
		for (int l = 0; l < lmax[spe] + 1; ++l)
		{
			const string key = spe + to_string(l);
			const dMatrix2 &psi = psi_nm[spe_idx][l];
			const size_t Mcut = psi.extent(1);
			for (int n = 0; n < nmax[key]; ++n)
			{
				if (projector_folded)
				{
					// The weights are in the PROJW projector (fold_salted_file): column n is C_n
					dMatrix1 c(psi.extent(0));
					for (size_t i = 0; i < psi.extent(0); i++) c(i) = psi(i, n);
					C[key + to_string(n)] = std::move(c);
					continue;
				}
				dMatrix1 weights_subset(Mcut);
				const double *w = weights.data() + weight_offset.at(key) + n * Mcut;
				std::copy(w, w + Mcut, weights_subset.data());
				C[key + to_string(n)] = dot(psi, weights_subset, false);
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

	// A BASIS block that does not describe what the model predicts is broken; name the disagreeing species
	// rather than fail later on a total count.
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
				// The short angular momentum is what has to be fixed in the model file.
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

vec SALTEDPredictor::gen_SALTED_densities()
{
	using namespace std;
	citations::cite(citations::Method::SALTED, std::cout);
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
