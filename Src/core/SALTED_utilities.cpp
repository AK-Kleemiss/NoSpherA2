#include "pch.h"
#include "tuning.h"
#include "SALTED_utilities.h"
#include "basis_set.h"
#include <set>
#include "integration_params.h"
#include "SALTED_io.h"
#include "constants.h"
#include "atoms.h"
#include "cube.h"
#include "wfn_class.h"
#include "metatensor.hpp"
#include "featomic.hpp"

cvec3 SALTED_Utils::complex_to_real_transformation(ivec sizes)
{
	const double sqrt_2 = sqrt(2.0);
	using namespace std;
	vector<cvec2> matrices{};
	for (int i = 0; i < sizes.size(); i++)
	{
		int lval = (sizes[i] - 1) / 2;
		int st = (lval & 1) ? 1 : -1;

		cvec2 transformed_matrix(sizes[i], cvec(sizes[i], 0.0));
		for (int j = 0; j < lval; j++)
		{
			transformed_matrix[j][j] = complex<double>(0.0, 1.0);
			transformed_matrix[j][sizes[i] - j - 1] = complex<double>(0.0, st);
			transformed_matrix[sizes[i] - j - 1][j] = complex<double>(1.0, 0.0);
			transformed_matrix[sizes[i] - j - 1][sizes[i] - j - 1] = complex<double>(-st, 0.0);
			st = -st;
		}
		transformed_matrix[lval][lval] = sqrt_2;
		// Divide each element by sqrt(2.0)
		for (auto& row : transformed_matrix)
		{
			for (auto& elem : row)
			{
				elem /= sqrt_2;
			}
		}
		matrices.push_back(transformed_matrix);
	}
	return matrices;
}

int SALTED_Utils::get_lmax_max(std::unordered_map<std::string, int>& lmax)
{
	int lmax_max = 0;
	for (auto& [key, value] : lmax)
	{
		if (value > lmax_max)
		{
			lmax_max = value;
		}
	}
	return lmax_max;
}

void SALTED_Utils::set_lmax_nmax(std::unordered_map<std::string, int>& lmax, std::unordered_map<std::string, int>& nmax, const BasisSet& basis_set, std::vector<std::string> species)
{
	// lmax = {"C": 5, "H":2,...} with the numbers beeing the maximum angular momentum (type) for the given atom
	// nmax = {C0: 10, C1: 7, ...} with the numbers beeing the maximum number of primitives for the given atom and type

	for (const auto& spe : species)
	{
		const int Z = constants::get_Z_from_label(spe.c_str());
		const auto& primitives = basis_set[Z];

		// shell_ids[l] contains all unique shells with angular momentum l.
		std::unordered_map<int, std::unordered_set<int>> shell_ids;

		for (const auto& primitive : primitives)
			shell_ids[primitive.type].insert(primitive.shell);

		const auto max_l = std::ranges::max(
			primitives | std::views::transform(&SimplePrimitive::type)
		);

		lmax[spe] = max_l;

		for (int l = 0; l <= max_l; ++l)
			nmax[spe + std::to_string(l)] = static_cast<int>(shell_ids[l].size());
	}
}

std::vector<char> SALTED_Utils::filter_input(WFN& wavy, options& opt, const SALTEDConfig& config) {
	// Two kinds of atom cannot be predicted, and both are handed to the spherical
	// Thakkar fill instead of guessed at:
	//
	//  1. a species the model was never trained on;
	//  2. an atom with nothing inside the descriptor cutoff.
	//
	// Removing them disturbs nobody: the test is symmetric, so an atom with
	// nothing within rcut is also nobody's neighbour within rcut, and every other
	// atom's environment is exactly what it was.
	//
	// The test here is purely geometric, while featomic additionally ignores
	// neighbours whose species is outside neighspe. An atom with neighbours but no
	// ALLOWED ones therefore still reaches equicomb - the zero guard there catches
	// it, leaves it spherical rather than NaN, and says so in the log.

	auto is_species_of_atom_defined_in_model = [&config](const std::string& atom_symbol) {
		return std::find(config.species.begin(), config.species.end(), atom_symbol) == config.species.end(); };

	const int ncen_in = wavy.get_ncen();
	std::vector<char> use_thakkar(ncen_in, 0);
	for (int a = 0; a < ncen_in; a++)
	{
		if (is_species_of_atom_defined_in_model(constants::atnr2letter(wavy.get_atom_charge(a))))
			use_thakkar[a] = 1;
	}
	const int n_unknown = static_cast<int>(std::count(use_thakkar.begin(), use_thakkar.end(), (char)1));

	const double rcut = std::min(config.rcut1, config.rcut2);
	// rcut is in Angstrom, the coordinates may not be
	const double rcut_internal = wavy.get_isBohr() ? constants::ang2bohr(rcut) : rcut;
	const double cut_sq = rcut_internal * rcut_internal;
	// Cells of edge rcut: a neighbour within rcut sits in the atom's own cell or one of the 26
	// around it, so the test grows with the atom count rather than its square. No periodic
	// images, because the featomic system has no cell either.
	std::vector<std::array<double, 3>> pos(ncen_in);
	std::array<double, 3> lo;
	lo.fill(std::numeric_limits<double>::max());
	for (int a = 0; a < ncen_in; a++)
		for (int ax = 0; ax < 3; ax++)
		{
			pos[a][ax] = wavy.get_atom_coordinate(a, ax);
			lo[ax] = std::min(lo[ax], pos[a][ax]);
		}
	auto cell = [&](const int a, const int ax) { return static_cast<long long>((pos[a][ax] - lo[ax]) / rcut_internal); };
	long long n_cell[3] = { 1, 1, 1 };
	for (int a = 0; a < ncen_in; a++)
		for (int ax = 0; ax < 3; ax++)
			n_cell[ax] = std::max(n_cell[ax], cell(a, ax) + 1);
	auto key = [&](const long long x, const long long y, const long long z) { return x + n_cell[0] * (y + n_cell[1] * z); };
	std::unordered_map<long long, ivec> cells;
	for (int a = 0; a < ncen_in; a++)
		cells[key(cell(a, 0), cell(a, 1), cell(a, 2))].push_back(a);

	int n_isolated = 0;
#pragma omp parallel for reduction(+ : n_isolated)
	for (int a = 0; a < ncen_in; a++)
	{
		if (use_thakkar[a]) continue;
		bool lonely = true;
		const long long cx = cell(a, 0), cy = cell(a, 1), cz = cell(a, 2);
		for (long long x = std::max(cx - 1, 0LL); x <= std::min(cx + 1, n_cell[0] - 1) && lonely; x++)
			for (long long y = std::max(cy - 1, 0LL); y <= std::min(cy + 1, n_cell[1] - 1) && lonely; y++)
				for (long long z = std::max(cz - 1, 0LL); z <= std::min(cz + 1, n_cell[2] - 1) && lonely; z++)
				{
					const auto it = cells.find(key(x, y, z));
					if (it == cells.end()) continue;
					for (const int b : it->second)
					{
						if (b == a) continue;
						double d_sq = 0.0;
						for (int ax = 0; ax < 3; ax++)
						{
							const double dx = pos[a][ax] - pos[b][ax];
							d_sq += dx * dx;
						}
						if (d_sq < cut_sq) { lonely = false; break; }
					}
				}
		if (lonely)
		{
			use_thakkar[a] = 1;   // distinct indices, and char so there is no bitfield to race on
			++n_isolated;
		}
	}

	if (n_unknown + n_isolated > 0)
	{
		if (n_unknown > 0)
		{
			std::cout << "WARNING: Not all species in the structure are known to the model. The following species are not known: ";
			for (int a = 0; a < ncen_in; a++)
			{
				if (is_species_of_atom_defined_in_model(constants::atnr2letter(wavy.get_atom_charge(a))))
				{
					std::cout << constants::atnr2letter(wavy.get_atom_charge(a)) << " ";
				}
			}
			std::cout << std::endl;
		}
		if (n_isolated > 0)
		{
			std::cout << "WARNING: " << n_isolated << " atom(s) have no neighbour within the "
				<< rcut << " A descriptor cutoff, so there is no environment to predict from."
				<< " Isolated solvent is the usual cause." << std::endl;
		}
		std::cout << "I will fill out these atoms using spherical Thakkar densities!\n";
		// make a copy of initial wavefunction, to leave the initial one untouched!
		for (int a = ncen_in - 1; a >= 0; a--)
		{
			if (use_thakkar[a])
			{
				wavy.erase_atom(a);
			}
		}
		// remove all known basis sets, to not get problems with newly loaded ones
		for (int a = 0; a < wavy.get_ncen(); a++)
		{
			wavy.clear_atom_basis_set(a);
		}
		//std::filesystem::path new_fn = wavy.get_path().parent_path() / "SALTED_temp.xyz"; //I think this is not actually neccecary....
		//wavy.write_xyz(new_fn);
		//wavy.set_path(new_fn);
		opt.needs_Thakkar_fill = true;
		return use_thakkar;
	}
	return {};
}

std::string SALTED_Utils::FeatomicHyperParameters::to_json() const
{
	std::ostringstream oss;
	oss << "{\n"
		<< "  \"cutoff\": {  \
					\"radius\": " << this->cutoff_radius << " , \"smoothing\": \
								  {\"type\": \"" << this->cutoff_function.type << "\", \"width\": " << this->cutoff_function.width << "} }, \n"
		<< "  \"density\": { \
					\"type\": \"Gaussian\", \"width\": " << this->atomic_gaussian_width << ", \"center_atom_weight\": " << this->center_atom_weight << "},\n"
		<< "  \"basis\": { \
					\"type\": \"TensorProduct\", \"max_angular\": " << this->max_angular << ", \"radial\": {\"type\":  \"" << this->radial_basis.type << "\", \"max_radial\": " << this->max_radial << "} , \"spline_accuracy\": " << this->radial_basis.spline_accuracy << "}\n"
		<< "}";
	return oss.str();
}


// Runs featomic's spherical expansion for the SALTED keys
static metatensor::TensorMap get_feats_projs(featomic::SimpleSystem featomic_system, const SALTED_Utils::FeatomicHyperParameters& parameters, const std::vector<int32_t>& centres)
{
	// size_t nspe1 = neighspe.size();
	std::vector<std::array<int32_t, 4>> keys_array;
	keys_array.reserve((parameters.max_angular + 1) * parameters.species.size() * parameters.neighspe.size());

	for (int l = 0; l < parameters.max_angular + 1; ++l)
	{
		for (const std::string& center_spe : parameters.species)
		{
			int32_t center_z = constants::get_Z_from_label(center_spe.c_str()) + 1;
			for (const std::string& neigh_spe : parameters.neighspe)
			{
				int32_t neigh_z = constants::get_Z_from_label(neigh_spe.c_str()) + 1;
				// Directly emplace back initializer_lists into keys_array
				keys_array.push_back({ l, 1, center_z, neigh_z });
			}
		}
	}

	// Assuming metatensor::Labels expects a flat sequence of integers for each label
	std::vector<int32_t> flattened_keys;
	for (const auto& subVector : keys_array)
	{
		flattened_keys.insert(flattened_keys.end(), subVector.begin(), subVector.end());
	}

	// Convert keys_array to rascaline::Labels
	std::vector<std::string> names = { "o3_lambda", "o3_sigma", "center_type", "neighbor_type" };
	metatensor::Labels keys_selection(names, flattened_keys.data(), flattened_keys.size() / names.size());


	//create the calculator with its name and parameters
	//Do not ask me, why Featomic expects the max_radial to be one less than the actual number of radial basis functions, but it does, so here we are
	SALTED_Utils::FeatomicHyperParameters modif_param = parameters;
	modif_param.max_radial -= 1;
	//featomic's rayon pool is built on first use from RAYON_NUM_THREADS and otherwise
	//takes every logical core, ignoring -cpus and OMP_NUM_THREADS
#ifdef _OPENMP
	if (std::getenv("RAYON_NUM_THREADS") == nullptr) { // Flawfinder: ignore
		const std::string n = std::to_string(omp_get_max_threads());
#ifdef _WIN32
		_putenv_s("RAYON_NUM_THREADS", n.c_str());
#else
		setenv("RAYON_NUM_THREADS", n.c_str(), 1);
#endif
	}
#endif
	const auto t0 = std::chrono::steady_clock::now();
	auto calculator = featomic::Calculator("spherical_expansion", modif_param.to_json().c_str());
	const auto t1 = std::chrono::steady_clock::now();

	featomic::CalculationOptions calc_opts;
	calc_opts.selected_keys = keys_selection;
	calc_opts.use_native_system = true;
	if (!centres.empty())
		calc_opts.selected_samples = featomic::LabelsSelection::subset(metatensor::Labels(std::vector<std::string>{ "atom" }, centres.data(), centres.size()));
	// run the calculation
	// one block per (o3_lambda, o3_sigma, center_type, neighbor_type); get_expansion_coeffs
	// densifies them itself
	metatensor::TensorMap out = calculator.compute(featomic_system, calc_opts);
	if (ProgressBar::report_counts)
		std::cout << "[stages] featomic calculator " << std::chrono::duration<double>(t1 - t0).count() << " s, compute "
				  << std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count() << " s" << std::endl;
	return out;
}

static size_t label_column(const metatensor::Labels& labels, const char* name)
{
	const auto& names = labels.names();
	for (size_t i = 0; i < names.size(); ++i)
		if (std::strcmp(names[i], name) == 0) return i;
	throw std::runtime_error(std::string("featomic labels lack the dimension ") + name);
}

// Packs featomic's raw blocks into omega exactly as keys_to_samples("center_type") followed by
// keys_to_properties("neighbor_type") (both sorting samples) would have laid them out, without
// building the merged TensorMap and without a save/load buffer round trip, which together
// cost more than featomic's own compute of the blocks:
//  - rows are the sorted union of the centre atoms of the lambda = 0 blocks
//  - channels of one lambda are the neighbour types in first-key order, each followed by its n
//  - a sample a neighbour block lacks stays zero
// Every output element sums its r contributions in the same order as the merged loop did, so
// the result is bit-identical; blocks of one lambda write disjoint (row, channel) cells.
static SALTEDDescriptors get_expansion_coeffs(metatensor::TensorMap& descriptor, const SALTED_Utils::FeatomicHyperParameters& parameters)
{
	const metatensor::Labels keys = descriptor.keys();
	const int32_t* key = keys.values().data();
	const size_t key_width = keys.size(), n_blocks = keys.count();
	const size_t k_lambda = label_column(keys, "o3_lambda"), k_neighbor = label_column(keys, "neighbor_type");

	std::vector<int> row_of;
	int nchannels = 0;
	{
		std::vector<int32_t> neighbors;
		for (size_t b = 0; b < n_blocks; ++b) {
			if (key[b * key_width + k_lambda] != 0) continue;
			metatensor::TensorBlock block = descriptor.block_by_id(b);
			const metatensor::Labels samples = block.samples();
			const int32_t* s = samples.values().data();
			const size_t width = samples.size(), atom = label_column(samples, "atom");
			for (size_t i = 0; i < samples.count(); ++i) {
				const int32_t a = s[i * width + atom];
				if (static_cast<size_t>(a) >= row_of.size()) row_of.resize(static_cast<size_t>(a) + 1, -1);
				row_of[a] = 0;
			}
			const int32_t neighbor = key[b * key_width + k_neighbor];
			if (std::find(neighbors.begin(), neighbors.end(), neighbor) == neighbors.end()) {
				neighbors.push_back(neighbor);
				nchannels += static_cast<int>(block.properties().count());
			}
		}
	}
	int n_atoms = 0;
	for (int& r : row_of)
		if (r >= 0) r = n_atoms++;

	SALTEDDescriptors omega(n_atoms, nchannels, parameters.max_angular);
	for (int l = 0; l <= parameters.max_angular; ++l)
	{
		const int m = 2 * l + 1;
		const cvec2 c2r = SALTED_Utils::complex_to_real_transformation({ m })[0];
		std::vector<std::pair<int32_t, int>> channel_of; // neighbour type -> first channel
		int next_channel = 0;
		for (size_t b = 0; b < n_blocks; ++b) {
			if (key[b * key_width + k_lambda] != l) continue;
			metatensor::TensorBlock block = descriptor.block_by_id(b);
			const metatensor::Labels samples = block.samples();
			const metatensor::NDArray<double> values = block.values();
			const int nprop = static_cast<int>(values.shape()[2]);
			const int32_t neighbor = key[b * key_width + k_neighbor];
			auto it = std::find_if(channel_of.begin(), channel_of.end(), [&](const auto& e) { return e.first == neighbor; });
			if (it == channel_of.end()) {
				channel_of.emplace_back(neighbor, next_channel);
				it = channel_of.end() - 1;
				next_channel += nprop;
			}
			const int first_channel = it->second;
			const int32_t* s = samples.values().data();
			const size_t width = samples.size(), atom = label_column(samples, "atom");
			const double* v = values.data();
			const int n_samples = static_cast<int>(samples.count());
#pragma omp parallel for schedule(static)
			for (int i = 0; i < n_samples; ++i) {
				const int row = row_of[s[i * width + atom]];
				for (int r = 0; r < m; ++r) {
					const double* value = v + (static_cast<size_t>(i) * m + r) * nprop;
					for (int d = 0; d < nprop; ++d) {
						cdouble* output = omega.block(row, first_channel + d, l);
						for (int c = 0; c < m; ++c)
							output[c] += conj(c2r[r][c]) * value[d];
					}
				}
			}
		}
	}

	return omega;
}


SALTEDDescriptors SALTED_Utils::calculate_SALTED_descriptors(const featomic::SimpleSystem& featomic_system, const SALTED_Utils::FeatomicHyperParameters& parameters, const std::vector<int32_t>& centres)
{
	metatensor::TensorMap descriptor = get_feats_projs(featomic_system, parameters, centres);
	return get_expansion_coeffs(descriptor, parameters);
}


namespace
{
	// Constructing a calculator splines the radial integral for every (n, l)
	// pair to `spline_accuracy`, which depends only on the hyperparameters, so
	// it is computed once and kept.
	//
	// Keyed on the parameter JSON, so a caller that changes any hyperparameter
	// gets a new calculator rather than a silently wrong one -- a descriptor of
	// the right shape computed with the wrong settings.
	//
	// `featomic::Calculator` is move-only, hence the indirection, and it is not
	// safe to use one calculator from several threads at once. The only caller
	// of `calculate_SOAP_Powerspectrum` is the serial descriptor path; SALTED
	// builds its own calculator in `get_feats_projs`.
	featomic::Calculator& cached_calculator(const std::string& name, const std::string& json)
	{
		static std::map<std::string, std::unique_ptr<featomic::Calculator>> cache;
		const std::string key = name + '\n' + json;
		auto found = cache.find(key);
		if (found == cache.end())
		{
			found = cache.emplace(key, std::make_unique<featomic::Calculator>(name, json)).first;
		}
		return *found->second;
	}
}

//FEATOMIC POWER Spectrum
metatensor::TensorMap SALTED_Utils::calculate_SOAP_Powerspectrum(featomic::SimpleSystem featomic_system, const SALTED_Utils::FeatomicHyperParameters& parameters) {
	// Phase timings, off unless NOSPHERA2_TIME_SOAP is set. Caching the
	// calculator does not remove the fixed per-call cost; it is in what follows.
	const bool time_phases = tuning("NOSPHERA2_TIME_SOAP") != nullptr;
	auto mark = std::chrono::steady_clock::now();
	auto lap = [&mark](const char* what, bool on) {
		if (!on) return;
		const auto now = std::chrono::steady_clock::now();
		std::cout << "  SOAP_PHASE " << what << " "
				  << std::chrono::duration<double>(now - mark).count() << std::endl;
		mark = now;
	};

	// Built once per set of hyperparameters and kept: see `cached_calculator`.
	auto& calculator = cached_calculator("soap_power_spectrum", parameters.to_json());
	lap("calculator", time_phases);

	std::vector<std::array<int32_t,3>> keys_array;
	//keys for centre types the system does not contain cost only metatensor bookkeeping, which for the full 726 is a large share of a call
	std::set<int32_t> present(featomic_system.types(), featomic_system.types() + featomic_system.size());
	for (const std::string& center_type : parameters.species)
	{
		int32_t z_center = constants::get_Z_from_label(center_type.c_str()) + 1;
		if (!present.count(z_center)) continue;

		for (size_t i = 0; i < parameters.neighspe.size(); ++i)
		{
			int32_t z1 = constants::get_Z_from_label(parameters.neighspe[i].c_str()) + 1;

			for (size_t j = i; j < parameters.neighspe.size(); ++j)
			{
				int32_t z2 = constants::get_Z_from_label(parameters.neighspe[j].c_str()) + 1;

				keys_array.push_back({ z_center, z1, z2 });
			}
		}
	}

	// Assuming metatensor::Labels expects a flat sequence of integers for each label
	std::vector<int32_t> flattened_keys;
	for (const auto& subVector : keys_array)
	{
		flattened_keys.insert(flattened_keys.end(), subVector.begin(), subVector.end());
	}

	// Convert keys_array to rascaline::Labels
	std::vector<std::string> names = {"center_type", "neighbor_1_type", "neighbor_2_type"};
	metatensor::Labels keys_selection(names, flattened_keys.data(), flattened_keys.size() / names.size());

	featomic::CalculationOptions calc_opts;
	calc_opts.use_native_system = true;
	calc_opts.selected_keys = keys_selection;
	// run the calculation
	// Initialize descriptor directly from computation result
	lap("keys", time_phases);
	metatensor::TensorMap descriptor = calculator.compute(featomic_system, calc_opts);
	lap("compute", time_phases);

	// The descriptor is a metatensor `TensorMap`, containing multiple blocks.
	// We can transform it to a single block containing a dense representation,
	// with one sample for each atom-centered environment.
	descriptor = descriptor.keys_to_samples("center_type");
	lap("keys_to_samples", time_phases);
	descriptor = descriptor.keys_to_properties(svec{ "neighbor_1_type" , "neighbor_2_type" });
	lap("keys_to_properties", time_phases);

	return descriptor;
}

void calc_cube_ML(const vec& data, WFN& dummy, cube& cube_data, const int& atom_nr)
{
	_time_point start = get_time();

	const int s1 = cube_data.get_size(0), s2 = cube_data.get_size(1), s3 = cube_data.get_size(2), total_size = s1 * s2 * s3;
	std::cout << "Lets go into the loop! There is " << total_size << " points" << std::endl;

	ProgressBar* progress = new ProgressBar(total_size, 60, "=", " ", "Calculating Values");
	vec v1{
		cube_data.get_vector(0, 0),
		cube_data.get_vector(1, 0),
		cube_data.get_vector(2, 0) },
		v2{
			cube_data.get_vector(0, 1),
			cube_data.get_vector(1, 1),
			cube_data.get_vector(2, 1) },
			v3{
				cube_data.get_vector(0, 2),
				cube_data.get_vector(1, 2),
				cube_data.get_vector(2, 2) },
				orig{
					cube_data.get_origin(0),
					cube_data.get_origin(1),
					cube_data.get_origin(2) };

	if (atom_nr != -1)
		std::cout << "Calculation for atom " << atom_nr << std::endl;

	const std::vector<atom> atoms = dummy.get_atoms();
	//atom_nr selects one atom: its own table and the slice of the coefficients that belongs to it
	const aux_density_table full(atoms);
	const aux_density_table t(atom_nr == -1 ? atoms : std::vector<atom>{ atoms[atom_nr] });
	vec coefs = data;
	if (atom_nr != -1) {
		const int off = full.coef_off[full.sh_start[atom_nr]];
		coefs.assign(data.begin() + off, data.begin() + off + t.n_coef);
	}
#pragma omp parallel for schedule(dynamic)
	for (int index = 0; index < total_size; index++)
	{
		int i = index / (s2 * s3);
		int j = (index / s3) % s2;
		int k = index % s3;

		cube_data.set_value(i, j, k, t(
			i * v1[0] + j * v2[0] + k * v3[0] + orig[0],
			i * v1[1] + j * v2[1] + k * v3[1] + orig[1],
			i * v1[2] + j * v2[2] + k * v3[2] + orig[2], coefs.data()));
		progress->update();
	}
	delete (progress);

	using namespace std;
	_time_point end = get_time();
	if (get_sec(start, end) < 60)
		std::cout << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) << " s" << endl;
	else if (get_sec(start, end) < 3600)
		std::cout << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 60 << " m " << get_sec(start, end) % 60 << " s" << endl;
	else
		std::cout << "Time to calculate Values: " << fixed << setprecision(0) << get_sec(start, end) / 3600 << " h " << (get_sec(start, end) % 3600) / 60 << " m" << endl;
	cube_data.calc_dv();
	std::cout << "Number of electrons: " << std::fixed << std::setprecision(4) << cube_data.sum() << std::endl;
};

cube calc_cube_ML(const vec& data, WFN& dummy, const int& atom_nr)
{
	properties_options opts;
	readxyzMinMax_fromWFN(dummy, opts);
	cube CubeRho(opts.NbSteps, dummy.get_ncen(), true);
	CubeRho.give_parent_wfn(dummy);

	for (int i = 0; i < 3; i++)
	{
		CubeRho.set_origin(i, opts.MinMax[i]);
		CubeRho.set_vector(i, i, (opts.MinMax[i + 3] - opts.MinMax[i]) / opts.NbSteps[i]);
	}
	CubeRho.set_comment1("Calculated density using NoSpherA2 from ML Data");
	CubeRho.set_comment2("from " + dummy.get_path().string());
	CubeRho.set_path((dummy.get_path().parent_path() / dummy.get_path().stem()).string() + "_RI_rho.cube");

	calc_cube_ML(data, dummy, CubeRho, atom_nr);

	return CubeRho;
};

#include "integrator.h"
#include "libCintMain.h"
#include "nos_math.h"
#include "npy.h"
void create_SALTED_training_data(const WFN& orbital, const WFN& aux, const options& opts) {
	std::cout << "Calculating density fitting coefficients..." << std::endl;
	DensityFitting::CONFIG config = DensityFitting::config_from_options(opts);
	config.analyze_quality = true;

	// SALTED learns one coefficient block per atom; keep all restrained moments on their centres.
	if (config.grid_higher_moments) {
		std::cout
			<< "SALTED training: restraining higher moments on their centres."
			<< std::endl;
		config.grid_higher_moments = false;
	}

	vec coefs = DensityFitting::density_fit(orbital, aux, config);

	vec overlap;
	Int_Params aux_basis(aux);
	compute2C<Overlap2C_SPH>(aux_basis, overlap);
	const int nao_max = aux_basis.get_nao();

	dMatrix1 coefs_vec(coefs.size());
	coefs_vec.container() = coefs;
	dMatrix2 overlap_mat(nao_max, nao_max);
	overlap_mat.container() = overlap;

	npy::write_npy("coefficients.npy",
		npy::npy_data<double>{
			coefs,
			{ static_cast<unsigned long>(coefs.size()) },
			false}
	);

	npy::write_npy("overlap.npy",
		npy::npy_data<double>{
		overlap,
		{ static_cast<unsigned long>(nao_max), static_cast<unsigned long>(nao_max) },
			false}
	);

}
