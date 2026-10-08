#pragma once
#include <filesystem>

#include "convenience.h"
#include "SALTED_io.h"
#include "SALTED_utilities.h"
#include "wfn_class.h"

//The SALTED prediction for a structure as a flat table with its coefficients, so rho and the ESP
//can be evaluated at points without a wavefunction (xyz input); opt.salted_model_dir names the model
class SALTEDPredictor
{
public:
	SALTEDPredictor(WFN wavy, options &opt_in);
	SALTEDPredictor() = default;

	const std::string get_dfbasis_name() const;
	vec gen_SALTED_densities();
	const std::filesystem::path get_salted_filename() const
	{
		return config.salted_filename;
	};
	WFN wavy;
	void shrink_intermediate_vectors();
	const bool basis_set_loaded() const { return bbasis_set_loaded; };
	//The basis this model predicts in: its own BASIS block if it has one, else the library set it names
	std::shared_ptr<BasisSet> get_model_basis() const;

private:
	bool bbasis_set_loaded = false;
	//An atom's coefficients mean something only with the basis they were trained on; kept for stitching models
	std::shared_ptr<BasisSet> model_basis{};
	// Several models: each element goes to the first model on the command line trained on it, and the
	// per-atom blocks are stitched back together. Empty for a single model.
	std::vector<std::unique_ptr<SALTEDPredictor>> sub_models{};
	// Per atom of the merged structure: which sub model predicts it, and its index there
	ivec atom_model{}, atom_in_model{};
	void build_merged(const WFN& wavy_in, options& opt_in);
	vec merge_predictions();
	// Estimate what the spherically filled atoms should carry, so the size of the
	// neutral-fill assumption is reported rather than hidden
	void estimate_fill_charges(const WFN& wavy_in, const std::vector<char>& use_thakkar, options& opt_in);
	// Set when atoms were moved to the spherical Thakkar fill. The charge
	// constraint needs it: with a mixed ML/Thakkar system the split of a net
	// charge between the two regions is undefined.
	bool spherical_fill_used = false;
	// EEQ estimate of the charge sitting on the spherically filled atoms, and
	// how many there are. Used to say how wrong the neutral-fill assumption is.
	double filled_eeq_charge = 0.0;
	// The part of that charge the spherical fill can actually carry (an ion
	// must be tabulated for the element). The ML target is shifted by exactly
	// this, so predicted + filled still sums to the right number of electrons.
	double applied_fill_charge = 0.0;
	int n_filled = 0;
	SALTEDConfig config;
	int natoms;
	std::filesystem::path SALTED_DIR;
	std::filesystem::path coef_file;
	bool debug;
	std::vector<std::string> atomic_symbols{};
	std::unordered_map<std::string, ivec> atom_idx{};

	std::unordered_map<std::string, int> natom_dict{}, lmax{}, nmax{};
	SALTEDDescriptors v1, v2;
	// Both hyperparameter sets identical: v2 is conj(v1), never filled, equicomb conjugates on read
	bool v2_is_conj_of_v1 = false;
	void setup_atomic_environment();

	vec weights{};
	std::unordered_map<std::string, dMatrix2> Vmat{};
	// Projector shapes for every species the model knows, present or not: the flat
	// weight vector is laid out over all of them, so absent widths shift the offsets
	std::unordered_map<std::string, std::array<size_t, 2>> proj_dims{};
	// Start of each present (species, l) block in that flat vector: n-major, one projector width per n
	std::unordered_map<std::string, size_t> weight_offset{};
	// A VERSION 4 file (PROJW) stores the projectors with the weights already folded in
	bool projector_folded = false;
	// Held open for the whole prediction, indexed by (species, lambda)
	std::unique_ptr<SALTED_BINARY_FILE> model_file{};
	std::unordered_map<std::string, SALTED_BINARY_FILE::block_ref> feat_index{};
	std::unordered_map<std::string, SALTED_BINARY_FILE::block_ref> proj_index{};
	std::unordered_set<std::string> model_species{};
	std::unordered_map<std::string, int> Mspe{};
	std::unordered_map<int, std::vector<int64_t>> vfps{};
	std::unordered_map<int, vec> wigner3j{};
	std::unordered_map<std::string, dMatrix2> power_env_sparse{};
	std::unordered_map<std::string, vec> av_coefs{};
	std::unordered_map<int, int> featsize{};

	featomic::SimpleSystem featomic_system;
	void read_model_data();
	// Drop the model matrices of one lambda; read_model_data() only indexes
	void free_model_lambda(const int lam);
	// (species+lambda key, projector, features); read and install split for the prefetch
	using lambda_blocks = std::vector<std::tuple<std::string, dMatrix2, dMatrix2>>;
	lambda_blocks read_model_lambda(const int lam);
	void install_model_lambda(lambda_blocks blocks);

	vec predict();
};

// Writes the model `in` to `out` as VERSION 4: PROJW holds what install_model_lambda() folds on every run
void fold_salted_file(const std::filesystem::path& in, const std::filesystem::path& out);

