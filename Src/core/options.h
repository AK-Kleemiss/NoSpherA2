#pragma once
#include "nos_types.h"
#include "properties_options.h"
#include "hkl.h"
#include "extinction.h"

//---------------- Object for handling all input options -------------------------------
struct options
	/** @brief All command line options and settings controlling a run. */
{
	std::ostream &log_file;
	double d_sfac_scan = 0.0;
	d3 sfac_diffuse = { 0.0, 0.0, 0.0 };
	double dmin = 99.0;
	double mem = 1000.0; // In MB
	//Set only when -mem was passed; only then is mem a budget the tsc block size and XCW's stored integrals are fitted to
	bool mem_given = false;
	double efield = 0.005;
	ivec2 groups;
	ivec2 hkl_min_max{ {-100, 100}, {-100, 100}, {-100, 100} };
	int hkl_scale = 1;
	vec2 twin_law;
	ivec2 combined_tsc_groups;
	pathvec combined_tsc_calc_files;
	pathvec combined_tsc_calc_cifs;
	std::vector<unsigned int> combined_tsc_calc_mult;
	ivec combined_tsc_calc_charge;
	ivec combined_tsc_calc_ECP;
	//Bounds-checked: the digesters read arguments[i + n] freely, so a flag that is last on the
	//line throws missing_argument (caught in digest_options with the flag's name) instead of
	//reading past the end
	struct missing_argument : std::out_of_range
	{
		using std::out_of_range::out_of_range;
	};
	struct checked_svec : svec
	{
		using svec::svec;
		std::string &operator[](size_t i) { return const_cast<std::string &>(std::as_const(*this)[i]); }
		const std::string &operator[](size_t i) const
		{
			if (i >= size())
				throw missing_argument("argument " + std::to_string(i));
			return svec::operator[](i);
		}
	};
	checked_svec arguments;
	pathvec combine_mo;
	svec Cations;
	svec Anions;
	pathvec pol_wfns;
	ivec cmo1;
	ivec cmo2;
	ivec ignore;
	std::filesystem::path salted_model_dir;
	//Every model given to -SALTED, in that order; salted_model_dir is the first of them
	pathvec salted_model_dirs;
	std::vector<std::shared_ptr<BasisSet>> aux_basis;
	std::filesystem::path wfn;
	std::filesystem::path wfn2;
	std::filesystem::path cube_density;
	std::filesystem::path fchk;
	std::string basis_set;
	std::filesystem::path hkl;
	std::filesystem::path cif;
	std::string method;
	std::filesystem::path xyz_file;
	std::filesystem::path coef_file;
	std::filesystem::path hirshfeld_surface;
	std::filesystem::path hirshfeld_surface2;
	std::filesystem::path fract_name;
	std::filesystem::path wavename;
	std::filesystem::path gaussian_path;
	std::filesystem::path turbomole_path;
	std::filesystem::path basis_set_path;
	std::filesystem::path mask_file;
	std::string occ;
	std::filesystem::path occ_toml_path;
	std::filesystem::path cwd;
	std::filesystem::path profiling_tests_root = "tests";
	pathvec promol_nci_xyz; //Two or more fragments; a grid point is intermolecular when no single fragment dominates
	//Geometry-aid jobs (-calc_featomic_descriptor(s), -classify_atoms(_list)): the flags queue, run_app_impl runs geometry_aid::run and quits
	bool calc_featomic_descriptor = false;
	std::filesystem::path classify_atoms_out;
	std::filesystem::path geometry_aid_model;
	pathvec featomic_structures;
	pathvec classify_structures;
	double geometry_aid_cutoff = 3.5;
	bool geometry_aid_metals = false;
	double geometry_aid_center_weight = 1.0;
	std::filesystem::path interaction_energies_job;
	std::filesystem::path xcw_settings_path;
	properties_options properties;
	bool debug = false;
	//Set by the one-shot options (-merge, -dipole_moments, -convert_to_47, ...) that used to
	//exit(0) inside the parser: that killed the host when run_app is a library call (Olex2, the
	//in-process tests), and no option after them was read. run_app_impl returns 0 instead.
	bool finished = false;
	bool all_charges = false;
	bool SALTED = false;
	bool Olex2_1_3_switch = false;
	bool iam_switch = false;
	bool read_k_pts = false;
	bool save_k_pts = false;
	bool combined_tsc_calc = false;
	bool binary_tsc = true;
	bool cif_based_combined_tsc_calc = false;
	bool no_date = false;
	bool gbw2wfn = false;
	bool old_tsc = false;
	bool write_tshb = false;
	bool label_tsc_output = false;
	bool write_CIF = false;
	bool test = false;
	bool electron_diffraction = false;
	bool ECP = false;
	bool RI_FIT = false;
	bool needs_Thakkar_fill = false;
	//Set around a spherical fill, where a disorder part already covered by an earlier one legitimately yields no atoms and must not read as a broken CIF
	bool allow_empty_asym = false;
	//Set while a spherical fill runs for somebody else: it must RETURN a block, not stream experimental.tscb out from under the table being built
	bool spherical_fill = false;
	// Per-atom EEQ charges for atoms handed to the spherical fill, as
	// {x, y, z, q} in the wavefunction's own coordinate units. Keyed by
	// POSITION rather than index because the fill rebuilds its wavefunction
	// from the original file, and that is how CIF and WFN atoms are matched
	// everywhere else here. Empty means "no charges known" -> neutral fill,
	// which is the previous behaviour.
	std::vector<std::array<double, 4>> spherical_fill_charges{};
	//Reflections per block when streaming the tsc; 0 restores the single allocation of the whole scatterers x reflections x 16 byte table
	//The exact block size is not performance-critical
	size_t tsc_block_size = 1000;
	//Set only when -tsc_block was passed, so an explicit block size wins over one derived from -mem
	bool tsc_block_given = false;
	//Reflections to hold at once for n_scat scatterers, 0 for the whole table; -tsc_block wins, then -mem, then the default
	//A block costs about three copies of n_scat * block * 16 bytes - producer, queue and writer - which is what the budget is spent against
	size_t tsc_block_for(const size_t n_refl, const size_t n_scat) const
	{
		if (tsc_block_given || !mem_given || mem <= 0.0)
			return tsc_block_size;
		const size_t item = 3 * (n_scat ? n_scat : 1) * sizeof(std::complex<double>);
		return items_within_budget(n_refl, item, static_cast<size_t>(mem * 1024.0 * 1024.0));
	}

	//set once a streamed run wrote the file itself, so the caller does not overwrite it with an empty one-shot block
	bool tsc_written_by_stream = false;
	//set when -tshb wrote experimental.tshb instead of the table, so the caller writes no (empty) table either
	bool tsc_replaced_by_tshb = false;
	//the run's .tshb while every part goes into it instead of the table; null otherwise
	tshb_block* tshb = nullptr;
	bool qct = false;
	// What structure_factors reads: the crystal model, the extinction and the grid precision
	struct SF_settings {
		// `i_sigma <x>`: only reflections with I/sigma(I) >= x enter chi^2 and the scale, as in
		// Tonto. Under the reader's sigma(F) = sigma(I)/2F that is F/sigma(F) >= 2x, so the
		// default 2 is SHELX's I > 2 sigma(I) and F > 4 sigma(F) at once
		double i_sigma_cutoff = 2;
		std::string basis_set_name;
		//`basis_overrides <file>`: per-element replacements of the orbital basis from a BSE JSON file
		std::filesystem::path basis_overrides;
		bool grown = false;
		bool mask = false;
		//`extinction <shelx|bc_gaussian|bc_lorentzian> [iso|aniso] [fixed] [start value]`; the start is one x,
		//or with aniso the six X11 X22 X33 X12 X13 X23 of the tensor
		extinction::model extinction_model = extinction::model::none;
		bool extinction_aniso = false;
		bool extinction_refine = true;
		vec extinction_start = { 1e-4 };
		//`wavelength <lambda>`, in Angstrom; overrides _diffrn_radiation_wavelength from the CIF
		double wavelength = 0.0;
		int n_params = 1;
		int refine_against = 1;
		int XWR_type = 1;
		// `weight [robust] sigma | unit | stl <a> | shelx <a> <b> | shelx_opt <a0> <b0> | chebyshev <A0..> |
		// chebyshev_opt <n> [c]` (CLI -xcw_weight): the crystallographic weight of every fit sum, see
		// structure_factors::update_cryst_weights. A new scheme is one more value here and one more case there;
		// weight_params holds its parameters in the order the keyword reads them (chebyshev_opt: n).
		// robust_c > 0 puts Tukey's biweight with that c on top of the scheme
		enum class weight_scheme { sigma, shelx, shelx_opt, unit, stl, chebyshev, chebyshev_opt };
		weight_scheme weight = weight_scheme::sigma;
		vec weight_params;
		double robust_c = 0.0;
		// `weight_rescale off` (CLI -xcw_weight_rescale 0) keeps unit/stl/chebyshev weights absolute; by default
		// they are multiplied once, at the first fill, by the k that makes the scan's criterion S = 1 (CRYSTALS REWEIGHT)
		bool weight_rescale = true;
		// `i_double` in the settings file: the XCW grid's AO values and block products in double precision
		bool grid_double = false;
	};

	// What the SCF solver reads: the convergence thresholds with their running state, the damping
	// and level shift schedule, the spin kind and the second-order switches
	struct SCF_settings {
		//The thresholds default to the `normal` preset of the XCW settings file
		double quant_diff = 1e-6;
		double current_quant_diff = 0;
		bool conv_quant_diff = false;
		double max_diis_error = 1e-5;
		double current_max_diis_error = 0;
		bool conv_max_diis_error = false;
		double gradient = 7e-5;
		double current_gradient = 0;
		bool conv_gradient = false;
		double RMSP_diff = 1e-6;
		double current_RMSP_diff = 0;
		bool conv_RMSP_diff = false;
		double MaxP_diff = 1e-5;
		double current_MaxP_diff = 0;
		bool conv_MaxP_diff = false;
		//The four damping/shift values default to off: `fast_conv` only clears the apply
		//flags, yet run_lambda copies alpha_damp and level_shift and SCF_iteration mixes the
		//density with alpha_damp regardless of the flag
		double diis_stop_damping = 0;
		//`check_hessian`: finite-difference check of the Hessian-vector product on the first
		//second-order step, reported in SCF.log; it enters TRAH once the DIIS error is below
		//trah_.start_threshold, so a cold start has a step to check
		bool check_hessian = false;
		bool apply_shift = true;
		bool method_apply_shift = true;
		double diis_stop_shift = 0;
		bool apply_damping = true;
		bool method_apply_damping = true;
		//`guess_basis <name>`: the first lambda starts from a Hartree-Fock converged in this
		//(smaller) basis by OCC's own driver, its density projected into the orbital basis
		std::string guess_basis_name;
		occ::qm::SpinorbitalKind hf_type = occ::qm::SpinorbitalKind::Restricted;
		double alpha_damp = 0;
		double level_shift = 0;
		int max_scf_iterations = 100;
		int charge = 0;
		int multiplicity = 1;
		// Not settings-file keys: filled by XCW_wrapper from the options when it builds the SCF solver
		bool incremental = true;
		bool gpu_eri = false;
		bool no_date = false;
		// `save_eri <path>` / `read_eri <path>`: the stored two-electron integrals, for a later run on the same geometry and basis
		std::filesystem::path eri_read_path, eri_save_path;

		// Clears the convergence flags
		void clear() {
			conv_quant_diff = false;
			conv_max_diis_error = false;
			conv_gradient = false;
			conv_RMSP_diff = false;
			conv_MaxP_diff = false;
			apply_shift = method_apply_shift;
			apply_damping = method_apply_damping;
		}

		// Performs the convergence check and sets flags accordingly
		bool convergence_check() const {
			if (conv_quant_diff == true && conv_max_diis_error == true && conv_gradient == true && conv_RMSP_diff == true && conv_MaxP_diff == true) {
				return true;
			}
			else {
				return false;
			}
		}

		// Updates the SCF routine in regards to damping and level shift
		void update(std::ostream& file, double& alpha_damp) {
			if (current_max_diis_error < diis_stop_damping && apply_damping == true) {
				apply_damping = false;
				alpha_damp = 0;
				print_centered_message("***Turned off damping***", 100, file);
			}
			if (current_max_diis_error < diis_stop_shift && apply_shift == true) {
				apply_shift = false;
				print_centered_message("***Turned off level shift***", 100, file);
			}
		}
	};

	// What XCW_wrapper reads itself; the two slices above go to structure_factors and the SCF solver
	struct XCW_settings : SF_settings, SCF_settings {
		//`df_basis <name>`: density fitting of the Fock build with this auxiliary basis
		std::string df_basis_name;
		//`hessian_df_basis <name>`: density fitting of the TRAH Hessian products only, iterations keep exact Fock builds.
		//Unset, def2-universal-jkfit is used by itself when the ERIs do not fit in memory; `none` keeps them exact
		std::string hessian_df_basis_name;
		//`scf_backend occ|dusk` (or -scf_backend): the SCF of the fit, occ's XCW_solver or molscf (dusk_scf.h)
		std::string scf_backend = "occ";
		double xcw_start_value;
		int num_xcw_steps;
		double xcw_step_size;
		bool read_first_guess;
		bool nbo_output = false;
		bool xcw_gaussian_halt = false;
		double xcw_strong_cutoff = 3.0;
	};

	XCW_settings xcw_settings;
	bool do_XCW = false;
	//-xcw_partition: how XCW shares rho among the atoms, "grid" (the -hirsh/-becke/-tfvc weights) or "ri" (an RI fit)
	std::string xcw_partition = "grid";
	bool rgbi = false;
	//-npa: NAO/NPA, run in-process
	bool npa = false;
	//per-NAO occupancy table beside the NPA; -npa_summary turns it off
	bool npa_orbitals = true;
	//-ibo: IAO charges and intrinsic bond orbitals, run in-process
	bool ibo = false;
	bool rgbi_no_sym = false;
	bool rgbi_EVs = false;
	bool rgbi_theta = false;
	//-rgbi_legacy_cutoff: atomic subspace by occupation threshold (1/6 NAO, 1/14 ANO) instead of the
	//free-atom orbital count; the projector rank, and so every bond index, jumps when an occupation
	//crosses it. Only for reproducing older numbers.
	bool rgbi_legacy_cutoff = false;
	RGBIOrbitalBasis rgbi_orbital_basis = RGBIOrbitalBasis::ANO;
	ivec3 rgbi_group_sets;
	bool fract = false;
	//GPU scattering factors when a device is present; -no_gpu forces the CPU loop
	bool use_gpu = true;
	//-gpu_fp64 keeps the double sincos on a card that would otherwise pick the fp32 one
	bool gpu_fp64 = false;
	//Single-precision AO values and block products on the XCW grid, as the device path runs; sgemm is twice
	//dgemm's rate. -no_cpu_itensor_fp32 keeps double.
	bool cpu_itensor_fp32 = true;
	//Seed each lambda step from the density extrapolated through the two previous steps
	//rather than the last one alone; the step is small and the trajectory smooth.
	bool xcw_extrapolate = true;
	//Build the two-electron Fock matrix from the change of the density between iterations,
	//rather than from the whole density every time: the stored integrals skip the segments
	//the difference cannot reach (stored_eri::JK), the direct kernel skips shell quartets.
	//The direct build only with xcw_int_precision 1e-12: at 1e-10 the increments' screening
	//error accumulates to a gradient floor of 3e-5 and the SCF never meets its 1e-5, and the
	//full build at 1e-10 is the faster of the two anyway.
	bool xcw_incremental = false;
	//-print_F: every converged lambda step writes its F_calc table, NA2_<lambda>_Fcalc.txt
	bool print_F = false;
	//Integral screening threshold of the XCW Fock build; OCC's own default is 1e-12.
	double xcw_int_precision = 1e-10;
	//-gflops reports achieved GFLOP/s per stage for the CPU and GPU paths at the end of a
	//run. The thresholds deciding what goes to the device were calibrated on one machine;
	//this is how they get re-derived on another.
	bool track_gflops = false;
	//-gpu_fp32 forces the reduced-argument single-precision sincos on a card that would
	//otherwise keep the double one. It exists for the test suite: without it the precision
	//a run uses depends on the card, so neither path can be pinned. -gpu_fp64 wins if both
	//are given, the accurate path being the safer thing to fall back to.
	bool gpu_fp32 = false;
	//-no_gpu_itensor keeps the XCW grid GEMMs and the stored integrals' Fock build on the CPU. On by
	//default: it moves the total energy only in the tenth significant figure. Read together with
	//use_gpu, so -no_gpu turns it off as well.
	bool gpu_itensor = true;
	//SALTED descriptor combination uses the device when one is available; -no_gpu_salted keeps it on the CPU.
	bool gpu_salted = true;
	//-no_gpu_grid keeps the Becke/TFVC integration weights on the CPU
	bool gpu_grid = true;
	//Owned by the caller; the scattering-factor grid is built in it instead of a local, so
	//a second table for the same geometry reuses the points and weights
	GridManager* grid_cache = nullptr;
	//-no_gpu_density keeps the fitted density of the Gordon-Kim repulsion grid and the spherical-atom grids on the CPU
	bool gpu_density = true;
	//-gpu_blas offers large dense GEMMs in nos_math to the device
	bool gpu_blas = false;
	//The XCW grid GEMMs go through cuBLAS when the machine has it, and through the
	//built-in CUTLASS path otherwise. cuBLAS is 1.65x faster on a V100 and level with
	//CUTLASS within measurement noise on consumer cards, so preferring it costs nothing
	//where it does not help. -no_gpu_cublas pins CUTLASS, which is what the reference tests
	//do: the two differ in the last digits, so a test left to pick would pass or fail on
	//whether a CUDA toolkit happened to be installed.
	bool gpu_cublas = true;
	//Standalone conceptual-DFT reactivity analysis (-fukui_analysis), run from run_app_impl rather than at parse time so its output survives
	bool fukui_analysis_run = false;
	//EQC decomposition (-eqc), run from run_app_impl like -fukui_analysis; fragments are 0-based atoms, charge, mult
	bool eqc = false;
	struct eqc_fragment
	{
		ivec atoms;
		int charge = 0;
		int mult = 1;
	};
	std::vector<eqc_fragment> eqc_frags;
	std::vector<std::filesystem::path> eqc_wfns;
	std::string eqc_method = "hf";  //occ's method name: hf or a functional
	std::string eqc_basis;          //a basis-library name; empty takes the basis out of the gbw
	bool eqc_cold = false;          //also converge everything from occ's own guess, to count what seeding saves
	//Basin analysis (-eli_analysis), run from run_app_impl for the same reason
	bool eli_analysis_run = false;
	//-fba: RGBI, native NBO/NPA with NRT, bondwise Laplacian, QTAIM and ELI-D on one wavefunction
	bool fba = false;
	//-partition: charges and atomic dipoles of the selected partition, written to <stem>.xyzq
	bool partition_run = false;
	bool profiling = false;
	bool promol_nci = false;
	bool get_g = false;
	int accuracy = 2;
	//-basin_grid <n>: the quadrature of the basin analysis pulled into the core, tightest
	//exponent sharpened n^2-fold, radial step divided by n, Lebedev order up n - 1 entries
	int basin_grid = 1;
	//-basin_cube: QTAIM and ELI-D basins on the cube instead of from the analytic critical points and field
	bool basin_cube = false;
	//-eli_basin_cube <A> <B>: with -eli_analysis, a cube of ELI-D inside the bond basins between atoms A and B (an element
	//like Zn or a label like Zn1) and zero outside, cropped to them, at the -eli_analysis resolution
	svec eli_basin_atoms;
	//-eli_only: -eli_analysis without the density critical points and QTAIM basins
	bool eli_only = false;
	//-qtaim <wfn>: -eli_analysis without the ELI-D basins
	bool qtaim_only = false;
	//-eli_contributions: with -eli_analysis, what each ELI-D basin's electrons come from (QTAIM atoms, PM/Boys IBOs, NLMOs)
	bool eli_contributions = false;
	//-eli_mask_cube: ELI-D on the nonzero voxels of this cube's grid, nothing else
	std::filesystem::path eli_mask_cube;
	//-no_spin_eli: skip the spin-resolved ELI-D basins of a spin-polarised wavefunction
	bool spin_eli = true;
	//-npa x.gbw: the wavefunctions written right after flags that otherwise read -wfn; each must be -wfn's
	pathvec flag_wfns;
	//a flag without its own file took -wfn from later on the line, so that -wfn is no ignored option
	bool wfn_read_ahead = false;
	int threads = -1;
	int pbc = 0;
	int charge = 0;
	int ECP_mode = 0;
	PartitionType partition_type = PartitionType::Hirshfeld;
	//-multipole_moments: the RI fit is restrained to this scheme's atomic moments up to this order, -1 = unrestrained
	int multipole_lmax = -1;
	MultipoleScheme multipole_scheme = MultipoleScheme::HIRSHFELD;
	double multipole_strength = 1.0;
	//-repulsion_overlap: exchange-repulsion of -interaction_energy as K * Int rhoA rhoB, 0 = not included
	double repulsion_overlap = 0.0;
	//-repulsion_exchange: exchange functional of the Gordon-Kim repulsion, 0 Dirac, 1 PBE, 2 B88, 3 r2SCAN-L
	int repulsion_exchange = 0;
	unsigned int mult = 0;
	hkl_list m_hkl_list;

	/** @brief Finds the debug flag, removes it from argc/argv and stores the rest internally. */
	void look_for_debug(int& argc, char** argv);
	/** @brief Digests the options; call only after look_for_debug(). */
	//file, format and conversion options
	bool digest_io_options(const std::string &temp, int &i);
	//resources, accuracy and general run control
	bool digest_run_options(const std::string &temp, int &i);
	//partitioning schemes, scattering factors and tsc tables
	bool digest_partition_options(const std::string &temp, int &i);
	//cube and property evaluation
	bool digest_property_options(const std::string &temp, int &i);
	//RI fitting, featomic descriptors and atom classification
	bool digest_ri_options(const std::string &temp, int &i);
	//X-ray constrained wavefunction fitting
	bool digest_xcw_options(const std::string &temp, int &i);
	//development and test-only switches
	bool digest_dev_options(const std::string &temp, int &i);
	void digest_options();
	//the wavefunction of a flag written without its own file: -wfn before it, else -wfn later on the line
	std::filesystem::path line_wfn(const std::string &flag, const int i);
	//-npa x.gbw: takes an existing wavefunction file right after a flag that otherwise reads -wfn
	void flag_wavefunction(int &i);
	void loadXCWsettings();
	/** @brief The error for a requested analysis with no wavefunction to run on, "" when runnable;
	 *  without -wfn/-occ run_app_impl skips its wavefunction branch silently. */
	std::string unrunnable_analysis() const;
	/** @brief Refuses -rgbi/-npa with an analysis that ends the run before them, and an -nbo_/-nrt_
	 *  option on a line with no NBO analysis. Called last in digest_options(), so order does not matter. */
	void refuse_unread_bonding_options();

	options() : log_file(std::cout)
	{
		groups.resize(1);
	};
	options(int& argc, char** argv, std::ostream& log) : log_file(log)
	{
		groups.resize(1);
		look_for_debug(argc, argv);
	};
};

/** @brief The analysis whose options begin with this flag's prefix, nullptr if none;
 *  digest_options refuses an unclaimed flag in such a family. */
const char *owning_analysis(const std::string &flag);

/** @brief The -nbo_/-nrt_ options the NBO handlers read from their own tokens, not through a digester.
 *  An option added to a handler's scan loop needs an entry here (CliRefusalTests asserts it). */
const std::set<std::string> &nbo_family_suboptions();

void convert_tonto_XCW_lambda_steps(const std::string& str, const std::string& lambda_step, bool debug, options& opt);
