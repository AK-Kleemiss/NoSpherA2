#pragma once
#include "convenience.h"
#include "cell.h"
#include "extinction.h"
#include "scattering_factors.h"
#include "wfn_class.h"
#include "GridManager.h"
#include <occ/qm/hf.h>

class structure_factors {

public:
	structure_factors() = default;
	structure_factors(options& opt_in);

	options* opt;

	// Store cristallographic quality criteria
	struct quality_criteria {
		double GooF1;
		double GooF2;
		double weighted_GooF1;
		double weighted_GooF2;
		// R1 = sum ||F_obs| - scale |F_calc|| / sum |F_obs|, over the fit set (Tonto's R1(gt))
		double R1;
		// The same numbers over every reflection, next to the fit-set ones above
		double R1_all;
		double GooF1_all;
		double GooF2_all;
		double weighted_GooF1_all;
		double weighted_GooF2_all;
		// wR2 = sqrt(Sum w (Fo^2 - Fc^2)^2 / Sum w Fo^4) with the crystallographic weight, fit set and all
		double wR2;
		double wR2_all;
		int refine_against; // 1 for F and 2 for F^2
		int goof_type; // 1 for traditional GooF and 2 for Coulomb weighted GooF
	};

	// Store information about the model (e.g. number of atoms, reflections)
	struct model_data {
		int ncen;
		int nr;
		int nr_enlarged;
		int nr_fit;
		int nmo;
		int n_params;
		double wavelength;
	};

	struct extinction_settings {
		extinction::model extinction_model;
		bool aniso;
		bool refine;
		double start;
	};

	// Data for contracted basis function
	struct ao_data {
		std::vector<primitive> prims;
		d3 pos;
		int m;
	};

	//A run of points of one atom grid with the AOs not negligible anywhere on it, ascending,
	//and their values row-major, one row of point_count per AO
	struct GridBlock {
		int point_start;
		int point_count;
		ivec active_aos;
		vec ao_values;
		std::vector<float> ao_values_f;
	};
	//XCW's integration grid, set up by setup_grid. An evaluation integrates rho = 2 sum_munu D chi_mu chi_nu block by
	//block and transforms it; the perturbation is the transpose of that chain. AO pairs in skip and AOs below the block
	//threshold count nowhere. Memory is the AO values of the blocks, in single precision when fp32
	struct xcw_grid {
		bool fp32 = true;
		std::vector<std::vector<GridBlock>> blocks;
		vec2 d1, d2, d3, w;
		ivec2 skip, lookup, inversion_link;
		ivec clean_of;
		bvec is_inverse;
		std::shared_ptr<struct xcw_fb> fb; //the planned Fourier-Bessel transform of every atom grid, null: the direct sum
		//grid_F's change-only update: a block whose rho can have moved by less than NOSPHERA2_XCW_GRID_DTOL
		//since it was last integrated keeps its rho. bound = Sum_p w_p (Sum_mu |chi_mu|)^2 per block, so
		//max|dD| bound bounds Sum_p |w drho| and with it |dF|; drift sums that over the calls it was kept
		occ::Mat D_last;
		vec2 rho_last;
		vec bound, drift;
		//the block stages run on the device (xcw_grid_gpu.h), atom grid g's points from gpu_off[g] there
		bool on_device = false;
		ivec gpu_off;
		//TFVC, MBIS and EMBIS shares follow rho (xcw_weights.cpp): w0 holds the integration weights alone and w = w0 s for the
		//current share s, which grid_rho recomputes from its own rho when refresh is set; rho_s is then its partitioned rho
		vec2 w0, rho_s;
		bool refresh = false;
		double tol = 0.0; //the MBIS/EMBIS refit's criterion when set, else NOSPHERA2_XCW_WEIGHT_TOL
		std::shared_ptr<struct xcw_weights> ws;
	} XG;
	//setup_grid's grid as TFVC, MBIS or EMBIS partition it; the Hirshfeld shares it was built with stand until rho is known
	void setup_grid_weights(const occ::qm::AOBasis& ao);
	//the scheme's shares from rho = w0 rho(p) of D, each fit starting from the last one
	void update_grid_weights(const vec2& rho, const dMatrix2& D);
	//the fits since the last call in one line, empty for the fixed shares
	std::string grid_weight_report();
	//Uploads the blocks for xcw_grid_gpu; false keeps them on the host
	bool hold_grid_gpu();
	//Plans XG.fb from setup_grid's layout of the grid points in their prototypes, see xcw_grid.cpp
	void plan_grid_fb(const struct sf_shell_layout& layout);
	//F_calc = anomalous + the transform of rho(D)
	void grid_F(const dMatrix2& D);
	//the per-atom transform of rho(D) at the cleaned entries, [g * nc + c]; grid_F's first half and the tscb's table.
	//The atoms' shares come from the grid weights, or from the RI fit when XR.on
	cvec grid_sfc(const dMatrix2& D);
	//out = Sum_r Re(pre_r dF_r/dD(mu, nu)) over the fit set, the transpose of grid_F: the walk the perturbation and its
	//response share; nmo x nmo, symmetric, without any prefactor
	void grid_contract(occ::Mat& out, const cvec& pre);
	//the weighted density 2 w_p rho(p) of every atom grid, grid_sfc's first half
	vec2& grid_rho(const dMatrix2& D);
	//out += Sum_p chi_mu V_p chi_nu over the blocks, the transpose of grid_rho for V_p = w_p v_p
	void grid_apply_V(const vec2& V, occ::Mat& out);

	//The atoms' shares of rho as an RI fit (-xcw_partition ri) instead of the grid weights: c = H^-1 (2 B D + R^T S^2 T(D)),
	//B = (mu nu|P) from occ's density-fitting engine, H = J + R^T S^2 R, R the analytic multipoles of the aux functions about
	//their atom and T(D) those of the grid partition (Hirshfeld, Becke, or TFVC, MBIS and EMBIS shares as they follow rho, fixed
	//within an evaluation), linear in D; f_g(k) = Sum_{P on g} c_P F_P(k) in closed form. ri_contract runs the chain backwards,
	//so the SCF gets the exact derivative at fixed shares. Nuclear and Sanderson charges are fixed targets (set_ri_targets,
	//frozen: left out of a response like the anomalous part), Mulliken's 2 Sum_{mu on a} (D S)_mumu linear in D as the grid's
	struct xcw_ri {
		bool on = false, frozen = false, response = false, mulliken = false, packed = false;
		int nbf = 0, naux = 0, lmax = -1, n_mom = 0;
		occ::Mat B;
		Eigen::LLT<occ::Mat> H;
		occ::Mat R, S;
		vec s2;
		ivec atom_start, bf_atom;
		std::vector<cvec> F;
		vec2 M;
		vec c, t0;
		std::vector<occ::Mat> A;
		occ::qm::AOBasis aux; //the fit's basis, pure, in occ's order and normalisation: XR.c's
		std::shared_ptr<BasisSet> aux_set; //the same basis as BasisSet, for a .tshb on it
	} XR;
	//moments [atom][l*l+l+m] in constants::spherical_harmonic's convention (DensityFitting's) as the fixed targets
	void set_ri_targets(const vec2& moments);
	void setup_ri(const occ::qm::AOBasis& ao, const occ::qm::AOBasis& aux, const int lmax, const double strength, const MultipoleScheme scheme);
	//the same in the -ri_fit basis, AutoAux on the atoms of ao without one, restrained as -multipole_moments says
	void setup_ri(const occ::qm::AOBasis& ao);
	//Re Sum_r pre_r (F_r(D) - anomalous_r) against 2 <D, P(pre)> at a random D and pre, printed: the perturbation has
	//to be the transpose of F_calc, whichever partition runs, with the shares the grid holds
	void check_adjoint();
	cvec ri_sfc(const dMatrix2& D);
	void ri_contract(const cvec& H, occ::Mat& out);

	//The entries of a rank-N symmetric tensor stored packed in the order i <= j <= k (<= l): the powers of x, y
	//and z each one multiplies and how often its index tuple occurs among all 3^N
	struct gc_term { int e[3]; double mult; };

	// Stores the ADP tensors
	vec3 ADPs;
	// Stores the Debye-Waller factors
	cvec2 DW_facts;
	// Stores the rotational phase factors
	cvec2 phase_facts;
	// Stores the translational phase factors
	cvec2 translation_phase_facts;
	// Store information about the model
	model_data model_data;
	// Store scattering data for each reflection
	::scatter_data scatter_data;
	// Store cristallographic quality criteria
	quality_criteria quality_criteria;
	// Store the unit cell
	cell unit_cell;
	// Store the atoms of the asymmetric unit (or the grown unit)
	std::vector<asym_atom> asym_atoms;
	// The asymmetric unit, the grown cluster with the operation behind each grown atom, the kept operations
	// and the sym_facts, for the SCF.log before the first cycle
	std::string model_summary;
	// The orbital basis, the library one or its per-element overrides
	std::shared_ptr<BasisSet> orbital_basis;
	// Store the k points
	vec2 k_pt;
	// The refined extinction coefficient, or the six Voigt components X11 X22 X33 X12 X13 X23
	// of the anisotropic tensor. Empty when no model is active, which is what every extinction
	// branch tests on.
	vec ext_p_;
	// Per-reflection geometry, built once: 0.001 lambda^3/sin(2 theta), cos(2 theta), and (for
	// the anisotropic models only) the nr x 6 coefficients a_{r,p}
	vec ext_c_, ext_cos2t_, ext_a_;
	// Per-iteration shape, see update_extinction. ext_dyc_ is dy/dt * c_r * |Fc_r|^2, so that
	// dy/dP_p = ext_dyc_[r] * a_{r,p}.
	vec ext_y_, ext_sqrt_y_, ext_g_, ext_m_, ext_dyc_;
	// Ordered vector of hkl
	std::vector<i3> hkl_ordered_;
	// 1/|H_r|^2 per reflection, see ensure_inv_H2_weights.
	vec inv_H2_;
	// The crystallographic weight of every fit sum (scale, extinction, criteria, perturbation), held as the
	// sigma(F) and sigma(F^2) it amounts to, sigma = 1/sqrt(w): the default scheme copies the read sigmas,
	// so every (x/sigma)^2 there stays bit for bit what it was. Filled by update_cryst_weights.
	vec sig_F_, sig_I_;
	// The scheme's current parameters (a, b of SHELX, a of stl, A_r of Chebyshev), from the settings until
	// shelx_opt / chebyshev_opt move them
	vec weight_ab_;
	// unit, stl, chebyshev: every w times weight_k_, set once at the first fill (ensure_cryst_weights) so that the
	// criterion S^2 of the base scheme is 1 there, then held for the whole scan: rescaling at every step would pin
	// S to 1 and take the stopping signal with it. Robust goes on top (its z = Delta sqrt(w) is scale-free)
	double weight_k_ = 1.0;
	bool rescale_pending_ = false;
	// where the rescale line goes besides the console (the XCW SCF.log)
	std::ostream* weight_log_ = nullptr;
	// "shelx_opt WGHT 0.1 0" / "robust(c 6) chebyshev_opt A 1 0.5", for the logs
	std::string weight_report() const;
	// Fills sig_F_ and sig_I_ from the scheme, the current F_calc and scale (so only after calc_F_calc).
	// The weights then stay fixed until the next call: no derivative wrt F_calc. The 1/|H|^2 of XWR_type 2
	// still multiplies on top at every site, so each scheme works with both criteria
	void update_cryst_weights();
	// The first fill, from eval_scale or calc_criteria: the read sigmas, and for another scheme its weights
	// at the scale of that 1/sigma^2 fit and the F_calc of the moment. No-op once filled
	void ensure_cryst_weights();
	// shelx_opt: a and b re-optimised on the fit set as cctbx's new_shelx_weighting.optimise_parameters
	// does, against the current F_calc and scale; chebyshev_opt: A_r fit to Delta^2 (fit_chebyshev). Refills the weights
	void optimise_cryst_weights();
	double fit_residual(int r) const;
	vec chebyshev_t() const;


	// The parameters the criteria divide by: the settings file's `params` plus the extinction
	// coefficients, but only while those are actually being refined
	int n_params() const {
		return model_data.n_params + static_cast<int>(extinction_settings.refine ? ext_p_.size() : 0);
	}
	// Reads the wavelength (settings file, else _diffrn_radiation_wavelength in the CIF),
	// sizes the coefficient vector and builds the per-reflection extinction geometry. No-op
	// when no extinction model was asked for. Called once from the constructor.
	void setup_extinction(const std::filesystem::path& cif);
	// Recomputes y_r, sqrt(y_r), dI/d|Fc|^2 and d|Fc_ext|/d|Fc| from the current F_calc and
	// the current coefficients. No-op when no model is active.
	void update_extinction();
	// One Levenberg-damped Gauss-Newton step on the coefficients against the same weighted
	// residual eval_scale minimises for the scale, with the scale held at its current value.
	// False when the step was negligible or had to be rejected.
	bool refine_extinction_step();
	// "ext 0.000123" or the six tensor components, for the per-lambda log
	std::string extinction_report() const;
	// The extinction shape of reflection r, 1 where no model is active
	double ext_y(const int r) const { return ext_y_.empty() ? 1.0 : ext_y_[r]; }
	double ext_sqrt_y(const int r) const { return ext_sqrt_y_.empty() ? 1.0 : ext_sqrt_y_[r]; }
	// dI/d|Fc|^2 and d(sqrt(y)|Fc|)/d|Fc|, the chain factors the gradient needs
	double ext_g(const int r) const { return ext_g_.empty() ? 1.0 : ext_g_[r]; }
	double ext_m(const int r) const { return ext_m_.empty() ? 1.0 : ext_m_[r]; }
	// a_{r,p} of x_r = sum_p a_{r,p} P_p; 1 for the isotropic models, which have one P
	double ext_a(const int r, const size_t p) const {
		return ext_a_.empty() ? 1.0 : ext_a_[static_cast<size_t>(r) * ext_p_.size() + p];
	}
	double ext_x(const int r) const {
		double x = 0.0;
		for (size_t p = 0; p < ext_p_.size(); p++) x += ext_a(r, p) * ext_p_[p];
		return x;
	}

	// Builds the (once-cached) ordered list of Miller indices matching the
	// index r used for F_obs[r]/F_calc[r] (see cell::generate_asym_lookup)
	void ensure_hkl_ordered();
	// Builds (once) the per-reflection 1/|H|^2 cache used by calc_criteria and the
	// perturbation when h2 weighting is set. No-op otherwise.
	void ensure_inv_H2_weights();

	//// Converts the ADP matrix (just U) from cif format into reciprocal space
	void U_cif2U_star();
	// Converts all ADP tensors from reciprocal space into real space
	void U_star2U_cart();

	// Each of the three setters computes its quantity and stores it; one never called leaves the
	// value the constructor put there (DW and phases 1, anomalous dispersion 0)
	// Evaluates and sets the Debye-Waller factors
	void set_DW();
	// Evaluates and sets the rotational and the translational phase factors
	void set_phases();
	// Evaluates and sets the direct corrections of the anomalous dispersion onto F_calc, from
	// whatever DW and phase factors are set at the time
	void set_anom();
	bool DW_set_ = false, phases_set_ = false, anom_set_ = false;

	// Creates primitive vectors from the basis set for calculating the XCW integrals
	void create_prims(std::vector<ao_data>& ao_data_shells, const occ::qm::AOBasis& occ_basis_set);
	// XCW's integration grid (XG) with the DW and phase factors as set: the symmetry links of the reflections, the
	// AO pair screen, the atom grids and the AO values on them in blocks, and the shares of the partition scheme
	void setup_grid(const occ::gto::AOBasis& aobasis);
	// The k-points of every reflection's symmetry images (XG.lookup) and which of them inversion makes conjugates
	// of others, so that only one of each pair is integrated (XG.inversion_link, clean_of, is_inverse)
	void link_symmetry_equivalents();
	// XG.skip: the AO pairs on different centres whose overlap bound is negligible
	void screen_ao_pairs(const std::vector<ao_data>& ao_data_shells);
	// The atom grids of the asymmetric unit in grid_manager, pruned on the promolecule, as XG.d1..d3 and XG.w, the dummy
	// wavefunction they were built for and where each point sits in its prototype grid
	void setup_atom_grids(GridManager& grid_manager, WFN& dummy_wave, struct sf_shell_layout& layout);
	// The AO values on the atom grids, Morton-ordered and cut into blocks that keep only the AOs not negligible
	// on them (XG.blocks)
	void build_grid_blocks(const std::vector<ao_data>& ao_data_shells, GridManager& grid_manager, const WFN& dummy_wave, struct sf_shell_layout& layout);
	size_t free_memory() const;
	size_t df_reserved_ = 0;

	// F_calc from the atomic scattering factors of every atom over hkl_enlarged
	void calc_F_calc(const cvec2& scattering_factors);
	// F_calc from a density matrix, on XCW's grid
	void calc_F_calc(const dMatrix2& D);
	// F_calc from the scattering factors of a tsc table
	void calc_F_calc(const itsc_block& tsc);
	// F_calc from wavefunctions, through the tsc table tsc_table computes for them
	void calc_F_calc(std::vector<WFN>& wavefunctions, GridManager* grids, std::ostream& log);
	// The tsc table of the wavefunctions at every entry of hkl_enlarged; grids keeps the integration grids for the next call
	itsc_block tsc_table(std::vector<WFN>& wavefunctions, GridManager* grids, std::ostream& log);

	// The perturbation matrix: the derivative of the fit criterion with respect to the density matrix, nmo rows per spin block
	void perturbation(occ::Mat& perturb, const bool unrestricted, std::ostream& log);
	// lambda times the change of the perturbation matrix for a change dD_eff of the density, the scale
	// moving along; leaves F_calc and the scale as they were
	occ::Mat perturbation_response(const dMatrix2& dD_eff, const double lambda);
	// The F_calc table of the current density, one reflection per line
	void write_F_calc(const std::string& path);
	// F_calc and scale as they are now, and back
	void save_F_calc();
	void restore_F_calc();
	cvec saved_F_calc_;
	double saved_scale_ = 0.0;

	// The tsc table of a converged density for the tscb output, see result_table in structure_factors.cpp;
	// refit_shares fits density-following shares to it first
	itsc_block result_table(const occ::Mat& D_wfn, std::vector<WFN>& wavefunctions, std::ostream& log, const bool refit_shares);
	// What result_table keeps between steps: the integration grids of the tsc route and, on the grid route,
	// the layout of the first table: its scatterers, hkl rows, the asym atom of each scatterer and the
	// hkl_enlarged entry of each row
	GridManager tsc_grids;
	ScattererLabels tscb_scat_;
	ivec2 tscb_hkl_;
	ivec tscb_atom_, tscb_entry_;

	// Evaluates the scaling factor for |F_calc| by least squares fitting, and with it the
	// extinction coefficients when a model is being refined
	void eval_scale();
	// The closed-form weighted least-squares scale alone, with the extinction shape as it is
	void solve_scale();
	// Calculates quality criteria like GooF and chi^2. When
	// h2 weighting is set, both are computed with an additional
	// 1/|H|^2 weighting (XCW_plan.md sec. 6.2, residual self-energy
	// criterion) instead of the traditional unweighted sums.
	void calc_criteria();
	// The criterion GooF^2 against F or F^2 (XWR_type x refine_against), over the fit set or over all
	double criterion(bool all) const;

	extinction_settings extinction_settings;


	// closing class
};

// SHELXL's weight on Fo^2 as cctbx's mainstream/new_shelx_weighting computes it, with Fc^2 unscaled and
// K the scale on it (Fo^2 ~ K Fc^2): P = (max(Fo^2, 0) + 2 K Fc^2)/3, w = 1/(sigma^2 + (a P)^2 + b K P)
double shelx_weight(double fo_sq, double sigma, double fc_sq, double K, double a, double b);
// cctbx new_shelx_weighting.optimise_parameters: a, b from the start values of SHELXL's formulae that
// make the GooF of 10 bins by Fc^2/max(Fc^2) closest to 1; {a, b}
vec optimise_shelx_weights(const vec& fo_sq, const vec& sigma, const vec& fc_sq, double K, int n_params);
// Carruthers & Watkin (1979): w = 1/Sum_r A_r T_r(t), 0 where the sum is not positive
double chebyshev_weight(double t, const vec& A);
// chebyshev_opt: the A_0..A_{n-1} whose Sum A_r T_r(t_i) fits d2_i (the squared residuals) in least squares
vec fit_chebyshev(const vec& t, const vec& d2, int n);
// Tukey's biweight on z_i = Delta_i sqrt(w_i): f_i = (1 - (z_i/(c s))^2)^2 for |z_i| < c s, else 0, with the robust
// scale s = median(|z|)/0.6745 over the reflections where in_scale is non-zero
vec tukey_biweight(const vec& z, const ivec& in_scale, double c);
