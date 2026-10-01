#pragma once
#include "convenience.h"
#include "cell.h"
#include "extinction.h"
#include "scattering_factors.h"
#include "wfn_class.h"
#include "i_tensor_stream.h"
#include <occ/qm/hf.h>
#include <thread>

class structure_factors {

public:
	structure_factors() = default;
	structure_factors(options& opt_in);

	options* opt;

//private:

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

	// The I tensor: held resident, or written to disk and read back a window
	// of reflections at a time. See decide_i_storage.
	struct I_tensor {
		// Held resident only while it fits i_tensor_max_mb; otherwise empty
		// and i_file_ carries the tensor. Read through i_block(r) either way.
		cvec I;
		//A tensor built in single precision is held, streamed and saved in single: half the
		//memory and half the traffic of the two walks per iteration, and the values are the
		//ones the GEMM produced either way. i_float / i_double in the settings override the
		//choice the build precision makes.
		std::vector<std::complex<float>> I32;
		bool i_float_ = false;
		//A copy of the resident tensor on the device does both SCF walks there
		bool i_on_device_ = false;
		i_tensor_file i_file_;
		bool i_streamed_ = false;
		int i_window_ = 0;
		//Only the (mu, nu) pairs the overlap screening kept are stored, i_compact_ of the
		//nmo (nmo + 1) / 2, in the order these two lists give
		ivec i_pair_mu_, i_pair_nu_;
		size_t i_compact_ = 0;
		// The packed (mu, nu) run of reflection r, from the loaded window or from the resident tensor
		const cdouble* i_block(const int r) const
		{
			return i_streamed_ ? i_file_.block(r) : I.data() + static_cast<size_t>(r) * i_compact_;
		}
		const std::complex<float>* i_block32(const int r) const
		{
			return i_streamed_ ? i_file_.block32(r) : I32.data() + static_cast<size_t>(r) * i_compact_;
		}
		bool read_tensor;
		std::filesystem::path i_tensor_file_path;
		std::filesystem::path i_tensor_save_path;
		bool tensor_single;
		bool tensor_double;
		size_t i_tensor_max_mb;
		std::string basis_set_name;
	};

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
	scatter_data scatter_data;
	// Store cristallographic quality criteria
	quality_criteria quality_criteria;
	// Store the unit cell
	cell unit_cell;
	// Store the atoms of the asymmetric unit (or the grown unit)
	std::vector<asym_atom> asym_atoms;
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
	// The I tensor, see eval_I_anom_disp
	I_tensor I_tens;
	// The background writer for `save <path>`. Joined, never detached: a thread still
	// running at exit is how the GPU warm-up bug of 939268f happened, and this one holds a
	// FILE* and reads the resident tensor.
	std::thread i_writer_;
	std::string i_writer_error_;
	static constexpr const char* i_tensor_default = "I_tensor_stream.bin";


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
	// index r used for F_obs[r]/F_calc[r] (see generate_asym_lookup)
	void ensure_hkl_ordered();
	// Builds (once) the per-reflection 1/|H|^2 cache used by calc_criteria and the
	// perturbation when h2 weighting is set. No-op otherwise.
	void ensure_inv_H2_weights();

	//// Converts the ADP matrix (just U) from cif format into reciprocal space
	void U_cif2U_star();
	// Converts all ADP tensors from reciprocal space into real space
	void U_star2U_cart();

	// Generates a list that links the symmetry operations to symmetry-generated reflexes for given reflex r
	ivec generate_asym_lookup(const int r);

	// Evaluates Debye-Waller factors
	void eval_DW();
	// Evaluates the rotational contribution to the phase factors
	void eval_phase();
	// Evaluates the translational contribution to the phase factors
	void eval_translation_phase();

	// Calculates direct corrections of the anomalous dispersion onto F_calc
	void eval_anom_disp();

	// Creates primitive vectors from the basis set for calculating the XCW integrals
	void create_prims(std::vector<ao_data>& ao_data_shells, occ::qm::AOBasis& occ_basis_set);
	// Helper function for flattening the I tensor
	size_t tri_index(int mu, int nu) const noexcept;
	// Combined method used to save memory, calculates (or reads) the I tensor and the correction for F_calc from anomalous dispersion
	I_tensor& eval_I_anom_disp(std::vector<ao_data>& ao_data_shells);
	// Evaluates the I tensor
	void eval_I(std::vector<ao_data>& ao_data_shells, double& time_taken, long long& screen_counter, long long& skipped_grids);
	void decide_i_storage();
	void open_i_stream_for_reading();
	std::filesystem::path i_tensor_path() const;
	void start_i_save();
	void finish_i_save();
	size_t i_budget(const char*& source, bool& automatic) const;

	// Calculates F_calc from the atomic scattering factors of every atom over hkl_enlarged
	void calc_F_calc(const cvec2& scattering_factors);
	// Calculates F_calc from the I tensor and a density matrix
	void calc_F_calc(I_tensor& I, const dMatrix2& D);

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
	// The criterion the SCF descends (XWR_type x refine_against), over the fit set or over all
	double criterion(bool all) const;

	extinction_settings extinction_settings;


	// closing class
};
