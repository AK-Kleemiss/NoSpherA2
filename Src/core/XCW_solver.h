#pragma once
#include "SCF_solver.h"

// The crystallographic quantities of the current density, as the SCF prints and descends them
struct fit_quality {
	double criterion, GooF2, R1, criterion_all, R1_all;
};

// What the SCF needs from the structure factors, implemented by XCW_wrapper so that the SCF never
// sees the structure factors, the grid or the reflection data
class SCF_coupling {

public:

	virtual ~SCF_coupling() = default;
	// F_calc, scale and the crystallographic refinement criteria of the effective density
	virtual void evaluate_cryst_crit(const dMatrix2& dm_eff) = 0;
	// The criteria evaluate_cryst_crit found
	virtual fit_quality cryst_crit() const = 0;
	// The perturbation matrix, the derivative of the criterion at the last evaluate_cryst_crit, nmo rows per spin block
	virtual void perturbation_matrix(occ::Mat& perturb, const bool unrestricted) = 0;
	// lambda times the change of the perturbation matrix for a change of the effective density,
	// the scale moving along; leaves F_calc and the scale as they were
	virtual occ::Mat perturbation_matrix_response(const dMatrix2& dD_eff, const double lambda) = 0;
	// F_calc and scale as they are now, and back
	virtual void save_F_calc() = 0;
	virtual void restore_F_calc() = 0;
	// Whether the grid blocks are on the device, for the timing record
	virtual bool grid_on_device() const = 0;
};

// The SCF of one lambda step of the X-ray constrained wavefunction fit: it minimises
// the target quantity E + lambda criterion (criterion = GooF^2 against F or F^2), with the criterion,
// its derivative (the perturbation matrix) and the response of that
// derivative taken from the coupling
class XCW_solver : public SCF_solver {

public:

	XCW_solver() = default;
	XCW_solver(const occ::core::Molecule& mol, const occ::qm::AOBasis& basis, const options::SCF_settings& settings, XCW_log_writer& writer, SCF_coupling& coupling);

	// The weight of the criterion, set before every run
	double lambda = 0.0;

protected:

	void evaluate_model_term(const scf_t& scf, const bool timed) override;
	double model_energy() const override;
	void add_model_fock(occ::Mat& F) const override;
	void add_model_response(const scf_t& scf, const occ::Mat& dD, occ::Mat& dF) override;
	bool model_active() const override { return lambda != 0; }
	void save_model_state() override;
	void restore_model_state() override;
	void begin_log() override { xcw_writer_->start_lambda(lambda); }
	void log_iteration(const scf_t& scf) override;
	void end_log(const scf_t& scf) override;
	std::string functional_name() const override { return "Target quantity"; }

private:

	SCF_coupling* coupling_ = nullptr;
	// The same writer as writer_, with the crystallographic columns
	XCW_log_writer* xcw_writer_ = nullptr;
	// The perturbation matrix and the fit quality at the density of the last evaluate_model_term, and their saved copies
	occ::Mat perturbation_, saved_perturbation_;
	fit_quality fit_{}, saved_fit_{};
	// The density matrix the structure factors are calculated from, the spin blocks summed
	void effective_dm(const scf_t& scf, const occ::Mat& D, dMatrix2& dm_eff) const;
};
