#include "pch.h"
#include "XCW_solver.h"

XCW_solver::XCW_solver(const occ::core::Molecule& mol, const occ::qm::AOBasis& basis, const options::SCF_settings& settings_in, SCF_log_writer& writer, SCF_coupling& coupling)
	: SCF_solver(mol, basis, settings_in, writer), coupling_(&coupling) {
	label_ = "XCW";
	step_name_ = "lambda step";
}

void XCW_solver::effective_dm(const scf_t& scf, const occ::Mat& D, dMatrix2& dm_eff) const {
	if (scf.ctx.mo.kind == occ::qm::SpinorbitalKind::Unrestricted) {
		for (int i = 0; i < nmo_; i++) {
			dm_eff(i, i) = D(i, i);
			dm_eff(i, i) += D(i + nmo_, i);
			for (int j = i + 1; j < nmo_; j++) {
				dm_eff(i, j) = D(i, j);
				dm_eff(i, j) += D(i + nmo_, j);
				dm_eff(j, i) = dm_eff(i, j);
			}
		}
	}
	else {
		for (int i = 0; i < nmo_; i++) {
			dm_eff(i, i) = D(i, i);
			for (int j = i + 1; j < nmo_; j++) {
				dm_eff(i, j) = D(i, j);
				dm_eff(j, i) = dm_eff(i, j);
			}
		}
	}
}

//Structure factors, scale and criteria of the current density and the perturbation matrix, the
//derivative of chi^2 with respect to the density
void XCW_solver::evaluate_model(const scf_t& scf, const bool timed) {
	const _time_point t0 = get_time();
	dMatrix2 dm_eff(nmo_, nmo_);
	effective_dm(scf, scf.ctx.mo.D, dm_eff);
	coupling_->evaluate(dm_eff);
	coupling_->perturbation(perturbation_, scf.ctx.mo.kind == occ::qm::SpinorbitalKind::Unrestricted);
	fit_ = coupling_->quality();
	if (timed) throughput::record_time("XCW structure factors + perturbation", coupling_->on_device(), get_msec(t0, get_time()));
}

//perturbation_ is the derivative of lambda * criterion^2 (the chi^2 of Jayatilaka's functional),
//so that, not the printed lambda * criterion, is what the SCF minimises
double XCW_solver::model_energy() const {
	return lambda * fit_.criterion * fit_.criterion;
}

void XCW_solver::add_model_fock(occ::Mat& F) const {
	F += perturbation_ * lambda;
}

//The change of the perturbation matrix for the density change dD, the same for both spin blocks
void XCW_solver::add_model_response(const scf_t& scf, const occ::Mat& dD, occ::Mat& dF) {
	if (lambda == 0) return;
	dMatrix2 dD_eff(nmo_, nmo_);
	effective_dm(scf, dD, dD_eff);
	const occ::Mat dP = coupling_->perturbation_response(dD_eff, lambda);
	const int spin_blocks = scf.ctx.mo.kind == occ::qm::SpinorbitalKind::Unrestricted ? 2 : 1;
	for (int b = 0; b < spin_blocks; b++) dF.middleRows(static_cast<Eigen::Index>(b) * nmo_, nmo_) += dP;
}

void XCW_solver::log_iteration(const scf_t& scf) {
	const double E = scf.ctx.energy.at("total"), penalty = fit_.criterion * lambda;
	quant_ = E + penalty;
	writer_->v = { scf.iter, lambda, fit_.criterion, fit_.GooF2, fit_.R1, E, penalty, quant_, fit_.criterion_all, fit_.R1_all, 0.0, false };
	writer_->iteration_line();
}

//What XCW_wrapper prints in the lambda table
void XCW_solver::end_log(const scf_t& scf) {
	writer_->v.lambda = lambda;
	writer_->v.E_total = scf.ctx.energy.at("total");
	writer_->v.quant = quant_;
}
