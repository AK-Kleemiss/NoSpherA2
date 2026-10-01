#include "pch.h"
#include "SCF_wrapper.h"
#include "XCW_solver.h"
#if defined(NOSPHERA2_USE_GPU) || defined(NOSPHERA2_USE_METAL)
#include "itensor_gpu.h"
#endif
#include "basis_set.h"
#include "citations.h"
#include <limits>
#include <random>

SCF_wrapper::SCF_wrapper(XCW_solver& xcw_in, structure_factors& sf_in) {
	xcw = &xcw_in;
	sf = &sf_in;
	opt = sf->opt;
	SCF_log.open("SCF.log");
}

//C = op(A) op(B) for column-major occ matrices through MKL, Src/core's Eigen being serial
static occ::Mat gemm(const occ::Mat& A, const occ::Mat& B, const bool ta = false, const bool tb = false)
{
	const int m = static_cast<int>(ta ? A.cols() : A.rows()), k = static_cast<int>(ta ? A.rows() : A.cols()), n = static_cast<int>(tb ? B.rows() : B.cols());
	occ::Mat C(m, n);
	cblas_dgemm(CblasColMajor, ta ? CblasTrans : CblasNoTrans, tb ? CblasTrans : CblasNoTrans, m, n, k, 1.0, A.data(), static_cast<int>(A.rows()), B.data(), static_cast<int>(B.rows()), 0.0, C.data(), m);
	return C;
}

void SCF_wrapper::setup_SCF_mol(occ::core::Molecule& mol) {
	double bohr2angstrom = constants::bohr2ang(1);
	std::ostringstream init_stream;

	init_stream << sf->model_data.ncen << "\n\n";

	for (int i = 0; i < sf->model_data.ncen; ++i) {
		init_stream
			<< constants::atnr2letter(sf->asym_atoms[i].type) << " "
			<< sf->asym_atoms[i].pos[0] * bohr2angstrom << " "
			<< sf->asym_atoms[i].pos[1] * bohr2angstrom << " "
			<< sf->asym_atoms[i].pos[2] * bohr2angstrom;

		if (i != sf->model_data.ncen - 1)
			init_stream << "\n";
	}

	std::string init = init_stream.str();
	mol = occ::io::molecule_from_xyz_string(init);
	mol.set_charge(opt->xcw_settings.charge);
	mol.set_multiplicity(opt->xcw_settings.multiplicity);
}

//The same walk as OCC's SOAD guess with a converged Hartree-Fock density in place of the
//tabulated atoms: F = H + G(D_small) through the mixed-basis build, the orbitals from its
//diagonalisation. The spin blocks are summed, the first iteration polarises them again.
void SCF_wrapper::small_basis_guess(occ::qm::SCF<occ::qm::HartreeFock>& scf) {
	occ::core::Molecule mol;
	setup_SCF_mol(mol);
	occ::qm::AOBasis small_bs;
	setup_basis(mol, opt->xcw_settings.guess_basis_name, small_bs);
	occ::qm::HartreeFock hf_small(small_bs);
	occ::qm::SCF scf_small(hf_small, opt->xcw_settings.hf_type);
	scf_small.set_charge_multiplicity(opt->xcw_settings.charge, opt->xcw_settings.multiplicity);
	scf_small.maxiter = opt->xcw_settings.max_scf_iterations;
	const double e_small = scf_small.compute_scf_energy();
	SCF_log << "XCW: initial guess from a " << opt->xcw_settings.guess_basis_name << " Hartree-Fock (" << small_bs.nbf()
		<< " functions), E = " << std::fixed << std::setprecision(8) << e_small << " Eh" << std::endl;

	occ::qm::MolecularOrbitals guess;
	guess.kind = scf.ctx.mo.kind;
	guess.n_ao = static_cast<int>(small_bs.nbf());
	guess.n_alpha = scf.n_alpha();
	guess.n_beta = scf.n_beta();
	guess.D = opt->xcw_settings.hf_type == occ::qm::SpinorbitalKind::Unrestricted
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

void SCF_wrapper::setup_basis(occ::core::Molecule& mol, std::string& basis_set_name, occ::qm::AOBasis& occ_basis_set) {
	std::shared_ptr<BasisSet> basis_set = BasisSetLibrary::get_basis_set(basis_set_name);
	occ_basis_set = basis_set->to_AOBasis(mol.atoms());
}

double SCF_wrapper::dynamic_damping(const occ::qm::SCF<occ::qm::HartreeFock>& scf, const double& current_alpha, const double& quant_diff, double& quant_diff_mem) {
	double new_alpha = current_alpha;
	if (quant_diff < quant_diff_mem / 10) {
		new_alpha *= 0.75;
		quant_diff_mem = quant_diff;
		if (quant_diff < 10 * scf.convergence_settings.energy_threshold) {
			print_centered_message("***Turned off damping***", 84, SCF_log);
			new_alpha = 0;
			opt->xcw_settings.apply_damping = false;
		}
		else {
			std::stringstream print_;
			print_ << "***Decreased damping to " << std::fixed << std::setprecision(3) << new_alpha << "***";
			print_centered_message(print_.str(), 84, SCF_log);
		}
	}
	return new_alpha;
	// closing function
}

void SCF_wrapper::apply_level_shift(const occ::Mat& C_old, const occ::qm::SCF<occ::qm::HartreeFock>& scf, occ::Mat& F_diis) {
	const int nocc = scf.ctx.mo.Cocc.cols();
	if (scf.ctx.mo.kind == occ::qm::SpinorbitalKind::Restricted) {
		const occ::Mat SC_virt = scf.ctx.S * C_old.rightCols(sf->model_data.nmo - nocc);
		F_diis.noalias() += opt->xcw_settings.level_shift * SC_virt * SC_virt.transpose();
	}
	else {
		const int nao = C_old.rows() / 2;
		const auto S_ao = scf.ctx.S.topRows(nao);
		//Cocc has max(n_alpha, n_beta) columns, the spin counts differ for open shells
		const occ::Mat SC_virt_a = S_ao * C_old.topRows(nao).rightCols(sf->model_data.nmo - scf.ctx.mo.n_alpha);
		const occ::Mat SC_virt_b = S_ao * C_old.bottomRows(nao).rightCols(sf->model_data.nmo - scf.ctx.mo.n_beta);
		F_diis.topRows(nao).noalias() += opt->xcw_settings.level_shift * SC_virt_a * SC_virt_a.transpose();
		F_diis.bottomRows(nao).noalias() += opt->xcw_settings.level_shift * SC_virt_b * SC_virt_b.transpose();
	}
}

void SCF_wrapper::build_effective_dm(const occ::qm::SCF<occ::qm::HartreeFock>& scf, dMatrix2& dm_ref, const occ::Mat& dm_old) {
	if (scf.ctx.mo.kind == occ::qm::SpinorbitalKind::Unrestricted) {
		for (int i = 0; i < sf->model_data.nmo; i++) {
			dm_ref(i, i) = dm_old(i, i);
			dm_ref(i, i) += dm_old(i + sf->model_data.nmo, i);
			for (int j = i + 1; j < sf->model_data.nmo; j++) {
				dm_ref(i, j) = dm_old(i, j);
				dm_ref(i, j) += dm_old(i + sf->model_data.nmo, j);
				dm_ref(j, i) = dm_ref(i, j);
			}
		}
	}
	else {
		for (int i = 0; i < sf->model_data.nmo; i++) {
			dm_ref(i, i) = dm_old(i, i);
			for (int j = i + 1; j < sf->model_data.nmo; j++) {
				dm_ref(i, j) = dm_old(i, j);
				dm_ref(j, i) = dm_ref(i, j);
			}
		}
	}
}

bool SCF_wrapper::do_SCF(const double& lambda, double& alpha, occ::qm::SCF<occ::qm::HartreeFock>& scf, occ::qm::Wavefunction& last_wfn, bool& has_guess, const bool write_result) {

	opt->xcw_settings.clear();
	diis_F_.clear();
	diis_E_.clear();
	adiis_.reset();
	ediis_.reset();
	best_quant_ = std::numeric_limits<double>::infinity();
	rescues_ = 0;
	soscf_reset();

	SCF_log << "Starting XCW SCF solver with lambda = " << std::fixed << std::setprecision(5) << lambda << "\n";
	SCF_log << "____________________________________________________________________________________\n";
	SCF_log << " Iteration\t\tCriterion\tGooF(F^2)\tR1(gt)\t\tTotal Energy\t\tPerturbation\tTarget quantity\n";
	SCF_log << "\t\t\t\t\t\t\t\t\t(Eh)\t\t\t(a. u.)\t\t(a. u.)\n";
	SCF_log << "____________________________________________________________________________________\n";

	// Compute first guess and update the energy according to this guess
	const _time_point guess_t0 = get_time();
	if (has_guess) {
		scf.set_initial_guess_from_wfn(last_wfn);
	}
	else {
		if (opt->xcw_settings.guess_basis_name.empty())
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
		SCF_log << "____________________________________________________________________________________\n";
		std::stringstream print_;
		print_ << "***SCF converged in " << scf.iter << " iterations***";
		print_centered_message(print_.str(), 84, SCF_log);

		//Before the summary line below so its A^2 can be appended as an extra column
		if (opt->xcw_settings.xcw_gaussian_halt) {
			xcw->evaluate_gaussian_halting(lambda);
		}


		const double current_criterion = sf->criterion(false);
		std::cout << std::fixed << std::setprecision(5) << lambda << "\t\t" << std::fixed << std::setprecision(4) << current_criterion << "\t\t" << sf->quality_criteria.GooF2 << "\t\t" << std::setprecision(5) << sf->quality_criteria.R1 << "\t\t" << std::fixed << std::setprecision(9) << scf.ctx.energy["total"] << "\t\t" << std::fixed << std::setprecision(3) << lambda * current_criterion << "\t\t" << std::fixed << std::setprecision(9) << quant
			<< "\t\t" << std::setprecision(4) << sf->criterion(true) << "\t\t" << std::setprecision(5) << sf->quality_criteria.R1_all;
		if (opt->xcw_settings.xcw_gaussian_halt && !xcw->gaussian_halt_history_.empty()) {
			std::cout << "\t\t" << std::setprecision(4) << xcw->gaussian_halt_history_.back().A2;
		}
		std::cout << std::endl;

		if (write_result) {
			const _time_point tscb_t0 = get_time();
			xcw->create_tscb(scf, lambda);
			throughput::record_time("XCW tscb", false, get_msec(tscb_t0, get_time()));
		}
	}
	else {
		SCF_log << "____________________________________________________________________________________\n";
		print_centered_message("***SCF did not converge***", 84, SCF_log);
		std::ostringstream perturbed_energy;
		perturbed_energy << " for perturbed energy: " << std::scientific << opt->xcw_settings.quant_diff << " (current: " << opt->xcw_settings.current_quant_diff << ") \n";
		std::ostringstream diis_error;
		diis_error << " for DIIS error: " << std::scientific << opt->xcw_settings.max_diis_error << " (current: " << opt->xcw_settings.current_max_diis_error << ") \n";
		std::ostringstream orbital_gradient;
		orbital_gradient << " for orbital gradient: " << std::scientific << opt->xcw_settings.gradient << " (current: " << opt->xcw_settings.current_gradient << ") \n";
		std::ostringstream max_density_diff;
		max_density_diff << " for maximum difference in density matrix: " << std::scientific << opt->xcw_settings.MaxP_diff << " (current: " << opt->xcw_settings.current_MaxP_diff << ") \n";
		std::ostringstream rmsd_density;
		rmsd_density << " for RMSD of density matrix: " << std::scientific << opt->xcw_settings.RMSP_diff << " (current: " << opt->xcw_settings.current_RMSP_diff << ") \n";
		if (opt->xcw_settings.conv_quant_diff) {
			SCF_log << "CONVERGED";
		}
		else {
			SCF_log << "NOT CONVERGED";
		}
		SCF_log << perturbed_energy.str();
		if (opt->xcw_settings.conv_max_diis_error) {
			SCF_log << "CONVERGED";
		}
		else {
			SCF_log << "NOT CONVERGED";
		}
		SCF_log << diis_error.str();
		if (opt->xcw_settings.conv_gradient) {
			SCF_log << "CONVERGED";
		}
		else {
			SCF_log << "NOT CONVERGED";
		}
		SCF_log << orbital_gradient.str();
		if (opt->xcw_settings.conv_MaxP_diff) {
			SCF_log << "CONVERGED";
		}
		else {
			SCF_log << "NOT CONVERGED";
		}
		SCF_log << max_density_diff.str();
		if (opt->xcw_settings.conv_RMSP_diff) {
			SCF_log << "CONVERGED";
		}
		else {
			SCF_log << "NOT CONVERGED";
		}
		SCF_log << rmsd_density.str();
	}
	return converged;
}

double SCF_wrapper::compute_orbital_gradient(const occ::qm::SCF<occ::qm::HartreeFock>& scf) {
	if (opt->xcw_settings.hf_type == occ::qm::SpinorbitalKind::Restricted) {
		const occ::Mat& C = scf.molecular_orbitals().C;
		const occ::Mat& Cocc = scf.molecular_orbitals().Cocc;
		const occ::Mat Cvir = C.rightCols(C.cols() - Cocc.cols());
		return 2.0 * gemm(Cvir, gemm(scf.ctx.F, Cocc), true).norm();
	}
	else if (opt->xcw_settings.hf_type == occ::qm::SpinorbitalKind::Unrestricted) {
		occ::Mat C_alpha = scf.molecular_orbitals().C.topRows(sf->model_data.nmo);
		occ::Mat C_beta = scf.molecular_orbitals().C.bottomRows(sf->model_data.nmo);
		occ::Mat Cocc_alpha = C_alpha.leftCols(scf.ctx.mo.n_alpha);
		occ::Mat Cocc_beta = C_beta.leftCols(scf.ctx.mo.n_beta);
		occ::Mat Cvir_alpha = C_alpha.rightCols(C_alpha.cols() - Cocc_alpha.cols());
		occ::Mat Cvir_beta = C_beta.rightCols(C_beta.cols() - Cocc_beta.cols());
		occ::Mat G_alpha = Cvir_alpha.transpose() * scf.ctx.F.topRows(sf->model_data.nmo) * Cocc_alpha;
		occ::Mat G_beta = Cvir_beta.transpose() * scf.ctx.F.bottomRows(sf->model_data.nmo) * Cocc_beta;
		return (std::hypot(G_alpha.norm(), G_beta.norm()));
	}
	err_not_impl_f("Orbital gradient for a general spinorbital kind", std::cout);
	return 0.0;
}

//F C = e S C through the orthogonaliser's X, one spin block at a time: X^T F X diagonalised
//by MKL, C = X C'. occ's MolecularOrbitals::update does the same with Eigen, which occ
//compiles serial; the occupation, smearing and density steps stay occ's.
void SCF_wrapper::solve_orbitals(occ::qm::SCF<occ::qm::HartreeFock>& scf, const occ::Mat& F) const {
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
occ::Mat SCF_wrapper::diis_update(occ::qm::SCF<occ::qm::HartreeFock>& scf) {
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
occ::Mat SCF_wrapper::cdiis_extrapolate(const occ::Mat& F, const occ::Mat& E) {
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
void SCF_wrapper::get_density_criteria(double& RMSP_diff, double& maxP_diff, const occ::Mat& dm, const occ::Mat& dm_last) {
	occ::Mat difference = dm - dm_last;
	RMSP_diff = std::sqrt(difference.squaredNorm() / difference.size());
	maxP_diff = difference.cwiseAbs().maxCoeff();
	// closing function
}

bool SCF_wrapper::SCF_iteration(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double& lambda, double& alpha, double& quant_diff_mem, double& quant, double& last_quant, occ::Mat& dm_last) {
	// Set up energy values & crystallographic information
	scf.iter++;
	const occ::Mat dm_old = scf.ctx.mo.D;
	dMatrix2 dm_eff(sf->model_data.nmo, sf->model_data.nmo);
	//This block is NoSpherA2 code, the Fock build below is OCC. Without the split the whole
	//remainder looks equally ours.
	const _time_point it_t0 = get_time();
	build_effective_dm(scf, dm_eff, dm_old);
	sf->calc_F_calc(*xcw->I_tens, dm_eff);
	sf->eval_scale();
	sf->calc_criteria();

	// Generates the perturbation matrix
	occ::Mat perturbation;
	xcw->calc_perturb(perturbation, scf);
	const _time_point it_t1 = get_time();
	throughput::record_time("XCW structure factors + perturbation", xcw->I_tens->i_on_device_, get_msec(it_t0, it_t1));

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

	const double current_criterion = sf->criterion(false);
	const double temp_penalty = current_criterion * lambda;
	quant = scf.ctx.energy["total"] + temp_penalty;

	scf.ctx.F += perturbation * lambda;

	// Prints output line for iteration
	SCF_log << "\t" << scf.iter << "\t\t" << std::fixed << std::setprecision(4) << current_criterion << "\t\t" << sf->quality_criteria.GooF2 << "\t\t" << std::setprecision(5) << sf->quality_criteria.R1 << "\t\t" << std::fixed << std::setprecision(9) << scf.ctx.energy["total"] << "\t\t" << std::fixed << std::setprecision(3) << temp_penalty << "\t\t" << std::fixed << std::setprecision(9) << quant << std::endl;

	//calc_perturb is the gradient of lambda * criterion^2 (the chi^2 of Jayatilaka's functional),
	//so that, not the printed lambda * criterion, is what the SCF descends and what the rescue ranks by
	const double phi = scf.ctx.energy["total"] + lambda * current_criterion * current_criterion;
	if (!soscf_ && rescue_scf(scf, phi, alpha)) return false;

	// DIIS extrapolation
	occ::Mat F_diis = diis_update(scf);
	opt->xcw_settings.current_max_diis_error = scf.diis_error;
	opt->xcw_settings.update(SCF_log, alpha);

	// Convergence check
	opt->xcw_settings.current_gradient = compute_orbital_gradient(scf);
	opt->xcw_settings.current_quant_diff = std::abs(quant - last_quant);
	if (SCF_convergence_check(scf, dm_last)) {
		return true;
	}
	last_quant = quant;

	if (!soscf_) {
		if (scf.iter == 1 || opt->xcw_settings.current_gradient < 0.5 * soscf_patience_grad_) {
			soscf_patience_grad_ = opt->xcw_settings.current_gradient;
			soscf_patience_iter_ = scf.iter;
		}
		//with soscf requested the DIIS stage only has to reach the quadratic region; Fe_phen HS lambda 0.08
		//oscillated for 30 iterations above the 1e-2 gate while TRAH converged in 6 from that very point
		const int patience = opt->xcw_settings.soscf ? trah_.patience_requested : trah_.patience;
		const bool stuck = scf.iter - soscf_patience_iter_ >= patience;
		if (stuck || (opt->xcw_settings.soscf && scf.diis_error < trah_.start_threshold)) {
			soscf_ = true;
			std::ostringstream what;
			what << "***" << (stuck ? "Orbital gradient not halved in " + std::to_string(patience) + " iterations" : "DIIS error below 1e-2")
				<< ": second-order steps on the orbital rotations from here***";
			print_centered_message(what.str(), 84, SCF_log);
			citations::cite(citations::Method::TRAH, SCF_log);
		}
	}
	if (soscf_) {
		soscf_step(scf, lambda, phi);
		dm_last = dm_old;
		return false;
	}

	// Apply level shift
	if (opt->xcw_settings.apply_shift) {
		const occ::Mat& C_old = scf.ctx.mo.C;
		apply_level_shift(C_old, scf, F_diis);
	}

	// Solves central eigenvalue problem
	solve_orbitals(scf, F_diis);

	// Apply damping
	if (opt->xcw_settings.apply_damping) {
		if (scf.iter == 2) {
			quant_diff_mem = opt->xcw_settings.current_quant_diff;
		}
		if (scf.iter > 2) {
			alpha = dynamic_damping(scf, alpha, opt->xcw_settings.current_quant_diff, quant_diff_mem);
		}
	}

	scf.ctx.mo.D *= (1 - alpha);
	scf.ctx.mo.D += alpha * dm_old;
	dm_last = dm_old;
	return false;

	//closing function
}

void SCF_wrapper::soscf_reset() {
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
occ::Vec SCF_wrapper::orbital_rotation_gradient(const occ::qm::SCF<occ::qm::HartreeFock>& scf, occ::Vec& diagonal_hessian) const {
	const occ::qm::MolecularOrbitals& mo = scf.ctx.mo;
	const int n = sf->model_data.nmo, nb = mo.kind == occ::qm::SpinorbitalKind::Unrestricted ? 2 : 1;
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
void SCF_wrapper::rotate_orbitals(occ::qm::SCF<occ::qm::HartreeFock>& scf, const occ::Mat& C_from, const occ::Vec& kappa) const {
	occ::qm::MolecularOrbitals& mo = scf.ctx.mo;
	const int n = sf->model_data.nmo, nb = mo.kind == occ::qm::SpinorbitalKind::Unrestricted ? 2 : 1;
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
void SCF_wrapper::soscf_step(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double lambda, const double phi) {
	const bool stepped = soscf_kappa_.size() > 0;
	if (stepped) {
		const double actual = phi - soscf_phi_;
		SCF_log << "\t\tTRAH: predicted " << std::scientific << std::setprecision(3) << soscf_pred_
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
			print_centered_message(give_up.str(), 84, SCF_log);
			//hand DIIS the orbitals the last accepted step left, not the rejected ones
			rotate_orbitals(scf, soscf_C_, occ::Vec::Zero(soscf_kappa_.size()));
			soscf_ = false;
			soscf_kappa_.resize(0);
			return;
		}
		std::ostringstream what;
		what << "***E + lambda chi^2 " << std::scientific << std::setprecision(1) << phi - soscf_phi_
			<< " Eh above the orbitals the step left: trust radius " << std::scientific << std::setprecision(2) << soscf_trust_ << ", step re-solved***";
		print_centered_message(what.str(), 84, SCF_log);
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
	if (opt->xcw_settings.check_hessian && !stepped) check_hessian(scf, lambda);
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
void SCF_wrapper::trah_solve(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double lambda, const bool extend) {
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
	SCF_log << "\t\tTRAH: " << micro << " micro-iterations (" << trah_micro_total_ << " in this lambda step), |g| " << std::scientific << std::setprecision(1) << gnorm
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
occ::Vec SCF_wrapper::hessian_vector(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double lambda, const occ::Vec& v) {
	occ::qm::MolecularOrbitals& mo = scf.ctx.mo;
	const int n = sf->model_data.nmo, nb = mo.kind == occ::qm::SpinorbitalKind::Unrestricted ? 2 : 1;
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
		const cvec F0 = sf->scatter_data.F_calc, fixed = sf->scatter_data.anom_correction;
		std::fill(sf->scatter_data.anom_correction.begin(), sf->scatter_data.anom_correction.end(), cdouble(0.0, 0.0));
		sf->calc_F_calc(*xcw->I_tens, dD_eff);
		const cvec dFr = sf->scatter_data.F_calc;
		sf->scatter_data.F_calc = F0;
		sf->scatter_data.anom_correction = fixed;
		//the scale's response from the weighted least-squares scale of eval_scale
		const bool against_F2 = opt->xcw_settings.refine_against == 2, weighted = opt->xcw_settings.XWR_type == 2;
		const double k = sf->scatter_data.scale, s = k * k;
		const int chunk = 128, nchunk = (sf->model_data.nr + chunk - 1) / chunk;
		vec dnum(nchunk), dden(nchunk), den(nchunk);
#pragma omp parallel for schedule(static)
		for (int ch = 0; ch < nchunk; ch++) {
			for (int r = ch * chunk; r < std::min((ch + 1) * chunk, sf->model_data.nr); r++) {
				if (!sf->scatter_data.hkl_mask[r]) continue;
				//ponytail: the extinction shape is frozen over the step (dy/d|Fc| dropped), so
				//sqrt(y)|Fc| and m d|Fc| reproduce I and dI exactly but their own curvature is
				//neglected. Only the Hessian's step proposal degrades - TRAH accepts on the
				//exact energy rebuild_at returns, and calc_perturb's gradient stays exact.
				const double w = weighted ? sf->inv_H2_[r] : 1.0, Fa = sf->ext_sqrt_y(r) * std::abs(F0[r]);
				const double dFa = std::abs(F0[r]) > 0 ? sf->ext_m(r) * std::real(std::conj(F0[r]) * dFr[r]) / std::abs(F0[r]) : 0.0;
				if (against_F2) {
					const double wi = w / (sf->scatter_data.sigma_obs2[r] * sf->scatter_data.sigma_obs2[r]);
					dnum[ch] += wi * 2.0 * Fa * dFa * sf->scatter_data.F_obs2[r];
					dden[ch] += wi * 4.0 * Fa * Fa * Fa * dFa;
					den[ch] += wi * Fa * Fa * Fa * Fa;
				}
				else {
					const double wi = w / (sf->scatter_data.sigma_obs[r] * sf->scatter_data.sigma_obs[r]);
					dnum[ch] += wi * dFa * sf->scatter_data.F_obs[r];
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
		cvec dq(sf->model_data.nr);
#pragma omp parallel for
		for (int r = 0; r < sf->model_data.nr; r++) {
			if (!sf->scatter_data.hkl_mask[r]) continue;
			const double w = weighted ? sf->inv_H2_[r] : 1.0, Fm = std::abs(F0[r]);
			if (Fm == 0) continue;
			//as above: the model amplitude is sqrt(y)|Fc| and its response m d|Fc|, with the
			//shape frozen. The carrier conj(F) keeps calc_perturb's chain factor.
			const double Fa = sf->ext_sqrt_y(r) * Fm, dFa = sf->ext_m(r) * std::real(std::conj(F0[r]) * dFr[r]) / Fm;
			const cdouble carrier = (against_F2 ? sf->ext_g(r) : sf->ext_m(r)) * std::conj(F0[r]);
			const cdouble dcarrier = (against_F2 ? sf->ext_g(r) : sf->ext_m(r)) * std::conj(dFr[r]);
			if (against_F2) {
				const double wi = w / (sf->scatter_data.sigma_obs2[r] * sf->scatter_data.sigma_obs2[r]), Fo2 = sf->scatter_data.F_obs2[r];
				dq[r] = wi * (dcarrier * (s * s * Fa * Fa - s * Fo2)
					+ carrier * (2.0 * s * dscale * Fa * Fa + 2.0 * s * s * Fa * dFa - dscale * Fo2));
			}
			else {
				const double wi = w / (sf->scatter_data.sigma_obs[r] * sf->scatter_data.sigma_obs[r]), Fo = sf->scatter_data.abs_F_obs[r];
				dq[r] = wi * (dcarrier * (k * k * sf->ext_sqrt_y(r) - k * Fo / Fm)
					+ carrier * (dscale * (2.0 * k * sf->ext_sqrt_y(r) - Fo / Fm) + k * Fo * std::real(std::conj(F0[r]) * dFr[r]) / (Fm * Fm * Fm)));
			}
		}
		occ::Mat dP;
		xcw->contract_I(dP, dq);
		dP *= (against_F2 ? 4.0 : 2.0) / (sf->model_data.nr_fit - sf->n_params()) * lambda;
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
double SCF_wrapper::rebuild_at(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double lambda, const occ::Mat& C_from, const occ::Vec& kappa) {
	if (kappa.size() > 0) rotate_orbitals(scf, C_from, kappa);
	else {
		scf.ctx.mo.update_occupied_orbitals();
		scf.ctx.mo.update_density_matrix();
	}
	dMatrix2 dm_eff(sf->model_data.nmo, sf->model_data.nmo);
	build_effective_dm(scf, dm_eff, scf.ctx.mo.D);
	sf->calc_F_calc(*xcw->I_tens, dm_eff);
	sf->eval_scale();
	sf->calc_criteria();
	occ::Mat perturbation;
	xcw->calc_perturb(perturbation, scf);
	scf.ctx.F = scf.ctx.H + (eri_ ? eri_fock(scf.ctx.mo, false) : scf.m_procedure.compute_fock(scf.ctx.mo, scf.ctx.K));
	scf.update_scf_energy(false);
	const double crit = sf->criterion(false);
	scf.ctx.F += lambda * perturbation;
	return scf.ctx.energy["total"] + lambda * crit * crit;
}

//orbital_rotation_gradient at C_from exp(kappa), rebuilt there and put back afterwards. Only
//the Hessian check pays for it.
occ::Vec SCF_wrapper::gradient_at(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double lambda, const occ::Mat& C_from, const occ::Vec& kappa) {
	const occ::qm::MolecularOrbitals mo_saved = scf.ctx.mo;
	const occ::Mat F_saved = scf.ctx.F;
	const auto energy_saved = scf.ctx.energy;
	const cvec F0_saved = sf->scatter_data.F_calc;
	const double scale_saved = sf->scatter_data.scale;
	rebuild_at(scf, lambda, C_from, kappa);
	occ::Vec h;
	const occ::Vec g = orbital_rotation_gradient(scf, h);
	scf.ctx.mo = mo_saved;
	scf.ctx.F = F_saved;
	scf.ctx.energy = energy_saved;
	sf->scatter_data.F_calc = F0_saved;
	sf->scatter_data.scale = scale_saved;
	return g;
}

//Central finite difference of the gradient along a fixed pseudo-random direction against
//hessian_vector, reported to XCW.log; `check_hessian` runs it on the first second-order step
void SCF_wrapper::check_hessian(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double lambda) {
	std::mt19937 rng(7);
	std::normal_distribution<double> gauss;
	occ::Vec v(soscf_grad_.size());
	for (Eigen::Index i = 0; i < v.size(); i++) v(i) = gauss(rng);
	v /= v.norm();
	const occ::Vec Hv = hessian_vector(scf, lambda, v);
	const double eps = 1e-4;
	const occ::Mat C0 = scf.ctx.mo.C;
	const occ::Vec fd = (gradient_at(scf, lambda, C0, eps * v) - gradient_at(scf, lambda, C0, -eps * v)) / (2.0 * eps);
	SCF_log << "\t\tHessian check: |Hv - FD| / |FD| = " << std::scientific << std::setprecision(2) << (Hv - fd).norm() / fd.norm()
		<< " (|Hv| " << Hv.norm() << ", |FD| " << fd.norm() << ", v.Hv " << v.dot(Hv) << ", v.FD " << v.dot(fd) << ")" << std::endl;
}

//The rescue of a lost SCF step, see best_mo_ in the header. True when it restarted the step:
//the orbitals are the best ones again, the iteration ends here and the next one rebuilds
//their Fock matrix in full. The stabiliser thresholds stay tightened for the later steps,
//which are at least as hard.
bool SCF_wrapper::rescue_scf(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double quant, double& alpha) {
	if (quant < best_quant_) {
		best_quant_ = quant;
		best_mo_ = scf.ctx.mo;
		return false;
	}
	if (quant - best_quant_ < rescue_rise_ || rescues_ >= 3) return false;
	rescues_++;
	opt->xcw_settings.diis_stop_shift /= 10;
	opt->xcw_settings.diis_stop_damping /= 10;
	std::ostringstream what;
	what << "***E + lambda chi^2 " << std::fixed << std::setprecision(3) << quant - best_quant_
		<< " Eh above its best: back to the best orbitals, level shift and damping on until DIIS error "
		<< std::scientific << std::setprecision(0) << opt->xcw_settings.diis_stop_shift << " / " << opt->xcw_settings.diis_stop_damping
		<< " (rescue " << rescues_ << "/3)***";
	print_centered_message(what.str(), 84, SCF_log);
	scf.ctx.mo = best_mo_;
	G_last_.resize(0, 0);
	diis_F_.clear();
	diis_E_.clear();
	adiis_.reset();
	ediis_.reset();
	opt->xcw_settings.apply_shift = opt->xcw_settings.method_apply_shift;
	opt->xcw_settings.apply_damping = opt->xcw_settings.method_apply_damping;
	alpha = opt->xcw_settings.alpha;
	return true;
}

bool SCF_wrapper::SCF_convergence_check(occ::qm::SCF<occ::qm::HartreeFock>& scf, occ::Mat& dm_last) {
	get_density_criteria(opt->xcw_settings.current_RMSP_diff, opt->xcw_settings.current_MaxP_diff, scf.ctx.mo.D, dm_last);
	opt->xcw_settings.conv_quant_diff = opt->xcw_settings.current_quant_diff < opt->xcw_settings.quant_diff;
	opt->xcw_settings.conv_max_diis_error = opt->xcw_settings.current_max_diis_error < opt->xcw_settings.max_diis_error;
	opt->xcw_settings.conv_gradient = opt->xcw_settings.current_gradient < opt->xcw_settings.gradient;
	opt->xcw_settings.conv_RMSP_diff = opt->xcw_settings.current_RMSP_diff < opt->xcw_settings.RMSP_diff;
	opt->xcw_settings.conv_MaxP_diff = opt->xcw_settings.current_MaxP_diff < opt->xcw_settings.MaxP_diff;
	return opt->xcw_settings.convergence_check();
	// closing function
}

//The stored integrals' two-electron Fock part, contracted on the device when it holds them
occ::Mat SCF_wrapper::eri_fock(const occ::qm::MolecularOrbitals& mo, bool screen) const {
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
