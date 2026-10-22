#include "pch.h"
#include "SCF_wrapper.h"
#if defined(NOSPHERA2_USE_GPU) || defined(NOSPHERA2_USE_METAL)
#include "itensor_gpu.h"
#endif
#include "basis_set.h"
#include "citations.h"
#include <limits>
#include <random>

SCF_wrapper::SCF_wrapper(const occ::core::Molecule& mol, const occ::qm::AOBasis& basis, const options::SCF_settings& settings_in, SCF_log_writer& writer) : settings(settings_in), mol_(mol), basis_(basis), nmo_(static_cast<int>(basis.nbf())), writer_(&writer) {}

//C = op(A) op(B) for column-major occ matrices through MKL, Src/core's Eigen being serial
static occ::Mat gemm(const occ::Mat& A, const occ::Mat& B, const bool ta = false, const bool tb = false)
{
	const int m = static_cast<int>(ta ? A.cols() : A.rows()), k = static_cast<int>(ta ? A.rows() : A.cols()), n = static_cast<int>(tb ? B.rows() : B.cols());
	occ::Mat C(m, n);
	cblas_dgemm(CblasColMajor, ta ? CblasTrans : CblasNoTrans, tb ? CblasTrans : CblasNoTrans, m, n, k, 1.0, A.data(), static_cast<int>(A.rows()), B.data(), static_cast<int>(B.rows()), 0.0, C.data(), m);
	return C;
}

occ::core::Molecule SCF_wrapper::setup_SCF_mol(const std::vector<asym_atom>& atoms, const int charge, const int multiplicity) {
	double bohr2angstrom = constants::bohr2ang(1);
	std::ostringstream init_stream;

	const int ncen = atoms.size();

	init_stream << ncen << "\n\n";

	for (int i = 0; i < ncen; ++i) {
		init_stream
			<< constants::atnr2letter(atoms[i].type) << " "
			<< atoms[i].pos[0] * bohr2angstrom << " "
			<< atoms[i].pos[1] * bohr2angstrom << " "
			<< atoms[i].pos[2] * bohr2angstrom;

		if (i != ncen - 1)
			init_stream << "\n";
	}

	occ::core::Molecule mol = occ::io::molecule_from_xyz_string(init_stream.str());
	mol.set_charge(charge);
	mol.set_multiplicity(multiplicity);
	return mol;
}

//Four fifths of what the I tensor left, under -mem when given, for the stored integrals
void SCF_wrapper::setup_eri(const occ::qm::HartreeFock& hf, const size_t avail) {
	eri_.clear();
	const _time_point eri_t0 = get_time();
	//At least 1: a budget of 0 means "no budget" to stored_eri
	const size_t budget = avail ? std::max<size_t>(avail / 5 * 4, 1) : 0;
	const bool read = !settings.eri_read_path.empty() && eri_.load(settings.eri_read_path, hf, budget, writer_->log);
	if (!read && eri_.build(hf, budget, writer_->log))
		throughput::record_time("XCW two-electron integrals", false, get_msec(eri_t0, get_time()));
	if (eri_ && !settings.eri_save_path.empty() && !eri_.save(settings.eri_save_path, hf))
		writer_->log << "XCW: could not write the two-electron integrals to " << settings.eri_save_path.string() << std::endl;
#if defined(NOSPHERA2_USE_GPU) || defined(NOSPHERA2_USE_METAL)
	if (eri_ && settings.gpu_eri) {
		const auto up_t0 = get_time();
		eri_on_device_ = eri_gpu_hold(eri_.data(), eri_.nbf(), eri_.npairs(), eri_.pair_a().data(), eri_.pair_b().data(),
			eri_.first_pair().data(), eri_.pair_index().data());
		if (eri_on_device_) throughput::record_time("XCW two-electron integrals upload", true, get_msec(up_t0, get_time()));
		if (!settings.no_date)
			std::cerr << "GPU in use: XCW Fock build from the stored integrals on "
			<< (eri_on_device_ ? "the device" : "the CPU - device unavailable, an APU without fast fp64, or the integrals too large") << std::endl;
	}
#endif
}

void SCF_wrapper::release_device() {
#if defined(NOSPHERA2_USE_GPU) || defined(NOSPHERA2_USE_METAL)
	eri_gpu_release();
#endif
	eri_on_device_ = false;
}

bool SCF_wrapper::solve(occ::qm::HartreeFock& hf, const double lambda, occ::qm::Wavefunction& wfn, const bool use_guess, SCF_coupling& coupling) {
	coupling_ = &coupling;
	//the stabilisers as configured, whatever the last lambda step's rescues made of them: divided
	//per rescue and never restored they reached 1e-33 by TEMPO's lambda 0.11
	if (base_.level_shift < 0)
		base_ = { settings.diis_stop_shift, settings.diis_stop_damping, settings.level_shift };
	settings.diis_stop_shift = base_.stop_shift;
	settings.diis_stop_damping = base_.stop_damping;
	settings.level_shift = base_.level_shift;
	occ::qm::SCF scf(hf, settings.hf_type);
	double alpha = settings.alpha;
	bool has_local_guess = use_guess;
	scf.set_charge_multiplicity(settings.charge, settings.multiplicity);
	scf.maxiter = settings.max_scf_iterations;
	scf.convergence_settings.level_shift = settings.level_shift;
	scf.convergence_settings.level_shift_threshold = 0;
	scf.update_occupied_orbital_count();
	const bool converged = do_SCF(lambda, alpha, scf, wfn, has_local_guess);
	wfn = scf.wavefunction();
	coupling_ = nullptr;
	return converged;
}

//The same walk as OCC's SOAD guess with a converged Hartree-Fock density in place of the
//tabulated atoms: F = H + G(D_small) through the mixed-basis build, the orbitals from its
//diagonalisation. The spin blocks are summed, the first iteration polarises them again.
void SCF_wrapper::small_basis_guess(occ::qm::SCF<occ::qm::HartreeFock>& scf) {
	const occ::qm::AOBasis small_bs = setup_basis(mol_, settings.guess_basis_name);
	occ::qm::HartreeFock hf_small(small_bs);
	occ::qm::SCF scf_small(hf_small, settings.hf_type);
	scf_small.set_charge_multiplicity(settings.charge, settings.multiplicity);
	scf_small.maxiter = settings.max_scf_iterations;
	const double e_small = scf_small.compute_scf_energy();
	writer_->log << "XCW: initial guess from a " << settings.guess_basis_name << " Hartree-Fock (" << small_bs.nbf()
		<< " functions), E = " << std::fixed << std::setprecision(8) << e_small << " Eh" << std::endl;

	occ::qm::MolecularOrbitals guess;
	guess.kind = scf.ctx.mo.kind;
	guess.n_ao = static_cast<int>(small_bs.nbf());
	guess.n_alpha = scf.n_alpha();
	guess.n_beta = scf.n_beta();
	guess.D = settings.hf_type == occ::qm::SpinorbitalKind::Unrestricted
		? occ::Mat(occ::qm::block::a(scf_small.ctx.mo.D) + occ::qm::block::b(scf_small.ctx.mo.D))
		: scf_small.ctx.mo.D;

	scf.update_occupied_orbital_count();
	scf.set_core_matrices();
	scf.ctx.F = scf.ctx.H;
	scf.set_conditioning_orthogonalizer();
	const bool same_basis = guess.kind == occ::qm::SpinorbitalKind::Restricted && small_bs.nbf() == scf.m_procedure.aobasis().nbf() && small_bs == scf.m_procedure.aobasis();
	scf.ctx.F += same_basis ? scf.m_procedure.compute_fock_from_density(guess) : scf.m_procedure.compute_fock_mixed_basis(guess, small_bs, false);
	scf.ctx.orthogonalizer.orthogonalize_molecular_orbitals(scf.ctx.mo, scf.ctx.F);
	scf.m_have_initial_guess = true;
}

occ::qm::AOBasis SCF_wrapper::setup_basis(const occ::core::Molecule& mol, const std::string& basis_set_name) {
	std::shared_ptr<BasisSet> basis_set = BasisSetLibrary::get_basis_set(basis_set_name);
	return basis_set->to_AOBasis(mol.atoms());
}

double SCF_wrapper::dynamic_damping(const occ::qm::SCF<occ::qm::HartreeFock>& scf, const double& current_alpha, const double& quant_diff, double& quant_diff_mem) {
	double new_alpha = current_alpha;
	if (quant_diff < quant_diff_mem / 10) {
		new_alpha *= 0.75;
		quant_diff_mem = quant_diff;
		if (quant_diff < 10 * scf.convergence_settings.energy_threshold) {
			print_centered_message("***Turned off damping***", 100, writer_->log);
			new_alpha = 0;
			settings.apply_damping = false;
		}
		else {
			std::stringstream print_;
			print_ << "***Decreased damping to " << std::fixed << std::setprecision(3) << new_alpha << "***";
			print_centered_message(print_.str(), 100, writer_->log);
		}
	}
	return new_alpha;
	// closing function
}

void SCF_wrapper::apply_level_shift(const occ::Mat& C_old, const occ::qm::SCF<occ::qm::HartreeFock>& scf, occ::Mat& F_diis) {
	const int nocc = scf.ctx.mo.Cocc.cols();
	if (scf.ctx.mo.kind == occ::qm::SpinorbitalKind::Restricted) {
		const occ::Mat SC_virt = scf.ctx.S * C_old.rightCols(nmo_ - nocc);
		F_diis.noalias() += settings.level_shift * SC_virt * SC_virt.transpose();
	}
	else {
		const int nao = C_old.rows() / 2;
		const auto S_ao = scf.ctx.S.topRows(nao);
		//Cocc has max(n_alpha, n_beta) columns, the spin counts differ for open shells
		const occ::Mat SC_virt_a = S_ao * C_old.topRows(nao).rightCols(nmo_ - scf.ctx.mo.n_alpha);
		const occ::Mat SC_virt_b = S_ao * C_old.bottomRows(nao).rightCols(nmo_ - scf.ctx.mo.n_beta);
		F_diis.topRows(nao).noalias() += settings.level_shift * SC_virt_a * SC_virt_a.transpose();
		F_diis.bottomRows(nao).noalias() += settings.level_shift * SC_virt_b * SC_virt_b.transpose();
	}
}

void SCF_wrapper::build_effective_dm(const occ::qm::SCF<occ::qm::HartreeFock>& scf, dMatrix2& dm_ref, const occ::Mat& dm_old) {
	if (scf.ctx.mo.kind == occ::qm::SpinorbitalKind::Unrestricted) {
		for (int i = 0; i < nmo_; i++) {
			dm_ref(i, i) = dm_old(i, i);
			dm_ref(i, i) += dm_old(i + nmo_, i);
			for (int j = i + 1; j < nmo_; j++) {
				dm_ref(i, j) = dm_old(i, j);
				dm_ref(i, j) += dm_old(i + nmo_, j);
				dm_ref(j, i) = dm_ref(i, j);
			}
		}
	}
	else {
		for (int i = 0; i < nmo_; i++) {
			dm_ref(i, i) = dm_old(i, i);
			for (int j = i + 1; j < nmo_; j++) {
				dm_ref(i, j) = dm_old(i, j);
				dm_ref(j, i) = dm_ref(i, j);
			}
		}
	}
}

bool SCF_wrapper::do_SCF(const double& lambda, double& alpha, occ::qm::SCF<occ::qm::HartreeFock>& scf, occ::qm::Wavefunction& last_wfn, bool& has_guess) {

#if !defined(__APPLE__) && !defined(NSA2_OPENBLAS)
	//The SCF's own MKL calls (DIIS, the diagonalisation) run serial on every path: they are small against
	//the Fock build, and threaded they wait for a time slice at every barrier once another job loads the
	//host (Radeon 780M box beside a 16-thread job: 170 s against 7 s for 14 iterations). A CPU or hybrid I tensor
	//has left this thread serial anyway, through mkl_set_num_threads_local in build_I.
	mkl_set_num_threads_local(1);
#endif
	settings.clear();
	diis_F_.clear();
	diis_E_.clear();
	adiis_.reset();
	ediis_.reset();
	best_quant_ = std::numeric_limits<double>::infinity();
	rescues_ = 0;
	builds_ = {};
	soscf_reset();
	warm_start_ = has_guess;

	writer_->start_lambda(lambda);

	// Compute first guess and update the energy according to this guess
	const _time_point guess_t0 = get_time();
	if (has_guess) {
		scf.set_initial_guess_from_wfn(last_wfn);
	}
	else {
		if (settings.guess_basis_name.empty())
			scf.compute_initial_guess();
		else
			small_basis_guess(scf);
		has_guess = true;
	}
	scf.ctx.K = scf.m_procedure.compute_schwarz_ints();
	scf.update_scf_energy(false);
	G_last_.resize(0, 0);
	reuse_best_G_ = false;
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
		writer_->log << "____________________________________________________________________________________________________\n";
		std::stringstream print_;
		print_ << "***SCF converged in " << scf.iter << " iterations***";
		print_centered_message(print_.str(), 100, writer_->log);
		writer_->log << "Fock builds: " << builds_.iter + builds_.hessian + builds_.rebuild << " (" << builds_.iter << " iterations, "
			<< builds_.hessian << " Hessian products, " << builds_.rebuild << " rebuilds)" << std::endl;
		writer_->v.lambda = lambda;
		writer_->v.E_total = scf.ctx.energy["total"];
		writer_->v.quant = quant;
	}
	else {
		writer_->log << "____________________________________________________________________________________________________\n";
		print_centered_message("***SCF did not converge***", 100, writer_->log);
		writer_->log << "Fock builds: " << builds_.iter + builds_.hessian + builds_.rebuild << " (" << builds_.iter << " iterations, "
			<< builds_.hessian << " Hessian products, " << builds_.rebuild << " rebuilds)" << std::endl;
		std::ostringstream perturbed_energy;
		perturbed_energy << " for perturbed energy: " << std::scientific << settings.quant_diff << " (current: " << settings.current_quant_diff << ") \n";
		std::ostringstream diis_error;
		diis_error << " for DIIS error: " << std::scientific << settings.max_diis_error << " (current: " << settings.current_max_diis_error << ") \n";
		std::ostringstream orbital_gradient;
		orbital_gradient << " for orbital gradient: " << std::scientific << settings.gradient << " (current: " << settings.current_gradient << ") \n";
		std::ostringstream max_density_diff;
		max_density_diff << " for maximum difference in density matrix: " << std::scientific << settings.MaxP_diff << " (current: " << settings.current_MaxP_diff << ") \n";
		std::ostringstream rmsd_density;
		rmsd_density << " for RMSD of density matrix: " << std::scientific << settings.RMSP_diff << " (current: " << settings.current_RMSP_diff << ") \n";
		if (settings.conv_quant_diff) {
			writer_->log << "CONVERGED";
		}
		else {
			writer_->log << "NOT CONVERGED";
		}
		writer_->log << perturbed_energy.str();
		if (settings.conv_max_diis_error) {
			writer_->log << "CONVERGED";
		}
		else {
			writer_->log << "NOT CONVERGED";
		}
		writer_->log << diis_error.str();
		if (settings.conv_gradient) {
			writer_->log << "CONVERGED";
		}
		else {
			writer_->log << "NOT CONVERGED";
		}
		writer_->log << orbital_gradient.str();
		if (settings.conv_MaxP_diff) {
			writer_->log << "CONVERGED";
		}
		else {
			writer_->log << "NOT CONVERGED";
		}
		writer_->log << max_density_diff.str();
		if (settings.conv_RMSP_diff) {
			writer_->log << "CONVERGED";
		}
		else {
			writer_->log << "NOT CONVERGED";
		}
		writer_->log << rmsd_density.str();
	}
	return converged;
}

double SCF_wrapper::compute_orbital_gradient(const occ::qm::SCF<occ::qm::HartreeFock>& scf) {
	if (settings.hf_type == occ::qm::SpinorbitalKind::Restricted) {
		const occ::Mat& C = scf.molecular_orbitals().C;
		const occ::Mat& Cocc = scf.molecular_orbitals().Cocc;
		const occ::Mat Cvir = C.rightCols(C.cols() - Cocc.cols());
		return 2.0 * gemm(Cvir, gemm(scf.ctx.F, Cocc), true).norm();
	}
	else if (settings.hf_type == occ::qm::SpinorbitalKind::Unrestricted) {
		occ::Mat C_alpha = scf.molecular_orbitals().C.topRows(nmo_);
		occ::Mat C_beta = scf.molecular_orbitals().C.bottomRows(nmo_);
		occ::Mat Cocc_alpha = C_alpha.leftCols(scf.ctx.mo.n_alpha);
		occ::Mat Cocc_beta = C_beta.leftCols(scf.ctx.mo.n_beta);
		occ::Mat Cvir_alpha = C_alpha.rightCols(C_alpha.cols() - Cocc_alpha.cols());
		occ::Mat Cvir_beta = C_beta.rightCols(C_beta.cols() - Cocc_beta.cols());
		occ::Mat G_alpha = Cvir_alpha.transpose() * scf.ctx.F.topRows(nmo_) * Cocc_alpha;
		occ::Mat G_beta = Cvir_beta.transpose() * scf.ctx.F.bottomRows(nmo_) * Cocc_beta;
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
	dMatrix2 dm_eff(nmo_, nmo_);
	//This block is NoSpherA2 code, the Fock build below is OCC. Without the split the whole
	//remainder looks equally ours.
	const _time_point it_t0 = get_time();
	build_effective_dm(scf, dm_eff, dm_old);
	coupling_->evaluate(dm_eff);

	// Generates the perturbation matrix
	occ::Mat perturbation;
	coupling_->perturbation(perturbation, scf.ctx.mo.kind == occ::qm::SpinorbitalKind::Unrestricted);
	const _time_point it_t1 = get_time();
	throughput::record_time("XCW structure factors + perturbation", coupling_->on_device(), get_msec(it_t0, it_t1));

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
	//against the value at the orbitals it left, and the noise puts that threshold at 1e-8 Eh
	//(1e-10 |phi| in heavy cells), far below the screening error of a difference build. Comparing one against the other made
	//good steps read as rises, and a rejection can only shrink the trust radius, so a single
	//artefact capped it for the rest of the lambda step. A rotation moves the density by a whole
	//step anyway, so there was little left to skip. OCC's own loop does the same (it forces a
	//full rebuild while the second-order step is active).
	const bool incremental = settings.incremental && !soscf_ && !eri_on_device_ && scf.m_procedure.fock_build_properties().density_screened && G_last_.size() > 0
		&& scf.iter - last_full_build_ < 8 && scf.diis_error > next_full_build_error_;
	//a rescue restored best_mo_ and its G: the density is the one G was built from
	const bool reused = std::exchange(reuse_best_G_, false);
	if (reused) {}
	else if (incremental) {
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
	if (!reused) builds_.iter++;
	if (!reused) throughput::record_time("XCW Fock build (OCC)", eri_on_device_, get_msec(fock_t0, get_time()));
	scf.update_scf_energy(false);

	const fit_quality fit = coupling_->quality();
	const double current_criterion = fit.criterion;
	const double temp_penalty = current_criterion * lambda;
	quant = scf.ctx.energy["total"] + temp_penalty;

	scf.ctx.F += perturbation * lambda;

	// Prints output line for iteration
	writer_->v = { scf.iter, lambda, current_criterion, fit.GooF2, fit.R1, scf.ctx.energy["total"], temp_penalty, quant, fit.criterion_all, fit.R1_all, 0.0, false };
	writer_->iteration_line();

	//perturbation is the gradient of lambda * criterion^2 (the chi^2 of Jayatilaka's functional),
	//so that, not the printed lambda * criterion, is what the SCF descends and what the rescue ranks by
	const double phi = scf.ctx.energy["total"] + lambda * current_criterion * current_criterion;
	//a lost step goes to TRAH from the best orbitals where a rescue would not help: a warm start
	//sits next to its solution, so DIIS overshot it, and a cold start that used its three rescues
	//has nothing left (TEMPO lambda 0.11: all three gone by iteration 12, TRAH after 30 more)
	const bool lost = !soscf_ && phi - best_quant_ >= rescue_rise_ && (warm_start_ || rescues_ >= 3);
	if (!soscf_ && !lost && rescue_scf(scf, phi, alpha)) return false;

	// DIIS extrapolation
	occ::Mat F_diis = diis_update(scf);
	settings.current_max_diis_error = scf.diis_error;
	settings.update(writer_->log, alpha);

	// Convergence check, the change on phi = E + lambda chi^2, the functional that is minimised and
	// stationary: quant (E + lambda chi) is not, so its change stalls above the threshold while
	// every other criterion is met (TEMPO lambda 0.01: 20 micro-iterations at |g| 1.6e-4)
	settings.current_gradient = compute_orbital_gradient(scf);
	settings.current_quant_diff = std::abs(phi - last_quant);
	//Under TRAH the density criteria judge the step still to come, the one trah_solve takes without
	//a Fock build once |g| is below the threshold, not the step just taken: judged on that, every
	//warm lambda step spent an iteration only to see its last step was small (TEMPO, P1). The
	//orbitals returned are the current ones, consistent with their Fock matrix and Fcalc.
	if (soscf_ && settings.current_gradient < settings.gradient) {
		occ::Vec h;
		const occ::Vec g = orbital_rotation_gradient(scf, h);
		std::swap(h, soscf_hdiag_);
		occ::Vec kappa = trah_precondition(scf, lambda, -g, 0);
		std::swap(h, soscf_hdiag_);
		if (kappa.norm() > soscf_trust_) kappa *= soscf_trust_ / kappa.norm();
		const occ::qm::MolecularOrbitals here = scf.ctx.mo;
		rotate_orbitals(scf, here.C, kappa);
		dm_last = scf.ctx.mo.D;
		scf.ctx.mo = here;
	}
	if (SCF_convergence_check(scf, dm_last)) {
		return true;
	}
	last_quant = phi;

	if (!soscf_) {
		if (scf.iter == 1 || settings.current_gradient < 0.5 * soscf_patience_grad_) {
			soscf_patience_grad_ = settings.current_gradient;
			soscf_patience_iter_ = scf.iter;
		}
		//A cold start runs DIIS: one Fock build per iteration beats a TRAH macro step far from the
		//solution (lambda 0: P1 14 builds against 26, TEMPO TZVPP 17 against 31). A warm start goes
		//straight to the second-order steps: it sits in the quadratic region, where DIIS lost its first
		//step on TEMPO (the gradient x4 to x6 at iteration 2) and crawled on P1 (x0.6 per iteration)
		//check_hessian needs a second-order step to check, so it enters TRAH at the DIIS error gate
		const bool second_order = settings.check_hessian;
		const int patience = second_order ? trah_.patience_requested : trah_.patience;
		const bool stuck = scf.iter - soscf_patience_iter_ >= patience;
		const bool warm = warm_start_;
		if (warm || lost || stuck || (second_order && scf.diis_error < trah_.start_threshold)) {
			soscf_ = true;
			occ::qm::MolecularOrbitals probe = scf.ctx.mo;
			probe.update_occupied_orbitals();
			probe.update_density_matrix();
			soscf_fresh_ = last_full_build_ == scf.iter && (probe.D - scf.ctx.mo.D).cwiseAbs().maxCoeff() < 1e-10;
			std::ostringstream what;
			what << "***" << (warm ? std::string("Warm start") : lost ? "E + lambda chi^2 " + std::to_string(phi - best_quant_) + " Eh above its best"
				: stuck ? "Orbital gradient not halved in " + std::to_string(patience) + " iterations" : "DIIS error below 1e-2")
				<< ": second-order steps on the orbital rotations from here***";
			print_centered_message(what.str(), 100, writer_->log);
			citations::cite(citations::Method::TRAH, writer_->log, 100);
			//TRAH from the best orbitals of this lambda step, not from wherever DIIS left them: TEMPO
			//lambda 0.11 (Cordula, AKL) switched at E -148 Eh, 330 Eh above the warm start, and spent
			//900 micro-iterations climbing back. The first step rebuilds the Fock matrix at them.
			if (phi > best_quant_ + trah_.noise) {
				writer_->log << "\t\tTRAH starts from the best orbitals of this step, " << std::scientific << std::setprecision(2)
					<< phi - best_quant_ << " Eh below the current ones" << std::endl;
				scf.ctx.mo = best_mo_;
				soscf_fresh_ = false;
			}
		}
	}
	if (soscf_) {
		soscf_step(scf, lambda, phi);
		dm_last = dm_old;
		return false;
	}

	// Apply level shift
	if (settings.apply_shift) {
		const occ::Mat& C_old = scf.ctx.mo.C;
		apply_level_shift(C_old, scf, F_diis);
	}

	// Solves central eigenvalue problem
	solve_orbitals(scf, F_diis);

	// Apply damping
	if (settings.apply_damping) {
		if (scf.iter == 2) {
			quant_diff_mem = settings.current_quant_diff;
		}
		if (scf.iter > 2) {
			alpha = dynamic_damping(scf, alpha, settings.current_quant_diff, quant_diff_mem);
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
	trah_pcg_total_ = 0;
}

//The gradient of E + lambda chi^2 with respect to the rotation kappa_ai of occupied orbital i
//into virtual a, C -> C exp(kappa), one spin block after the other: 4 F_ai for the restricted
//density 2 C_occ C_occ^T, 2 F_ai per spin block otherwise, F the perturbed Fock matrix in the
//MO basis. diagonal_hessian gets the usual approximation of its diagonal, 4 (e_a - e_i) or
//2 (e_a - e_i) over the diagonal of that F, floored so a near-degenerate pair cannot blow up
//the step.
occ::Vec SCF_wrapper::orbital_rotation_gradient(const occ::qm::SCF<occ::qm::HartreeFock>& scf, occ::Vec& diagonal_hessian) const {
	const occ::qm::MolecularOrbitals& mo = scf.ctx.mo;
	const int n = nmo_, nb = mo.kind == occ::qm::SpinorbitalKind::Unrestricted ? 2 : 1;
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
	const int n = nmo_, nb = mo.kind == occ::qm::SpinorbitalKind::Unrestricted ? 2 : 1;
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
	//A step the model expects to gain less than the noise of a Fock build cannot be judged by what
	//the functional does: P1 lambda 0.05, predicted -1e-8 Eh, actual +3e-8 at -1962 Eh, rejected,
	//and the radius it fell to (1e-4) cost the next lambda step 8 macro steps doubling it back.
	//It keeps the radius and is taken unless the rise is beyond any noise.
	//The noise grows with the functional: P1 lambda 0.04 at -1962 Eh, predicted -1.8e-8, actual +3.4e-8,
	//judged at the fixed 1e-8, and the radius it fell to cost lambda 0.05 seven extra macro steps.
	const double noise = std::max(trah_.noise, 1e-10 * std::abs(phi));
	const bool judged = -soscf_pred_ > noise;
	if (stepped) {
		const double actual = phi - soscf_phi_;
		writer_->log << "\t\tTRAH: predicted " << std::scientific << std::setprecision(3) << soscf_pred_
			<< ", actual " << actual << " Eh, rho " << std::fixed << std::setprecision(2)
			<< (soscf_pred_ < 0 ? actual / soscf_pred_ : 0.0) << ", model error "
			<< std::scientific << std::setprecision(1) << std::abs(actual - soscf_pred_) << std::endl;
		if (judged) soscf_trust_ = occ::qm::trust_radius_update(soscf_trust_, soscf_kappa_.norm(), soscf_pred_, actual, soscf_boundary_, trah_);
	}
	if (stepped && phi > soscf_phi_ + (judged ? noise : unjudged_rise_) && soscf_kappa_.cwiseAbs().maxCoeff() > 1e-6) {
		//A radius this small is the model saying it is worthless at these orbitals, not that the
		//step should be shorter again. Two in a row and DIIS gets the orbitals back, with the
		//rescue path available again, rather than the lambda step grinding out max_iter on steps
		//too short to move anything.
		soscf_floored_ = soscf_trust_ <= trah_.trust_min * 1.000001 ? soscf_floored_ + 1 : 0;
		if (soscf_floored_ >= 2) {
			std::ostringstream give_up;
			give_up << "***E + lambda chi^2 still rising at the smallest trust radius: back to DIIS***";
			print_centered_message(give_up.str(), 100, writer_->log);
			//hand DIIS the orbitals the last accepted step left, not the rejected ones
			rotate_orbitals(scf, soscf_C_, occ::Vec::Zero(soscf_kappa_.size()));
			soscf_ = false;
			soscf_kappa_.resize(0);
			return;
		}
		std::ostringstream what;
		what << "***E + lambda chi^2 " << std::scientific << std::setprecision(1) << phi - soscf_phi_
			<< " Eh above the orbitals the step left: trust radius " << std::scientific << std::setprecision(2) << soscf_trust_ << ", step re-solved***";
		print_centered_message(what.str(), 100, writer_->log);
		trah_solve(scf, lambda, false);
		rotate_orbitals(scf, soscf_C_, soscf_kappa_);
		return;
	}
	soscf_floored_ = 0;
	//The Roothaan iterations damp the density, so the Fock matrix of the iteration that hands
	//over belongs to a mix of densities, not to the orbitals: rebuilt from them once, so the
	//first gradient, Hessian and reference functional are consistent
	soscf_phi_ = stepped || std::exchange(soscf_fresh_, false) ? phi : rebuild_at(scf, lambda, scf.ctx.mo.C, occ::Vec());
	soscf_grad_ = orbital_rotation_gradient(scf, soscf_hdiag_);
	trah_B_.clear();
	trah_HB_.clear();
	soscf_C_ = scf.ctx.mo.C;
	if (settings.check_hessian && !stepped) check_hessian(scf, lambda);
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
	//never below the convergence threshold on the orbital gradient: a residual under it buys nothing,
	//and the last macro step spent 20 products on TEMPO lambda 0.01 getting there
	const double tol = std::max(gnorm * std::min(0.1, std::max(1e-3, std::sqrt(gnorm))), settings.gradient);
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
	//a gradient already below the convergence threshold only needs the step that settles the density:
	//the one from the Hessian without its two-electron part, no Fock build (P1: one per lambda step)
	if (extend && trah_B_.empty() && gnorm < tol) {
		soscf_kappa_ = trah_precondition(scf, lambda, -g, 0);
		if (soscf_kappa_.norm() > soscf_trust_) soscf_kappa_ *= soscf_trust_ / soscf_kappa_.norm();
		soscf_pred_ = g.dot(soscf_kappa_) + 0.5 * soscf_kappa_.dot(hessian_vector(scf, lambda, soscf_kappa_, false));
		soscf_boundary_ = false;
		writer_->log << "\t\tTRAH: |g| " << std::scientific << std::setprecision(1) << gnorm << " below the threshold, step without a Fock build, predicted " << soscf_pred_ << " Eh" << std::endl;
		return;
	}
	if (extend && trah_B_.empty()) add(trah_precondition(scf, lambda, -g, 0));
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
		if (!add(trah_precondition(scf, lambda, -r, theta))) break;
		micro++;
	}
	trah_micro_total_ += micro;
	soscf_kappa_ = kappa;
	soscf_pred_ = g.dot(kappa) + 0.5 * kappa.dot(Hkappa);
	writer_->log << "\t\tTRAH: " << micro << " micro-iterations (" << trah_micro_total_ << " in this lambda step, " << trah_pcg_total_ << " PCG), |g| " << std::scientific << std::setprecision(1) << gnorm
		<< ", residual " << rnorm << ", |kappa| " << knorm << (soscf_boundary_ ? " on" : " within") << " the trust radius " << std::fixed << std::setprecision(3) << soscf_trust_
		<< ", shift " << std::scientific << std::setprecision(1) << theta << ", predicted " << soscf_pred_ << " Eh" << std::endl;
}

//H v for the rotation direction v, the derivative of orbital_rotation_gradient along it:
//dC_occ = C_vir V, dC_vir = -C_occ V^T from C exp(kappa); dD from the density's own
//convention (C_occ C_occ^T, halved per spin block for UHF); dF = G(dD) by the same Fock
//build the SCF uses, plus lambda times the response of the perturbation, which is the
//derivative of perturbation's per-reflection scalar with the scale moving along (the scale
//minimises chi^2, so the gradient does not see it but the Hessian does); and
//H v = fac (dC_vir^T F C_occ + C_vir^T dF C_occ + C_vir^T F dC_occ) with F the perturbed
//Fock matrix of the current orbitals. The occupied-occupied and virtual-virtual parts of
//the second-order orbital change leave the functional alone, so this is the exact Hessian.
occ::Vec SCF_wrapper::hessian_vector(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double lambda, const occ::Vec& v, const bool two_electron) {
	occ::qm::MolecularOrbitals& mo = scf.ctx.mo;
	const int n = nmo_, nb = mo.kind == occ::qm::SpinorbitalKind::Unrestricted ? 2 : 1;
	const double fac = nb == 1 ? 4.0 : 2.0, dfac = nb == 1 ? 1.0 : 0.5;
	occ::Mat dC = occ::Mat::Zero(mo.C.rows(), n), dD = occ::Mat::Zero(mo.D.rows(), n), dCo = occ::Mat::Zero(mo.Cocc.rows(), mo.Cocc.cols());
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
		dCo.block(row, 0, n, nocc) = dCocc;
	}
	occ::Mat dF;
	const occ::qm::HartreeFock* df = hessian_df_ ? hessian_df_.get() : !eri_ && scf.m_procedure.using_density_fitting() ? &scf.m_procedure : nullptr;
	if (two_electron && df) {
		//occ's fitted exchange reads Cocc, not D: swapping in dD would give the exchange of the
		//current density. G is linear in the density and dD = (D(X + tY) - D(X - tY)) / 2t for
		//X = Cocc, Y = dCocc in the density's own convention, so two builds; t keeps X +- tY of
		//one size, the difference then loses no digits
		const double y = dCo.norm();
		dF = occ::Mat::Zero(mo.D.rows(), n);
		if (y > 0) {
			const double t = mo.Cocc.norm() / y;
			occ::qm::MolecularOrbitals m = mo;
			m.Cocc = mo.Cocc + t * dCo;
			m.update_density_matrix();
			dF = df->compute_fock(m);
			m.Cocc = mo.Cocc - t * dCo;
			m.update_density_matrix();
			dF -= df->compute_fock(m);
			dF /= 2 * t;
		}
		builds_.hessian++;
	}
	else if (two_electron) {
		//G(dD): the Fock build reads mo.D, and a screen on the density would drop the small dD
		std::swap(mo.D, dD);
		dF = eri_ ? eri_fock(mo, false) : scf.m_procedure.compute_fock(mo, scf.ctx.K);
		builds_.hessian++;
		std::swap(mo.D, dD);
	}
	else dF = occ::Mat::Zero(mo.D.rows(), n);
	if (lambda != 0) {
		dMatrix2 dD_eff(n, n);
		build_effective_dm(scf, dD_eff, dD);
		const occ::Mat dP = coupling_->perturbation_response(dD_eff, lambda);
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

//(H - theta)^-1 rhs for the Davidson subspace: the diagonal, refined at lambda > 0 by up to trah_pcg_
//conjugate-gradient steps on H without its two-electron part, whose products cost a structure
//factor evaluation instead of a Fock build. That part is what the diagonal misses: chi^2 is
//stiff along few directions. Davidson takes the exact product of whatever comes back. TEMPO def2-SVP,
//lambda 0.01-0.05: 70-106 Fock builds per step on the diagonal alone, 21-40 with it.
occ::Vec SCF_wrapper::trah_precondition(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double lambda, const occ::Vec& rhs, const double theta) {
	const occ::Vec d = (soscf_hdiag_.array() - theta).max(1e-2).matrix();
	occ::Vec x = rhs.cwiseQuotient(d);
	if (lambda == 0) return x;
	occ::Vec r = rhs - (hessian_vector(scf, lambda, x, false) - theta * x);
	occ::Vec z = r.cwiseQuotient(d), p = z;
	double rz = r.dot(z);
	const double stop = trah_pcg_tol_ * rhs.norm();
	for (int k = 0; k < trah_pcg_ && r.norm() > stop; k++) {
		const occ::Vec Ap = hessian_vector(scf, lambda, p, false) - theta * p;
		const double pAp = p.dot(Ap);
		if (!(pAp > 0)) break;
		const double a = rz / pAp;
		x += a * p;
		r -= a * Ap;
		z = r.cwiseQuotient(d);
		const double rzn = r.dot(z);
		p = z + (rzn / rz) * p;
		rz = rzn;
		trah_pcg_total_++;
	}
	return x;
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
	dMatrix2 dm_eff(nmo_, nmo_);
	build_effective_dm(scf, dm_eff, scf.ctx.mo.D);
	coupling_->evaluate(dm_eff);
	occ::Mat perturbation;
	coupling_->perturbation(perturbation, scf.ctx.mo.kind == occ::qm::SpinorbitalKind::Unrestricted);
	builds_.rebuild++;
	scf.ctx.F = scf.ctx.H + (eri_ ? eri_fock(scf.ctx.mo, false) : scf.m_procedure.compute_fock(scf.ctx.mo, scf.ctx.K));
	scf.update_scf_energy(false);
	const double crit = coupling_->quality().criterion;
	scf.ctx.F += lambda * perturbation;
	return scf.ctx.energy["total"] + lambda * crit * crit;
}

//orbital_rotation_gradient at C_from exp(kappa), rebuilt there and put back afterwards. Only
//the Hessian check pays for it.
occ::Vec SCF_wrapper::gradient_at(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double lambda, const occ::Mat& C_from, const occ::Vec& kappa) {
	const occ::qm::MolecularOrbitals mo_saved = scf.ctx.mo;
	const occ::Mat F_saved = scf.ctx.F;
	const auto energy_saved = scf.ctx.energy;
	coupling_->save_state();
	rebuild_at(scf, lambda, C_from, kappa);
	occ::Vec h;
	const occ::Vec g = orbital_rotation_gradient(scf, h);
	scf.ctx.mo = mo_saved;
	scf.ctx.F = F_saved;
	scf.ctx.energy = energy_saved;
	coupling_->restore_state();
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
	writer_->log << "\t\tHessian check: |Hv - FD| / |FD| = " << std::scientific << std::setprecision(2) << (Hv - fd).norm() / fd.norm()
		<< " (|Hv| " << Hv.norm() << ", |FD| " << fd.norm() << ", v.Hv " << v.dot(Hv) << ", v.FD " << v.dot(fd) << ")" << std::endl;
}

//The rescue of a lost SCF step, see best_mo_ in the header. True when it restarted the step:
//the orbitals are the best ones again, the iteration ends here and the next one takes their
//two-electron part from best_G_ instead of building it again. Each rescue doubles the level
//shift, halves the undamped share of the density and holds both ten times longer: restoring the start state with the start
//flags replayed the failure digit for digit (TEMPO, Cordula, three identical rescues). All
//of it is this lambda step's only, do_SCF resets it.
bool SCF_wrapper::rescue_scf(occ::qm::SCF<occ::qm::HartreeFock>& scf, const double quant, double& alpha) {
	if (quant < best_quant_) {
		best_quant_ = quant;
		best_mo_ = scf.ctx.mo;
		best_G_ = G_last_;
		return false;
	}
	if (quant - best_quant_ < rescue_rise_ || rescues_ >= 3) return false;
	rescues_++;
	const double tighten = std::pow(10.0, rescues_);
	settings.diis_stop_shift = (base_.stop_shift > 0 ? base_.stop_shift : 1e-2) / tighten;
	settings.diis_stop_damping = (base_.stop_damping > 0 ? base_.stop_damping : 1e-3) / tighten;
	settings.level_shift = std::max(base_.level_shift, 0.25) * std::pow(2.0, rescues_);
	alpha = 1.0 - (1.0 - settings.alpha) / std::pow(2.0, rescues_);
	std::ostringstream what1, what2, what3;
	what1 << "***E + lambda chi^2 " << std::fixed << std::setprecision(3) << quant - best_quant_
		<< " Eh above its best: back to the best orbitals";
	what2 << "level shift " << std::setprecision(2) << settings.level_shift
		<< ", damping " << alpha << " until DIIS error "
		<< std::scientific << std::setprecision(0) << settings.diis_stop_shift << " / " << settings.diis_stop_damping;
	what3 << "(rescue " << rescues_ << "/3)***";
	print_centered_message(what1.str(), 100, writer_->log);
	print_centered_message(what2.str(), 100, writer_->log);
	print_centered_message(what3.str(), 100, writer_->log);
	scf.ctx.mo = best_mo_;
	G_last_ = best_G_;
	reuse_best_G_ = true;
	diis_F_.clear();
	diis_E_.clear();
	adiis_.reset();
	ediis_.reset();
	settings.apply_shift = true;
	settings.apply_damping = true;
	return true;
}

bool SCF_wrapper::SCF_convergence_check(occ::qm::SCF<occ::qm::HartreeFock>& scf, occ::Mat& dm_last) {
	get_density_criteria(settings.current_RMSP_diff, settings.current_MaxP_diff, scf.ctx.mo.D, dm_last);
	settings.conv_quant_diff = settings.current_quant_diff < settings.quant_diff;
	settings.conv_max_diis_error = settings.current_max_diis_error < settings.max_diis_error;
	settings.conv_gradient = settings.current_gradient < settings.gradient;
	settings.conv_RMSP_diff = settings.current_RMSP_diff < settings.RMSP_diff;
	settings.conv_MaxP_diff = settings.current_MaxP_diff < settings.MaxP_diff;
	return settings.convergence_check();
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
