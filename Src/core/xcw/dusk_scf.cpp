#include "pch.h"
#include "dusk_scf.h"
#ifdef NOSPHERA2_USE_DUSK
#include "basis_set.h"

dusk::BasisInput dusk_occ::basis(const BasisSet& set, const std::vector<occ::core::Atom>& atoms) {
	dusk::BasisInput b;
	for (int a = 0; a < static_cast<int>(atoms.size()); a++) {
		const int Z = atoms[a].atomic_number;
		if (b.elements.count(Z)) continue;
		std::vector<dusk::BasisShell>& shells = b.elements[Z];
		const std::span<const SimplePrimitive> p = set[Z - 1];
		for (int i = 0; i < static_cast<int>(p.size()); i++) {
			if (i == 0 || p[i].shell != p[i - 1].shell)
				shells.push_back({ p[i].type, {}, {}, 1 });
			shells.back().exponents.push_back(p[i].exp);
			shells.back().coefficients.push_back(p[i].coefficient);
		}
	}
	return b;
}

dusk::Input dusk_occ::input(const occ::core::Molecule& mol, const BasisSet& orbital) {
	dusk::Input in;
	for (const occ::core::Atom& a : mol.atoms()) {
		in.atomic_numbers.push_back(a.atomic_number);
		in.positions.push_back({ a.x, a.y, a.z });
	}
	in.charge = mol.charge();
	in.multiplicity = mol.multiplicity();
	in.basis = basis(orbital, mol.atoms());
	return in;
}

dusk::Matrix dusk_occ::to_dusk(const occ::Mat& A) {
	dusk::Matrix M(A.rows(), A.cols());
	Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(M.data(), A.rows(), A.cols()) = A;
	return M;
}

occ::Mat dusk_occ::to_occ(const dusk::Matrix& A) {
	return Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(A.data(), A.rows, A.cols);
}

dusk::Wavefunction dusk_occ::to_dusk(const occ::qm::Wavefunction& wfn) {
	const occ::qm::MolecularOrbitals& mo = wfn.mo;
	const Eigen::Index n = static_cast<Eigen::Index>(mo.n_ao), m = mo.C.cols();
	dusk::Wavefunction w;
	w.nbf = n;
	w.nmo = m;
	w.nalpha = static_cast<int>(mo.n_alpha);
	w.nbeta = static_cast<int>(mo.n_beta);
	w.unrestricted = mo.kind == occ::qm::SpinorbitalKind::Unrestricted;
	for (int s = 0; s < (w.unrestricted ? 2 : 1); s++) {
		w.C[s] = to_dusk(mo.C.middleRows(s * n, n));
		w.eps[s].assign(mo.energies.data() + s * m, mo.energies.data() + (s + 1) * m);
	}
	return w;
}

void dusk_occ::to_occ(const dusk::Result& r, occ::qm::Wavefunction& wfn) {
	const dusk::Wavefunction& w = r.wfn;
	const Eigen::Index n = static_cast<Eigen::Index>(w.nbf), m = static_cast<Eigen::Index>(w.nmo);
	const int ns = w.unrestricted ? 2 : 1;
	err_checkf(m == n, "scf_backend dusk: the basis is linearly dependent, which occ's wavefunction cannot hold", std::cout);
	occ::qm::MolecularOrbitals& mo = wfn.mo;
	mo.kind = w.unrestricted ? occ::qm::SpinorbitalKind::Unrestricted : occ::qm::SpinorbitalKind::Restricted;
	mo.n_alpha = w.nalpha;
	mo.n_beta = w.nbeta;
	mo.n_ao = n;
	mo.C.resize(ns * n, m);
	mo.energies.resize(ns * m);
	for (int s = 0; s < ns; s++) {
		mo.C.middleRows(s * n, n) = to_occ(w.C[s]);
		for (Eigen::Index j = 0; j < m; j++) mo.energies(s * m + j) = w.eps[s][j];
	}
	mo.update_occupied_orbitals();
	mo.update_density_matrix();
	wfn.nbf = n;
	wfn.num_electrons = w.nalpha + w.nbeta;
	wfn.have_energies = true;
	wfn.energy.total = r.energy - r.components.at("perturbation");
	wfn.converged = r.converged;
}

dusk::Options dusk_occ::scf_options(const options::SCF_settings& s, std::ostream& log) {
	dusk::Options o;
	o.verbosity = 0;
	o.convergence = { s.quant_diff, s.max_diis_error, s.gradient, s.RMSP_diff, s.MaxP_diff, true, true };
	o.max_iterations = s.max_scf_iterations;
	o.reference = s.rohf ? dusk::Reference::RestrictedOpen
		: s.hf_type == occ::qm::SpinorbitalKind::Unrestricted ? dusk::Reference::Unrestricted : dusk::Reference::Restricted;
	//`nodiis`: one stored Fock matrix, Roothaan steps that dusk shortens to the minimum of the energy's quadratic model
	o.diis_size = s.use_diis ? s.diis_size : 1;
	o.soscf = s.use_soscf;
	o.second_order = s.use_trah ? dusk::SecondOrder::Auto : dusk::SecondOrder::Off;
	//`nososcf`: TRAH takes SOSCF's place at the same DIIS error
	if (!s.use_soscf && s.use_trah) o.second_order_start = o.soscf_start;
	o.rescues = 3;
	o.confirm_convergence = true;
	o.stability_check = s.stability_check;
	o.broken_symmetry = s.broken_symmetry;
	o.check_hessian = s.check_hessian;
	o.damping = s.method_apply_damping ? s.alpha_damp : 0.0;
	o.damping_off = s.diis_stop_damping;
	o.level_shift = s.method_apply_shift ? s.level_shift : 0.0;
	o.level_shift_off = s.diis_stop_shift;
	o.incremental_fock = s.incremental;
	o.full_rebuild_drop = 10;
	if (!s.guess_basis_name.empty()) {
		o.guess = dusk::Guess::SmallBasis;
		o.guess_basis = s.guess_basis_name;
	}
	if (s.rohf) log << "dusk: rohf is HIGHLY EXPERIMENTAL, not for use yet (dusk::Reference::RestrictedOpen)\n";
	return o;
}

double xcw_perturbation::evaluate(const dusk::Matrix& D, dusk::Matrix& F) {
	const int n = static_cast<int>(D.rows);
	dMatrix2 dm(n, n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) dm(i, j) = 0.5 * D(i, j);
	coupling_->evaluate_cryst_crit(dm, incremental);
	fit = coupling_->cryst_crit();
	evaluations++;
	if (lambda != 0) {
		occ::Mat P;
		coupling_->perturbation_matrix(P, false, lambda, incremental);
		F = dusk_occ::to_dusk(P);
	}
	return fit.criterion;
}

bool xcw_perturbation::response(const dusk::Matrix& D, const dusk::Matrix& dD, dusk::Matrix& F1) {
	const int n = static_cast<int>(dD.rows);
	dMatrix2 dm(n, n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) dm(i, j) = 0.5 * dD(i, j);
	F1 = dusk_occ::to_dusk(coupling_->perturbation_matrix_response(dm, 1.0));
	return true;
}

//the input with the criterion as its first term (index 0), for the member initialiser of scf
static dusk::Input with_term(dusk::Input in, const std::shared_ptr<xcw_perturbation>& xcw) {
	in.perturbations.insert(in.perturbations.begin(), { xcw, 0.0 });
	return in;
}

dusk_solver::dusk_solver(dusk::Input in, const occ::core::Molecule& mol, const occ::qm::AOBasis& occ_basis, XCW_log_writer& writer, SCF_coupling& coupling)
	: xcw(std::make_shared<xcw_perturbation>(coupling)), scf(dusk::make_scf(with_term(std::move(in), xcw))),
	basis_(occ_basis), atoms_(mol.atoms()), writer_(&writer) {
	//XCW_solver::log_iteration's row: the fit of the density this Fock build evaluated, E_total without lambda * criterion
	scf.set_iteration_callback([this](const dusk::IterationInfo& i) {
		const fit_quality& f = xcw->fit;
		const double penalty = xcw->lambda * i.perturbation_energies[0];
		writer_->v = { i.iteration, xcw->lambda, f.criterion, f.GooF2, f.R1, i.energy - penalty, penalty, i.energy, f.criterion_all, f.R1_all, 0.0, false, f.wR2, f.wR2_all };
		writer_->iteration_line();
		});
}

bool dusk_solver::run(const double lambda, occ::qm::Wavefunction& wfn, const bool use_guess) {
	writer_->start_lambda(lambda);
	scf.set_strength(index_, lambda);
	if (use_guess) scf.set_guess(dusk_occ::to_dusk(wfn));
	const int evaluations = xcw->evaluations;
	const dusk::Result r = scf.run();
	wfn.basis = basis_;
	wfn.atoms = atoms_;
	dusk_occ::to_occ(r, wfn);
	std::ostringstream s;
	s << "***SCF " << (r.converged ? "converged" : "did not converge") << " in " << r.iterations << " iterations (dusk)***";
	print_centered_message(s.str(), 100, writer_->log);
	writer_->log << "dusk: E " << std::fixed << std::setprecision(10) << wfn.energy.total << " Eh, lambda criterion " << r.components.at("perturbation")
		<< ", criterion " << xcw->fit.criterion << ", " << xcw->evaluations - evaluations << " perturbation evaluations, "
		<< r.fock_builds << " Fock builds, " << r.response_builds << " response builds";
	if (!std::isnan(r.lowest_hessian_eigenvalue)) writer_->log << ", lowest Hessian eigenvalue " << std::scientific << r.lowest_hessian_eigenvalue;
	writer_->log << std::defaultfloat << std::endl;
	writer_->v.lambda = lambda;
	writer_->v.E_total = wfn.energy.total;
	writer_->v.quant = r.energy;
	return r.converged;
}

void dusk_solver::fd_check(const occ::qm::Wavefunction& wfn, std::ostream& out) {
	const Eigen::Index n = static_cast<Eigen::Index>(wfn.mo.n_ao);
	const occ::Mat Dt = 2.0 * (wfn.mo.kind == occ::qm::SpinorbitalKind::Unrestricted ? occ::Mat(wfn.mo.D.topRows(n) + wfn.mo.D.bottomRows(n)) : wfn.mo.D);
	const dusk::Matrix D = dusk_occ::to_dusk(Dt);
	const double keep = xcw->lambda, h = 1e-4;
	xcw->incremental = false;
	xcw->lambda = keep == 0 ? 1.0 : keep;
	std::mt19937 gen(7);
	std::normal_distribution<double> nd;
	out << std::scientific << std::setprecision(3);
	for (int k = 0; k < 4; k++) {
		occ::Mat X(n, n);
		for (Eigen::Index i = 0; i < n; i++)
			for (Eigen::Index j = 0; j <= i; j++) X(i, j) = X(j, i) = nd(gen);
		X /= X.cwiseAbs().maxCoeff();
		dusk::Matrix F0(n, n), F1(n, n), Fp(n, n), Fm(n, n);
		xcw->evaluate(D, F0);
		xcw->response(D, dusk_occ::to_dusk(X), F1);
		const double Ep = xcw->evaluate(dusk_occ::to_dusk(Dt + h * X), Fp), Em = xcw->evaluate(dusk_occ::to_dusk(Dt - h * X), Fm);
		const double fd = (Ep - Em) / (2 * h), an = (dusk_occ::to_occ(F0).array() * X.array()).sum();
		const occ::Mat dF = (dusk_occ::to_occ(Fp) - dusk_occ::to_occ(Fm)) / (2 * h);
		out << "dusk FD check " << k << ": dE_p/dD.X " << an << " central difference " << fd << " rel " << std::abs(an - fd) / std::abs(fd)
			<< "; response max|F1 - dF| / max|dF| " << (dusk_occ::to_occ(F1) - dF).cwiseAbs().maxCoeff() / dF.cwiseAbs().maxCoeff() << std::endl;
	}
	dusk::Matrix F(n, n);
	xcw->evaluate(D, F);
	xcw->lambda = keep;
	out << std::defaultfloat;
}
#endif
