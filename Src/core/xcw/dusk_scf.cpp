#include "pch.h"
#include "dusk_scf.h"
#ifdef NOSPHERA2_USE_DUSK
#include "basis_set.h"

molscf::Molecule dusk::molecule(const occ::core::Molecule& mol) {
	molscf::Molecule m;
	for (const occ::core::Atom& a : mol.atoms())
		m.atoms.push_back({ a.atomic_number, { a.x, a.y, a.z } });
	m.charge = mol.charge();
	m.multiplicity = mol.multiplicity();
	return m;
}

gtoint::Basis dusk::basis(const BasisSet& set, const std::vector<occ::core::Atom>& atoms) {
	std::vector<gtoint::Shell> shells;
	for (int a = 0; a < static_cast<int>(atoms.size()); a++) {
		const std::span<const SimplePrimitive> p = set[atoms[a].atomic_number - 1];
		for (int i = 0; i < static_cast<int>(p.size()); i++) {
			if (i == 0 || p[i].shell != p[i - 1].shell)
				shells.push_back({ p[i].type, { atoms[a].x, atoms[a].y, atoms[a].z }, {}, {}, 1, a });
			shells.back().exponents.push_back(p[i].exp);
			shells.back().coefficients.push_back(p[i].coefficient);
		}
	}
	return gtoint::Basis(std::move(shells), gtoint::Angular::Spherical);
}

molscf::Matrix dusk::to_dusk(const occ::Mat& A) {
	molscf::Matrix M(A.rows(), A.cols());
	Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(M.data(), A.rows(), A.cols()) = A;
	return M;
}

occ::Mat dusk::to_occ(const molscf::Matrix& A) {
	return Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(A.data(), A.rows, A.cols);
}

molscf::Wavefunction dusk::to_dusk(const occ::qm::Wavefunction& wfn) {
	const occ::qm::MolecularOrbitals& mo = wfn.mo;
	const Eigen::Index n = static_cast<Eigen::Index>(mo.n_ao), m = mo.C.cols();
	molscf::Wavefunction w;
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

void dusk::to_occ(const molscf::Result& r, occ::qm::Wavefunction& wfn) {
	const molscf::Wavefunction& w = r.wfn;
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

molscf::Options dusk::scf_options(const options::SCF_settings& s) {
	molscf::Options o;
	o.reference = s.hf_type == occ::qm::SpinorbitalKind::Unrestricted ? molscf::Reference::Unrestricted : molscf::Reference::Restricted;
	o.convergence = { s.quant_diff, s.max_diis_error, s.gradient, s.RMSP_diff, s.MaxP_diff };
	o.max_iterations = s.max_scf_iterations;
	if (!s.guess_basis_name.empty()) {
		o.guess = molscf::Guess::SmallBasis;
		o.guess_basis = s.guess_basis_name;
	}
	o.verbosity = 0;
	return o;
}

double xcw_perturbation::evaluate(const molscf::Matrix& D, molscf::Matrix& F) {
	const int n = static_cast<int>(D.rows);
	dMatrix2 dm(n, n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) dm(i, j) = 0.5 * D(i, j);
	coupling_->evaluate_cryst_crit(dm, false);
	fit = coupling_->cryst_crit();
	evaluations++;
	if (lambda != 0) {
		occ::Mat P;
		coupling_->perturbation_matrix(P, false, lambda, false);
		F = dusk::to_dusk(P);
	}
	return fit.criterion;
}

bool xcw_perturbation::response(const molscf::Matrix& D, const molscf::Matrix& dD, molscf::Matrix& F1) {
	const int n = static_cast<int>(dD.rows);
	dMatrix2 dm(n, n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) dm(i, j) = 0.5 * dD(i, j);
	F1 = dusk::to_dusk(coupling_->perturbation_matrix_response(dm, 1.0));
	return true;
}

dusk_solver::dusk_solver(const occ::core::Molecule& mol, const occ::qm::AOBasis& occ_basis, const BasisSet& orbital, const molscf::Options& o,
	XCW_log_writer& writer, SCF_coupling& coupling)
	: scf(dusk::molecule(mol), dusk::basis(orbital, mol.atoms()), o), xcw(std::make_shared<xcw_perturbation>(coupling)),
	basis_(occ_basis), atoms_(mol.atoms()), writer_(&writer) {
	index_ = scf.add_perturbation(xcw, 0.0);
}

bool dusk_solver::run(const double lambda, occ::qm::Wavefunction& wfn, const bool use_guess) {
	writer_->start_lambda(lambda);
	xcw->lambda = lambda;
	scf.set_strength(index_, lambda);
	if (use_guess) scf.set_guess(dusk::to_dusk(wfn));
	const int evaluations = xcw->evaluations;
	const molscf::Result r = scf.run();
	wfn.basis = basis_;
	wfn.atoms = atoms_;
	dusk::to_occ(r, wfn);
	std::ostringstream s;
	s << "***SCF " << (r.converged ? "converged" : "did not converge") << " in " << r.iterations << " iterations (molscf)***";
	print_centered_message(s.str(), 100, writer_->log);
	writer_->log << "molscf: E " << std::fixed << std::setprecision(10) << wfn.energy.total << " Eh, lambda criterion " << r.components.at("perturbation")
		<< ", criterion " << xcw->fit.criterion << ", " << xcw->evaluations - evaluations << " perturbation evaluations";
	if (r.components.count("lowest_hessian_eigenvalue")) writer_->log << ", lowest Hessian eigenvalue " << std::scientific << r.components.at("lowest_hessian_eigenvalue");
	writer_->log << std::defaultfloat << std::endl;
	writer_->v.lambda = lambda;
	writer_->v.E_total = wfn.energy.total;
	writer_->v.quant = r.energy;
	return r.converged;
}

void dusk_solver::fd_check(const occ::qm::Wavefunction& wfn, std::ostream& out) {
	const Eigen::Index n = static_cast<Eigen::Index>(wfn.mo.n_ao);
	const occ::Mat Dt = 2.0 * (wfn.mo.kind == occ::qm::SpinorbitalKind::Unrestricted ? occ::Mat(wfn.mo.D.topRows(n) + wfn.mo.D.bottomRows(n)) : wfn.mo.D);
	const molscf::Matrix D = dusk::to_dusk(Dt);
	const double keep = xcw->lambda, h = 1e-4;
	xcw->lambda = keep == 0 ? 1.0 : keep;
	std::mt19937 gen(7);
	std::normal_distribution<double> nd;
	out << std::scientific << std::setprecision(3);
	for (int k = 0; k < 4; k++) {
		occ::Mat X(n, n);
		for (Eigen::Index i = 0; i < n; i++)
			for (Eigen::Index j = 0; j <= i; j++) X(i, j) = X(j, i) = nd(gen);
		X /= X.cwiseAbs().maxCoeff();
		molscf::Matrix F0(n, n), F1(n, n), Fp(n, n), Fm(n, n);
		xcw->evaluate(D, F0);
		xcw->response(D, dusk::to_dusk(X), F1);
		const double Ep = xcw->evaluate(dusk::to_dusk(Dt + h * X), Fp), Em = xcw->evaluate(dusk::to_dusk(Dt - h * X), Fm);
		const double fd = (Ep - Em) / (2 * h), an = (dusk::to_occ(F0).array() * X.array()).sum();
		const occ::Mat dF = (dusk::to_occ(Fp) - dusk::to_occ(Fm)) / (2 * h);
		out << "dusk FD check " << k << ": dE_p/dD.X " << an << " central difference " << fd << " rel " << std::abs(an - fd) / std::abs(fd)
			<< "; response max|F1 - dF| / max|dF| " << (dusk::to_occ(F1) - dF).cwiseAbs().maxCoeff() / dF.cwiseAbs().maxCoeff() << std::endl;
	}
	molscf::Matrix F(n, n);
	xcw->evaluate(D, F);
	xcw->lambda = keep;
	out << std::defaultfloat;
}
#endif
