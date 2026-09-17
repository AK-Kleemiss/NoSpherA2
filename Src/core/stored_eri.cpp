#include "pch.h"
#include "stored_eri.h"

namespace {
	//libcint mallocs its scratch per quartet unless handed one; the size query is the same call
	//without an output
	std::array<int, 4> quartet(occ::qm::cint::IntegralEnvironment& env, const bool sph, std::array<int, 4> sh,
		occ::qm::cint::Optimizer& opt, std::vector<double>& buffer, std::vector<double>& cache) {
		const size_t need = sph
			? libcint::int2e_sph(nullptr, nullptr, sh.data(), env.atom_data_ptr(), env.num_atoms(), env.basis_data_ptr(), env.num_basis(), env.env_data_ptr(), nullptr, nullptr)
			: libcint::int2e_cart(nullptr, nullptr, sh.data(), env.atom_data_ptr(), env.num_atoms(), env.basis_data_ptr(), env.num_basis(), env.env_data_ptr(), nullptr, nullptr);
		if (need > cache.size()) cache.resize(need);
		return sph
			? env.four_center_helper<occ::qm::cint::Operator::coulomb, occ::qm::Shell::Kind::Spherical>(sh, opt.optimizer_ptr(), buffer.data(), cache.data())
			: env.four_center_helper<occ::qm::cint::Operator::coulomb, occ::qm::Shell::Kind::Cartesian>(sh, opt.optimizer_ptr(), buffer.data(), cache.data());
	}

	inline size_t packed(const int a, const int b) {
		return a >= b ? (size_t)a * (a + 1) / 2 + b : (size_t)b * (b + 1) / 2 + a;
	}
}

void stored_eri::clear() {
	v_.reset();
	nbf_ = 0;
	pa_.clear(), pb_.clear(), first_.clear(), idx_.clear(), q_.clear(), qseg_.clear();
	qmax_ = 0.0;
}

//Two passes over the significant shell pairs: the diagonal quartets (pq|pq) give the Schwarz
//bound of every basis-function pair, which decides the kept set, then every unique quartet
//of kept pairs is computed and each slot written by the one quartet that holds it
bool stored_eri::build(const occ::qm::HartreeFock& hf, const size_t budget_bytes, std::ostream& log) {
	clear();
	if (!hf.fock_build_properties().density_screened) return false;
	occ::qm::IntegralEngine engine(hf.aobasis());
	const int nbf = static_cast<int>(engine.nbf()), nsh = static_cast<int>(engine.nsh()), npq = nsh * (nsh + 1) / 2;
	const size_t npair = (size_t)nbf * (nbf + 1) / 2;
	const auto& shellpairs = engine.shellpairs();
	const auto& first_bf = engine.first_bf();
	const bool sph = engine.is_spherical();
	auto& env = engine.env();
	const auto significant = [&](const int p, const int q) {
		return std::binary_search(shellpairs[p].begin(), shellpairs[p].end(), (size_t)q);
	};
	//Bound of every pair, -1 where the shell pair is not significant
	std::vector<double> q(npair, -1.0);
#pragma omp parallel
	{
		occ::qm::IntegralEngine local_engine(hf.aobasis());
		auto& local_env = local_engine.env();
		occ::qm::cint::Optimizer opt(local_env, occ::qm::cint::Operator::coulomb, 4);
		std::vector<double> buffer(local_env.buffer_size_2e()), cache;
#pragma omp for schedule(static)
		for (int pq = 0; pq < npq; pq++) {
			const int p = static_cast<int>((std::sqrt(8.0 * pq + 1.0) - 1.0) / 2.0), qs = pq - p * (p + 1) / 2;
			if (!significant(p, qs)) continue;
			const std::array<int, 4> dims = quartet(local_env, sph, { p, qs, p, qs }, opt, buffer, cache);
			//dims[0] is -1 for an all-zero quartet; dims[2] is the same shell's size
			const int d0 = dims[2], d1 = dims[1];
			for (int f1 = 0; f1 < d1; f1++)
				for (int f0 = 0; f0 < d0; f0++)
					q[packed(first_bf[p] + f0, first_bf[qs] + f1)] = dims[0] < 0 ? 0.0
						: std::sqrt(std::abs(buffer[f0 + d0 * (f1 + d1 * (f0 + d0 * f1))]));
		}
	}
	qmax_ = *std::max_element(q.begin(), q.end());
	idx_.assign(npair, -1);
	first_.assign(nbf + 1, 0);
	for (int a = 0, ab = 0; a < nbf; a++) {
		first_[a] = npairs();
		for (int b = 0; b <= a; b++, ab++) {
			if (q[ab] < 0.0 || q[ab] * qmax_ < threshold) continue;
			idx_[ab] = npairs();
			pa_.push_back(a), pb_.push_back(b);
		}
	}
	const int npk = npairs();
	first_[nbf] = npk;
	const size_t nint = (size_t)npk * (npk + 1) / 2;
	if (budget_bytes != 0 && nint * sizeof(double) > budget_bytes) {
		log << "XCW: " << nint * sizeof(double) / 1048576.0 << " MB of two-electron integrals over " << npk << " of " << npair
			<< " pairs do not fit in memory, direct Fock builds" << std::endl;
		clear();
		return false;
	}
	nbf_ = nbf;
	v_.reset(new double[nint]);
#pragma omp parallel
	{
		occ::qm::IntegralEngine local_engine(hf.aobasis());
		auto& local_env = local_engine.env();
		occ::qm::cint::Optimizer opt(local_env, occ::qm::cint::Operator::coulomb, 4);
		std::vector<double> buffer(local_env.buffer_size_2e()), cache;
#pragma omp for schedule(static)
		for (long long i = 0; i < static_cast<long long>(nint); i++) v_[i] = 0.0;
#pragma omp for schedule(static)
		for (int pq = 0; pq < npq; pq++) {
			const int p = static_cast<int>((std::sqrt(8.0 * pq + 1.0) - 1.0) / 2.0), qs = pq - p * (p + 1) / 2;
			if (!significant(p, qs)) continue;
			for (int r = 0; r <= p; r++) {
				const int s_max = p == r ? qs : r;
				for (const size_t s : shellpairs[r]) {
					if (static_cast<int>(s) > s_max) break;
					const std::array<int, 4> dims = quartet(local_env, sph, { p, qs, r, static_cast<int>(s) }, opt, buffer, cache);
					if (dims[0] < 0) continue;
					const double* v = buffer.data();
					for (int f3 = 0; f3 < dims[3]; f3++) {
						const int d = first_bf[s] + f3;
						for (int f2 = 0; f2 < dims[2]; f2++) {
							const int cd = idx_[packed(first_bf[r] + f2, d)];
							for (int f1 = 0; f1 < dims[1]; f1++) {
								const int b = first_bf[qs] + f1;
								for (int f0 = 0; f0 < dims[0]; f0++, v++) {
									const int ab = idx_[packed(first_bf[p] + f0, b)];
									if (ab < 0 || cd < 0) continue;
									v_[ab >= cd ? (size_t)ab * (ab + 1) / 2 + cd : (size_t)cd * (cd + 1) / 2 + ab] = *v;
								}
							}
						}
					}
				}
			}
		}
	}
	q_.resize(npk);
	for (int k = 0; k < npk; k++) q_[k] = std::sqrt(std::abs(v_[(size_t)k * (k + 1) / 2 + k]));
	qseg_.assign(nbf, 0.0);
	for (int c = 0; c < nbf; c++)
		for (int k = first_[c]; k < first_[c + 1]; k++) qseg_[c] = std::max(qseg_[c], q_[k]);
	log << "XCW: " << nint * sizeof(double) / 1048576.0 << " MB of two-electron integrals over " << npk << " of " << npair
		<< " basis-function pairs held in memory" << std::endl;
	return true;
}

//With the off-diagonal pairs of D doubled J is the symmetric packed matrix times that vector,
//one dot and one axpy per segment. K's four scatters per integral run along the segment c of a
//row, where the second index is pair_b(): two dots against columns of the symmetric D and two
//axpys into columns of the (symmetrised) K, on the run d = 0..c directly and through the index
//list otherwise. The diagonal cd == ab carries half the weight and is taken back with its
//segment. Per-thread partials merged in thread order, so the result does not depend on the
//schedule.
void stored_eri::JK(const occ::Mat& D, occ::Mat& J, occ::Mat& K, const bool screen) const {
	const int n = nbf_, npk = npairs(), nthr = omp_get_max_threads();
	const int *pa = pa_.data(), *pb = pb_.data(), *first = first_.data();
	std::vector<double> dp(npk), fp(npk);
	for (int k = 0; k < npk; k++) fp[k] = pa[k] == pb[k] ? 1.0 : 2.0, dp[k] = fp[k] * D(pa[k], pb[k]);
	//Largest element of every row of the difference, and overall
	std::vector<double> rm;
	double rmax = 0.0;
	if (screen) {
		rm.resize(n);
		for (int a = 0; a < n; a++) rmax = std::max(rmax, rm[a] = D.row(a).cwiseAbs().maxCoeff());
	}
	const double* Dd = D.data();
	std::vector<std::vector<double>> Jp(nthr);
	std::vector<occ::Mat> Kp(nthr);
#pragma omp parallel
	{
		std::vector<double> Jl(npk, 0.0), t(n);
		occ::Mat Kl = occ::Mat::Zero(n, n);
		double* Kd = Kl.data();
#pragma omp for schedule(static)
		for (int k = 0; k < npk; k++) {
			const int a = pa[k], b = pb[k];
			if (screen && q_[k] * qmax_ * rmax < threshold) continue;
			const double* v = v_.get() + (size_t)k * (k + 1) / 2;
			const double fab = 2.0 * fp[k], dab = fp[k] * Dd[a + b * n];
			const double *Da = Dd + a * n, *Db = Dd + b * n;
			double *Ka = Kd + a * n, *Kb = Kd + b * n;
			for (int c = 0; c <= a; c++) {
				const int o = first[c], m = (c < a ? first[c + 1] : k + 1) - o;
				if (m <= 0) continue;
				if (screen && q_[k] * qseg_[c] * std::max({ rm[a], rm[b], rm[c] }) < threshold) continue;
				const double* vc = v + o;
				const Eigen::Map<const Eigen::VectorXd> vm(vc, m);
				//The row's own slot enters J once from the column side and once from the row side
				Eigen::Map<Eigen::VectorXd>(Jl.data() + o, c < a ? m : m - 1) += dab * vm.head(c < a ? m : m - 1);
				Jl[k] += vm.dot(Eigen::Map<const Eigen::VectorXd>(dp.data() + o, m));
				const double kca = fab * Da[c], kcb = fab * Db[c];
				Eigen::Map<Eigen::VectorXd> tm(t.data(), m);
				tm = vm.cwiseProduct(Eigen::Map<const Eigen::VectorXd>(fp.data() + o, m));
				if (pb[o] == 0 && pb[o + m - 1] == m - 1) {
					Kd[c + a * n] += fab * tm.dot(Eigen::Map<const Eigen::VectorXd>(Db, m));
					Kd[c + b * n] += fab * tm.dot(Eigen::Map<const Eigen::VectorXd>(Da, m));
					Eigen::Map<Eigen::VectorXd>(Ka, m) += kcb * tm;
					Eigen::Map<Eigen::VectorXd>(Kb, m) += kca * tm;
				}
				else {
					double sa = 0.0, sb = 0.0;
					for (int i = 0; i < m; i++) {
						const int d = pb[o + i];
						sa += t[i] * Db[d], sb += t[i] * Da[d];
						Ka[d] += kcb * t[i], Kb[d] += kca * t[i];
					}
					Kd[c + a * n] += fab * sa, Kd[c + b * n] += fab * sb;
				}
				if (c == a) {
					const double h = fp[k] * fp[k] * v[k];
					Kd[a + a * n] -= h * Db[b], Kd[a + b * n] -= h * Da[b];
					Ka[b] -= h * Da[b], Kb[b] -= h * Da[a];
				}
			}
		}
		Jp[omp_get_thread_num()].swap(Jl);
		Kp[omp_get_thread_num()].swap(Kl);
	}
	J = occ::Mat::Zero(n, n);
	K.resize(n, n);
	//The partials outweigh the matrices many times over, so their sum is parallel too
#pragma omp parallel for schedule(static)
	for (int a = 0; a < n; a++) {
		for (int b = 0; b <= a; b++) {
			double sk = 0.0;
			for (int t = 0; t < nthr; t++)
				if (!Jp[t].empty()) sk += Kp[t](a, b) + Kp[t](b, a);
			K(a, b) = K(b, a) = 0.125 * sk;
		}
		for (int k = first[a]; k < first[a + 1]; k++) {
			double sj = 0.0;
			for (int t = 0; t < nthr; t++)
				if (!Jp[t].empty()) sj += Jp[t][k];
			J(a, pb[k]) = J(pb[k], a) = sj;
		}
	}
}

occ::Mat stored_eri::fock(const occ::qm::MolecularOrbitals& mo, const bool screen, const jk_fn& jk) const {
	const auto JK_ = [&](const occ::Mat& D, occ::Mat& J, occ::Mat& K) { if (jk) jk(D, J, K); else JK(D, J, K, screen); };
	occ::Mat J, K;
	if (mo.kind == occ::qm::SpinorbitalKind::Restricted) {
		JK_(mo.D, J, K);
		return 2.0 * J - K;
	}
	occ::Mat F = occ::Mat::Zero(mo.D.rows(), mo.D.cols()), Jb, Kb;
	JK_(occ::qm::block::a(mo.D), J, K);
	JK_(occ::qm::block::b(mo.D), Jb, Kb);
	occ::qm::block::a(F) = (J + Jb) - K;
	occ::qm::block::b(F) = (J + Jb) - Kb;
	return F;
}
