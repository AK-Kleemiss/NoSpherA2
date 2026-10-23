#include "pch.h"
#include "stored_eri.h"

namespace {
	//libcint mallocs scratch per quartet unless handed one; the size query is the same call without output
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

	//FNV-1a over what the integrals depend on: each shell's kind, l, centre, exponents and coefficients
	uint64_t basis_key(const occ::qm::AOBasis& basis) {
		uint64_t h = 1469598103934665603ull;
		const auto mix = [&h](const void* p, const size_t n) {
			const unsigned char* c = static_cast<const unsigned char*>(p);
			for (size_t i = 0; i < n; i++) h = (h ^ c[i]) * 1099511628211ull;
		};
		for (const occ::qm::Shell& s : basis.shells()) {
			const int kl[2] = { static_cast<int>(s.kind), static_cast<int>(s.l) };
			mix(kl, sizeof(kl));
			mix(s.origin.data(), 3 * sizeof(double));
			mix(s.exponents.data(), s.exponents.size() * sizeof(double));
			mix(s.contraction_coefficients.data(), s.contraction_coefficients.size() * sizeof(double));
		}
		const double t = stored_eri::threshold;
		mix(&t, sizeof(t));
		return h;
	}

	constexpr char eri_magic[8] = { 'N', 'O', 'S', 'E', 'R', 'I', '0', '1' };

	template <typename T> void put(std::ofstream& f, const T* p, const size_t n) {
		f.write(reinterpret_cast<const char*>(p), static_cast<std::streamsize>(n * sizeof(T)));
	}
	template <typename T> void get(std::ifstream& f, T* p, const size_t n) {
		f.read(reinterpret_cast<char*>(p), static_cast<std::streamsize>(n * sizeof(T)));
	}
}

void stored_eri::clear() {
	v_.reset();
	nbf_ = 0;
	pa_.clear(), pb_.clear(), first_.clear(), idx_.clear(), q_.clear(), qseg_.clear();
	qmax_ = 0.0;
}

bool stored_eri::save(const std::filesystem::path& path, const occ::qm::HartreeFock& hf) const {
	if (!v_) return false;
	std::ofstream f(path, std::ios::binary);
	const uint64_t key = basis_key(hf.aobasis());
	const int npk = npairs();
	f.write(eri_magic, sizeof(eri_magic));
	put(f, &key, 1), put(f, &nbf_, 1), put(f, &npk, 1), put(f, &qmax_, 1);
	put(f, pa_.data(), npk), put(f, pb_.data(), npk), put(f, first_.data(), first_.size()), put(f, idx_.data(), idx_.size());
	put(f, q_.data(), npk), put(f, qseg_.data(), qseg_.size()), put(f, v_.get(), bytes() / sizeof(double));
	return static_cast<bool>(f);
}

bool stored_eri::load(const std::filesystem::path& path, const occ::qm::HartreeFock& hf, const size_t budget_bytes, std::ostream& log) {
	clear();
	if (!hf.fock_build_properties().density_screened) return false;
	if (!std::filesystem::exists(path)) {
		log << "XCW: " << path.string() << " not found, building the two-electron integrals" << std::endl;
		return false;
	}
	std::ifstream f(path, std::ios::binary);
	char magic[sizeof(eri_magic)] = {};
	uint64_t key = 0;
	int nbf = 0, npk = 0;
	f.read(magic, sizeof(magic));
	get(f, &key, 1), get(f, &nbf, 1), get(f, &npk, 1);
	if (!f || !std::equal(magic, magic + sizeof(magic), eri_magic) || key != basis_key(hf.aobasis()) || nbf != static_cast<int>(hf.aobasis().nbf())
		|| npk < 0) {
		log << "XCW: " << path.string() << " holds no two-electron integrals of this basis and geometry, building them" << std::endl;
		return false;
	}
	const size_t npair = (size_t)nbf * (nbf + 1) / 2, nint = (size_t)npk * (npk + 1) / 2;
	if (budget_bytes != 0 && nint * sizeof(double) > budget_bytes) return false;
	get(f, &qmax_, 1);
	pa_.resize(npk), pb_.resize(npk), first_.resize(nbf + 1), idx_.resize(npair), q_.resize(npk), qseg_.resize(nbf);
	get(f, pa_.data(), npk), get(f, pb_.data(), npk), get(f, first_.data(), first_.size()), get(f, idx_.data(), npair);
	get(f, q_.data(), npk), get(f, qseg_.data(), nbf);
	v_.reset(new double[nint]);
	get(f, v_.get(), nint);
	if (!f || f.peek() != std::char_traits<char>::eof()) {
		clear();
		log << "XCW: " << path.string() << " is cut short or too long, building the two-electron integrals" << std::endl;
		return false;
	}
	nbf_ = nbf;
	log << "XCW: " << nint * sizeof(double) / 1048576.0 << " MB of two-electron integrals over " << npk << " of " << npair
		<< " basis-function pairs read from " << path.string() << std::endl;
	return true;
}

//Diagonal quartets (pq|pq) give each pair's Schwarz bound and so the kept set; then each slot is
//written by the one unique quartet of kept pairs that holds it
bool stored_eri::build(const occ::qm::HartreeFock& hf, const size_t budget_bytes, std::ostream& log) {
	clear();
	// Returns false for density-fitted calculations
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
		//r runs to p, so the work per pq grows along the loop: static chunks left the
		//last thread with most of it. Every
		//integral is written by exactly one pq, so the schedule does not touch the result
#pragma omp for schedule(dynamic, 8)
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

//With off-diagonal D pairs doubled, J is the packed symmetric matrix times that vector, a dot and an
//axpy per segment. K's four scatters per integral run along segment c (second index pair_b()): two
//dots against columns of D, two axpys into columns of the symmetrised K, directly on d = 0..c, else
//through the index list. The diagonal cd == ab has half weight and is taken back with its segment.
//With two densities (UHF) one pass over the integrals gives J of their sum and K of each, so the
//bandwidth-bound stream is read once per Fock build, not once per spin.
//Per-thread partials merge in thread order, so the result is schedule-independent.
void stored_eri::jk(const occ::Mat* const* D, const int nd, occ::Mat& J, occ::Mat* K, const bool screen) const {
	const int n = nbf_, npk = npairs(), nthr = omp_get_max_threads();
	const int *pa = pa_.data(), *pb = pb_.data(), *first = first_.data();
	const double* Dd[2] = { D[0]->data(), D[nd - 1]->data() };
	std::vector<double> dp(npk), fp(npk);
	for (int k = 0; k < npk; k++) {
		fp[k] = pa[k] == pb[k] ? 1.0 : 2.0;
		double s = 0.0;
		for (int i = 0; i < nd; i++) s += Dd[i][pa[k] + pb[k] * n];
		dp[k] = fp[k] * s;
	}
	//Largest element of every row of the difference (over both spins), and overall
	std::vector<double> rm;
	double rmax = 0.0;
	if (screen) {
		rm.assign(n, 0.0);
		for (int i = 0; i < nd; i++)
			for (int a = 0; a < n; a++) rm[a] = std::max(rm[a], D[i]->row(a).cwiseAbs().maxCoeff());
		for (int a = 0; a < n; a++) rmax = std::max(rmax, rm[a]);
	}
	std::vector<std::vector<double>> Jp(nthr);
	std::vector<occ::Mat> Kp(nthr * nd);
#pragma omp parallel
	{
		std::vector<double> Jl(npk, 0.0), t(n);
		occ::Mat Kl[2];
		double* Kd[2];
		for (int i = 0; i < nd; i++) Kl[i] = occ::Mat::Zero(n, n), Kd[i] = Kl[i].data();
		//c runs to a: cyclic rows balance the triangle and, unlike a dynamic schedule,
		//keep the per-thread partials and so the summation order fixed for a thread count
#pragma omp for schedule(static, 1)
		for (int k = 0; k < npk; k++) {
			const int a = pa[k], b = pb[k];
			if (screen && q_[k] * qmax_ * rmax < threshold) continue;
			const double* v = v_.get() + (size_t)k * (k + 1) / 2;
			const double fab = 2.0 * fp[k], dab = dp[k];
			for (int c = 0; c <= a; c++) {
				const int o = first[c], m = (c < a ? first[c + 1] : k + 1) - o;
				if (m <= 0) continue;
				if (screen && q_[k] * qseg_[c] * std::max({ rm[a], rm[b], rm[c] }) < threshold) continue;
				const double* vc = v + o;
				const Eigen::Map<const Eigen::VectorXd> vm(vc, m);
				//The row's own slot enters J once from the column side and once from the row side
				Eigen::Map<Eigen::VectorXd>(Jl.data() + o, c < a ? m : m - 1) += dab * vm.head(c < a ? m : m - 1);
				Jl[k] += vm.dot(Eigen::Map<const Eigen::VectorXd>(dp.data() + o, m));
				Eigen::Map<Eigen::VectorXd> tm(t.data(), m);
				tm = vm.cwiseProduct(Eigen::Map<const Eigen::VectorXd>(fp.data() + o, m));
				const bool direct = pb[o] == 0 && pb[o + m - 1] == m - 1;
				for (int s = 0; s < nd; s++) {
					const double *Da = Dd[s] + a * n, *Db = Dd[s] + b * n;
					double *Ks = Kd[s], *Ka = Ks + a * n, *Kb = Ks + b * n;
					const double kca = fab * Da[c], kcb = fab * Db[c];
					if (direct) {
						Ks[c + a * n] += fab * tm.dot(Eigen::Map<const Eigen::VectorXd>(Db, m));
						Ks[c + b * n] += fab * tm.dot(Eigen::Map<const Eigen::VectorXd>(Da, m));
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
						Ks[c + a * n] += fab * sa, Ks[c + b * n] += fab * sb;
					}
					if (c == a) {
						const double h = fp[k] * fp[k] * v[k];
						Ks[a + a * n] -= h * Db[b], Ks[a + b * n] -= h * Da[b];
						Ka[b] -= h * Da[b], Kb[b] -= h * Da[a];
					}
				}
			}
		}
		const int me = omp_get_thread_num();
		Jp[me].swap(Jl);
		for (int i = 0; i < nd; i++) Kp[me * nd + i].swap(Kl[i]);
	}
	J = occ::Mat::Zero(n, n);
	for (int i = 0; i < nd; i++) K[i].resize(n, n);
	//The partials outweigh the matrices many times, so their sum is parallel too
#pragma omp parallel for schedule(static)
	for (int a = 0; a < n; a++) {
		for (int i = 0; i < nd; i++)
			for (int b = 0; b <= a; b++) {
				double sk = 0.0;
				for (int t = 0; t < nthr; t++)
					if (!Jp[t].empty()) sk += Kp[t * nd + i](a, b) + Kp[t * nd + i](b, a);
				K[i](a, b) = K[i](b, a) = 0.125 * sk;
			}
		for (int k = first[a]; k < first[a + 1]; k++) {
			double sj = 0.0;
			for (int t = 0; t < nthr; t++)
				if (!Jp[t].empty()) sj += Jp[t][k];
			J(a, pb[k]) = J(pb[k], a) = sj;
		}
	}
}

void stored_eri::JK(const occ::Mat& D, occ::Mat& J, occ::Mat& K, const bool screen) const {
	const occ::Mat* d[1] = { &D };
	jk(d, 1, J, &K, screen);
}

void stored_eri::JK(const occ::Mat& Da, const occ::Mat& Db, occ::Mat& J, occ::Mat& Ka, occ::Mat& Kb, const bool screen) const {
	const occ::Mat* d[2] = { &Da, &Db };
	occ::Mat K[2];
	jk(d, 2, J, K, screen);
	Ka.swap(K[0]), Kb.swap(K[1]);
}

occ::Mat stored_eri::fock(const occ::qm::MolecularOrbitals& mo, const bool screen, const jk_fn& jk) const {
	const auto JK_ = [&](const occ::Mat& D, occ::Mat& J, occ::Mat& K) { if (jk) jk(D, J, K); else JK(D, J, K, screen); };
	occ::Mat J, K;
	if (mo.kind == occ::qm::SpinorbitalKind::Restricted) {
		JK_(mo.D, J, K);
		return 2.0 * J - K;
	}
	occ::Mat F = occ::Mat::Zero(mo.D.rows(), mo.D.cols()), Kb;
	const occ::Mat Da = occ::qm::block::a(mo.D), Db = occ::qm::block::b(mo.D);
	if (jk) {
		occ::Mat Jb;
		jk(Da, J, K);
		jk(Db, Jb, Kb);
		J += Jb;
	}
	else
		JK(Da, Db, J, K, Kb, screen);
	occ::qm::block::a(F) = 2.0 * J - 2.0 * K;
	occ::qm::block::b(F) = 2.0 * J - 2.0 * Kb;
	return F;
}
