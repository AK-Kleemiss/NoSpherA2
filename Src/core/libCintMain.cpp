#include "pch.h"
#include "constants.h"
#include "libCintMain.h"
#include "libCintKernels.h"
#include "tuning.h"
//
#if defined(__APPLE__)
// On macOS we are using Accelerate for BLAS/LAPACK
#include <Accelerate/Accelerate.h>
#elif !defined(NSA2_OPENBLAS)
// Linux/Windows with oneMKL
#include <mkl.h>
#endif

// Function to compute three-center two-electron integrals (eri3c)
template <typename Kernel>
void computeEri3c(Int_Params &param1,
	Int_Params &param2,
	vec &eri3c)
{
	int nQM = param1.get_nbas();
	int nAux = param2.get_nbas();

	Int_Params combined(param1, param2);
	// combined.print_data("combined");

	ivec bas = combined.get_bas();
	ivec atm = combined.get_atm();
	vec env = combined.get_env();

	ivec shl_slice = {
		0,
		nQM,
		0,
		nQM,
		nQM,
		nQM + nAux,
	};

	int nat = combined.get_natoms();
	int nbas = combined.get_nbas();

	assert(shl_slice[1] <= nbas);
	assert(shl_slice[3] <= nbas);
	assert(shl_slice[5] <= nbas);

	ivec aoloc = Kernel::gen_loc(bas, nbas);
	int naoi = aoloc[shl_slice[1]] - aoloc[shl_slice[0]];
	int naoj = aoloc[shl_slice[3]] - aoloc[shl_slice[2]];
	int naok = aoloc[shl_slice[5]] - aoloc[shl_slice[4]];

	libcint::CINTOpt* opty = nullptr;
	Kernel::optimizer(opty, atm.data(), nat, bas.data(), nbas, env.data());

	// Compute integrals
	vec res((size_t)naoi * (size_t)naoj * (size_t)naok, 0.0);
	eri3c.resize((size_t)naoi * (size_t)naoj * (size_t)naok, 0.0);

	Kernel::drv(res.data(),
		1,
		shl_slice.data(),
		aoloc.data(),
		opty,
		atm.data(), nat,
		bas.data(), nbas,
		env.data());

	// FOR TESTING PURPOSES!!!!
	// GTOnr3c_drv(int3c2e_sph, res.data(), 1, shl_slice.data(), aoloc.data(), NULL, atm.data(), nat, bas.data(), nbas, env.data());

	// res is in fortran order, write the result in regular ordering
	for (int k = 0; k < naok; k++)
	{
		for (int j = 0; j < naoj; j++)
		{
			for (int i = 0; i < naoi; i++)
			{
				std::size_t idx_F = i + j * (size_t)naoi + k * ((size_t)naoi * (size_t)naoj);
				std::size_t idx_C = i * ((size_t)naoj * (size_t)naok) + j * (size_t)naok + k;
				eri3c[idx_C] = res[idx_F];
			}
		}
	}
}
template void computeEri3c<Coulomb3C_SPH>(Int_Params &param1,
	Int_Params &param2,
	vec &eri3c);

template <typename Kernel>
void compute2C(Int_Params &params, vec &ret) {
	ivec bas = params.get_bas();
	ivec atm = params.get_atm();
	vec env = params.get_env();

	int nbas = params.get_nbas();
	int nat = params.get_natoms();

	ivec shl_slice = { 0, nbas, 0, nbas };
	ivec aoloc = Kernel::gen_loc(bas, nbas);

	int naoi = aoloc[shl_slice[1]] - aoloc[shl_slice[0]];
	int naoj = aoloc[shl_slice[3]] - aoloc[shl_slice[2]];

	libcint::CINTOpt* opty = nullptr;
	Kernel::optimizer(opty, atm.data(), nat, bas.data(), nbas, env.data());

	// Compute integrals
	vec res((size_t)naoi * (size_t)naoj, 0.0);
	ret.resize((size_t)naoi * (size_t)naoj, 0.0);
	Kernel::drv(res.data(), 1, shl_slice.data(), aoloc.data(), opty, atm.data(), nat, bas.data(), nbas, env.data());

	// res is in fortran order, write the result in regular ordering
	for (int i = 0; i < naoi; i++)
	{
		for (int j = 0; j < naoj; j++)
		{
			ret[(size_t)j * (size_t)naoi + i] = res[(size_t)i * (size_t)naoj + j];
		}
	}
}
template void compute2C<Coulomb2C_SPH>(Int_Params &params, vec &ret);
template void compute2C<Coulomb2C_CRT>(Int_Params &params, vec &ret);
template void compute2C<Overlap2C_SPH>(Int_Params &params, vec &ret);
template void compute2C<Overlap2C_CRT>(Int_Params &params, vec &ret);


extern "C" libcint::CINTIntegralFunction int3c1e_sph;

//Per-shell integrals behind each 3C kernel. pair and aux give the Schwarz factors of
//|(ab|P)| <= sqrt((ab|ab)) sqrt((P|P)); nullptr where the metric has none. The overlap bound would need
//int4c1e, which this libcint build may lack, so the overlap metric is not screened.
template <typename K> struct Shell3C;
//far_field: a pure Coulomb aux shell outside the orbital density sees it only through one local-expansion term (see computeRho);
//a Cartesian d shell mixes in r^2 s, an overlap kernel has no far field.
template <> struct Shell3C<Coulomb3C_SPH> {
	static constexpr libcint::CINTIntegralFunction *three = libcint::int3c2e_sph, *pair = libcint::int2e_sph, *aux = libcint::int2c2e_sph;
	static constexpr bool far_field = true;
};
template <> struct Shell3C<Coulomb3C_CRT> {
	static constexpr libcint::CINTIntegralFunction *three = libcint::int3c2e_cart, *pair = libcint::int2e_cart, *aux = libcint::int2c2e_cart;
	static constexpr bool far_field = false;
};
template <> struct Shell3C<Overlap3C_SPH> {
	static constexpr libcint::CINTIntegralFunction *three = int3c1e_sph, *pair = nullptr, *aux = nullptr;
	static constexpr bool far_field = false;
};

//rho_P = sum_ab w D_ab (ab|P) over orbital shell pairs a >= b (w = 2 off the diagonal), one shell triplet at a time,
//skipped when w max|D_ab| Q_ab Q_P sigma_P < thr. sigma_p = max_A |(J^-1 n_A)_p| is the population a unit error in
//rho_p moves onto atom A, so thr is in electrons; without sigma nothing is screened.
template <typename Kernel>
void computeRho(
	const Int_Params &normal_basis,
	const Int_Params &aux_basis,
	const dMatrix2 &dm,
	vec &rho,
	const std::optional<ivec> asym_atm_list,
	const vec *sensitivity)
{
	using S = Shell3C<Kernel>;
	Int_Params combined(normal_basis, aux_basis);

	ivec bas = combined.get_bas();
	ivec atm = combined.get_atm();
	vec  env = combined.get_env();

	const int nQM = normal_basis.get_nbas();
	const int nAux = aux_basis.get_nbas();
	const int nat = combined.get_natoms();
	const int nbas = combined.get_nbas();

	ivec aoloc = Kernel::gen_loc(bas, nbas);
	const int aux0 = aoloc[nQM];
	const int naux = aoloc[nQM + nAux] - aux0;
	rho.assign(naux, 0.0);
	double thr = constants::ri_screen_threshold;
	if (const char *e = tuning("NOS_RI_SCREEN")) thr = std::atof(e);
	const bool bound = S::pair != nullptr && thr > 0.0 && sensitivity != nullptr;
	int dmax_orb = 0, dmax_aux = 0;
	for (int s = 0; s < nbas; s++) {
		int &d = s < nQM ? dmax_orb : dmax_aux;
		d = std::max(d, aoloc[s + 1] - aoloc[s]);
	}
	//Far field: once every orbital primitive pair (exponent p, centre R_p) overlaps an aux shell on atom C by less than
	//erfc(x) < eps, x^2 = pq/(p+q) |R_p - C|^2 with q the shell's most diffuse exponent, the shell R(r) Y_lm acts as
	//4pi/(2l+1) m Y_lm/r_C^(l+1) with m = sum_i c_i int r^(2l+2) e^(-q_i r^2) dr ~ sum_i c_i q_i^-(l+3/2). Every far shell of
	//one l on C is then the same integral times its m: one single-primitive reference shell per (C, l), appended after
	//the aux shells at the tightest exponent (so it is far whenever any of them is), is computed and scaled by m/m_ref.
	double far_eps = constants::ri_far_erfc;
	if (const char *e = tuning("NOS_RI_FAR")) far_eps = std::atof(e);
	const bool use_far = S::far_field && far_eps > 0.0;
	double X2 = 0.0;
	ivec ref_of(nAux, -1);
	vec ratio(nAux, 0.0), qmin(nAux, 0.0);
	int nbx = nbas, maxprim = 1;
	if (use_far) {
		double lo = 0.0, hi = 30.0; //erfc(X) = eps
		for (int it = 0; it < 100; it++) {
			const double mid = 0.5 * (lo + hi);
			(std::erfc(mid) > far_eps ? lo : hi) = mid;
		}
		X2 = hi * hi;
		for (int s = 0; s < nQM; s++) maxprim = std::max(maxprim, bas[s * BAS_SLOTS + NPRIM_OF]);
		std::map<std::pair<int, int>, double> qref;
		for (int P = 0; P < nAux; P++) {
			const int *s = &bas[(nQM + P) * BAS_SLOTS];
			if (s[NCTR_OF] != 1) continue;
			const double *ex = env.data() + s[PTR_EXP];
			double &q = qref[{ s[ATOM_OF], s[ANG_OF] }];
			q = std::max(q, *std::max_element(ex, ex + s[NPRIM_OF]));
		}
		std::map<std::pair<int, int>, int> ref_shell;
		for (const auto &[key, q] : qref) {
			ref_shell[key] = nbx++;
			const int ptr = static_cast<int>(env.size());
			env.push_back(q);
			env.push_back(1.0);
			const int row[BAS_SLOTS] = { key.first, key.second, 1, 1, 0, ptr, ptr + 1, 0 };
			bas.insert(bas.end(), row, row + BAS_SLOTS);
		}
		for (int P = 0; P < nAux; P++) {
			const int *s = &bas[(nQM + P) * BAS_SLOTS];
			const double *ex = env.data() + s[PTR_EXP], *co = env.data() + s[PTR_COEFF];
			const int l = s[ANG_OF];
			qmin[P] = *std::min_element(ex, ex + s[NPRIM_OF]);
			if (s[NCTR_OF] != 1) continue;
			const std::pair<int, int> key{ s[ATOM_OF], l };
			double m = 0.0;
			for (int i = 0; i < s[NPRIM_OF]; i++) m += co[i] * std::pow(ex[i], -(l + 1.5));
			ref_of[P] = ref_shell[key];
			ratio[P] = m / std::pow(qref[key], -(l + 1.5));
		}
	}
	//libcint mallocs scratch per call unless handed one; the size query is the call without output.
	//{s,s,s,s} over every shell is the bound pyscf's GTOmax_cache_size uses
	CACHE_SIZE_T ncache = 0;
	for (int s = 0; s < nbx; s++) {
		int shls[4] = { s, s, s, s };
		ncache = std::max(ncache, S::three(nullptr, nullptr, shls, atm.data(), nat, bas.data(), nbx, env.data(), nullptr, nullptr));
		if (bound) ncache = std::max(ncache, s < nQM
			? S::pair(nullptr, nullptr, shls, atm.data(), nat, bas.data(), nQM, env.data(), nullptr, nullptr)
			: S::aux(nullptr, nullptr, shls, atm.data(), nat, bas.data(), nbas, env.data(), nullptr, nullptr));
	}
	//Q_P = sqrt(max_p (p|p)) max_p sigma_p per aux shell
	vec qaux(nAux, 1.0);
	double qaux_max = 1.0;
	if (bound) {
		vec buf(static_cast<size_t>(dmax_aux) * dmax_aux), cache(ncache);
		qaux_max = 0.0;
		for (int P = 0; P < nAux; P++) {
			int shls[2] = { nQM + P, nQM + P };
			const int dp = aoloc[nQM + P + 1] - aoloc[nQM + P];
			double m = 0.0, s = 0.0;
			if (S::aux(buf.data(), nullptr, shls, atm.data(), nat, bas.data(), nbas, env.data(), nullptr, cache.data()))
				for (int f = 0; f < dp; f++) m = std::max(m, std::abs(buf[f + dp * f]));
			for (int f = 0; f < dp; f++) s = std::max(s, std::abs((*sensitivity)[aoloc[nQM + P] - aux0 + f]));
			qaux[P] = std::sqrt(m) * s;
			qaux_max = std::max(qaux_max, qaux[P]);
		}
	}

	libcint::CINTOpt *opty = nullptr, *opt2 = nullptr;
	Kernel::optimizer(opty, atm.data(), nat, bas.data(), nbx, env.data());
	//(ab|ab) only touches orbital shells, the first nQM, so its optimizer skips the aux pairs
	if (bound) libcint::int2e_optimizer(&opt2, atm.data(), nat, bas.data(), nQM, env.data());

	//Every shell a sums into its own vector, added to rho in a fixed order once all before it are done: the
	//ill-conditioned metric turns summation-order noise into tsc noise, so rho must not depend on the thread schedule
	std::vector<vec> part(nQM);
	int merged = 0;
	long long done = 0, nfar = 0, nref = 0;
	const double far_lnk = use_far ? -std::log(far_eps) : 0.0;
	{ //the bar ends its line when it goes out of scope, before the summary below
		ProgressBar pb(nQM, 60, "#", " ", "Calculating Eri3c Matrix");
#pragma omp parallel reduction(+:done, nfar, nref)
		{
			vec cache(ncache), dblk(static_cast<size_t>(dmax_orb) * dmax_orb);
			vec buf3(static_cast<size_t>(dmax_orb) * dmax_orb * dmax_aux), buf2(static_cast<size_t>(dmax_orb) * dmax_orb * dmax_orb * dmax_orb);
			//far field per thread: primitive pairs (p, centre) of the current a, b; the contracted reference integrals,
			//valid while refstamp matches the pair's stamp
			vec pp(static_cast<size_t>(4) * maxprim * maxprim), refv(static_cast<size_t>(nbx - nbas) * dmax_aux);
			std::vector<long long> refstamp(nbx - nbas, -1);
			long long stamp = 0;
			auto contract = [&](int k, int dab) {
				const double *col = buf3.data() + static_cast<size_t>(dab) * k;
				double sum = 0.0;
				for (int ij = 0; ij < dab; ij++) sum += col[ij] * dblk[ij];
				return sum;
			};
#pragma omp for schedule(dynamic)
			for (int a = 0; a < nQM; a++) {
				vec loc(naux, 0.0);
				const int a0 = aoloc[a], da = aoloc[a + 1] - a0;
				for (int b = 0; b <= a; b++) {
					const int b0 = aoloc[b], db = aoloc[b + 1] - b0, dab = da * db;
					const double w = a == b ? 1.0 : 2.0;
					//same element order as the integral block: i (shell a) fastest, then j (shell b)
					double dm_max = 0.0;
					for (int j = 0; j < db; j++)
						for (int i = 0; i < da; i++) {
							dblk[i + da * j] = w * dm(b0 + j, a0 + i);
							dm_max = std::max(dm_max, std::abs(dblk[i + da * j]));
						}
					if (dm_max == 0.0) continue;
					double qab = 1.0;
					if (bound) {
						int shls[4] = { a, b, a, b };
						double m = 0.0;
						if (S::pair(buf2.data(), nullptr, shls, atm.data(), nat, bas.data(), nQM, env.data(), opt2, cache.data()))
							for (int j = 0; j < db; j++)
								for (int i = 0; i < da; i++) m = std::max(m, std::abs(buf2[i + da * j + dab * (i + da * j)]));
						qab = std::sqrt(m);
						if (dm_max * qab * qaux_max < thr) continue;
					}
					int npp = 0, atom_c = -1;
					double qstar = 0.0;
					stamp++;
					if (use_far) { //primitive pairs with exp(-ab/p |AB|^2) < eps carry no density
						const int *sa = &bas[a * BAS_SLOTS], *sb = &bas[b * BAS_SLOTS];
						const double *ea = env.data() + sa[PTR_EXP], *eb = env.data() + sb[PTR_EXP];
						const double *A = env.data() + atm[PTR_COORD + ATM_SLOTS * sa[ATOM_OF]], *B = env.data() + atm[PTR_COORD + ATM_SLOTS * sb[ATOM_OF]];
						const double AB2 = (A[0] - B[0]) * (A[0] - B[0]) + (A[1] - B[1]) * (A[1] - B[1]) + (A[2] - B[2]) * (A[2] - B[2]);
						for (int i = 0; i < sa[NPRIM_OF]; i++)
							for (int j = 0; j < sb[NPRIM_OF]; j++) {
								const double p = ea[i] + eb[j];
								if (ea[i] * eb[j] / p * AB2 > far_lnk) continue;
								double *q = &pp[4 * npp++];
								q[0] = p;
								for (int d = 0; d < 3; d++) q[d + 1] = (ea[i] * A[d] + eb[j] * B[d]) / p;
							}
					}
					for (int P = 0; P < nAux; P++) {
						if (bound && dm_max * qab * qaux[P] < thr) continue;
						done++;
						const int p0 = aoloc[nQM + P] - aux0, dp = aoloc[nQM + P + 1] - aoloc[nQM + P];
						if (ref_of[P] >= 0) {
							const int c = bas[(nQM + P) * BAS_SLOTS + ATOM_OF];
							if (c != atom_c) { //smallest q that is far from every primitive pair: p q R^2 / (p + q) >= X^2
								atom_c = c;
								const double *C = env.data() + atm[PTR_COORD + ATM_SLOTS * c];
								qstar = 0.0;
								for (int k = 0; k < npp; k++) {
									const double *q = &pp[4 * k];
									const double R2 = (q[1] - C[0]) * (q[1] - C[0]) + (q[2] - C[1]) * (q[2] - C[1]) + (q[3] - C[2]) * (q[3] - C[2]);
									const double d = q[0] * R2 - X2;
									if (d <= 0.0) { qstar = HUGE_VAL; break; }
									qstar = std::max(qstar, X2 * q[0] / d);
								}
							}
							if (qmin[P] >= qstar) {
								const int r = ref_of[P] - nbas;
								double *v = refv.data() + static_cast<size_t>(r) * dmax_aux;
								if (refstamp[r] != stamp) {
									refstamp[r] = stamp;
									nref++;
									int shls[3] = { a, b, ref_of[P] };
									const bool nz = S::three(buf3.data(), nullptr, shls, atm.data(), nat, bas.data(), nbx, env.data(), opty, cache.data());
									for (int k = 0; k < dp; k++) v[k] = nz ? contract(k, dab) : 0.0;
								}
								for (int k = 0; k < dp; k++) loc[p0 + k] += ratio[P] * v[k];
								nfar++;
								continue;
							}
						}
						int shls[3] = { a, b, nQM + P };
						if (!S::three(buf3.data(), nullptr, shls, atm.data(), nat, bas.data(), nbx, env.data(), opty, cache.data())) continue;
						for (int k = 0; k < dp; k++) loc[p0 + k] += contract(k, dab);
					}
				}
				pb.update();
#pragma omp critical
				{
					part[a] = std::move(loc);
					for (; merged < nQM && !part[merged].empty(); merged++) {
						for (int p = 0; p < naux; p++) rho[p] += part[merged][p];
						vec().swap(part[merged]);
					}
				}
			}
		}
	}
	libcint::CINTdel_optimizer(&opty);
	if (opt2) libcint::CINTdel_optimizer(&opt2);
	const long long total = static_cast<long long>(nQM) * (nQM + 1) / 2 * nAux;
	if (bound) std::cout << "Charge-weighted Schwarz screening (" << thr << " e) kept " << done << " of " << total << " shell triplets" << std::endl;
	if (use_far) std::cout << "Far field (erfc < " << far_eps << ") took " << nfar << " of " << done << " triplets from " << nref << " reference integrals" << std::endl;
}
template void computeRho<Coulomb3C_SPH>(
	const Int_Params &normal_basis,
	const Int_Params &aux_basis,
	const dMatrix2 &dm,
	vec &rho,
	const std::optional<ivec> asym_atm_list,
	const vec *sensitivity);
template void computeRho<Coulomb3C_CRT>(
	const Int_Params &normal_basis,
	const Int_Params &aux_basis,
	const dMatrix2 &dm,
	vec &rho,
	const std::optional<ivec> asym_atm_list,
	const vec *sensitivity);
template void computeRho<Overlap3C_SPH>(
	const Int_Params &normal_basis,
	const Int_Params &aux_basis,
	const dMatrix2 &dm,
	vec &rho,
	const std::optional<ivec> asym_atm_list,
	const vec *sensitivity);


template <typename Kernel>
void compute3C(Int_Params &param1,
	Int_Params &param2,
	vec &eri3c) {
	int nQM = param1.get_nbas();
	int nAux = param2.get_nbas();
	Int_Params combined(param1, param2);
	ivec bas = combined.get_bas();
	ivec atm = combined.get_atm();
	vec env = combined.get_env();
	int nat = combined.get_natoms();
	int nbas = combined.get_nbas();
	//ivec aoloc = make_loc(bas, nbas);
	ivec aoloc = Kernel::gen_loc(bas, nbas);
	unsigned long long int naoi = aoloc[nQM] - aoloc[0];
	unsigned long long int naoj = aoloc[nQM] - aoloc[0];
	unsigned long long int naok = aoloc[nQM + nAux] - aoloc[nQM];
	eri3c.resize(naoi * naoj * naok, 0.0);
	libcint::CINTOpt* opty = nullptr;
	Kernel::optimizer(opty, atm.data(), nat, bas.data(), nbas, env.data());
	ivec shl_slice = { 0, nQM, 0, nQM, nQM, nQM + nAux };
	Kernel::drv(eri3c.data(), 1, shl_slice.data(), aoloc.data(), opty, atm.data(), nat, bas.data(), nbas, env.data());

}
template void compute3C<Coulomb3C_SPH>(Int_Params &param1,
	Int_Params &param2,
	vec &eri3c);
template void compute3C<Overlap3C_SPH>(Int_Params &param1,
	Int_Params &param2,
	vec &eri3c);



//Cartesian to real spherical transformation matrix
//
//Kwargs :
//normalized:
//How the Cartesian GTOs are normalized.  'sp' means the s and p
//functions are normalized(this is the convention used by libcint
//    library).
dMatrix2 cart2sph(const int l, const bool normalized) {
	int n_cart = (l + 1) * (l + 2) / 2;

	dMatrix2 c_tensor(n_cart, n_cart);
	for (int i = 0; i < n_cart; i++) {
		c_tensor(i, i) = 1.0;
	}

	if (l == 0 || l == 1) { //For s and p functions, the transformation is trivial
		if (l == 1) { //libcint is built with PYPZPX: spherical p comes as py, pz, px
			std::fill(c_tensor.container().begin(), c_tensor.container().end(), 0.0);
			c_tensor(1, 0) = c_tensor(2, 1) = c_tensor(0, 2) = 1.0;
		}
		if (normalized) {
			return c_tensor;
		}
		else {
			double norm_factor = (l == 0) ? 0.282094791773878143 : 0.488602511902919921;
			for (auto &val : c_tensor.container()) {
				val *= norm_factor;
			}
			return c_tensor;
		}
	}

	err_checkf(l <= 15, "cart2sph_matrix: l must be <= 15", std::cout);

	int n_sph = 2 * l + 1;
	vec c_sph(static_cast<size_t>(n_sph) * n_cart, 0.0);

	libcint::CINTc2s_ket_sph(c_sph.data(), n_cart, c_tensor.data(), l);
	//Transform back to row-major order
	dMatrix2 c_sph_RM(n_cart, n_sph);
	for (int i = 0; i < n_cart; i++) {
		for (int j = 0; j < n_sph; j++) {
			c_sph_RM(i, j) = c_sph[j * n_cart + i];
		}
	}
	return c_sph_RM;
}


//Returns a n_cart*n_sph matrix used in transforming a cartesian DM to a spherical DM
dMatrix2 get_cart2sph_matrix(const WFN &cart_wfn, const bool normalized) {
	//First collect the complete number of spherical and cartesian functions used in the wavefunction
	int max_l = 0;
	int n_cart = 0, n_sph = 0;
	for (const atom &a : cart_wfn.get_atoms()) {
		int prim = 0;
		for (int shell = 0; shell < a.get_shellcount_size(); shell++) {
			const int type = a.get_basis_set_type(prim) - 1;

			n_cart += ((type + 1) * (type + 2)) / 2;
			n_sph += 2 * type + 1;
			if (type > max_l) {
				max_l = type;
			}

			prim += a.get_shellcount(shell);
		}
	}

	std::vector<dMatrix2> conversion_matrices(max_l + 1);
	for (int l = 0; l <= max_l; l++) {
		conversion_matrices[l] = cart2sph(l, normalized);
	}

	std::cout << "Number of cartesian functions: " << n_cart << ", number of spherical functions: " << n_sph << std::endl;
	dMatrix2 c_dm(n_cart, n_sph);
	int cart_idx = 0, sph_idx = 0;
	for (const atom &a : cart_wfn.get_atoms()) {
		int prim = 0;
		for (int shell = 0; shell < a.get_shellcount_size(); shell++) {
			const int type = a.get_basis_set_type(prim) - 1;
			const dMatrix2 &c_sph = conversion_matrices[type];
			const int n_cart_shell = ((type + 1) * (type + 2)) / 2;
			const int n_sph_shell = 2 * type + 1;
			//Fill in the appropriate block in the c_dm matrix
			for (int i = 0; i < n_cart_shell; i++) {
				for (int j = 0; j < n_sph_shell; j++) {
					c_dm(cart_idx + i, sph_idx + j) = c_sph(i, j);
				}
			}
			cart_idx += n_cart_shell;
			sph_idx += n_sph_shell;
			prim += a.get_shellcount(shell);
		}
	}
	////print all of c_dm
	//for (int i = 0; i < c_dm.extent(0); i++) {
	//    for (int j = 0; j < c_dm.extent(1); j++) {
	//        std::cout << c_dm(i, j) << " ";
	//    }
	//    std::cout << std::endl;
	//}

	return c_dm;
}


vec eval_GTO_sph(Int_Params& params, vec2& grid, ivec& shl_slice) {
	ivec bas = params.get_bas();
	ivec atm = params.get_atm();
	vec env = params.get_env();


	//grid = numpy.asarray(grid, dtype = numpy.double, order = 'F'); grid holds the three coordinate rows
	const int ngrid = (int)grid[0].size();
	vec fortran_grid((size_t)ngrid * 3);
	for (int i = 0; i < ngrid; i++) {
		for (int j = 0; j < 3; j++) {
			fortran_grid[(size_t)i + (size_t)j * ngrid] = grid[j][i];
		}
	}

	int nbas = params.get_nbas();
	int nat = params.get_natoms();

	if (shl_slice.size() == 0) {
		shl_slice = { 0, nbas};
	}
	//ivec aoloc = Kernel::gen_loc(bas, nbas);
	ivec aoloc = make_loc<COORDINATE_TYPE::SPH>(bas, nbas);
	int nao = aoloc[shl_slice[1]] - aoloc[shl_slice[0]];

	//non0tab = numpy.ones(((ngrids+BLKSIZE-1)//BLKSIZE,nbas),dtype = numpy.uint8)
	std::vector<uint8_t> non0table(static_cast<size_t>((ngrid + 56 - 1) / 56) * nbas, 1);

	// Compute integrals
	vec res((size_t)nao * (size_t)ngrid, 0.0);
	GTOval_sph(ngrid, shl_slice.data(), aoloc.data(), res.data(), fortran_grid.data(), non0table.data(), atm.data(), nat, bas.data(), nbas, env.data());

	vec ret((size_t)nao * (size_t)ngrid, 0.0);
	// res is in fortran order, write the result in regular ordering
	for (int i = 0; i < nao; i++)
	{
		for (int j = 0; j < ngrid; j++)
		{
			ret[(size_t)j * (size_t)nao + i] = res[(size_t)i * (size_t)ngrid + j];
		}
	}
	return ret;
}