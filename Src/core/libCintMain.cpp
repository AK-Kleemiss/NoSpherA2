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

//d^t/dX^t d^u/dY^u d^v/dZ^v 1/|D| for t+u+v <= L into R[(t S + u) S + v]: McMurchie-Davidson with the Boys function of an
//infinite exponent, R^(n)_000 = (1/D d/dD)^n 1/D = (-1)^n (2n-1)!!/D^(2n+1), R^(n)_(t+1)uv = t R^(n+1)_(t-1)uv + X R^(n+1)_tuv.
//The levels n alternate between R and the scratch T, so n = 0 lands in R
static void inverse_r_derivatives(const int L, const int S, const double *D, double *R, double *T) {
	const double r2 = D[0] * D[0] + D[1] * D[1] + D[2] * D[2];
	double g[64];
	g[0] = 1.0 / std::sqrt(r2);
	for (int n = 0; n < L; n++) g[n + 1] = -(2 * n + 1) * g[n] / r2;
	const int S2 = S * S;
	for (int n = L; n >= 0; n--) {
		double *cur = (n & 1) ? T : R;
		const double *prev = (n & 1) ? R : T;
		cur[0] = g[n];
		for (int s = 1; s <= L - n; s++)
			for (int t = s; t >= 0; t--)
				for (int u = s - t; u >= 0; u--) {
					const int v = s - t - u, i = (t * S + u) * S + v;
					if (t > 0) cur[i] = D[0] * prev[i - S2] + (t > 1 ? (t - 1) * prev[i - 2 * S2] : 0.0);
					else if (u > 0) cur[i] = D[1] * prev[i - S] + (u > 1 ? (u - 1) * prev[i - 2 * S] : 0.0);
					else cur[i] = D[2] * prev[i - 1] + (v > 1 ? (v - 1) * prev[i - 2] : 0.0);
				}
	}
}

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
	//Far field: once every orbital primitive pair (exponent p, centre P_k) overlaps an aux shell on atom C by less than
	//erfc(x) < eps, x^2 = pq/(p+q) |P_k - C|^2 with q the shell's most diffuse exponent, the shell R(r) S_lm(r) (S_lm = r^l Y_lm)
	//acts as the point multipole 4pi/(2l+1) M Y_lm/r_C^(l+1), M = sum_i c_i int r^(2l+2) e^(-q_i r^2) dr, and a Gaussian
	//e^(-p|r-P|^2) sees that harmonic potential only at its centre, times (pi/p)^(3/2). With Hobson's
	//Y_lm/r^(l+1) = (-1)^l/(2l-1)!! S_lm(grad) 1/r and the pair density in Hermite Gaussians (McMurchie-Davidson,
	//sum_k sum_tuv h^k_tuv d^tuv/dP^tuv e^(-p_k|r-P_k|^2)), (ab|P_m) = M_P U_lm with
	//U_lm = 4pi/(2l+1) (-1)^l/(2l-1)!! sum_k (pi/p_k)^(3/2) sum_tuv h^k_tuv S_lm(grad) d^tuv/dD^tuv 1/|D|, D = P_k - C,
	//one sum per orbital pair and aux atom for every l on it.
	double far_eps = constants::ri_far_erfc;
	if (const char *e = tuning("NOS_RI_FAR")) far_eps = std::atof(e);
	bool use_far = S::far_field && far_eps > 0.0;
	double X2 = 0.0;
	vec mom(nAux, 0.0), qmin(nAux, 0.0);
	ivec far_l(nAux, -1), lmax_at(nat, -1); //far_l: -1 = never far (general contraction)
	int maxprim = 1, lorb = 0, laux = 0;
	for (int s = 0; s < nQM; s++) {
		const int *sh = &bas[s * BAS_SLOTS];
		use_far = use_far && sh[NCTR_OF] == 1; //the Hermite expansion takes one contraction per orbital shell
		maxprim = std::max(maxprim, sh[NPRIM_OF]);
		lorb = std::max(lorb, sh[ANG_OF]);
	}
	if (use_far) {
		double lo = 0.0, hi = 30.0; //erfc(X) = eps
		for (int it = 0; it < 100; it++) {
			const double mid = 0.5 * (lo + hi);
			(std::erfc(mid) > far_eps ? lo : hi) = mid;
		}
		X2 = hi * hi;
		for (int P = 0; P < nAux; P++) {
			const int *sh = &bas[(nQM + P) * BAS_SLOTS];
			if (sh[NCTR_OF] != 1) continue;
			const double *ex = env.data() + sh[PTR_EXP], *co = env.data() + sh[PTR_COEFF];
			const int l = sh[ANG_OF];
			qmin[P] = *std::min_element(ex, ex + sh[NPRIM_OF]);
			for (int i = 0; i < sh[NPRIM_OF]; i++) mom[P] += co[i] * std::tgamma(l + 1.5) / (2.0 * std::pow(ex[i], l + 1.5));
			far_l[P] = l;
			lmax_at[sh[ATOM_OF]] = std::max(lmax_at[sh[ATOM_OF]], l);
			laux = std::max(laux, l);
		}
	}
	//Hermite/Cartesian triples (t,u,v) by total degree, libcint's Cartesian order within a degree, so degree l starts at
	//nh(l - 1); hoff is a triple's offset in a dense hs^3 array, where offsets add like the triples, and hidx maps back
	const int hs = 2 * lorb + laux + 1;
	auto nh = [](const int n) { return (n + 1) * (n + 2) * (n + 3) / 6; };
	ivec hoff, htuv, hidx;
	std::vector<vec> c2s(std::max(lorb, laux) + 1); //S_lm = sum_x c2s[l][x (2l+1) + m] x^lx y^ly z^lz, s/p factors included
	vec pref(laux + 1);
	if (use_far) {
		hidx.assign(static_cast<size_t>(hs) * hs * hs, -1);
		for (int s = 0; s < hs; s++)
			for (int t = s; t >= 0; t--)
				for (int u = s - t; u >= 0; u--) {
					const int o = (t * hs + u) * hs + s - t - u;
					hidx[o] = static_cast<int>(hoff.size());
					hoff.push_back(o);
					htuv.insert(htuv.end(), { t, u, s - t - u });
				}
		for (int l = 0; l < static_cast<int>(c2s.size()); l++) {
			dMatrix2 m = cart2sph(l, false);
			const int ns = 2 * l + 1;
			c2s[l].resize(static_cast<size_t>(nh(l) - nh(l - 1)) * ns);
			for (int x = 0; x < nh(l) - nh(l - 1); x++)
				for (int j = 0; j < ns; j++) c2s[l][x * ns + j] = m(x, j);
		}
		double dfact = 1.0; //(2l-1)!!
		for (int l = 0; l <= laux; l++) {
			pref[l] = constants::FOUR_PI / (2 * l + 1) * (l % 2 ? -1.0 : 1.0) / dfact;
			dfact *= 2 * l + 1;
		}
	}
	//libcint mallocs scratch per call unless handed one; the size query is the call without output.
	//{s,s,s,s} over every shell is the bound pyscf's GTOmax_cache_size uses
	CACHE_SIZE_T ncache = 0;
	for (int s = 0; s < nbas; s++) {
		int shls[4] = { s, s, s, s };
		ncache = std::max(ncache, S::three(nullptr, nullptr, shls, atm.data(), nat, bas.data(), nbas, env.data(), nullptr, nullptr));
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
	Kernel::optimizer(opty, atm.data(), nat, bas.data(), nbas, env.data());
	//(ab|ab) only touches orbital shells, the first nQM, so its optimizer skips the aux pairs
	if (bound) libcint::int2e_optimizer(&opt2, atm.data(), nat, bas.data(), nQM, env.data());

	//Every shell a sums into its own vector, added to rho in a fixed order once all before it are done: the
	//ill-conditioned metric turns summation-order noise into tsc noise, so rho must not depend on the thread schedule
	std::vector<vec> part(nQM);
	int merged = 0;
	long long done = 0, nfar = 0, nsum = 0;
	const double far_lnk = use_far ? -std::log(far_eps) : 0.0;
	{ //the bar ends its line when it goes out of scope, before the summary below
		ProgressBar pb(nQM, 60, "#", " ", "Calculating Eri3c Matrix");
#pragma omp parallel reduction(+:done, nfar, nsum)
		{
			vec cache(ncache), dblk(static_cast<size_t>(dmax_orb) * dmax_orb);
			vec buf3(static_cast<size_t>(dmax_orb) * dmax_orb * dmax_aux), buf2(static_cast<size_t>(dmax_orb) * dmax_orb * dmax_orb * dmax_orb);
			//far field per thread: primitive pairs {p, P_k, c_a c_b e^(-ab/p AB^2) (pi/p)^(3/2)} of the current a, b, their
			//Hermite coefficients h^k (one summed set when a and b share an atom: every P_k is that atom), the Cartesian
			//density block, E^ij_t per axis, d^tuv 1/|D| and its scratch, and G_xyz / U_lm of the last aux atom
			const int ncart = nh(lorb) - nh(lorb - 1), ne = (lorb + 1) * (lorb + 1) * (2 * lorb + 1);
			vec pp(static_cast<size_t>(5) * maxprim * maxprim), hk(static_cast<size_t>(maxprim) * maxprim * nh(2 * lorb));
			vec dc(static_cast<size_t>(ncart) * ncart), tc(static_cast<size_t>(ncart) * dmax_orb), E(3 * static_cast<size_t>(ne));
			vec Rb(static_cast<size_t>(hs) * hs * hs), Tb(Rb.size()), G(nh(laux)), U(static_cast<size_t>(laux + 1) * (laux + 1));
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
					const int *sa = &bas[a * BAS_SLOTS], *sb = &bas[b * BAS_SLOTS];
					const double *A = env.data() + atm[PTR_COORD + ATM_SLOTS * sa[ATOM_OF]], *B = env.data() + atm[PTR_COORD + ATM_SLOTS * sb[ATOM_OF]];
					const int la = sa[ANG_OF], lb = sb[ANG_OF], lab = la + lb, nhab = nh(lab);
					const bool same = sa[ATOM_OF] == sb[ATOM_OF];
					int npp = 0, nk = -1, atom_c = -1, atom_u = -1; //nk < 0: h^k not built yet for this pair
					double qstar = 0.0;
					const double *C = nullptr;
					if (use_far) { //primitive pairs with exp(-ab/p |AB|^2) < eps carry no density
						const double *ea = env.data() + sa[PTR_EXP], *eb = env.data() + sb[PTR_EXP];
						const double *ca = env.data() + sa[PTR_COEFF], *cb = env.data() + sb[PTR_COEFF];
						const double AB2 = (A[0] - B[0]) * (A[0] - B[0]) + (A[1] - B[1]) * (A[1] - B[1]) + (A[2] - B[2]) * (A[2] - B[2]);
						for (int i = 0; i < sa[NPRIM_OF]; i++)
							for (int j = 0; j < sb[NPRIM_OF]; j++) {
								const double p = ea[i] + eb[j], mu = ea[i] * eb[j] / p * AB2;
								if (mu > far_lnk) continue;
								double *q = &pp[5 * npp++];
								q[0] = p;
								for (int d = 0; d < 3; d++) q[d + 1] = (ea[i] * A[d] + eb[j] * B[d]) / p;
								q[4] = ca[i] * cb[j] * std::exp(-mu) * constants::PI3_2 / (p * std::sqrt(p));
							}
					}
					for (int P = 0; P < nAux; P++) {
						if (bound && dm_max * qab * qaux[P] < thr) continue;
						done++;
						const int p0 = aoloc[nQM + P] - aux0, dp = aoloc[nQM + P + 1] - aoloc[nQM + P];
						if (far_l[P] >= 0) {
							const int c = bas[(nQM + P) * BAS_SLOTS + ATOM_OF];
							if (c != atom_c) { //smallest q that is far from every primitive pair: p q R^2 / (p + q) >= X^2
								atom_c = c;
								C = env.data() + atm[PTR_COORD + ATM_SLOTS * c];
								qstar = 0.0;
								for (int k = 0; k < npp; k++) {
									const double *q = &pp[5 * k];
									const double R2 = (q[1] - C[0]) * (q[1] - C[0]) + (q[2] - C[1]) * (q[2] - C[1]) + (q[3] - C[2]) * (q[3] - C[2]);
									const double d = q[0] * R2 - X2;
									if (d <= 0.0) { qstar = HUGE_VAL; break; }
									qstar = std::max(qstar, X2 * q[0] / d);
								}
							}
							if (qmin[P] >= qstar) {
								if (nk < 0) { //h^k_tuv = pre_k sum_xy Dcart_xy E^(ax bx)_t E^(ay by)_u E^(az bz)_v, Dcart = C_a dblk C_b^T
									const int na = nh(la) - nh(la - 1), nb = nh(lb) - nh(lb - 1), nt = lab + 1, nij = (lb + 1) * nt;
									for (int x = 0; x < na; x++)
										for (int j = 0; j < db; j++) {
											double s = 0.0;
											for (int i = 0; i < da; i++) s += c2s[la][x * da + i] * dblk[i + da * j];
											tc[x * db + j] = s;
										}
									for (int x = 0; x < na; x++)
										for (int y = 0; y < nb; y++) {
											double s = 0.0;
											for (int j = 0; j < db; j++) s += tc[x * db + j] * c2s[lb][y * db + j];
											dc[x * nb + y] = s;
										}
									nk = same ? 1 : npp;
									std::fill(hk.begin(), hk.begin() + static_cast<size_t>(nk) * nhab, 0.0);
									for (int k = 0; k < npp; k++) {
										const double *q = &pp[5 * k];
										for (int d = 0; d < 3; d++) { //E^(i+1)j_t = E^ij_(t-1)/2p + X_PA E^ij_t + (t+1) E^ij_(t+1), likewise j with X_PB
											double *e = &E[d * static_cast<size_t>(ne)];
											std::fill(e, e + (la + 1) * nij, 0.0);
											e[0] = 1.0;
											const double xa = q[d + 1] - A[d], xb = q[d + 1] - B[d], h2p = 0.5 / q[0];
											auto step = [&](const double *src, double *dst, const double x, const int tmax) {
												for (int t = 0; t <= tmax; t++)
													dst[t] = (t > 0 ? h2p * src[t - 1] : 0.0) + x * src[t] + (t < lab ? (t + 1) * src[t + 1] : 0.0);
											};
											for (int i = 0; i < la; i++) step(e + i * nij, e + (i + 1) * nij, xa, i + 1);
											for (int j = 0; j < lb; j++)
												for (int i = 0; i <= la; i++) step(e + i * nij + j * nt, e + i * nij + (j + 1) * nt, xb, i + j + 1);
										}
										double *h = &hk[static_cast<size_t>(same ? 0 : k) * nhab];
										for (int x = 0; x < na; x++)
											for (int y = 0; y < nb; y++) {
												const double w0 = q[4] * dc[x * nb + y];
												if (w0 == 0.0) continue;
												const int *ia = &htuv[3 * (nh(la - 1) + x)], *ib = &htuv[3 * (nh(lb - 1) + y)];
												const double *ex = &E[(ia[0] * (lb + 1) + ib[0]) * nt], *ey = &E[ne + (ia[1] * (lb + 1) + ib[1]) * nt],
													*ez = &E[2 * static_cast<size_t>(ne) + (ia[2] * (lb + 1) + ib[2]) * nt];
												for (int t = 0; t <= ia[0] + ib[0]; t++)
													for (int u = 0; u <= ia[1] + ib[1]; u++) {
														const double wtu = w0 * ex[t] * ey[u];
														for (int v = 0; v <= ia[2] + ib[2]; v++) h[hidx[(t * hs + u) * hs + v]] += wtu * ez[v];
													}
											}
									}
								}
								if (c != atom_u) { //G_xyz = sum_k sum_tuv h^k_tuv d^(tuv+xyz) 1/|P_k - C|, then U_lm for every l on C
									atom_u = c;
									nsum++;
									const int lc = lmax_at[c], nhc = nh(lc);
									std::fill(G.begin(), G.begin() + nhc, 0.0);
									for (int k = 0; k < nk; k++) {
										const double *Pk = same ? A : &pp[5 * k + 1];
										const double D[3] = { Pk[0] - C[0], Pk[1] - C[1], Pk[2] - C[2] };
										inverse_r_derivatives(lab + lc, hs, D, Rb.data(), Tb.data());
										const double *h = &hk[static_cast<size_t>(k) * nhab];
										for (int g = 0; g < nhc; g++) {
											const double *r = Rb.data() + hoff[g];
											double s = 0.0;
											for (int i = 0; i < nhab; i++) s += h[i] * r[hoff[i]];
											G[g] += s;
										}
									}
									for (int l = 0; l <= lc; l++)
										for (int m = 0; m < 2 * l + 1; m++) {
											double s = 0.0;
											for (int x = 0; x < nh(l) - nh(l - 1); x++) s += c2s[l][x * (2 * l + 1) + m] * G[nh(l - 1) + x];
											U[l * l + m] = pref[l] * s;
										}
								}
								const double *u = &U[far_l[P] * far_l[P]];
								for (int k = 0; k < dp; k++) loc[p0 + k] += mom[P] * u[k];
								nfar++;
								continue;
							}
						}
						int shls[3] = { a, b, nQM + P };
						if (!S::three(buf3.data(), nullptr, shls, atm.data(), nat, bas.data(), nbas, env.data(), opty, cache.data())) continue;
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
	if (use_far) std::cout << "Far field (erfc < " << far_eps << ") took " << nfar << " of " << done << " triplets as point multipoles from " << nsum << " pair-atom sums" << std::endl;
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