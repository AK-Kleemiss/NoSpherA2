#include "pch.h"
#include "SALTED_equicomb.h"

#if defined(__APPLE__)
// On macOS we�re using Accelerate for BLAS/LAPACK
#include <Accelerate/Accelerate.h>
#elif !defined(NSA2_OPENBLAS)
// Linux/Windows with oneMKL
#include <mkl.h>
#endif

#include "constants.h"

static bool g_equicomb_use_gpu = false;
void equicomb_set_gpu(bool on) { g_equicomb_use_gpu = on; }
bool equicomb_gpu_enabled() { return g_equicomb_use_gpu; }
#ifdef NOSPHERA2_USE_GPU
#include "salted_gpu.h"
#endif

// Start of the l block in one atom's matrices: sum_{k<l} (2k+1)^2
static size_t dm_offset(const int l) { return static_cast<size_t>(l * (2 * l - 1) * (2 * l + 1) / 3); }

// One atom's A[m,m'] = sum_n x_m conj(x_m') for l = 0..lmax into m ([are, aim][(2l+1)^2 blocks],
// dsz = dm_offset(lmax + 1) apart), real and imaginary parts apart: plain double loops
// vectorise, std::complex products without -ffast-math do not. The sums over n run over
// xr/xi[m * nch + n] ((2 lmax + 1) * nch each); for a real density A[-m',-m] =
// (-1)^(m+m') A[m,m'] and A is Hermitian, so only m' <= m, m + m' <= 0 are summed.
static void atom_density(const SALTEDDescriptors &d, const int nch, const int lmax, const int iat, const size_t dsz, double *m,
	double *__restrict xr, double *__restrict xi)
{
	for (int l = 0; l <= lmax; ++l)
	{
		const int nm = 2 * l + 1;
		const size_t o = dm_offset(l);
		double *__restrict are = m + o, *__restrict aim = m + dsz + o;
		for (int n = 0; n < nch; ++n)
		{
			const cdouble *x = d.block(iat, n, l);
			for (int a = 0; a < nm; ++a) { xr[a * nch + n] = x[a].real(); xi[a * nch + n] = x[a].imag(); }
		}
		for (int i = 0; i <= l; ++i)
			for (int j = i; j < nm - i; ++j)
			{
				const double *rj = xr + j * nch, *ij = xi + j * nch, *ri = xr + i * nch, *ii = xi + i * nch;
				double sr = 0.0, si = 0.0;
				for (int n = 0; n < nch; ++n)
				{
					sr += rj[n] * ri[n] + ij[n] * ii[n];
					si += ij[n] * ri[n] - rj[n] * ii[n];
				}
				const double s = (i + j) % 2 ? -1.0 : 1.0;
				const int jm = nm - 1 - j, im = nm - 1 - i;
				are[j * nm + i] = sr; aim[j * nm + i] = si;
				are[i * nm + j] = sr; aim[i * nm + j] = -si;
				are[im * nm + jm] = s * sr; aim[im * nm + jm] = s * si;
				are[jm * nm + im] = s * sr; aim[jm * nm + im] = -s * si;
			}
	}
}

// Which (im1, im2) pairs contribute is set by |im1 - mu| <= l2, which depends only on
// (il, imu, lam), never on n1 or n2. That test selects a CONTIGUOUS run of im1, and
// im2 = im1 - mu + l2 advances in lockstep with it, so a first index, a length and an
// offset into w3j describe a group completely. w3j is consumed in exactly (il, imu, im1)
// order, so a group's weights are a contiguous slice of it and need no copy.
struct w3j_run { int im1_begin, im2_begin, count, w_off; };

// Everything of one lambda the walk needs that does not depend on the atoms
struct lambda_plan
{
	int lam = 0, l21 = 1, llmax = 0, lmax1 = 0, lmax2 = 0;
	const vec *w3j = nullptr;
	const ivec2 *llvec = nullptr;
	std::vector<w3j_run> runs;
	size_t total_terms = 0;
	vec cpm_re, cpm_im;
};

static lambda_plan make_plan(const vec &w3j, const ivec2 &llvec, const int lam, const cvec2 &c2r)
{
	lambda_plan q;
	q.lam = lam; q.l21 = 2 * lam + 1; q.llmax = static_cast<int>(llvec[0].size());
	q.w3j = &w3j; q.llvec = &llvec;
	const int l21 = q.l21;
	q.runs.assign(static_cast<size_t>(q.llmax) * l21, w3j_run{0, 0, 0, 0});
	int w_idx = 0;
	for (int til = 0; til < q.llmax; ++til)
	{
		const int tl1 = llvec[0][til], tl2 = llvec[1][til];
		err_checkf(tl1 >= 0 && tl2 >= 0, "equicomb: negative angular momentum in llvec", std::cout);
		err_checkf((tl1 + tl2 + lam) % 2 == 0, "equicomb: l1 + l2 + lambda odd in llvec", std::cout);
		for (int timu = 0; timu < l21; ++timu)
		{
			const int tmu = timu - lam + tl1;
			const int lo = std::max(0, tmu - tl2);
			const int hi = std::min(2 * tl1, tmu + tl2);
			const int cnt = (hi >= lo) ? (hi - lo + 1) : 0;
			q.runs[static_cast<size_t>(til) * l21 + timu] = { lo, lo - tmu + tl2, cnt, w_idx };
			w_idx += cnt;
		}
	}
	q.total_terms = static_cast<size_t>(w_idx);
	// w3j is consumed once per (n1,n2) shell pair, one entry per surviving m2, which is what the runs add up to
	err_checkf(w3j.size() >= q.total_terms, "equicomb: w3j holds " + std::to_string(w3j.size()) +
		" entries, the shell loop consumes " + std::to_string(q.total_terms), std::cout);
	q.lmax1 = *std::max_element(llvec[0].begin(), llvec[0].end());
	q.lmax2 = *std::max_element(llvec[1].begin(), llvec[1].end());

	// The complex-to-real matrix is a mirror-pair transform: row i couples only column i and
	// column l21-1-i. The structure is read off c2r rather than assumed; anything else is an error.
	for (int i2 = 0; i2 < l21; ++i2)
		for (int j2 = 0; j2 < l21; ++j2)
			err_checkf((c2r[i2][j2].real() == 0.0 && c2r[i2][j2].imag() == 0.0) || j2 == i2 || j2 == l21 - 1 - i2,
				"equicomb: complex-to-real row " + std::to_string(i2) + " is not a mirror pair", std::cout);
	// pc[-mu] = (-1)^mu conj(pc[mu]), so row imu of c2r gives Re(cpm[imu] pc[|mu|]) with
	// cpm = c2r[imu][lam + mu] + (-1)^mu conj(c2r[imu][lam - mu])
	q.cpm_re.assign(l21, 0.0);
	q.cpm_im.assign(l21, 0.0);
	q.cpm_re[lam] = c2r[lam][lam].real();
	for (int mu = 1; mu <= lam; ++mu)
		for (const int imu : { lam - mu, lam + mu })
		{
			const cdouble c = c2r[imu][lam + mu] + (mu % 2 ? -1.0 : 1.0) * std::conj(c2r[imu][lam - mu]);
			q.cpm_re[imu] = c.real();
			q.cpm_im[imu] = c.imag();
		}
	return q;
}

// The l21 real features of channels (n1, n2) and (l1, l2) = llvec[.][il] of one atom into
// out[imu * stride]: Re(c2r pc), pc[mu] = sum_m1 w v1[n1,l1,m1] u[n2,l2,m1-mu] with u = conj(v1)
// or v2 as stored. For real densities and l1 + l2 + lam even pc[-mu] = (-1)^mu conj(pc[mu]),
// so only mu >= 0 is summed; in pc[0] the terms of m1 and -m1 are complex conjugates.
static inline void shell_features(const lambda_plan &q, const SALTEDDescriptors &v1, const SALTEDDescriptors &u, const bool conj,
	const int iat, const int n1, const int n2, const int il, double *out, const size_t stride)
{
	const int lam = q.lam, l21 = q.l21;
	const cdouble *__restrict a = v1.block(iat, n1, (*q.llvec)[0][il]);
	const cdouble *__restrict v2_ptr = u.block(iat, n2, (*q.llvec)[1][il]);
	{
		const w3j_run &run = q.runs[static_cast<size_t>(il) * l21 + lam];
		const double *__restrict w = q.w3j->data() + run.w_off;
		const cdouble *__restrict av = a + run.im1_begin;
		const cdouble *__restrict b = v2_ptr + run.im2_begin;
		const int h = run.count / 2;
		const double sg = conj ? 1.0 : -1.0;
		double x = 0.0;
		for (int k = 0; k < h; ++k)
			x += w[k] * (av[k].real() * b[k].real() + sg * av[k].imag() * b[k].imag());
		out[static_cast<size_t>(lam) * stride] = q.cpm_re[lam] * (2.0 * x + w[h] * av[h].real() * b[h].real());
	}
	for (int mu = 1; mu <= lam; ++mu)
	{
		double sr = 0.0, si = 0.0;
		const w3j_run &run = q.runs[static_cast<size_t>(il) * l21 + lam + mu];
		const double *__restrict w = q.w3j->data() + run.w_off;
		const cdouble *__restrict av = a + run.im1_begin;
		const cdouble *__restrict b = v2_ptr + run.im2_begin;
		if (conj) [[likely]]
		{
			for (int k = 0; k < run.count; ++k)
			{
				const double ar = w[k] * av[k].real(), ai = w[k] * av[k].imag();
				sr += ar * b[k].real() + ai * b[k].imag();
				si += ai * b[k].real() - ar * b[k].imag();
			}
		}
		else
		{
			for (int k = 0; k < run.count; ++k)
			{
				const double ar = w[k] * av[k].real(), ai = w[k] * av[k].imag();
				sr += ar * b[k].real() - ai * b[k].imag();
				si += ar * b[k].imag() + ai * b[k].real();
			}
		}
		out[static_cast<size_t>(lam + mu) * stride] = q.cpm_re[lam + mu] * sr - q.cpm_im[lam + mu] * si;
		out[static_cast<size_t>(lam - mu) * stride] = q.cpm_re[lam - mu] * sr - q.cpm_im[lam - mu] * si;
	}
}

// sum of the squared features of one atom for one lambda, feature by feature
static double atom_norm2_direct(const lambda_plan &q, const SALTEDDescriptors &v1, const SALTEDDescriptors &u, const bool conj,
	const int iat, const int nrad1, const int nrad2)
{
	vec f(q.l21);
	double inner = 0.0;
	for (int n1 = 0; n1 < nrad1; ++n1)
		for (int n2 = 0; n2 < nrad2; ++n2)
			for (int il = 0; il < q.llmax; ++il)
			{
				shell_features(q, v1, u, conj, iat, n1, n2, il, f.data(), 1);
				for (int i = 0; i < q.l21; ++i)
					inner += f[i] * f[i];
			}
	return inner;
}

// sum of the squared features of one atom for one lambda, from its density matrices. For real
// densities, v[n,l,-m] = (-1)^m conj(v[n,l,m]), and l1 + l2 + lam even every feature is real, so
// with c2r unitary the sum is sum_mu |pc_mu|^2 = |pc_0|^2 + 2 sum_mu>0 |pc_mu|^2, and summed over
// n1 and n2 that is sum_mu sum_m1,m1' w w' A1[m1,m1'] A2[m1-mu,m1'-mu], A2 of u = conj(v1)
// (s2 = -1) or of v2
static double atom_norm2(const lambda_plan &q, const double *m1, const size_t dsz1, const double *m2, const size_t dsz2, const double s2)
{
	const int l21 = q.l21;
	const vec &w3j = *q.w3j;
	const ivec2 &llvec = *q.llvec;
	double inner = 0.0;
	for (int il = 0; il < q.llmax; ++il)
	{
		const int l1 = llvec[0][il], l2 = llvec[1][il];
		const int nm1 = 2 * l1 + 1, nm2 = 2 * l2 + 1;
		const size_t o1 = dm_offset(l1), o2 = dm_offset(l2);
		const double *a1r = m1 + o1, *a1i = m1 + dsz1 + o1;
		const double *a2r = m2 + o2, *a2i = m2 + dsz2 + o2;
		for (int mu = q.lam; mu < l21; ++mu)
		{
			const w3j_run &r = q.runs[static_cast<size_t>(il) * l21 + mu];
			double P = 0.0;
			for (int x = 0; x < r.count; ++x)
			{
				const double w = w3j[static_cast<size_t>(r.w_off) + x];
				const int i1 = (r.im1_begin + x) * nm1 + r.im1_begin;
				const int i2 = (r.im2_begin + x) * nm2 + r.im2_begin;
				for (int y = 0; y < r.count; ++y)
					P += w * w3j[static_cast<size_t>(r.w_off) + y] * (a1r[i1 + y] * a2r[i2 + y] - s2 * a1i[i1 + y] * a2i[i2 + y]);
			}
			inner += (mu == q.lam ? 1.0 : 2.0) * P;
		}
	}
	return inner;
}

// Said once per run. An atom with no neighbours still has an l = 0 descriptor, its own density
// being spherically symmetric; only the equivariant lam >= 1 parts vanish, so zeroing them
// leaves the atom spherical, which is the right answer for it.
static void warn_empty(const int empty_environments)
{
	static std::atomic<unsigned> warned_empty_environment{0};
	if (empty_environments > 0 && constants::first_this_run(warned_empty_environment))
	{
		std::cout << "WARNING: " << empty_environments << " atom(s) have no neighbour"
				  << " inside the descriptor cutoff.\n"
				  << "         Their environment singles out no direction, so their"
				  << " predicted density stays spherical.\n"
				  << "         Isolated solvent is the usual cause."
				  << std::endl;
	}
}

// normfact[k][atom] for every plan. The density matrices do not depend on lambda, so each
// atom's are built once into a per-thread buffer and contracted for all plans, instead of
// being stored for all atoms (natoms * 2 * sum (2l+1)^2 doubles): no faster, but that memory
// is never held. The work per atom is fixed, so the norm scales linearly with the atom count.
// An empty environment keeps 0, and so does an atom whose squared norm is below noise2 times
// tr(A1) tr(A2) / (2 lam + 1) over all l.
static vec2 plan_norms(const int natoms, const int nrad1, const int nrad2, const SALTEDDescriptors &v1, const SALTEDDescriptors &v2,
	const bool v2_is_conj_of_v1, const lambda_plan *const *plans, const size_t nplans, const double noise2)
{
	// u = conj(v1) has conj(A1) as its matrix; with nrad2 == nrad1 that is v1's
	// own (share) and the pair sum flips the sign of its imaginary part
	const SALTEDDescriptors &u = v2_is_conj_of_v1 ? v1 : v2;
	const bool share = v2_is_conj_of_v1 && nrad2 == nrad1;
	const double s2 = v2_is_conj_of_v1 ? -1.0 : 1.0;
	int lmax1 = 0, lmax2 = 0;
	for (size_t k = 0; k < nplans; ++k) { lmax1 = std::max(lmax1, plans[k]->lmax1); lmax2 = std::max(lmax2, plans[k]->lmax2); }
	if (share) lmax1 = lmax2 = std::max(lmax1, lmax2);
	err_checkf(2 * std::max(lmax1, lmax2) + 1 <= 64, "equicomb: l above 31 in the descriptors", std::cout);
	err_checkf(lmax1 < static_cast<int>(v1.offsets().size()) && lmax2 < static_cast<int>(u.offsets().size()),
		"equicomb: the shells need an l the descriptors do not hold", std::cout);
	const size_t dsz1 = dm_offset(lmax1 + 1), dsz2 = dm_offset(lmax2 + 1);
	vec2 normfact(nplans, vec(natoms, 0.0));
	ivec direct(nplans, 0), zeroed(nplans, 0);
	int empty = 0;
#pragma omp parallel
	{
		vec m1(2 * dsz1), m2(share ? 0 : 2 * dsz2);
		vec x(2 * static_cast<size_t>(2 * std::max(lmax1, lmax2) + 1) * std::max(nrad1, nrad2));
		double *xr = x.data(), *xi = xr + x.size() / 2;
#pragma omp for schedule(dynamic, 8)
		for (int iat = 0; iat < natoms; ++iat)
		{
			atom_density(v1, nrad1, lmax1, iat, dsz1, m1.data(), xr, xi);
			if (!share)
				atom_density(u, nrad2, lmax2, iat, dsz2, m2.data(), xr, xi);
			const double *M2 = share ? m1.data() : m2.data();
			double tr1[32], tr2[32], tot1 = 0.0, tot2 = 0.0;
			for (int l = 0; l <= lmax1; ++l)
			{
				tr1[l] = 0.0;
				for (int a = 0; a <= 2 * l; ++a) tr1[l] += m1[dm_offset(l) + a * (2 * l + 1) + a];
				tot1 += tr1[l];
			}
			for (int l = 0; l <= lmax2; ++l)
			{
				tr2[l] = 0.0;
				for (int a = 0; a <= 2 * l; ++a) tr2[l] += M2[dm_offset(l) + a * (2 * l + 1) + a];
				tot2 += tr2[l];
			}
			for (size_t k = 0; k < nplans; ++k)
			{
				const lambda_plan &q = *plans[k];
				double inner = atom_norm2(q, m1.data(), dsz1, M2, dsz2, s2);
				// sum_il tr(A1[l1]) tr(A2[l2]) / (2 lam + 1) bounds inner from above and is the size of
				// the terms the density-matrix sum adds; where they cancel to a small fraction of it,
				// as for features that nearly vanish by site symmetry, the features are summed directly
				double bound = 0.0;
				for (int il = 0; il < q.llmax; ++il)
					bound += tr1[(*q.llvec)[0][il]] * tr2[(*q.llvec)[1][il]];
				if (inner < 1e-5 * bound / q.l21)
				{
					inner = atom_norm2_direct(q, v1, u, v2_is_conj_of_v1, iat, nrad1, nrad2);
#pragma omp atomic
					++direct[k];
				}
				// Features that vanish by site symmetry come out as rounding noise of the positions
				// and the descriptors; normalised they would point anywhere, so they are zeroed.
				// An empty environment gives an all-zero descriptor, so inner is 0 and
				// 1/sqrt(inner) is +inf, making every feature NaN. Zero is the meaningful
				// answer: the kernel contributes nothing and the atom keeps the species
				// average the model adds separately.
				if (inner > 0.0 && inner <= noise2 * tot1 * tot2 / q.l21)
				{
#pragma omp atomic
					++zeroed[k];
				}
				else if (inner > 0.0) [[likely]]
					normfact[k][iat] = 1.0 / sqrt(inner);
			}
			// no neighbour: only the atom's own spherical density, every l >= 1 trace exactly zero
			if (tot1 == tr1[0] && tot2 == tr2[0])
			{
#pragma omp atomic
				++empty;
			}
		}
	}
	if (ProgressBar::report_counts)
	{
		std::cout << "[equicomb] direct norms of " << natoms << " atoms, per lambda:";
		for (size_t k = 0; k < nplans; ++k)
			std::cout << " " << plans[k]->lam << ":" << direct[k];
		std::cout << "; zeroed as rounding noise:";
		for (size_t k = 0; k < nplans; ++k)
			std::cout << " " << plans[k]->lam << ":" << zeroed[k];
		std::cout << std::endl;
	}
	warn_empty(empty);
	return normfact;
}

vec2 equicomb_norms(int natoms, int nrad1, int nrad2,
	const SALTEDDescriptors &v1, const SALTEDDescriptors &v2,
	const std::vector<const vec *> &w3j, const std::vector<ivec2> &llvec, const std::vector<cvec2> &c2r,
	bool v2_is_conj_of_v1, double noise2)
{
	err_checkf(w3j.size() == llvec.size() && c2r.size() == llvec.size(), "equicomb_norms: one w3j, llvec and c2r per lambda", std::cout);
	std::vector<lambda_plan> plans;
	plans.reserve(llvec.size());
	for (size_t lam = 0; lam < llvec.size(); ++lam)
		plans.push_back(make_plan(*w3j[lam], llvec[lam], static_cast<int>(lam), c2r[lam]));
	std::vector<const lambda_plan *> ptr;
	for (const lambda_plan &q : plans) ptr.push_back(&q);
	return plan_norms(natoms, nrad1, nrad2, v1, v2, v2_is_conj_of_v1, ptr.data(), ptr.size(), noise2);
}

// BE AWARE, THAT V2 IS ALREADY ASSUMED TO BE CONJUGATED!!!!!
void equicomb(int natoms, int nrad1, int nrad2,
			  const SALTEDDescriptors &v1,
			  const SALTEDDescriptors &v2,
			  const vec &w3j,
			  const ivec2 &llvec, const int &lam,
			  const cvec2 &c2r, const int &featsize,
			  const int &nfps, const std::vector<int64_t> &vfps,
			  vec &p,
			  bool v2_is_conj_of_v1,
			  const double *norms)
{
	if (natoms < 0 || nrad1 < 0 || nrad2 < 0 || lam < 0 || featsize < 0 || nfps < 0)
	{
		throw std::invalid_argument("equicomb: negative dimensions are not allowed");
	}
	if (natoms == 0 || nrad1 == 0 || nrad2 == 0 || featsize == 0 || nfps == 0)
	{
		return;
	}

	const long long l21_ll = 2LL * static_cast<long long>(lam) + 1LL;
	if (l21_ll <= 0LL || l21_ll > static_cast<long long>(std::numeric_limits<int>::max()))
	{
		throw std::overflow_error("equicomb: invalid lam leads to invalid 2*lam+1");
	}
	const int l21 = static_cast<int>(l21_ll);
	const int llmax = (int)llvec[0].size();

	const size_t required_p = static_cast<size_t>(natoms) * static_cast<size_t>(l21) * static_cast<size_t>(nfps);
	if (p.size() < required_p)
	{
		throw std::out_of_range("equicomb: output buffer p is smaller than required size");
	}

	// Sizes below come from the model file, not the structure: ifeat advances
	// nrad1*nrad2*llmax into a featsize buffer, w3j is walked unbounded and vfps
	// indexes ptemp. A structure missing a species the model knows makes those
	// disagree, so check here rather than run off the end inside the threads.
	// v2 is not filled when the two descriptor sets match, so read v2_src, never v2.
	const SALTEDDescriptors &v2_src = v2_is_conj_of_v1 ? v1 : v2;
	const size_t shells = static_cast<size_t>(nrad1) * nrad2 * llmax;
	err_checkf(shells <= static_cast<size_t>(featsize), "equicomb: featsize " + std::to_string(featsize) +
		" is smaller than nrad1*nrad2*llmax " + std::to_string(shells), std::cout);
	for (int chk = 0; chk < nfps; ++chk)
		err_checkf(vfps[chk] >= 0 && static_cast<size_t>(vfps[chk]) < static_cast<size_t>(featsize),
			"equicomb: vfps entry " + std::to_string(chk) + " is outside featsize", std::cout);

	const lambda_plan plan = make_plan(w3j, llvec, lam, c2r);
	const std::vector<w3j_run> &runs = plan.runs;
	const size_t total_terms = plan.total_terms;
	if (ProgressBar::report_counts)
	{
		std::cout << "[equicomb] lam " << lam
				  << ", " << total_terms << " wigner terms"
				  << ", natoms " << natoms << ", nrad1 " << nrad1 << ", nrad2 " << nrad2
				  << ", llmax " << llmax << ", l21 " << l21
				  << ", featsize " << featsize << ", nfps " << nfps << std::endl;
	}

	//Timed from here so both throughput rows include the norm; counted as the nfps
	//features actually built
	const _time_point eq_t0 = get_time();
	vec normfact(natoms, 0.0);
	// Features past nrad1*nrad2*llmax stay zero, as they were when ptemp was built in full
	const int shells_i = static_cast<int>(shells);

	// The predictor hands in all lambda's norms from one pass (equicomb_norms); a caller
	// without them gets this lambda's alone
	if (norms)
		std::copy(norms, norms + natoms, normfact.begin());
	else
	{
		const lambda_plan *one = &plan;
		normfact = plan_norms(natoms, nrad1, nrad2, v1, v2, v2_is_conj_of_v1, &one, 1, 0.0)[0];
	}
	if (ProgressBar::report_counts)
		std::cout << "[equicomb] lam " << lam << ": norm " << get_msec(eq_t0, get_time()) << " ms" << std::endl;

#ifdef NOSPHERA2_USE_GPU
	//The device builds the selected features with the norm above; it falls through to the
	//CPU loop below if no device is present or it will not fit.
	if (g_equicomb_use_gpu && salted_gpu_available())
	{
		ivec flat_runs(static_cast<size_t>(llmax) * l21 * 4);
		for (size_t r = 0; r < runs.size(); ++r) {
			flat_runs[r * 4 + 0] = runs[r].im1_begin;
			flat_runs[r * 4 + 1] = runs[r].im2_begin;
			flat_runs[r * 4 + 2] = runs[r].count;
			flat_runs[r * 4 + 3] = runs[r].w_off;
		}
		ivec fps(vfps.begin(), vfps.begin() + nfps);
		const SALTEDDescriptors &v2_gpu = v2_is_conj_of_v1 ? v1 : v2;
		salted_gpu_problem q;
		q.natoms = natoms; q.nrad1 = nrad1; q.nrad2 = nrad2; q.llmax = llmax;
		q.lam = lam; q.l21 = l21; q.shells = shells_i; q.nfps = nfps;
		q.v2_is_conj_of_v1 = v2_is_conj_of_v1;
		q.v1_values = reinterpret_cast<const double *>(v1.values().data());
		q.v1_offsets = v1.offsets().data();
		q.v1_nchannels = v1.nchannels();
		q.v1_noff = static_cast<int>(v1.offsets().size());
		q.v1_len_doubles = static_cast<long long>(v1.values().size()) * 2;
		q.v2_values = reinterpret_cast<const double *>(v2_gpu.values().data());
		q.v2_offsets = v2_gpu.offsets().data();
		q.v2_nchannels = v2_gpu.nchannels();
		q.v2_noff = static_cast<int>(v2_gpu.offsets().size());
		q.v2_len_doubles = static_cast<long long>(v2_gpu.values().size()) * 2;
		q.w3j = w3j.data(); q.w3j_len = static_cast<long long>(w3j.size());
		q.llvec0 = llvec[0].data(); q.llvec1 = llvec[1].data();
		q.runs = flat_runs.data(); q.cpm_re = plan.cpm_re.data(); q.cpm_im = plan.cpm_im.data();
		q.vfps = fps.data(); q.normfact = normfact.data(); q.p = p.data();
		const bool gpu_ok = salted_gpu_equicomb(q);
		if (gpu_ok)
			throughput::record("SALTED equicomb", true,
				throughput::flops_equicomb(natoms, nfps, 1, 1, l21),
				get_msec(eq_t0, get_time()));
		//Once per run, not once per lambda: nine identical lines say nothing extra.
		//Printed before the progress bar exists, whose carriage returns would eat it.
		static std::atomic<unsigned> announced{0};
		if (!constants::hide_gpu_notes && constants::first_this_run(announced)) {
			std::cout << "GPU in use: SALTED descriptors on "
					  << (gpu_ok ? "the device (double precision)" : "the CPU - device unavailable") << std::endl;
		}
		if (ProgressBar::report_counts)
			std::cout << "[equicomb] lam " << lam << ": GPU " << (gpu_ok ? "used" : "refused, using the CPU loop")
					  << ", " << get_msec(eq_t0, get_time()) << " ms with the norm" << std::endl;
		if (gpu_ok)
			return;
	}
#endif

	// Only the nfps selected features are built, as in shell_features; only normfact comes
	// from the norm above. They are grouped by il so that each (il, mu) run is set up once
	// for all of them: features start[il]..start[il + 1] have channels sn1, sn2 and go to
	// column scol of p. The columns pad past nrad1*nrad2*llmax get zeros, so every element of p is written.
	ivec start(llmax + 1, 0);
	for (int i = 0; i < nfps; ++i)
		if (vfps[i] < shells_i) ++start[vfps[i] % llmax + 1];
	for (int il = 0; il < llmax; ++il) start[il + 1] += start[il];
	ivec sn1(start[llmax]), sn2(start[llmax]), scol(start[llmax]), next(start.begin(), start.end() - 1), pad;
	for (int i = 0; i < nfps; ++i)
	{
		const int f = static_cast<int>(vfps[i]);
		if (f >= shells_i)
		{
			pad.push_back(i);
			continue;
		}
		const int k = next[f % llmax]++;
		sn1[k] = f / (llmax * nrad2); sn2[k] = (f / llmax) % nrad2; scol[k] = i;
	}
	const double sg = v2_is_conj_of_v1 ? 1.0 : -1.0;
	// Scoped: the bar rewinds to its own line when it is destroyed, so anything
	// printed after the loop but before that is silently overwritten.
	{
	ProgressBar pb(natoms, 60, "#", " ", "Calculating descriptors for l = " + toString(lam));
#pragma omp parallel for schedule(dynamic, 1)
	for (int iat = 0; iat < natoms; ++iat)
	{
		double *pa = p.data() + static_cast<size_t>(iat) * l21 * nfps;
		const double nf = normfact[iat];
		for (int imu = 0; imu < l21; ++imu)
			for (size_t k = 0; k < pad.size(); ++k)
				pa[static_cast<size_t>(imu) * nfps + pad[k]] = 0.0;
		for (int il = 0; il < llmax; ++il)
		{
			const int l1 = llvec[0][il], l2 = llvec[1][il];
			{
				const w3j_run &run = runs[static_cast<size_t>(il) * l21 + lam];
				const double *__restrict w = w3j.data() + run.w_off;
				const int h = run.count / 2;
				const double c = plan.cpm_re[lam] * nf;
				double *__restrict out = pa + static_cast<size_t>(lam) * nfps;
				for (int k = start[il]; k < start[il + 1]; ++k)
				{
					const cdouble *__restrict av = v1.block(iat, sn1[k], l1) + run.im1_begin;
					const cdouble *__restrict b = v2_src.block(iat, sn2[k], l2) + run.im2_begin;
					double x = 0.0;
					for (int t = 0; t < h; ++t)
						x += w[t] * (av[t].real() * b[t].real() + sg * av[t].imag() * b[t].imag());
					out[scol[k]] = c * (2.0 * x + w[h] * av[h].real() * b[h].real());
				}
			}
			for (int mu = 1; mu <= lam; ++mu)
			{
				const w3j_run &run = runs[static_cast<size_t>(il) * l21 + lam + mu];
				const double *__restrict w = w3j.data() + run.w_off;
				const double cpr = plan.cpm_re[lam + mu] * nf, cpi = plan.cpm_im[lam + mu] * nf;
				const double cmr = plan.cpm_re[lam - mu] * nf, cmi = plan.cpm_im[lam - mu] * nf;
				double *__restrict op = pa + static_cast<size_t>(lam + mu) * nfps, *__restrict om = pa + static_cast<size_t>(lam - mu) * nfps;
				for (int k = start[il]; k < start[il + 1]; ++k)
				{
					const cdouble *__restrict av = v1.block(iat, sn1[k], l1) + run.im1_begin;
					const cdouble *__restrict b = v2_src.block(iat, sn2[k], l2) + run.im2_begin;
					double sr = 0.0, si = 0.0;
					for (int t = 0; t < run.count; ++t)
					{
						const double ar = w[t] * av[t].real(), ai = w[t] * av[t].imag();
						sr += ar * b[t].real() + sg * ai * b[t].imag();
						si += ai * b[t].real() - sg * ar * b[t].imag();
					}
					op[scol[k]] = cpr * sr - cpi * si;
					om[scol[k]] = cmr * sr - cmi * si;
				}
			}
		}
	}
	}

	throughput::record("SALTED equicomb", false,
		throughput::flops_equicomb(natoms, nfps, 1, 1, l21),
		get_msec(eq_t0, get_time()));
}

void equicomb(int natoms, int nrad1, int nrad2,
			  const SALTEDDescriptors &v1,
			  const SALTEDDescriptors &v2,
			  vec &w3j, int llmax,
			  ivec2 &llvec, int lam,
			  cvec2 &c2r, int featsize,
			  vec &p,
			  bool v2_is_conj_of_v1)
{
	if (natoms < 0 || nrad1 < 0 || nrad2 < 0 || llmax < 0 || lam < 0 || featsize < 0)
	{
		throw std::invalid_argument("equicomb: negative dimensions are not allowed");
	}
	if (natoms == 0 || nrad1 == 0 || nrad2 == 0 || llmax == 0 || featsize == 0)
	{
		return;
	}

	const long long l21_ll = 2LL * static_cast<long long>(lam) + 1LL;
	if (l21_ll <= 0LL || l21_ll > static_cast<long long>(std::numeric_limits<int>::max()))
	{
		throw std::overflow_error("equicomb: invalid lam leads to invalid 2*lam+1");
	}
	const int l21 = static_cast<int>(l21_ll);

	const size_t required_p = static_cast<size_t>(natoms) * static_cast<size_t>(l21) * static_cast<size_t>(featsize);
	if (p.size() < required_p)
	{
		throw std::out_of_range("equicomb: output buffer p is smaller than required size");
	}

	int iat, n1, n2, il, imu, im1, im2, i, j, ifeat, iwig, l1, l2, mu, m1, m2;
	double inner, normfact;

	// default(none) means every name the region touches must be listed, read-only
	// parameters included. clang enforces that; MSVC's OpenMP 2.0 does not, so a
	// missing name builds clean on Windows and breaks the macOS and Linux jobs.

#pragma omp parallel for private(iat, n1, n2, il, imu, im1, im2, i, j, ifeat, iwig, l1, l2, mu, m1, m2, inner, normfact) default(none) shared(natoms, nrad1, nrad2, v1, v2, v2_is_conj_of_v1, w3j, llmax, llvec, lam, l21, c2r, p, featsize, constants::cnull)
	for (iat = 0; iat < natoms; ++iat)
	{
		vec2 ptemp(l21, vec(featsize, 0.0));
		cvec pcmplx(l21, constants::cnull);
		vec preal(l21, 0.0);
		inner = 0.0;
		ifeat = 0;
		for (n1 = 0; n1 < nrad1; ++n1)
		{
			for (n2 = 0; n2 < nrad2; ++n2)
			{
				iwig = 0;
				for (il = 0; il < llmax; ++il)
				{
					l1 = llvec[0][il];
					l2 = llvec[1][il];

					fill(pcmplx.begin(), pcmplx.end(), constants::cnull);

					for (imu = 0; imu < l21; ++imu)
					{
						mu = imu - lam;
						for (im1 = 0; im1 < 2 * l1 + 1; ++im1)
						{
							m1 = im1 - l1;
							m2 = m1 - mu;
							if (abs(m2) <= l2)
							{
								im2 = m2 + l2;
								// v2 is conj(v1) elementwise when the two descriptor sets match
								pcmplx[imu] += w3j[iwig] * v1.block(iat, n1, l1)[im1] *
									(v2_is_conj_of_v1 ? std::conj(v1.block(iat, n2, l2)[im2])
													  : v2.block(iat, n2, l2)[im2]);
								iwig++;
							}
						}
					}

					fill(preal.begin(), preal.end(), 0.0);
					for (i = 0; i < l21; ++i)
					{
						for (j = 0; j < l21; ++j)
						{
							preal[i] += real(c2r[i][j] * pcmplx[j]);
						}
						inner += preal[i] * preal[i];
						ptemp[i][ifeat] = preal[i];
					}
					ifeat++;
				}
			}
		}
		// See the sparsified overload: an empty environment makes this zero and every feature NaN
		normfact = sqrt(inner);
		const double inv_normfact = (normfact > 0.0) ? (1.0 / normfact) : 0.0;
		for (ifeat = 0; ifeat < featsize; ++ifeat)
		{
			for (imu = 0; imu < l21; ++imu)
			{
				const size_t out_idx = static_cast<size_t>(iat) * static_cast<size_t>(l21) * static_cast<size_t>(featsize)
					+ static_cast<size_t>(imu) * static_cast<size_t>(featsize)
					+ static_cast<size_t>(ifeat);
				p[out_idx] = ptemp[imu][ifeat] * inv_normfact;
			}
		}
	}
}
