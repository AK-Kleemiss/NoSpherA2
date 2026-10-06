#include "pch.h"
#include "SALTED_equicomb.h"

#if defined(__APPLE__)
// On macOS we�re using Accelerate for BLAS/LAPACK
#include <Accelerate/Accelerate.h>
#else
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

// BE AWARE, THAT V2 IS ALREADY ASSUMED TO BE CONJUGATED!!!!!
void equicomb(int natoms, int nrad1, int nrad2,
			  const SALTEDDescriptors &v1,
			  const SALTEDDescriptors &v2,
			  const vec &w3j,
			  const ivec2 &llvec, const int &lam,
			  const cvec2 &c2r, const int &featsize,
			  const int &nfps, const std::vector<int64_t> &vfps,
			  vec &p,
			  bool v2_is_conj_of_v1)
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

	std::fill(p.begin(), p.begin() + static_cast<std::ptrdiff_t>(required_p), 0.0);

	// Which (im1, im2) pairs contribute is set by |im1 - mu| <= l2, which depends
	// only on (il, imu, lam), never on n1 or n2. That test selects a CONTIGUOUS run
	// of im1, and im2 = im1 - mu + l2 advances in lockstep with it, so a first
	// index, a length and an offset into w3j describe a group completely. w3j is
	// consumed in exactly (il, imu, im1) order, so a group's weights are a
	// contiguous slice of it and need no copy here.
	struct w3j_run { int im1_begin, im2_begin, count, w_off; };
	std::vector<w3j_run> runs(static_cast<size_t>(llmax) * l21, w3j_run{0, 0, 0, 0});
	size_t total_terms = 0;
	{
		int w_idx = 0;
		for (int til = 0; til < llmax; ++til)
		{
			const int tl1 = llvec[0][til], tl2 = llvec[1][til];
			err_checkf(tl1 >= 0 && tl2 >= 0, "equicomb: negative angular momentum in llvec", std::cout);
			for (int timu = 0; timu < l21; ++timu)
			{
				const int tmu = timu - lam + tl1;
				const int lo = std::max(0, tmu - tl2);
				const int hi = std::min(2 * tl1, tmu + tl2);
				const int cnt = (hi >= lo) ? (hi - lo + 1) : 0;
				runs[static_cast<size_t>(til) * l21 + timu] = { lo, lo - tmu + tl2, cnt, w_idx };
				w_idx += cnt;
			}
		}
		total_terms = static_cast<size_t>(w_idx);
	}
	// w3j is consumed once per (n1,n2) shell pair, one entry per surviving m2, which is what the runs add up to
	err_checkf(w3j.size() >= total_terms, "equicomb: w3j holds " + std::to_string(w3j.size()) +
		" entries, the shell loop consumes " + std::to_string(total_terms), std::cout);

	// The complex-to-real matrix is a mirror-pair transform: row i couples only
	// column i and column l21-1-i, so every row holds exactly two nonzeros (the
	// middle row holds one). Dropping the exact zeros leaves a finite sum bit for
	// bit the same as long as the survivors stay in ascending column order, which
	// is the order the dense loop added them in. The structure is read off c2r
	// rather than assumed; anything not two-per-row is an error.
	struct c2r_entry { int j; double re, im; };
	std::vector<c2r_entry> c2r_nz(static_cast<size_t>(l21) * 2, c2r_entry{0, 0.0, 0.0});
	ivec c2r_cnt(l21, 0);
	for (int i2 = 0; i2 < l21; ++i2)
	{
		int cnt = 0;
		for (int j2 = 0; j2 < l21; ++j2)
		{
			const cdouble &e = c2r[i2][j2];
			if (e.real() == 0.0 && e.imag() == 0.0) continue;
			err_checkf(cnt < 2, "equicomb: complex-to-real row " + std::to_string(i2) + " has more than two non-zeros", std::cout);
			c2r_nz[static_cast<size_t>(i2) * 2 + cnt] = { j2, e.real(), e.imag() };
			++cnt;
		}
		c2r_cnt[i2] = cnt;
	}
	if (ProgressBar::report_counts)
	{
		std::cout << "[equicomb] lam " << lam
				  << ", " << total_terms << " wigner terms"
				  << ", natoms " << natoms << ", nrad1 " << nrad1 << ", nrad2 " << nrad2
				  << ", llmax " << llmax << ", l21 " << l21
				  << ", featsize " << featsize << ", nfps " << nfps << std::endl;
	}

	// Every atom is normalised over all featsize features, but only nfps of them are
	// kept, and the norm needs none of the others. Feature (n1, n2, il) is Re(q) with
	// q = c2r pc and pc[mu] = sum_m1 w v1[n1,l1,m1] u[n2,l2,m1-mu], where u is the factor
	// the loop multiplies by (conj(v1), or v2 as stored). Then
	//   sum_i Re(q_i)^2 = (q^H q + Re q^T q) / 2
	//                   = Re(sum K[mu,mu'] pc_mu conj(pc_mu') + sum G[mu,mu'] pc_mu pc_mu') / 2
	// with K = sum_i c2r[i,mu] conj(c2r[i,mu']) and G = sum_i c2r[i,mu] c2r[i,mu'], and
	// summed over n1 and n2 both products factor into per-atom, per-l density matrices
	// A[m,m'] = sum_n x_m conj(x_m') and B[m,m'] = sum_n x_m x_m' of v1 and of u. That is
	// one pass over m pairs per shell instead of nrad1*nrad2 feature evaluations; the
	// value is the same up to rounding (SALTED's equicombsparse_numba, a7fbc3a).
	const int lmax1 = *std::max_element(llvec[0].begin(), llvec[0].end());
	const int lmax2 = *std::max_element(llvec[1].begin(), llvec[1].end());
	vec K_re(static_cast<size_t>(l21) * l21, 0.0), K_im(K_re), G_re(K_re), G_im(K_re);
	for (int a = 0; a < l21; ++a)
		for (int b = 0; b < l21; ++b)
		{
			cdouble k = constants::cnull, g = constants::cnull;
			for (int i2 = 0; i2 < l21; ++i2)
			{
				k += c2r[i2][a] * std::conj(c2r[i2][b]);
				g += c2r[i2][a] * c2r[i2][b];
			}
			const size_t ab = static_cast<size_t>(a) * l21 + b;
			K_re[ab] = k.real(); K_im[ab] = k.imag(); G_re[ab] = g.real(); G_im[ab] = g.imag();
		}
	// A and B of one atom for l = 0..lmax, packed per l as (2l+1)^2 blocks, real and
	// imaginary parts apart: plain double loops vectorise, std::complex products without
	// -ffast-math do not. When u = conj(v1) its matrices are the conjugates of v1's, so the
	// pair sum flips the sign of their imaginary parts, and with nrad2 == nrad1 they are
	// v1's own and are not built twice.
	auto dm_offset = [](int l) { size_t o = 0; for (int k = 0; k < l; ++k) o += static_cast<size_t>(2 * k + 1) * (2 * k + 1); return o; };
	struct dmat { vec are, aim, bre, bim; };
	auto density_matrices = [&](const SALTEDDescriptors &d, int nch, int lmax, int iat, dmat &m) {
		std::fill(m.are.begin(), m.are.end(), 0.0);
		std::fill(m.aim.begin(), m.aim.end(), 0.0);
		std::fill(m.bre.begin(), m.bre.end(), 0.0);
		std::fill(m.bim.begin(), m.bim.end(), 0.0);
		double xr[64], xi[64];
		for (int l = 0; l <= lmax; ++l)
		{
			const int nm = 2 * l + 1;
			const size_t o = dm_offset(l);
			double *__restrict are = m.are.data() + o, *__restrict aim = m.aim.data() + o;
			double *__restrict bre = m.bre.data() + o, *__restrict bim = m.bim.data() + o;
			for (int n = 0; n < nch; ++n)
			{
				const cdouble *x = d.block(iat, n, l);
				for (int a = 0; a < nm; ++a) { xr[a] = x[a].real(); xi[a] = x[a].imag(); }
				for (int a = 0; a < nm; ++a)
					for (int b = 0; b < nm; ++b)
					{
						are[a * nm + b] += xr[a] * xr[b] + xi[a] * xi[b];
						aim[a * nm + b] += xi[a] * xr[b] - xr[a] * xi[b];
						bre[a * nm + b] += xr[a] * xr[b] - xi[a] * xi[b];
						bim[a * nm + b] += xi[a] * xr[b] + xr[a] * xi[b];
					}
			}
		}
	};
	err_checkf(2 * std::max(lmax1, lmax2) + 1 <= 64, "equicomb: l above 31 in llvec", std::cout);
	//Timed from here so both throughput rows include the norm; counted as the nfps
	//features actually built
	const _time_point eq_t0 = get_time();
	// Shared, v1's matrices stand in for u's, so they then reach lmax2 as well
	const bool share = v2_is_conj_of_v1 && nrad2 == nrad1;
	const int lm1 = share ? std::max(lmax1, lmax2) : lmax1;
	const size_t dm1_size = dm_offset(lm1 + 1), dm2_size = dm_offset(lmax2 + 1);
	const double s2 = v2_is_conj_of_v1 ? -1.0 : 1.0;
	vec normfact(natoms, 0.0);
	int empty_environments = 0;
#pragma omp parallel
	{
		dmat m1{vec(dm1_size), vec(dm1_size), vec(dm1_size), vec(dm1_size)};
		dmat m2_own;
		if (!share)
			m2_own = dmat{vec(dm2_size), vec(dm2_size), vec(dm2_size), vec(dm2_size)};
		const dmat &m2 = share ? m1 : m2_own;
#pragma omp for schedule(dynamic, 1) reduction(+ : empty_environments)
		for (int iat = 0; iat < natoms; ++iat)
		{
			density_matrices(v1, nrad1, lm1, iat, m1);
			if (!share)
				density_matrices(v2_src, nrad2, lmax2, iat, m2_own);
			double inner = 0.0;
			for (int il = 0; il < llmax; ++il)
			{
				const int l1 = llvec[0][il], l2 = llvec[1][il];
				const int nm1 = 2 * l1 + 1, nm2 = 2 * l2 + 1;
				const size_t o1 = dm_offset(l1), o2 = dm_offset(l2);
				const double *a1r = m1.are.data() + o1, *a1i = m1.aim.data() + o1, *b1r = m1.bre.data() + o1, *b1i = m1.bim.data() + o1;
				const double *a2r = m2.are.data() + o2, *a2i = m2.aim.data() + o2, *b2r = m2.bre.data() + o2, *b2i = m2.bim.data() + o2;
				double t = 0.0;
				for (int mu = 0; mu < l21; ++mu)
				{
					const w3j_run &r = runs[static_cast<size_t>(il) * l21 + mu];
					if (r.count == 0) continue;
					for (int mu2 = 0; mu2 < l21; ++mu2)
					{
						const size_t mm = static_cast<size_t>(mu) * l21 + mu2;
						const w3j_run &r2 = runs[static_cast<size_t>(il) * l21 + mu2];
						if ((K_re[mm] == 0.0 && K_im[mm] == 0.0 && G_re[mm] == 0.0 && G_im[mm] == 0.0) || r2.count == 0) continue;
						double Pr = 0.0, Pi = 0.0, Qr = 0.0, Qi = 0.0;
						for (int x = 0; x < r.count; ++x)
						{
							const double w = w3j[static_cast<size_t>(r.w_off) + x];
							const int i1 = (r.im1_begin + x) * nm1 + r2.im1_begin;
							const int i2 = (r.im2_begin + x) * nm2 + r2.im2_begin;
							for (int y = 0; y < r2.count; ++y)
							{
								const double ww = w * w3j[static_cast<size_t>(r2.w_off) + y];
								const double ar = a2r[i2 + y], ai = s2 * a2i[i2 + y];
								const double br = b2r[i2 + y], bi = s2 * b2i[i2 + y];
								Pr += ww * (a1r[i1 + y] * ar - a1i[i1 + y] * ai);
								Pi += ww * (a1r[i1 + y] * ai + a1i[i1 + y] * ar);
								Qr += ww * (b1r[i1 + y] * br - b1i[i1 + y] * bi);
								Qi += ww * (b1r[i1 + y] * bi + b1i[i1 + y] * br);
							}
						}
						t += K_re[mm] * Pr - K_im[mm] * Pi + G_re[mm] * Qr - G_im[mm] * Qi;
					}
				}
				inner += 0.5 * t;
			}
			// An empty environment gives an all-zero descriptor, so inner is 0 and
			// 1/sqrt(inner) is +inf, making every feature NaN. Zero is the meaningful
			// answer: the kernel contributes nothing and the atom keeps the species
			// average the model adds separately.
			if (inner > 0.0) [[likely]]
				normfact[iat] = 1.0 / sqrt(inner);
			else
				++empty_environments;
		}
	}

	// Said once per run, on the first lambda that sees it - NOT gated on lam == 0.
	// An atom with no neighbours still has an l = 0 descriptor, its own density being
	// spherically symmetric; only the equivariant lam >= 1 parts vanish, so zeroing
	// them leaves the atom spherical, which is the right answer for it.
	static bool warned_empty_environment = false;
	if (empty_environments > 0 && !warned_empty_environment)
	{
		warned_empty_environment = true;
		std::cout << "WARNING: " << empty_environments << " atom(s) have no neighbour"
				  << " inside the descriptor cutoff.\n"
				  << "         Their environment singles out no direction, so their"
				  << " predicted density stays spherical.\n"
				  << "         Isolated solvent is the usual cause."
				  << std::endl;
	}
	if (ProgressBar::report_counts)
		std::cout << "[equicomb] lam " << lam << ": norm " << get_msec(eq_t0, get_time()) << " ms" << std::endl;
	// Features past nrad1*nrad2*llmax stay zero, as they were when ptemp was built in full
	const int shells_i = static_cast<int>(shells);

#ifdef NOSPHERA2_USE_GPU
	//The device reproduces this walk exactly; it falls through to the CPU loop below if
	//no device is present or it will not fit.
	if (g_equicomb_use_gpu)
	{
		ivec flat_runs(static_cast<size_t>(llmax) * l21 * 4);
		for (size_t r = 0; r < runs.size(); ++r) {
			flat_runs[r * 4 + 0] = runs[r].im1_begin;
			flat_runs[r * 4 + 1] = runs[r].im2_begin;
			flat_runs[r * 4 + 2] = runs[r].count;
			flat_runs[r * 4 + 3] = runs[r].w_off;
		}
		ivec cols(static_cast<size_t>(l21) * 2, 0);
		vec cre(static_cast<size_t>(l21) * 2, 0.0), cim(static_cast<size_t>(l21) * 2, 0.0);
		for (int i2 = 0; i2 < l21; ++i2)
			for (int k2 = 0; k2 < 2; ++k2) {
				cols[i2 * 2 + k2] = c2r_nz[static_cast<size_t>(i2) * 2 + k2].j;
				cre[i2 * 2 + k2] = c2r_nz[static_cast<size_t>(i2) * 2 + k2].re;
				cim[i2 * 2 + k2] = c2r_nz[static_cast<size_t>(i2) * 2 + k2].im;
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
		q.runs = flat_runs.data(); q.c2r_cols = cols.data();
		q.c2r_re = cre.data(); q.c2r_im = cim.data(); q.c2r_cnt = c2r_cnt.data();
		q.vfps = fps.data(); q.normfact = normfact.data(); q.p = p.data();
		const bool gpu_ok = salted_gpu_equicomb(q);
		if (gpu_ok)
			throughput::record("SALTED equicomb", true,
				throughput::flops_equicomb(natoms, nfps, 1, 1, l21),
				get_msec(eq_t0, get_time()));
		//Once per run, not once per lambda: nine identical lines say nothing extra.
		//Printed before the progress bar exists, whose carriage returns would eat it.
		static bool announced = false;
		if (!announced && !constants::hide_gpu_notes) {
			announced = true;
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

	// Only the nfps selected features are built. Each is the same arithmetic the full
	// walk did (w * v1 rounded to a double first, then the run sum, then c2r), so the
	// features are bit-identical to it; only normfact comes from the norm above.
	// Scoped: the bar rewinds to its own line when it is destroyed, so anything
	// printed after the loop but before that is silently overwritten.
	{
	ProgressBar pb(natoms, 60, "#", " ", "Calculating descriptors for l = " + toString(lam));
#pragma omp parallel
	{
		vec pre(l21), pim(l21);
#pragma omp for schedule(dynamic, 1)
		for (int iat = 0; iat < natoms; ++iat)
		{
			const size_t offset = static_cast<size_t>(iat) * l21 * nfps;
			for (int i = 0; i < nfps; ++i)
			{
				const int f = static_cast<int>(vfps[i]);
				if (f >= shells_i) continue;
				const int il = f % llmax, n2 = (f / llmax) % nrad2, n1 = f / (llmax * nrad2);
				const cdouble *__restrict a = v1.block(iat, n1, llvec[0][il]);
				// v2 is conj(v1) when the two descriptor sets are the same; the sign is
				// applied by choosing between the two forms of the inner loop
				const cdouble *__restrict v2_ptr = v2_src.block(iat, n2, llvec[1][il]);
				for (int imu = 0; imu < l21; ++imu)
				{
					double sr = 0.0, si = 0.0;
					const w3j_run &run = runs[static_cast<size_t>(il) * l21 + imu];
					const double *__restrict w = w3j.data() + run.w_off;
					const cdouble *__restrict av = a + run.im1_begin;
					const cdouble *__restrict b = v2_ptr + run.im2_begin;
					if (v2_is_conj_of_v1) [[likely]]
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
					pre[imu] = sr;
					pim[imu] = si;
				}
				for (int imu = 0; imu < l21; ++imu)
				{
					double preal = 0.0;
					const c2r_entry *__restrict row = &c2r_nz[static_cast<size_t>(imu) * 2];
					for (int k = 0; k < c2r_cnt[imu]; ++k)
						preal += row[k].re * pre[row[k].j] - row[k].im * pim[row[k].j];
					p[offset + i + static_cast<size_t>(imu) * nfps] = preal * normfact[iat];
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
