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

    int iat, n1, n2, il, imu, i, j, ifeat, l2;
    double inner, normfact, preal;
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
    // rather than assumed; anything not two-per-row falls back to the dense walk.
    struct c2r_entry { int j; double re, im; };
    std::vector<c2r_entry> c2r_nz(static_cast<size_t>(l21) * 2, c2r_entry{0, 0.0, 0.0});
    std::vector<int> c2r_cnt(l21, 0);
    bool c2r_is_sparse = (c2r.size() >= static_cast<size_t>(l21));
    for (int i2 = 0; i2 < l21 && c2r_is_sparse; ++i2)
    {
        if (c2r[i2].size() < static_cast<size_t>(l21)) { c2r_is_sparse = false; break; }
        int cnt = 0;
        for (int j2 = 0; j2 < l21; ++j2)
        {
            const cdouble &e = c2r[i2][j2];
            if (e.real() == 0.0 && e.imag() == 0.0) continue;
            if (cnt == 2) { c2r_is_sparse = false; break; }
            c2r_nz[static_cast<size_t>(i2) * 2 + cnt] = { j2, e.real(), e.imag() };
            ++cnt;
        }
        c2r_cnt[i2] = cnt;
    }
    if (ProgressBar::report_counts)
    {
        std::cout << "[equicomb] lam " << lam << ": c2r "
                  << (c2r_is_sparse ? "sparse (2 per row)" : "dense fallback")
                  << ", " << total_terms << " wigner terms"
                  << ", natoms " << natoms << ", nrad1 " << nrad1 << ", nrad2 " << nrad2
                  << ", llmax " << llmax << ", l21 " << l21
                  << ", featsize " << featsize << ", nfps " << nfps << std::endl;
    }

    int empty_environments = 0;

#ifdef NOSPHERA2_USE_GPU
    //The device reproduces this walk exactly; it falls through to the CPU loop below if
    //no device is present, the transform is not the two-per-row form, or it will not fit.
    if (g_equicomb_use_gpu && c2r_is_sparse)
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
        //Reverse of vfps: which output slot, if any, each shell triple feeds
        ivec sel(static_cast<size_t>(featsize), -1);
        for (int i2 = 0; i2 < nfps; ++i2) sel[static_cast<size_t>(vfps[i2])] = i2;
        const SALTEDDescriptors &v2_gpu = v2_is_conj_of_v1 ? v1 : v2;
        salted_gpu_problem q;
        q.natoms = natoms; q.nrad1 = nrad1; q.nrad2 = nrad2; q.llmax = llmax;
        q.lam = lam; q.l21 = l21; q.featsize = featsize; q.nfps = nfps;
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
        q.sel = sel.data(); q.p = p.data();
        const bool gpu_ok = salted_gpu_equicomb(q, &empty_environments);
        //Once per run, not once per lambda: nine identical lines say nothing extra.
        //Printed before the progress bar exists, whose carriage returns would eat it.
        static bool announced = false;
        if (!announced) {
            announced = true;
            std::cout << "GPU in use: SALTED descriptors on "
                      << (gpu_ok ? "the device (double precision)" : "the CPU - device unavailable") << std::endl;
        }
        if (ProgressBar::report_counts)
            std::cout << "[equicomb] lam " << lam << ": GPU " << (gpu_ok ? "used" : "refused, using the CPU loop") << std::endl;
        if (gpu_ok)
            return;
    }
#endif


    // Scoped: the bar rewinds to its own line when it is destroyed, so anything
    // printed after the loop but before that is silently overwritten.
    {
    ProgressBar pb(natoms, 60, "#", " ", "Calculating descriptors for l = " + toString(lam));
#pragma omp parallel 
    {
        vec ptemp(static_cast<size_t>(l21) * featsize, 0.0);
        vec pcmplx_real(l21);
        vec pcmplx_imag(l21);
        // w * v1 depends on n1 but not on n2, so build it once per (atom, n1); real
        // and imaginary parts live in separate arrays so the inner loop reads plain
        // double streams rather than picking fields out of a complex.
        // The conjugation sign is deliberately NOT folded in here: conj(v2) puts that
        // sign on v1_i*v2_i in the real accumulator and on v1_r*v2_i in the imaginary
        // one, so one signed copy of w*v1 cannot serve both. It is applied by choosing
        // between two forms of the inner loop instead.
        vec wv1_re(total_terms, 0.0);
        vec wv1_im(total_terms, 0.0);
        const double *wigner_ptr = NULL;
        int limit_l1 = 0;
#pragma omp for private(iat, n1, n2, il, imu, i, j, ifeat, l2, inner, normfact, preal) schedule(dynamic, 1)
        for (iat = 0; iat < natoms; ++iat)
        {
            inner = 0.0;
            ifeat = 0;
            for (n1 = 0; n1 < nrad1; ++n1)
            {
                for (int fl = 0; fl < llmax; ++fl)
                {
                    const cdouble *__restrict v1_fill = v1.block(iat, n1, llvec[0][fl]);
                    for (int fmu = 0; fmu < l21; ++fmu)
                    {
                        const w3j_run &fr = runs[static_cast<size_t>(fl) * l21 + fmu];
                        for (int k = 0; k < fr.count; ++k)
                        {
                            const double wk = w3j[static_cast<size_t>(fr.w_off) + k];
                            const cdouble &av = v1_fill[fr.im1_begin + k];
                            wv1_re[static_cast<size_t>(fr.w_off) + k] = wk * av.real();
                            wv1_im[static_cast<size_t>(fr.w_off) + k] = wk * av.imag();
                        }
                    }
                }
                for (n2 = 0; n2 < nrad2; ++n2)
                {
                    for (il = 0; il < llmax; ++il)
                    {
                        l2 = llvec[1][il];

                        // v2 is conj(v1) when the two descriptor sets are the same
                        const cdouble *v2_ptr = v2_src.block(iat, n2, l2);

                        for (imu = 0; imu < l21; imu++)
                        {
                            double acc_real = 0.0;
                            double acc_imag = 0.0;

                            const w3j_run &run = runs[static_cast<size_t>(il) * l21 + imu];
                            const double *__restrict ar = wv1_re.data() + run.w_off;
                            const double *__restrict ai = wv1_im.data() + run.w_off;
                            const cdouble *__restrict b = v2_ptr + run.im2_begin;
                            if (v2_is_conj_of_v1)
                            {
                                for (int k = 0; k < run.count; ++k)
                                {
                                    const double v2_r = b[k].real();
                                    const double v2_i = b[k].imag();
                                    acc_real += ar[k] * v2_r + ai[k] * v2_i;
                                    acc_imag += ai[k] * v2_r - ar[k] * v2_i;
                                }
                            }
                            else
                            {
                                for (int k = 0; k < run.count; ++k)
                                {
                                    const double v2_r = b[k].real();
                                    const double v2_i = b[k].imag();
                                    acc_real += ar[k] * v2_r - ai[k] * v2_i;
                                    acc_imag += ar[k] * v2_i + ai[k] * v2_r;
                                }
                            }
                            pcmplx_real[imu] = acc_real;
                            pcmplx_imag[imu] = acc_imag;
                        }
                        //recycling this variable
                        limit_l1 = l21 * ifeat;
                        const double *__restrict pvec_real_ptr = pcmplx_real.data();
                        const double *__restrict pvec_imag_ptr = pcmplx_imag.data();
                        if (c2r_is_sparse)
                        {
                            for (i = 0; i < l21; ++i)
                            {
                                preal = 0.0;
                                const c2r_entry *__restrict row = &c2r_nz[static_cast<size_t>(i) * 2];
                                const int nz = c2r_cnt[i];
                                for (int k = 0; k < nz; ++k)
                                {
                                    preal += row[k].re * pvec_real_ptr[row[k].j] - row[k].im * pvec_imag_ptr[row[k].j];
                                }
                                inner += preal * preal;
                                ptemp[i + limit_l1] = preal;
                            }
                        }
                        else
                        {
                            for (i = 0; i < l21; ++i)
                            {
                                preal = 0.0;
                                const cdouble *__restrict cvec_ptr = c2r[i].data();
#if defined(_MSC_VER)
#pragma loop(ivdep)
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC ivdep
#endif
                                for (j = 0; j < l21; ++j)
                                {
                                    const cdouble &c2r_ih = cvec_ptr[j];
                                    preal += c2r_ih.real() * pvec_real_ptr[j] - c2r_ih.imag() * pvec_imag_ptr[j];
                                }
                                inner += preal * preal;
                                ptemp[i + limit_l1] = preal;
                            }
                        }
                        ifeat++;
                    }
                }
            }

            // An empty environment gives an all-zero descriptor, so inner is 0 and
            // 1/sqrt(inner) is +inf, making every feature NaN. Zero is the meaningful
            // answer: the kernel contributes nothing and the atom keeps the species
            // average the model adds separately.
            if (inner > 0.0)
            {
                normfact = 1.0 / sqrt(inner);
            }
            else
            {
                normfact = 0.0;
#pragma omp atomic
                ++empty_environments;
            }
            const int offset = iat * l21 * nfps;
            for (i = 0; i < nfps; ++i)
            {
                const int off_i = i + offset;
                const int feat_i = static_cast<int>(vfps[i]) * l21;
                for (imu = 0; imu < l21; ++imu)
                {
                    const size_t out_idx = static_cast<size_t>(off_i + (imu * nfps));
                    const size_t feat_idx = static_cast<size_t>(imu + feat_i);
                    p[out_idx] = ptemp[feat_idx] * normfact;
                }
            }
            //pb.update(std::cout);
        }
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
