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

// Guard failures have to survive: the program log is buffered and the progress
// bar seeks back over it, and an uncaught throw leaves no message at all.
static void equicomb_bail(const std::string &msg)
{
    std::ofstream f("equicomb_error.txt", std::ios::app);
    f << msg << std::endl;
    f.close();
    std::cerr << msg << std::endl;
    std::cerr.flush();
    throw std::out_of_range(msg);
}

// BE AWARE, THAT V2 IS ALREADY ASSUMED TO BE CONJUGATED!!!!!
void equicomb(int natoms, int nrad1, int nrad2,
              const cvec4 &v1,
              const cvec4 &v2,
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

    // The walks below trust sizes taken from the model file, not from the
    // structure: ifeat advances nrad1*nrad2*llmax times into a buffer sized by
    // featsize, wigner_ptr runs through w3j with no bound, and vfps indexes
    // ptemp. A structure missing a species the model knows makes those
    // disagree, so check here rather than run off the end inside the threads.
    // When the two descriptor sets are identical v2 is not stored at all and
    // conj(v1) is read instead, so every check below has to look at the array
    // actually used - not at v2, which is deliberately empty in that case.
    const cvec4 &v2_src = v2_is_conj_of_v1 ? v1 : v2;
    const size_t shells = static_cast<size_t>(nrad1) * nrad2 * llmax;
    if (shells > static_cast<size_t>(featsize))
    {
        equicomb_bail("equicomb: featsize " + std::to_string(featsize) +
            " is smaller than nrad1*nrad2*llmax " + std::to_string(shells));
    }
    if (v1.size() < static_cast<size_t>(natoms) || v2_src.size() < static_cast<size_t>(natoms))
    {
        equicomb_bail("equicomb: expansion coefficients hold fewer atoms than requested");
    }
    for (int chk = 0; chk < natoms; ++chk)
    {
        if (v1[chk].size() < static_cast<size_t>(nrad1) || v2_src[chk].size() < static_cast<size_t>(nrad2))
        {
            equicomb_bail("equicomb: atom " + std::to_string(chk) +
                " carries fewer radial channels than the model expects");
        }
    }
    // w3j is consumed once per (n1,n2) shell pair, one entry per surviving m2
    size_t w3j_needed = 0;
    for (int chk = 0; chk < llmax; ++chk)
    {
        const int cl1 = llvec[0][chk], cl2 = llvec[1][chk];
        if (cl1 < 0 || cl2 < 0)
        {
            equicomb_bail("equicomb: negative angular momentum in llvec");
        }
        if (v1[0][0].size() <= static_cast<size_t>(cl1) || v2_src[0][0].size() <= static_cast<size_t>(cl2))
        {
            equicomb_bail("equicomb: llvec asks for l beyond the expansion coefficients");
        }
        for (int cmu = 0; cmu < l21; ++cmu)
        {
            const int cm = cmu - lam + cl1;
            for (int cm1 = 0; cm1 < 2 * cl1 + 1; ++cm1)
            {
                if (abs(cm1 - cm) <= cl2) ++w3j_needed;
            }
        }
    }
    if (w3j.size() < w3j_needed)
    {
        equicomb_bail("equicomb: w3j holds " + std::to_string(w3j.size()) +
            " entries, the shell loop consumes " + std::to_string(w3j_needed));
    }
    for (int chk = 0; chk < nfps; ++chk)
    {
        if (vfps[chk] < 0 || static_cast<size_t>(vfps[chk]) >= static_cast<size_t>(featsize))
        {
            equicomb_bail("equicomb: vfps entry " + std::to_string(chk) + " is outside featsize");
        }
    }

    // Initialize p with zeros
    std::fill(p.begin(), p.begin() + static_cast<std::ptrdiff_t>(required_p), 0.0);

    // Declare variables at the beginning
    int iat, n1, n2, il, imu, im1, im2, i, j, ifeat, l1, l2, mu, m2;
    double inner, normfact, preal;
    // Which (im1, im2) pairs contribute is decided by |im1 - mu| <= l2, and that
    // depends only on (il, imu, lam) - never on n1 or n2. The test was being
    // re-evaluated inside the n1 x n2 loop, which runs nrad1*nrad2 times per atom
    // per lambda, and the iterations that fail it do no work.
    //
    // The survivors are more than a subset. |im1 - mu| <= l2 selects a CONTIGUOUS
    // run of im1, and im2 = im1 - mu + l2 advances in lockstep with it, so a group
    // needs no index list at all: a first index, a length and an offset into w3j
    // describe it completely. The inner loop then reads three unit-stride streams
    // instead of gathering two indices out of a term list for every term.
    //
    // w3j is consumed in exactly (il, imu, im1) order, which is why the weights of
    // one group are a contiguous slice of it and no copy of them is needed here.
    //
    // Same values, same operations, same summation order: the result is unchanged.
    struct w3j_run { int im1_begin, im2_begin, count, w_off; };
    std::vector<w3j_run> runs(static_cast<size_t>(llmax) * l21, w3j_run{0, 0, 0, 0});
    size_t total_terms = 0;
    {
        int w_idx = 0;
        for (int til = 0; til < llmax; ++til)
        {
            const int tl1 = llvec[0][til], tl2 = llvec[1][til];
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

    // The complex-to-real matrix is a mirror-pair transform: row i couples only
    // column i and column l21-1-i, so every row holds exactly two nonzeros (the
    // middle row holds one). The transform below was walking all l21 columns of
    // every row, which for lam = 5 multiplies by 119 exact zeros out of 121.
    //
    // Adding a term that is exactly zero does not change a finite sum, so dropping
    // those terms - while keeping the survivors in ascending column order, which is
    // the order the dense loop added them in - leaves the result bit for bit the
    // same. The structure is read off c2r rather than assumed, and anything that
    // does not fit two-per-row falls back to the dense walk, so a future transform
    // cannot silently lose terms here.
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
                  << ", " << total_terms << " wigner terms" << std::endl;
    }

    ProgressBar pb(natoms, 60, "#", " ", "Calculating descriptors for l = " + toString(lam));
#pragma omp parallel 
    {
        vec ptemp(l21 * featsize, 0.0);
        vec pcmplx_real(l21);
        vec pcmplx_imag(l21);
        // w * v1 depends on n1 but not on n2, and the n2 loop below runs nrad2
        // times, so that product was being formed nrad2 times over. Build it once
        // per (atom, n1) instead. Real and imaginary parts live in separate arrays
        // so the inner loop reads plain double streams rather than picking fields
        // out of a complex.
        //
        // The conjugation sign is deliberately NOT folded in here. It looks like it
        // could be, since multiplying by -1 is exact, but conj(v2) puts that sign on
        // v1_i*v2_i in the real accumulator and on v1_r*v2_i in the imaginary one,
        // so one signed copy of w*v1 cannot serve both. Folding it anyway moved the
        // form factors by 5% of the largest one - not a rounding effect, a wrong
        // answer. It is applied by choosing between two forms of the inner loop,
        // which costs one branch per (il, imu) rather than a multiply per term.
        vec wv1_re(total_terms, 0.0);
        vec wv1_im(total_terms, 0.0);
        const double *wigner_ptr = NULL;
        int limit_l1 = 0;
#pragma omp for private(iat, n1, n2, il, imu, im1, im2, i, j, ifeat, l1, l2, mu, m2, inner, normfact, preal) schedule(dynamic, 1)
        for (iat = 0; iat < natoms; ++iat)
        {
            inner = 0.0;
            ifeat = 0;
            for (n1 = 0; n1 < nrad1; ++n1)
            {
                for (int fl = 0; fl < llmax; ++fl)
                {
                    const cdouble *__restrict v1_fill = v1[iat][n1][llvec[0][fl]].data();
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
                        l1 = llvec[0][il];
                        l2 = llvec[1][il];

                        // v2 is conj(v1) when the two descriptor sets are the same
                        const cdouble *v2_ptr =
                            v2_src[iat][n2][l2].data();

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

            normfact = 1.0 / sqrt(inner);
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

void equicomb(int natoms, int nrad1, int nrad2,
              cvec4 &v1,
              cvec4 &v2,
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

    // Declare variables at the beginning
    int iat, n1, n2, il, imu, im1, im2, i, j, ifeat, iwig, l1, l2, mu, m1, m2;
    double inner, normfact;

#pragma omp parallel for private(iat, n1, n2, il, imu, im1, im2, i, j, ifeat, iwig, l1, l2, mu, m1, m2, inner, normfact) default(none) shared(natoms, nrad1, nrad2, v1, v2, w3j, llmax, llvec, lam, l21, c2r, p, featsize, constants::cnull)
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
                                // v2 is conj(v1) elementwise when the two descriptor
                                // sets match, so the same value is available from v1
                                // whatever the index order happens to be here
                                pcmplx[imu] += w3j[iwig] * v1[l1][iat][im1][n1] *
                                    (v2_is_conj_of_v1 ? std::conj(v1[l2][iat][im2][n2])
                                                      : v2[l2][iat][im2][n2]);
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
        normfact = sqrt(inner);
        for (ifeat = 0; ifeat < featsize; ++ifeat)
        {
            for (imu = 0; imu < l21; ++imu)
            {
                const size_t out_idx = static_cast<size_t>(iat) * static_cast<size_t>(l21) * static_cast<size_t>(featsize)
                    + static_cast<size_t>(imu) * static_cast<size_t>(featsize)
                    + static_cast<size_t>(ifeat);
                p[out_idx] = ptemp[imu][ifeat] / normfact;
            }
        }
    }
}
