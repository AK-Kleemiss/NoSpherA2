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

// BE AWARE, THAT V2 IS ALREADY ASSUMED TO BE CONJUGATED!!!!!
void equicomb(int natoms, int nrad1, int nrad2,
              const cvec4 &v1,
              const cvec4 &v2,
              const vec &w3j,
              const ivec2 &llvec, const int &lam,
              const cvec2 &c2r, const int &featsize,
              const int &nfps, const std::vector<int64_t> &vfps,
              vec &p)
{
    const int l21 = 2 * lam + 1;
    const int llmax = (int)llvec[0].size();

    // Initialize p with zeros
    std::memset(p.data(), 0, natoms * l21 * nfps * sizeof(double));
    const vec f_vec(featsize, 0.0);

    // Declare variables at the beginning
    int iat, n1, n2, il, imu, im1, im2, i, j, ifeat, l1, l2, mu, m2;
    double inner, normfact, preal;
    ProgressBar pb(natoms, 60, "#", " ", "Calculating descriptors for l = " + toString(lam));
#pragma omp parallel 
    {
        vec ptemp(l21 * featsize, 0.0);
        vec pcmplx_real(l21);
        vec pcmplx_imag(l21);
        const double *wigner_ptr = NULL;
        int limit_l1 = 0;
#pragma omp for private(iat, n1, n2, il, imu, im1, im2, i, j, ifeat, l1, l2, mu, m2, inner, normfact, preal) schedule(dynamic, 1)
        for (iat = 0; iat < natoms; ++iat)
        {
            inner = 0.0;
            ifeat = 0;
            for (n1 = 0; n1 < nrad1; ++n1)
            {
                for (n2 = 0; n2 < nrad2; ++n2)
                {
                    wigner_ptr = w3j.data();
                    for (il = 0; il < llmax; ++il)
                    {
                        l1 = llvec[0][il];
                        l2 = llvec[1][il];
                        limit_l1 = 2 * l1 + 1;

                        const cdouble *v1_ptr = v1[iat][n1][l1].data();
                        const cdouble *v2_ptr = v2[iat][n2][l2].data();

                        for (imu = 0; imu < l21; imu++)
                        {
                            mu = imu - lam + l1;
                            double acc_real = 0.0;
                            double acc_imag = 0.0;

                            for (im1 = 0; im1 < limit_l1; ++im1)
                            {
                                m2 = im1 - mu;
                                if (abs(m2) <= l2)
                                {
                                    im2 = m2 + l2;
                                    const double wigner_val = *wigner_ptr;
                                    const cdouble &v1_val = v1_ptr[im1];
                                    const cdouble &v2_val = v2_ptr[im2];

                                    const double& v1_r = v1_val.real();
                                    const double& v1_i = v1_val.imag();
                                    const double& v2_r = v2_val.real();
                                    const double& v2_i = v2_val.imag();

                                    acc_real += wigner_val * (v1_r * v2_r - v1_i * v2_i);
                                    acc_imag += wigner_val * (v1_r * v2_i + v1_i * v2_r);
                                    wigner_ptr++;
                                }
                            }
                            pcmplx_real[imu] = acc_real;
                            pcmplx_imag[imu] = acc_imag;
                        }
                        //recycling this variable
                        limit_l1 = l21 * ifeat;
                        for (i = 0; i < l21; ++i)
                        {
                            preal = 0.0;
                            const cdouble *__restrict cvec_ptr = c2r[i].data();
                            const double *__restrict pvec_real_ptr = pcmplx_real.data();
                            const double *__restrict pvec_imag_ptr = pcmplx_imag.data();
#pragma ivdep
                            for (j = 0; j < l21; ++j)
                            {
                                const cdouble &c2r_ih = cvec_ptr[j];
                                preal += c2r_ih.real() * pvec_real_ptr[j] - c2r_ih.imag() * pvec_imag_ptr[j];
                            }
                            inner += preal * preal;
                            ptemp[i + limit_l1] = preal;
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
                const int feat_i = vfps[i] * l21;
                for (imu = 0; imu < l21; ++imu)
                {
                    p[off_i + (imu * nfps)] = ptemp[imu + feat_i] * normfact;
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
              vec &p)
{
    // Declare variables at the beginning
    int iat, n1, n2, il, imu, im1, im2, i, j, ifeat, iwig, l1, l2, mu, m1, m2;
    double inner, normfact;

#pragma omp parallel for private(iat, n1, n2, il, imu, im1, im2, i, j, ifeat, iwig, l1, l2, mu, m1, m2, inner, normfact) default(none) shared(natoms, nrad1, nrad2, v1, v2, w3j, llmax, llvec, lam, c2r, p, featsize, constants::cnull)
    for (iat = 0; iat < natoms; ++iat)
    {
        const int l21 = 2 * lam + 1;
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
                                pcmplx[imu] += w3j[iwig] * v1[l1][iat][im1][n1] * v2[l2][iat][im2][n2];
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
            for (imu = 0; imu < 2 * lam + 1; ++imu)
            {
                p[iat * (2 * lam + 1) * featsize + (imu * featsize) + ifeat] = ptemp[imu][ifeat] / normfact;
            }
        }
    }
}
