#include "SALTED_equicomb.h"

#if has_RAS == 1
#include "cblas.h"
#endif

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
    p.assign(natoms * l21 * nfps, 0.0);
    const vec f_vec(featsize, 0.0);

    // Declare variables at the beginning
    int iat, n1, n2, il, imu, im1, im2, i, j, ifeat, iwig, l1, l2, mu, m1, m2;
    double inner, normfact;
    const cdouble null(0.0, 0.0);
    ProgressBar pb(natoms, 60, "#", " ", "Calculating descriptors for l = " + toString(lam));
#pragma omp parallel for private(iat, n1, n2, il, imu, im1, im2, i, j, ifeat, iwig, l1, l2, mu, m1, m2, inner, normfact) default(none) shared(pb, natoms, nrad1, nrad2, v1, v2, w3j, llmax, llvec, lam, c2r, nfps, vfps, p, featsize, l21, null, f_vec, std::cout)
    for (iat = 0; iat < natoms; ++iat)
    {
        vec2 ptemp(l21, f_vec);
        cvec pcmplx(l21, null);
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

                    //cvec2 *v1_ptr = (cvec2 *)&v1[l1][iat];
                    //cvec2 *v2_ptr = (cvec2 *)&v2[l2][iat];

                    cvec *v1_ptr = (cvec *)&v1[iat][n1][l1];
                    cvec *v2_ptr = (cvec *)&v2[iat][n2][l2];

                    fill(pcmplx.begin(), pcmplx.end(), null);

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
                                //pcmplx[imu] += w3j[iwig] * (*v1_ptr)[im1][n1] * (*v2_ptr)[im2][n2];
                                pcmplx[imu] += w3j[iwig] * (*v1_ptr)[im1] * (*v2_ptr)[im2];
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
        int offset = iat * l21 * nfps;
        for (int n = 0; n < nfps; ++n)
        {
            for (imu = 0; imu < l21; ++imu)
            {
                p[offset + (imu * nfps)] = ptemp[imu][vfps[n]] / normfact;
            }
            offset++;
        }
        pb.update(std::cout);
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
    // Initialize p with zeros
    p.assign(natoms * (2 * lam + 1) * featsize, 0.0);

    // Declare variables at the beginning
    int iat, n1, n2, il, imu, im1, im2, i, j, ifeat, iwig, l1, l2, mu, m1, m2;
    double inner, normfact;
    const cdouble null(0.0, 0.0);

#pragma omp parallel for private(iat, n1, n2, il, imu, im1, im2, i, j, ifeat, iwig, l1, l2, mu, m1, m2, inner, normfact) default(none) shared(natoms, nrad1, nrad2, v1, v2, w3j, llmax, llvec, lam, c2r, p, featsize, null)
    for (iat = 0; iat < natoms; ++iat)
    {
        const int l21 = 2 * lam + 1;
        vec2 ptemp(l21, vec(featsize, 0.0));
        cvec pcmplx(l21, null);
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

                    fill(pcmplx.begin(), pcmplx.end(), null);

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
