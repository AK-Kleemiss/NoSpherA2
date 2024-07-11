#include "SALTED_equicomb.h"

using namespace std;

void equicomb(int natoms, int nang1, int nang2, int nrad1, int nrad2,
    const cvec4& v1,
    const cvec4& v2,
    const vec& w3j, int llmax,
    const ivec2& llvec, int lam,
    const cvec2& c2r, int featsize,
    int nfps, const vector<int64_t>& vfps,
    vec3& p)
{
    // Initialize p with zeros
    p.assign(2 * lam + 1, vec2(nfps, vec(natoms, 0.0)));

    // Declare variables at the beginning
    int iat, n1, n2, il, imu, im1, im2, i, j, ifeat, iwig, l1, l2, mu, m1, m2;
    double inner, normfact;

#pragma omp parallel for private(iat, n1, n2, il, imu, im1, im2, i, j, ifeat, iwig, l1, l2, mu, m1, m2, inner, normfact) default(none) shared(natoms, nang1, nang2, nrad1, nrad2, v1, v2, w3j, llmax, llvec, lam, c2r, nfps, vfps, p, featsize)
    for (iat = 0; iat < natoms; ++iat)
    {
        const int l21 = 2 * lam + 1;
        vec2 ptemp(l21, vec(featsize, 0.0));
        cvec pcmplx(l21, complex<double>(0.0, 0.0));
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

                    fill(pcmplx.begin(), pcmplx.end(), complex<double>(0.0, 0.0));

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
                                pcmplx[imu] += w3j[iwig] * v1[im1][l1][n1][iat] * conj(v2[im2][l2][n2][iat]);
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
        for (int n = 0; n < nfps; ++n)
        {
            for (imu = 0; imu < l21; ++imu)
            {
                p[imu][n][iat] = ptemp[imu][vfps[n]] / normfact;
            }
        }
    }
}

void equicomb(int natoms, int nang1, int nang2, int nrad1, int nrad2,
    cvec4& v1,
    cvec4& v2,
    vec& w3j, int llmax,
    ivec2& llvec, int lam,
    cvec2& c2r, int featsize,
    vec3& p)
{

    cout << "WARNING EQUICOMB IS HAS NOT BEEN TESTED, PROCEED WITH CAUTION!!!" << endl;
    // Initialize p with zeros
    p = vec3(2 * lam + 1, vec2(featsize, vec(natoms, 0.0)));

    // Parallel region
#pragma omp parallel for default(none) shared(natoms, nang1, nang2, nrad1, nrad2, v1, v2, w3j, llmax, llvec, lam, c2r, featsize, p)
    for (int iat = 0; iat < natoms; ++iat)
    {
        double inner = 0.0;
        vec2 ptemp(2 * lam + 1, vec(featsize, 0.0));
        int ifeat = 0;
        for (int n1 = 0; n1 < nrad1; ++n1)
        {
            for (int n2 = 0; n2 < nrad2; ++n2)
            {
                int iwig = 0;
                for (int il = 0; il < llmax; ++il)
                {
                    int l1 = llvec[0][il];
                    int l2 = llvec[1][il];
                    cvec pcmplx(2 * lam + 1, complex<double>(0.0, 0.0));
                    for (int imu = 0; imu < 2 * lam + 1; ++imu)
                    {
                        int mu = imu - lam;
                        for (int im1 = 0; im1 < 2 * l1 + 1; ++im1)
                        {
                            int m1 = im1 - l1;
                            int m2 = m1 - mu;
                            if (std::abs(m2) <= l2)
                            {
                                int im2 = m2 + l2;
                                pcmplx[imu] += w3j[iwig] * v1[im1][l1][n1][iat] * std::conj(v2[im2][l2][n2][iat]);
                                iwig++;
                            }
                        }
                    }
                    vec preal(2 * lam + 1);
                    // Matrix multiplication c2r * pcmplx
                    for (int i = 0; i < 2 * lam + 1; ++i)
                    {
                        cdouble sum = 0.0;
                        for (int j = 0; j < 2 * lam + 1; ++j)
                        {
                            sum += c2r[i][j] * pcmplx[j];
                        }
                        preal[i] = std::real(sum);
                        inner += std::norm(sum);
                        ptemp[i][ifeat] = preal[i];
                    }
                    ifeat++;
                }
            }
        }
        const double normfact = std::sqrt(inner);
        for (ifeat = 0; ifeat < featsize; ++ifeat)
        {
            for (int imu = 0; imu < 2 * lam + 1; ++imu)
            {
                p[imu][ifeat][iat] = ptemp[imu][ifeat] / normfact;
            }
        }
    }
}
