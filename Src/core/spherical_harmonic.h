#pragma once
#include <cmath>

//The real spherical harmonics of one l contracted with their 2l+1 coefficients, x, y, z
//normalised; the l <= 8 switch of constants::spherical_harmonic(l, d, coefs) in constants.cpp,
//which forwards here. Templated on the number type so aux_density.h differentiates it with
//dual, and free of project includes so aux_density_gpu.cu compiles it under nvcc and hipcc:
//those see the host-device attribute, every other compiler plain inline functions.
#if defined(__CUDACC__) || defined(__HIPCC__)
#define AUX_HD __host__ __device__
#else
#define AUX_HD
#endif

namespace constants
{
    template <class T>
    AUX_HD inline T spherical_harmonic(const int l, const T x, const T y, const T z, const double* coefs)
    {
        constexpr double PI = 3.1415926535897932384626433832795028;
        switch (l)
        {
        case 0:
            return sqrt(1.0 / (4 * PI)) * coefs[0];
        case 1:
            return sqrt(3.0 / (4 * PI)) * (coefs[0] * y + coefs[1] * z + coefs[2] * x);
        case 2:
            return sqrt(15.0 / (4 * PI)) * (y * x * coefs[0] + y * z * coefs[1] + x * z * coefs[3]) + sqrt(5.0 / (16.0 * PI)) * (3 * z * z - 1.0) * coefs[2] + sqrt(15.0 / (16.0 * PI)) * (x * x - y * y) * coefs[4];
        case 3:
        {
            const T y2 = y * y, x2 = x * x, z2 = z * z;
            return sqrt(35.0 / (32.0 * PI)) * (y * (3 * x2 - y2) * coefs[0] + x * (x2 - 3 * y2) * coefs[6]) +
                sqrt(105.0 / (4 * PI)) * x * y * z * coefs[1] +
                sqrt(21.0 / (32.0 * PI)) * (y * (5 * z2 - 1.0) * coefs[2] + x * (5 * z2 - 1.0) * coefs[4]) +
                sqrt(7.0 / (16.0 * PI)) * (5 * z2 * z - 3 * z) * coefs[3] +
                sqrt(105.0 / (16.0 * PI)) * ((x2 - y2) * z) * coefs[5];
        }
        case 4:
        {
            const T x2 = x * x, y2 = y * y, z2 = z * z;
            return sqrt(315.0 / (16.0 * PI)) * x * y * (x2 - y2) * coefs[0] +
                sqrt(315.0 / (32.0 * PI)) * (y * (3 * x2 - y2) * z * coefs[1] + x * (x2 - 3 * y2) * z * coefs[7]) +
                sqrt(45.0 / (16.0 * PI)) * x * y * (7 * z2 - 1.0) * coefs[2] +
                sqrt(45.0 / (32.0 * PI)) * (y * (7 * z2 * z - 3 * z) * coefs[3] + x * (7 * z2 * z - 3 * z) * coefs[5]) +
                sqrt(9.0 / (256.0 * PI)) * (35 * z2 * z2 - 30 * z2 + 3.0) * coefs[4] +
                sqrt(45.0 / (64.0 * PI)) * (x2 - y2) * (7 * z2 - 1.0) * coefs[6] +
                sqrt(315.0 / (256.0 * PI)) * ((x2 * (x2 - 3 * y2)) - (y2 * (3 * x2 - y2))) * coefs[8];
        }
        case 5:
        {
            const T x2 = x * x, y2 = y * y, z2 = z * z;
            return sqrt(693.0 / (2048.0 * PI)) * (2 * y2 * y2 * y - 20 * x2 * y2 * y + 10 * y * x2 * x2) * coefs[0] +
                -sqrt(3465.0 / (256.0 * PI)) * z * ((4 * x * y2 * y - 4 * x2 * x * y) * coefs[1] - (x2 * x2 - 6 * x2 * y2 + y2 * y2) * coefs[9]) +
                sqrt(385.0 / (512.0 * PI)) * (9 * z2 - 1.0) * (y * (3 * x2 - y2) * coefs[2] + x * (x2 - 3 * y2) * coefs[8]) +
                sqrt(1155.0 / (64.0 * PI)) * (3 * z2 * z - z) * (2 * x * y * coefs[3] + (x2 - y2) * coefs[7]) +
                sqrt(165.0 / (256.0 * PI)) * (21 * z2 * z2 - 14 * z2 + 1.0) * (y * coefs[4] + x * coefs[6]) +
                sqrt(11.0 / (256.0 * PI)) * (63 * z2 * z2 * z - 70 * z2 * z + 15 * z) * coefs[5] +
                sqrt(693.0 / (2048.0 * PI)) * (2 * x2 * x2 * x - 20 * x2 * x * y2 + 10 * x * y2 * y2) * coefs[10];
        }
        case 6:
        {
            const T x2 = x * x, y2 = y * y, z2 = z * z;
            const T x4 = x2 * x2, y4 = y2 * y2, z4 = z2 * z2;
            return (1.0 / 32.0) * sqrt(13.0 / PI) * (z2 * (z2 * (231.0 * z2 - 315.0) + 105.0) - 5.0) * coefs[6] +
                (1.0 / 16.0) * sqrt(273.0 / (PI)) * (33 * z4 - 30 * z2 + 5) * (x * z * coefs[7] + y * z * coefs[5]) +
                (1.0 / 32.0) * sqrt(1365.0 / (2 * PI)) * (33 * z4 - 18 * z2 + 1) * ((x2 - y2) * coefs[8] + 2 * x * y * coefs[4]) +
                (1.0 / 16.0) * sqrt(1365.0 / (2 * PI)) * z * (11 * z2 - 3) * (x * (x2 - 3 * y2) * coefs[9] + y * (3 * x2 - y2) * coefs[3]) +
                (3.0 / 32.0) * sqrt(91.0 / (PI)) * (11 * z2 - 1) * ((x2 * (x2 - 6 * y2) + y4) * coefs[10] + x * y * 4 * (x2 - y2) * coefs[2]) +
                (3.0 / 16.0) * sqrt(1001.0 / (2 * PI)) * z * ((x4 * x - 10 * x2 * x * y2 + 5 * x * y4) * coefs[11] + (5 * x4 * y - 10 * x2 * y2 * y + y4 * y) * coefs[1]) +
                (1.0 / 32.0) * sqrt(3003.0 / (2 * PI)) * (-y4 * y2 + 15 * y4 * x2 - 15 * y2 * x4 + x4 * x2) * coefs[12] +
                (1.0 / 32.0) * sqrt(3003.0 / (2 * PI)) * (6 * x4 * x * y - 20 * x2 * x * y2 * y + 6 * x * y4 * y) * coefs[0];
        }
        case 7:
        {
            const T x2 = x * x, y2 = y * y, z2 = z * z;
            const T x4 = x2 * x2, y4 = y2 * y2;
            return (1.0 / 32.0) * sqrt(15.0 / PI) * z * (z2 * (z2 * (429.0 * z2 - 693.0) + 315.0) - 35.0) * coefs[7] +
                (1.0 / 64.0) * sqrt(105.0 / PI) * (z2 * (z2 * (429.0 * z2 - 495.0) + 135.0) - 5.0) * (x * coefs[8] + y * coefs[6]) +
                (3.0 / 32.0) * sqrt(35.0 / (2 * PI)) * z * (z2 * (143.0 * z2 - 110.0) + 15.0) * ((x2 - y2) * coefs[9] + 2 * x * y * coefs[5]) +
                (3.0 / 64.0) * sqrt(35.0 / PI) * (z2 * (143.0 * z2 - 66.0) + 3.0) * (x * (x2 - 3 * y2) * coefs[10] + y * (3 * x2 - y2) * coefs[4]) +
                (3.0 / 32.0) * sqrt(385.0 / PI) * z * (13.0 * z2 - 3.0) * ((x2 * (x2 - 6 * y2) + y4) * coefs[11] + x * y * 4 * (x2 - y2) * coefs[3]) +
                (3.0 / 64.0) * sqrt(385.0 / PI) * (13.0 * z2 - 1.0) * (x * (x4 + 5 * y2 * (-2 * x2 + y2)) * coefs[12] + y * (5 * x2 * (x2 - 2 * y2) + y4) * coefs[2]) +
                (3.0 / 32.0) * sqrt(5005.0 / (2 * PI)) * z * ((y2 * (y2 * (-y2 + 15 * x2) - 15 * x4) + x4 * x2) * coefs[13] + x * y * (x2 * (6 * x2 - 20 * y2) + 6 * y4) * coefs[1]) +
                (3.0 / 64.0) * sqrt(715.0 / PI) * x * (x4 * x2 + 7 * y2 * (-3 * x4 + y2 * (5 * x2 - y2))) * coefs[14] +
                (3.0 / 64.0) * sqrt(715.0 / PI) * y * (x2 * (x2 * (7 * x2 - 35 * y2) + 21 * y4) - y4 * y2) * coefs[0];
        }
        case 8:
        {
            const T x2 = x * x, y2 = y * y, z2 = z * z;
            const T x4 = x2 * x2, y4 = y2 * y2, z4 = z2 * z2;
            return (1.0 / 256.0) * sqrt(17.0 / PI) * (z2 * (z2 * (z2 * (6435.0 * z2 - 12012.0) + 6930.0) - 1260.0) + 35.0) * coefs[8] +
                (3.0 / 64.0) * sqrt(17.0 / PI) * z * (715.0 * z4 * z2 - 1001.0 * z4 + 385.0 * z2 - 35.0) * (x * coefs[9] + y * coefs[7]) +
                (3.0 / 64.0) * sqrt(595.0 / (2 * PI)) * (143 * z4 * z2 - 143 * z4 + 33 * z2 - 1) * ((x2 - y2) * coefs[10] + 2 * x * y * coefs[6]) +
                (1.0 / 64.0) * sqrt(19635.0 / PI) * z * (39 * z4 - 26 * z2 + 3) * ((x2 * x - 3 * x * y2) * coefs[11] + (3 * x2 * y - y2 * y) * coefs[5]) +
                (3.0 / 128.0) * sqrt(1309.0 / PI) * (65 * z4 - 26 * z2 + 1) * ((x2 * (x2 - 6 * y2) + y4) * coefs[12] + x * y * 4 * (x2 - y2) * coefs[4]) +
                (3.0 / 64.0) * sqrt(17017.0 / PI) * z * (5 * z2 - 1) * (x * (x4 + 5 * y2 * (-2 * x2 + y2)) * coefs[13] + y * (5 * x2 * (x2 - 2 * y2) + y4) * coefs[3]) +
                (1.0 / 64.0) * sqrt(7293.0 / (2 * PI)) * (15 * z2 - 1) * ((y2 * (y2 * (-y2 + 15 * x2) - 15 * x4) + x4 * x2) * coefs[14] + x * y * (x2 * (6 * x2 - 20 * y2) + 6 * y4) * coefs[2]) +
                (3.0 / 64.0) * sqrt(12155.0 / PI) * z * (x * (x4 * x2 + 7 * y2 * (-3 * x4 + y2 * (5 * x2 - y2))) * coefs[15] + y * (x2 * (x2 * (7 * x2 - 35 * y2) + 21 * y4) - y4 * y2) * coefs[1]) +
                (3.0 / 256.0) * sqrt(12155.0 / PI) * (y2 * (y2 * (y4 - 28 * x2 * y2 + 70 * x4) - 28 * x4 * x2) + x4 * x4) * coefs[16] +
                (3.0 / 256.0) * sqrt(12155.0 / PI) * y * x * (y2 * (y2 * (-8 * y2 + 56 * x2) - 56 * x4) + 8 * x4 * x2) * coefs[0];
        }
        }
        return T(0.0);
    }
}
