#pragma once

namespace constants
{
#include <limits>
#include <complex>
#include <cmath>
#include <map>

    double constexpr sqrtNewtonRaphson(double x, double curr, double prev)
    {
        return curr == prev
            ? curr
            : sqrtNewtonRaphson(x, 0.5 * (curr + x / curr), curr);
    }

    /*
     * Constexpr version of the square root
     * Return value:
     *   - For a finite and non-negative value of "x", returns an approximation for the square root of "x"
     *   - Otherwise, returns NaN
     * Taken from https://stackoverflow.com/questions/8622256/in-c11-is-sqrt-defined-as-constexpr
     */
    double constexpr sqrt(double x)
    {
        return x >= 0 && x < std::numeric_limits<double>::infinity()
            ? sqrtNewtonRaphson(x, x, 0)
            : std::numeric_limits<double>::quiet_NaN();
    }
    // Constants for later use
    constexpr double SQRT2 = sqrt(2.0);
    constexpr double SQRT3 = sqrt(3.0);
    constexpr double SQRT5 = sqrt(5.0);
    constexpr double INV_SQRT2 = 1.0 / SQRT2;
    constexpr double INV_SQRT3 = 1.0 / SQRT3;
    constexpr double INV_SQRT5 = 1.0 / SQRT5;
    constexpr int hardness = 3;
    constexpr double cutoff = 1.0e-20;
    constexpr double PI = 3.14159265358979323846;
    constexpr double INV_PI = 1.0 / PI;
    constexpr double PI_2 = PI / 2.0;
    constexpr double TWO_PI = 2 * PI;
    constexpr double FOUR_PI = 4 * PI;
    constexpr double C0 = SQRT2 * FOUR_PI;
    const double sqr_pi = sqrt(PI);
    constexpr double PI2 = PI * PI;
    constexpr double PI3 = PI2 * PI;
    constexpr double PI_180 = PI / 180.0;
    const double TG32 = tgamma(3.0 / 2.0);
    constexpr double ED_fact = 0.023934;
    constexpr int max_LT = 33;
    constexpr int MAG = 5810;
    //                       3,     5     7,    9,    11,   13,   15,   17
    //                      19,    21
    constexpr int lebedev_table[33] = { 6, 14, 26, 38, 50, 74, 86, 110,
                                       146, 170, 194, 230, 266, 302, 350, 434,
                                       590, 770, 974, 1202, 1454, 1730, 2030, 2354,
                                       2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810 };
    constexpr long long int ft[21]{ 1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200, 1307674368000, 20922789888000, 355687428096000, 6402373705728000, 121645100408832000, 2432902008176640000 };
    constexpr double alpha_coef = 0.1616204596739954813316614;
    constexpr double c_13 = 1.0 / 3.0;
    constexpr double c_43 = 4.0 / 3.0;
    constexpr double c_38 = 3.0 / 8.0;
    constexpr double c_16 = 1.0 / 6.0;
    constexpr double c_23 = 2.0 / 3.0;
    constexpr double c_53 = 5.0 / 3.0;
    constexpr double c_1_21 = 1.0 / 21.0;
    constexpr double c_1_30 = 1.0 / 30.0;
    constexpr double c_1_35 = 1.0 / 35.0;
    constexpr double c_1_15 = 1.0 / 15.0;
    constexpr double c_3_40 = 3.0 / 40.0;
    constexpr double c_1_79 = 1.0 / 79.0;
    constexpr double c_4_105 = 4.0 / 105.0;
    constexpr double c_1_105 = 1.0 / 105.0;
    constexpr double c_9_280 = 9.0 / 280.0;
    constexpr double c_m43 = -4.0 / 3.0;
    constexpr double c_m53 = -5.0 / 3.0;
    constexpr double barnsbohr = 2.80028520539078E+7;
    constexpr double fine_struct = 7.2973525693E-3;
    constexpr double inv_fine_struct = 1 / fine_struct;
    constexpr double fine_pi = inv_fine_struct / TWO_PI / PI;
    constexpr double inv_fine_mod = inv_fine_struct / FOUR_PI;
    constexpr double keV_per_hartree = 0.027211386245988;
    constexpr double angstrom2eV = 1.23984193 * 10000;
    constexpr double angstrom2keV = 12.3984193;
    constexpr double f_to_mu = 4208.031548;
    constexpr double barns_to_electrons = 1.43110541E-8;
    constexpr double SI2Debye = 3.33564E-30;               // in C*m
    constexpr double e_A2Debye = 0.2081943;                // in e*Angstrom
    constexpr double a0 = 0.529177210903E-10;              // in m
    constexpr double h = 6.62607015E-34 / 1.602176634E-19; // in eV*s
    constexpr double Ryd_ener = 13.6056923;                // in eV
    constexpr double alpha = 0.0072973525693;              // Sommerfeld fine structure constant
    constexpr double el_mass = 9.1093837015E-31;           // in kg
    constexpr double el_charge = 1.602176634E-19;          // in C
    constexpr double speed_of_light = 2.99792458E8;        // m/s
    constexpr double null = 0.0;
    constexpr std::complex<double> cnull = std::complex<double>(0.0, 0.0);
    constexpr std::complex<double> c_i = std::complex<double>(0.0, 1.0);

    const double ctelf = 10 * pow(2, -2.0 / 3.0) * pow(3, c_m53) * pow(PI, -c_43);
    constexpr double c_1_4p = sqrt(1.0 / (FOUR_PI));
    constexpr double c_3_4p = sqrt(3.0 / (FOUR_PI));
    constexpr double c_5_16p = sqrt(5.0 / (16.0 * PI));
    constexpr double c_7_16p = sqrt(7.0 / (16.0 * PI));
    constexpr double c_9_256p = sqrt(9.0 / (256.0 * PI));
    constexpr double c_11_256p = sqrt(11.0 / (256.0 * PI));
    constexpr double c_13_1024p = sqrt(13.0 / (1024.0 * PI));
    constexpr double c_15_4p = sqrt(15.0 / (FOUR_PI));
    constexpr double c_15_16p = sqrt(15.0 / (16.0 * PI));
    constexpr double c_21_32p = sqrt(21.0 / (32.0 * PI));
    constexpr double c_35_32p = sqrt(35.0 / (32.0 * PI));
    constexpr double c_45_16p = sqrt(45.0 / (16.0 * PI));
    constexpr double c_45_32p = sqrt(45.0 / (32.0 * PI));
    constexpr double c_45_64p = sqrt(45.0 / (64.0 * PI));
    constexpr double c_105_4p = sqrt(105.0 / (FOUR_PI));
    constexpr double c_105_16p = sqrt(105.0 / (16.0 * PI));
    constexpr double c_165_256p = sqrt(165.0 / (256.0 * PI));
    constexpr double c_273_256p = sqrt(273.0 / (256.0 * PI));
    constexpr double c_315_16p = sqrt(315.0 / (16.0 * PI));
    constexpr double c_315_32p = sqrt(315.0 / (32.0 * PI));
    constexpr double c_315_256p = sqrt(315.0 / (256.0 * PI));
    constexpr double c_385_512p = sqrt(385.0 / (512.0 * PI));
    constexpr double c_693_2048p = sqrt(693.0 / (2048.0 * PI));
    constexpr double c_1155_64p = sqrt(1155.0 / (64.0 * PI));
    constexpr double c_3465_256p = sqrt(3465.0 / (256.0 * PI));

    constexpr size_t sod = sizeof(double);
    constexpr size_t soi = sizeof(int);
    constexpr size_t soc = sizeof(char);
    constexpr size_t sob = sizeof(bool);
    constexpr size_t socd = sizeof(std::complex<double>);
    constexpr size_t soli = sizeof(long int);

    constexpr long long int ft_fun(const int& nr)
    {
        if (nr >= 0 && nr <= 20)
            return ft[nr];
        else if (nr < 0)
            return 0;
        else
            return ft_fun(nr - 1) * nr;
    }

    constexpr const double bohr2ang(const double& inp)
    {
        return inp * a0 * 1E10;
    }

    constexpr const double bohr2ang_p(const double& inp, const int& p)
    {
        if (p == 0)
            return 1.0;
        else if (p == 1)
            return bohr2ang(inp);
        else
            return bohr2ang_p(bohr2ang(inp), p - 1);
    }

    constexpr const double ang2bohr(const double& inp)
    {
        return inp / a0 * 1E-10;
    }
    constexpr const double ang2bohr_p(const double& inp, const int& p)
    {
        if (p == 0)
            return 1.0;
        else if (p == 1)
            return ang2bohr(inp);
        else
            return ang2bohr_p(ang2bohr(inp), p - 1);
    }

    constexpr const double cubic_ang2bohr(const double& inp)
    {
        return inp / ((a0 * a0 * a0) * 1E30);
    }

    constexpr const double cubic_bohr2ang(const double& inp)
    {
        return inp * (a0 * a0 * a0) * 1E30;
    }

    //------------------general functions for easy use of terminal input--------------------
    constexpr double bragg_angstrom[114]{
        0.00, // DUMMY LINE
        0.35, 0.35,
        1.45, 1.05, 0.85, 0.70, 0.65, 0.60, 0.50, 0.45,
        1.80, 1.50, 1.25, 1.10, 1.00, 1.00, 1.00, 1.00,
        2.20, 1.80, 1.60, 1.40, 1.35, 1.40, 1.40, 1.40, 1.35, 1.35, 1.35, 1.35, 1.30, 1.25, 1.15, 1.15, 1.15, 1.10,
        2.35, 2.00, 1.80, 1.55, 1.45, 1.45, 1.35, 1.30, 1.35, 1.40, 1.60, 1.55, 1.55, 1.45, 1.45, 1.40, 1.40, 1.40,
        2.60, 2.15, 1.95, 1.85, 1.85, 1.85, 1.85, 1.85, 1.85, 1.80, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.55, 1.45, 1.35, 1.30, 1.30, 1.35, 1.35, 1.35, 1.50, 1.90, 1.75, 1.60, 1.90, 1.50, 1.50,
        2.80, 2.35, 2.15, 2.05, 2.05, 2.05, 2.05, 2.05, 2.05, 2.00, 1.95, 1.95, 1.95, 1.95, 1.95, 1.95 };

    // Covalent Radii according to the CSD
    constexpr double covalent_radii[114]{
        0.0,
        0.23, 1.5,
        1.28, 0.96, 0.83, 0.68, 0.68, 0.68, 0.64, 1.5,
        1.66, 1.41, 1.21, 1.2, 1.05, 1.02, 0.99, 1.51,
        2.03, 1.76, 1.7, 1.6, 1.53, 1.39, 1.61, 1.52, 1.26, 1.24, 1.32, 1.22, 1.22, 1.17, 1.21, 1.22, 1.21, 1.5,
        2.2, 1.95, 1.9, 1.75, 1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.54, 1.42, 1.39, 1.39, 1.47, 1.4, 1.5,
        2.44, 2.15, 2.07, 2.04, 2.03, 2.01, 1.99, 1.98, 1.98, 1.96, 1.94, 1.92, 1.92, 1.89, 1.9, 1.87, 1.87, 1.75, 1.7, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32, 1.45, 1.46, 1.48, 1.4, 1.21, 1.5,
        2.6, 2.21, 2.15, 2.06, 2.00, 1.96, 1.9, 1.87, 1.8, 1.69, 1.54, 1.83, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5 };

    // Integer atom masses
    constexpr unsigned int integer_masses[]{
        1, 4,
        7, 9, 11, 12, 14, 16, 19, 20,
        23, 24, 27, 28, 31, 32, 35, 40,
        39, 40, 45, 48, 51, 52, 55, 56, 59, 58, 63, 64, 69, 74, 75, 80, 79, 84,
        85, 87, 88, 91, 92, 96, 98, 101, 103, 106, 108, 112, 115, 119, 122, 128, 127, 131,
        132, 137, 139, 140, 141, 144, 145, 150, 152, 157, 159, 163, 165, 167, 169, 173, 175, 178, 181, 184, 186, 190, 192, 195, 197, 201, 204, 207, 209, 209, 210, 222 };

    constexpr double real_masses[]{
        1.0079, 4.0026,
        6.941, 9.0122, 10.811, 12.011, 14.007, 15.999, 18.998, 20.18,
        22.99, 24.305, 26.986, 28.086, 30.974, 32.065, 35.453, 39.948,
        39.098, 40.078, 44.956, 47.867, 50.942, 51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.38, 69.723, 72.64, 74.922, 78.96, 79.904, 83.798,
        85.468, 87.62, 88.906, 91.224, 92.906, 95.96, 97.90, 101.07, 102.91, 106.42, 107.87, 112.41, 114.82, 118.71, 121.76, 127.6, 126.9, 131.29,
        132.91, 137.33, 139.91, 140.12, 140.91, 144.24, 144.9, 150.36, 151.96, 157.25, 158.93, 162.5, 164.93, 167.26, 168.93, 173.05, 174.97, 178.49, 180.95, 183.84, 186.21, 190.23, 192.22, 195.08, 196.97, 200.59, 204.38, 207.2, 208.98, 208.9, 209.9, 222.0 };

    const std::map<std::string, std::string> SPACE_GROUPS_HM = { {"1", "P1"}, {"2", "P-1"}, {"3:b", "P121"}, {"3:c", "P112"}, {"3:a", "P211"}, {"4:b", "P1211"}, {"4:c", "P1121"}, {"4:a", "P2111"}, {"5:b1", "C121"}, {"5:b2", "A121"}, {"5:b3", "I121"}, {"5:c1", "A112"}, {"5:c2", "B112"}, {"5:c3", "I112"}, {"5:a1", "B211"}, {"5:a2", "C211"}, {"5:a3", "I211"}, {"6:b", "P1m1"}, {"6:c", "P11m"}, {"6:a", "Pm11"}, {"7:b1", "P1c1"}, {"7:b2", "P1n1"}, {"7:b3", "P1a1"}, {"7:c1", "P11a"}, {"7:c2", "P11n"}, {"7:c3", "P11b"}, {"7:a1", "Pb11"}, {"7:a2", "Pn11"}, {"7:a3", "Pc11"}, {"8:b1", "C1m1"}, {"8:b2", "A1m1"}, {"8:b3", "I1m1"}, {"8:c1", "A11m"}, {"8:c2", "B11m"}, {"8:c3", "I11m"}, {"8:a1", "Bm11"}, {"8:a2", "Cm11"}, {"8:a3", "Im11"}, {"9:b1", "C1c1"}, {"9:b2", "A1n1"}, {"9:b3", "I1a1"}, {"9:-b1", "A1a1"}, {"9:-b2", "C1n1"}, {"9:-b3", "I1c1"}, {"9:c1", "A11a"}, {"9:c2", "B11n"}, {"9:c3", "I11b"}, {"9:-c1", "B11b"}, {"9:-c2", "A11n"}, {"9:-c3", "I11a"}, {"9:a1", "Bb11"}, {"9:a2", "Cn11"}, {"9:a3", "Ic11"}, {"9:-a1", "Cc11"}, {"9:-a2", "Bn11"}, {"9:-a3", "Ib11"}, {"10:b", "P12/m1"}, {"10:c", "P112/m"}, {"10:a", "P2/m11"}, {"11:b", "P121/m1"}, {"11:c", "P1121/m"}, {"11:a", "P21/m11"}, {"12:b1", "C12/m1"}, {"12:b2", "A12/m1"}, {"12:b3", "I12/m1"}, {"12:c1", "A112/m"}, {"12:c2", "B112/m"}, {"12:c3", "I112/m"}, {"12:a1", "B2/m11"}, {"12:a2", "C2/m11"}, {"12:a3", "I2/m11"}, {"13:b1", "P12/c1"}, {"13:b2", "P12/n1"}, {"13:b3", "P12/a1"}, {"13:c1", "P112/a"}, {"13:c2", "P112/n"}, {"13:c3", "P112/b"}, {"13:a1", "P2/b11"}, {"13:a2", "P2/n11"}, {"13:a3", "P2/c11"}, {"14:b1", "P121/c1"}, {"14:b2", "P121/n1"}, {"14:b3", "P121/a1"}, {"14:c1", "P1121/a"}, {"14:c2", "P1121/n"}, {"14:c3", "P1121/b"}, {"14:a1", "P21/b11"}, {"14:a2", "P21/n11"}, {"14:a3", "P21/c11"}, {"15:b1", "C12/c1"}, {"15:b2", "A12/n1"}, {"15:b3", "I12/a1"}, {"15:-b1", "A12/a1"}, {"15:-b2", "C12/n1"}, {"15:-b3", "I12/c1"}, {"15:c1", "A112/a"}, {"15:c2", "B112/n"}, {"15:c3", "I112/b"}, {"15:-c1", "B112/b"}, {"15:-c2", "A112/n"}, {"15:-c3", "I112/a"}, {"15:a1", "B2/b11"}, {"15:a2", "C2/n11"}, {"15:a3", "I2/c11"}, {"15:-a1", "C2/c11"}, {"15:-a2", "B2/n11"}, {"15:-a3", "I2/b11"}, {"16", "P222"}, {"17", "P2221"}, {"17:cab", "P2122"}, {"17:bca", "P2212"}, {"18", "P21212"}, {"18:cab", "P22121"}, {"18:bca", "P21221"}, {"19", "P212121"}, {"20", "C2221"}, {"20:cab", "A2122"}, {"20:bca", "B2212"}, {"21", "C222"}, {"21:cab", "A222"}, {"21:bca", "B222"}, {"22", "F222"}, {"23", "I222"}, {"24", "I212121"}, {"25", "Pmm2"}, {"25:cab", "P2mm"}, {"25:bca", "Pm2m"}, {"26", "Pmc21"}, {"26:ba-", "Pcm21"}, {"26:cab", "P21ma"}, {"26:-cb", "P21am"}, {"26:bca", "Pb21m"}, {"26:a-c", "Pm21b"}, {"27", "Pcc2"}, {"27:cab", "P2aa"}, {"27:bca", "Pb2b"}, {"28", "Pma2"}, {"28:ba-", "Pbm2"}, {"28:cab", "P2mb"}, {"28:-cb", "P2cm"}, {"28:bca", "Pc2m"}, {"28:a-c", "Pm2a"}, {"29", "Pca21"}, {"29:ba-", "Pbc21"}, {"29:cab", "P21ab"}, {"29:-cb", "P21ca"}, {"29:bca", "Pc21b"}, {"29:a-c", "Pb21a"}, {"30", "Pnc2"}, {"30:ba-", "Pcn2"}, {"30:cab", "P2na"}, {"30:-cb", "P2an"}, {"30:bca", "Pb2n"}, {"30:a-c", "Pn2b"}, {"31", "Pmn21"}, {"31:ba-", "Pnm21"}, {"31:cab", "P21mn"}, {"31:-cb", "P21nm"}, {"31:bca", "Pn21m"}, {"31:a-c", "Pm21n"}, {"32", "Pba2"}, {"32:cab", "P2cb"}, {"32:bca", "Pc2a"}, {"33", "Pna21"}, {"33:ba-", "Pbn21"}, {"33:cab", "P21nb"}, {"33:-cb", "P21cn"}, {"33:bca", "Pc21n"}, {"33:a-c", "Pn21a"}, {"34", "Pnn2"}, {"34:cab", "P2nn"}, {"34:bca", "Pn2n"}, {"35", "Cmm2"}, {"35:cab", "A2mm"}, {"35:bca", "Bm2m"}, {"36", "Cmc21"}, {"36:ba-", "Ccm21"}, {"36:cab", "A21ma"}, {"36:-cb", "A21am"}, {"36:bca", "Bb21m"}, {"36:a-c", "Bm21b"}, {"37", "Ccc2"}, {"37:cab", "A2aa"}, {"37:bca", "Bb2b"}, {"38", "Amm2"}, {"38:ba-", "Bmm2"}, {"38:cab", "B2mm"}, {"38:-cb", "C2mm"}, {"38:bca", "Cm2m"}, {"38:a-c", "Am2m"}, {"39", "Abm2"}, {"39:ba-", "Bma2"}, {"39:cab", "B2cm"}, {"39:-cb", "C2mb"}, {"39:bca", "Cm2a"}, {"39:a-c", "Ac2m"}, {"40", "Ama2"}, {"40:ba-", "Bbm2"}, {"40:cab", "B2mb"}, {"40:-cb", "C2cm"}, {"40:bca", "Cc2m"}, {"40:a-c", "Am2a"}, {"41", "Aba2"}, {"41:ba-", "Bba2"}, {"41:cab", "B2cb"}, {"41:-cb", "C2cb"}, {"41:bca", "Cc2a"}, {"41:a-c", "Ac2a"}, {"42", "Fmm2"}, {"42:cab", "F2mm"}, {"42:bca", "Fm2m"}, {"43", "Fdd2"}, {"43:cab", "F2dd"}, {"43:bca", "Fd2d"}, {"44", "Imm2"}, {"44:cab", "I2mm"}, {"44:bca", "Im2m"}, {"45", "Iba2"}, {"45:cab", "I2cb"}, {"45:bca", "Ic2a"}, {"46", "Ima2"}, {"46:ba-", "Ibm2"}, {"46:cab", "I2mb"}, {"46:-cb", "I2cm"}, {"46:bca", "Ic2m"}, {"46:a-c", "Im2a"}, {"47", "Pmmm"}, {"48:1", "Pnnn:1"}, {"48:2", "Pnnn:2"}, {"49", "Pccm"}, {"49:cab", "Pmaa"}, {"49:bca", "Pbmb"}, {"50:1", "Pban:1"}, {"50:2", "Pban:2"}, {"50:1ca", "Pncb:1"}, {"50:2ca", "Pncb:2"}, {"50:1bc", "Pcna:1"}, {"50:2bc", "Pcna:2"}, {"51", "Pmma"}, {"51:ba-", "Pmmb"}, {"51:cab", "Pbmm"}, {"51:-cb", "Pcmm"}, {"51:bca", "Pmcm"}, {"51:a-c", "Pmam"}, {"52", "Pnna"}, {"52:ba-", "Pnnb"}, {"52:cab", "Pbnn"}, {"52:-cb", "Pcnn"}, {"52:bca", "Pncn"}, {"52:a-c", "Pnan"}, {"53", "Pmna"}, {"53:ba-", "Pnmb"}, {"53:cab", "Pbmn"}, {"53:-cb", "Pcnm"}, {"53:bca", "Pncm"}, {"53:a-c", "Pman"}, {"54", "Pcca"}, {"54:ba-", "Pccb"}, {"54:cab", "Pbaa"}, {"54:-cb", "Pcaa"}, {"54:bca", "Pbcb"}, {"54:a-c", "Pbab"}, {"55", "Pbam"}, {"55:cab", "Pmcb"}, {"55:bca", "Pcma"}, {"56", "Pccn"}, {"56:cab", "Pnaa"}, {"56:bca", "Pbnb"}, {"57", "Pbcm"}, {"57:ba-", "Pcam"}, {"57:cab", "Pmca"}, {"57:-cb", "Pmab"}, {"57:bca", "Pbma"}, {"57:a-c", "Pcmb"}, {"58", "Pnnm"}, {"58:cab", "Pmnn"}, {"58:bca", "Pnmn"}, {"59:1", "Pmmn:1"}, {"59:2", "Pmmn:2"}, {"59:1ca", "Pnmm:1"}, {"59:2ca", "Pnmm:2"}, {"59:1bc", "Pmnm:1"}, {"59:2bc", "Pmnm:2"}, {"60", "Pbcn"}, {"60:ba-", "Pcan"}, {"60:cab", "Pnca"}, {"60:-cb", "Pnab"}, {"60:bca", "Pbna"}, {"60:a-c", "Pcnb"}, {"61", "Pbca"}, {"61:ba-", "Pcab"}, {"62", "Pnma"}, {"62:ba-", "Pmnb"}, {"62:cab", "Pbnm"}, {"62:-cb", "Pcmn"}, {"62:bca", "Pmcn"}, {"62:a-c", "Pnam"}, {"63", "Cmcm"}, {"63:ba-", "Ccmm"}, {"63:cab", "Amma"}, {"63:-cb", "Amam"}, {"63:bca", "Bbmm"}, {"63:a-c", "Bmmb"}, {"64", "Cmca"}, {"64:ba-", "Ccmb"}, {"64:cab", "Abma"}, {"64:-cb", "Acam"}, {"64:bca", "Bbcm"}, {"64:a-c", "Bmab"}, {"65", "Cmmm"}, {"65:cab", "Ammm"}, {"65:bca", "Bmmm"}, {"66", "Cccm"}, {"66:cab", "Amaa"}, {"66:bca", "Bbmb"}, {"67", "Cmma"}, {"67:ba-", "Cmmb"}, {"67:cab", "Abmm"}, {"67:-cb", "Acmm"}, {"67:bca", "Bmcm"}, {"67:a-c", "Bmam"}, {"68:1", "Ccca:1"}, {"68:2", "Ccca:2"}, {"68:1ba", "cCccb:1"}, {"68:2ba", "cCccb:2"}, {"68:1ca", "Abaa:1"}, {"68:2ca", "Abaa:2"}, {"68:1-c", "aAcaa:1"}, {"68:2-c", "aAcaa:2"}, {"68:1bc", "Bbcb:1"}, {"68:2bc", "Bbcb:2"}, {"68:1a-", "bBbab:1"}, {"68:2a-", "bBbab:2"}, {"69", "Fmmm"}, {"70:1", "Fddd:1"}, {"70:2", "Fddd:2"}, {"71", "Immm"}, {"72", "Ibam"}, {"72:cab", "Imcb"}, {"72:bca", "Icma"}, {"73", "Ibca"}, {"73:ba-", "Icab"}, {"74", "Imma"}, {"74:ba-", "Immb"}, {"74:cab", "Ibmm"}, {"74:-cb", "Icmm"}, {"74:bca", "Imcm"}, {"74:a-c", "Imam"}, {"75", "P4"}, {"76", "P41"}, {"77", "P42"}, {"78", "P43"}, {"79", "I4"}, {"80", "I41"}, {"81", "P-4"}, {"82", "I-4"}, {"83", "P4/m"}, {"84", "P42/m"}, {"85:1", "P4/n:1"}, {"85:2", "P4/n:2"}, {"86:1", "P42/n:1"}, {"86:2", "P42/n:2"}, {"87", "I4/m"}, {"88:1", "I41/a:1"}, {"88:2", "I41/a:2"}, {"89", "P422"}, {"90", "P4212"}, {"91", "P4122"}, {"92", "P41212"}, {"93", "P4222"}, {"94", "P42212"}, {"95", "P4322"}, {"96", "P43212"}, {"97", "I422"}, {"98", "I4122"}, {"99", "P4mm"}, {"100", "P4bm"}, {"101", "P42cm"}, {"102", "P42nm"}, {"103", "P4cc"}, {"104", "P4nc"}, {"105", "P42mc"}, {"106", "P42bc"}, {"107", "I4mm"}, {"108", "I4cm"}, {"109", "I41md"}, {"110", "I41cd"}, {"111", "P-42m"}, {"112", "P-42c"}, {"113", "P-421m"}, {"114", "P-421c"}, {"115", "P-4m2"}, {"116", "P-4c2"}, {"117", "P-4b2"}, {"118", "P-4n2"}, {"119", "I-4m2"}, {"120", "I-4c2"}, {"121", "I-42m"}, {"122", "I-42d"}, {"123", "P4/mmm"}, {"124", "P4/mcc"}, {"125:1", "P4/nbm:1"}, {"125:2", "P4/nbm:2"}, {"126:1", "P4/nnc:1"}, {"126:2", "P4/nnc:2"}, {"127", "P4/mbm"}, {"128", "P4/mnc"}, {"129:1", "P4/nmm:1"}, {"129:2", "P4/nmm:2"}, {"130:1", "P4/ncc:1"}, {"130:2", "P4/ncc:2"}, {"131", "P42/mmc"}, {"132", "P42/mcm"}, {"133:1", "P42/nbc:1"}, {"133:2", "P42/nbc:2"}, {"134:1", "P42/nnm:1"}, {"134:2", "P42/nnm:2"}, {"135", "P42/mbc"}, {"136", "P42/mnm"}, {"137:1", "P42/nmc:1"}, {"137:2", "P42/nmc:2"}, {"138:1", "P42/ncm:1"}, {"138:2", "P42/ncm:2"}, {"139", "I4/mmm"}, {"140", "I4/mcm"}, {"141:1", "I41/amd:1"}, {"141:2", "I41/amd:2"}, {"142:1", "I41/acd:1"}, {"142:2", "I41/acd:2"}, {"143", "P3"}, {"144", "P31"}, {"145", "P32"}, {"146:H", "R3:H"}, {"146:R", "R3:R"}, {"147", "P-3"}, {"148:H", "R-3:H"}, {"148:R", "R-3:R"}, {"149", "P312"}, {"150", "P321"}, {"151", "P3112"}, {"152", "P3121"}, {"153", "P3212"}, {"154", "P3221"}, {"155:H", "R32:H"}, {"155:R", "R32:R"}, {"156", "P3m1"}, {"157", "P31m"}, {"158", "P3c1"}, {"159", "P31c"}, {"160:H", "R3m:H"}, {"160:R", "R3m:R"}, {"161:H", "R3c:H"}, {"161:R", "R3c:R"}, {"162", "P-31m"}, {"163", "P-31c"}, {"164", "P-3m1"}, {"165", "P-3c1"}, {"166:H", "R-3m:H"}, {"166:R", "R-3m:R"}, {"167:H", "R-3c:H"}, {"167:R", "R-3c:R"}, {"168", "P6"}, {"169", "P61"}, {"170", "P65"}, {"171", "P62"}, {"172", "P64"}, {"173", "P63"}, {"174", "P-6"}, {"175", "P6/m"}, {"176", "P63/m"}, {"177", "P622"}, {"178", "P6122"}, {"179", "P6522"}, {"180", "P6222"}, {"181", "P6422"}, {"182", "P6322"}, {"183", "P6mm"}, {"184", "P6cc"}, {"185", "P63cm"}, {"186", "P63mc"}, {"187", "P-6m2"}, {"188", "P-6c2"}, {"189", "P-62m"}, {"190", "P-62c"}, {"191", "P6/mmm"}, {"192", "P6/mcc"}, {"193", "P63/mcm"}, {"194", "P63/mmc"}, {"195", "P23"}, {"196", "F23"}, {"197", "I23"}, {"198", "P213"}, {"199", "I213"}, {"200", "Pm-3"}, {"201:1", "Pn-3:1"}, {"201:2", "Pn-3:2"}, {"202", "Fm-3"}, {"203:1", "Fd-3:1"}, {"203:2", "Fd-3:2"}, {"204", "Im-3"}, {"205", "Pa-3"}, {"206", "Ia-3"}, {"207", "P432"}, {"208", "P4232"}, {"209", "F432"}, {"210", "F4132"}, {"211", "I432"}, {"212", "P4332"}, {"213", "P4132"}, {"214", "I4132"}, {"215", "P-43m"}, {"216", "F-43m"}, {"217", "I-43m"}, {"218", "P-43n"}, {"219", "F-43c"}, {"220", "I-43d"}, {"221", "Pm-3m"}, {"222:1", "Pn-3n:1"}, {"222:2", "Pn-3n:2"}, {"223", "Pm-3n"}, {"224:1", "Pn-3m:1"}, {"224:2", "Pn-3m:2"}, {"225", "Fm-3m"}, {"226", "Fm-3c"}, {"227:1", "Fd-3m:1"}, {"227:2", "Fd-3m:2"}, {"228:1", "Fd-3c:1"}, {"228:2", "Fd-3c:2"}, {"229", "Im-3m"}, {"230", "Ia-3d"} };

    constexpr int ECP_electrons[] = { 0, 0, 0,
                             0, 0, 0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                             28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28,
                             46, 46, 46, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60 };

    constexpr int ECP_electrons_xTB[] = { 0, 0, 0,
                                     2, 2, 2, 2, 2, 2, 2, 2,
                                     10, 10, 10, 10, 10, 10, 10, 10,
                                     18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 28, 28, 28, 28, 28, 28, 28,
                                     36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 46, 46, 46, 46, 46, 46, 46,
                                     54, 54, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 68, 68, 68, 68, 68, 68, 68, 68, 78, 78, 78, 78, 78, 78, 78 };

    constexpr int ECP_electrons_pTB[] = { 0, 0, 0,
                                     0, 0, 2, 2, 2, 2, 2, 2,
                                     2, 2, 10, 10, 10, 10, 10, 10,
                                     10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 28, 28, 28, 28, 28, 28,
                                     28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 46, 46, 46, 46, 46, 46,
                                     46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 60, 60, 60, 60, 60, 60, 60, 60, 60, 78, 78, 78, 78, 78, 78 };

    constexpr const char* Labels[] = { "DM", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr" };
    constexpr const char* atnr2letter(const int& nr)
    {
        if (nr == 0)
        {
            // Exception for Q peaks in residual maps
            return "Q";
        }
        if (nr > 103 || nr < 0)
        {
            if (nr == 119)
            {
                // Exception for Q in ECPs from ORCA
                return "Q";
            }
            std::cout << "Only yet implemented from H-Lr, ask Florian for improvements or give a reasonable number between 1-103!" << std::endl;
            return ("PROBLEM");
        }
        else
            return Labels[nr];
    };

    constexpr int get_Z_from_label(const char* tmp)
    {
        return
            (tmp[0] == 'H' && tmp[1] == '\0') ? 0 :
            (tmp[0] == 'D' && tmp[1] == '\0') ? 0 :
            (tmp[0] == 'T' && tmp[1] == '\0') ? 0 :
            (tmp[0] == 'H' && tmp[1] == 'e' && tmp[2] == '\0') ? 1 :
            (tmp[0] == 'L' && tmp[1] == 'i' && tmp[2] == '\0') ? 2 :
            (tmp[0] == 'B' && tmp[1] == 'e' && tmp[2] == '\0') ? 3 :
            (tmp[0] == 'B' && tmp[1] == '\0') ? 4 :
            (tmp[0] == 'C' && tmp[1] == '\0') ? 5 :
            (tmp[0] == 'N' && tmp[1] == '\0') ? 6 :
            (tmp[0] == 'O' && tmp[1] == '\0') ? 7 :
            (tmp[0] == 'F' && tmp[1] == '\0') ? 8 :
            (tmp[0] == 'N' && tmp[1] == 'e' && tmp[2] == '\0') ? 9 :
            (tmp[0] == 'N' && tmp[1] == 'a' && tmp[2] == '\0') ? 10 :
            (tmp[0] == 'M' && tmp[1] == 'g' && tmp[2] == '\0') ? 11 :
            (tmp[0] == 'A' && tmp[1] == 'l' && tmp[2] == '\0') ? 12 :
            (tmp[0] == 'S' && tmp[1] == 'i' && tmp[2] == '\0') ? 13 :
            (tmp[0] == 'P' && tmp[1] == '\0') ? 14 :
            (tmp[0] == 'S' && tmp[1] == '\0') ? 15 :
            (tmp[0] == 'C' && tmp[1] == 'l' && tmp[2] == '\0') ? 16 :
            (tmp[0] == 'A' && tmp[1] == 'r' && tmp[2] == '\0') ? 17 :
            (tmp[0] == 'K' && tmp[1] == '\0') ? 18 :
            (tmp[0] == 'C' && tmp[1] == 'a' && tmp[2] == '\0') ? 19 :
            (tmp[0] == 'S' && tmp[1] == 'c' && tmp[2] == '\0') ? 20 :
            (tmp[0] == 'T' && tmp[1] == 'i' && tmp[2] == '\0') ? 21 :
            (tmp[0] == 'V' && tmp[1] == '\0') ? 22 :
            (tmp[0] == 'C' && tmp[1] == 'r' && tmp[2] == '\0') ? 23 :
            (tmp[0] == 'M' && tmp[1] == 'n' && tmp[2] == '\0') ? 24 :
            (tmp[0] == 'F' && tmp[1] == 'e' && tmp[2] == '\0') ? 25 :
            (tmp[0] == 'C' && tmp[1] == 'o' && tmp[2] == '\0') ? 26 :
            (tmp[0] == 'N' && tmp[1] == 'i' && tmp[2] == '\0') ? 27 :
            (tmp[0] == 'C' && tmp[1] == 'u' && tmp[2] == '\0') ? 28 :
            (tmp[0] == 'Z' && tmp[1] == 'n' && tmp[2] == '\0') ? 29 :
            (tmp[0] == 'G' && tmp[1] == 'a' && tmp[2] == '\0') ? 30 :
            (tmp[0] == 'G' && tmp[1] == 'e' && tmp[2] == '\0') ? 31 :
            (tmp[0] == 'A' && tmp[1] == 's' && tmp[2] == '\0') ? 32 :
            (tmp[0] == 'S' && tmp[1] == 'e' && tmp[2] == '\0') ? 33 :
            (tmp[0] == 'B' && tmp[1] == 'r' && tmp[2] == '\0') ? 34 :
            (tmp[0] == 'K' && tmp[1] == 'r' && tmp[2] == '\0') ? 35 :
            (tmp[0] == 'R' && tmp[1] == 'b' && tmp[2] == '\0') ? 36 :
            (tmp[0] == 'S' && tmp[1] == 'r' && tmp[2] == '\0') ? 37 :
            (tmp[0] == 'Y' && tmp[1] == '\0') ? 38 :
            (tmp[0] == 'Z' && tmp[1] == 'r' && tmp[2] == '\0') ? 39 :
            (tmp[0] == 'N' && tmp[1] == 'b' && tmp[2] == '\0') ? 40 :
            (tmp[0] == 'M' && tmp[1] == 'o' && tmp[2] == '\0') ? 41 :
            (tmp[0] == 'T' && tmp[1] == 'c' && tmp[2] == '\0') ? 42 :
            (tmp[0] == 'R' && tmp[1] == 'u' && tmp[2] == '\0') ? 43 :
            (tmp[0] == 'R' && tmp[1] == 'h' && tmp[2] == '\0') ? 44 :
            (tmp[0] == 'P' && tmp[1] == 'd' && tmp[2] == '\0') ? 45 :
            (tmp[0] == 'A' && tmp[1] == 'g' && tmp[2] == '\0') ? 46 :
            (tmp[0] == 'C' && tmp[1] == 'd' && tmp[2] == '\0') ? 47 :
            (tmp[0] == 'I' && tmp[1] == 'n' && tmp[2] == '\0') ? 48 :
            (tmp[0] == 'S' && tmp[1] == 'n' && tmp[2] == '\0') ? 49 :
            (tmp[0] == 'S' && tmp[1] == 'b' && tmp[2] == '\0') ? 50 :
            (tmp[0] == 'T' && tmp[1] == 'e' && tmp[2] == '\0') ? 51 :
            (tmp[0] == 'I' && tmp[1] == '\0') ? 52 :
            (tmp[0] == 'X' && tmp[1] == 'e' && tmp[2] == '\0') ? 53 :
            (tmp[0] == 'C' && tmp[1] == 's' && tmp[2] == '\0') ? 54 :
            (tmp[0] == 'B' && tmp[1] == 'a' && tmp[2] == '\0') ? 55 :
            (tmp[0] == 'L' && tmp[1] == 'a' && tmp[2] == '\0') ? 56 :
            (tmp[0] == 'C' && tmp[1] == 'e' && tmp[2] == '\0') ? 57 :
            (tmp[0] == 'P' && tmp[1] == 'r' && tmp[2] == '\0') ? 58 :
            (tmp[0] == 'N' && tmp[1] == 'd' && tmp[2] == '\0') ? 59 :
            (tmp[0] == 'P' && tmp[1] == 'm' && tmp[2] == '\0') ? 60 :
            (tmp[0] == 'S' && tmp[1] == 'm' && tmp[2] == '\0') ? 61 :
            (tmp[0] == 'E' && tmp[1] == 'u' && tmp[2] == '\0') ? 62 :
            (tmp[0] == 'G' && tmp[1] == 'd' && tmp[2] == '\0') ? 63 :
            (tmp[0] == 'T' && tmp[1] == 'b' && tmp[2] == '\0') ? 64 :
            (tmp[0] == 'D' && tmp[1] == 'y' && tmp[2] == '\0') ? 65 :
            (tmp[0] == 'H' && tmp[1] == 'o' && tmp[2] == '\0') ? 66 :
            (tmp[0] == 'E' && tmp[1] == 'r' && tmp[2] == '\0') ? 67 :
            (tmp[0] == 'T' && tmp[1] == 'm' && tmp[2] == '\0') ? 68 :
            (tmp[0] == 'Y' && tmp[1] == 'b' && tmp[2] == '\0') ? 69 :
            (tmp[0] == 'L' && tmp[1] == 'u' && tmp[2] == '\0') ? 70 :
            (tmp[0] == 'H' && tmp[1] == 'f' && tmp[2] == '\0') ? 71 :
            (tmp[0] == 'T' && tmp[1] == 'a' && tmp[2] == '\0') ? 72 :
            (tmp[0] == 'W' && tmp[1] == '\0') ? 73 :
            (tmp[0] == 'R' && tmp[1] == 'e' && tmp[2] == '\0') ? 74 :
            (tmp[0] == 'O' && tmp[1] == 's' && tmp[2] == '\0') ? 75 :
            (tmp[0] == 'I' && tmp[1] == 'r' && tmp[2] == '\0') ? 76 :
            (tmp[0] == 'P' && tmp[1] == 't' && tmp[2] == '\0') ? 77 :
            (tmp[0] == 'A' && tmp[1] == 'u' && tmp[2] == '\0') ? 78 :
            (tmp[0] == 'H' && tmp[1] == 'g' && tmp[2] == '\0') ? 79 :
            (tmp[0] == 'T' && tmp[1] == 'l' && tmp[2] == '\0') ? 80 :
            (tmp[0] == 'P' && tmp[1] == 'b' && tmp[2] == '\0') ? 81 :
            (tmp[0] == 'B' && tmp[1] == 'i' && tmp[2] == '\0') ? 82 :
            (tmp[0] == 'P' && tmp[1] == 'o' && tmp[2] == '\0') ? 83 :
            (tmp[0] == 'A' && tmp[1] == 't' && tmp[2] == '\0') ? 85 :
            (tmp[0] == 'R' && tmp[1] == 'n' && tmp[2] == '\0') ? 86 :
            (tmp[0] == 'F' && tmp[1] == 'r' && tmp[2] == '\0') ? 87 :
            (tmp[0] == 'R' && tmp[1] == 'a' && tmp[2] == '\0') ? 88 :
            (tmp[0] == 'A' && tmp[1] == 'c' && tmp[2] == '\0') ? 89 :
            (tmp[0] == 'T' && tmp[1] == 'h' && tmp[2] == '\0') ? 90 :
            (tmp[0] == 'P' && tmp[1] == 'a' && tmp[2] == '\0') ? 91 :
            (tmp[0] == 'U' && tmp[1] == '\0') ? 92 :
            (tmp[0] == 'N' && tmp[1] == 'p' && tmp[2] == '\0') ? 93 :
            (tmp[0] == 'P' && tmp[1] == 'u' && tmp[2] == '\0') ? 94 :
            (tmp[0] == 'A' && tmp[1] == 'm' && tmp[2] == '\0') ? 95 :
            (tmp[0] == 'C' && tmp[1] == 'm' && tmp[2] == '\0') ? 96 :
            (tmp[0] == 'B' && tmp[1] == 'k' && tmp[2] == '\0') ? 97 :
            (tmp[0] == 'C' && tmp[1] == 'f' && tmp[2] == '\0') ? 98 :
            (tmp[0] == 'E' && tmp[1] == 's' && tmp[2] == '\0') ? 99 :
            (tmp[0] == 'F' && tmp[1] == 'm' && tmp[2] == '\0') ? 100 :
            (tmp[0] == 'M' && tmp[1] == 'd' && tmp[2] == '\0') ? 101 :
            (tmp[0] == 'N' && tmp[1] == 'o' && tmp[2] == '\0') ? 102 :
            (tmp[0] == 'L' && tmp[1] == 'r' && tmp[2] == '\0') ? 103 :
            (tmp[0] == 'R' && tmp[1] == 'f' && tmp[2] == '\0') ? 104 :
            (tmp[0] == 'D' && tmp[1] == 'b' && tmp[2] == '\0') ? 105 :
            (tmp[0] == 'S' && tmp[1] == 'g' && tmp[2] == '\0') ? 106 :
            (tmp[0] == 'B' && tmp[1] == 'h' && tmp[2] == '\0') ? 107 :
            (tmp[0] == 'H' && tmp[1] == 's' && tmp[2] == '\0') ? 108 :
            (tmp[0] == 'M' && tmp[1] == 't' && tmp[2] == '\0') ? 109 :
            (tmp[0] == 'D' && tmp[1] == 's' && tmp[2] == '\0') ? 110 :
            (tmp[0] == 'R' && tmp[1] == 'g' && tmp[2] == '\0') ? 111 :
            -1;
    };

    const int type_vector[168]{
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    0, 0, 1,
    2, 0, 0,
    0, 2, 0,
    0, 0, 2,
    1, 1, 0,
    1, 0, 1,
    0, 1, 1,
    3, 0, 0,
    0, 3, 0,
    0, 0, 3,
    2, 1, 0,
    2, 0, 1,
    0, 2, 1,
    1, 2, 0,
    1, 0, 2,
    0, 1, 2,
    1, 1, 1,
    0, 0, 4,
    0, 1, 3,
    0, 2, 2,
    0, 3, 1,
    0, 4, 0,
    1, 0, 3,
    1, 1, 2,
    1, 2, 1,
    1, 3, 0,
    2, 0, 2,
    2, 1, 1,
    2, 2, 0,
    3, 0, 1,
    3, 1, 0,
    4, 0, 0,
    0, 0, 5,
    0, 1, 4,
    0, 2, 3,
    0, 3, 2,
    0, 4, 1,
    0, 5, 0,
    1, 0, 4,
    1, 1, 3,
    1, 2, 2,
    1, 3, 1,
    1, 4, 0,
    2, 0, 3,
    2, 1, 2,
    2, 2, 1,
    2, 3, 0,
    3, 0, 2,
    3, 1, 1,
    3, 2, 0,
    4, 0, 1,
    4, 1, 0,
    5, 0, 0 };

    constexpr void type2vector(int index, int* vector)
    {
        if (index < 1 || index > 56)
        {
            vector[0] = -1;
            vector[1] = -1;
            vector[2] = -1;
            return;
        }
        const int temp = index - 1;
        vector[0] = constants::type_vector[temp * 3];
        vector[1] = constants::type_vector[temp * 3 + 1];
        vector[2] = constants::type_vector[temp * 3 + 2];
    };

    constexpr double normgauss(const int& type, const double& exp)
    {
        int t[3]{ 0,0,0 };
        if (type > 0)
        {
            constants::type2vector(type, t);
            err_checkf(t[0] != -1, "Problem with type2vector!", std::cout);
            err_checkf(t[1] != -1, "Problem with type2vector!", std::cout);
            err_checkf(t[2] != -1, "Problem with type2vector!", std::cout);
        }
        else
            t[0] = t[1] = t[2] = 0;
        long long int temp = constants::ft[t[0]] * constants::ft[t[1]] * constants::ft[t[2]];
        int t1 = 2 * t[0], t2 = 2 * t[1], t3 = 2 * t[2];
        long long int temp2 = constants::ft[t1] * constants::ft[t2] * constants::ft[t3];
        double temp1 = 2 * exp / constants::PI;
        temp1 = temp1 * temp1 * temp1;
        temp1 = constants::sqrt(constants::sqrt(temp1));
        const int exponent = t[0] + t[1] + t[2];
        double temp_e = 1.0;
        for (int i = 0; i < exponent; i++)
            temp_e *= 8 * exp;
        return temp1 * constants::sqrt(temp_e * temp / temp2);
    };

    constexpr double spherical_harmonic(const int& l, const int& m, const double* d)
    {
        /*Here d[0] = x
                     d[1] = y
                     d[2] = z
                     d[3] = r^2 IGNORED
                     d[4] = r   IGNORED
                     */
                     // Will need extension for up to l=8
                     // calc spherical harmonic
        double SH = 0, x = d[0], y = d[1], z = d[2];
        switch (l)
        {
        case 0: // S
            SH = constants::c_1_4p;
            break;
        case 1:
            switch (m)
            {
            case 0: // P 0 Z
                SH = constants::c_3_4p * z;
                break;
            case 1: // P 1 X
                SH = constants::c_3_4p * x;
                break;
            case -1: // P -1 Y
                SH = constants::c_3_4p * y;
                break;
            default:
                err_not_impl_f("Wrong spherical harmonic called!", std::cout);
            }
            break;
        case 2:
            switch (m)
            {
            case 0: // D 0 Z2
                SH = constants::c_5_16p * (3 * z*z - 1.0);
                break;
            case 1: // D 1 XZ
                SH = constants::c_15_4p * x * z;
                break;
            case -1: // D -1 YZ
                SH = constants::c_15_4p * y * z;
                break;
            case 2: // D 2 X2-Y2
                SH = constants::c_15_16p * (x*x - y*y);
                break;
            case -2: // D -2 XY
                SH = constants::c_15_4p * y * x;
                break;
            default:
                err_not_impl_f("Wrong spherical harmonic called!", std::cout);
            }
            break;
        case 3:
            switch (m)
            {
            case 0: // F 0 Z3
                SH = constants::c_7_16p * (5 * z*z*z - 3 * z);
                break;
            case 1: // F 1 XZZ
                SH = constants::c_21_32p * x * (5 * z*z - 1.0);
                break;
            case -1: // F -1 YZZ
                SH = constants::c_21_32p * y * (5 * z*z - 1.0);
                break;
            case 2: // F 2 Z(X2-Y2)
                SH = constants::c_105_16p * ((x*x - y*y) * z);
                break;
            case -2: // F -2 XYZ
                SH = constants::c_105_4p * x * y * z;
                break;
            case 3: // F 3 X(X^2-3Y^2)
                SH = constants::c_35_32p * x * (x*x - 3 * y*y);
                break;
            case -3: // F -3 Y(3X^2-Y^2)
                SH = constants::c_35_32p * y * (3 * x*x - y*y);
                break;
            default:
                err_not_impl_f("Wrong spherical harmonic called!", std::cout);
            }
            break;
        case 4:
            switch (m)
            {
            case 0: // G 0 Z^4
                SH = constants::c_9_256p * (35 * z*z*z*z - 30 * z*z + 3.0);
                break;
            case 1: // G 1 X(7Z^3-3ZR^2)
                SH = constants::c_45_32p * x * (7 * z*z*z - 3 * z);
                break;
            case -1: // G -1 Y(7Z^2-3ZR^2)
                SH = constants::c_45_32p * y * (7 * z*z*z - 3 * z);
                break;
            case 2: // G 2
                SH = constants::c_45_64p * (x*x - y*y) * (7 * z*z - 1.0);
                break;
            case -2: // G -2
                SH = constants::c_45_16p * x * y * (7 * z*z - 1.0);
                break;
            case 3: // G 3 XZ(X^2-3Y^2)
                SH = constants::c_315_32p * x * (x*x - 3 * y*y) * z;
                break;
            case -3: // G -3 XZ(3X^2-Y^2)
                SH = constants::c_315_32p * y * (3 * x*x - y*y) * z;
                break;
            case 4: // G 4 X^2(X^-3Y^2)-Y^2(3X^2-Y^2)
                SH = constants::c_315_256p * ((x*x * (x*x - 3 * y*y)) -
                    (y*y * (3 * x*x - y*y)));
                break;
            case -4: // G -4 XY(X^2-Y^2)
                SH = constants::c_315_16p * x * y * (x*x - y*y);
                break;
            default:
                err_not_impl_f("Wrong spherical harmonic called!", std::cout);
            }
            break;
        case 5:
            switch (m)
            {
            case 0: // H Z^5
                SH = constants::c_11_256p * (63 * z*z*z*z*z - 70 * z*z*z + 15 * z);
                break;
            case 1:
                SH = constants::c_165_256p * x * (21 * z*z*z*z - 14 * z*z + 1.0);
                break;
            case -1:
                SH = constants::c_165_256p * y * (21 * z*z*z*z - 14 * z*z + 1.0);
                break;
            case 2:
                SH = constants::c_1155_64p * (x*x - y*y) * (3 * z*z*z - z);
                break;
            case -2:
                SH = constants::c_1155_64p * 2 * x * y * (3 * z*z*z - z);
                break;
            case 3:
                SH = constants::c_385_512p * x * (x*x - 3 * y*y) * (9 * z*z - 1.0);
                break;
            case -3:
                SH = constants::c_385_512p * y * (3 * x*x - y*y) * (9 * z*z - 1.0);
                break;
            case 4:
                SH = constants::c_3465_256p * (x*x*x*x - 6 * x * x * y * y + y*y*y*y) * z;
                break;
            case -4:
                SH = -constants::c_3465_256p * (4 * x * y*y*y - 4 * x*x*x * y) * z;
                break;
            case 5:
                SH = constants::c_693_2048p * (2 * x*x*x*x*x - 20 * x*x*x * y*y + 10 * x * y*y*y*y);
                break;
            case -5:
                SH = constants::c_693_2048p * (2 * y*y*y*y*y - 20 * x*x * y*y*y + 10 * y * x*x*x*x);
                break;
            default:
                err_not_impl_f("Wrong spherical harmonic called!", std::cout);
            }
            break;
        case 6:
            switch (m)
            {
            case 0: // I Z^6
                SH = constants::c_13_1024p * (231 * z*z*z*z*z*z - 315 * z*z*z*z + 105 * z*z - 5);
                break;
            case 1:
                SH = constants::c_273_256p * x * (21 * z*z*z*z - 14 * z*z + 1.0);
                break;
            case -1:
                SH = constants::c_165_256p * 2 * y * (21 * z*z*z*z - 14 * z*z + 1.0);
                break;
            case 2:
                SH = constants::c_1155_64p * (x*x - y*y) * (3 * z*z*z - z);
                break;
            case -2:
                SH = constants::c_1155_64p * 2 * x * y * (3 * z*z*z - z);
                break;
            case 3:
                SH = constants::c_385_512p * x * (x*x - 3 * y*y) * (9 * z*z - 1.0);
                break;
            case -3:
                SH = constants::c_385_512p * y * (3 * x*x - y*y) * (9 * z*z - 1.0);
                break;
            case 4:
                SH = constants::c_3465_256p * (x*x*x*x - 6 * x * x * y * y + y*y*y*y) * z;
                break;
            case -4:
                SH = -constants::c_3465_256p * (4 * x * y*y*y - 4 * x*x*x * y) * z;
                break;
            case 5:
                SH = constants::c_693_2048p * (2 * x*x*x*x*x - 20 * x*x*x * y*y + 10 * x * y*y*y*y);
                break;
            case -5:
                SH = constants::c_693_2048p * (2 * y*y*y*y*y - 20 * x*x * y*y*y + 10 * y * x*x*x*x);
                break;
            default:
                err_not_impl_f("Wrong spherical harmonic called!", std::cout);
            }
            break;
        default:
            err_not_impl_f("Higehr than l=6 not done for spherical harmonic!", std::cout);
        }
        return SH;
    }
}