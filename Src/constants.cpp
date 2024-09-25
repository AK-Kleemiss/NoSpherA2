#include "constants.h"
#include "convenience.h"

namespace constants {
    const char* atnr2letter(const int& nr)
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
    }

    const double normgauss(const int& type, const double& exp)
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

    const double spherical_harmonic(const int& l, const int& m, const double* d)
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
                SH = constants::c_5_16p * (3 * z * z - 1.0);
                break;
            case 1: // D 1 XZ
                SH = constants::c_15_4p * x * z;
                break;
            case -1: // D -1 YZ
                SH = constants::c_15_4p * y * z;
                break;
            case 2: // D 2 X2-Y2
                SH = constants::c_15_16p * (x * x - y * y);
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
                SH = constants::c_7_16p * (5 * z * z * z - 3 * z);
                break;
            case 1: // F 1 XZZ
                SH = constants::c_21_32p * x * (5 * z * z - 1.0);
                break;
            case -1: // F -1 YZZ
                SH = constants::c_21_32p * y * (5 * z * z - 1.0);
                break;
            case 2: // F 2 Z(X2-Y2)
                SH = constants::c_105_16p * ((x * x - y * y) * z);
                break;
            case -2: // F -2 XYZ
                SH = constants::c_105_4p * x * y * z;
                break;
            case 3: // F 3 X(X^2-3Y^2)
                SH = constants::c_35_32p * x * (x * x - 3 * y * y);
                break;
            case -3: // F -3 Y(3X^2-Y^2)
                SH = constants::c_35_32p * y * (3 * x * x - y * y);
                break;
            default:
                err_not_impl_f("Wrong spherical harmonic called!", std::cout);
            }
            break;
        case 4:
            switch (m)
            {
            case 0: // G 0 Z^4
                SH = constants::c_9_256p * (35 * z * z * z * z - 30 * z * z + 3.0);
                break;
            case 1: // G 1 X(7Z^3-3ZR^2)
                SH = constants::c_45_32p * x * (7 * z * z * z - 3 * z);
                break;
            case -1: // G -1 Y(7Z^2-3ZR^2)
                SH = constants::c_45_32p * y * (7 * z * z * z - 3 * z);
                break;
            case 2: // G 2
                SH = constants::c_45_64p * (x * x - y * y) * (7 * z * z - 1.0);
                break;
            case -2: // G -2
                SH = constants::c_45_16p * x * y * (7 * z * z - 1.0);
                break;
            case 3: // G 3 XZ(X^2-3Y^2)
                SH = constants::c_315_32p * x * (x * x - 3 * y * y) * z;
                break;
            case -3: // G -3 XZ(3X^2-Y^2)
                SH = constants::c_315_32p * y * (3 * x * x - y * y) * z;
                break;
            case 4: // G 4 X^2(X^-3Y^2)-Y^2(3X^2-Y^2)
                SH = constants::c_315_256p * ((x * x * (x * x - 3 * y * y)) -
                    (y * y * (3 * x * x - y * y)));
                break;
            case -4: // G -4 XY(X^2-Y^2)
                SH = constants::c_315_16p * x * y * (x * x - y * y);
                break;
            default:
                err_not_impl_f("Wrong spherical harmonic called!", std::cout);
            }
            break;
        case 5:
            switch (m)
            {
            case 0: // H Z^5
                SH = constants::c_11_256p * (63 * z * z * z * z * z - 70 * z * z * z + 15 * z);
                break;
            case 1:
                SH = constants::c_165_256p * x * (21 * z * z * z * z - 14 * z * z + 1.0);
                break;
            case -1:
                SH = constants::c_165_256p * y * (21 * z * z * z * z - 14 * z * z + 1.0);
                break;
            case 2:
                SH = constants::c_1155_64p * (x * x - y * y) * (3 * z * z * z - z);
                break;
            case -2:
                SH = constants::c_1155_64p * 2 * x * y * (3 * z * z * z - z);
                break;
            case 3:
                SH = constants::c_385_512p * x * (x * x - 3 * y * y) * (9 * z * z - 1.0);
                break;
            case -3:
                SH = constants::c_385_512p * y * (3 * x * x - y * y) * (9 * z * z - 1.0);
                break;
            case 4:
                SH = constants::c_3465_256p * (x * x * x * x - 6 * x * x * y * y + y * y * y * y) * z;
                break;
            case -4:
                SH = -constants::c_3465_256p * (4 * x * y * y * y - 4 * x * x * x * y) * z;
                break;
            case 5:
                SH = constants::c_693_2048p * (2 * x * x * x * x * x - 20 * x * x * x * y * y + 10 * x * y * y * y * y);
                break;
            case -5:
                SH = constants::c_693_2048p * (2 * y * y * y * y * y - 20 * x * x * y * y * y + 10 * y * x * x * x * x);
                break;
            default:
                err_not_impl_f("Wrong spherical harmonic called!", std::cout);
            }
            break;
        case 6:
            switch (m)
            {
            case 0: // I Z^6
                SH = constants::c_13_1024p * (231 * z * z * z * z * z * z - 315 * z * z * z * z + 105 * z * z - 5);
                break;
            case 1:
                SH = constants::c_273_256p * x * (21 * z * z * z * z - 14 * z * z + 1.0);
                break;
            case -1:
                SH = constants::c_165_256p * 2 * y * (21 * z * z * z * z - 14 * z * z + 1.0);
                break;
            case 2:
                SH = constants::c_1155_64p * (x * x - y * y) * (3 * z * z * z - z);
                break;
            case -2:
                SH = constants::c_1155_64p * 2 * x * y * (3 * z * z * z - z);
                break;
            case 3:
                SH = constants::c_385_512p * x * (x * x - 3 * y * y) * (9 * z * z - 1.0);
                break;
            case -3:
                SH = constants::c_385_512p * y * (3 * x * x - y * y) * (9 * z * z - 1.0);
                break;
            case 4:
                SH = constants::c_3465_256p * (x * x * x * x - 6 * x * x * y * y + y * y * y * y) * z;
                break;
            case -4:
                SH = -constants::c_3465_256p * (4 * x * y * y * y - 4 * x * x * x * y) * z;
                break;
            case 5:
                SH = constants::c_693_2048p * (2 * x * x * x * x * x - 20 * x * x * x * y * y + 10 * x * y * y * y * y);
                break;
            case -5:
                SH = constants::c_693_2048p * (2 * y * y * y * y * y - 20 * x * x * y * y * y + 10 * y * x * x * x * x);
                break;
            case 6:
                SH = constants::c_3003_2048p * (x * x * x * x * x * x - 15 * x * x * x * x * y * y + 15 * x * x * y * y * y * y - y * y * y * y * y);  //THIS WAS GUESS BY GITHUB MAKE SURE THIS WORKS
                break;
            case -6:
                SH = constants::c_3003_2048p * (6 * x * x * x * x * y - 20 * x * x * y * y * y + 6 * y * y * y * y * y);  //THIS WAS GUESS BY GITHUB MAKE SURE THIS WORKS
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

    double associated_legendre_polynomial(const int &l, const int &m, const double &x) {
        switch (l) {
        case 0:
            return 1.0;
            break;
        case 1:
            switch (m) {
            case 0:
                return x;
                break;
            case 1:
                return sqrt(1 - x * x);
                break;
            case -1:
                return -0.5 * sqrt(1 - x * x);
                break;
            default:
                err_checkf(false, "This is impossible!", std::cout);
            }
            break;
        case 2:
            switch (m) {
            case 0:
                return 0.5 * (3 * x * x - 1);
                break;
            case 1:
                return 3.0 * x * sqrt(1 - x * x);
                break;
            case 2:
                return -3.0 * (x * x - 1);
                break;
            case -1:
                return -0.5 * x * sqrt(1 - x * x);
                break;
            case -2:
                return -0.125 * (1 - x * x);
                break;
            default:
                err_checkf(false, "This is impossible!", std::cout);
            }
            break;
        case 3:
            switch (m) {
            case 0:
                return 0.5 * x * (5 * x * x - 3);
                break;
            case 1:
                return 1.5 * (5 * x * x - 1) * sqrt(1 - x * x);
                break;
            case 2:
                return -15.0 * x * (x * x - 1);
                break;
            case 3:
                return 15.0 * pow(1 - x * x, 1.5);
                break;
            case -1:
                return -0.125 * (5 * x * x - 1) * sqrt(1 - x * x);
                break;
            case -2:
                return 0.125 * x * (x * x - 1);
                break;
            case -3:
                return -(1.0 / 48.0) * pow(1 - x * x, 1.5);
                break;
            default:
                err_checkf(false, "This is impossible!", std::cout);
            }
            break;
        case 4:
            switch (m) {
            case 0:
                return 4.375 * x * x * x * x - 3.75 * x * x + 0.375;
                break;
            case 1:
                return 2.5 * sqrt(1 - x * x) * (7 * x * x * x - 3 * x);
                break;
            case 2:
                return -7.5 * (x * x - 1) * (7 * x * x - 1);
                break;
            case 3:
                return 105.0 * x * pow(1 - x * x, 1.5);
                break;
            case 4:
                return 105.0 * pow(x * x - 1, 2.0);
                break;
            case -1:
                return -0.125 * sqrt(1 - x * x) * (7 * x * x * x - 3 * x);
                break;
            case -2:
                return (1.0 / 48.0) * (x * x - 1) * (7 * x * x - 1);
                break;
            case -3:
                return -(1.0 / 48.0) * pow(1 - x, 1.5) * x * pow(x + 1, 1.5);
                break;
            case -4:
                return -(1.0 / 384.0) * pow(x * x - 1, 2.0);
                break;
            default:
                err_checkf(false, "This is impossible!", std::cout);
            }
            break;
        case 5:
            switch (m) {
            case 0:
                return 7.875 * x * x * x * x * x - 8.75 * x * x * x + 1.875 * x;
                break;
            case 1:
                return (15.0 / 8.0) * sqrt(1 - x * x) * (21 * x * x * x * x - 14 * x * x + 1);
                break;
            case 2:
                return -(105.0 / 2.0) * (x * x - 1) * (3 * x * x * x - x);
                break;
            case 3:
                return (105.0 / 2.0) * pow(1 - x * x, 1.5) * (9 * x * x - 1);
                break;
            case 4:
                return 945.0 * x * pow(x * x - 1, 2);
                break;
            case 5:
                return 945.0 * pow(1 - x * x, 2.5);
                break;
            case -1:
                return -sqrt(1 - x * x) * (1.3125 * x * x * x * x - 0.875 * x * x + 0.0625);
                break;
            case -2:
                return (1.0 / 16.0) * (x * x - 1) * (3 * x * x * x - x);
                break;
            case -3:
                return -(1.0 / 384.0) * pow(1 - x * x, 1.5) * (9 * x * x - 1);
                break;
            case -4:
                return -(1.0 / 384.0) * x * pow(x * x - 1, 2.0);
                break;
            case -5:
                return -(1.0 / 3840.0) * pow(1 - x * x, 2.5);
                break;
            default:
                err_checkf(false, "This is impossible!", std::cout);
            }
            break;
        case 6:  //THIS IS WILDLY UNTESTED BECAUSE I DONT TRUST constants::spherical_harmonic!!!! (WORKS FOR m=0)
            switch (m) {
            case 0:
                return 14.4375 * x * x * x * x * x * x - 19.6875 * x * x * x * x + 6.5625 * x * x - 0.3125;
                break;
            case 1:
                return -2.625 * x * sqrt(1 - x * x) * (33 * x * x * x * x - 30 * x * x + 5);
                break;
            case 2:
                return -13.125 * (x * x - 1) * (33 * x * x * x * x - 18 * x * x + 1);
                break;
            case 3:
                return -157.5 * pow(1 - x * x, 1.5) * x * (11 * x * x - 3);
                break;
            case 4:
                return 472.5 * pow(x * x - 1, 2.0) * (11 * x * x - 1);
                break;
            case 5:
                return -10395.0 * x * pow(1 - x * x, 2.5);
                break;
            case 6:
                return -10395.0 * pow(x * x - 1, 3.0);
                break;
            case -1:
                return 0.0625 * x * sqrt(1 - x * x) * (33 * x * x * x * x - 30 * x * x + 5);
                break;
            case -2:
                return -0.0078125 * (x * x - 1) * (33 * x * x * x * x - 18 * x * x + 1);
                break;
            case -3:
                return -(1.0 / 384.0) * pow(1 - x * x, 1.5) * x * (11 * x * x - 3);
                break;
            case -4:
                return (1.0 / 3840.0) * pow(x * x - 1, 2.0) * (11 * x * x - 1);
                break;
            case -5:
                return (1.0 / 3840.0) * x * pow(1 - x * x, 2.5);
                break;
            case -6:
                return -(1.0 / 46080.0) * pow(x * x - 1, 3.0);
                break;
            default:
                err_checkf(false, "This is impossible!", std::cout);
            }
            break;
        default:
            err_not_impl_f("associated_legendre_polynomial l > 6 ", std::cout);
            return -1000.;
        }
        return -1000.;
    };

    //Returns the spherical coordinates of a given cartesian vector
    //the output is a vector with the (radius, theta and phi)
    vec cartesian_to_spherical(const double& x, const double& y, const double& z) {
        double r = sqrt(x * x + y * y + z * z);
        if (r == 0) {
            return { r, 0., 0. };
        }
        else {
            return { r, acos(z / r), atan2(y, x) };
        }
    }


    //Original implementation after P. Coppens DOI: 10.1107/97809553602060000759 Eq. 1.2.7.2b
    //I omitted the abs(m) in the factorial as most other sources do not include it
    double real_spherical(const int& l, const int& m, const double& theta, const double& phi) {
        double N;
        m == 0 ? N = sqrt((2 * l + 1) / constants::FOUR_PI) : N = sqrt(((2 * l + 1) / constants::TWO_PI) * double(constants::ft[l - (m)]) / double(constants::ft[l + (m)]));
        if (m >= 0) {
            return N * associated_legendre_polynomial(l, m, cos(theta)) * cos(m * phi);
        }
        else {
            return N * associated_legendre_polynomial(l, m, cos(theta)) * sin(m * phi);
        }
    }
}