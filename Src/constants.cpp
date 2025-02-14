#include "pch.h"
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

    const std::complex<double> complex_spherical_harmonic(const int& l, const int& m, double& theta, double& phi){
        if (l < 6) {
			std::cout << "Nopyt Nope, spherical harmonic" << std::endl;
            exit(1);
        }
        std::complex<double> res = 0;
        switch (l) {
        case 6:
            switch (m) {
            case -6:
				res = 1.0 / 64.0 * sqrt(3003.0 / PI) * pow(sin(theta), 6) * exp(-6.0 * cone * phi);
				break;
			case -5:
				res = 1.0 / 32.0 * sqrt(1001.0 / PI) * pow(sin(theta), 5) * cos(theta) * exp(-5.0 * cone * phi);
                break;
			case -4:
				res = 3.0 / 32.0 * sqrt(91.0 / (2.0 * PI)) * pow(sin(theta), 4) * (11.0 * cos(theta) * cos(theta) - 1.0) * exp(-4.0 * cone * phi);
				break;
            case -3:
				res = 1.0 / 32.0 * sqrt(1365.0 / PI) * pow(sin(theta), 3) * cos(theta) * (11.0 * cos(theta) * cos(theta) - 3.0) * exp(-3.0 * cone * phi);
                break;
			case -2:
				res = 1.0 / 64.0 * sqrt(1365.0 / PI) * pow(sin(theta), 2) * (33.0 *  pow(cos(theta), 4.0) - 18.0 * cos(theta) * cos(theta) + 1) * exp(-2.0 * cone * phi);
                break;
			case -1:
				res = 1.0 / 16.0 * sqrt(273.0 / (2.0 * PI)) * sin(theta) * cos(theta) * (33.0 * pow(cos(theta), 4.0) - 30.0 * cos(theta) * cos(theta) + 5) * exp(-1.0 * cone * phi);
				break;
            case 0:
				res = 1.0 / 32.0 * sqrt(13.0 / PI) * (231.0 * pow(cos(theta), 6) - 315.0 * pow(cos(theta), 4) + 105.0 * pow(cos(theta), 2) - 5.0);
				break;
			case 1:
				res = -1.0 / 16.0 * sqrt(273.0 / (2.0 * PI)) * sin(theta) * cos(theta) * (33.0 * pow(cos(theta), 4.0) - 30.0 * cos(theta) * cos(theta) + 5) * exp(1.0 * cone * phi);
                break;
			case 2:
				res = 1.0 / 64.0 * sqrt(1365.0 / PI) * pow(sin(theta), 2) * (33.0 * pow(cos(theta), 4.0) - 18.0 * cos(theta) * cos(theta) + 1) * exp(2.0 * cone * phi);
                break;
			case 3:
				res = -1.0 / 32.0 * sqrt(1365.0 / PI) * pow(sin(theta), 3) * cos(theta) * (11.0 * cos(theta) * cos(theta) - 3.0) * exp(3.0 * cone * phi);
				break;
			case 4:
				res = 3.0 / 32.0 * sqrt(91.0 / (2.0 * PI)) * pow(sin(theta), 4) * (11.0 * cos(theta) * cos(theta) - 1.0) * exp(4.0 * cone * phi);
				break;
			case 5:
				res = -3 / 32.0 * sqrt(1001.0 / PI) * pow(sin(theta), 5) * cos(theta) * exp(5.0 * cone * phi);
				break;
			case 6:
				res = 1.0 / 64.0 * sqrt(3003.0 / PI) * pow(sin(theta), 6) * exp(6.0 * cone * phi);
                break;
            }
        }
        return res;
    }

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
        double leng;
        std::pair<double, double> spherical;
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
        case 6:  // Generated via Wolfram Alpha using: (SphericalHarmonicY[6,-m,theta,phi] + SphericalHarmonicY[6,m,theta,phi])  //BEWARE THE SIGN!!!!
            leng = sqrt(x * x + y * y + z * z);
            spherical = constants::norm_cartesian_to_spherical(x / leng, y / leng, z / leng);
            switch (m)
            {
            case 0: // I Z^6
                //1/32 sqrt(13/PI) (231 cos^6(theta) - 315 cos^4(theta) + 105 cos^2(theta) - 5)
				SH = (1.0 / 32.0) * sqrt(13.0 / PI) * (231.0 * pow(cos(spherical.first), 6.0) - 315.0 * pow(cos(spherical.first), 4) + 105.0 * pow(cos(spherical.first), 2.0) - 5.0);
                break;
            case 1:
				SH = (1.0 / 16.0) * sqrt(273.0 / (PI)) * sin(spherical.first) * cos(spherical.first) * (33 * pow(cos(spherical.first), 4) - 30 * pow(cos(spherical.first), 2) + 5) * cos(spherical.second);
                break;
            case -1:
                //-1/256 i sqrt(273/(2 PI)) (5 sin(2 theta) + 12 sin(4 theta) + 33 sin(6 theta)) sin(phi)
				SH = (1.0 / 256.0) * sqrt(273.0 / (4 * PI)) * (5 * sin(2 * spherical.first) + 12 * sin(4 * spherical.first) + 33 * sin(6 * spherical.first)) * sin(spherical.second);
                break;
            case 2:
				SH = (1.0 / 32.0) * sqrt(1365.0 / (2 * PI)) * pow(sin(spherical.first), 2) * (33 * pow(cos(spherical.first), 4) - 18 * pow(cos(spherical.first), 2) + 1) * cos(spherical.second * 2);
                break;
            case -2:
                // -1/256 i sqrt(1365/PI) sin^2(theta) (60 cos(2 theta) + 33 cos(4 theta) + 35) sin(2 phi)
                SH = (1.0 / 256.0) * sqrt(1365.0 / (2*PI)) * pow(sin(spherical.first), 2) * (60 * cos(spherical.first * 2) + 33 * cos(spherical.first * 4) + 35) * sin(spherical.second * 2);
                break;
            case 3:
				SH = (1.0 / 16.0) * sqrt(1365.0 / (2 * PI)) * pow(sin(spherical.first), 3) * cos(spherical.first) * (11 * pow(cos(spherical.first), 2) - 3) * cos(spherical.second * 3);
                break;
            case -3:
                //-1/64 i sqrt(1365/PI) sin^3(theta) (21 cos(theta) + 11 cos(3 theta)) sin(3 phi)
                SH = (1.0 / 64.0) * sqrt(1365.0 / (2*PI)) * pow(sin(spherical.first), 3) * (21 * cos(spherical.first) + 11 * cos(spherical.first * 3)) * sin(spherical.second * 3);
                break;
            case 4:
				SH = (3.0 / 64.0) * sqrt(91.0 / (PI)) * pow(sin(spherical.first), 4) * (11 * cos(2*spherical.first) + 9) * cos(spherical.second * 4);
                break;
            case -4:
                //-3/32 i sqrt(91/(2 PI)) sin^4(theta) (11 cos(2 theta) + 9) sin(4 phi)
                SH = (3.0 / 32.0) * sqrt(91.0 / (4 * PI)) * pow(sin(spherical.first), 4) * (11 * cos(spherical.first * 2) + 9) * sin(spherical.second * 4);
                break;
            case 5:
				SH = (3.0 / 16.0) * sqrt(1001.0 / (2 * PI)) * pow(sin(spherical.first), 5) * cos(spherical.first) * cos(spherical.second * 5);
                break;
            case -5:
                SH = (3.0 / 16.0) * sqrt(1001.0 / (2 * PI)) * pow(sin(spherical.first), 5) * cos(spherical.first) * sin(spherical.second * 5);
                break;
            case 6:
				SH = (1.0 / 32.0) * sqrt(3003.0 / (2 * PI)) * pow(sin(spherical.first), 6) * cos(spherical.second * 6);
                break;
            case -6:
                SH = (1.0 / 32.0) * sqrt(3003.0 / (2 * PI)) * pow(sin(spherical.first), 6.0) * sin(spherical.second * 6);
                break;
            default:
                err_not_impl_f("Wrong spherical harmonic called!", std::cout);
            }
            break;
		case 7:// Generated via Wolfram Alpha using: (SphericalHarmonicY[7,-m,theta,phi] + SphericalHarmonicY[7,m,theta,phi])  //BEWARE THE SIGN!!!!
            leng = sqrt(x * x + y * y + z * z);
            spherical = constants::norm_cartesian_to_spherical(x / leng, y / leng, z / leng);
			switch (m)
			{
			case 0:
				SH = (1.0 / 32.0) * sqrt(15.0 / PI) * (429.0 * pow(cos(spherical.first), 7) - 693.0 * pow(cos(spherical.first), 5) + 315.0 * pow(cos(spherical.first), 3) - 35.0 * cos(spherical.first));
                break;
			case 1:
				SH = (1.0 / 64.0) * sqrt(105.0 / PI) * sin(spherical.first) * (429.0 * pow(cos(spherical.first), 6) - 495.0 * pow(cos(spherical.first), 4) + 135.0 * pow(cos(spherical.first), 2) - 5.0) * cos(spherical.second);
				break;
			case -1:
                SH = (1.0 / 64.0) * sqrt(105.0 / PI) * sin(spherical.first) * (429.0 * pow(cos(spherical.first), 6) - 495.0 * pow(cos(spherical.first), 4) + 135.0 * pow(cos(spherical.first), 2) - 5.0) * sin(spherical.second);
				break;
			case 2:
				SH = (3.0 / 32.0) * sqrt(35.0 / (2 * PI)) * pow(sin(spherical.first), 2) * cos(spherical.first) * (143.0 * pow(cos(spherical.first), 4.0) - 110.0 * pow(cos(spherical.first), 2.0) + 15.0) * cos(spherical.second * 2);
				break;
			case -2:
                SH = (3.0 / 128.0) * sqrt(35.0 / (2 * PI)) * pow(sin(spherical.first), 2) * cos(spherical.first) * (132.0 * cos(spherical.first * 2) + 143.0 * cos(spherical.first * 4) + 109.0) * sin(spherical.second) * cos(spherical.second);
				break;
			case 3:
				SH = (3.0 / 64.0) * sqrt(35.0 / PI) * pow(sin(spherical.first), 3) * (143.0 * pow(cos(spherical.first), 4.0) - 66.0 * pow(cos(spherical.first), 2.0) + 3.0) * cos(spherical.second * 3);
				break;
			case -3:
                SH = (3.0 / 64.0) * sqrt(35.0 / PI) * pow(sin(spherical.first), 3) * (143.0 * pow(cos(spherical.first), 4.0) - 66.0 * pow(cos(spherical.first), 2.0) + 3.0) * sin(spherical.second * 3);
				break;
			case 4:
				SH = (3.0 / 32.0) * sqrt(385.0 / PI) * pow(sin(spherical.first), 4) * cos(spherical.first) * (13.0 * pow(cos(spherical.first), 2.0) - 3.0) * cos(spherical.second * 4);
				break;
			case -4:
                SH = (3.0 / 32.0) * sqrt(385.0 / PI) * pow(sin(spherical.first), 4) * cos(spherical.first) * (13.0 * pow(cos(spherical.first), 2.0) - 3.0) * sin(spherical.second * 4);
				break;
			case 5:
				SH = (3.0 / 64.0) * sqrt(385.0 / PI) * pow(sin(spherical.first), 5) * (13.0 * pow(cos(spherical.first), 2.0) - 1.0) * cos(spherical.second * 5);
				break;
			case -5:
				SH = (3.0 / 64.0) * sqrt(385.0 / PI) * pow(sin(spherical.first), 5) * (13.0 * pow(cos(spherical.first), 2.0) - 1.0) * sin(spherical.second * 5);
				break;
			case 6:
                SH = (3.0 / 32.0) * sqrt(5005.0 / (2 * PI)) * pow(sin(spherical.first), 6) * cos(spherical.first) * cos(spherical.second * 6);
				break;
			case -6:
				SH = (3.0 / 32.0) * sqrt(5005.0 / (2 * PI)) * pow(sin(spherical.first), 6) * cos(spherical.first) * sin(spherical.second * 6);
				break;
			case 7:
				SH = (3.0 / 64.0) * sqrt(715.0 / PI) * pow(sin(spherical.first), 7) * cos(spherical.second * 7);
				break;
			case -7:
				SH = (3.0 / 64.0) * sqrt(715.0 / PI) * pow(sin(spherical.first), 7) * sin(spherical.second * 7);
				break;
            default:
                err_not_impl_f("Wrong spherical harmonic called!", std::cout);
			}
            break;
        case 8:
            leng = sqrt(x * x + y * y + z * z);
            spherical = constants::norm_cartesian_to_spherical(x / leng, y / leng, z / leng);
			switch (m)
			{
			case 0:
				SH = (1.0 / 256.0) * sqrt(17.0 / PI) * (6435.0 * pow(cos(spherical.first), 8) - 12012.0 * pow(cos(spherical.first), 6) + 6930.0 * pow(cos(spherical.first), 4) - 1260.0 * pow(cos(spherical.first), 2) + 35.0);
                break;
			case 1:
				SH = (3.0 / 64.0) * sqrt(17.0 / PI) * sin(spherical.first) * cos(spherical.first) * (715.0 * pow(cos(spherical.first), 6) - 1001.0 * pow(cos(spherical.first), 4) + 385.0 * pow(cos(spherical.first), 2) - 35.0) * cos(spherical.second);
				break;
			case -1:
				SH = (3.0 / 64.0) * sqrt(17.0 / PI) * sin(spherical.first) * cos(spherical.first) * (715.0 * pow(cos(spherical.first), 6) - 1001.0 * pow(cos(spherical.first), 4) + 385.0 * pow(cos(spherical.first), 2) - 35.0) * sin(spherical.second);
				break;
			case 2:
				SH = (3.0 / 64.0) * sqrt(595.0 / (2 * PI)) * pow(sin(spherical.first), 2) * (143 * pow(cos(spherical.first), 6) - 143 * pow(cos(spherical.first), 4) + 33 * pow(cos(spherical.first), 2) - 1) * cos(spherical.second * 2);
				break;
			case -2:
				SH = (3.0 / 1024.0) * sqrt(595.0 / (2 * PI)) * pow(sin(spherical.first), 2) * (385 * cos(spherical.first * 2) + 286 * cos(spherical.first * 4) + 143 * cos(spherical.first * 6) + 210) * sin(spherical.second) * cos(spherical.second);
				break;
			case 3:
				SH = (1.0 / 64.0) * sqrt(19635.0 / PI) * pow(sin(spherical.first), 3) * cos(spherical.first) * (39 * pow(cos(spherical.first), 4) - 26 * pow(cos(spherical.first), 2) + 3) * cos(spherical.second * 3);
				break;
			case -3:
				SH = (1.0 / 64.0) * sqrt(19635.0 / PI) * pow(sin(spherical.first), 3) * cos(spherical.first) * (39 * pow(cos(spherical.first), 4) - 26 * pow(cos(spherical.first), 2) + 3) * sin(spherical.second * 3);
				break;
			case 4:
				SH = (3.0 / 128.0) * sqrt(1309.0 / PI) * pow(sin(spherical.first), 4) * (65 * pow(cos(spherical.first), 4) - 26 * pow(cos(spherical.first), 2) + 1) * cos(spherical.second * 4); //TEST
				break;
			case -4:
				SH = (3.0 / 128.0) * sqrt(1309.0 / PI) * pow(sin(spherical.first), 4) * (65 * pow(cos(spherical.first), 4) - 26 * pow(cos(spherical.first), 2) + 1) * sin(spherical.second * 4);
				break;
			case 5:
				SH = (3.0 / 64.0) * sqrt(17017.0 / PI) * pow(sin(spherical.first), 5) * cos(spherical.first) * (5 * pow(cos(spherical.first), 2) - 1) * cos(spherical.second * 5);
				break;
			case -5:
				SH = (3.0 / 64.0) * sqrt(17017.0 / PI) * pow(sin(spherical.first), 5) * cos(spherical.first) * (5 * pow(cos(spherical.first), 2) - 1) * sin(spherical.second * 5);
				break;
			case 6:
                SH = (1.0 / 64.0) * sqrt(7293.0 / (2 * PI)) * pow(sin(spherical.first), 6) * (15 * pow(cos(spherical.first), 2) - 1) * cos(spherical.second * 6);
				break;
			case -6:
				SH = (1.0 / 64.0) * sqrt(7293.0 / (2 * PI)) * pow(sin(spherical.first), 6) * (15 * pow(cos(spherical.first), 2) - 1) * sin(spherical.second * 6);
				break;
			case 7:
				SH = (3.0 / 64.0) * sqrt(12155.0 / PI) * pow(sin(spherical.first), 7) * cos(spherical.first) * cos(spherical.second * 7);
				break;
			case -7:
				SH = (3.0 / 64.0) * sqrt(12155.0 / PI) * pow(sin(spherical.first), 7) * cos(spherical.first) * sin(spherical.second * 7);
				break;
			case 8:
				SH = (3.0 / 256.0) * sqrt(12155.0 / PI) * pow(sin(spherical.first), 8) * cos(spherical.second * 8);
				break;
			case -8:
				SH = (3.0 / 256.0) * sqrt(12155.0 / PI) * pow(sin(spherical.first), 8) * sin(spherical.second * 8);
				break;
			default:
				err_not_impl_f("Wrong spherical harmonic called!", std::cout);
			}
			break;
        case 9:
            leng = sqrt(x * x + y * y + z * z);
            spherical = constants::norm_cartesian_to_spherical(x / leng, y / leng, z / leng);
            switch (m) {
                case 0:
                    // 1/256 sqrt(19/PI) (12155 cos^9(theta) - 25740 cos^7(theta) + 18018 cos^5(theta) - 4620 cos^3(theta) + 315 cos(theta))
                    SH = (1.0 / 256.0) * sqrt(19.0 / PI) * (12155 * pow(cos(spherical.first), 9) - 25740 * pow(cos(spherical.first), 7) + 18018 * pow(cos(spherical.first), 5) - 4620 * pow(cos(spherical.first), 3) + 315 * cos(spherical.first));
                    break;
                case 1:
                    // 3/256 sqrt(95/PI) sin(theta) (2431 cos^8(theta) - 4004 cos^6(theta) + 2002 cos^4(theta) - 308 cos^2(theta) + 7) cos(phi) UIUIUIUI
					SH = (3.0 / 256.0) * sqrt(95.0 / PI) * sin(spherical.first) * (2431 * pow(cos(spherical.first), 8) - 4004 * pow(cos(spherical.first), 6) + 2002 * pow(cos(spherical.first), 4) - 308 * pow(cos(spherical.first), 2) + 7) * cos(spherical.second);
                    break;
                case -1:
                    SH = (3.0 / 256.0) * sqrt(95.0 / PI) * sin(spherical.first) * (2431 * pow(cos(spherical.first), 8) - 4004 * pow(cos(spherical.first), 6) + 2002 * pow(cos(spherical.first), 4) - 308 * pow(cos(spherical.first), 2) + 7) * sin(spherical.second);
                    break;
                case 2:
                    // 3/64 sqrt(1045/(2 PI)) sin^2(theta) cos(theta) (221 cos^6(theta) - 273 cos^4(theta) + 91 cos^2(theta) - 7) cos(2 phi) UIUIUIUI
					SH = (3.0 / 64.0) * sqrt(1045.0 / (2 * PI)) * pow(sin(spherical.first), 2) * cos(spherical.first) * (221 * pow(cos(spherical.first), 6) - 273 * pow(cos(spherical.first), 4) + 91 * pow(cos(spherical.first), 2) - 7) * cos(spherical.second * 2);
                    break;
                case -2:
                    SH = (3.0 / 64.0) * sqrt(1045.0 / (2 * PI)) * pow(sin(spherical.first), 2) * cos(spherical.first) * (221 * pow(cos(spherical.first), 6) - 273 * pow(cos(spherical.first), 4) + 91 * pow(cos(spherical.first), 2) - 7) * sin(spherical.second * 2);
                    break;
                case 3:
                    // 1/128 sqrt(21945/(2 PI)) sin^3(theta) (221 cos^6(theta) - 195 cos^4(theta) + 39 cos^2(theta) - 1) cos(3 phi)  UIUIUIUI
					SH = (1.0 / 128.0) * sqrt(21945.0 / (2 * PI)) * pow(sin(spherical.first), 3) * (221 * pow(cos(spherical.first), 6) - 195 * pow(cos(spherical.first), 4) + 39 * pow(cos(spherical.first), 2) - 1) * cos(spherical.second * 3);
                    break;
                case -3:
                    SH = (1.0 / 128.0) * sqrt(21945.0 / (2 * PI)) * pow(sin(spherical.first), 3) * (221 * pow(cos(spherical.first), 6) - 195 * pow(cos(spherical.first), 4) + 39 * pow(cos(spherical.first), 2) - 1) * sin(spherical.second * 3);
                    break;
                case 4:
                    // 3/128 sqrt(95095/PI) sin^4(theta) cos(theta) (17 cos^4(theta) - 10 cos^2(theta) + 1) cos(4 phi)
					SH = (3.0 / 128.0) * sqrt(95095.0 / PI) * pow(sin(spherical.first), 4) * cos(spherical.first) * (17 * pow(cos(spherical.first), 4) - 10 * pow(cos(spherical.first), 2) + 1) * cos(spherical.second * 4);
                    break;
                case -4:
                    SH = (3.0 / 128.0) * sqrt(95095.0 / PI) * pow(sin(spherical.first), 4) * cos(spherical.first) * (17 * pow(cos(spherical.first), 4) - 10 * pow(cos(spherical.first), 2) + 1) * sin(spherical.second * 4);
                    break;
                case 5:
                    // 3/128 sqrt(2717/(2 PI)) sin^5(theta) (85 cos^4(theta) - 30 cos^2(theta) + 1) cos(5 phi)
					SH = (3.0 / 128.0) * sqrt(2717.0 / (2 * PI)) * pow(sin(spherical.first), 5) * (85 * pow(cos(spherical.first), 4) - 30 * pow(cos(spherical.first), 2) + 1) * cos(spherical.second * 5);
                    break;
                case -5:
                    SH = (3.0 / 128.0) * sqrt(2717.0 / (2 * PI)) * pow(sin(spherical.first), 5) * (85 * pow(cos(spherical.first), 4) - 30 * pow(cos(spherical.first), 2) + 1) * sin(spherical.second * 5);
                    break;
                case 6:
                    // 1/64 sqrt(40755/(2 PI)) sin^6(theta) cos(theta) (17 cos^2(theta) - 3) cos(6 phi)
					SH = (1.0 / 64.0) * sqrt(40755.0 / (2 * PI)) * pow(sin(spherical.first), 6) * cos(spherical.first) * (17 * pow(cos(spherical.first), 2) - 3) * cos(spherical.second * 6);
                    break;
                case -6:
                    SH = (1.0 / 64.0) * sqrt(40755.0 / (2 * PI)) * pow(sin(spherical.first), 6) * cos(spherical.first) * (17 * pow(cos(spherical.first), 2) - 3) * sin(spherical.second * 6);
                    break;
                case 7:
                    // 3/256 sqrt(13585/(2 PI)) sin^7(theta) (17 cos^2(theta) - 1) cos(7 phi)
					SH = (3.0 / 256.0) * sqrt(13585.0 / (2 * PI)) * pow(sin(spherical.first), 7) * (17 * pow(cos(spherical.first), 2) - 1) * cos(spherical.second * 7);
                    break;
                case -7:
                    SH = (3.0 / 256.0) * sqrt(13585.0 / (2 * PI)) * pow(sin(spherical.first), 7) * (17 * pow(cos(spherical.first), 2) - 1) * sin(spherical.second * 7);
                    break;
                case 8:
                    // 3/256 sqrt(230945/PI) sin^8(theta) cos(theta) cos(8 phi)
					SH = (3.0 / 256.0) * sqrt(230945.0 / PI) * pow(sin(spherical.first), 8) * cos(spherical.first) * cos(spherical.second * 8);
                    break;
                case -8:
                    SH = (3.0 / 256.0) * sqrt(230945.0 / PI) * pow(sin(spherical.first), 8) * cos(spherical.first) * sin(spherical.second * 8);
                    break;
                case 9:
                    // 1/256 sqrt(230945/(2 PI)) sin^9(theta) cos(9 phi)
					SH = (1.0 / 256.0) * sqrt(230945.0 / (2 * PI)) * pow(sin(spherical.first), 9) * cos(spherical.second * 9);
                    break;
                case -9:
                    // 1/256 sqrt(230945/(2 PI)) sin^9(theta) sin(9 phi)
                    SH = (1.0 / 256.0) * sqrt(230945.0 / (2 * PI)) * pow(sin(spherical.first), 9) * sin(spherical.second * 9);
                    break;
                default:
                    err_not_impl_f("Wrong spherical harmonic called!", std::cout);
            }
            break;
        default:
            err_not_impl_f("Higehr than l=7 not done for spherical harmonic!", std::cout);
        }
        return SH;
    }

    // https://www.wolframalpha.com/input?i=LegendreP%5B6%2C6%2Cx%5D watch the signs!
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
        case 6: 
            switch (m) {
            case 0:
                return (1.0 / 16.0) * (231.0 * x * x * x * x * x * x - 315.0 * x * x * x * x + 105 * x * x - 5);
                break;
            case 1:
                return (21.0 / 8.0) * sqrt(1 - x * x) * (33 * x * x * x * x * x - 30 * x * x * x + 5 * x);
                break;
            case 2:
                return -(105.0 / 8.0) * (x * x - 1) * (33 * x * x * x * x - 18 * x * x + 1);
                break;
            case 3:
                return (315.0 / 2.0) * pow(1 - x * x, 1.5) * (11 * x * x * x - 3 * x);
                break;
            case 4:
                return (945.0 / 2.0) * pow(x * x - 1.0, 2.0) * (11 * x * x - 1);
                break;
            case 5:
                return 10395.0 * x * pow(1 - x * x, 2.5);
                break;
            case 6:
                return -10395.0 * pow(x * x - 1, 3.0);
                break;
            case -1:
                return -(1.0 / 16.0) * sqrt(1 - x * x) * (33 * x * x * x * x * x - 30 * x * x * x + 5 * x);
                break;
            case -2:
                return  (1.0 / 128.0) * (x * x - 1) * (33 * x * x * x * x - 18 * x * x + 1);
                break;
            case -3:
                return -(1.0 / 384.0) * pow(1 - x * x, 1.5) * (11 * x * x * x - 3 * x);
                break;
            case -4:
                return -(1.0 / 3840.0) * pow(x * x - 1, 2.0) * (11 * x * x - 1);
                break;
            case -5:
                return -(1.0 / 3840.0) * x * pow(1 - x * x, 2.5);
                break;
            case -6:
                return (1.0 / 46080.0) * pow(x * x - 1, 3.0);
                break;
            default:
                err_checkf(false, "This is impossible!", std::cout);
            }
            break;
        case 7: //Values are correct, only sign may be wrong
            switch (m)
            {
            case 0:
                return (1.0 / 16.0) * (429 * x * x * x * x * x * x * x - 693 * x * x * x * x * x + 315 * x * x * x - 35 * x);
                break;
            case 1:
                return (7.0 / 16.0) * sqrt(1 - x * x) * (429 * x * x * x * x * x * x - 495 * x * x * x * x + 135 * x * x - 5);
                break;
            case 2:
                return -(63.0 / 8.0) * (x * x - 1) * (143 * x * x * x * x * x - 110 * x * x * x + 15 * x);
                break;
            case 3:
                return (315.0 / 8.0) * pow(1 - x * x, 1.5) * (143 * x * x * x * x - 66 * x * x + 3);
                break;
            case 4:
                return (3465.0 / 2.0) * pow(x * x - 1, 2.0) * (13 * x * x * x - 3 * x);
                break;
            case 5:
                return (10395.0 / 2.0) * pow(1 - x * x, 2.5) * (13 * x * x - 1);
                break;
            case 6:
                return -135135.0 * pow(x * x - 1, 3.0) * x;
                break;
            case 7:
                return 135135.0 * pow(1 - x * x, 3.5);
                break;
            case -1:
                return -(1.0 / 128.0) * sqrt(1 - x * x) * (429 * x * x * x * x * x * x - 495 * x * x * x * x + 135 * x * x - 5);
                break;
            case -2:
                return (1.0 / 384.0) * (x * x - 1) * (143 * x * x * x * x * x - 110 * x * x * x + 15 * x);
                break;
            case -3:
                return -(1.0 / 3840.0) * pow(1 - x * x, 1.5) * (143 * x * x * x * x - 66 * x * x + 3);
                break;
            case -4:
                return -(1.0 / 3840.0) * pow(x * x - 1, 2.0) * (13 * x * x * x - 3 * x);
                break;
            case -5:
                return -(1.0 / 46080.0) * pow(1 - x * x, 2.5) * (13 * x * x - 1);
                break;
            case -6:
                return (1.0 / 46080.0) * pow(x * x - 1.0, 3.0) * x;
                break;
            case -7:
                return -(1.0 / 645120.0) * pow(1 - x * x, 3.5);
                break;
            default:
                err_checkf(false, "This is impossible!", std::cout);
            }
        case 8:
			switch (m)
			{
			case 0:
				return (1.0 / 128.0) * (6435 * x * x * x * x * x * x * x * x - 12012 * x * x * x * x * x * x + 6930 * x * x * x * x - 1260 * x * x + 35);
			case 1:
				return (9.0 / 16.0) * sqrt(1 - x * x) * (715 * x * x * x * x * x * x * x - 1001 * x * x * x * x * x + 385 * x * x * x - 35 * x);
			case 2:
				return -(315.0 / 16.0) * (x * x - 1) * (143 * x * x * x * x * x * x - 143 * x * x * x * x + 33 * x * x - 1);
			case 3:
				return (3465.0 / 8.0) * pow(1 - x * x, 1.5) * (39 * x * x * x * x * x - 26 * x * x * x + 3 * x);
			case 4:
				return (10395.0 / 8.0) * pow(x * x - 1, 2.0) * (65 * x * x * x * x - 26 * x * x + 1);
			case 5:
				return (135135.0 / 2.0) * pow(1 - x * x, 2.5) * (5 * x * x * x - x);
			case 6:
				return -(135135.0 / 2.0) * pow(x * x - 1, 3.0) * (15 * x * x - 1);
			case 7:
				return 2027025.0 * x * pow(1 - x * x, 3.5);
			case 8:
				return 2027025.0 * pow(x * x - 1, 4.0);
			case -1:
				return -(1.0 / 128.0) * sqrt(1 - x * x) * (715 * x * x * x * x * x * x * x - 1001 * x * x * x * x * x + 385 * x * x * x - 35 * x);
			case -2:
				return (1.0 / 256.0) * (x * x - 1) * (143 * x * x * x * x * x * x - 143 * x * x * x * x + 33 * x * x - 1);
			case -3:
				return -(1.0 / 768.0) * pow(1 - x * x, 1.5) * (39 * x * x * x * x * x - 26 * x * x * x + 3 * x);
			case -4:
				return -(1.0 / 15360.0) * pow(x * x - 1, 2.0) * (65 * x * x * x * x - 26 * x * x + 1);
			case -5:
				return -(1.0 / 15360.0) * pow(1 - x * x, 2.5) * (5 * x * x * x - x);
			case -6:
				return (1.0 / 645120.0) * pow(x * x - 1, 3.0) * (15 * x * x - 1);
			case -7:
				return -(1.0 / 645120.0) * x * pow(1 - x * x, 3.5);
			case -8:
				return -(1.0 / 10321920.0) * pow(x * x - 1, 4.0);
            default:
                err_checkf(false, "This is impossible!", std::cout);
			}
        case 9:
			switch (m)
			{
			case 0:
                // 1/128 (12155 x^9 - 25740 x^7 + 18018 x^5 - 4620 x^3 + 315 x)
				return (1.0 / 128.0) * (12155 * pow(x, 9) - 25740 * pow(x, 7) + 18018 * pow(x, 5) - 4620 * pow(x, 3) + 315 * x);
			case 1:
                // -45/128 sqrt(1 - x^2) (2431 x^8 - 4004 x^6 + 2002 x^4 - 308 x^2 + 7)
				return (45.0 / 128.0) * sqrt(1 - x * x) * (2431 * pow(x, 8) - 4004 * pow(x, 6) + 2002 * pow(x, 4) - 308 * pow(x, 2) + 7);
			case 2:
                // -495/16 (x^2 - 1) (221 x^7 - 273 x^5 + 91 x^3 - 7 x)
				return -(495.0 / 16.0) * (x * x - 1) * (221 * pow(x, 7) - 273 * pow(x, 5) + 91 * pow(x, 3) - 7 * x);
			case 3:
                // -3465/16 (1 - x^2)^(3/2) (221 x^6 - 195 x^4 + 39 x^2 - 1)
				return (3465.0 / 16.0) * pow(1 - x * x, 1.5) * (221 * pow(x, 6) - 195 * pow(x, 4) + 39 * pow(x, 2) - 1);
			case 4:
                // 135135/8 (x^2 - 1)^2 (17 x^5 - 10 x^3 + x)
				return (135135.0 / 8.0) * (x * x - 1) * (x * x - 1) * (17 * pow(x, 5) - 10 * pow(x, 3) + x);
			case 5:
                // -135135/8 (1 - x^2)^(5/2) (85 x^4 - 30 x^2 + 1)
				return (135135.0 / 8.0) * pow(1 - x * x, 2.5) * (85 * pow(x, 4) - 30 * pow(x, 2) + 1);
			case 6:
                // -675675/2 (x^2 - 1)^3 (17 x^3 - 3 x)
				return -(675675.0 / 2.0) * pow(x * x - 1, 3) * (17 * pow(x, 3) - 3 * x);
			case 7:
                // -2027025/2 (1 - x^2)^(7/2) (17 x^2 - 1)
				return (2027025.0 / 2.0) * pow(1 - x * x, 3.5) * (17 * pow(x, 2) - 1);
			case 8:
                // 34459425 x (x^2 - 1)^4
				return 34459425 * x * pow(x * x - 1, 4);
			case 9:
                // -34459425 (1 - x^2)^(9/2)
				return 34459425 * pow(1 - x * x, 4.5);
			case -1:
                // 1/256 sqrt(1 - x^2) (2431 x^8 - 4004 x^6 + 2002 x^4 - 308 x^2 + 7)
				return -(1.0 / 256.0) * sqrt(1 - x * x) * (2431 * pow(x, 8) - 4004 * pow(x, 6) + 2002 * pow(x, 4) - 308 * pow(x, 2) + 7);
			case -2:
                // -1/256 (x^2 - 1) (221 x^7 - 273 x^5 + 91 x^3 - 7 x)
				return (1.0 / 256.0) * (x * x - 1) * (221 * pow(x, 7) - 273 * pow(x, 5) + 91 * pow(x, 3) - 7 * x);
			case -3:
                // ((1 - x^2)^(3/2) (221 x^6 - 195 x^4 + 39 x^2 - 1))/3072
				return -(1.0 / 3072.0) * pow(1 - x * x, 1.5) * (221 * pow(x, 6) - 195 * pow(x, 4) + 39 * pow(x, 2) - 1);
			case -4:
                // ((x^2 - 1)^2 (17 x^5 - 10 x^3 + x))/3072
				return -(1.0 / 3072.0) * (x * x - 1) * (x * x - 1) * (17 * pow(x, 5) - 10 * pow(x, 3) + x);
			case -5:
                // ((1 - x^2)^(5/2) (85 x^4 - 30 x^2 + 1))/215040
				return -(1.0 / 215040.0) * pow(1 - x * x, 2.5) * (85 * pow(x, 4) - 30 * pow(x, 2) + 1);
			case -6:
                // -((x^2 - 1)^3 (17 x^3 - 3 x))/645120
				return (1.0 / 645120.0) * pow(x * x - 1, 3) * (17 * pow(x, 3) - 3 * x);
			case -7:
                // ((1 - x^2)^(7/2) (17 x^2 - 1))/10321920
				return -(1.0 / 10321920.0) * pow(1 - x * x, 3.5) * (17 * pow(x, 2) - 1);
			case -8:
                // (x (x^2 - 1)^4)/10321920
				return -x * pow(x * x - 1, 4) / 10321920;
			case -9:
                // (1 - x^2)^(9/2)/185794560
				return -pow(1 - x * x, 4.5) / 185794560;
			default:
				err_checkf(false, "This is impossible!", std::cout);
			}
        default:
            err_not_impl_f("associated_legendre_polynomial l > 9 ", std::cout);
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

    //Returns the spherical coordinates of a given cartesian vector
    //the output is a vector with the (radius, theta and phi)
    std::pair<double,double> norm_cartesian_to_spherical(const double& x, const double& y, const double& z) {
        return { acos(z), atan2(y, x) };
    }

    const vec2 spherical_norms{
        {sqrt(1 / constants::FOUR_PI)}, // l=0
        {sqrt(3 / constants::TWO_PI * double(constants::ft[2]) / double(constants::ft[0])),
         sqrt(3 / constants::FOUR_PI),
         sqrt(3 / constants::TWO_PI * double(constants::ft[0]) / double(constants::ft[2]))}, // l=1 m=-1, m=0, m=1
        {sqrt(5 / constants::TWO_PI * double(constants::ft[4]) / double(constants::ft[0])),
         sqrt(5 / constants::TWO_PI * double(constants::ft[3]) / double(constants::ft[1])),
         sqrt(5 / constants::FOUR_PI) ,
         sqrt(5 / constants::TWO_PI * double(constants::ft[1]) / double(constants::ft[3])),
         sqrt(5 / constants::TWO_PI * double(constants::ft[0]) / double(constants::ft[4]))}, // l=2 m=-2, m=-1, m=0, m=1, m=2
        {sqrt(7 / constants::TWO_PI * double(constants::ft[6]) / double(constants::ft[0])),
         sqrt(7 / constants::TWO_PI * double(constants::ft[5]) / double(constants::ft[1])),
         sqrt(7 / constants::TWO_PI * double(constants::ft[4]) / double(constants::ft[2])),
         sqrt(7 / constants::FOUR_PI),
         sqrt(7 / constants::TWO_PI * double(constants::ft[2]) / double(constants::ft[4])),
         sqrt(7 / constants::TWO_PI * double(constants::ft[1]) / double(constants::ft[5])),
         sqrt(7 / constants::TWO_PI * double(constants::ft[0]) / double(constants::ft[6]))}, // l=3 m=-3, m=-2, m=-1, m=0, m=1, m=2, m=3
        {sqrt(9 / constants::TWO_PI * double(constants::ft[8]) / double(constants::ft[0])),
         sqrt(9 / constants::TWO_PI * double(constants::ft[7]) / double(constants::ft[1])),
         sqrt(9 / constants::TWO_PI * double(constants::ft[6]) / double(constants::ft[2])),
         sqrt(9 / constants::TWO_PI * double(constants::ft[5]) / double(constants::ft[3])),
         sqrt(9 / constants::FOUR_PI),
         sqrt(9 / constants::TWO_PI * double(constants::ft[3]) / double(constants::ft[5])),
         sqrt(9 / constants::TWO_PI * double(constants::ft[2]) / double(constants::ft[6])),
         sqrt(9 / constants::TWO_PI * double(constants::ft[1]) / double(constants::ft[7])),
         sqrt(9 / constants::TWO_PI * double(constants::ft[0]) / double(constants::ft[8]))}, // l=4 m=-4, m=-3, m=-2, m=-1, m=0, m=1, m=2, m=3, m=4
        {sqrt(11 / constants::TWO_PI * double(constants::ft[10]) / double(constants::ft[0])),
         sqrt(11 / constants::TWO_PI * double(constants::ft[9]) / double(constants::ft[1])),
         sqrt(11 / constants::TWO_PI * double(constants::ft[8]) / double(constants::ft[2])),
         sqrt(11 / constants::TWO_PI * double(constants::ft[7]) / double(constants::ft[3])),
         sqrt(11 / constants::TWO_PI * double(constants::ft[6]) / double(constants::ft[4])),
         sqrt(11 / constants::FOUR_PI),
         sqrt(11 / constants::TWO_PI * double(constants::ft[4]) / double(constants::ft[6])),
         sqrt(11 / constants::TWO_PI * double(constants::ft[3]) / double(constants::ft[7])),
         sqrt(11 / constants::TWO_PI * double(constants::ft[2]) / double(constants::ft[8])),
         sqrt(11 / constants::TWO_PI * double(constants::ft[1]) / double(constants::ft[9])),
         sqrt(11 / constants::TWO_PI * double(constants::ft[0]) / double(constants::ft[10]))}, // l=5 m=-5, m=-4, m=-3, m=-2, m=-1, m=0, m=1, m=2, m=3, m=4, m=5
        {sqrt(13 / constants::TWO_PI * double(constants::ft[12]) / double(constants::ft[0])),
         sqrt(13 / constants::TWO_PI * double(constants::ft[11]) / double(constants::ft[1])),
         sqrt(13 / constants::TWO_PI * double(constants::ft[10]) / double(constants::ft[2])),
         sqrt(13 / constants::TWO_PI * double(constants::ft[9]) / double(constants::ft[3])),
         sqrt(13 / constants::TWO_PI * double(constants::ft[8]) / double(constants::ft[4])),
         sqrt(13 / constants::TWO_PI * double(constants::ft[7]) / double(constants::ft[5])),
         sqrt(13 / constants::FOUR_PI),
         sqrt(13 / constants::TWO_PI * double(constants::ft[5]) / double(constants::ft[7])),
         sqrt(13 / constants::TWO_PI * double(constants::ft[4]) / double(constants::ft[8])),
         sqrt(13 / constants::TWO_PI * double(constants::ft[3]) / double(constants::ft[9])),
         sqrt(13 / constants::TWO_PI * double(constants::ft[2]) / double(constants::ft[10])),
         sqrt(13 / constants::TWO_PI * double(constants::ft[1]) / double(constants::ft[11])),
         sqrt(13 / constants::TWO_PI * double(constants::ft[0]) / double(constants::ft[12]))}, // l=6 m=-6, m=-5, m=-4, m=-3, m=-2, m=-1, m=0, m=1, m=2, m=3, m=4, m=5, m=6
        {sqrt(15 / constants::TWO_PI * double(constants::ft[14]) / double(constants::ft[0])),
		sqrt(15 / constants::TWO_PI * double(constants::ft[13]) / double(constants::ft[1])),
		sqrt(15 / constants::TWO_PI * double(constants::ft[12]) / double(constants::ft[2])),
		sqrt(15 / constants::TWO_PI * double(constants::ft[11]) / double(constants::ft[3])),
		sqrt(15 / constants::TWO_PI * double(constants::ft[10]) / double(constants::ft[4])),
		sqrt(15 / constants::TWO_PI * double(constants::ft[9]) / double(constants::ft[5])),
		sqrt(15 / constants::TWO_PI * double(constants::ft[8]) / double(constants::ft[6])),
		sqrt(15 / constants::FOUR_PI),
		sqrt(15 / constants::TWO_PI * double(constants::ft[6]) / double(constants::ft[8])),
		sqrt(15 / constants::TWO_PI * double(constants::ft[5]) / double(constants::ft[9])),
		sqrt(15 / constants::TWO_PI * double(constants::ft[4]) / double(constants::ft[10])),
		sqrt(15 / constants::TWO_PI * double(constants::ft[3]) / double(constants::ft[11])),
		sqrt(15 / constants::TWO_PI * double(constants::ft[2]) / double(constants::ft[12])),
		sqrt(15 / constants::TWO_PI * double(constants::ft[1]) / double(constants::ft[13])),
		sqrt(15 / constants::TWO_PI * double(constants::ft[0]) / double(constants::ft[14]))}, // l=7 m=-7, m=-6, m=-5, m=-4, m=-3, m=-2, m=-1, m=0, m=1, m=2, m=3, m=4, m=5, m=6, m=7
		{sqrt(17 / constants::TWO_PI * double(constants::ft[16]) / double(constants::ft[0])),
		sqrt(17 / constants::TWO_PI * double(constants::ft[15]) / double(constants::ft[1])),
		sqrt(17 / constants::TWO_PI * double(constants::ft[14]) / double(constants::ft[2])),
		sqrt(17 / constants::TWO_PI * double(constants::ft[13]) / double(constants::ft[3])),
		sqrt(17 / constants::TWO_PI * double(constants::ft[12]) / double(constants::ft[4])),
		sqrt(17 / constants::TWO_PI * double(constants::ft[11]) / double(constants::ft[5])),
		sqrt(17 / constants::TWO_PI * double(constants::ft[10]) / double(constants::ft[6])),
		sqrt(17 / constants::TWO_PI * double(constants::ft[9]) / double(constants::ft[7])),
		sqrt(17 / constants::FOUR_PI),
		sqrt(17 / constants::TWO_PI * double(constants::ft[7]) / double(constants::ft[9])),
		sqrt(17 / constants::TWO_PI * double(constants::ft[6]) / double(constants::ft[10])),
		sqrt(17 / constants::TWO_PI * double(constants::ft[5]) / double(constants::ft[11])),
		sqrt(17 / constants::TWO_PI * double(constants::ft[4]) / double(constants::ft[12])),
		sqrt(17 / constants::TWO_PI * double(constants::ft[3]) / double(constants::ft[13])),
		sqrt(17 / constants::TWO_PI * double(constants::ft[2]) / double(constants::ft[14])),
		sqrt(17 / constants::TWO_PI * double(constants::ft[1]) / double(constants::ft[15])),
		sqrt(17 / constants::TWO_PI * double(constants::ft[0]) / double(constants::ft[16]))}, // l=8 m=-8, m=-7, m=-6, m=-5, m=-4, m=-3, m=-2, m=-1, m=0, m=1, m=2, m=3, m=4, m=5, m=6, m=7, m=8
		{sqrt(19 / constants::TWO_PI * double(constants::ft[18]) / double(constants::ft[0])),
		sqrt(19 / constants::TWO_PI * double(constants::ft[17]) / double(constants::ft[1])),
		sqrt(19 / constants::TWO_PI * double(constants::ft[16]) / double(constants::ft[2])),
		sqrt(19 / constants::TWO_PI * double(constants::ft[15]) / double(constants::ft[3])),
		sqrt(19 / constants::TWO_PI * double(constants::ft[14]) / double(constants::ft[4])),
		sqrt(19 / constants::TWO_PI * double(constants::ft[13]) / double(constants::ft[5])),
		sqrt(19 / constants::TWO_PI * double(constants::ft[12]) / double(constants::ft[6])),
		sqrt(19 / constants::TWO_PI * double(constants::ft[11]) / double(constants::ft[7])),
		sqrt(19 / constants::TWO_PI * double(constants::ft[10]) / double(constants::ft[8])),
		sqrt(19 / constants::FOUR_PI),
		sqrt(19 / constants::TWO_PI * double(constants::ft[8]) / double(constants::ft[10])),
		sqrt(19 / constants::TWO_PI * double(constants::ft[7]) / double(constants::ft[11])),
		sqrt(19 / constants::TWO_PI * double(constants::ft[6]) / double(constants::ft[12])),
		sqrt(19 / constants::TWO_PI * double(constants::ft[5]) / double(constants::ft[13])),
		sqrt(19 / constants::TWO_PI * double(constants::ft[4]) / double(constants::ft[14])),
		sqrt(19 / constants::TWO_PI * double(constants::ft[3]) / double(constants::ft[15])),
		sqrt(19 / constants::TWO_PI * double(constants::ft[2]) / double(constants::ft[16])),
		sqrt(19 / constants::TWO_PI * double(constants::ft[1]) / double(constants::ft[17])),
		sqrt(19 / constants::TWO_PI * double(constants::ft[0]) / double(constants::ft[18]))}, // l=9 m=-9, m=-8, m=-7, m=-6, m=-5, m=-4, m=-3, m=-2, m=-1, m=0, m=1, m=2, m=3, m=4, m=5, m=6, m=7, m=8, m=9
    };

    //Original implementation after P. Coppens DOI: 10.1107/97809553602060000759 Eq. 1.2.7.2b
    //I omitted the abs(m) in the factorial as most other sources do not include it
    double real_spherical(const int& l, const int& m, const double& theta, const double& phi) {
        //double N;
        //m == 0 ? N = sqrt((2 * l + 1) / constants::FOUR_PI) : N = sqrt(((2 * l + 1) / constants::TWO_PI) * double(constants::ft[l - (m)]) / double(constants::ft[l + (m)]));
        if (m >= 0) {
            return spherical_norms[l][l + m] * associated_legendre_polynomial(l, m, cos(theta)) * cos(m * phi);
        }
        else {
            return spherical_norms[l][l + m] * associated_legendre_polynomial(l, m, cos(theta)) * sin(m * phi);
        }
    }

    int get_closest_num_angular(const int& n)
    {
        int m = 0;

        for (int i = 0; i < max_LT; i++) {
            m = lebedev_table[i];
            if (m >= n)
                return m;
        }
        return -1;
    };

    int get_angular_order(const int& n)
    {
        for (int i = 0; i < max_LT; i++) {
            if (lebedev_table[i] == n)
                return i;
        }
        return -1;
    };
}