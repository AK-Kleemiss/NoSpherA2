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
                SH = constants::c_315_256p * ((x * x * (x * x - 3 * y * y)) - (y * y * (3 * x * x - y * y)));
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
            switch (m)
            {
            case 0: // I Z^6
                //1/32 sqrt(13/PI) (231 cos^6(theta) - 315 cos^4(theta) + 105 cos^2(theta) - 5)
                SH = (1.0 / 32.0) * std::sqrt(13.0 / PI) * (z * z * (z*z*(231.0 * z * z - 315.0) + 105.0) - 5.0);
                break;
            case 1:
                SH = (1.0 / 16.0) * std::sqrt(273.0 / (PI)) * x * z * (33 * z * z * z * z - 30 * z * z + 5);
                break;
            case -1:
                //-1/256 i sqrt(273/(2 PI)) (5 sin(2 theta) + 12 sin(4 theta) + 33 sin(6 theta)) sin(phi)
                SH = (1.0 / 16.0) * std::sqrt(273.0 / (PI)) * y * z * (33 * z * z * z * z - 30 * z * z + 5);
                break;
            case 2:
                SH = (1.0 / 32.0) * std::sqrt(1365.0 / (2 * PI)) * (33 * z * z * z * z - 18 * z * z + 1) * (x * x - y * y);
                break;
            case -2:
                // -1/256 i sqrt(1365/PI) sin^2(theta) (60 cos(2 theta) + 33 cos(4 theta) + 35) sin(2 phi) 
                SH = (1.0 / 256.0) * std::sqrt(1365.0 / (2 * PI)) * (60 * (2 * z * z - 1) + 33 * (8 * z * z * z * z - 8 * z * z + 1) + 35) * 2 * x * y;
                break;
            case 3:
				//SH = (1.0 / 16.0) * std::sqrt(1365.0 / (2 * PI)) * pow(sin(spherical.first), 3) * z * (11 * z*z - 3) * cos(spherical.second * 3);
                SH = (1.0 / 16.0) * std::sqrt(1365.0 / (2 * PI)) * z * (11 * z * z - 3) * (4 * x * x * x - 3 * x * (x * x + y * y));
                break;
            case -3:
                //-1/64 i sqrt(1365/PI) sin^3(theta) (21 cos(theta) + 11 cos(3 theta)) sin(3 phi)
                SH = (1.0 / 64.0) * std::sqrt(1365.0 / (2 * PI)) * (21 * z + 11 * (4 * z * z * z - 3 * z)) * (3 * y * (x * x + y * y) - 4 * y * y * y);
                break;
            case 4:
                SH = (3.0 / 64.0) * std::sqrt(91.0 / (PI)) * (2 * (11 * z * z - 1) * (x * x * x * x - 6 * x * x * y * y + y * y * y * y));
                break;
            case -4:
                //-3/32 i sqrt(91/(2 PI)) sin^4(theta) (11 cos(2 theta) + 9) sin(4 phi)
                SH = (3.0 / 32.0) * std::sqrt(91.0 / (4 * PI)) * (11 * (2 * z * z - 1) + 9) * 4 * x * y * (x * x - y * y);
                break;
            case 5:
                SH = (3.0 / 16.0) * std::sqrt(1001.0 / (2 * PI)) * z * (x * x * x * x * x - 10 * x * x * x * y * y + 5 * x * y * y * y * y);
                break;
            case -5:
                SH = (3.0 / 16.0) * std::sqrt(1001.0 / (2 * PI)) * z * (5 * x * x * x * x * y - 10 * x * x * y * y * y + y * y * y * y * y);
                break;
            case 6:
                SH = (1.0 / 32.0) * std::sqrt(3003.0 / (2 * PI)) * (-y * y * y * y * y * y + 15 * y * y * y * y * x * x - 15 * y * y * x * x * x * x + x * x * x * x * x * x);
                break;
            case -6:
                SH = (1.0 / 32.0) * std::sqrt(3003.0 / (2 * PI)) * (6 * x * x * x * x * x * y - 20 * x * x * x * y * y * y + 6 * x * y * y * y * y * y);
                break;
            default:
                err_not_impl_f("Wrong spherical harmonic called!", std::cout);
            }
            break;
		case 7:// Generated via Wolfram Alpha using: (SphericalHarmonicY[7,-m,theta,phi] + SphericalHarmonicY[7,m,theta,phi])  //BEWARE THE SIGN!!!!
			switch (m)
			{
			case 0:
                SH = (1.0 / 32.0) * std::sqrt(15.0 / PI) * z * (z * z * (z * z * (429.0 * z * z - 693.0) + 315.0) - 35.0);
                break;
			case 1:
                SH = (1.0 / 64.0) * std::sqrt(105.0 / PI) * x * (z * z * (z * z * (429.0 * z * z - 495.0) + 135.0) - 5.0);
				break;
			case -1:
                SH = (1.0 / 64.0) * std::sqrt(105.0 / PI) * y * (z * z * (z * z * (429.0 * z * z - 495.0) + 135.0) - 5.0);
				break;
			case 2:
                SH = (3.0 / 32.0) * std::sqrt(35.0 / (2 * PI)) * (x * x - y * y) * z * (z * z * (143.0 * z * z - 110.0) + 15.0);
				break;
			case -2:
                SH = (3.0 / 32.0) * std::sqrt(35.0 / (2 * PI)) * 2 * x * y * z * (z * z * (143.0 * z * z - 110.0) + 15.0);
				break;
			case 3:
                SH = (3.0 / 64.0) * std::sqrt(35.0 / PI) * x * (x * x - 3 * y * y) * (z * z * (143.0 * z * z - 66.0) + 3.0);
				break;
			case -3:
                SH = (3.0 / 64.0) * std::sqrt(35.0 / PI) * y * (3 * x * x - y * y) * (z * z * (143.0 * z * z - 66.0) + 3.0);
				break;
			case 4:
                SH = (3.0 / 32.0) * std::sqrt(385.0 / PI) * (x * x * (x * x - 6 * y * y) + y * y * y * y) * z * (13.0 * z * z - 3.0);
				break;
			case -4:
                SH = (3.0 / 32.0) * std::sqrt(385.0 / PI) * x * y * 4 * (x * x - y * y) * z * (13.0 * z * z - 3.0);
				break;
			case 5:
                SH = (3.0 / 64.0) * std::sqrt(385.0 / PI) * x * (x * x * x * x + 5 * y * y * (-2 * x * x + y * y)) * (13.0 * z * z - 1.0);
				break;
			case -5:
                SH = (3.0 / 64.0) * std::sqrt(385.0 / PI) * y * (5 * x * x * (x * x - 2 * y * y) + y * y * y * y) * (13.0 * z * z - 1.0);
				break;
			case 6:
                SH = (3.0 / 32.0) * std::sqrt(5005.0 / (2 * PI)) * (y * y * (y * y * (-y * y + 15 * x * x) - 15 * x * x * x * x) + x * x * x * x * x * x) * z;
				break;
			case -6:
                SH = (3.0 / 32.0) * std::sqrt(5005.0 / (2 * PI)) * x * y * (x * x * (6 * x * x - 20 * y * y) + 6 * y * y * y * y) * z;
				break;
			case 7:
                SH = (3.0 / 64.0) * std::sqrt(715.0 / PI) * x * (x * x * x * x * x * x + 7 * y * y * (-3 * x * x * x * x + y * y * (5 * x * x - y * y)));
				break;
			case -7:
                SH = (3.0 / 64.0) * std::sqrt(715.0 / PI) * y * (x * x * (x * x * (7 * x * x - 35 * y * y) + 21 * y * y * y * y) - y * y * y * y * y * y);
				break;
            default:
                err_not_impl_f("Wrong spherical harmonic called!", std::cout);
			}
            break;
        case 8:
			switch (m)
			{
			case 0:
                SH = (1.0 / 256.0) * std::sqrt(17.0 / PI) * (z * z * (z * z * (z * z * (6435.0 * z * z - 12012.0) + 6930.0) - 1260.0) + 35.0);
                break;
			case 1:
				SH = (3.0 / 64.0) * std::sqrt(17.0 / PI) * x * z * (715.0 * z*z*z*z*z*z - 1001.0 * z*z*z*z + 385.0 * z*z - 35.0);
				break;
			case -1:
				SH = (3.0 / 64.0) * std::sqrt(17.0 / PI)* y * z * (715.0 * z*z*z*z*z*z - 1001.0 * z*z*z*z + 385.0 * z*z - 35.0);
				break;
			case 2:
				SH = (3.0 / 64.0) * std::sqrt(595.0 / (2 * PI)) * (x * x - y * y) * (143 * z*z*z*z*z*z - 143 * z*z*z*z + 33 * z*z - 1);
				break;
			case -2:
                SH = (3.0 / 64.0) * std::sqrt(595.0 / (2 * PI)) * 2 * x * y * (143 * z * z * z * z * z * z - 143 * z * z * z * z + 33 * z * z - 1);
				break;
			case 3:
				SH = (1.0 / 64.0) * std::sqrt(19635.0 / PI) * (x * x * x - 3 * x * y * y) * z * (39 * z*z*z*z - 26 * z*z + 3);
				break;
			case -3:
				SH = (1.0 / 64.0) * std::sqrt(19635.0 / PI) * (3 * x * x * y - y * y * y) *z * (39 * z*z*z*z - 26 * z*z + 3);
				break;
			case 4:
				SH = (3.0 / 128.0) * std::sqrt(1309.0 / PI) * (x * x * (x * x - 6 * y * y) + y * y * y * y) * (65 * z*z*z*z - 26 * z*z + 1); //TEST
				break;
			case -4:
				SH = (3.0 / 128.0) * std::sqrt(1309.0 / PI) * x * y * 4 * (x * x - y * y) * (65 * z*z*z*z - 26 * z*z + 1) ;
				break;
			case 5:
				SH = (3.0 / 64.0) * std::sqrt(17017.0 / PI) * x * (x * x * x * x + 5 * y * y * (-2 * x * x + y * y)) *z * (5 * z*z - 1) ;
				break;
			case -5:
				SH = (3.0 / 64.0) * std::sqrt(17017.0 / PI) * y * (5 * x * x * (x * x - 2 * y * y) + y * y * y * y) *z * (5 * z*z - 1) ;
				break;
			case 6:
                SH = (1.0 / 64.0) * std::sqrt(7293.0 / (2 * PI)) * (y * y * (y * y * (-y * y + 15 * x * x) - 15 * x * x * x * x) + x * x * x * x * x * x) * (15 * z*z - 1);
				break;
			case -6:
				SH = (1.0 / 64.0) * std::sqrt(7293.0 / (2 * PI)) * x * y * (x * x * (6 * x * x - 20 * y * y) + 6 * y * y * y * y) * (15 * z*z - 1);
				break;
			case 7:
				SH = (3.0 / 64.0) * std::sqrt(12155.0 / PI) * x * (x * x * x * x * x * x + 7 * y * y * (-3 * x * x * x * x + y * y * (5 * x * x - y * y))) *z ;
				break;
			case -7:
				SH = (3.0 / 64.0) * std::sqrt(12155.0 / PI) * y * (x * x * (x * x * (7 * x * x - 35 * y * y) + 21 * y * y * y * y) - y * y * y * y * y * y) *z ;
				break;
			case 8:
                SH = (3.0 / 256.0) * std::sqrt(12155.0 / PI) * (y * y * (y * y * (y * y * y * y - 28 * x * x * y * y + 70 * x * x * x * x) - 28 * x * x * x * x * x * x) + x * x * x * x * x * x * x * x);
				break;
			case -8:
                SH = (3.0 / 256.0) * std::sqrt(12155.0 / PI) * y * x * (y * y * (y * y * (-8 * y * y + 56 * x * x) - 56 * x * x * x * x) + 8 * x * x * x * x * x * x);
				break;
			default:
				err_not_impl_f("Wrong spherical harmonic called!", std::cout);
			}
			break;
        case 9:
            switch (m) {
                case 0:
                    // 1/256 sqrt(19/PI) cos(theta) (12155 cos^8(theta) - 25740 cos^6(theta) + 18018 cos^4(theta) - 4620 cos^2(theta) + 315)
                    SH = (1.0 / 256.0) * std::sqrt(19.0 / PI) *z * (12155 * z*z*z*z*z*z*z*z - 25740 * z*z*z*z*z*z + 18018 * z*z*z*z - 4620 * z*z + 315 );
                    break;
                case 1:
                    // 3/256 sqrt(95/PI) sin(theta) (2431 cos^8(theta) - 4004 cos^6(theta) + 2002 cos^4(theta) - 308 cos^2(theta) + 7) cos(phi) UIUIUIUI
					SH = (3.0 / 256.0) * std::sqrt(95.0 / PI) * x * (2431 * z*z*z*z*z*z*z*z - 4004 * z*z*z*z*z*z + 2002 * z*z*z*z - 308 * z*z + 7) ;
                    break;
                case -1:
                    SH = (3.0 / 256.0) * std::sqrt(95.0 / PI) * y * (2431 * z*z*z*z*z*z*z*z - 4004 * z*z*z*z*z*z + 2002 * z*z*z*z - 308 * z*z + 7) ;
                    break;
                case 2:
                    // 3/64 sqrt(1045/(2 PI)) sin^2(theta) cos(theta) (221 cos^6(theta) - 273 cos^4(theta) + 91 cos^2(theta) - 7) cos(2 phi) UIUIUIUI
					SH = (3.0 / 64.0) * std::sqrt(1045.0 / (2 * PI)) * (x * x - y * y) *z * (221 * z*z*z*z*z*z - 273 * z*z*z*z + 91 * z*z - 7);
                    break;
                case -2:
                    SH = (3.0 / 64.0) * std::sqrt(1045.0 / (2 * PI)) * 2 * x * y *z * (221 * z*z*z*z*z*z - 273 * z*z*z*z + 91 * z*z - 7);
                    break;
                case 3:
                    // 1/128 sqrt(21945/(2 PI)) sin^3(theta) (221 cos^6(theta) - 195 cos^4(theta) + 39 cos^2(theta) - 1) cos(3 phi)  UIUIUIUI
					SH = (1.0 / 128.0) * std::sqrt(21945.0 / (2 * PI)) * x * (x * x - 3 * y * y) * (221 * z*z*z*z*z*z - 195 * z*z*z*z + 39 * z*z - 1);
                    break;
                case -3:
                    SH = (1.0 / 128.0) * std::sqrt(21945.0 / (2 * PI)) * y * (3 * x * x - y * y) * (221 * z*z*z*z*z*z - 195 * z*z*z*z + 39 * z*z - 1);
                    break;
                case 4:
                    // 3/128 sqrt(95095/PI) sin^4(theta) cos(theta) (17 cos^4(theta) - 10 cos^2(theta) + 1) cos(4 phi)
					SH = (3.0 / 128.0) * std::sqrt(95095.0 / PI) * (x * x * (x * x - 6 * y * y) + y * y * y * y) *z * (17 * z*z*z*z - 10 * z*z + 1);
                    break;
                case -4:
                    SH = (3.0 / 128.0) * std::sqrt(95095.0 / PI) * x * y * 4 * (x * x - y * y) *z * (17 * z*z*z*z - 10 * z*z + 1);
                    break;
                case 5:
                    // 3/128 sqrt(2717/(2 PI)) sin^5(theta) (85 cos^4(theta) - 30 cos^2(theta) + 1) cos(5 phi)
					SH = (3.0 / 128.0) * std::sqrt(2717.0 / (2 * PI)) * x * (x * x * x * x + 5 * y * y * (-2 * x * x + y * y)) * (85 * z*z*z*z - 30 * z*z + 1);
                    break;
                case -5:
                    SH = (3.0 / 128.0) * std::sqrt(2717.0 / (2 * PI)) * y * (5 * x * x * (x * x - 2 * y * y) + y * y * y * y) * (85 * z*z*z*z - 30 * z*z + 1);
                    break;
                case 6:
                    // 1/64 sqrt(40755/(2 PI)) sin^6(theta) cos(theta) (17 cos^2(theta) - 3) cos(6 phi)
					SH = (1.0 / 64.0) * std::sqrt(40755.0 / (2 * PI)) * (y * y * (y * y * (-y * y + 15 * x * x) - 15 * x * x * x * x) + x * x * x * x * x * x) *z * (17 * z*z - 3);
                    break;
                case -6:
                    SH = (1.0 / 64.0) * std::sqrt(40755.0 / (2 * PI)) * x * y * (x * x * (6 * x * x - 20 * y * y) + 6 * y * y * y * y) *z * (17 * z*z - 3);
                    break;
                case 7:
                    // 3/256 sqrt(13585/(2 PI)) sin^7(theta) (17 cos^2(theta) - 1) cos(7 phi)
					SH = (3.0 / 256.0) * std::sqrt(13585.0 / (2 * PI)) * x * (x * x * x * x * x * x + 7 * y * y * (-3 * x * x * x * x + y * y * (5 * x * x - y * y))) * (17 * z*z - 1);
                    break;
                case -7:
                    SH = (3.0 / 256.0) * std::sqrt(13585.0 / (2 * PI)) * y * (x * x * (x * x * (7 * x * x - 35 * y * y) + 21 * y * y * y * y) - y * y * y * y * y * y) * (17 * z*z - 1);
                    break;
                case 8:
                    // 3/256 sqrt(230945/PI) sin^8(theta) cos(theta) cos(8 phi)
					SH = (3.0 / 256.0) * std::sqrt(230945.0 / PI) * (y * y * (y * y * (y * y * y * y - 28 * x * x * y * y + 70 * x * x * x * x) - 28 * x * x * x * x * x * x) + x * x * x * x * x * x * x * x) *z;
                    break;
                case -8:
                    SH = (3.0 / 256.0) * std::sqrt(230945.0 / PI) * y * x * (y * y * (y * y * (-8 * y * y + 56 * x * x) - 56 * x * x * x * x) + 8 * x * x * x * x * x * x) *z;
                    break;
                case 9:
                    // 1/256 sqrt(230945/(2 PI)) sin^9(theta) cos(9 phi)
                    SH = (1.0 / 256.0) * std::sqrt(230945.0 / (2 * PI)) * x * (y*y * (y * y * (y * y *(9 * y * y - 84 * x * x) + 126 * x * x * x * x) - 36 * x * x * x * x * x * x) + x * x * x * x * x * x * x * x);
                    break;
                case -9:
                    // 1/256 sqrt(230945/(2 PI)) sin^9(theta) sin(9 phi)
                    SH = (1.0 / 256.0) * std::sqrt(230945.0 / (2 * PI)) * y * (9 * x * x * x * x * x * x * x * x + y * y * (-84 * x * x * x * x * x * x + y * y * (126 * x * x * x * x + y * y * (-36 * x * x + y * y))));
                    break;
                default:
                    err_not_impl_f("Wrong spherical harmonic called!", std::cout);
            }
            break;
        default:
        {
			vec d = cartesian_to_spherical(x, y, z);
            SH = real_spherical(l, m, d[1], d[2]);
        }
            
        }
        return SH;
    }

   // // https://www.wolframalpha.com/input?i=LegendreP%5B6%2C6%2Cx%5D watch the signs!
   // double associated_legendre_polynomial(const int &l, const int &m, const double &x) {
   //     switch (l) {
   //     case 0:
   //         return 1.0;
   //         break;
   //     case 1:
   //         switch (m) {
   //         case 0:
   //             return x;
   //             break;
   //         case 1:
   //             return std::sqrt(1 - x * x);
   //             break;
   //         case -1:
   //             return -0.5 * std::sqrt(1 - x * x);
   //             break;
   //         default:
   //             err_checkf(false, "This is impossible!", std::cout);
   //         }
   //         break;
   //     case 2:
   //         switch (m) {
   //         case 0:
   //             return 0.5 * (3 * x * x - 1);
   //             break;
   //         case 1:
   //             return 3.0 * x * std::sqrt(1 - x * x);
   //             break;
   //         case 2:
   //             return -3.0 * (x * x - 1);
   //             break;
   //         case -1:
   //             return -0.5 * x * std::sqrt(1 - x * x);
   //             break;
   //         case -2:
   //             return -0.125 * (1 - x * x);
   //             break;
   //         default:
   //             err_checkf(false, "This is impossible!", std::cout);
   //         }
   //         break;
   //     case 3:
   //         switch (m) {
   //         case 0:
   //             return 0.5 * x * (5 * x * x - 3);
   //             break;
   //         case 1:
   //             return 1.5 * (5 * x * x - 1) * std::sqrt(1 - x * x);
   //             break;
   //         case 2:
   //             return -15.0 * x * (x * x - 1);
   //             break;
   //         case 3:
   //             return 15.0 * pow(1 - x * x, 1.5);
   //             break;
   //         case -1:
   //             return -0.125 * (5 * x * x - 1) * std::sqrt(1 - x * x);
   //             break;
   //         case -2:
   //             return 0.125 * x * (x * x - 1);
   //             break;
   //         case -3:
   //             return -(1.0 / 48.0) * pow(1 - x * x, 1.5);
   //             break;
   //         default:
   //             err_checkf(false, "This is impossible!", std::cout);
   //         }
   //         break;
   //     case 4:
   //         switch (m) {
   //         case 0:
   //             return 4.375 * x * x * x * x - 3.75 * x * x + 0.375;
   //             break;
   //         case 1:
   //             return 2.5 * std::sqrt(1 - x * x) * (7 * x * x * x - 3 * x);
   //             break;
   //         case 2:
   //             return -7.5 * (x * x - 1) * (7 * x * x - 1);
   //             break;
   //         case 3:
   //             return 105.0 * x * pow(1 - x * x, 1.5);
   //             break;
   //         case 4:
   //             return 105.0 * pow(x * x - 1, 2.0);
   //             break;
   //         case -1:
   //             return -0.125 * std::sqrt(1 - x * x) * (7 * x * x * x - 3 * x);
   //             break;
   //         case -2:
   //             return (1.0 / 48.0) * (x * x - 1) * (7 * x * x - 1);
   //             break;
   //         case -3:
   //             return -(1.0 / 48.0) * pow(1 - x, 1.5) * x * pow(x + 1, 1.5);
   //             break;
   //         case -4:
   //             return -(1.0 / 384.0) * pow(x * x - 1, 2.0);
   //             break;
   //         default:
   //             err_checkf(false, "This is impossible!", std::cout);
   //         }
   //         break;
   //     case 5:
   //         switch (m) {
   //         case 0:
   //             return 7.875 * x * x * x * x * x - 8.75 * x * x * x + 1.875 * x;
   //             break;
   //         case 1:
   //             return (15.0 / 8.0) * std::sqrt(1 - x * x) * (21 * x * x * x * x - 14 * x * x + 1);
   //             break;
   //         case 2:
   //             return -(105.0 / 2.0) * (x * x - 1) * (3 * x * x * x - x);
   //             break;
   //         case 3:
   //             return (105.0 / 2.0) * pow(1 - x * x, 1.5) * (9 * x * x - 1);
   //             break;
   //         case 4:
   //             return 945.0 * x * pow(x * x - 1, 2);
   //             break;
   //         case 5:
   //             return 945.0 * pow(1 - x * x, 2.5);
   //             break;
   //         case -1:
   //             return -std::sqrt(1 - x * x) * (1.3125 * x * x * x * x - 0.875 * x * x + 0.0625);
   //             break;
   //         case -2:
   //             return (1.0 / 16.0) * (x * x - 1) * (3 * x * x * x - x);
   //             break;
   //         case -3:
   //             return -(1.0 / 384.0) * pow(1 - x * x, 1.5) * (9 * x * x - 1);
   //             break;
   //         case -4:
   //             return -(1.0 / 384.0) * x * pow(x * x - 1, 2.0);
   //             break;
   //         case -5:
   //             return -(1.0 / 3840.0) * pow(1 - x * x, 2.5);
   //             break;
   //         default:
   //             err_checkf(false, "This is impossible!", std::cout);
   //         }
   //         break;
   //     case 6: 
   //         switch (m) {
   //         case 0:
   //             return (1.0 / 16.0) * (231.0 * x * x * x * x * x * x - 315.0 * x * x * x * x + 105 * x * x - 5);
   //             break;
   //         case 1:
   //             return (21.0 / 8.0) * std::sqrt(1 - x * x) * (33 * x * x * x * x * x - 30 * x * x * x + 5 * x);
   //             break;
   //         case 2:
   //             return -(105.0 / 8.0) * (x * x - 1) * (33 * x * x * x * x - 18 * x * x + 1);
   //             break;
   //         case 3:
   //             return (315.0 / 2.0) * pow(1 - x * x, 1.5) * (11 * x * x * x - 3 * x);
   //             break;
   //         case 4:
   //             return (945.0 / 2.0) * pow(x * x - 1.0, 2.0) * (11 * x * x - 1);
   //             break;
   //         case 5:
   //             return 10395.0 * x * pow(1 - x * x, 2.5);
   //             break;
   //         case 6:
   //             return -10395.0 * pow(x * x - 1, 3.0);
   //             break;
   //         case -1:
   //             return -(1.0 / 16.0) * std::sqrt(1 - x * x) * (33 * x * x * x * x * x - 30 * x * x * x + 5 * x);
   //             break;
   //         case -2:
   //             return  (1.0 / 128.0) * (x * x - 1) * (33 * x * x * x * x - 18 * x * x + 1);
   //             break;
   //         case -3:
   //             return -(1.0 / 384.0) * pow(1 - x * x, 1.5) * (11 * x * x * x - 3 * x);
   //             break;
   //         case -4:
   //             return -(1.0 / 3840.0) * pow(x * x - 1, 2.0) * (11 * x * x - 1);
   //             break;
   //         case -5:
   //             return -(1.0 / 3840.0) * x * pow(1 - x * x, 2.5);
   //             break;
   //         case -6:
   //             return (1.0 / 46080.0) * pow(x * x - 1, 3.0);
   //             break;
   //         default:
   //             err_checkf(false, "This is impossible!", std::cout);
   //         }
   //         break;
   //     case 7: 
   //         switch (m)
   //         {
   //         case 0:
   //             return (1.0 / 16.0) * (429 * x * x * x * x * x * x * x - 693 * x * x * x * x * x + 315 * x * x * x - 35 * x);
   //             break;
   //         case 1:
   //             return (7.0 / 16.0) * std::sqrt(1 - x * x) * (429 * x * x * x * x * x * x - 495 * x * x * x * x + 135 * x * x - 5);
   //             break;
   //         case 2:
   //             return -(63.0 / 8.0) * (x * x - 1) * (143 * x * x * x * x * x - 110 * x * x * x + 15 * x);
   //             break;
   //         case 3:
   //             return (315.0 / 8.0) * pow(1 - x * x, 1.5) * (143 * x * x * x * x - 66 * x * x + 3);
   //             break;
   //         case 4:
   //             return (3465.0 / 2.0) * pow(x * x - 1, 2.0) * (13 * x * x * x - 3 * x);
   //             break;
   //         case 5:
   //             return (10395.0 / 2.0) * pow(1 - x * x, 2.5) * (13 * x * x - 1);
   //             break;
   //         case 6:
   //             return -135135.0 * pow(x * x - 1, 3.0) * x;
   //             break;
   //         case 7:
   //             return 135135.0 * pow(1 - x * x, 3.5);
   //             break;
   //         case -1:
   //             return -(1.0 / 128.0) * std::sqrt(1 - x * x) * (429 * x * x * x * x * x * x - 495 * x * x * x * x + 135 * x * x - 5);
   //             break;
   //         case -2:
   //             return (1.0 / 384.0) * (x * x - 1) * (143 * x * x * x * x * x - 110 * x * x * x + 15 * x);
   //             break;
   //         case -3:
   //             return -(1.0 / 3840.0) * pow(1 - x * x, 1.5) * (143 * x * x * x * x - 66 * x * x + 3);
   //             break;
   //         case -4:
   //             return -(1.0 / 3840.0) * pow(x * x - 1, 2.0) * (13 * x * x * x - 3 * x);
   //             break;
   //         case -5:
   //             return -(1.0 / 46080.0) * pow(1 - x * x, 2.5) * (13 * x * x - 1);
   //             break;
   //         case -6:
   //             return (1.0 / 46080.0) * pow(x * x - 1.0, 3.0) * x;
   //             break;
   //         case -7:
   //             return -(1.0 / 645120.0) * pow(1 - x * x, 3.5);
   //             break;
   //         default:
   //             err_checkf(false, "This is impossible!", std::cout);
   //         }
   //         break;
   //     case 8:
			//switch (m)
			//{
			//case 0:
			//	return (1.0 / 128.0) * (6435 * x * x * x * x * x * x * x * x - 12012 * x * x * x * x * x * x + 6930 * x * x * x * x - 1260 * x * x + 35);
			//case 1:
			//	return (9.0 / 16.0) * std::sqrt(1 - x * x) * (715 * x * x * x * x * x * x * x - 1001 * x * x * x * x * x + 385 * x * x * x - 35 * x);
			//case 2:
			//	return -(315.0 / 16.0) * (x * x - 1) * (143 * x * x * x * x * x * x - 143 * x * x * x * x + 33 * x * x - 1);
			//case 3:
			//	return (3465.0 / 8.0) * pow(1 - x * x, 1.5) * (39 * x * x * x * x * x - 26 * x * x * x + 3 * x);
			//case 4:
			//	return (10395.0 / 8.0) * pow(x * x - 1, 2.0) * (65 * x * x * x * x - 26 * x * x + 1);
			//case 5:
			//	return (135135.0 / 2.0) * pow(1 - x * x, 2.5) * (5 * x * x * x - x);
			//case 6:
			//	return -(135135.0 / 2.0) * pow(x * x - 1, 3.0) * (15 * x * x - 1);
			//case 7:
			//	return 2027025.0 * x * pow(1 - x * x, 3.5);
			//case 8:
			//	return 2027025.0 * pow(x * x - 1, 4.0);
			//case -1:
			//	return -(1.0 / 128.0) * std::sqrt(1 - x * x) * (715 * x * x * x * x * x * x * x - 1001 * x * x * x * x * x + 385 * x * x * x - 35 * x);
			//case -2:
			//	return (1.0 / 256.0) * (x * x - 1) * (143 * x * x * x * x * x * x - 143 * x * x * x * x + 33 * x * x - 1);
			//case -3:
			//	return -(1.0 / 768.0) * pow(1 - x * x, 1.5) * (39 * x * x * x * x * x - 26 * x * x * x + 3 * x);
			//case -4:
			//	return -(1.0 / 15360.0) * pow(x * x - 1, 2.0) * (65 * x * x * x * x - 26 * x * x + 1);
			//case -5:
			//	return -(1.0 / 15360.0) * pow(1 - x * x, 2.5) * (5 * x * x * x - x);
			//case -6:
			//	return (1.0 / 645120.0) * pow(x * x - 1, 3.0) * (15 * x * x - 1);
			//case -7:
			//	return -(1.0 / 645120.0) * x * pow(1 - x * x, 3.5);
			//case -8:
			//	return -(1.0 / 10321920.0) * pow(x * x - 1, 4.0);
   //         default:
   //             err_checkf(false, "This is impossible!", std::cout);
			//}
   //         break;
   //     default:
   //         return std::assoc_legendre(l, m, x) * ((m < 0) ? ((m % 2 == 0) ? -1 : 1) : 1);
   //     }
   // };

    //Returns the spherical coordinates of a given cartesian vector
    //the output is a vector with the (radius, theta and phi)
    vec cartesian_to_spherical(const double& x, const double& y, const double& z) {
        double r = std::sqrt(x * x + y * y + z * z);
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