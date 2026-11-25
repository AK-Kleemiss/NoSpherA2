#pragma once
#include "convenience.h"
#include "atoms.h"

namespace
{
	//Per-thread scratch: K components per function summed over its primitives, then one multiply into
	//the MOs; only functions the exponent cutoff kept are visited and reset
	template <int K>
	struct ao_scratch {
		vec v;
		std::vector<char> hit;
		ivec touched;
		void begin(const size_t nao)
		{
			if (hit.size() < nao) {
				hit.resize(nao, 0);
				v.resize(K * nao);
			}
			touched.clear();
		}
		void add(const int a, const double s, const double* chi)
		{
			double* t = v.data() + (size_t)K * a;
			if (!hit[a]) {
				hit[a] = 1;
				touched.push_back(a);
				for (int k = 0; k < K; k++) t[k] = s * chi[k];
			}
			else
				for (int k = 0; k < K; k++) t[k] += s * chi[k];
		}
		//phi[k * nmo + mo] += sum_A v[A][k] C[A][mo] over the occupied MOs only: every caller weighs phi by the occupation,
		//so the virtuals stay at the caller's zero. The occ != 0 runs are rescanned per call, O(nmo), so a changed occupation
		//never goes stale; the per-MO sum order is unchanged, so the occupied values are bit-identical to the full loop
		ivec runs;
		template <class Occ>
		void to_mo(const double* coef_ao, const int nmo, double* phi, Occ occ)
		{
			runs.clear();
			for (int mo = 0; mo < nmo; mo++)
				if (occ(mo) != 0.0) {
					if (runs.empty() || runs.back() != mo) { runs.push_back(mo); runs.push_back(mo + 1); }
					else runs.back() = mo + 1;
				}
			//Eight functions per pass over phi, added one after the other: every element sees the products in the
			//order the one-function loop adds them, with an eighth of the loads and stores of phi. The coefficients
			//of an MO are loaded once for all K components, not once per component
			const int nt = static_cast<int>(touched.size());
			int i = 0;
			for (; i + 8 <= nt; i += 8) {
				const double* c[8];
				double s[K][8];
				for (int j = 0; j < 8; j++) {
					c[j] = coef_ao + (size_t)touched[i + j] * nmo;
					const double* t = v.data() + (size_t)K * touched[i + j];
					for (int k = 0; k < K; k++) s[k][j] = t[k];
				}
				for (size_t r = 0; r < runs.size(); r += 2)
#pragma omp simd
					for (int mo = runs[r]; mo < runs[r + 1]; mo++) {
						const double c0 = c[0][mo], c1 = c[1][mo], c2 = c[2][mo], c3 = c[3][mo], c4 = c[4][mo], c5 = c[5][mo], c6 = c[6][mo], c7 = c[7][mo];
						for (int k = 0; k < K; k++) {
							double x = phi[(size_t)k * nmo + mo];
							x += s[k][0] * c0;
							x += s[k][1] * c1;
							x += s[k][2] * c2;
							x += s[k][3] * c3;
							x += s[k][4] * c4;
							x += s[k][5] * c5;
							x += s[k][6] * c6;
							x += s[k][7] * c7;
							phi[(size_t)k * nmo + mo] = x;
						}
					}
			}
			for (; i < nt; i++) {
				const double* c = coef_ao + (size_t)touched[i] * nmo;
				const double* t = v.data() + (size_t)K * touched[i];
				for (size_t r = 0; r < runs.size(); r += 2)
#pragma omp simd
					for (int mo = runs[r]; mo < runs[r + 1]; mo++) {
						const double cm = c[mo];
						for (int k = 0; k < K; k++) phi[(size_t)k * nmo + mo] += t[k] * cm;
					}
			}
			for (const int a : touched) hit[a] = 0;
		}
	};
	//The Cartesian monomial x^l y^m z^n of a primitive type, from the powers table d_ of compute_dens_cartesian:
	//0-2 dx dy dz, 3 r^2, 4-6 squares, 7-9 cubes, 10-12 fourth, 13-15 fifth powers
	inline double cart_monomial(const int type, const double* d_)
	{
		switch (type)
		{
		case 0: return 1.0;
		case 1: return 1.0;// 0, 0, 0,
		case 2: return d_[0];// 1, 0, 0,
		case 3: return d_[1];// 0, 1, 0,
		case 4: return d_[2];// 0, 0, 1,
		case 5: return d_[4];// 2, 0, 0,
		case 6: return d_[5];// 0, 2, 0,
		case 7: return d_[6];// 0, 0, 2,
		case 8: return d_[0] * d_[1];// 1, 1, 0,
		case 9: return d_[0] * d_[2];// 1, 0, 1,
		case 10: return d_[1] * d_[2];// 0, 1, 1,
		case 11: return d_[7];// 3, 0, 0,
		case 12: return d_[8];// 0, 3, 0,
		case 13: return d_[9];// 0, 0, 3,
		case 14: return d_[4] * d_[1];// 2, 1, 0,
		case 15: return d_[4] * d_[2];// 2, 0, 1,
		case 16: return d_[5] * d_[2];// 0, 2, 1,
		case 17: return d_[0] * d_[5];// 1, 2, 0,
		case 18: return d_[0] * d_[6];// 1, 0, 2,
		case 19: return d_[1] * d_[6];// 0, 1, 2,
		case 20: return d_[0] * d_[1] * d_[2];// 1, 1, 1,
		case 21: return d_[12];// 0, 0, 4,
		case 22: return d_[1] * d_[9];// 0, 1, 3,
		case 23: return d_[5] * d_[6];// 0, 2, 2,
		case 24: return d_[8] * d_[2];// 0, 3, 1,
		case 25: return d_[11];// 0, 4, 0,
		case 26: return d_[0] * d_[9];// 1, 0, 3,
		case 27: return d_[0] * d_[1] * d_[6];// 1, 1, 2,
		case 28: return d_[0] * d_[5] * d_[2];// 1, 2, 1,
		case 29: return d_[0] * d_[8];// 1, 3, 0,
		case 30: return d_[4] * d_[6];// 2, 0, 2,
		case 31: return d_[4] * d_[1] * d_[2];// 2, 1, 1,
		case 32: return d_[4] * d_[5];// 2, 2, 0,
		case 33: return d_[7] * d_[2];// 3, 0, 1,
		case 34: return d_[7] * d_[1];// 3, 1, 0,
		case 35: return d_[10];// 4, 0, 0,
		case 36: return d_[15];// 0, 0, 5,
		case 37: return d_[1] * d_[12];// 0, 1, 4,
		case 38: return d_[5] * d_[9];// 0, 2, 3,
		case 39: return d_[8] * d_[6];// 0, 3, 2,
		case 40: return d_[11] * d_[2];// 0, 4, 1,
		case 41: return d_[14];// 0, 5, 0,
		case 42: return d_[0] * d_[12];// 1, 0, 4,
		case 43: return d_[0] * d_[1] * d_[9];// 1, 1, 3,
		case 44: return d_[0] * d_[5] * d_[6];// 1, 2, 2,
		case 45: return d_[0] * d_[8] * d_[2];// 1, 3, 1,
		case 46: return d_[0] * d_[11];// 1, 4, 0,
		case 47: return d_[4] * d_[9];// 2, 0, 3,
		case 48: return d_[4] * d_[1] * d_[6];// 2, 1, 2,
		case 49: return d_[4] * d_[5] * d_[2];// 2, 2, 1,
		case 50: return d_[4] * d_[8];// 2, 3, 0,
		case 51: return d_[7] * d_[6];// 3, 0, 2,
		case 52: return d_[7] * d_[1] * d_[2];// 3, 1, 1,
		case 53: return d_[7] * d_[5];// 3, 2, 0,
		case 54: return d_[10] * d_[2];// 4, 0, 1,
		case 55: return d_[10] * d_[1];// 4, 1, 0,
		case 56: return d_[13];// 5, 0, 0,
		case 57: return d_[15] * d_[2];// 0, 0, 6,
		case 58: return d_[1] * d_[15];// 0, 1, 5,
		case 59: return d_[5] * d_[12];// 0, 2, 4,
		case 60: return d_[8] * d_[9];// 0, 3, 3,
		case 61: return d_[11] * d_[6];// 0, 4, 2,
		case 62: return d_[14] * d_[2];// 0, 5, 1,
		case 63: return d_[14] * d_[1];// 0, 6, 0,
		case 64: return d_[0] * d_[15];// 1, 0, 5,
		case 65: return d_[0] * d_[1] * d_[12];// 1, 1, 4,
		case 66: return d_[0] * d_[5] * d_[9];// 1, 2, 3,
		case 67: return d_[0] * d_[8] * d_[6];// 1, 3, 2,
		case 68: return d_[0] * d_[11] * d_[2];// 1, 4, 1,
		case 69: return d_[0] * d_[14];// 1, 5, 0,
		case 70: return d_[4] * d_[12];// 2, 0, 4,
		case 71: return d_[4] * d_[1] * d_[9];// 2, 1, 3,
		case 72: return d_[4] * d_[5] * d_[6];// 2, 2, 2,
		case 73: return d_[4] * d_[8] * d_[2];// 2, 3, 1,
		case 74: return d_[4] * d_[11];// 2, 4, 0,
		case 75: return d_[7] * d_[9];// 3, 0, 3,
		case 76: return d_[7] * d_[1] * d_[6];// 3, 1, 2,
		case 77: return d_[7] * d_[5] * d_[2];// 3, 2, 1,
		case 78: return d_[7] * d_[8];// 3, 3, 0,
		case 79: return d_[10] * d_[6];// 4, 0, 2,
		case 80: return d_[10] * d_[1] * d_[2];// 4, 1, 1,
		case 81: return d_[10] * d_[5];// 4, 2, 0,
		case 82: return d_[13] * d_[2];// 5, 0, 1,
		case 83: return d_[13] * d_[1];// 5, 1, 0,
		case 84: return d_[13] * d_[0];// 6, 0, 0,
		case 85: return d_[15] * d_[6];// 0, 0, 7,
		case 86: return d_[1] * d_[15] * d_[2];// 0, 1, 6,
		case 87: return d_[5] * d_[15];// 0, 2, 5,
		case 88: return d_[8] * d_[12];// 0, 3, 4,
		case 89: return d_[11] * d_[9];// 0, 4, 3,
		case 90: return d_[14] * d_[6];// 0, 5, 2,
		case 91: return d_[14] * d_[1] * d_[2];// 0, 6, 1,
		case 92: return d_[14] * d_[5];// 0, 7, 0,
		case 93: return d_[0] * d_[15] * d_[2];// 1, 0, 6,
		case 94: return d_[0] * d_[1] * d_[15];// 1, 1, 5,
		case 95: return d_[0] * d_[5] * d_[12];// 1, 2, 4,
		case 96: return d_[0] * d_[8] * d_[9];// 1, 3, 3,
		case 97: return d_[0] * d_[11] * d_[6];// 1, 4, 2,
		case 98: return d_[0] * d_[14] * d_[2];// 1, 5, 1,
		case 99: return d_[0] * d_[14] * d_[1];// 1, 6, 0,
		case 100: return d_[4] * d_[15];// 2, 0, 5,
		case 101: return d_[4] * d_[1] * d_[12];// 2, 1, 4,
		case 102: return d_[4] * d_[5] * d_[9];// 2, 2, 3,
		case 103: return d_[4] * d_[8] * d_[6];// 2, 3, 2,
		case 104: return d_[4] * d_[11] * d_[2];// 2, 4, 1,
		case 105: return d_[4] * d_[14];// 2, 5, 0,
		case 106: return d_[7] * d_[12];// 3, 0, 4,
		case 107: return d_[7] * d_[1] * d_[9];// 3, 1, 3,
		case 108: return d_[7] * d_[5] * d_[6];// 3, 2, 2,
		case 109: return d_[7] * d_[8] * d_[2];// 3, 3, 1,
		case 110: return d_[7] * d_[11];// 3, 4, 0,
		case 111: return d_[10] * d_[9];// 4, 0, 3,
		case 112: return d_[10] * d_[1] * d_[6];// 4, 1, 2,
		case 113: return d_[10] * d_[5] * d_[2];// 4, 2, 1,
		case 114: return d_[10] * d_[8];// 4, 3, 0,
		case 115: return d_[13] * d_[6];// 5, 0, 2,
		case 116: return d_[13] * d_[1] * d_[2];// 5, 1, 1,
		case 117: return d_[13] * d_[5];// 5, 2, 0,
		case 118: return d_[13] * d_[0] * d_[2];// 6, 0, 1,
		case 119: return d_[13] * d_[0] * d_[1];// 6, 1, 0,
		case 120: return d_[13] * d_[4];// 7, 0, 0,
		case 121: return d_[15] * d_[9];// 0, 0, 8,
		case 122: return d_[1] * d_[15] * d_[6];// 0, 1, 7,
		case 123: return d_[5] * d_[15] * d_[2];// 0, 2, 6,
		case 124: return d_[8] * d_[15];// 0, 3, 5,
		case 125: return d_[11] * d_[12];// 0, 4, 4,
		case 126: return d_[14] * d_[9];// 0, 5, 3,
		case 127: return d_[14] * d_[1] * d_[6];// 0, 6, 2,
		case 128: return d_[14] * d_[5] * d_[2];// 0, 7, 1,
		case 129: return d_[14] * d_[8];// 0, 8, 0,
		case 130: return d_[0] * d_[15] * d_[6];// 1, 0, 7,
		case 131: return d_[0] * d_[1] * d_[15] * d_[2];// 1, 1, 6,
		case 132: return d_[0] * d_[5] * d_[15];// 1, 2, 5,
		case 133: return d_[0] * d_[8] * d_[12];// 1, 3, 4,
		case 134: return d_[0] * d_[11] * d_[9];// 1, 4, 3,
		case 135: return d_[0] * d_[14] * d_[6];// 1, 5, 2,
		case 136: return d_[0] * d_[14] * d_[1] * d_[2];// 1, 6, 1,
		case 137: return d_[0] * d_[14] * d_[5];// 1, 7, 0,
		case 138: return d_[4] * d_[15] * d_[2];// 2, 0, 6,
		case 139: return d_[4] * d_[1] * d_[15];// 2, 1, 5,
		case 140: return d_[4] * d_[5] * d_[12];// 2, 2, 4,
		case 141: return d_[4] * d_[8] * d_[9];// 2, 3, 3,
		case 142: return d_[4] * d_[11] * d_[6];// 2, 4, 2,
		case 143: return d_[4] * d_[14] * d_[2];// 2, 5, 1,
		case 144: return d_[4] * d_[14] * d_[1];// 2, 6, 0,
		case 145: return d_[7] * d_[15];// 3, 0, 5,
		case 146: return d_[7] * d_[1] * d_[12];// 3, 1, 4,
		case 147: return d_[7] * d_[5] * d_[9];// 3, 2, 3,
		case 148: return d_[7] * d_[8] * d_[6];// 3, 3, 2,
		case 149: return d_[7] * d_[11] * d_[2];// 3, 4, 1,
		case 150: return d_[7] * d_[14];// 3, 5, 0,
		case 151: return d_[10] * d_[12];// 4, 0, 4,
		case 152: return d_[10] * d_[1] * d_[9];// 4, 1, 3,
		case 153: return d_[10] * d_[5] * d_[6];// 4, 2, 2,
		case 154: return d_[10] * d_[8] * d_[2];// 4, 3, 1,
		case 155: return d_[10] * d_[11];// 4, 4, 0,
		case 156: return d_[13] * d_[9];// 5, 0, 3,
		case 157: return d_[13] * d_[1] * d_[6];// 5, 1, 2,
		case 158: return d_[13] * d_[5] * d_[2];// 5, 2, 1,
		case 159: return d_[13] * d_[8];// 5, 3, 0,
		case 160: return d_[13] * d_[0] * d_[6];// 6, 0, 2,
		case 161: return d_[13] * d_[0] * d_[1] * d_[2];// 6, 1, 1,
		case 162: return d_[13] * d_[0] * d_[5];// 6, 2, 0,
		case 163: return d_[13] * d_[4] * d_[2];// 7, 0, 1,
		case 164: return d_[13] * d_[4] * d_[1];// 7, 1, 0,
		case 165: return d_[13] * d_[7];// 8, 0, 0,
		case 166: return d_[15] * d_[12];// 0, 0, 9,
		case 167: return d_[1] * d_[15] * d_[9];// 0, 1, 8,
		case 168: return d_[5] * d_[15] * d_[6];// 0, 2, 7,
		case 169: return d_[8] * d_[15] * d_[2];// 0, 3, 6,
		case 170: return d_[11] * d_[15];// 0, 4, 5,
		case 171: return d_[14] * d_[12];// 0, 5, 4,
		case 172: return d_[14] * d_[1] * d_[9];// 0, 6, 3,
		case 173: return d_[14] * d_[5] * d_[6];// 0, 7, 2,
		case 174: return d_[14] * d_[8] * d_[2];// 0, 8, 1,
		case 175: return d_[14] * d_[11];// 0, 9, 0,
		case 176: return d_[0] * d_[15] * d_[9];// 1, 0, 8,
		case 177: return d_[0] * d_[1] * d_[15] * d_[6];// 1, 1, 7,
		case 178: return d_[0] * d_[5] * d_[15] * d_[2];// 1, 2, 6,
		case 179: return d_[0] * d_[8] * d_[15];// 1, 3, 5,
		case 180: return d_[0] * d_[11] * d_[12];// 1, 4, 4,
		case 181: return d_[0] * d_[14] * d_[9];// 1, 5, 3,
		case 182: return d_[0] * d_[14] * d_[1] * d_[6];// 1, 6, 2,
		case 183: return d_[0] * d_[14] * d_[5] * d_[2];// 1, 7, 1,
		case 184: return d_[0] * d_[14] * d_[8];// 1, 8, 0,
		case 185: return d_[4] * d_[15] * d_[6];// 2, 0, 7,
		case 186: return d_[4] * d_[1] * d_[15] * d_[2];// 2, 1, 6,
		case 187: return d_[4] * d_[5] * d_[15];// 2, 2, 5,
		case 188: return d_[4] * d_[8] * d_[12];// 2, 3, 4,
		case 189: return d_[4] * d_[11] * d_[9];// 2, 4, 3,
		case 190: return d_[4] * d_[14] * d_[6];// 2, 5, 2,
		case 191: return d_[4] * d_[14] * d_[1] * d_[2];// 2, 6, 1,
		case 192: return d_[4] * d_[14] * d_[5];// 2, 7, 0,
		case 193: return d_[7] * d_[15] * d_[2];// 3, 0, 6,
		case 194: return d_[7] * d_[1] * d_[15];// 3, 1, 5,
		case 195: return d_[7] * d_[5] * d_[12];// 3, 2, 4,
		case 196: return d_[7] * d_[8] * d_[9];// 3, 3, 3,
		case 197: return d_[7] * d_[11] * d_[6];// 3, 4, 2,
		case 198: return d_[7] * d_[14] * d_[2];// 3, 5, 1,
		case 199: return d_[7] * d_[14] * d_[1];// 3, 6, 0,
		case 200: return d_[10] * d_[15];// 4, 0, 5,
		case 201: return d_[10] * d_[1] * d_[12];// 4, 1, 4,
		case 202: return d_[10] * d_[5] * d_[9];// 4, 2, 3,
		case 203: return d_[10] * d_[8] * d_[6];// 4, 3, 2,
		case 204: return d_[10] * d_[11] * d_[2];// 4, 4, 1,
		case 205: return d_[10] * d_[14];// 4, 5, 0,
		case 206: return d_[13] * d_[12];// 5, 0, 4,
		case 207: return d_[13] * d_[1] * d_[9];// 5, 1, 3,
		case 208: return d_[13] * d_[5] * d_[6];// 5, 2, 2,
		case 209: return d_[13] * d_[8] * d_[2];// 5, 3, 1,
		case 210: return d_[13] * d_[11];// 5, 4, 0,
		case 211: return d_[13] * d_[0] * d_[9];// 6, 0, 3,
		case 212: return d_[13] * d_[0] * d_[1] * d_[6];// 6, 1, 2,
		case 213: return d_[13] * d_[0] * d_[5] * d_[2];// 6, 2, 1,
		case 214: return d_[13] * d_[0] * d_[8];// 6, 3, 0,
		case 215: return d_[13] * d_[4] * d_[6];// 7, 0, 2,
		case 216: return d_[13] * d_[4] * d_[1] * d_[2];// 7, 1, 1,
		case 217: return d_[13] * d_[4] * d_[5];// 7, 2, 0,
		case 218: return d_[13] * d_[7] * d_[2];// 8, 0, 1,
		case 219: return d_[13] * d_[7] * d_[1];// 8, 1, 0,
		case 220: return d_[13] * d_[10];// 9, 0, 0,
		case 221: return d_[15] * d_[15];// 0, 0, 10,
		case 222: return d_[1] * d_[15] * d_[12];// 0, 1, 9,
		case 223: return d_[5] * d_[15] * d_[9];// 0, 2, 8,
		case 224: return d_[8] * d_[15] * d_[6];// 0, 3, 7,
		case 225: return d_[11] * d_[15] * d_[2];// 0, 4, 6,
		case 226: return d_[14] * d_[15];// 0, 5, 5,
		case 227: return d_[14] * d_[1] * d_[12];// 0, 6, 4,
		case 228: return d_[14] * d_[5] * d_[9];// 0, 7, 3,
		case 229: return d_[14] * d_[8] * d_[6];// 0, 8, 2,
		case 230: return d_[14] * d_[11] * d_[2];// 0, 9, 1,
		case 231: return d_[14] * d_[14];// 0, 10, 0,
		case 232: return d_[0] * d_[15] * d_[12];// 1, 0, 9,
		case 233: return d_[0] * d_[1] * d_[15] * d_[9];// 1, 1, 8,
		case 234: return d_[0] * d_[5] * d_[15] * d_[6];// 1, 2, 7,
		case 235: return d_[0] * d_[8] * d_[15] * d_[2];// 1, 3, 6,
		case 236: return d_[0] * d_[11] * d_[15];// 1, 4, 5,
		case 237: return d_[0] * d_[14] * d_[12];// 1, 5, 4,
		case 238: return d_[0] * d_[14] * d_[1] * d_[9];// 1, 6, 3,
		case 239: return d_[0] * d_[14] * d_[5] * d_[6];// 1, 7, 2,
		case 240: return d_[0] * d_[14] * d_[8] * d_[2];// 1, 8, 1,
		case 241: return d_[0] * d_[14] * d_[11];// 1, 9, 0,
		case 242: return d_[4] * d_[15] * d_[9];// 2, 0, 8,
		case 243: return d_[4] * d_[1] * d_[15] * d_[6];// 2, 1, 7,
		case 244: return d_[4] * d_[5] * d_[15] * d_[2];// 2, 2, 6,
		case 245: return d_[4] * d_[8] * d_[15];// 2, 3, 5,
		case 246: return d_[4] * d_[11] * d_[12];// 2, 4, 4,
		case 247: return d_[4] * d_[14] * d_[9];// 2, 5, 3,
		case 248: return d_[4] * d_[14] * d_[1] * d_[6];// 2, 6, 2,
		case 249: return d_[4] * d_[14] * d_[5] * d_[2];// 2, 7, 1,
		case 250: return d_[4] * d_[14] * d_[8];// 2, 8, 0,
		case 251: return d_[7] * d_[15] * d_[6];// 3, 0, 7,
		case 252: return d_[7] * d_[1] * d_[15] * d_[2];// 3, 1, 6,
		case 253: return d_[7] * d_[5] * d_[15];// 3, 2, 5,
		case 254: return d_[7] * d_[8] * d_[12];// 3, 3, 4,
		case 255: return d_[7] * d_[11] * d_[9];// 3, 4, 3,
		case 256: return d_[7] * d_[14] * d_[6];// 3, 5, 2,
		case 257: return d_[7] * d_[14] * d_[1] * d_[2];// 3, 6, 1,
		case 258: return d_[7] * d_[14] * d_[5];// 3, 7, 0,
		case 259: return d_[10] * d_[15] * d_[2];// 4, 0, 6,
		case 260: return d_[10] * d_[1] * d_[15];// 4, 1, 5,
		case 261: return d_[10] * d_[5] * d_[12];// 4, 2, 4,
		case 262: return d_[10] * d_[8] * d_[9];// 4, 3, 3,
		case 263: return d_[10] * d_[11] * d_[6];// 4, 4, 2,
		case 264: return d_[10] * d_[14] * d_[2];// 4, 5, 1,
		case 265: return d_[10] * d_[14] * d_[1];// 4, 6, 0,
		case 266: return d_[13] * d_[15];// 5, 0, 5,
		case 267: return d_[13] * d_[1] * d_[12];// 5, 1, 4,
		case 268: return d_[13] * d_[5] * d_[9];// 5, 2, 3,
		case 269: return d_[13] * d_[8] * d_[6];// 5, 3, 2,
		case 270: return d_[13] * d_[11] * d_[2];// 5, 4, 1,
		case 271: return d_[13] * d_[14];// 5, 5, 0,
		case 272: return d_[13] * d_[0] * d_[12];// 6, 0, 4,
		case 273: return d_[13] * d_[0] * d_[1] * d_[9];// 6, 1, 3,
		case 274: return d_[13] * d_[0] * d_[5] * d_[6];// 6, 2, 2,
		case 275: return d_[13] * d_[0] * d_[8] * d_[2];// 6, 3, 1,
		case 276: return d_[13] * d_[0] * d_[11];// 6, 4, 0,
		case 277: return d_[13] * d_[4] * d_[9];// 7, 0, 3,
		case 278: return d_[13] * d_[4] * d_[1] * d_[6];// 7, 1, 2,
		case 279: return d_[13] * d_[4] * d_[5] * d_[2];// 7, 2, 1,
		case 280: return d_[13] * d_[4] * d_[8];// 7, 3, 0,
		case 281: return d_[13] * d_[7] * d_[6];// 8, 0, 2,
		case 282: return d_[13] * d_[7] * d_[1] * d_[2];// 8, 1, 1,
		case 283: return d_[13] * d_[7] * d_[5];// 8, 2, 0,
		case 284: return d_[13] * d_[10] * d_[2];// 9, 0, 1,
		case 285: return d_[13] * d_[10] * d_[1];// 9, 1, 0,
		case 286: return d_[13] * d_[13];// 10, 0, 0,
		default: return 1.0;
		}
	}
	//The table cart_monomial reads, for the point Pos relative to the centre c
	inline void cart_powers(double* d_, const d3& Pos, const atom& c)
	{
		d_[0] = Pos[0] - c.get_coordinate(0);
		d_[1] = Pos[1] - c.get_coordinate(1);
		d_[2] = Pos[2] - c.get_coordinate(2);
		d_[4] = d_[0] * d_[0];
		d_[5] = d_[1] * d_[1];
		d_[6] = d_[2] * d_[2];
		d_[3] = d_[4] + d_[5] + d_[6];
		d_[7] = d_[0] * d_[4];
		d_[8] = d_[1] * d_[5];
		d_[9] = d_[2] * d_[6];
		d_[10] = d_[0] * d_[7];
		d_[11] = d_[1] * d_[8];
		d_[12] = d_[2] * d_[9];
		d_[13] = d_[0] * d_[10];
		d_[14] = d_[1] * d_[11];
		d_[15] = d_[2] * d_[12];
	}
}
