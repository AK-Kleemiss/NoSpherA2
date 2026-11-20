#include "pch.h"

#include "core/convenience.h"
#include "core/constants.h"
#include "core/atoms.h"
#include "core/wfn_class.h"
#include "core/mo_class.h"

#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <random>

namespace
{
//Verbatim copies of master d624d8c4 (Src/core/wfn_density.cpp) with the members they read bound to locals
double compute_g_cartesian_master_ref(
	WFN& w,
	const d3 &Pos,
	vec2 &d,
	vec &phi)
{
	const int ncen = w.get_ncen(), nex = w.get_nex(), nmo = w.get_nmo();
	const std::vector<atom>& atoms = *w.get_atoms_ptr();
	const ivec& centers = w.get_centers();
	const ivec& types = w.get_types();
	const vec& exponents = w.get_exponents();
	std::vector<MO> MOs;
	for (int m = 0; m < nmo; m++) MOs.push_back(w.get_MO(m));
	auto get_coef_primitive_major = [&w]() { return w.get_coef_primitive_major(); };
	std::fill(phi.begin(), phi.end(), 0.0);
	double g = 0.0;
	int j;
	double ex, *d_;

	// precalculate some distances and powers of distances for faster computation
	for (j = 0; j < ncen; j++)
	{
		const atom &a = atoms[j];
		d_ = d[j].data();
		d_[0] = Pos[0] - a.get_coordinate(0);
		d_[1] = Pos[1] - a.get_coordinate(1);
		d_[2] = Pos[2] - a.get_coordinate(2);
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

	// Pre-cache frequently accessed data
	const int *centers_data = centers.data();
	const int *types_data = types.data();
	const double *exponents_data = exponents.data();
	double *phi_data = phi.data();
	const double exp_cutoff = constants::exp_cutoff;

	for (j = 0; j < nex; j++)
	{
		d_ = d[centers_data[j] - 1].data();
		ex = -exponents_data[j] * d_[3];
		if (ex < exp_cutoff)
		{ // corresponds to cutoff of maximum density contribution of 1E-5
			continue;
		}
		ex = exp(ex);
		switch (types_data[j])
		{
		case 0:  break;
		case 1:  break;// 0, 0, 0,
		case 2:  ex *= d_[0]; break;// 1, 0, 0,
		case 3:  ex *= d_[1]; break;// 0, 1, 0,
		case 4:  ex *= d_[2]; break;// 0, 0, 1,
		case 5:  ex *= d_[4]; break;// 2, 0, 0,
		case 6:  ex *= d_[5]; break;// 0, 2, 0,
		case 7:  ex *= d_[6]; break;// 0, 0, 2,
		case 8:  ex *= d_[0] * d_[1]; break;// 1, 1, 0,
		case 9:  ex *= d_[0] * d_[2]; break;// 1, 0, 1,
		case 10: ex *= d_[1] * d_[2]; break;// 0, 1, 1,
		case 11: ex *= d_[7]; break;// 3, 0, 0,
		case 12: ex *= d_[8]; break;// 0, 3, 0,
		case 13: ex *= d_[9]; break;// 0, 0, 3,
		case 14: ex *= d_[4] * d_[1]; break;// 2, 1, 0,
		case 15: ex *= d_[4] * d_[2]; break;// 2, 0, 1,
		case 16: ex *= d_[5] * d_[2]; break;// 0, 2, 1,
		case 17: ex *= d_[0] * d_[5]; break;// 1, 2, 0,
		case 18: ex *= d_[0] * d_[6]; break;// 1, 0, 2,
		case 19: ex *= d_[1] * d_[6]; break;// 0, 1, 2,
		case 20: ex *= d_[0] * d_[1] * d_[2]; break;// 1, 1, 1,
		case 21: ex *= d_[12]; break;// 0, 0, 4,
		case 22: ex *= d_[1] * d_[9]; break;// 0, 1, 3,
		case 23: ex *= d_[5] * d_[6]; break;// 0, 2, 2,
		case 24: ex *= d_[8] * d_[2]; break;// 0, 3, 1,
		case 25: ex *= d_[11]; break;// 0, 4, 0,
		case 26: ex *= d_[0] * d_[9]; break;// 1, 0, 3,
		case 27: ex *= d_[0] * d_[1] * d_[6]; break;// 1, 1, 2,
		case 28: ex *= d_[0] * d_[5] * d_[2]; break;// 1, 2, 1,
		case 29: ex *= d_[0] * d_[8]; break;// 1, 3, 0,
		case 30: ex *= d_[4] * d_[6]; break;// 2, 0, 2,
		case 31: ex *= d_[4] * d_[1] * d_[2]; break;// 2, 1, 1,
		case 32: ex *= d_[4] * d_[5]; break;// 2, 2, 0,
		case 33: ex *= d_[7] * d_[2]; break;// 3, 0, 1,
		case 34: ex *= d_[7] * d_[1]; break;// 3, 1, 0,
		case 35: ex *= d_[10]; break;// 4, 0, 0,
		case 36: ex *= d_[15]; break;// 0, 0, 5,
		case 37: ex *= d_[1] * d_[12]; break;// 0, 1, 4,
		case 38: ex *= d_[5] * d_[9]; break;// 0, 2, 3,
		case 39: ex *= d_[8] * d_[6]; break;// 0, 3, 2,
		case 40: ex *= d_[11] * d_[2]; break;// 0, 4, 1,
		case 41: ex *= d_[14]; break;// 0, 5, 0,
		case 42: ex *= d_[0] * d_[12]; break;// 1, 0, 4,
		case 43: ex *= d_[0] * d_[1] * d_[9]; break;// 1, 1, 3,
		case 44: ex *= d_[0] * d_[5] * d_[6]; break;// 1, 2, 2,
		case 45: ex *= d_[0] * d_[8] * d_[2]; break;// 1, 3, 1,
		case 46: ex *= d_[0] * d_[11]; break;// 1, 4, 0,
		case 47: ex *= d_[4] * d_[9]; break;// 2, 0, 3,
		case 48: ex *= d_[4] * d_[1] * d_[6]; break;// 2, 1, 2,
		case 49: ex *= d_[4] * d_[5] * d_[2]; break;// 2, 2, 1,
		case 50: ex *= d_[4] * d_[8]; break;// 2, 3, 0,
		case 51: ex *= d_[7] * d_[6]; break;// 3, 0, 2,
		case 52: ex *= d_[7] * d_[1] * d_[2]; break;// 3, 1, 1,
		case 53: ex *= d_[7] * d_[5]; break;// 3, 2, 0,
		case 54: ex *= d_[10] * d_[2]; break;// 4, 0, 1,
		case 55: ex *= d_[10] * d_[1]; break;// 4, 1, 0,
		case 56: ex *= d_[13]; break;// 5, 0, 0,
		case 57: ex *= d_[15] * d_[2]; break;// 0, 0, 6,
		case 58: ex *= d_[1] * d_[15]; break;// 0, 1, 5,
		case 59: ex *= d_[5] * d_[12]; break;// 0, 2, 4,
		case 60: ex *= d_[8] * d_[9]; break;// 0, 3, 3,
		case 61: ex *= d_[11] * d_[6]; break;// 0, 4, 2,
		case 62: ex *= d_[14] * d_[2]; break;// 0, 5, 1,
		case 63: ex *= d_[14] * d_[1]; break;// 0, 6, 0,
		case 64: ex *= d_[0] * d_[15]; break;// 1, 0, 5,
		case 65: ex *= d_[0] * d_[1] * d_[12]; break;// 1, 1, 4,
		case 66: ex *= d_[0] * d_[5] * d_[9]; break;// 1, 2, 3,
		case 67: ex *= d_[0] * d_[8] * d_[6]; break;// 1, 3, 2,
		case 68: ex *= d_[0] * d_[11] * d_[2]; break;// 1, 4, 1,
		case 69: ex *= d_[0] * d_[14]; break;// 1, 5, 0,
		case 70: ex *= d_[4] * d_[12]; break;// 2, 0, 4,
		case 71: ex *= d_[4] * d_[1] * d_[9]; break;// 2, 1, 3,
		case 72: ex *= d_[4] * d_[5] * d_[6]; break;// 2, 2, 2,
		case 73: ex *= d_[4] * d_[8] * d_[2]; break;// 2, 3, 1,
		case 74: ex *= d_[4] * d_[11]; break;// 2, 4, 0,
		case 75: ex *= d_[7] * d_[9]; break;// 3, 0, 3,
		case 76: ex *= d_[7] * d_[1] * d_[6]; break;// 3, 1, 2,
		case 77: ex *= d_[7] * d_[5] * d_[2]; break;// 3, 2, 1,
		case 78: ex *= d_[7] * d_[8]; break;// 3, 3, 0,
		case 79: ex *= d_[10] * d_[6]; break;// 4, 0, 2,
		case 80: ex *= d_[10] * d_[1] * d_[2]; break;// 4, 1, 1,
		case 81: ex *= d_[10] * d_[5]; break;// 4, 2, 0,
		case 82: ex *= d_[13] * d_[2]; break;// 5, 0, 1,
		case 83: ex *= d_[13] * d_[1]; break;// 5, 1, 0,
		case 84: ex *= d_[13] * d_[0]; break;// 6, 0, 0,
		case 85: ex *= d_[15] * d_[6]; break;// 0, 0, 7,
		case 86: ex *= d_[1] * d_[15] * d_[2]; break;// 0, 1, 6,
		case 87: ex *= d_[5] * d_[15]; break;// 0, 2, 5,
		case 88: ex *= d_[8] * d_[12]; break;// 0, 3, 4,
		case 89: ex *= d_[11] * d_[9]; break;// 0, 4, 3,
		case 90: ex *= d_[14] * d_[6]; break;// 0, 5, 2,
		case 91: ex *= d_[14] * d_[1] * d_[2]; break;// 0, 6, 1,
		case 92: ex *= d_[14] * d_[5]; break;// 0, 7, 0,
		case 93: ex *= d_[0] * d_[15] * d_[2]; break;// 1, 0, 6,
		case 94: ex *= d_[0] * d_[1] * d_[15]; break;// 1, 1, 5,
		case 95: ex *= d_[0] * d_[5] * d_[12]; break;// 1, 2, 4,
		case 96: ex *= d_[0] * d_[8] * d_[9]; break;// 1, 3, 3,
		case 97: ex *= d_[0] * d_[11] * d_[6]; break;// 1, 4, 2,
		case 98: ex *= d_[0] * d_[14] * d_[2]; break;// 1, 5, 1,
		case 99: ex *= d_[0] * d_[14] * d_[1]; break;// 1, 6, 0,
		case 100: ex *= d_[4] * d_[15]; break;// 2, 0, 5,
		case 101: ex *= d_[4] * d_[1] * d_[12]; break;// 2, 1, 4,
		case 102: ex *= d_[4] * d_[5] * d_[9]; break;// 2, 2, 3,
		case 103: ex *= d_[4] * d_[8] * d_[6]; break;// 2, 3, 2,
		case 104: ex *= d_[4] * d_[11] * d_[2]; break;// 2, 4, 1,
		case 105: ex *= d_[4] * d_[14]; break;// 2, 5, 0,
		case 106: ex *= d_[7] * d_[12]; break;// 3, 0, 4,
		case 107: ex *= d_[7] * d_[1] * d_[9]; break;// 3, 1, 3,
		case 108: ex *= d_[7] * d_[5] * d_[6]; break;// 3, 2, 2,
		case 109: ex *= d_[7] * d_[8] * d_[2]; break;// 3, 3, 1,
		case 110: ex *= d_[7] * d_[11]; break;// 3, 4, 0,
		case 111: ex *= d_[10] * d_[9]; break;// 4, 0, 3,
		case 112: ex *= d_[10] * d_[1] * d_[6]; break;// 4, 1, 2,
		case 113: ex *= d_[10] * d_[5] * d_[2]; break;// 4, 2, 1,
		case 114: ex *= d_[10] * d_[8]; break;// 4, 3, 0,
		case 115: ex *= d_[13] * d_[6]; break;// 5, 0, 2,
		case 116: ex *= d_[13] * d_[1] * d_[2]; break;// 5, 1, 1,
		case 117: ex *= d_[13] * d_[5]; break;// 5, 2, 0,
		case 118: ex *= d_[13] * d_[0] * d_[2]; break;// 6, 0, 1,
		case 119: ex *= d_[13] * d_[0] * d_[1]; break;// 6, 1, 0,
		case 120: ex *= d_[13] * d_[4]; break;// 7, 0, 0,
		case 121: ex *= d_[15] * d_[9]; break;// 0, 0, 8,
		case 122: ex *= d_[1] * d_[15] * d_[6]; break;// 0, 1, 7,
		case 123: ex *= d_[5] * d_[15] * d_[2]; break;// 0, 2, 6,
		case 124: ex *= d_[8] * d_[15]; break;// 0, 3, 5,
		case 125: ex *= d_[11] * d_[12]; break;// 0, 4, 4,
		case 126: ex *= d_[14] * d_[9]; break;// 0, 5, 3,
		case 127: ex *= d_[14] * d_[1] * d_[6]; break;// 0, 6, 2,
		case 128: ex *= d_[14] * d_[5] * d_[2]; break;// 0, 7, 1,
		case 129: ex *= d_[14] * d_[8]; break;// 0, 8, 0,
		case 130: ex *= d_[0] * d_[15] * d_[6]; break;// 1, 0, 7,
		case 131: ex *= d_[0] * d_[1] * d_[15] * d_[2]; break;// 1, 1, 6,
		case 132: ex *= d_[0] * d_[5] * d_[15]; break;// 1, 2, 5,
		case 133: ex *= d_[0] * d_[8] * d_[12]; break;// 1, 3, 4,
		case 134: ex *= d_[0] * d_[11] * d_[9]; break;// 1, 4, 3,
		case 135: ex *= d_[0] * d_[14] * d_[6]; break;// 1, 5, 2,
		case 136: ex *= d_[0] * d_[14] * d_[1] * d_[2]; break;// 1, 6, 1,
		case 137: ex *= d_[0] * d_[14] * d_[5]; break;// 1, 7, 0,
		case 138: ex *= d_[4] * d_[15] * d_[2]; break;// 2, 0, 6,
		case 139: ex *= d_[4] * d_[1] * d_[15]; break;// 2, 1, 5,
		case 140: ex *= d_[4] * d_[5] * d_[12]; break;// 2, 2, 4,
		case 141: ex *= d_[4] * d_[8] * d_[9]; break;// 2, 3, 3,
		case 142: ex *= d_[4] * d_[11] * d_[6]; break;// 2, 4, 2,
		case 143: ex *= d_[4] * d_[14] * d_[2]; break;// 2, 5, 1,
		case 144: ex *= d_[4] * d_[14] * d_[1]; break;// 2, 6, 0,
		case 145: ex *= d_[7] * d_[15]; break;// 3, 0, 5,
		case 146: ex *= d_[7] * d_[1] * d_[12]; break;// 3, 1, 4,
		case 147: ex *= d_[7] * d_[5] * d_[9]; break;// 3, 2, 3,
		case 148: ex *= d_[7] * d_[8] * d_[6]; break;// 3, 3, 2,
		case 149: ex *= d_[7] * d_[11] * d_[2]; break;// 3, 4, 1,
		case 150: ex *= d_[7] * d_[14]; break;// 3, 5, 0,
		case 151: ex *= d_[10] * d_[12]; break;// 4, 0, 4,
		case 152: ex *= d_[10] * d_[1] * d_[9]; break;// 4, 1, 3,
		case 153: ex *= d_[10] * d_[5] * d_[6]; break;// 4, 2, 2,
		case 154: ex *= d_[10] * d_[8] * d_[2]; break;// 4, 3, 1,
		case 155: ex *= d_[10] * d_[11]; break;// 4, 4, 0,
		case 156: ex *= d_[13] * d_[9]; break;// 5, 0, 3,
		case 157: ex *= d_[13] * d_[1] * d_[6]; break;// 5, 1, 2,
		case 158: ex *= d_[13] * d_[5] * d_[2]; break;// 5, 2, 1,
		case 159: ex *= d_[13] * d_[8]; break;// 5, 3, 0,
		case 160: ex *= d_[13] * d_[0] * d_[6]; break;// 6, 0, 2,
		case 161: ex *= d_[13] * d_[0] * d_[1] * d_[2]; break;// 6, 1, 1,
		case 162: ex *= d_[13] * d_[0] * d_[5]; break;// 6, 2, 0,
		case 163: ex *= d_[13] * d_[4] * d_[2]; break;// 7, 0, 1,
		case 164: ex *= d_[13] * d_[4] * d_[1]; break;// 7, 1, 0,
		case 165: ex *= d_[13] * d_[7]; break;// 8, 0, 0,
		case 166: ex *= d_[15] * d_[12]; break;// 0, 0, 9,
		case 167: ex *= d_[1] * d_[15] * d_[9]; break;// 0, 1, 8,
		case 168: ex *= d_[5] * d_[15] * d_[6]; break;// 0, 2, 7,
		case 169: ex *= d_[8] * d_[15] * d_[2]; break;// 0, 3, 6,
		case 170: ex *= d_[11] * d_[15]; break;// 0, 4, 5,
		case 171: ex *= d_[14] * d_[12]; break;// 0, 5, 4,
		case 172: ex *= d_[14] * d_[1] * d_[9]; break;// 0, 6, 3,
		case 173: ex *= d_[14] * d_[5] * d_[6]; break;// 0, 7, 2,
		case 174: ex *= d_[14] * d_[8] * d_[2]; break;// 0, 8, 1,
		case 175: ex *= d_[14] * d_[11]; break;// 0, 9, 0,
		case 176: ex *= d_[0] * d_[15] * d_[9]; break;// 1, 0, 8,
		case 177: ex *= d_[0] * d_[1] * d_[15] * d_[6]; break;// 1, 1, 7,
		case 178: ex *= d_[0] * d_[5] * d_[15] * d_[2]; break;// 1, 2, 6,
		case 179: ex *= d_[0] * d_[8] * d_[15]; break;// 1, 3, 5,
		case 180: ex *= d_[0] * d_[11] * d_[12]; break;// 1, 4, 4,
		case 181: ex *= d_[0] * d_[14] * d_[9]; break;// 1, 5, 3,
		case 182: ex *= d_[0] * d_[14] * d_[1] * d_[6]; break;// 1, 6, 2,
		case 183: ex *= d_[0] * d_[14] * d_[5] * d_[2]; break;// 1, 7, 1,
		case 184: ex *= d_[0] * d_[14] * d_[8]; break;// 1, 8, 0,
		case 185: ex *= d_[4] * d_[15] * d_[6]; break;// 2, 0, 7,
		case 186: ex *= d_[4] * d_[1] * d_[15] * d_[2]; break;// 2, 1, 6,
		case 187: ex *= d_[4] * d_[5] * d_[15]; break;// 2, 2, 5,
		case 188: ex *= d_[4] * d_[8] * d_[12]; break;// 2, 3, 4,
		case 189: ex *= d_[4] * d_[11] * d_[9]; break;// 2, 4, 3,
		case 190: ex *= d_[4] * d_[14] * d_[6]; break;// 2, 5, 2,
		case 191: ex *= d_[4] * d_[14] * d_[1] * d_[2]; break;// 2, 6, 1,
		case 192: ex *= d_[4] * d_[14] * d_[5]; break;// 2, 7, 0,
		case 193: ex *= d_[7] * d_[15] * d_[2]; break;// 3, 0, 6,
		case 194: ex *= d_[7] * d_[1] * d_[15]; break;// 3, 1, 5,
		case 195: ex *= d_[7] * d_[5] * d_[12]; break;// 3, 2, 4,
		case 196: ex *= d_[7] * d_[8] * d_[9]; break;// 3, 3, 3,
		case 197: ex *= d_[7] * d_[11] * d_[6]; break;// 3, 4, 2,
		case 198: ex *= d_[7] * d_[14] * d_[2]; break;// 3, 5, 1,
		case 199: ex *= d_[7] * d_[14] * d_[1]; break;// 3, 6, 0,
		case 200: ex *= d_[10] * d_[15]; break;// 4, 0, 5,
		case 201: ex *= d_[10] * d_[1] * d_[12]; break;// 4, 1, 4,
		case 202: ex *= d_[10] * d_[5] * d_[9]; break;// 4, 2, 3,
		case 203: ex *= d_[10] * d_[8] * d_[6]; break;// 4, 3, 2,
		case 204: ex *= d_[10] * d_[11] * d_[2]; break;// 4, 4, 1,
		case 205: ex *= d_[10] * d_[14]; break;// 4, 5, 0,
		case 206: ex *= d_[13] * d_[12]; break;// 5, 0, 4,
		case 207: ex *= d_[13] * d_[1] * d_[9]; break;// 5, 1, 3,
		case 208: ex *= d_[13] * d_[5] * d_[6]; break;// 5, 2, 2,
		case 209: ex *= d_[13] * d_[8] * d_[2]; break;// 5, 3, 1,
		case 210: ex *= d_[13] * d_[11]; break;// 5, 4, 0,
		case 211: ex *= d_[13] * d_[0] * d_[9]; break;// 6, 0, 3,
		case 212: ex *= d_[13] * d_[0] * d_[1] * d_[6]; break;// 6, 1, 2,
		case 213: ex *= d_[13] * d_[0] * d_[5] * d_[2]; break;// 6, 2, 1,
		case 214: ex *= d_[13] * d_[0] * d_[8]; break;// 6, 3, 0,
		case 215: ex *= d_[13] * d_[4] * d_[6]; break;// 7, 0, 2,
		case 216: ex *= d_[13] * d_[4] * d_[1] * d_[2]; break;// 7, 1, 1,
		case 217: ex *= d_[13] * d_[4] * d_[5]; break;// 7, 2, 0,
		case 218: ex *= d_[13] * d_[7] * d_[2]; break;// 8, 0, 1,
		case 219: ex *= d_[13] * d_[7] * d_[1]; break;// 8, 1, 0,
		case 220: ex *= d_[13] * d_[10]; break;// 9, 0, 0,
		case 221: ex *= d_[15] * d_[15]; break;// 0, 0, 10,
		case 222: ex *= d_[1] * d_[15] * d_[12]; break;// 0, 1, 9,
		case 223: ex *= d_[5] * d_[15] * d_[9]; break;// 0, 2, 8,
		case 224: ex *= d_[8] * d_[15] * d_[6]; break;// 0, 3, 7,
		case 225: ex *= d_[11] * d_[15] * d_[2]; break;// 0, 4, 6,
		case 226: ex *= d_[14] * d_[15]; break;// 0, 5, 5,
		case 227: ex *= d_[14] * d_[1] * d_[12]; break;// 0, 6, 4,
		case 228: ex *= d_[14] * d_[5] * d_[9]; break;// 0, 7, 3,
		case 229: ex *= d_[14] * d_[8] * d_[6]; break;// 0, 8, 2,
		case 230: ex *= d_[14] * d_[11] * d_[2]; break;// 0, 9, 1,
		case 231: ex *= d_[14] * d_[14]; break;// 0, 10, 0,
		case 232: ex *= d_[0] * d_[15] * d_[12]; break;// 1, 0, 9,
		case 233: ex *= d_[0] * d_[1] * d_[15] * d_[9]; break;// 1, 1, 8,
		case 234: ex *= d_[0] * d_[5] * d_[15] * d_[6]; break;// 1, 2, 7,
		case 235: ex *= d_[0] * d_[8] * d_[15] * d_[2]; break;// 1, 3, 6,
		case 236: ex *= d_[0] * d_[11] * d_[15]; break;// 1, 4, 5,
		case 237: ex *= d_[0] * d_[14] * d_[12]; break;// 1, 5, 4,
		case 238: ex *= d_[0] * d_[14] * d_[1] * d_[9]; break;// 1, 6, 3,
		case 239: ex *= d_[0] * d_[14] * d_[5] * d_[6]; break;// 1, 7, 2,
		case 240: ex *= d_[0] * d_[14] * d_[8] * d_[2]; break;// 1, 8, 1,
		case 241: ex *= d_[0] * d_[14] * d_[11]; break;// 1, 9, 0,
		case 242: ex *= d_[4] * d_[15] * d_[9]; break;// 2, 0, 8,
		case 243: ex *= d_[4] * d_[1] * d_[15] * d_[6]; break;// 2, 1, 7,
		case 244: ex *= d_[4] * d_[5] * d_[15] * d_[2]; break;// 2, 2, 6,
		case 245: ex *= d_[4] * d_[8] * d_[15]; break;// 2, 3, 5,
		case 246: ex *= d_[4] * d_[11] * d_[12]; break;// 2, 4, 4,
		case 247: ex *= d_[4] * d_[14] * d_[9]; break;// 2, 5, 3,
		case 248: ex *= d_[4] * d_[14] * d_[1] * d_[6]; break;// 2, 6, 2,
		case 249: ex *= d_[4] * d_[14] * d_[5] * d_[2]; break;// 2, 7, 1,
		case 250: ex *= d_[4] * d_[14] * d_[8]; break;// 2, 8, 0,
		case 251: ex *= d_[7] * d_[15] * d_[6]; break;// 3, 0, 7,
		case 252: ex *= d_[7] * d_[1] * d_[15] * d_[2]; break;// 3, 1, 6,
		case 253: ex *= d_[7] * d_[5] * d_[15]; break;// 3, 2, 5,
		case 254: ex *= d_[7] * d_[8] * d_[12]; break;// 3, 3, 4,
		case 255: ex *= d_[7] * d_[11] * d_[9]; break;// 3, 4, 3,
		case 256: ex *= d_[7] * d_[14] * d_[6]; break;// 3, 5, 2,
		case 257: ex *= d_[7] * d_[14] * d_[1] * d_[2]; break;// 3, 6, 1,
		case 258: ex *= d_[7] * d_[14] * d_[5]; break;// 3, 7, 0,
		case 259: ex *= d_[10] * d_[15] * d_[2]; break;// 4, 0, 6,
		case 260: ex *= d_[10] * d_[1] * d_[15]; break;// 4, 1, 5,
		case 261: ex *= d_[10] * d_[5] * d_[12]; break;// 4, 2, 4,
		case 262: ex *= d_[10] * d_[8] * d_[9]; break;// 4, 3, 3,
		case 263: ex *= d_[10] * d_[11] * d_[6]; break;// 4, 4, 2,
		case 264: ex *= d_[10] * d_[14] * d_[2]; break;// 4, 5, 1,
		case 265: ex *= d_[10] * d_[14] * d_[1]; break;// 4, 6, 0,
		case 266: ex *= d_[13] * d_[15]; break;// 5, 0, 5,
		case 267: ex *= d_[13] * d_[1] * d_[12]; break;// 5, 1, 4,
		case 268: ex *= d_[13] * d_[5] * d_[9]; break;// 5, 2, 3,
		case 269: ex *= d_[13] * d_[8] * d_[6]; break;// 5, 3, 2,
		case 270: ex *= d_[13] * d_[11] * d_[2]; break;// 5, 4, 1,
		case 271: ex *= d_[13] * d_[14]; break;// 5, 5, 0,
		case 272: ex *= d_[13] * d_[0] * d_[12]; break;// 6, 0, 4,
		case 273: ex *= d_[13] * d_[0] * d_[1] * d_[9]; break;// 6, 1, 3,
		case 274: ex *= d_[13] * d_[0] * d_[5] * d_[6]; break;// 6, 2, 2,
		case 275: ex *= d_[13] * d_[0] * d_[8] * d_[2]; break;// 6, 3, 1,
		case 276: ex *= d_[13] * d_[0] * d_[11]; break;// 6, 4, 0,
		case 277: ex *= d_[13] * d_[4] * d_[9]; break;// 7, 0, 3,
		case 278: ex *= d_[13] * d_[4] * d_[1] * d_[6]; break;// 7, 1, 2,
		case 279: ex *= d_[13] * d_[4] * d_[5] * d_[2]; break;// 7, 2, 1,
		case 280: ex *= d_[13] * d_[4] * d_[8]; break;// 7, 3, 0,
		case 281: ex *= d_[13] * d_[7] * d_[6]; break;// 8, 0, 2,
		case 282: ex *= d_[13] * d_[7] * d_[1] * d_[2]; break;// 8, 1, 1,
		case 283: ex *= d_[13] * d_[7] * d_[5]; break;// 8, 2, 0,
		case 284: ex *= d_[13] * d_[10] * d_[2]; break;// 9, 0, 1,
		case 285: ex *= d_[13] * d_[10] * d_[1]; break;// 9, 1, 0,
		case 286: ex *= d_[13] * d_[13]; break;// 10, 0, 0,
		default: break;
		}

		// use pointer arithmetic and cache coefficient pointer
		// This avoids repeated virtual function calls to get_coefficient_f
		double *phi_ptr = phi_data;
		const double *phi_end = phi_ptr + nmo;

		// Cache the coefficient pointer for this primitive across all MOs
		for (; phi_ptr != phi_end; ++phi_ptr)
		{
			*phi_ptr += ex;
		}
	}

	// use pointer arithmetic and minimize overhead
	const double *phi_ptr = phi_data;
	const double *phi_end = phi_ptr + nmo;

	for (; phi_ptr != phi_end; ++phi_ptr)
	{
		const double &phi_val = *phi_ptr;
		g += phi_val * phi_val;
	}

	return g;
}

double compute_spin_dens_cartesian_master_ref(
	WFN& w,
	const d3 &Pos,
	vec2 &d,
	vec &phi)
{
	const int ncen = w.get_ncen(), nex = w.get_nex(), nmo = w.get_nmo();
	const std::vector<atom>& atoms = *w.get_atoms_ptr();
	const ivec& centers = w.get_centers();
	const ivec& types = w.get_types();
	const vec& exponents = w.get_exponents();
	std::vector<MO> MOs;
	for (int m = 0; m < nmo; m++) MOs.push_back(w.get_MO(m));
	auto get_coef_primitive_major = [&w]() { return w.get_coef_primitive_major(); };
	std::fill(phi.begin(), phi.end(), 0.0);
	double alpha = 0.0, beta = 0.0, ex, *d_;
	int j;

	for (j = 0; j < ncen; j++)
	{
		d_ = d[j].data();
		d_[0] = Pos[0] - atoms[j].get_coordinate(0);
		d_[1] = Pos[1] - atoms[j].get_coordinate(1);
		d_[2] = Pos[2] - atoms[j].get_coordinate(2);
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

	const int *centers_data = centers.data();
	const int *types_data = types.data();
	const double *exponents_data = exponents.data();
	const MO *MOs_data = MOs.data();
	double *phi_data = phi.data();
	//Primitive-major coefficients: every MO for one primitive is contiguous here, where
	//MOs keeps them nex apart. Identical arithmetic in identical order - one sequential
	//read per primitive instead of nmo scattered ones, and no per-point heap allocation.
	const double *const coefs = get_coef_primitive_major();

	for (j = 0; j < nex; j++)
	{
		d_ = d[centers_data[j] - 1].data();
		ex = -exponents_data[j] * d_[3];
		if (ex < constants::exp_cutoff)
		{ // corresponds to cutoff of ex ~< 1E-20
			continue;
		}
		ex = exp(ex);
		switch (types_data[j])
		{
		case 0:  break;
		case 1:  break;// 0, 0, 0,
		case 2:  ex *= d_[0]; break;// 1, 0, 0,
		case 3:  ex *= d_[1]; break;// 0, 1, 0,
		case 4:  ex *= d_[2]; break;// 0, 0, 1,
		case 5:  ex *= d_[4]; break;// 2, 0, 0,
		case 6:  ex *= d_[5]; break;// 0, 2, 0,
		case 7:  ex *= d_[6]; break;// 0, 0, 2,
		case 8:  ex *= d_[0] * d_[1]; break;// 1, 1, 0,
		case 9:  ex *= d_[0] * d_[2]; break;// 1, 0, 1,
		case 10: ex *= d_[1] * d_[2]; break;// 0, 1, 1,
		case 11: ex *= d_[7]; break;// 3, 0, 0,
		case 12: ex *= d_[8]; break;// 0, 3, 0,
		case 13: ex *= d_[9]; break;// 0, 0, 3,
		case 14: ex *= d_[4] * d_[1]; break;// 2, 1, 0,
		case 15: ex *= d_[4] * d_[2]; break;// 2, 0, 1,
		case 16: ex *= d_[5] * d_[2]; break;// 0, 2, 1,
		case 17: ex *= d_[0] * d_[5]; break;// 1, 2, 0,
		case 18: ex *= d_[0] * d_[6]; break;// 1, 0, 2,
		case 19: ex *= d_[1] * d_[6]; break;// 0, 1, 2,
		case 20: ex *= d_[0] * d_[1] * d_[2]; break;// 1, 1, 1,
		case 21: ex *= d_[12]; break;// 0, 0, 4,
		case 22: ex *= d_[1] * d_[9]; break;// 0, 1, 3,
		case 23: ex *= d_[5] * d_[6]; break;// 0, 2, 2,
		case 24: ex *= d_[8] * d_[2]; break;// 0, 3, 1,
		case 25: ex *= d_[11]; break;// 0, 4, 0,
		case 26: ex *= d_[0] * d_[9]; break;// 1, 0, 3,
		case 27: ex *= d_[0] * d_[1] * d_[6]; break;// 1, 1, 2,
		case 28: ex *= d_[0] * d_[5] * d_[2]; break;// 1, 2, 1,
		case 29: ex *= d_[0] * d_[8]; break;// 1, 3, 0,
		case 30: ex *= d_[4] * d_[6]; break;// 2, 0, 2,
		case 31: ex *= d_[4] * d_[1] * d_[2]; break;// 2, 1, 1,
		case 32: ex *= d_[4] * d_[5]; break;// 2, 2, 0,
		case 33: ex *= d_[7] * d_[2]; break;// 3, 0, 1,
		case 34: ex *= d_[7] * d_[1]; break;// 3, 1, 0,
		case 35: ex *= d_[10]; break;// 4, 0, 0,
		case 36: ex *= d_[15]; break;// 0, 0, 5,
		case 37: ex *= d_[1] * d_[12]; break;// 0, 1, 4,
		case 38: ex *= d_[5] * d_[9]; break;// 0, 2, 3,
		case 39: ex *= d_[8] * d_[6]; break;// 0, 3, 2,
		case 40: ex *= d_[11] * d_[2]; break;// 0, 4, 1,
		case 41: ex *= d_[14]; break;// 0, 5, 0,
		case 42: ex *= d_[0] * d_[12]; break;// 1, 0, 4,
		case 43: ex *= d_[0] * d_[1] * d_[9]; break;// 1, 1, 3,
		case 44: ex *= d_[0] * d_[5] * d_[6]; break;// 1, 2, 2,
		case 45: ex *= d_[0] * d_[8] * d_[2]; break;// 1, 3, 1,
		case 46: ex *= d_[0] * d_[11]; break;// 1, 4, 0,
		case 47: ex *= d_[4] * d_[9]; break;// 2, 0, 3,
		case 48: ex *= d_[4] * d_[1] * d_[6]; break;// 2, 1, 2,
		case 49: ex *= d_[4] * d_[5] * d_[2]; break;// 2, 2, 1,
		case 50: ex *= d_[4] * d_[8]; break;// 2, 3, 0,
		case 51: ex *= d_[7] * d_[6]; break;// 3, 0, 2,
		case 52: ex *= d_[7] * d_[1] * d_[2]; break;// 3, 1, 1,
		case 53: ex *= d_[7] * d_[5]; break;// 3, 2, 0,
		case 54: ex *= d_[10] * d_[2]; break;// 4, 0, 1,
		case 55: ex *= d_[10] * d_[1]; break;// 4, 1, 0,
		case 56: ex *= d_[13]; break;// 5, 0, 0,
		case 57: ex *= d_[15] * d_[2]; break;// 0, 0, 6,
		case 58: ex *= d_[1] * d_[15]; break;// 0, 1, 5,
		case 59: ex *= d_[5] * d_[12]; break;// 0, 2, 4,
		case 60: ex *= d_[8] * d_[9]; break;// 0, 3, 3,
		case 61: ex *= d_[11] * d_[6]; break;// 0, 4, 2,
		case 62: ex *= d_[14] * d_[2]; break;// 0, 5, 1,
		case 63: ex *= d_[14] * d_[1]; break;// 0, 6, 0,
		case 64: ex *= d_[0] * d_[15]; break;// 1, 0, 5,
		case 65: ex *= d_[0] * d_[1] * d_[12]; break;// 1, 1, 4,
		case 66: ex *= d_[0] * d_[5] * d_[9]; break;// 1, 2, 3,
		case 67: ex *= d_[0] * d_[8] * d_[6]; break;// 1, 3, 2,
		case 68: ex *= d_[0] * d_[11] * d_[2]; break;// 1, 4, 1,
		case 69: ex *= d_[0] * d_[14]; break;// 1, 5, 0,
		case 70: ex *= d_[4] * d_[12]; break;// 2, 0, 4,
		case 71: ex *= d_[4] * d_[1] * d_[9]; break;// 2, 1, 3,
		case 72: ex *= d_[4] * d_[5] * d_[6]; break;// 2, 2, 2,
		case 73: ex *= d_[4] * d_[8] * d_[2]; break;// 2, 3, 1,
		case 74: ex *= d_[4] * d_[11]; break;// 2, 4, 0,
		case 75: ex *= d_[7] * d_[9]; break;// 3, 0, 3,
		case 76: ex *= d_[7] * d_[1] * d_[6]; break;// 3, 1, 2,
		case 77: ex *= d_[7] * d_[5] * d_[2]; break;// 3, 2, 1,
		case 78: ex *= d_[7] * d_[8]; break;// 3, 3, 0,
		case 79: ex *= d_[10] * d_[6]; break;// 4, 0, 2,
		case 80: ex *= d_[10] * d_[1] * d_[2]; break;// 4, 1, 1,
		case 81: ex *= d_[10] * d_[5]; break;// 4, 2, 0,
		case 82: ex *= d_[13] * d_[2]; break;// 5, 0, 1,
		case 83: ex *= d_[13] * d_[1]; break;// 5, 1, 0,
		case 84: ex *= d_[13] * d_[0]; break;// 6, 0, 0,
		case 85: ex *= d_[15] * d_[6]; break;// 0, 0, 7,
		case 86: ex *= d_[1] * d_[15] * d_[2]; break;// 0, 1, 6,
		case 87: ex *= d_[5] * d_[15]; break;// 0, 2, 5,
		case 88: ex *= d_[8] * d_[12]; break;// 0, 3, 4,
		case 89: ex *= d_[11] * d_[9]; break;// 0, 4, 3,
		case 90: ex *= d_[14] * d_[6]; break;// 0, 5, 2,
		case 91: ex *= d_[14] * d_[1] * d_[2]; break;// 0, 6, 1,
		case 92: ex *= d_[14] * d_[5]; break;// 0, 7, 0,
		case 93: ex *= d_[0] * d_[15] * d_[2]; break;// 1, 0, 6,
		case 94: ex *= d_[0] * d_[1] * d_[15]; break;// 1, 1, 5,
		case 95: ex *= d_[0] * d_[5] * d_[12]; break;// 1, 2, 4,
		case 96: ex *= d_[0] * d_[8] * d_[9]; break;// 1, 3, 3,
		case 97: ex *= d_[0] * d_[11] * d_[6]; break;// 1, 4, 2,
		case 98: ex *= d_[0] * d_[14] * d_[2]; break;// 1, 5, 1,
		case 99: ex *= d_[0] * d_[14] * d_[1]; break;// 1, 6, 0,
		case 100: ex *= d_[4] * d_[15]; break;// 2, 0, 5,
		case 101: ex *= d_[4] * d_[1] * d_[12]; break;// 2, 1, 4,
		case 102: ex *= d_[4] * d_[5] * d_[9]; break;// 2, 2, 3,
		case 103: ex *= d_[4] * d_[8] * d_[6]; break;// 2, 3, 2,
		case 104: ex *= d_[4] * d_[11] * d_[2]; break;// 2, 4, 1,
		case 105: ex *= d_[4] * d_[14]; break;// 2, 5, 0,
		case 106: ex *= d_[7] * d_[12]; break;// 3, 0, 4,
		case 107: ex *= d_[7] * d_[1] * d_[9]; break;// 3, 1, 3,
		case 108: ex *= d_[7] * d_[5] * d_[6]; break;// 3, 2, 2,
		case 109: ex *= d_[7] * d_[8] * d_[2]; break;// 3, 3, 1,
		case 110: ex *= d_[7] * d_[11]; break;// 3, 4, 0,
		case 111: ex *= d_[10] * d_[9]; break;// 4, 0, 3,
		case 112: ex *= d_[10] * d_[1] * d_[6]; break;// 4, 1, 2,
		case 113: ex *= d_[10] * d_[5] * d_[2]; break;// 4, 2, 1,
		case 114: ex *= d_[10] * d_[8]; break;// 4, 3, 0,
		case 115: ex *= d_[13] * d_[6]; break;// 5, 0, 2,
		case 116: ex *= d_[13] * d_[1] * d_[2]; break;// 5, 1, 1,
		case 117: ex *= d_[13] * d_[5]; break;// 5, 2, 0,
		case 118: ex *= d_[13] * d_[0] * d_[2]; break;// 6, 0, 1,
		case 119: ex *= d_[13] * d_[0] * d_[1]; break;// 6, 1, 0,
		case 120: ex *= d_[13] * d_[4]; break;// 7, 0, 0,
		case 121: ex *= d_[15] * d_[9]; break;// 0, 0, 8,
		case 122: ex *= d_[1] * d_[15] * d_[6]; break;// 0, 1, 7,
		case 123: ex *= d_[5] * d_[15] * d_[2]; break;// 0, 2, 6,
		case 124: ex *= d_[8] * d_[15]; break;// 0, 3, 5,
		case 125: ex *= d_[11] * d_[12]; break;// 0, 4, 4,
		case 126: ex *= d_[14] * d_[9]; break;// 0, 5, 3,
		case 127: ex *= d_[14] * d_[1] * d_[6]; break;// 0, 6, 2,
		case 128: ex *= d_[14] * d_[5] * d_[2]; break;// 0, 7, 1,
		case 129: ex *= d_[14] * d_[8]; break;// 0, 8, 0,
		case 130: ex *= d_[0] * d_[15] * d_[6]; break;// 1, 0, 7,
		case 131: ex *= d_[0] * d_[1] * d_[15] * d_[2]; break;// 1, 1, 6,
		case 132: ex *= d_[0] * d_[5] * d_[15]; break;// 1, 2, 5,
		case 133: ex *= d_[0] * d_[8] * d_[12]; break;// 1, 3, 4,
		case 134: ex *= d_[0] * d_[11] * d_[9]; break;// 1, 4, 3,
		case 135: ex *= d_[0] * d_[14] * d_[6]; break;// 1, 5, 2,
		case 136: ex *= d_[0] * d_[14] * d_[1] * d_[2]; break;// 1, 6, 1,
		case 137: ex *= d_[0] * d_[14] * d_[5]; break;// 1, 7, 0,
		case 138: ex *= d_[4] * d_[15] * d_[2]; break;// 2, 0, 6,
		case 139: ex *= d_[4] * d_[1] * d_[15]; break;// 2, 1, 5,
		case 140: ex *= d_[4] * d_[5] * d_[12]; break;// 2, 2, 4,
		case 141: ex *= d_[4] * d_[8] * d_[9]; break;// 2, 3, 3,
		case 142: ex *= d_[4] * d_[11] * d_[6]; break;// 2, 4, 2,
		case 143: ex *= d_[4] * d_[14] * d_[2]; break;// 2, 5, 1,
		case 144: ex *= d_[4] * d_[14] * d_[1]; break;// 2, 6, 0,
		case 145: ex *= d_[7] * d_[15]; break;// 3, 0, 5,
		case 146: ex *= d_[7] * d_[1] * d_[12]; break;// 3, 1, 4,
		case 147: ex *= d_[7] * d_[5] * d_[9]; break;// 3, 2, 3,
		case 148: ex *= d_[7] * d_[8] * d_[6]; break;// 3, 3, 2,
		case 149: ex *= d_[7] * d_[11] * d_[2]; break;// 3, 4, 1,
		case 150: ex *= d_[7] * d_[14]; break;// 3, 5, 0,
		case 151: ex *= d_[10] * d_[12]; break;// 4, 0, 4,
		case 152: ex *= d_[10] * d_[1] * d_[9]; break;// 4, 1, 3,
		case 153: ex *= d_[10] * d_[5] * d_[6]; break;// 4, 2, 2,
		case 154: ex *= d_[10] * d_[8] * d_[2]; break;// 4, 3, 1,
		case 155: ex *= d_[10] * d_[11]; break;// 4, 4, 0,
		case 156: ex *= d_[13] * d_[9]; break;// 5, 0, 3,
		case 157: ex *= d_[13] * d_[1] * d_[6]; break;// 5, 1, 2,
		case 158: ex *= d_[13] * d_[5] * d_[2]; break;// 5, 2, 1,
		case 159: ex *= d_[13] * d_[8]; break;// 5, 3, 0,
		case 160: ex *= d_[13] * d_[0] * d_[6]; break;// 6, 0, 2,
		case 161: ex *= d_[13] * d_[0] * d_[1] * d_[2]; break;// 6, 1, 1,
		case 162: ex *= d_[13] * d_[0] * d_[5]; break;// 6, 2, 0,
		case 163: ex *= d_[13] * d_[4] * d_[2]; break;// 7, 0, 1,
		case 164: ex *= d_[13] * d_[4] * d_[1]; break;// 7, 1, 0,
		case 165: ex *= d_[13] * d_[7]; break;// 8, 0, 0,
		case 166: ex *= d_[15] * d_[12]; break;// 0, 0, 9,
		case 167: ex *= d_[1] * d_[15] * d_[9]; break;// 0, 1, 8,
		case 168: ex *= d_[5] * d_[15] * d_[6]; break;// 0, 2, 7,
		case 169: ex *= d_[8] * d_[15] * d_[2]; break;// 0, 3, 6,
		case 170: ex *= d_[11] * d_[15]; break;// 0, 4, 5,
		case 171: ex *= d_[14] * d_[12]; break;// 0, 5, 4,
		case 172: ex *= d_[14] * d_[1] * d_[9]; break;// 0, 6, 3,
		case 173: ex *= d_[14] * d_[5] * d_[6]; break;// 0, 7, 2,
		case 174: ex *= d_[14] * d_[8] * d_[2]; break;// 0, 8, 1,
		case 175: ex *= d_[14] * d_[11]; break;// 0, 9, 0,
		case 176: ex *= d_[0] * d_[15] * d_[9]; break;// 1, 0, 8,
		case 177: ex *= d_[0] * d_[1] * d_[15] * d_[6]; break;// 1, 1, 7,
		case 178: ex *= d_[0] * d_[5] * d_[15] * d_[2]; break;// 1, 2, 6,
		case 179: ex *= d_[0] * d_[8] * d_[15]; break;// 1, 3, 5,
		case 180: ex *= d_[0] * d_[11] * d_[12]; break;// 1, 4, 4,
		case 181: ex *= d_[0] * d_[14] * d_[9]; break;// 1, 5, 3,
		case 182: ex *= d_[0] * d_[14] * d_[1] * d_[6]; break;// 1, 6, 2,
		case 183: ex *= d_[0] * d_[14] * d_[5] * d_[2]; break;// 1, 7, 1,
		case 184: ex *= d_[0] * d_[14] * d_[8]; break;// 1, 8, 0,
		case 185: ex *= d_[4] * d_[15] * d_[6]; break;// 2, 0, 7,
		case 186: ex *= d_[4] * d_[1] * d_[15] * d_[2]; break;// 2, 1, 6,
		case 187: ex *= d_[4] * d_[5] * d_[15]; break;// 2, 2, 5,
		case 188: ex *= d_[4] * d_[8] * d_[12]; break;// 2, 3, 4,
		case 189: ex *= d_[4] * d_[11] * d_[9]; break;// 2, 4, 3,
		case 190: ex *= d_[4] * d_[14] * d_[6]; break;// 2, 5, 2,
		case 191: ex *= d_[4] * d_[14] * d_[1] * d_[2]; break;// 2, 6, 1,
		case 192: ex *= d_[4] * d_[14] * d_[5]; break;// 2, 7, 0,
		case 193: ex *= d_[7] * d_[15] * d_[2]; break;// 3, 0, 6,
		case 194: ex *= d_[7] * d_[1] * d_[15]; break;// 3, 1, 5,
		case 195: ex *= d_[7] * d_[5] * d_[12]; break;// 3, 2, 4,
		case 196: ex *= d_[7] * d_[8] * d_[9]; break;// 3, 3, 3,
		case 197: ex *= d_[7] * d_[11] * d_[6]; break;// 3, 4, 2,
		case 198: ex *= d_[7] * d_[14] * d_[2]; break;// 3, 5, 1,
		case 199: ex *= d_[7] * d_[14] * d_[1]; break;// 3, 6, 0,
		case 200: ex *= d_[10] * d_[15]; break;// 4, 0, 5,
		case 201: ex *= d_[10] * d_[1] * d_[12]; break;// 4, 1, 4,
		case 202: ex *= d_[10] * d_[5] * d_[9]; break;// 4, 2, 3,
		case 203: ex *= d_[10] * d_[8] * d_[6]; break;// 4, 3, 2,
		case 204: ex *= d_[10] * d_[11] * d_[2]; break;// 4, 4, 1,
		case 205: ex *= d_[10] * d_[14]; break;// 4, 5, 0,
		case 206: ex *= d_[13] * d_[12]; break;// 5, 0, 4,
		case 207: ex *= d_[13] * d_[1] * d_[9]; break;// 5, 1, 3,
		case 208: ex *= d_[13] * d_[5] * d_[6]; break;// 5, 2, 2,
		case 209: ex *= d_[13] * d_[8] * d_[2]; break;// 5, 3, 1,
		case 210: ex *= d_[13] * d_[11]; break;// 5, 4, 0,
		case 211: ex *= d_[13] * d_[0] * d_[9]; break;// 6, 0, 3,
		case 212: ex *= d_[13] * d_[0] * d_[1] * d_[6]; break;// 6, 1, 2,
		case 213: ex *= d_[13] * d_[0] * d_[5] * d_[2]; break;// 6, 2, 1,
		case 214: ex *= d_[13] * d_[0] * d_[8]; break;// 6, 3, 0,
		case 215: ex *= d_[13] * d_[4] * d_[6]; break;// 7, 0, 2,
		case 216: ex *= d_[13] * d_[4] * d_[1] * d_[2]; break;// 7, 1, 1,
		case 217: ex *= d_[13] * d_[4] * d_[5]; break;// 7, 2, 0,
		case 218: ex *= d_[13] * d_[7] * d_[2]; break;// 8, 0, 1,
		case 219: ex *= d_[13] * d_[7] * d_[1]; break;// 8, 1, 0,
		case 220: ex *= d_[13] * d_[10]; break;// 9, 0, 0,
		case 221: ex *= d_[15] * d_[15]; break;// 0, 0, 10,
		case 222: ex *= d_[1] * d_[15] * d_[12]; break;// 0, 1, 9,
		case 223: ex *= d_[5] * d_[15] * d_[9]; break;// 0, 2, 8,
		case 224: ex *= d_[8] * d_[15] * d_[6]; break;// 0, 3, 7,
		case 225: ex *= d_[11] * d_[15] * d_[2]; break;// 0, 4, 6,
		case 226: ex *= d_[14] * d_[15]; break;// 0, 5, 5,
		case 227: ex *= d_[14] * d_[1] * d_[12]; break;// 0, 6, 4,
		case 228: ex *= d_[14] * d_[5] * d_[9]; break;// 0, 7, 3,
		case 229: ex *= d_[14] * d_[8] * d_[6]; break;// 0, 8, 2,
		case 230: ex *= d_[14] * d_[11] * d_[2]; break;// 0, 9, 1,
		case 231: ex *= d_[14] * d_[14]; break;// 0, 10, 0,
		case 232: ex *= d_[0] * d_[15] * d_[12]; break;// 1, 0, 9,
		case 233: ex *= d_[0] * d_[1] * d_[15] * d_[9]; break;// 1, 1, 8,
		case 234: ex *= d_[0] * d_[5] * d_[15] * d_[6]; break;// 1, 2, 7,
		case 235: ex *= d_[0] * d_[8] * d_[15] * d_[2]; break;// 1, 3, 6,
		case 236: ex *= d_[0] * d_[11] * d_[15]; break;// 1, 4, 5,
		case 237: ex *= d_[0] * d_[14] * d_[12]; break;// 1, 5, 4,
		case 238: ex *= d_[0] * d_[14] * d_[1] * d_[9]; break;// 1, 6, 3,
		case 239: ex *= d_[0] * d_[14] * d_[5] * d_[6]; break;// 1, 7, 2,
		case 240: ex *= d_[0] * d_[14] * d_[8] * d_[2]; break;// 1, 8, 1,
		case 241: ex *= d_[0] * d_[14] * d_[11]; break;// 1, 9, 0,
		case 242: ex *= d_[4] * d_[15] * d_[9]; break;// 2, 0, 8,
		case 243: ex *= d_[4] * d_[1] * d_[15] * d_[6]; break;// 2, 1, 7,
		case 244: ex *= d_[4] * d_[5] * d_[15] * d_[2]; break;// 2, 2, 6,
		case 245: ex *= d_[4] * d_[8] * d_[15]; break;// 2, 3, 5,
		case 246: ex *= d_[4] * d_[11] * d_[12]; break;// 2, 4, 4,
		case 247: ex *= d_[4] * d_[14] * d_[9]; break;// 2, 5, 3,
		case 248: ex *= d_[4] * d_[14] * d_[1] * d_[6]; break;// 2, 6, 2,
		case 249: ex *= d_[4] * d_[14] * d_[5] * d_[2]; break;// 2, 7, 1,
		case 250: ex *= d_[4] * d_[14] * d_[8]; break;// 2, 8, 0,
		case 251: ex *= d_[7] * d_[15] * d_[6]; break;// 3, 0, 7,
		case 252: ex *= d_[7] * d_[1] * d_[15] * d_[2]; break;// 3, 1, 6,
		case 253: ex *= d_[7] * d_[5] * d_[15]; break;// 3, 2, 5,
		case 254: ex *= d_[7] * d_[8] * d_[12]; break;// 3, 3, 4,
		case 255: ex *= d_[7] * d_[11] * d_[9]; break;// 3, 4, 3,
		case 256: ex *= d_[7] * d_[14] * d_[6]; break;// 3, 5, 2,
		case 257: ex *= d_[7] * d_[14] * d_[1] * d_[2]; break;// 3, 6, 1,
		case 258: ex *= d_[7] * d_[14] * d_[5]; break;// 3, 7, 0,
		case 259: ex *= d_[10] * d_[15] * d_[2]; break;// 4, 0, 6,
		case 260: ex *= d_[10] * d_[1] * d_[15]; break;// 4, 1, 5,
		case 261: ex *= d_[10] * d_[5] * d_[12]; break;// 4, 2, 4,
		case 262: ex *= d_[10] * d_[8] * d_[9]; break;// 4, 3, 3,
		case 263: ex *= d_[10] * d_[11] * d_[6]; break;// 4, 4, 2,
		case 264: ex *= d_[10] * d_[14] * d_[2]; break;// 4, 5, 1,
		case 265: ex *= d_[10] * d_[14] * d_[1]; break;// 4, 6, 0,
		case 266: ex *= d_[13] * d_[15]; break;// 5, 0, 5,
		case 267: ex *= d_[13] * d_[1] * d_[12]; break;// 5, 1, 4,
		case 268: ex *= d_[13] * d_[5] * d_[9]; break;// 5, 2, 3,
		case 269: ex *= d_[13] * d_[8] * d_[6]; break;// 5, 3, 2,
		case 270: ex *= d_[13] * d_[11] * d_[2]; break;// 5, 4, 1,
		case 271: ex *= d_[13] * d_[14]; break;// 5, 5, 0,
		case 272: ex *= d_[13] * d_[0] * d_[12]; break;// 6, 0, 4,
		case 273: ex *= d_[13] * d_[0] * d_[1] * d_[9]; break;// 6, 1, 3,
		case 274: ex *= d_[13] * d_[0] * d_[5] * d_[6]; break;// 6, 2, 2,
		case 275: ex *= d_[13] * d_[0] * d_[8] * d_[2]; break;// 6, 3, 1,
		case 276: ex *= d_[13] * d_[0] * d_[11]; break;// 6, 4, 0,
		case 277: ex *= d_[13] * d_[4] * d_[9]; break;// 7, 0, 3,
		case 278: ex *= d_[13] * d_[4] * d_[1] * d_[6]; break;// 7, 1, 2,
		case 279: ex *= d_[13] * d_[4] * d_[5] * d_[2]; break;// 7, 2, 1,
		case 280: ex *= d_[13] * d_[4] * d_[8]; break;// 7, 3, 0,
		case 281: ex *= d_[13] * d_[7] * d_[6]; break;// 8, 0, 2,
		case 282: ex *= d_[13] * d_[7] * d_[1] * d_[2]; break;// 8, 1, 1,
		case 283: ex *= d_[13] * d_[7] * d_[5]; break;// 8, 2, 0,
		case 284: ex *= d_[13] * d_[10] * d_[2]; break;// 9, 0, 1,
		case 285: ex *= d_[13] * d_[10] * d_[1]; break;// 9, 1, 0,
		case 286: ex *= d_[13] * d_[13]; break;// 10, 0, 0,
		default: break;
		}
		// use pointer arithmetic and cache coefficient pointer
		// This avoids repeated virtual function calls to get_coefficient_f
		const double *c_row = coefs + (size_t)j * nmo;
		double *phi_ptr = phi_data;
		for (int mo = 0; mo < nmo; ++mo, ++phi_ptr)
		{
			*phi_ptr += c_row[mo] * ex;
		}
	}

	// use pointer arithmetic and minimize overhead
	const double *phi_ptr = phi_data;
	const double *phi_end = phi_ptr + nmo;
	const MO *mo_ptr = MOs_data;

	for (; phi_ptr != phi_end; ++phi_ptr, ++mo_ptr)
	{
		const double &phi_val = *phi_ptr;
		if (mo_ptr->get_op())
			beta += mo_ptr->get_occ() * phi_val * phi_val;
		else
			alpha += mo_ptr->get_occ() * phi_val * phi_val;
	}

	return alpha - beta;
}

	uint64_t bits(const double x) { return std::bit_cast<uint64_t>(x); }

	//Three centres and one primitive of every type 0..287 (0 and 287 take the default branch), spread
	//over the centres with exponents 0.15-1.2 so l = 10 still contributes a few bohr out.
	//spin 0: every MO alpha; 1: alternating alpha/beta; 2: all beta. Occupations include 0 and fractions.
	WFN make_all_types(const int spin)
	{
		WFN w(e_origin::NOT_YET_DEFINED);
		w.push_back_atom("C", 0.0, 0.0, 0.0, 6);
		w.push_back_atom("O", 1.3, -0.4, 0.9, 8);
		w.push_back_atom("H", -1.1, 1.7, -0.6, 1);
		const double occ[5] = { 2.0, 1.0, 0.5, 0.0, 1.5 };
		for (int m = 0; m < 5; m++)
			w.push_back_MO(m + 1, occ[m], -1.0 + 0.3 * m, spin == 0 ? 0 : spin == 2 ? 1 : m % 2);
		std::mt19937 gen(4711);
		std::uniform_real_distribution<double> c(-1.0, 1.0);
		for (int t = 0; t <= 287; t++)
		{
			double v[5];
			for (int m = 0; m < 5; m++) v[m] = c(gen);
			w.add_primitive(1 + t % 3, t, 0.15 + 0.05 * (t % 22), v);
		}
		w.set_exp_cutoff();
		return w;
	}

	std::vector<d3> probe_points()
	{
		std::vector<d3> p{ { 0.0, 0.0, 0.0 }, { 1.3, -0.4, 0.9 }, { -1.1, 1.7, -0.6 }, { -0.0, 0.0, -0.0 },
			{ 1e-9, -1e-12, 3e-15 }, { 0.5, 0.5, 0.5 }, { -2.5, -3.0, 2.0 }, { 40.0, -35.0, 50.0 }, { 1e3, 1e3, -1e3 },
			{ 6.0, 0.1, -0.2 } };
		std::mt19937 gen(12345);
		std::uniform_real_distribution<double> u(-5.0, 5.0);
		for (int i = 0; i < 300; i++) p.push_back({ u(gen), u(gen), u(gen) });
		return p;
	}

	//|a - b| against a scale: phi and the spin density are sums of signed terms
	int n_not_bitwise = 0;
	bool close(const double a, const double b, const double scale)
	{
		n_not_bitwise += bits(a) != bits(b);
		return std::fabs(a - b) <= 1e-12 * scale;
	}

	//compute_g_cartesian, both compute_spin_dens and compute_dens over make_all_types at every probe point, then every
	//type alone at random points; the bits of each value and phi go to out. The master copies above are held to 1e-12
	//only: compiled in this file under /fp:fast they differ from wfn_density.cpp in the last bits already on master.
	//The bitwise comparison is two builds' dumps: NOS_WFN_KERNEL_DUMP writes one, NOS_WFN_KERNEL_REF compares with one.
	//computeRhoELI and computeGrad are only dumped
	int sweep(std::vector<uint64_t>& out)
	{
		int bad = 0;
		for (int spin = 0; spin < 3; spin++)
		{
			WFN w = make_all_types(spin);
			const int nmo = w.get_nmo();
			vec2 d(w.get_ncen(), vec(16, 0.0)), d_ref(w.get_ncen(), vec(16, 0.0));
			vec phi(nmo, 0.0), phi_ref(nmo, 0.0);
			for (const d3& P : probe_points())
			{
				const double g = w.compute_g_cartesian(P, d, phi);
				const double g_ref = compute_g_cartesian_master_ref(w, P, d_ref, phi_ref);
				bad += !close(g, g_ref, g_ref);
				out.push_back(bits(g));
				for (int m = 0; m < nmo; m++)
				{
					bad += !close(phi[m], phi_ref[m], std::sqrt(g_ref));
					out.push_back(bits(phi[m]));
				}
				const double rho = w.compute_dens(P);
				const double s = w.compute_spin_dens(P, d, phi);
				const double s_ref = compute_spin_dens_cartesian_master_ref(w, P, d_ref, phi_ref);
				bad += !close(s, s_ref, rho);
				out.push_back(bits(s));
				out.push_back(bits(rho));
				for (int m = 0; m < nmo; m++)
				{
					bad += !close(phi[m], phi_ref[m], std::sqrt(g_ref));
					out.push_back(bits(phi[m]));
				}
				const double s1 = w.compute_spin_dens(P);
				bad += bits(s1) != bits(s);
				out.push_back(bits(s1));
				double r_eli = 0.0, eli = 0.0, r_grad = 0.0;
				d3 grad{ 0.0, 0.0, 0.0 };
				w.computeRhoELI(P, r_eli, eli);
				w.computeGrad(P, grad, &r_grad);
				for (const double v : { r_eli, eli, grad[0], grad[1], grad[2], r_grad }) out.push_back(bits(v));
			}
		}
		std::mt19937 gen(99);
		std::uniform_real_distribution<double> u(-2.0, 2.0);
		for (int t = 0; t <= 287; t++)
		{
			WFN w(e_origin::NOT_YET_DEFINED);
			w.push_back_atom("C", 0.1, -0.2, 0.3, 6);
			w.push_back_MO(1, 1.0, -1.0, 0);
			w.push_back_MO(2, 1.0, -0.5, 1);
			double v[2] = { 0.8, -0.3 };
			w.add_primitive(1, t, 0.4, v);
			w.set_exp_cutoff();
			vec2 d(1, vec(16, 0.0)), d_ref(1, vec(16, 0.0));
			vec phi(2, 0.0), phi_ref(2, 0.0);
			for (int i = 0; i < 20; i++)
			{
				const d3 P{ u(gen), u(gen), u(gen) };
				const double g = w.compute_g_cartesian(P, d, phi), g_ref = compute_g_cartesian_master_ref(w, P, d_ref, phi_ref);
				bad += !close(g, g_ref, g_ref) || !close(phi[0], phi_ref[0], std::sqrt(g_ref));
				out.push_back(bits(g));
				out.push_back(bits(phi[0]));
				const double s = w.compute_spin_dens(P, d, phi), s_ref = compute_spin_dens_cartesian_master_ref(w, P, d_ref, phi_ref);
				bad += !close(s, s_ref, g_ref) || !close(phi[0], phi_ref[0], std::sqrt(g_ref)) || !close(phi[1], phi_ref[1], std::sqrt(g_ref));
				out.push_back(bits(s));
				out.push_back(bits(phi[0]));
				out.push_back(bits(phi[1]));
			}
		}
		return bad;
	}

	TEST(WfnCleanupTests, CartesianKernelsMatchMaster)
	{
		std::vector<uint64_t> out;
		EXPECT_EQ(sweep(out), 0);
		std::cout << "values not bitwise equal to the master copies: " << n_not_bitwise << "\n";
		int nonzero = 0;
		for (const uint64_t b : out) nonzero += (b << 1) != 0;
		EXPECT_GT(nonzero, (int)out.size() / 2);
		if (const char* f = std::getenv("NOS_WFN_KERNEL_DUMP"))
		{
			std::ofstream o(f, std::ios::binary);
			o.write(reinterpret_cast<const char*>(out.data()), out.size() * sizeof(uint64_t));
			EXPECT_TRUE(o.good());
		}
		//A dump of another build (master's): every value bit for bit
		if (const char* f = std::getenv("NOS_WFN_KERNEL_REF"))
		{
			std::ifstream i(f, std::ios::binary);
			std::vector<uint64_t> ref(out.size() + 1, 0);
			i.read(reinterpret_cast<char*>(ref.data()), ref.size() * sizeof(uint64_t));
			ASSERT_EQ((size_t)i.gcount(), out.size() * sizeof(uint64_t));
			int diff = 0;
			for (size_t k = 0; k < out.size(); k++) diff += ref[k] != out[k];
			EXPECT_EQ(diff, 0);
		}
	}
}
