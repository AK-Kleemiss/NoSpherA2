#include "pch.h"
#include "wfn_class.h"
#include "convenience.h"
#include "mo_class.h"
#include "cube.h"
#include "constants.h"
#include "fchk.h"
#include "basis_set.h"
#include "nos_math.h"
#include "libCintMain.h"
#include "integrator.h"
#include "cell.h"


const double WFN::compute_dens(
	const d3 &Pos,
	vec2 &d,
	vec &phi) const
{
	//if (d_f_switch)
	//{
	//    err_checkf(d.size() >= 5, "d is too small!", std::cout);
	//    err_checkf(phi.size() >= get_nmo(true), "phi is too small!", std::cout);
	//    return compute_dens_spherical(Pos1, Pos2, Pos3, d, phi);
	//}
	//else
	//{
	err_checkf(d.size() >= ncen, "d is too small!", std::cout);
	//the kernel accumulates every MO, virtuals included; callers size phi by the occupied count
	if (phi.size() < (size_t)nmo) phi.resize(nmo);
	return compute_dens_cartesian(Pos, d, phi);
	//}
};

const double WFN::compute_dens(
	const d3 &Pos)
	const
{
	//Called per point of the basin lookup: the ncen + 1 allocations per call were
	//most of its allocator time, so the buffers persist per thread (compute_dens_cartesian zeroes phi)
	thread_local vec2 d;
	thread_local vec phi;
	if (phi.size() < (size_t)nmo) phi.resize(nmo);
	if (d.size() < (size_t)ncen) d.resize(ncen, vec(16, 0.0));
	return compute_dens_cartesian(Pos, d, phi);
};

const double WFN::compute_spin_dens(
	const d3 &Pos,
	vec2 &d,
	vec &phi) const
{
	//if (d_f_switch)
	//{
	//    err_checkf(d.size() >= 5, "d is too small!", std::cout);
	//    err_checkf(phi.size() >= get_nmo(true), "phi is too small!", std::cout);
	//    err_not_impl_f("Nah.. not yet implemented correctly", std::cout);
	//    return compute_dens_spherical(Pos1, Pos2, Pos3, d, phi);
	//}
	//else
	//{
	err_checkf(d.size() >= ncen, "d is too small!", std::cout);
	if (phi.size() < (size_t)nmo) phi.resize(nmo);
	return compute_spin_dens_cartesian(Pos, d, phi);
	//}
};

const double WFN::compute_spin_dens(
	const d3 &Pos) const
{
	vec2 d;
	vec phi(nmo, 0.0);

	//if (d_f_switch)
	//{
	//    d.resize(5);
	//    for (int i = 0; i < 5; i++)
	//        d[i].resize(ncen, 0.0);
	//    err_not_impl_f("Nah.. not yet implemented correctly", std::cout);
	//    return compute_dens_spherical(Pos1, Pos2, Pos3, d, phi);
	//}
	//else
	//{
	d.resize(ncen);
	for (int i = 0; i < ncen; i++)
		d[i].resize(16, 0.0);
	return compute_spin_dens_cartesian(Pos, d, phi);
	//}
};

const double WFN::compute_dens_cartesian(
	const d3 &Pos,
	vec2 &d,
	vec &phi) const
{
	std::fill(phi.begin(), phi.end(), 0.0);
	double Rho = 0.0;
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
	const MO *MOs_data = MOs.data();
	double *phi_data = phi.data();
	//Primitive-major: every MO for one primitive is contiguous, where MOs keeps them nex
	//apart. Same arithmetic in the same order, and no per-point allocation.
	const double *const coefs = get_coef_primitive_major();
	thread_local vec exps;
	if (exps.size() < group_exponent.size()) exps.resize(group_exponent.size());
	exp_table([&d](const int c) { return d[c][3]; }, exps.data());
	const int *group = prim_exp_group.data();

	for (j = 0; j < nex; j++)
	{
		ex = exps[group[j]];
		if (ex == 0.0)
		{ // corresponds to cutoff of maximum density contribution of 1E-5
			continue;
		}
		d_ = d[centers_data[j] - 1].data();
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
		Rho += mo_ptr->get_occ() * phi_val * phi_val;
	}

	return Rho;
}

const double WFN::eval_ao(
	std::array<double, 4>& d,
	const std::vector<primitive>& prims,
	const int &m
) const
{

	// normalize distances for spherical harmonic
	const int type = prims[0].get_type();
	const double scale = 1 / d[3];
	d[0] *= scale;
	d[1] *= scale;
	d[2] *= scale;
	const double rl = std::pow(d[3], type);
	d[3] *= d[3];

	double radial = 0;
	const primitive* p = prims.data();
	const primitive* const p_end = p + prims.size();

	for (; p != p_end; p++) {
		radial += p->eval_gaussian_unnormalized(rl, d[3]);
	}

	return radial * constants::spherical_harmonic(type, m, d.data());
	// err_checkf(coef_counter == exp_coefs, "WRONG NUMBER OF COEFFICIENTS! " + std::to_string(coef_counter) + " vs. " + std::to_string(exp_coefs), std::cout);
}

const double WFN::compute_g_cartesian(
	const d3 &Pos,
	vec2 &d,
	vec &phi) const
{
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


const double WFN::compute_spin_dens_cartesian(
	const d3 &Pos,
	vec2 &d,
	vec &phi) const
{
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

const double WFN::compute_MO_spherical(
	const d3 &Pos,
	const int &MO) const
{
	err_not_impl_f("This one is not tested an will most likely not work, therefore aborting!", std::cout);
	return 0.0;
	(void)MO;
	//err_checkf(d_f_switch, "Only works for spheriacl wavefunctions!", std::cout);
	//int iat;
	//int l = 0;
	//// ex will carry information about radial function
	//double ex;
	//vec2 d(5);
	//for (int i = 0; i < 5; i++)
	//    d[i].resize(ncen);
	//double phi(0.0);

	//for (iat = 0; iat < ncen; iat++)
	//{
	//    d[0][iat] = Pos1 - atoms[iat].x;
	//    d[1][iat] = Pos2 - atoms[iat].y;
	//    d[2][iat] = Pos3 - atoms[iat].z;
	//    d[3][iat] = d[0][iat] * d[0][iat] + d[1][iat] * d[1][iat] + d[2][iat] * d[2][iat];
	//    d[4][iat] = sqrt(d[3][iat]);
	//}
	///*Here d[0] = x
	//             d[1] = y
	//             d[2] = z
	//             d[3] = r^2
	//             d[4] = r
	//             */

	//for (int j = 0; j < nex; j++)
	//{
	//    iat = centers[j] - 1;
	//    ex = -exponents[j] * d[3][iat];
	//    if (ex < -46.0517)
	//    { // corresponds to cutoff of ex ~< 1E-20
	//        continue;
	//    }
	//    // apply radial function:
	//    ex = exp(ex);
	//    int lam = 0, m = 0;
	//    if (l > 1)
	//    {
	//        if (l <= 4)
	//        {
	//            ex *= d[4][iat];
	//            lam = 1;
	//            if (l == 2)
	//                m = 0;
	//            else if (l == 3)
	//                m = 1;
	//            else if (l == 4)
	//                m = -1;
	//        }
	//        else if (l <= 9)
	//        {
	//            ex *= d[3][iat];
	//            lam = 2;
	//            if (l == 5)
	//                m = 0;
	//            else if (l == 6)
	//                m = 1;
	//            else if (l == 7)
	//                m = -1;
	//            else if (l == 8)
	//                m = 2;
	//            else if (l == 9)
	//                m = -2;
	//        }
	//        else if (l <= 16)
	//        {
	//            ex *= d[3][iat] * d[4][iat];
	//            lam = 3;
	//            if (l == 10)
	//                m = 0;
	//            else if (l == 11)
	//                m = 1;
	//            else if (l == 12)
	//                m = -1;
	//            else if (l == 13)
	//                m = 2;
	//            else if (l == 14)
	//                m = -2;
	//            else if (l == 15)
	//                m = 3;
	//            else if (l == 16)
	//                m = -3;
	//        }
	//        else if (l <= 25)
	//        {
	//            ex *= d[3][iat] * d[3][iat];
	//            lam = 4;
	//            if (l == 17)
	//                m = 0;
	//            else if (l == 18)
	//                m = 1;
	//            else if (l == 19)
	//                m = -1;
	//            else if (l == 20)
	//                m = 2;
	//            else if (l == 21)
	//                m = -2;
	//            else if (l == 22)
	//                m = 3;
	//            else if (l == 23)
	//                m = -3;
	//            else if (l == 24)
	//                m = 4;
	//            else if (l == 25)
	//                m = -4;
	//        }
	//        else if (l <= 36)
	//        {
	//            ex *= pow(d[4][iat], 5);
	//            lam = 5;
	//            if (l == 26)
	//                m = 0;
	//            else if (l == 27)
	//                m = 1;
	//            else if (l == 28)
	//                m = -1;
	//            else if (l == 29)
	//                m = 2;
	//            else if (l == 30)
	//                m = -2;
	//            else if (l == 31)
	//                m = 3;
	//            else if (l == 32)
	//                m = -3;
	//            else if (l == 33)
	//                m = 4;
	//            else if (l == 34)
	//                m = -4;
	//            else if (l == 35)
	//                m = 5;
	//            else if (l == 36)
	//                m = -5;
	//        }
	//    }
	//    // calc spherical harmonic
	//    double d_t[]{d[0][iat], d[1][iat], d[2][iat], d[3][iat], d[4][iat]};
	//    double SH = spherical_harmonic(lam, m, d_t);
	//    SH *= ex;                                 // multiply radial part with spherical harmonic
	//    phi += MOs[MO].get_coefficient_f(j) * SH; // build MO values at this point
	//}
	//shrink_vector<vec>(d);

	//return phi;
}

const double WFN::compute_dens_spherical(
	const d3 &Pos,
	vec2 &d,
	vec &phi) const
{
	err_not_impl_f("This one is not tested an will most likely not work, therefore aborting!", std::cout);
	return 0.0;
	(void)d;
	(void)phi;
	//err_checkf(d_f_switch, "Only works for spheriacl wavefunctions!", std::cout);
	//std::fill(phi.begin(), phi.end(), 0.0);
	//double Rho = 0.0;
	//int iat;
	//int l;
	//// ex will carry information about radial function
	//double ex;
	//int mo;

	//for (iat = 0; iat < ncen; iat++)
	//{
	//    d[0][iat] = Pos1 - atoms[iat].x;
	//    d[1][iat] = Pos2 - atoms[iat].y;
	//    d[2][iat] = Pos3 - atoms[iat].z;
	//    d[3][iat] = d[0][iat] * d[0][iat] + d[1][iat] * d[1][iat] + d[2][iat] * d[2][iat];
	//    d[4][iat] = sqrt(d[3][iat]);
	//}
	///*Here d[0] = x
	//             d[1] = y
	//             d[2] = z
	//             d[3] = r^2
	//             d[4] = r
	//             */
	//for (int j = 0; j < nex; j++)
	//{
	//    iat = centers[j] - 1;
	//    ex = -exponents[j] * d[3][iat];
	//    if (ex < -46.0517)
	//    { // corresponds to cutoff of ex ~< 1E-20
	//        continue;
	//    }
	//    ex = exp(ex);
	//    // apply radial function:
	//    l = types[j];
	//    int lam = 0, m = 0;
	//    if (l > 1)
	//    {
	//        if (l <= 4)
	//        {
	//            ex *= d[4][iat];
	//            lam = 1;
	//            if (l == 2)
	//                m = 0;
	//            else if (l == 3)
	//                m = 1;
	//            else if (l == 4)
	//                m = -1;
	//        }
	//        else if (l <= 9)
	//        {
	//            ex *= d[3][iat];
	//            lam = 2;
	//            if (l == 5)
	//                m = 0;
	//            else if (l == 6)
	//                m = 1;
	//            else if (l == 7)
	//                m = -1;
	//            else if (l == 8)
	//                m = 2;
	//            else if (l == 9)
	//                m = -2;
	//        }
	//        else if (l <= 16)
	//        {
	//            ex *= d[3][iat] * d[4][iat];
	//            lam = 3;
	//            if (l == 10)
	//                m = 0;
	//            else if (l == 11)
	//                m = 1;
	//            else if (l == 12)
	//                m = -1;
	//            else if (l == 13)
	//                m = 2;
	//            else if (l == 14)
	//                m = -2;
	//            else if (l == 15)
	//                m = 3;
	//            else if (l == 16)
	//                m = -3;
	//        }
	//        else if (l <= 25)
	//        {
	//            ex *= d[3][iat] * d[3][iat];
	//            lam = 4;
	//            if (l == 17)
	//                m = 0;
	//            else if (l == 18)
	//                m = 1;
	//            else if (l == 19)
	//                m = -1;
	//            else if (l == 20)
	//                m = 2;
	//            else if (l == 21)
	//                m = -2;
	//            else if (l == 22)
	//                m = 3;
	//            else if (l == 23)
	//                m = -3;
	//            else if (l == 24)
	//                m = 4;
	//            else if (l == 25)
	//                m = -4;
	//        }
	//        else if (l <= 36)
	//        {
	//            ex *= pow(d[4][iat], 5);
	//            lam = 5;
	//            if (l == 26)
	//                m = 0;
	//            else if (l == 27)
	//                m = 1;
	//            else if (l == 28)
	//                m = -1;
	//            else if (l == 29)
	//                m = 2;
	//            else if (l == 30)
	//                m = -2;
	//            else if (l == 31)
	//                m = 3;
	//            else if (l == 32)
	//                m = -3;
	//            else if (l == 33)
	//                m = 4;
	//            else if (l == 34)
	//                m = -4;
	//            else if (l == 35)
	//                m = 5;
	//            else if (l == 36)
	//                m = -5;
	//        }
	//    }
	//    double d_t[]{d[0][iat], d[1][iat], d[2][iat], d[3][iat], d[4][iat]};
	//    double SH = spherical_harmonic(lam, m, d_t);
	//    SH *= ex; // multiply radial part with spherical harmonic
	//    auto run = phi.data();
	//    auto run2 = MOs.data();
	//    for (mo = 0; mo < phi.size(); mo++)
	//    {
	//        *run += (*run2).get_coefficient_f(j) * SH; // build MO values at this point
	//        run++, run2++;
	//    }
	//}

	//auto run = phi.data();
	//auto run2 = MOs.data();
	//for (mo = 0; mo < phi.size(); mo++)
	//{
	//    Rho += (*run2).get_occ() * pow(*run, 2);
	//    run++, run2++;
	//}

	//return Rho;
}

//Transposed MO coefficients, [primitive * nmo + mo], built once and reused by every point.
//The first caller is usually several threads at once: make_chi parallelises over atom
//pairs, and compute_dens underneath it is what asks for this. Unlocked, two threads both
//find it cold and one reallocates under the other. The valid flag is read through an
//atomic_ref so the warm path takes no lock (it is called per grid point from every
//density evaluator); the member stays a plain bool so WFN stays copyable.
const double* WFN::get_coef_primitive_major() const
{
	const int _nmo = get_nmo(false);
	if (_nmo <= 0 || nex <= 0) return nullptr;
	std::atomic_ref<bool> valid(coef_primitive_major_valid);
	if (valid.load(std::memory_order_acquire)
		&& coef_primitive_major.size() == (size_t)nex * (size_t)_nmo)
		return coef_primitive_major.data();
	static std::mutex coef_cache_mutex;
	std::lock_guard<std::mutex> lock(coef_cache_mutex);
	if (valid.load(std::memory_order_acquire)
		&& coef_primitive_major.size() == (size_t)nex * (size_t)_nmo)
		return coef_primitive_major.data();
	valid.store(false, std::memory_order_release);
	coef_primitive_major.assign((size_t)nex * (size_t)_nmo, 0.0);
	for (int mo = 0; mo < _nmo; mo++) {
		const double* src = MOs[mo].get_coefficient_ptr();
		if (!src) continue;
		for (int j = 0; j < nex; j++)
			coef_primitive_major[(size_t)j * _nmo + mo] = src[j];
	}
	build_exp_groups();
	valid.store(true, std::memory_order_release);
	return coef_primitive_major.data();
}

void WFN::build_exp_groups() const
{
	std::vector<vec> per_center(ncen);
	ivec local(nex);
	for (int j = 0; j < nex; j++) {
		vec &e = per_center[centers[j] - 1];
		const size_t g = std::find(e.begin(), e.end(), exponents[j]) - e.begin();
		if (g == e.size()) e.push_back(exponents[j]);
		local[j] = (int)g;
	}
	group_exponent.clear();
	center_group_start.assign((size_t)ncen + 1, 0);
	center_min_exponent.assign(ncen, 0.0);
	for (int c = 0; c < ncen; c++) {
		center_group_start[c] = (int)group_exponent.size();
		group_exponent.insert(group_exponent.end(), per_center[c].begin(), per_center[c].end());
		if (!per_center[c].empty()) center_min_exponent[c] = *std::min_element(per_center[c].begin(), per_center[c].end());
	}
	center_group_start[ncen] = (int)group_exponent.size();
	prim_exp_group.resize(nex);
	for (int j = 0; j < nex; j++) prim_exp_group[j] = center_group_start[centers[j] - 1] + local[j];
}

const void WFN::computeValues(
	const d3 &PosGrid, // [3] vector with current position on te grid
	double &Rho,           // Value of Electron Density
	double &normGrad,      // Gradiant Vector
	double *Hess,          // Hessian Matrix, later used to determine lambda2
	double &Elf,           // Value of the ELF
	double &Eli,           // Value of the ELI
	double &Lap            // Value for the Laplacian
) const
{
	d3 Grad;
	double tau;
	computeValues(PosGrid, Rho, Grad, Hess, tau);
	Elf = 0;
	if (Rho > 0)
	{
		normGrad = constants::alpha_coef * sqrt(Grad[0] * Grad[0] + Grad[1] * Grad[1] + Grad[2] * Grad[2]) / pow(Rho, constants::c_43);
		Elf = 1 / (1 + pow(constants::ctelf * pow(Rho, constants::c_m53) * (tau * 0.5 - 0.125 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2)) / Rho), 2));
		Eli = 0.5 * Rho * pow(48 / (Rho * tau - 0.25 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2))), constants::c_38);
	}
	Lap = Hess[0] + Hess[4] + Hess[8];
};

void WFN::computeValues(const d3 &PosGrid, double &Rho, d3 &Grad, double *Hess, double &tau) const
{
	const int _nmo = get_nmo(false);
	vec phi(10 * _nmo, 0.0);
	double *phi_temp;
	double chi[10]{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	double d[4]{ 0, 0, 0, 0 };
	int iat = 0, k, j;
	int l[3]{ 0, 0, 0 };
	double ex = 0;
	double xl[3][3]{ {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
	Grad = { 0, 0, 0 };
	tau = 0;

	Rho = 0;
	Hess[0] = 0;
	Hess[1] = 0;
	Hess[2] = 0;
	Hess[8] = 0;
	Hess[4] = 0;
	Hess[5] = 0;

	const MO *MOs_data = MOs.data();
	//Primitive-major coefficients: every MO for one primitive is contiguous here, where
	//MOs keeps them nex apart. Identical arithmetic in identical order - one sequential
	//read per primitive instead of nmo scattered ones, and no per-point heap allocation.
	const double *const coefs = get_coef_primitive_major();

	for (j = 0; j < nex; j++)
	{
		iat = get_center(j) - 1;

		constants::type2vector(get_type(j), l);
		d[0] = PosGrid[0] - atoms[iat].get_coordinate(0);
		d[1] = PosGrid[1] - atoms[iat].get_coordinate(1);
		d[2] = PosGrid[2] - atoms[iat].get_coordinate(2);
		d[3] = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
		double temp = -get_exponent(j) * (d[3]);
		if (temp < constants::exp_cutoff)
			continue;
		ex = exp(temp);
		for (k = 0; k < 3; k++)
		{
			if (l[k] == 0)
			{
				xl[k][0] = 1.0;
				xl[k][1] = 0.0;
				xl[k][2] = 0.0;
			}
			else if (l[k] == 1)
			{
				xl[k][0] = d[k];
				xl[k][1] = 1.0;
				xl[k][2] = 0.0;
			}
			else if (l[k] == 2)
			{
				xl[k][0] = d[k] * d[k];
				xl[k][1] = 2 * d[k];
				xl[k][2] = 2;
			}
			else if (l[k] == 3)
			{
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d[k];
				xl[k][1] = 3 * d2;
				xl[k][2] = 6 * d[k];
			}
			else if (l[k] == 4)
			{
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d2;
				xl[k][1] = 4 * d2 * d[k];
				xl[k][2] = 12 * d2;
			}
			else if (l[k] == 5)
			{
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d2 * d[k];
				xl[k][1] = 5 * d2 * d2;
				xl[k][2] = 20 * d2 * d[k];
			}
			else if (l[k] == 6)
			{
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d2 * d2;
				xl[k][1] = 6 * d2 * d2 * d[k];
				xl[k][2] = 30 * d2 * d2;
			}
			else
			{
				return;
			}
		}
		const double ex2 = 2 * get_exponent(j);
		chi[0] = xl[0][0] * xl[1][0] * xl[2][0] * ex;
		chi[1] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * xl[1][0] * xl[2][0] * ex;
		chi[2] = (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[0][0] * xl[2][0] * ex;
		chi[3] = (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * xl[0][0] * xl[1][0] * ex;
		const double temp_ex = pow(ex2, 2);
		chi[4] = (xl[0][2] - ex2 * (2 * l[0] + 1) * xl[0][0] + temp_ex * pow(d[0], l[0] + 2)) * xl[1][0] * xl[2][0] * ex;
		chi[5] = (xl[1][2] - ex2 * (2 * l[1] + 1) * xl[1][0] + temp_ex * pow(d[1], l[1] + 2)) * xl[2][0] * xl[0][0] * ex;
		chi[6] = (xl[2][2] - ex2 * (2 * l[2] + 1) * xl[2][0] + temp_ex * pow(d[2], l[2] + 2)) * xl[0][0] * xl[1][0] * ex;
		chi[7] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[2][0] * ex;
		chi[8] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * xl[1][0] * ex;
		chi[9] = (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[0][0] * ex;

		const double *c_row = coefs + (size_t)j * nmo;
		double *phi_ptr = phi.data();
		for (int mo = 0; mo < nmo; ++mo, phi_ptr += 10)
		{
			const double c = c_row[mo];
			for (k = 0; k < 10; k++)
				phi_ptr[k] += c * chi[k];
		}
	}
	for (int mo = 0; mo < _nmo; mo++)
	{
		const double occ = get_MO_occ(mo);
		const double docc = 2 * occ;
		if (occ != 0)
		{
			phi_temp = &phi[mo * 10];
			Rho += occ * pow(*phi_temp, 2);
			Grad[0] += docc * *phi_temp * phi_temp[1];
			Grad[1] += docc * *phi_temp * phi_temp[2];
			Grad[2] += docc * *phi_temp * phi_temp[3];
			Hess[0] += docc * (*phi_temp * phi_temp[4] + pow(phi_temp[1], 2));
			Hess[4] += docc * (*phi_temp * phi_temp[5] + pow(phi_temp[2], 2));
			Hess[8] += docc * (*phi_temp * phi_temp[6] + pow(phi_temp[3], 2));
			Hess[1] += docc * (*phi_temp * phi_temp[7] + phi_temp[1] * phi_temp[2]);
			Hess[2] += docc * (*phi_temp * phi_temp[8] + phi_temp[1] * phi_temp[3]);
			Hess[5] += docc * (*phi_temp * phi_temp[9] + phi_temp[2] * phi_temp[3]);
			tau += occ * (pow(phi_temp[1], 2) + pow(phi_temp[2], 2) + pow(phi_temp[3], 2));
		}
	}

	Hess[3] = Hess[1];
	Hess[6] = Hess[2];
	Hess[7] = Hess[5];
};

const void WFN::computeELIELF(
	const d3 &PosGrid, // [3] vector with current position on te grid
	double &Elf,           // Value of the ELF
	double &Eli            // Value of the ELI
) const
{
	const int _nmo = get_nmo(false);
	vec phi(4 * _nmo, 0.0);
	double *phi_temp;
	double chi[4]{ 0, 0, 0, 0 };
	double d[3]{ 0, 0, 0 };
	int iat = 0, k, j;
	int l[3]{ 0, 0, 0 };
	double ex = 0;
	double xl[3][3]{ {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };

	const MO *MOs_data = MOs.data();
	//Primitive-major coefficients: every MO for one primitive is contiguous here, where
	//MOs keeps them nex apart. Identical arithmetic in identical order - one sequential
	//read per primitive instead of nmo scattered ones, and no per-point heap allocation.
	const double *const coefs = get_coef_primitive_major();

	for (j = 0; j < nex; j++)
	{
		iat = centers[j] - 1;

		constants::type2vector(get_type(j), l);
		d[0] = PosGrid[0] - atoms[iat].get_coordinate(0);
		d[1] = PosGrid[1] - atoms[iat].get_coordinate(1);
		d[2] = PosGrid[2] - atoms[iat].get_coordinate(2);
		double temp = -get_exponent(j) * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
		if (temp < constants::exp_cutoff)
			continue;
		ex = exp(temp);
		for (int k = 0; k < 3; k++)
		{
			if (l[k] == 0)
			{
				xl[k][0] = 1.0;
				xl[k][1] = 0.0;
				xl[k][2] = 0.0;
			}
			else if (l[k] == 1)
			{
				xl[k][0] = d[k];
				xl[k][1] = 1.0;
				xl[k][2] = 0.0;
			}
			else if (l[k] == 2)
			{
				xl[k][0] = d[k] * d[k];
				xl[k][1] = 2 * d[k];
				xl[k][2] = 2;
			}
			else if (l[k] == 3)
			{
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d[k];
				xl[k][1] = 3 * d2;
				xl[k][2] = 6 * d[k];
			}
			else if (l[k] == 4)
			{
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d2;
				xl[k][1] = 4 * d2 * d[k];
				xl[k][2] = 12 * d2;
			}
			else if (l[k] == 5)
			{
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d2 * d[k];
				xl[k][1] = 5 * d2 * d2;
				xl[k][2] = 20 * d2 * d[k];
			}
			else if (l[k] == 6)
			{
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d2 * d2;
				xl[k][1] = 6 * d2 * d2 * d[k];
				xl[k][2] = 30 * d2 * d2;
			}
			else
			{
				return;
			}
		}
		double ex2 = 2 * get_exponent(j);
		chi[0] = xl[0][0] * xl[1][0] * xl[2][0] * ex;
		chi[1] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * xl[1][0] * xl[2][0] * ex;
		chi[2] = (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[0][0] * xl[2][0] * ex;
		chi[3] = (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * xl[0][0] * xl[1][0] * ex;

		const double *c_row = coefs + (size_t)j * nmo;
		double *phi_ptr = phi.data();
		for (int mo = 0; mo < nmo; ++mo, phi_ptr += 4)
		{
			const double c = c_row[mo];
			for (k = 0; k < 4; k++)
				phi_ptr[k] += c * chi[k];
		}
	}

	double Grad[3]{ 0, 0, 0 };
	double tau = 0;
	double Rho = 0;

	for (int mo = 0; mo < _nmo; mo++)
	{
		const double occ = get_MO_occ(mo);
		const double docc = 2 * occ;
		if (occ != 0)
		{
			phi_temp = &phi[mo * 4];
			Rho += occ * pow(*phi_temp, 2);
			Grad[0] += docc * *phi_temp * phi_temp[1];
			Grad[1] += docc * *phi_temp * phi_temp[2];
			Grad[2] += docc * *phi_temp * phi_temp[3];
			tau += occ * (pow(phi_temp[1], 2) + pow(phi_temp[2], 2) + pow(phi_temp[3], 2));
		}
	}
	if (Rho > 0)
	{
		Elf = 1 / (1 + pow(constants::ctelf * pow(Rho, constants::c_m53) * (tau * 0.5 - 0.125 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2)) / Rho), 2));
		Eli = 0.5 * Rho * pow(48 / (Rho * tau - 0.25 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2))), constants::c_38);
	}
};

const double WFN::computeELI(
	const d3 &PosGrid // [3] vector with current position on te grid
) const
{
	const int _nmo = get_nmo(false);
	vec phi(4 * _nmo, 0.0);
	double *phi_temp;
	double chi[4]{ 0, 0, 0, 0 };
	double d[3]{ 0, 0, 0 };
	int iat = 0, k, j;
	int l[3]{ 0, 0, 0 };
	double ex = 0;
	double xl[3][3]{ {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };

	const MO *MOs_data = MOs.data();
	//Primitive-major coefficients: every MO for one primitive is contiguous here, where
	//MOs keeps them nex apart. Identical arithmetic in identical order - one sequential
	//read per primitive instead of nmo scattered ones, and no per-point heap allocation.
	const double *const coefs = get_coef_primitive_major();

	for (j = 0; j < nex; j++)
	{
		iat = centers[j] - 1;

		constants::type2vector(types[j], l);
		d[0] = PosGrid[0] - atoms[iat].get_coordinate(0);
		d[1] = PosGrid[1] - atoms[iat].get_coordinate(1);
		d[2] = PosGrid[2] - atoms[iat].get_coordinate(2);
		double temp = -exponents[j] * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
		if (temp < constants::exp_cutoff)
			continue;
		ex = exp(temp);
		for (int k = 0; k < 3; k++)
		{
			if (l[k] == 0)
			{
				xl[k][0] = 1.0;
				xl[k][1] = 0.0;
				xl[k][2] = 0.0;
			}
			else if (l[k] == 1)
			{
				xl[k][0] = d[k];
				xl[k][1] = 1.0;
				xl[k][2] = 0.0;
			}
			else if (l[k] == 2)
			{
				xl[k][0] = d[k] * d[k];
				xl[k][1] = 2 * d[k];
				xl[k][2] = 2;
			}
			else if (l[k] == 3)
			{
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d[k];
				xl[k][1] = 3 * d2;
				xl[k][2] = 6 * d[k];
			}
			else if (l[k] == 4)
			{
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d2;
				xl[k][1] = 4 * d2 * d[k];
				xl[k][2] = 12 * d2;
			}
			else if (l[k] == 5)
			{
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d2 * d[k];
				xl[k][1] = 5 * d2 * d2;
				xl[k][2] = 20 * d2 * d[k];
			}
			else if (l[k] == 6)
			{
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d2 * d2;
				xl[k][1] = 6 * d2 * d2 * d[k];
				xl[k][2] = 30 * d2 * d2;
			}
			else
			{
				err_not_impl_f("Higher angular momentum of cartesian function in ELI computation", std::cout);
			}
		}
		double ex2 = 2 * get_exponent(j);
		chi[0] = xl[0][0] * xl[1][0] * xl[2][0] * ex;
		chi[1] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * xl[1][0] * xl[2][0] * ex;
		chi[2] = (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[0][0] * xl[2][0] * ex;
		chi[3] = (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * xl[0][0] * xl[1][0] * ex;

		const double *c_row = coefs + (size_t)j * nmo;
		double *phi_ptr = phi.data();
		for (int mo = 0; mo < nmo; ++mo, phi_ptr += 4)
		{
			const double c = c_row[mo];
			for (k = 0; k < 4; k++)
				phi_ptr[k] += c * chi[k];
		}
	}

	double Grad[3]{ 0, 0, 0 };
	double tau = 0;
	double Rho = 0;

	for (int mo = 0; mo < _nmo; mo++)
	{
		const double occ = get_MO_occ(mo);
		const double docc = 2 * occ;
		if (occ != 0)
		{
			phi_temp = &phi[mo * 4];
			Rho += occ * pow(*phi_temp, 2);
			Grad[0] += docc * *phi_temp * phi_temp[1];
			Grad[1] += docc * *phi_temp * phi_temp[2];
			Grad[2] += docc * *phi_temp * phi_temp[3];
			tau += occ * (pow(phi_temp[1], 2) + pow(phi_temp[2], 2) + pow(phi_temp[3], 2));
		}
	}
	return 0.5 * Rho * pow(48 / (Rho * tau - 0.25 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2))), constants::c_38);
};

void WFN::computeRhoELI(
	const d3 &PosGrid, // [3] vector with current position on te grid
	double& out_Rho,
	double& out_Eli
) const
{
	const int _nmo = get_nmo(false);
	//Called per point of the basin quadrature (the ELI climb uses computeELIGrad);
	//the two buffers were a third of its allocator time, so they persist per thread
	thread_local vec phi, d;
	phi.assign(4 * _nmo, 0.0);
	if (d.size() < 16 * (size_t)ncen) d.resize(16 * (size_t)ncen);
	double *phi_temp;
	double chi[4]{ 0, 0, 0, 0 };
	int k, j;
	double ex = 0;

	for (j = 0; j < ncen; j++)
	{
		const atom &a = atoms[j];
		double *d_ = d.data() + 16 * j;
		d_[0] = PosGrid[0] - a.get_coordinate(0);
		d_[1] = PosGrid[1] - a.get_coordinate(1);
		d_[2] = PosGrid[2] - a.get_coordinate(2);
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
	//Primitive-major coefficients: every MO for one primitive is contiguous here, where
	//MOs keeps them nex apart. Identical arithmetic in identical order - one sequential
	//read per primitive instead of nmo scattered ones, and no per-point heap allocation.
	const double *const coefs = get_coef_primitive_major();
	thread_local vec exps;
	if (exps.size() < group_exponent.size()) exps.resize(group_exponent.size());
	exp_table([r2 = d.data()](const int c) { return r2[16 * (size_t)c + 3]; }, exps.data());
	const int *group = prim_exp_group.data();

	for (j = 0; j < nex; j++)
	{
		ex = exps[group[j]];
		if (ex == 0.0)
			continue;
		const double *d_ = d.data() + 16 * (centers_data[j] - 1);
		const int type = types_data[j];
		const int type_index = (type - 1) * 3;
		int lx = 0;
		int ly = 0;
		int lz = 0;
		if (type > 1 && type <= 286) {
			lx = constants::type_vector[type_index];
			ly = constants::type_vector[type_index + 1];
			lz = constants::type_vector[type_index + 2];
		}

		double x0 = 1.0, x1 = 0.0, xnext = d_[0];
		double y0 = 1.0, y1 = 0.0, ynext = d_[1];
		double z0 = 1.0, z1 = 0.0, znext = d_[2];

		switch (lx) {
		case 0: x0 = 1.0;   x1 = 0.0;    xnext = d_[0];  break;
		case 1: x0 = d_[0]; x1 = 1.0;    xnext = d_[4];  break;
		case 2: x0 = d_[4]; x1 = 2*d_[0]; xnext = d_[7];  break;
		case 3: x0 = d_[7]; x1 = 3*d_[4]; xnext = d_[10]; break;
		case 4: x0 = d_[10]; x1 = 4*d_[7]; xnext = d_[13]; break;
		case 5: x0 = d_[13]; x1 = 5*d_[10]; xnext = d_[13]*d_[0]; break;
		case 6: x0 = d_[13]*d_[0]; x1 = 6*d_[13]; xnext = d_[13]*d_[4]; break;
		case 7: x0 = d_[13] * d_[4]; x1 = 7*d_[13] * d_[0]; xnext = d_[13] * d_[7]; break;
		case 8: x0 = d_[13] * d_[7]; x1 = 8*d_[13] * d_[4]; xnext = d_[13] * d_[10]; break;
		case 9: x0 = d_[13] * d_[10]; x1 = 9*d_[13] * d_[7]; xnext = d_[13] * d_[13]; break;
		case 10: x0 = d_[13] * d_[13]; x1 = 10*d_[13] * d_[10]; xnext = d_[13] * d_[13] * d_[0]; break;
		default:
			err_not_impl_f("Higher angular momentum of cartesian function in ELI computation", std::cout);
			break;
		}
		switch (ly) {
		case 0: y0 = 1.0;   y1 = 0.0;    ynext = d_[1];  break;
		case 1: y0 = d_[1]; y1 = 1.0;    ynext = d_[5];  break;
		case 2: y0 = d_[5]; y1 = 2*d_[1]; ynext = d_[8];  break;
		case 3: y0 = d_[8]; y1 = 3*d_[5]; ynext = d_[11]; break;
		case 4: y0 = d_[11]; y1 = 4*d_[8]; ynext = d_[14]; break;
		case 5: y0 = d_[14]; y1 = 5*d_[11]; ynext = d_[14]*d_[1]; break;
		case 6: y0 = d_[14]*d_[1]; y1 = 6*d_[14]; ynext = d_[14]*d_[5]; break;
		case 7: y0 = d_[14] * d_[5]; y1 = 7*d_[14] * d_[1]; ynext = d_[14] * d_[8]; break;
		case 8: y0 = d_[14] * d_[8]; y1 = 8*d_[14] * d_[5]; ynext = d_[14] * d_[11]; break;
		case 9: y0 = d_[14] * d_[11]; y1 = 9*d_[14] * d_[8]; ynext = d_[14] * d_[14]; break;
		case 10: y0 = d_[14] * d_[14]; y1 = 10*d_[14] * d_[11]; ynext = d_[14] * d_[14] * d_[1]; break;
		default:
			err_not_impl_f("Higher angular momentum of cartesian function in ELI computation", std::cout);
			break;
		}
		switch (lz) {
		case 0: z0 = 1.0;   z1 = 0.0;    znext = d_[2];  break;
		case 1: z0 = d_[2]; z1 = 1.0;    znext = d_[6];  break;
		case 2: z0 = d_[6]; z1 = 2*d_[2]; znext = d_[9];  break;
		case 3: z0 = d_[9]; z1 = 3*d_[6]; znext = d_[12]; break;
		case 4: z0 = d_[12]; z1 = 4*d_[9]; znext = d_[15]; break;
		case 5: z0 = d_[15]; z1 = 5*d_[12]; znext = d_[15]*d_[2]; break;
		case 6: z0 = d_[15]*d_[2]; z1 = 6*d_[15]; znext = d_[15]*d_[6]; break;
		case 7: z0 = d_[15] * d_[6]; z1 = 7*d_[15] * d_[2]; znext = d_[15] * d_[9]; break;
		case 8: z0 = d_[15] * d_[9]; z1 = 8*d_[15] * d_[6]; znext = d_[15] * d_[12]; break;
		case 9: z0 = d_[15] * d_[12]; z1 = 9*d_[15] * d_[9]; znext = d_[15] * d_[15]; break;
		case 10: z0 = d_[15] * d_[15]; z1 = 10*d_[15] * d_[12]; znext = d_[15] * d_[15] * d_[2]; break;
		default:
			err_not_impl_f("Higher angular momentum of cartesian function in ELI computation", std::cout);
			break;
		}

		const double ex2 = 2 * exponents_data[j];
		chi[0] = x0 * y0 * z0 * ex;
		chi[1] = (x1 - ex2 * xnext) * y0 * z0 * ex;
		chi[2] = (y1 - ex2 * ynext) * x0 * z0 * ex;
		chi[3] = (z1 - ex2 * znext) * x0 * y0 * ex;

		const double *c_row = coefs + (size_t)j * _nmo;
		double *phi_ptr = phi.data();
		for (int mo = 0; mo < _nmo; ++mo, phi_ptr += 4)
		{
			const double c = c_row[mo];
			for (k = 0; k < 4; k++)
				phi_ptr[k] += c * chi[k];
		}
	}

	double Grad[3]{ 0, 0, 0 };
	double tau = 0;
	double Rho = 0;

	for (int mo = 0; mo < _nmo; mo++)
	{
		const double occ = get_MO_occ(mo);
		const double docc = 2 * occ;
		if (occ != 0)
		{
			phi_temp = &phi[mo * 4];
			Rho += occ * pow(*phi_temp, 2);
			Grad[0] += docc * *phi_temp * phi_temp[1];
			Grad[1] += docc * *phi_temp * phi_temp[2];
			Grad[2] += docc * *phi_temp * phi_temp[3];
			tau += occ * (pow(phi_temp[1], 2) + pow(phi_temp[2], 2) + pow(phi_temp[3], 2));
		}
	}
	//ELI-D of one spin channel of a closed shell, Kohout's definition and DGrid's alpha-alpha
	//field: rho_s (12 / g_s)^(3/8) with g_s = rho_s tau_s - |grad rho_s|^2 / 4 and every
	//sigma quantity half the total, which is where the 1/2 and the 48 come from.
	//g vanishes where a single orbital carries the density (also where the exponent cutoff
	//has dropped every other one); ELI-D has no finite value there and the grid holds 0,
	//as the basin climb in computeELIGrad does, instead of an inf or NaN no cube reader parses
	const double g = Rho * tau - 0.25 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2));
	out_Eli = g > 0 ? 0.5 * Rho * pow(48 / g, constants::c_38) : 0.0;
	out_Rho = Rho;
};

void WFN::computeELIGrad(
	const d3 &PosGrid,
	double& out_Eli,
	d3& out_grad
) const
{
	//ELI-D Y = rho/2 (48/g)^(3/8) with g = rho tau - |grad rho|^2 / 4, so
	//grad Y = (48/g)^(3/8) / 2 (grad rho - 3/8 rho grad g / g) with
	//grad g = tau grad rho + rho grad tau - H grad rho / 2, H the density Hessian.
	//One pass with orbital values, gradients and Hessians replaces the six ELI
	//evaluations of the central difference the basin climb used before
	const int _nmo = get_nmo(false);
	thread_local vec phi, d;
	phi.assign(10 * _nmo, 0.0);
	if (d.size() < 4 * (size_t)ncen) d.resize(4 * (size_t)ncen);
	for (int j = 0; j < ncen; j++)
	{
		double *d_ = d.data() + 4 * j;
		for (int k = 0; k < 3; k++) d_[k] = PosGrid[k] - atoms[j].get_coordinate(k);
		d_[3] = d_[0] * d_[0] + d_[1] * d_[1] + d_[2] * d_[2];
	}
	const double *const coefs = get_coef_primitive_major();
	thread_local vec exps;
	if (exps.size() < group_exponent.size()) exps.resize(group_exponent.size());
	exp_table([r2 = d.data()](const int c) { return r2[4 * (size_t)c + 3]; }, exps.data());
	const int *group = prim_exp_group.data();
	for (int j = 0; j < nex; j++)
	{
		const double ex = exps[group[j]];
		if (ex == 0.0)
			continue;
		const double *d_ = d.data() + 4 * (centers[j] - 1);
		const double ex2 = 2 * exponents[j];
		int l[3]{ 0, 0, 0 };
		constants::type2vector(types[j], l);
		if (l[0] < 0) continue; //type outside 1..286 (l > 10)
		//per axis f = x^l e^(-a x^2) without the exponential: f, f', f''
		double f[3], f1[3], f2[3];
		for (int k = 0; k < 3; k++)
		{
			double p[13];
			p[0] = 1.0;
			for (int n = 1; n <= l[k] + 2; n++) p[n] = p[n - 1] * d_[k];
			f[k] = p[l[k]];
			f1[k] = (l[k] ? l[k] * p[l[k] - 1] : 0.0) - ex2 * p[l[k] + 1];
			f2[k] = (l[k] > 1 ? l[k] * (l[k] - 1) * p[l[k] - 2] : 0.0) - ex2 * (2 * l[k] + 1) * p[l[k]] + ex2 * ex2 * p[l[k] + 2];
		}
		//0 value, 1-3 x y z, 4-6 xx yy zz, 7 xy, 8 xz, 9 yz
		const double chi[10] = {
			f[0] * f[1] * f[2] * ex,
			f1[0] * f[1] * f[2] * ex, f[0] * f1[1] * f[2] * ex, f[0] * f[1] * f1[2] * ex,
			f2[0] * f[1] * f[2] * ex, f[0] * f2[1] * f[2] * ex, f[0] * f[1] * f2[2] * ex,
			f1[0] * f1[1] * f[2] * ex, f1[0] * f[1] * f1[2] * ex, f[0] * f1[1] * f1[2] * ex };
		const double *c_row = coefs + (size_t)j * _nmo;
		double *phi_ptr = phi.data();
		for (int mo = 0; mo < _nmo; ++mo, phi_ptr += 10)
		{
			const double c = c_row[mo];
			for (int k = 0; k < 10; k++)
				phi_ptr[k] += c * chi[k];
		}
	}
	static constexpr int hidx[3][3] = { {4, 7, 8}, {7, 5, 9}, {8, 9, 6} };
	double rho = 0, tau = 0, G[3]{ 0, 0, 0 }, T[3]{ 0, 0, 0 }, H[3][3]{ {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
	for (int mo = 0; mo < _nmo; mo++)
	{
		const double occ = get_MO_occ(mo);
		if (occ == 0) continue;
		const double *p = &phi[mo * 10], docc = 2 * occ;
		rho += occ * p[0] * p[0];
		for (int i = 0; i < 3; i++)
		{
			tau += occ * p[1 + i] * p[1 + i];
			G[i] += docc * p[0] * p[1 + i];
			for (int k = 0; k < 3; k++)
			{
				H[i][k] += docc * (p[0] * p[hidx[i][k]] + p[1 + i] * p[1 + k]);
				T[k] += docc * p[1 + i] * p[hidx[i][k]];
			}
		}
	}
	const double g = rho * tau - 0.25 * (G[0] * G[0] + G[1] * G[1] + G[2] * G[2]);
	if (!(g > 0))
	{
		//a single occupied orbital: g vanishes and ELI-D has no finite value there
		out_Eli = 0;
		out_grad = { 0, 0, 0 };
		return;
	}
	const double f = pow(48 / g, constants::c_38);
	out_Eli = 0.5 * rho * f;
	for (int k = 0; k < 3; k++)
	{
		const double dg = G[k] * tau + rho * T[k] - 0.5 * (G[0] * H[0][k] + G[1] * H[1][k] + G[2] * H[2][k]);
		out_grad[k] = 0.5 * f * (G[k] - 0.375 * rho * dg / g);
	}
};

void WFN::computeGrad(
	const d3 &PosGrid, // [3] vector with current position on te grid
	d3& gradient,
	double *rho
) const
{
	const int _nmo = get_nmo(false);
	//Called per point of the basin quadrature (the ELI climb uses computeELIGrad);
	//the two buffers were a third of its allocator time, so they persist per thread
	thread_local vec phi, d;
	phi.assign(4 * _nmo, 0.0);
	if (d.size() < 16 * (size_t)ncen) d.resize(16 * (size_t)ncen);
	double *phi_temp;
	double chi[4]{ 0, 0, 0, 0 };
	int k, j;
	double ex = 0;

	for (j = 0; j < ncen; j++)
	{
		const atom &a = atoms[j];
		double *d_ = d.data() + 16 * j;
		d_[0] = PosGrid[0] - a.get_coordinate(0);
		d_[1] = PosGrid[1] - a.get_coordinate(1);
		d_[2] = PosGrid[2] - a.get_coordinate(2);
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
	//Primitive-major coefficients: every MO for one primitive is contiguous here, where
	//MOs keeps them nex apart. Identical arithmetic in identical order - one sequential
	//read per primitive instead of nmo scattered ones, and no per-point heap allocation.
	const double *const coefs = get_coef_primitive_major();
	thread_local vec exps;
	if (exps.size() < group_exponent.size()) exps.resize(group_exponent.size());
	exp_table([r2 = d.data()](const int c) { return r2[16 * (size_t)c + 3]; }, exps.data());
	const int *group = prim_exp_group.data();

	for (j = 0; j < nex; j++)
	{
		ex = exps[group[j]];
		if (ex == 0.0)
			continue;
		const double *d_ = d.data() + 16 * (centers_data[j] - 1);
		const int type = types_data[j];
		const int type_index = (type - 1) * 3;
		int lx = 0;
		int ly = 0;
		int lz = 0;
		if (type > 1 && type <= 286) {
			lx = constants::type_vector[type_index];
			ly = constants::type_vector[type_index + 1];
			lz = constants::type_vector[type_index + 2];
		}

		double x0 = 1.0, x1 = 0.0, xnext = d_[0];
		double y0 = 1.0, y1 = 0.0, ynext = d_[1];
		double z0 = 1.0, z1 = 0.0, znext = d_[2];

		switch (lx) {
		case 0: x0 = 1.0;   x1 = 0.0;    xnext = d_[0];  break;
		case 1: x0 = d_[0]; x1 = 1.0;    xnext = d_[4];  break;
		case 2: x0 = d_[4]; x1 = 2*d_[0]; xnext = d_[7];  break;
		case 3: x0 = d_[7]; x1 = 3*d_[4]; xnext = d_[10]; break;
		case 4: x0 = d_[10]; x1 = 4*d_[7]; xnext = d_[13]; break;
		case 5: x0 = d_[13]; x1 = 5*d_[10]; xnext = d_[13]*d_[0]; break;
		case 6: x0 = d_[13]*d_[0]; x1 = 6*d_[13]; xnext = d_[13]*d_[4]; break;
		case 7: x0 = d_[13] * d_[4]; x1 = 7*d_[13] * d_[0]; xnext = d_[13] * d_[7]; break;
		case 8: x0 = d_[13] * d_[7]; x1 = 8*d_[13] * d_[4]; xnext = d_[13] * d_[10]; break;
		case 9: x0 = d_[13] * d_[10]; x1 = 9*d_[13] * d_[7]; xnext = d_[13] * d_[13]; break;
		case 10: x0 = d_[13] * d_[13]; x1 = 10*d_[13] * d_[10]; xnext = d_[13] * d_[13] * d_[0]; break;
		default:
			err_not_impl_f("Higher angular momentum of cartesian function in ELI computation", std::cout);
			break;
		}
		switch (ly) {
		case 0: y0 = 1.0;   y1 = 0.0;    ynext = d_[1];  break;
		case 1: y0 = d_[1]; y1 = 1.0;    ynext = d_[5];  break;
		case 2: y0 = d_[5]; y1 = 2*d_[1]; ynext = d_[8];  break;
		case 3: y0 = d_[8]; y1 = 3*d_[5]; ynext = d_[11]; break;
		case 4: y0 = d_[11]; y1 = 4*d_[8]; ynext = d_[14]; break;
		case 5: y0 = d_[14]; y1 = 5*d_[11]; ynext = d_[14]*d_[1]; break;
		case 6: y0 = d_[14]*d_[1]; y1 = 6*d_[14]; ynext = d_[14]*d_[5]; break;
		case 7: y0 = d_[14] * d_[5]; y1 = 7*d_[14] * d_[1]; ynext = d_[14] * d_[8]; break;
		case 8: y0 = d_[14] * d_[8]; y1 = 8*d_[14] * d_[5]; ynext = d_[14] * d_[11]; break;
		case 9: y0 = d_[14] * d_[11]; y1 = 9*d_[14] * d_[8]; ynext = d_[14] * d_[14]; break;
		case 10: y0 = d_[14] * d_[14]; y1 = 10*d_[14] * d_[11]; ynext = d_[14] * d_[14] * d_[1]; break;
		default:
			err_not_impl_f("Higher angular momentum of cartesian function in ELI computation", std::cout);
			break;
		}
		switch (lz) {
		case 0: z0 = 1.0;   z1 = 0.0;    znext = d_[2];  break;
		case 1: z0 = d_[2]; z1 = 1.0;    znext = d_[6];  break;
		case 2: z0 = d_[6]; z1 = 2*d_[2]; znext = d_[9];  break;
		case 3: z0 = d_[9]; z1 = 3*d_[6]; znext = d_[12]; break;
		case 4: z0 = d_[12]; z1 = 4*d_[9]; znext = d_[15]; break;
		case 5: z0 = d_[15]; z1 = 5*d_[12]; znext = d_[15]*d_[2]; break;
		case 6: z0 = d_[15]*d_[2]; z1 = 6*d_[15]; znext = d_[15]*d_[6]; break;
		case 7: z0 = d_[15] * d_[6]; z1 = 7*d_[15] * d_[2]; znext = d_[15] * d_[9]; break;
		case 8: z0 = d_[15] * d_[9]; z1 = 8*d_[15] * d_[6]; znext = d_[15] * d_[12]; break;
		case 9: z0 = d_[15] * d_[12]; z1 = 9*d_[15] * d_[9]; znext = d_[15] * d_[15]; break;
		case 10: z0 = d_[15] * d_[15]; z1 = 10*d_[15] * d_[12]; znext = d_[15] * d_[15] * d_[2]; break;
		default:
			err_not_impl_f("Higher angular momentum of cartesian function in ELI computation", std::cout);
			break;
		}

		const double ex2 = 2 * exponents_data[j];
		chi[0] = x0 * y0 * z0 * ex;
		chi[1] = (x1 - ex2 * xnext) * y0 * z0 * ex;
		chi[2] = (y1 - ex2 * ynext) * x0 * z0 * ex;
		chi[3] = (z1 - ex2 * znext) * x0 * y0 * ex;

		const double *c_row = coefs + (size_t)j * _nmo;
		double *phi_ptr = phi.data();
		for (int mo = 0; mo < _nmo; ++mo, phi_ptr += 4)
		{
			const double c = c_row[mo];
			for (k = 0; k < 4; k++)
				phi_ptr[k] += c * chi[k];
		}
	}

	double Grad[3]{ 0, 0, 0 }, Rho = 0.0;

	for (int mo = 0; mo < _nmo; mo++)
	{
		const double occ = get_MO_occ(mo);
		const double docc = 2 * occ;
		if (occ != 0)
		{
			phi_temp = &phi[mo * 4];
			Grad[0] += docc * *phi_temp * phi_temp[1];
			Grad[1] += docc * *phi_temp * phi_temp[2];
			Grad[2] += docc * *phi_temp * phi_temp[3];
			if (rho) Rho += occ * *phi_temp * *phi_temp;
		}
	}
	gradient[0] = Grad[0];
	gradient[1] = Grad[1];
	gradient[2] = Grad[2];
	if (rho) *rho = Rho;

};

const double WFN::computeELF(
	const d3 &PosGrid // [3] vector with current position on te grid
) const
{
	const int _nmo = get_nmo(false);
	vec phi(4 * _nmo, 0.0);
	double *phi_temp;
	double chi[4]{ 0, 0, 0, 0 };
	double d[3]{ 0, 0, 0 };
	int iat = 0;
	int l[3]{ 0, 0, 0 };
	double ex = 0;
	double xl[3][3]{ {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };

	for (int j = 0; j < nex; j++)
	{
		iat = get_center(j) - 1;

		constants::type2vector(get_type(j), l);
		d[0] = PosGrid[0] - atoms[iat].get_coordinate(0);
		d[1] = PosGrid[1] - atoms[iat].get_coordinate(1);
		d[2] = PosGrid[2] - atoms[iat].get_coordinate(2);
		double temp = -get_exponent(j) * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
		if (temp < constants::exp_cutoff)
			continue;
		ex = exp(temp);
		for (int k = 0; k < 3; k++)
		{
			if (l[k] == 0)
			{
				xl[k][0] = 1.0;
				xl[k][1] = 0.0;
				xl[k][2] = 0.0;
			}
			else if (l[k] == 1)
			{
				xl[k][0] = d[k];
				xl[k][1] = 1.0;
				xl[k][2] = 0.0;
			}
			else if (l[k] == 2)
			{
				xl[k][0] = d[k] * d[k];
				xl[k][1] = 2 * d[k];
				xl[k][2] = 2;
			}
			else if (l[k] == 3)
			{
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d[k];
				xl[k][1] = 3 * d2;
				xl[k][2] = 6 * d[k];
			}
			else if (l[k] == 4)
			{
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d2;
				xl[k][1] = 4 * d2 * d[k];
				xl[k][2] = 12 * d2;
			}
			else if (l[k] == 5)
			{
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d2 * d[k];
				xl[k][1] = 5 * d2 * d2;
				xl[k][2] = 20 * d2 * d[k];
			}
			else if (l[k] == 6)
			{
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d2 * d2;
				xl[k][1] = 6 * d2 * d2 * d[k];
				xl[k][2] = 30 * d2 * d2;
			}
			else
			{
				err_not_impl_f("Higher angular momentum of cartesian function in ELF computation", std::cout);
			}
		}
		double ex2 = 2 * get_exponent(j);
		chi[0] = xl[0][0] * xl[1][0] * xl[2][0] * ex;
		chi[1] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * xl[1][0] * xl[2][0] * ex;
		chi[2] = (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[0][0] * xl[2][0] * ex;
		chi[3] = (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * xl[0][0] * xl[1][0] * ex;
		for (int mo = 0; mo < _nmo; mo++)
		{
			phi_temp = &phi[mo * 4];
			for (int i = 0; i < 4; i++)
				//                if( abs(chi[i]) * coefficients[nprim*maxc[j]+j] > pow(10.0,-10) )
				phi_temp[i] += MOs[mo].get_coefficient_f(j) * chi[i]; // build MO values at this point
		}
	}

	double Grad[3]{ 0, 0, 0 };
	double tau = 0;
	double Rho = 0;

	for (int mo = 0; mo < _nmo; mo++)
	{
		const double occ = get_MO_occ(mo);
		const double docc = 2 * occ;
		if (occ != 0)
		{
			phi_temp = &phi[mo * 4];
			Rho += occ * pow(*phi_temp, 2);
			Grad[0] += docc * *phi_temp * phi_temp[1];
			Grad[1] += docc * *phi_temp * phi_temp[2];
			Grad[2] += docc * *phi_temp * phi_temp[3];
			tau += occ * (pow(phi_temp[1], 2) + pow(phi_temp[2], 2) + pow(phi_temp[3], 2));
		}
	}
	return 1 / (1 + pow(constants::ctelf * pow(Rho, constants::c_m53) * (tau * 0.5 - 0.125 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2)) / Rho), 2));
};

const void WFN::computeLapELIELF(
	const d3 &PosGrid, // [3] vector with current position on te grid
	double &Elf,           // Value of the ELF
	double &Eli,           // Value of the ELI
	double &Lap            // Value for the Laplacian
) const
{
	const int _nmo = get_nmo(false);
	vec phi(7 * _nmo, 0.0);
	double *phi_temp;
	double chi[7]{ 0, 0, 0, 0, 0, 0, 0 };
	double d[3]{ 0, 0, 0 };
	int iat = 0, k, j;
	int l[3]{ 0, 0, 0 };
	double ex = 0;
	double xl[3][3]{ {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };

	const MO *MOs_data = MOs.data();
	//Primitive-major coefficients: every MO for one primitive is contiguous here, where
	//MOs keeps them nex apart. Identical arithmetic in identical order - one sequential
	//read per primitive instead of nmo scattered ones, and no per-point heap allocation.
	const double *const coefs = get_coef_primitive_major();

	for (j = 0; j < nex; j++)
	{
		iat = centers[j] - 1;

		constants::type2vector(types[j], l);
		d[0] = PosGrid[0] - atoms[iat].get_coordinate(0);
		d[1] = PosGrid[1] - atoms[iat].get_coordinate(1);
		d[2] = PosGrid[2] - atoms[iat].get_coordinate(2);
		double temp = -exponents[j] * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
		if (temp < constants::exp_cutoff)
			continue;
		ex = exp(temp);
		for (int k = 0; k < 3; k++)
		{
			if (l[k] == 0)
			{
				xl[k][0] = 1.0;
				xl[k][1] = 0.0;
				xl[k][2] = 0.0;
			}
			else if (l[k] == 1)
			{
				xl[k][0] = d[k];
				xl[k][1] = 1.0;
				xl[k][2] = 0.0;
			}
			else if (l[k] == 2)
			{
				xl[k][0] = d[k] * d[k];
				xl[k][1] = 2 * d[k];
				xl[k][2] = 2;
			}
			else if (l[k] == 3)
			{
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d[k];
				xl[k][1] = 3 * d2;
				xl[k][2] = 6 * d[k];
			}
			else if (l[k] == 4)
			{
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d2;
				xl[k][1] = 4 * d2 * d[k];
				xl[k][2] = 12 * d2;
			}
			else if (l[k] == 5)
			{
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d2 * d[k];
				xl[k][1] = 5 * d2 * d2;
				xl[k][2] = 20 * d2 * d[k];
			}
			else if (l[k] == 6)
			{
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d2 * d2;
				xl[k][1] = 6 * d2 * d2 * d[k];
				xl[k][2] = 30 * d2 * d2;
			}
			else
			{
				return;
			}
		}
		double ex2 = 2 * get_exponent(j);
		chi[0] = xl[0][0] * xl[1][0] * xl[2][0] * ex;
		chi[1] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * xl[1][0] * xl[2][0] * ex;
		chi[2] = (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[0][0] * xl[2][0] * ex;
		chi[3] = (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * xl[0][0] * xl[1][0] * ex;
		double temp_ex = pow(ex2, 2);
		chi[4] = (xl[0][2] - ex2 * (2 * l[0] + 1) * xl[0][0] + temp_ex * pow(d[0], l[0] + 2)) * xl[1][0] * xl[2][0] * ex;
		chi[5] = (xl[1][2] - ex2 * (2 * l[1] + 1) * xl[1][0] + temp_ex * pow(d[1], l[1] + 2)) * xl[2][0] * xl[0][0] * ex;
		chi[6] = (xl[2][2] - ex2 * (2 * l[2] + 1) * xl[2][0] + temp_ex * pow(d[2], l[2] + 2)) * xl[0][0] * xl[1][0] * ex;

		const double *c_row = coefs + (size_t)j * nmo;
		double *phi_ptr = phi.data();
		for (int mo = 0; mo < nmo; ++mo, phi_ptr += 7)
		{
			const double c = c_row[mo];
			for (k = 0; k < 7; k++)
				phi_ptr[k] += c * chi[k];
		}
	}

	double Grad[3]{ 0, 0, 0 };
	double Hess[3]{ 0, 0, 0 };
	double tau = 0;
	double Rho = 0;

	Elf = 0;

	for (int mo = 0; mo < _nmo; mo++)
	{
		const double occ = get_MO_occ(mo);
		const double docc = 2 * occ;
		if (occ != 0)
		{
			phi_temp = &phi[mo * 7];
			Rho += occ * pow(*phi_temp, 2);
			Grad[0] += docc * *phi_temp * phi_temp[1];
			Grad[1] += docc * *phi_temp * phi_temp[2];
			Grad[2] += docc * *phi_temp * phi_temp[3];
			Hess[0] += docc * (*phi_temp * phi_temp[4] + pow(phi_temp[1], 2));
			Hess[1] += docc * (*phi_temp * phi_temp[5] + pow(phi_temp[2], 2));
			Hess[2] += docc * (*phi_temp * phi_temp[6] + pow(phi_temp[3], 2));
			tau += occ * (pow(phi_temp[1], 2) + pow(phi_temp[2], 2) + pow(phi_temp[3], 2));
		}
	}
	Elf = 1 / (1 + pow(constants::ctelf * pow(Rho, constants::c_m53) * (tau * 0.5 - 0.125 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2)) / Rho), 2));
	Eli = 0.5 * Rho * pow(48 / (Rho * tau - 0.25 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2))), constants::c_38);
	Lap = Hess[0] + Hess[1] + Hess[2];
};

const void WFN::computeLapELI(
	const d3 &PosGrid, // [3] vector with current position on te grid
	double &Eli,           // Value of the ELI
	double &Lap            // Value for the Laplacian
) const
{
	const int _nmo = get_nmo(false);
	vec phi(7 * _nmo, 0.0);
	double *phi_temp;
	double chi[7]{ 0, 0, 0, 0, 0, 0, 0 };
	double d[3]{ 0, 0, 0 };
	int k, j;
	int l[3]{ 0, 0, 0 };
	double ex = 0;
	double xl[3][3]{ {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };

	const MO *MOs_data = MOs.data();
	double *phi_data = phi.data();
	//Primitive-major coefficients: every MO for one primitive is contiguous here, where
	//MOs keeps them nex apart. Identical arithmetic in identical order - one sequential
	//read per primitive instead of nmo scattered ones, and no per-point heap allocation.
	const double *const coefs = get_coef_primitive_major();

	for (j = 0; j < nex; j++)
	{
		k = centers[j] - 1;

		constants::type2vector(get_type(j), l);
		d[0] = PosGrid[0] - atoms[k].get_coordinate(0);
		d[1] = PosGrid[1] - atoms[k].get_coordinate(1);
		d[2] = PosGrid[2] - atoms[k].get_coordinate(2);
		const double temp = -get_exponent(j) * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
		if (temp < constants::exp_cutoff)
			continue;
		ex = exp(temp);
		for (int k = 0; k < 3; k++)
		{
			switch (l[k])
			{
			case 0:
			{
				xl[k][0] = 1.0;
				xl[k][1] = 0.0;
				xl[k][2] = 0.0;
				break;
			}
			case 1:
			{
				xl[k][0] = d[k];
				xl[k][1] = 1.0;
				xl[k][2] = 0.0;
				break;
			}
			case 2:
			{
				xl[k][0] = d[k] * d[k];
				xl[k][1] = 2 * d[k];
				xl[k][2] = 2;
				break;
			}
			case 3:
			{
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d[k];
				xl[k][1] = 3 * d2;
				xl[k][2] = 6 * d[k];
				break;
			}
			case 4:
			{
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d2;
				xl[k][1] = 4 * d2 * d[k];
				xl[k][2] = 12 * d2;
				break;
			}
			case 5:
			{
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d2 * d[k];
				xl[k][1] = 5 * d2 * d2;
				xl[k][2] = 20 * d2 * d[k];
				break;
			}
			case 6:
			{
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d2 * d2;
				xl[k][1] = 6 * d2 * d2 * d[k];
				xl[k][2] = 30 * d2 * d2;
				break;
			}
			default:
			{
				break;
			}
			}
		}
		const double ex2 = 2 * get_exponent(j);
		chi[0] = xl[0][0] * xl[1][0] * xl[2][0] * ex;
		chi[1] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * xl[1][0] * xl[2][0] * ex;
		chi[2] = (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[0][0] * xl[2][0] * ex;
		chi[3] = (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * xl[0][0] * xl[1][0] * ex;
		const double temp_ex = pow(ex2, 2);
		chi[4] = (xl[0][2] - ex2 * (2 * l[0] + 1) * xl[0][0] + temp_ex * pow(d[0], l[0] + 2)) * xl[1][0] * xl[2][0] * ex;
		chi[5] = (xl[1][2] - ex2 * (2 * l[1] + 1) * xl[1][0] + temp_ex * pow(d[1], l[1] + 2)) * xl[2][0] * xl[0][0] * ex;
		chi[6] = (xl[2][2] - ex2 * (2 * l[2] + 1) * xl[2][0] + temp_ex * pow(d[2], l[2] + 2)) * xl[0][0] * xl[1][0] * ex;

		const double *c_row = coefs + (size_t)j * nmo;
		double *phi_ptr = phi.data();
		for (int mo = 0; mo < nmo; ++mo, phi_ptr += 7)
		{
			const double c = c_row[mo];
			phi_ptr[0] += c * chi[0];
			phi_ptr[1] += c * chi[1];
			phi_ptr[2] += c * chi[2];
			phi_ptr[3] += c * chi[3];
			phi_ptr[4] += c * chi[4];
			phi_ptr[5] += c * chi[5];
			phi_ptr[6] += c * chi[6];
		}
	}

	double Grad[3]{ 0, 0, 0 };
	double Hess[3]{ 0, 0, 0 };
	double tau = 0;
	double Rho = 0;

	for (int mo = 0; mo < _nmo; mo++)
	{
		const double occ = get_MO_occ(mo);
		const double docc = 2 * occ;
		if (occ != 0)
		{
			phi_temp = &phi[mo * 7];
			Rho += occ * pow(*phi_temp, 2);
			Grad[0] += docc * *phi_temp * phi_temp[1];
			Grad[1] += docc * *phi_temp * phi_temp[2];
			Grad[2] += docc * *phi_temp * phi_temp[3];
			Hess[0] += docc * (*phi_temp * phi_temp[4] + pow(phi_temp[1], 2));
			Hess[1] += docc * (*phi_temp * phi_temp[5] + pow(phi_temp[2], 2));
			Hess[2] += docc * (*phi_temp * phi_temp[6] + pow(phi_temp[3], 2));
			tau += occ * (pow(phi_temp[1], 2) + pow(phi_temp[2], 2) + pow(phi_temp[3], 2));
		}
	}
	Eli = 0.5 * Rho * pow(48 / (Rho * tau - 0.25 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2))), constants::c_38);
	Lap = Hess[0] + Hess[1] + Hess[2];
};

//WTF?! Not having this leads to undefined behaviour on Rocky apparently but not Ubuntu/Windows?!
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC push_options
#pragma GCC optimize ("O3")
#endif
const double WFN::computeLap(
	const d3 &PosGrid // [3] vector with current position on te grid
) const
{
	const int _nmo = get_nmo(false);
	vec phi(7 * _nmo, 0.0);
	double chi[7]{ 0, 0, 0, 0, 0, 0, 0 };
	double d[3]{ 0, 0, 0 };
	int j, k;
	int l[3]{ 0, 0, 0 };
	double ex = 0;
	double xl[9]{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	const MO *MOs_data = MOs.data();
	//Primitive-major coefficients: every MO for one primitive is contiguous here, where
	//MOs keeps them nex apart. Identical arithmetic in identical order - one sequential
	//read per primitive instead of nmo scattered ones, and no per-point heap allocation.
	const double *const coefs = get_coef_primitive_major();

	for (j = 0; j < nex; j++)
	{
		k = centers[j] - 1;

		constants::type2vector(get_type(j), l);
		d[0] = PosGrid[0] - atoms[k].get_coordinate(0);
		d[1] = PosGrid[1] - atoms[k].get_coordinate(1);
		d[2] = PosGrid[2] - atoms[k].get_coordinate(2);
		const double temp = -exponents[j] * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
		if (temp < constants::exp_cutoff)
			continue;
		ex = exp(temp);
		for (k = 0; k < 3; k++)
		{
			switch (l[k])
			{
			case 0:
			{
				xl[3 * k + 0] = 1.0;
				xl[3 * k + 1] = 0.0;
				xl[3 * k + 2] = 0.0;
				break;
			}
			case 1:
			{
				xl[3 * k + 0] = d[k];
				xl[3 * k + 1] = 1.0;
				xl[3 * k + 2] = 0.0;
				break;
			}
			case 2:
			{
				xl[3 * k + 0] = d[k] * d[k];
				xl[3 * k + 1] = 2 * d[k];
				xl[3 * k + 2] = 2;
				break;
			}
			case 3:
			{
				const double d2 = d[k] * d[k];
				xl[3 * k + 0] = d2 * d[k];
				xl[3 * k + 1] = 3 * d2;
				xl[3 * k + 2] = 6 * d[k];
				break;
			}
			case 4:
			{
				const double d2 = d[k] * d[k];
				xl[3 * k + 0] = d2 * d2;
				xl[3 * k + 1] = 4 * d2 * d[k];
				xl[3 * k + 2] = 12 * d2;
				break;
			}
			case 5:
			{
				const double d2 = d[k] * d[k];
				xl[3 * k + 0] = d2 * d2 * d[k];
				xl[3 * k + 1] = 5 * d2 * d2;
				xl[3 * k + 2] = 20 * d2 * d[k];
				break;
			}
			case 6:
			{
				const double d2 = d[k] * d[k];
				xl[3 * k + 0] = d2 * d2 * d2;
				xl[3 * k + 1] = 6 * d2 * d2 * d[k];
				xl[3 * k + 2] = 30 * d2 * d2;
				break;
			}
			default:
			{
				return -100;
				break;
			}
			}
		}
		const double ex2 = 2 * exponents[j];
		chi[0] = xl[0] * xl[3] * xl[6] * ex;
		chi[1] = (xl[1] - ex2 * pow(d[0], l[0] + 1)) * xl[3] * xl[6] * ex;
		chi[2] = (xl[4] - ex2 * pow(d[1], l[1] + 1)) * xl[0] * xl[6] * ex;
		chi[3] = (xl[7] - ex2 * pow(d[2], l[2] + 1)) * xl[0] * xl[3] * ex;
		const double ex4 = ex2 * ex2;
		chi[4] = (xl[2] - ex2 * (2 * l[0] + 1) * xl[0] + ex4 * pow(d[0], l[0] + 2)) * xl[3] * xl[6] * ex;
		chi[5] = (xl[5] - ex2 * (2 * l[1] + 1) * xl[3] + ex4 * pow(d[1], l[1] + 2)) * xl[6] * xl[0] * ex;
		chi[6] = (xl[8] - ex2 * (2 * l[2] + 1) * xl[6] + ex4 * pow(d[2], l[2] + 2)) * xl[0] * xl[3] * ex;
		const double *c_row = coefs + (size_t)j * nmo;
		double *phi_ptr = phi.data();
		for (int mo = 0; mo < nmo; ++mo, phi_ptr += 7)
		{
			const double c = c_row[mo];
			phi_ptr[0] += c * chi[0];
			phi_ptr[1] += c * chi[1];
			phi_ptr[2] += c * chi[2];
			phi_ptr[3] += c * chi[3];
			phi_ptr[4] += c * chi[4];
			phi_ptr[5] += c * chi[5];
			phi_ptr[6] += c * chi[6];
		}

	}

	double Hess[3]{ 0, 0, 0 };

	for (int mo = 0; mo < _nmo; mo++)
	{
		const double occ = 2 * get_MO_occ(mo);
		if (occ != 0)
		{
			const double *phi_temp = &phi[mo * 7];
			Hess[0] += occ * (*phi_temp * phi_temp[4] + pow(phi_temp[1], 2));
			Hess[1] += occ * (*phi_temp * phi_temp[5] + pow(phi_temp[2], 2));
			Hess[2] += occ * (*phi_temp * phi_temp[6] + pow(phi_temp[3], 2));
		}
	}
	return Hess[0] + Hess[1] + Hess[2];
};
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC pop_options
#endif

const double WFN::computeMO(
	const d3 &PosGrid, // [3] array with current position on the grid
	const int &mo) const
{
	double result = 0.0;
	int l[3]{ 0, 0, 0 };
	double ex = 0, *d_;
	// x, y, z, r^2 and the powers 2..5 laid out as in compute_dens_cartesian; per-thread scratch
	// because this is called per grid point per orbital from parallel loops
	thread_local vec d;
	if (d.size() < (size_t)16 * ncen) d.resize((size_t)16 * ncen);
	for (int iat = 0; iat < ncen; iat++)
	{
		d_ = d.data() + (size_t)16 * iat;
		d_[0] = PosGrid[0] - atoms[iat].get_coordinate(0);
		d_[1] = PosGrid[1] - atoms[iat].get_coordinate(1);
		d_[2] = PosGrid[2] - atoms[iat].get_coordinate(2);
		d_[4] = d_[0] * d_[0];
		d_[5] = d_[1] * d_[1];
		d_[6] = d_[2] * d_[2];
		d_[3] = d_[4] + d_[5] + d_[6];
		for (int k = 7; k < 16; k++)
			d_[k] = d_[k - 3] * d_[(k - 7) % 3];
	}
	const double *c = MOs[mo].get_coefficient_ptr();
	for (int j = 0; j < nex; j++)
	{
		d_ = d.data() + (size_t)16 * (get_center(j) - 1);
		ex = -get_exponent(j) * d_[3];
		if (ex < constants::exp_cutoff)
			continue;
		ex = exp(ex);
		constants::type2vector(get_type(j), l);
		// power p of coordinate k sits at k for p = 1 and at 3p - 2 + k above
		for (int k = 0; k < 3; k++)
			if (l[k] == 6) ex *= d_[13 + k] * d_[k];
			else if (l[k]) ex *= d_[(l[k] == 1 ? 0 : 3 * l[k] - 2) + k];
		result += c[j] * ex;
	}
	return result;
}

// Boys function F_m(T) by a 6-term Taylor expansion around a tabulated grid (step 0.1 up to
// T = 30, error < 1E-10), asymptotic beyond (1E-14); expn = exp(-T), which the caller has anyway
static double boys(const int m, const double T, const double expn)
{
	constexpr int mmax = 24 + 6, nT = 301;
	constexpr double step = 0.1;
	static const vec table = []()
	{
		vec F((size_t)nT * (mmax + 1));
		for (int i = 0; i < nT; i++)
		{
			const double T0 = i * step, e = exp(-T0);
			// F_mmax by its series e^-T sum_k (2T)^k / ((2m+1)...(2m+2k+1)), then downward
			double term = 1.0 / (2 * mmax + 1), sum = term;
			for (int k = 1; term > 1E-17 * sum; k++)
				term *= 2 * T0 / (2 * mmax + 2 * k + 1), sum += term;
			double *row = F.data() + (size_t)i * (mmax + 1);
			row[mmax] = e * sum;
			for (int n = mmax - 1; n >= 0; n--)
				row[n] = (2 * T0 * row[n + 1] + e) / (2 * n + 1);
		}
		return F;
	}();
	if (T >= (nT - 1) * step)
	{
		double f = 0.5 * sqrt(constants::PI / T); // F_0, then upward, stable at this T
		for (int n = 1; n <= m; n++)
			f = ((2 * n - 1) * f - expn) / (2 * T);
		return f;
	}
	const int i = (int)(T / step + 0.5);
	const double dT = i * step - T, *row = table.data() + (size_t)i * (mmax + 1) + m;
	double f = 0, pw = 1;
	for (int k = 0; k <= 5; k++, pw *= dT / k)
		f += row[k] * pw;
	return f;
}
const double WFN::fj(int &j, int &l, int &m, double &aa, double &bb) const
{
	double temp = 0.0;
	double temp2 = 0.0;
	int a = 0, b = 0;
	for (int i = std::max(0, j - m); i <= std::min(j, l); i++)
	{
		// pre = factorial[l] / factorial[l - i] / factorial[i] * factorial[m] / factorial[m - j + i] / factorial[j - i];
		temp2 = static_cast<double>(pre[j][l][m][i]);
		a = l - i;
		b = m + i - j;
		if (a != 0)
			temp2 *= pow(aa, a);
		if (b != 0)
			temp2 *= pow(bb, b);
		temp += temp2;
	}
	return temp;
};

const double WFN::Afac(int &l, int &r, int &i, double &PC, double &gamma, double &fjtmp) const
{
	double temp = fjtmp * pow(0.25 / gamma, r + i) / Afac_pre[l][r][i];
	const int num = l - 2 * r - 2 * i;
	if (num != 0)
		temp *= pow(PC, num);
	if (i % 2 == 1)
		return -temp;
	else
		return temp;
}

// number of (l, r, s) terms of one axis with l_i + l_j = L
static int esp_axis_terms(const int L)
{
	int n = 0;
	for (int l = 0; l <= L; l++)
		for (int r = 0; r <= l / 2; r++)
			n += (l - 2 * r) / 2 + 1;
	return n;
}

WFN::ESP_pairs WFN::build_ESP_pairs() const
{
	ESP_pairs t;
	const int MO = get_nmo(true), stride = get_nmo(false);
	const int nprim = get_nex();
	const double *coef = get_coef_primitive_major(); // [prim][mo]
	int l_i[3], l_j[3];
	double Pi[3], Pj[3];
	t.off.push_back(0);
	for (int iprim = 0; iprim < nprim; iprim++)
	{
		const int iat = get_center(iprim) - 1;
		constants::type2vector(get_type(iprim), l_i);
		const double iex = get_exponent(iprim);
		for (int jprim = iprim; jprim < nprim; jprim++)
		{
			const int jat = get_center(jprim) - 1;
			const double jex = get_exponent(jprim);
			const double ex_sum = iex + jex;
			double sqd = 0;
			for (int k = 0; k < 3; k++)
				sqd += pow(atoms[iat].get_coordinate(k) - atoms[jat].get_coordinate(k), 2);
			const double prefac = constants::TWO_PI / ex_sum * exp(-iex * jex * sqd / ex_sum);
			if (prefac < 1E-10)
				continue;
			double dij = 0;
			for (int mo = 0; mo < MO; mo++)
				dij += get_MO_occ(mo) * coef[iprim * stride + mo] * coef[jprim * stride + mo];
			// screening: the dropped pairs move the far-field ESP of sucrose by 3E-6 at most (18 Sep 2026)
			if (abs(prefac * dij) < 1E-8)
				continue;
			constants::type2vector(get_type(jprim), l_j);
			t.ex_sum.push_back(ex_sum);
			t.weight.push_back(prefac * dij * (iprim != jprim ? 2.0 : 1.0));
			d3 P;
			for (int k = 0; k < 3; k++)
			{
				P[k] = (atoms[iat].get_coordinate(k) * iex + atoms[jat].get_coordinate(k) * jex) / ex_sum;
				Pi[k] = P[k] - atoms[iat].get_coordinate(k);
				Pj[k] = P[k] - atoms[jat].get_coordinate(k);
			}
			t.P.push_back(P);
			t.L.push_back({ l_i[0] + l_j[0], l_i[1] + l_j[1], l_i[2] + l_j[2] });
			// Afac without its PC power: sign * fj * (1/4g)^(r+s) / (r! s! (l-2r-2s)!)
			for (int k = 0; k < 3; k++)
				for (int l = 0; l <= l_i[k] + l_j[k]; l++)
				{
					const double fjtmp = (l % 2 ? -1.0 : 1.0) * fj(l, l_i[k], l_j[k], Pi[k], Pj[k]);
					for (int r = 0; r <= l / 2; r++)
						for (int s = 0; s <= (l - 2 * r) / 2; s++)
						{
							t.coef.push_back((s % 2 ? -fjtmp : fjtmp) * pow(0.25 / ex_sum, r + s) / Afac_pre[l][r][s]);
							t.pc_pow.push_back((unsigned char)(l - 2 * r - 2 * s));
							t.fn_idx.push_back((unsigned char)(l - 2 * r - s));
						}
				}
			t.off.push_back((int)t.coef.size());
		}
	}
	return t;
}

const double WFN::computeESP(const d3 &PosGrid, const ESP_pairs &t) const
{
	double ESP = 0;
	for (int iat = 0; iat < get_ncen(); iat++)
	{
		double r2 = 0;
		for (int k = 0; k < 3; k++)
			r2 += pow(PosGrid[k] - atoms[iat].get_coordinate(k), 2);
		ESP += (get_atom_charge(iat) - atoms[iat].get_ECP_electrons()) / sqrt(r2); // ECP/xTB/pTB: only the valence electrons are in the MOs, so the core must not count as nuclear charge
	}

	double Fn[25], pcp[3][9], Al[506], Am[506], An[506]; // l_i, l_j <= 4 per axis (pre tables)
	int mapl[506], mapm[506], mapn[506];
	const int npairs = (int)t.weight.size();
	for (int p = 0; p < npairs; p++)
	{
		const double ex_sum = t.ex_sum[p];
		const std::array<int, 3> &L = t.L[p];
		double sqpc = 0;
		for (int k = 0; k < 3; k++)
		{
			const double PC = t.P[p][k] - PosGrid[k];
			sqpc += PC * PC;
			pcp[k][0] = 1.0;
			for (int n = 1; n <= L[k]; n++)
				pcp[k][n] = pcp[k][n - 1] * PC;
		}
		double expc = exp(-ex_sum * sqpc);
		int MaxFn = L[0] + L[1] + L[2];
		Fn[MaxFn] = boys(MaxFn, ex_sum * sqpc, expc);
		const double twoexpc = 2 * ex_sum * sqpc;
		for (int nu = MaxFn - 1; nu >= 0; nu--)
			Fn[nu] = (expc + twoexpc * Fn[nu + 1]) / (2 * (nu + 1) - 1);

		int c = t.off[p];
		const int nl = esp_axis_terms(L[0]), nm = esp_axis_terms(L[1]), nn = esp_axis_terms(L[2]);
		for (int l = 0; l < nl; l++, c++)
			Al[l] = t.coef[c] * pcp[0][t.pc_pow[c]], mapl[l] = t.fn_idx[c];
		for (int m = 0; m < nm; m++, c++)
			Am[m] = t.coef[c] * pcp[1][t.pc_pow[c]], mapm[m] = t.fn_idx[c];
		for (int n = 0; n < nn; n++, c++)
			An[n] = t.coef[c] * pcp[2][t.pc_pow[c]], mapn[n] = t.fn_idx[c];

		double term = 0.0;
		for (int l = 0; l < nl; l++)
			for (int m = 0; m < nm; m++)
			{
				const double lm = Al[l] * Am[m];
				for (int n = 0; n < nn; n++)
					term += lm * An[n] * Fn[mapl[l] + mapm[m] + mapn[n]];
			}
		ESP -= t.weight[p] * term;
	}
	return ESP;
};
