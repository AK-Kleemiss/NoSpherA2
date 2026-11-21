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
#include "eli_family.h"
#include "integrator.h"
#include "cell.h"
#include "wfn_kernels.h"
#ifdef NOSPHERA2_USE_GPU
#include "aux_density_gpu.h"
#include "esp_gpu.h"
#include "basin_field_gpu.h"
#endif

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

bool WFN::eli_orbital_pass(const d3 &PosGrid, vec &phi) const
{
	const int _nmo = get_nmo(false);
	thread_local vec d;
	phi.assign(10 * _nmo, 0.0);
	if (d.size() < 4 * (size_t)ncen) d.resize(4 * (size_t)ncen);
	for (int j = 0; j < ncen; j++)
	{
		double *d_ = d.data() + 4 * j;
		for (int k = 0; k < 3; k++) d_[k] = PosGrid[k] - atoms[j].get_coordinate(k);
		d_[3] = d_[0] * d_[0] + d_[1] * d_[1] + d_[2] * d_[2];
	}
	if (!get_coef_primitive_major())
		return false;
	thread_local vec exps;
	if (exps.size() < group_exponent.size()) exps.resize(group_exponent.size());
	exp_table([r2 = d.data()](const int c) { return r2[4 * (size_t)c + 3]; }, exps.data());
	const int *group = prim_exp_group.data();
	thread_local ao_scratch<10> ao;
	ao.begin(coef_ao_major.size() / _nmo);
	for (int j = 0; j < nex; j++)
	{
		const double ex = exps[group[j]];
		if (ex == 0.0 || prim_ao[j] < 0)
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
		ao.add(prim_ao[j], prim_ao_scale[j], chi);
	}
	//Component-major, as in computeGrad: ten unit-stride accumulations per function
	ao.to_mo(coef_ao_major.data(), _nmo, phi.data(), [this](const int mo) { return MOs[mo].get_occ(); });
	return true;
}

void WFN::computeELIGrad(
	const d3 &PosGrid,
	double& out_Eli,
	d3& out_grad,
	double *out_rho
) const
{
	//ELI-D Y = rho/2 (48/g)^(3/8) with g = rho tau - |grad rho|^2 / 4, so
	//grad Y = (48/g)^(3/8) / 2 (grad rho - 3/8 rho grad g / g) with
	//grad g = tau grad rho + rho grad tau - H grad rho / 2, H the density Hessian
	const int _nmo = get_nmo(false);
	thread_local vec phi;
	if (!eli_orbital_pass(PosGrid, phi))
	{
		out_Eli = 0;
		out_grad = { 0, 0, 0 };
		if (out_rho) *out_rho = 0.0;
		return;
	}
	static constexpr int hidx[3][3] = { {4, 7, 8}, {7, 5, 9}, {8, 9, 6} };
	double rho = 0, tau = 0, G[3]{ 0, 0, 0 }, T[3]{ 0, 0, 0 }, H[3][3]{ {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };
	for (int mo = 0; mo < _nmo; mo++)
	{
		const double occ = get_MO_occ(mo);
		if (occ == 0) continue;
		double p[10];
		for (int k = 0; k < 10; k++) p[k] = phi[(size_t)k * _nmo + mo];
		const double docc = 2 * occ;
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
	if (out_rho) *out_rho = rho;
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

void WFN::computeELISpinGrad(
	const d3 &PosGrid,
	const int field,
	const double triplet_factor,
	double& out_Eli,
	d3& out_grad,
	double *aux
) const
{
	//Per spin channel s the same ingredients as computeELIGrad: rho_s, T_s = sum n |grad phi|^2,
	//G_s = grad rho_s, dT_s = grad T_s, H_s the Hessian of rho_s.  Same spin (eli_family.h):
	//g = rho_s T_s - |G_s|^2/4, grad g = G_s T_s + rho_s dT_s - H_s G_s / 2; triplet:
	//g = g_a + g_b + (rho_a T_b + rho_b T_a)/2 - G_a.G_b/4, differentiated term by term.
	//Y = r (12/g)^(3/8), r = rho_s or triplet_factor * rho, grad Y = (12/g)^(3/8) (grad r - 3/8 r grad g / g)
	const int _nmo = get_nmo(false);
	thread_local vec phi;
	out_Eli = 0;
	out_grad = { 0, 0, 0 };
	if (aux) for (int k = 0; k < 4; k++) aux[k] = 0.0;
	if (!eli_orbital_pass(PosGrid, phi)) return;
	static constexpr int hidx[3][3] = { {4, 7, 8}, {7, 5, 9}, {8, 9, 6} };
	//spin_split scans serially: get_MO_op_count opens an OpenMP region and this runs per climb point
	const eli_family::SpinSplit how = eli_family::spin_split(*this);
	const bool halves = how == eli_family::SpinSplit::halves;
	double rho[2]{ 0, 0 }, tau[2]{ 0, 0 }, G[2][3]{}, T[2][3]{}, H[2][3][3]{};
	for (int mo = 0; mo < _nmo; mo++)
	{
		const double occ = get_MO_occ(mo);
		if (occ == 0) continue;
		double p[10];
		for (int k = 0; k < 10; k++) p[k] = phi[(size_t)k * _nmo + mo];
		//halves: half of every occupation in each channel, accumulated once into alpha and copied
		double nc[2];
		eli_family::mo_spin_occupations(occ, get_MO_op(mo), how, nc);
		for (int s = 0; s < (halves ? 1 : 2); s++)
		{
			const double n = nc[s], dn = 2 * n;
			if (n == 0.0) continue;
			rho[s] += n * p[0] * p[0];
			for (int i = 0; i < 3; i++)
			{
				tau[s] += n * p[1 + i] * p[1 + i];
				G[s][i] += dn * p[0] * p[1 + i];
				for (int k = 0; k < 3; k++)
				{
					H[s][i][k] += dn * (p[0] * p[hidx[i][k]] + p[1 + i] * p[1 + k]);
					T[s][k] += dn * p[1 + i] * p[hidx[i][k]];
				}
			}
		}
	}
	if (halves)
	{
		rho[1] = rho[0]; tau[1] = tau[0];
		for (int i = 0; i < 3; i++) { G[1][i] = G[0][i]; T[1][i] = T[0][i]; for (int k = 0; k < 3; k++) H[1][i][k] = H[0][i][k]; }
	}
	//g_s and grad g_s of one channel
	auto same = [&](const int s, double &g, double *dg) {
		g = rho[s] * tau[s] - 0.25 * (G[s][0] * G[s][0] + G[s][1] * G[s][1] + G[s][2] * G[s][2]);
		for (int k = 0; k < 3; k++)
			dg[k] = G[s][k] * tau[s] + rho[s] * T[s][k] - 0.5 * (G[s][0] * H[s][0][k] + G[s][1] * H[s][1][k] + G[s][2] * H[s][2][k]);
	};
	double g, dg[3], r, dr[3];
	if (field == 2)
	{
		double ga, gb, dga[3], dgb[3];
		same(0, ga, dga);
		same(1, gb, dgb);
		g = ga + gb + 0.5 * (rho[0] * tau[1] + rho[1] * tau[0]) - 0.25 * (G[0][0] * G[1][0] + G[0][1] * G[1][1] + G[0][2] * G[1][2]);
		for (int k = 0; k < 3; k++)
		{
			dg[k] = dga[k] + dgb[k] + 0.5 * (G[0][k] * tau[1] + rho[0] * T[1][k] + G[1][k] * tau[0] + rho[1] * T[0][k])
				- 0.25 * (H[0][k][0] * G[1][0] + H[0][k][1] * G[1][1] + H[0][k][2] * G[1][2] + H[1][k][0] * G[0][0] + H[1][k][1] * G[0][1] + H[1][k][2] * G[0][2]);
			dr[k] = triplet_factor * (G[0][k] + G[1][k]);
		}
		r = triplet_factor * (rho[0] + rho[1]);
	}
	else
	{
		same(field, g, dg);
		r = rho[field];
		for (int k = 0; k < 3; k++) dr[k] = G[field][k];
	}
	if (aux)
	{
		aux[0] = rho[0] + rho[1];
		aux[1] = rho[0];
		aux[2] = rho[1];
		//rho_s Y_q = g_s / (12 rho_s^(5/3)), so the basin average of ELI-q is this integral over N_s
		aux[3] = field != 2 && r > 0.0 && g > 0.0 ? g / (12.0 * std::pow(r, 5.0 / 3.0)) : 0.0;
	}
	//an empty channel, or a single orbital carrying it: no pair, no finite ELI-D
	if (!(g > 0) || !(r > 0)) return;
	const double f = pow(12 / g, constants::c_38);
	out_Eli = r * f;
	for (int k = 0; k < 3; k++) out_grad[k] = f * (dr[k] - 0.375 * r * dg[k] / g);
}

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
	double chi[4]{ 0, 0, 0, 0 };
	int k, j;
	double ex = 0;

	for (j = 0; j < ncen; j++)
		cart_powers(d.data() + 16 * j, PosGrid, atoms[j]);

	const int *centers_data = centers.data();
	const int *types_data = types.data();
	const double *exponents_data = exponents.data();

	if (!get_coef_primitive_major())
	{
		gradient = { 0, 0, 0 };
		if (rho) *rho = 0.0;
		return;
	}
	thread_local vec exps;
	if (exps.size() < group_exponent.size()) exps.resize(group_exponent.size());
	exp_table([r2 = d.data()](const int c) { return r2[16 * (size_t)c + 3]; }, exps.data());
	const int *group = prim_exp_group.data();
	thread_local ao_scratch<4> ao;
	ao.begin(coef_ao_major.size() / _nmo);

	for (j = 0; j < nex; j++)
	{
		ex = exps[group[j]];
		if (ex == 0.0 || prim_ao[j] < 0)
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

		ao.add(prim_ao[j], prim_ao_scale[j], chi);
	}
	//Component-major: phi[k * nmo + mo], unit stride against one broadcast scalar per function
	ao.to_mo(coef_ao_major.data(), _nmo, phi.data(), [this](const int mo) { return MOs[mo].get_occ(); });

	double Grad[3]{ 0, 0, 0 }, Rho = 0.0;

	for (int mo = 0; mo < _nmo; mo++)
	{
		const double occ = get_MO_occ(mo);
		const double docc = 2 * occ;
		if (occ != 0)
		{
			const double v = phi[mo], gx = phi[(size_t)_nmo + mo],
				gy = phi[(size_t)2 * _nmo + mo], gz = phi[(size_t)3 * _nmo + mo];
			Grad[0] += docc * v * gx;
			Grad[1] += docc * v * gy;
			Grad[2] += docc * v * gz;
			if (rho) Rho += occ * v * v;
		}
	}
	gradient[0] = Grad[0];
	gradient[1] = Grad[1];
	gradient[2] = Grad[2];
	if (rho) *rho = Rho;

};

bool WFN::field_grad_gpu(const bool eli, const int np, const double *pts, double *val, double *grad, double *rho) const
{
	void* ctx = field_gpu_open(eli, np);
	const bool ok = field_gpu_run(ctx, np, pts, val, grad, rho);
	field_gpu_close(ctx);
	return ok;
}

bool WFN::field_gpu_run(void* ctx, const int np, const double *pts, double *val, double *grad, double *rho)
{
#ifdef NOSPHERA2_USE_GPU
	return ctx && basin_field_gpu_run(ctx, np, pts, val, grad, rho);
#else
	(void)ctx; (void)np; (void)pts; (void)val; (void)grad; (void)rho;
	return false;
#endif
}

void WFN::field_gpu_close(void* ctx)
{
#ifdef NOSPHERA2_USE_GPU
	if (ctx) basin_field_gpu_close(ctx);
#else
	(void)ctx;
#endif
}

void* WFN::field_gpu_open(const bool eli, const int max_points) const
{
#ifdef NOSPHERA2_USE_GPU
	const int _nmo = get_nmo(false);
	if (max_points <= 0 || !get_coef_primitive_major())
		return nullptr;
	const int nao = (int)(coef_ao_major.size() / _nmo);
	//The primitives the host evaluators visit (a function, a type up to l = 10), grouped by
	//function in ascending wfn order, so every function sums its primitives as ao.add does
	ivec start(nao + 1, 0), lj(3 * (size_t)nex, -1);
	for (int j = 0; j < nex; j++)
	{
		if (prim_ao[j] < 0) continue;
		constants::type2vector(types[j], &lj[3 * (size_t)j]);
		if (lj[3 * (size_t)j] >= 0) start[prim_ao[j] + 1]++;
	}
	for (int a = 0; a < nao; a++) start[a + 1] += start[a];
	const int nprim = start[nao];
	ivec fill(start.begin(), start.end() - 1), pc(nprim), pl(3 * (size_t)nprim);
	vec pe(nprim), ps(nprim);
	for (int j = 0; j < nex; j++)
	{
		if (prim_ao[j] < 0 || lj[3 * (size_t)j] < 0) continue;
		const int i = fill[prim_ao[j]]++;
		pc[i] = centers[j] - 1;
		std::copy_n(&lj[3 * (size_t)j], 3, &pl[3 * (size_t)i]);
		pe[i] = exponents[j];
		ps[i] = prim_ao_scale[j];
	}
	//only the occupied MOs, in order: the host reduction skips the others
	ivec mos;
	vec occ;
	for (int mo = 0; mo < _nmo; mo++)
		if (get_MO_occ(mo) != 0) { mos.push_back(mo); occ.push_back(get_MO_occ(mo)); }
	const int nocc = (int)mos.size();
	vec coef((size_t)nao * nocc), cxyz(3 * (size_t)ncen);
	for (int a = 0; a < nao; a++)
		for (int i = 0; i < nocc; i++) coef[(size_t)a * nocc + i] = coef_ao_major[(size_t)a * _nmo + mos[i]];
	for (int c = 0; c < ncen; c++)
		for (int k = 0; k < 3; k++) cxyz[3 * (size_t)c + k] = atoms[c].get_coordinate(k);
	if (nocc == 0) return nullptr;
	return basin_field_gpu_open(eli ? 10 : 4, ncen, cxyz.data(), center_min_exponent.data(), constants::exp_cutoff,
		nao, start.data(), pc.data(), pl.data(), pe.data(), ps.data(), nocc, coef.data(), occ.data(), max_points);
#else
	(void)eli; (void)max_points;
	return nullptr;
#endif
}

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
