#pragma once

#include "convenience.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char* file, int line, bool abort = true)
{
	if (code != cudaSuccess)
	{
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}

__constant__ int gpu_nex[1];
__constant__ int gpu_nmo[1];
__constant__ int gpu_ncen[1];
__constant__ int gpu_MaxGrid[1];
__constant__ float gpu_log_incr[1];
__constant__ float gpu_start_radial_dens[1];
__constant__ float gpu_bragg_angstrom[114]{
0.00,
	0.35, 0.35,
	1.45, 1.05,																																					0.85, 0.70, 0.65, 0.60, 0.50, 0.45,
	1.80, 1.50,																																					1.25, 1.10, 1.00, 1.00, 1.00, 1.00,
	2.20, 1.80,																						1.60, 1.40, 1.35, 1.40, 1.40, 1.40, 1.35, 1.35, 1.35, 1.35, 1.30, 1.25, 1.15, 1.15, 1.15, 1.10,
	2.35, 2.00,																						1.80, 1.55, 1.45, 1.45, 1.35, 1.30, 1.35, 1.40, 1.60, 1.55, 1.55, 1.45, 1.45, 1.40, 1.40, 1.40,
	2.60, 2.15, 1.95, 1.85, 1.85, 1.85, 1.85, 1.85, 1.85, 1.80, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.55, 1.45, 1.35, 1.30, 1.30, 1.35, 1.35, 1.35, 1.50, 1.90, 1.75, 1.60, 1.90, 1.50, 1.50,
	2.80, 2.35, 2.15, 2.05, 2.05, 2.05, 2.05, 2.05, 2.05, 2.00, 1.95, 1.95, 1.95, 1.95, 1.95, 1.95 };

__device__ void gpu_ptype(
	int& vector0, int& vector1, int& vector2,
	const int& l)
{
	if (l == 2)
		vector0 = 1, vector1 = 0, vector2 = 0;
	else if (l == 3)
		vector0 = 0, vector1 = 1, vector2 = 0;
	else if (l == 4)
		vector0 = 0, vector1 = 0, vector2 = 1;
}

__device__ float gpu_ptype(
	const float& vector0, const float& vector1, const float& vector2,
	const int& l)
{
	if (l == 2)
		return vector0;
	else if (l == 3)
		return vector1;
	else if (l == 4)
		return vector2;
	else return -1;
}

__device__ void gpu_dtype(
	int& vector0, int& vector1, int& vector2,
	const int& l)
{
	if (l == 5)
		vector0 = 2, vector1 = 0, vector2 = 0;
	else if (l == 6)
		vector0 = 0, vector1 = 2, vector2 = 0;
	else if (l == 7)
		vector0 = 0, vector1 = 0, vector2 = 2;
	else if (l == 8)
		vector0 = 1, vector1 = 1, vector2 = 0;
	else if (l == 9)
		vector0 = 1, vector1 = 0, vector2 = 1;
	else if (l == 10)
		vector0 = 0, vector1 = 1, vector2 = 1;
}

__device__ float gpu_dtype(
	const float& vector0, const float& vector1, const float& vector2,
	const int& l)
{
	if (l == 5)
		return vector0 * vector0;
	else if (l == 6)
		return vector1 * vector1;
	else if (l == 7)
		return vector2 * vector2;
	else if (l == 8)
		return vector0 * vector1;
	else if (l == 9)
		return vector0 * vector2;
	else if (l == 10)
		return vector1 * vector2;
	else return -1;
}

__device__ void gpu_ftype(
	int& vector0, int& vector1, int& vector2,
	const int& l)
{
	if (l == 11)
		vector0 = 3, vector1 = 0, vector2 = 0;
	else if (l == 12)
		vector0 = 0, vector1 = 3, vector2 = 0;
	else if (l == 13)
		vector0 = 0, vector1 = 0, vector2 = 3;
	else if (l == 14)
		vector0 = 2, vector1 = 1, vector2 = 0;
	else if (l == 15)
		vector0 = 2, vector1 = 0, vector2 = 1;
	else if (l == 16)
		vector0 = 0, vector1 = 2, vector2 = 1;
	else if (l == 17)
		vector0 = 1, vector1 = 2, vector2 = 0;
	else if (l == 18)
		vector0 = 1, vector1 = 0, vector2 = 2;
	else if (l == 19)
		vector0 = 0, vector1 = 1, vector2 = 1;
	else if (l == 20)
		vector0 = 1, vector1 = 1, vector2 = 1;
}

__device__ float gpu_ftype(
	const float& vector0, const float& vector1, const float& vector2,
	const int& l)
{
	if (l == 11)
		return vector0 * vector0 * vector0;
	else if (l == 12)
		return vector1 * vector1 * vector1;
	else if (l == 13)
		return vector2 * vector2 * vector2;
	else if (l == 14)
		return vector0 * vector0 * vector1;
	else if (l == 15)
		return vector0 * vector0 * vector2;
	else if (l == 16)
		return vector1 * vector1 * vector2;
	else if (l == 17)
		return vector0 * vector1 * vector1;
	else if (l == 18)
		return vector0 * vector2 * vector2;
	else if (l == 19)
		return vector1 * vector2 * vector2;
	else if (l == 20)
		return vector0 * vector1 * vector2;
	else return -1;
}

__device__ void gpu_gtype(
	int& vector0, int& vector1, int& vector2,
	const int& l)
{
	if (l == 21)
		vector0 = 0, vector1 = 0, vector2 = 4;
	else if (l == 22)
		vector0 = 0, vector1 = 1, vector2 = 3;
	else if (l == 23)
		vector0 = 0, vector1 = 2, vector2 = 2;
	else if (l == 24)
		vector0 = 0, vector1 = 3, vector2 = 1;
	else if (l == 25)
		vector0 = 0, vector1 = 4, vector2 = 0;
	else if (l == 26)
		vector0 = 1, vector1 = 0, vector2 = 3;
	else if (l == 27)
		vector0 = 1, vector1 = 1, vector2 = 2;
	else if (l == 28)
		vector0 = 1, vector1 = 2, vector2 = 1;
	else if (l == 29)
		vector0 = 1, vector1 = 3, vector2 = 0;
	else if (l == 30)
		vector0 = 2, vector1 = 0, vector2 = 2;
	else if (l == 31)
		vector0 = 2, vector1 = 1, vector2 = 1;
	else if (l == 32)
		vector0 = 2, vector1 = 2, vector2 = 0;
	else if (l == 33)
		vector0 = 3, vector1 = 0, vector2 = 1;
	else if (l == 34)
		vector0 = 3, vector1 = 1, vector2 = 0;
	else if (l == 35)
		vector0 = 4, vector1 = 0, vector2 = 0;
}

__device__ float gpu_gtype(
	const float& vector0, const float& vector1, const float& vector2,
	const int& l)
{
	if (l == 21)
		return vector2 * vector2 * vector2 * vector2;
	else if (l == 22)
		return vector1 * vector2 * vector2 * vector2;
	else if (l == 23)
		return vector1 * vector1 * vector2 * vector2;
	else if (l == 24)
		return vector1 * vector1 * vector1 * vector2;
	else if (l == 25)
		return vector1 * vector1 * vector1 * vector1;
	else if (l == 26)
		return vector0 * vector2 * vector2 * vector2;
	else if (l == 27)
		return vector0 * vector1 * vector2 * vector2;
	else if (l == 28)
		return vector0 * vector1 * vector1 * vector2;
	else if (l == 29)
		return vector0 * vector1 * vector1 * vector1;
	else if (l == 30)
		return vector0 * vector0 * vector2 * vector2;
	else if (l == 31)
		return vector0 * vector0 * vector1 * vector2;
	else if (l == 32)
		return vector0 * vector0 * vector1 * vector1;
	else if (l == 33)
		return vector0 * vector0 * vector0 * vector2;
	else if (l == 34)
		return vector0 * vector0 * vector0 * vector1;
	else if (l == 35)
		return vector0 * vector0 * vector0 * vector0;
	else return -1;
}

__device__ void gpu_htype(
	int& vector0, int& vector1, int& vector2,
	const int& l)
{
	if (l == 36)
		vector0 = 0, vector1 = 0, vector2 = 5;
	else if (l == 37)
		vector0 = 0, vector1 = 1, vector2 = 4;
	else if (l == 38)
		vector0 = 0, vector1 = 2, vector2 = 3;
	else if (l == 39)
		vector0 = 0, vector1 = 3, vector2 = 2;
	else if (l == 40)
		vector0 = 0, vector1 = 4, vector2 = 1;
	else if (l == 41)
		vector0 = 0, vector1 = 5, vector2 = 0;
	else if (l == 42)
		vector0 = 1, vector1 = 0, vector2 = 4;
	else if (l == 43)
		vector0 = 1, vector1 = 1, vector2 = 3;
	else if (l == 44)
		vector0 = 1, vector1 = 2, vector2 = 2;
	else if (l == 45)
		vector0 = 1, vector1 = 3, vector2 = 1;
	else if (l == 46)
		vector0 = 1, vector1 = 4, vector2 = 0;
	else if (l == 47)
		vector0 = 2, vector1 = 0, vector2 = 3;
	else if (l == 48)
		vector0 = 2, vector1 = 1, vector2 = 2;
	else if (l == 49)
		vector0 = 2, vector1 = 2, vector2 = 1;
	else if (l == 50)
		vector0 = 2, vector1 = 3, vector2 = 0;
	else if (l == 51)
		vector0 = 3, vector1 = 0, vector2 = 2;
	else if (l == 52)
		vector0 = 3, vector1 = 1, vector2 = 1;
	else if (l == 53)
		vector0 = 3, vector1 = 2, vector2 = 0;
	else if (l == 54)
		vector0 = 4, vector1 = 0, vector2 = 1;
	else if (l == 55)
		vector0 = 4, vector1 = 1, vector2 = 0;
	else if (l == 56)
		vector0 = 5, vector1 = 0, vector2 = 0;
}

__device__ float gpu_htype(
	const float& vector0, const float& vector1, const float& vector2,
	const int& l)
{
	if (l == 36)
		return vector2 * vector2 * vector2 * vector2 * vector2;
	else if (l == 37)
		return vector1 * vector2 * vector2 * vector2 * vector2;
	else if (l == 38)
		return vector1 * vector1 * vector2 * vector2 * vector2;
	else if (l == 39)
		return vector1 * vector1 * vector1 * vector2 * vector2;
	else if (l == 40)
		return vector1 * vector1 * vector1 * vector1 * vector2;
	else if (l == 41)
		return vector1 * vector1 * vector1 * vector1 * vector1;
	else if (l == 42)
		return vector0 * vector2 * vector2 * vector2 * vector2;
	else if (l == 43)
		return vector0 * vector1 * vector2 * vector2 * vector2;
	else if (l == 44)
		return vector0 * vector1 * vector1 * vector2 * vector2;
	else if (l == 45)
		return vector0 * vector1 * vector1 * vector1 * vector2;
	else if (l == 46)
		return vector0 * vector1 * vector1 * vector1 * vector1;
	else if (l == 47)
		return vector0 * vector0 * vector2 * vector2 * vector2;
	else if (l == 48)
		return vector0 * vector0 * vector1 * vector2 * vector2;
	else if (l == 49)
		return vector0 * vector0 * vector1 * vector1 * vector2;
	else if (l == 50)
		return vector0 * vector0 * vector1 * vector1 * vector1;
	else if (l == 51)
		return vector0 * vector0 * vector0 * vector2 * vector2;
	else if (l == 52)
		return vector0 * vector0 * vector0 * vector1 * vector2;
	else if (l == 53)
		return vector0 * vector0 * vector0 * vector1 * vector1;
	else if (l == 54)
		return vector0 * vector0 * vector0 * vector0 * vector2;
	else if (l == 55)
		return vector0 * vector0 * vector0 * vector0 * vector1;
	else if (l == 56)
		return vector0 * vector0 * vector0 * vector0 * vector0;
	else return -1;
}

__device__ void gpu_type2vector(
	const int& index,
	int& vector0, int& vector1, int& vector2)
{
	if (index == 1) {
		vector0 = 0;
		vector1 = 0;
		vector2 = 0;
	}
	else if (index < 5)
		gpu_ptype(vector0, vector1, vector2, index);
	else if (index < 11)
		gpu_dtype(vector0, vector1, vector2, index);
	else if (index < 21)
		gpu_ftype(vector0, vector1, vector2, index);
	else if (index < 36)
		gpu_gtype(vector0, vector1, vector2, index);
	else
		gpu_htype(vector0, vector1, vector2, index);
}

__device__ float gpu_type2mult(
	const int& index,
	float& vector0, float& vector1, float& vector2)
{
	if (index == 1)
		return 1.0;
	else if (index == 2)
		return vector0;
	else if (index == 3)
		return vector1;
	else if (index == 4)
		return vector2;
	else if (index < 11)
		return gpu_dtype(vector0, vector1, vector2, index);
	else if (index < 21)
		return gpu_ftype(vector0, vector1, vector2, index);
	else if (index < 36)
		return gpu_gtype(vector0, vector1, vector2, index);
	else
		return gpu_htype(vector0, vector1, vector2, index);
}

__device__ int gpu_type_to_order(const int index)
{
	if (index < 1 || index > 35) return -1;
	else if (index == 1) return 0;
	else if (index < 5)  return 1;
	else if (index < 11) return 2;
	else if (index < 21) return 3;
	else if (index < 36) return 4;
	else return 5;
}

// JCP 88, 2547 (1988), eq. 20
__device__ float gpu_f3(const float x)
{
	float f = x;
	for (int i = 0; i < 3; i++) {
		f *= (1.5 - 0.5 * f * f);
	}
	return f;
}

__device__ float get_becke_w(const int num_centers,
	const int* proton_charges,
	const float* x_atoms,
	const float* y_atoms,
	const float* z_atoms,
	const int center_index,
	const float x,
	const float y,
	const float z)
{
	float R_a, R_b;
	float u_ab, a_ab, mu_ab, nu_ab;
	float f, chi;
	float dist_a, dist_b, dist_ab;
	float vx, vy, vz;

	float* pa = (float*)malloc(sizeof(float) * num_centers);

	for (int a = 0; a < num_centers; a++)
		pa[a] = 1.0;

	for (int a = 0; a < num_centers; a++)
	{
		vx = x_atoms[a] - x;
		vy = y_atoms[a] - y;
		vz = z_atoms[a] - z;
		dist_a = vx * vx + vy * vy + vz * vz;
		dist_a = sqrt(dist_a);

		R_a = gpu_bragg_angstrom[proton_charges[a]];

		for (int b = 0; b < a; b++) {
			vx = x_atoms[b] - x;
			vy = y_atoms[b] - y;
			vz = z_atoms[b] - z;
			dist_b = vx * vx + vy * vy + vz * vz;
			dist_b = sqrt(dist_b);

			R_b = gpu_bragg_angstrom[proton_charges[b]];

			vx = x_atoms[b] - x_atoms[a];
			vy = y_atoms[b] - y_atoms[a];
			vz = z_atoms[b] - z_atoms[a];
			dist_ab = vx * vx + vy * vy + vz * vz;
			dist_ab = sqrt(dist_ab);

			// JCP 88, 2547 (1988), eq. 11
			mu_ab = (dist_a - dist_b) / dist_ab;

			if (abs(R_a - R_b) > 1.0e-14) {
				chi = R_a / R_b;
				u_ab = (chi - 1) / (chi + 1);
				a_ab = u_ab / (u_ab * u_ab - 1.0);

				// JCP 88, 2547 (1988), eq. A3
				if (a_ab > 0.5)
					a_ab = 0.5;
				if (a_ab < -0.5)
					a_ab = -0.5;

				nu_ab = mu_ab + a_ab * (1.0 - mu_ab * mu_ab);
			}
			else
				nu_ab = mu_ab;

			f = gpu_f3(nu_ab);

			if (abs(1.0 - f) < 1.0e-14)
				// if f == 1.0 we need to take care
				// otherwise we can get numerical problems
				pa[a] = 0.0;
			else {
				if (pa[a] > 1E-20 || pa[a] < -1E-20)
					pa[a] *= 0.5 * (1.0 - f);
				else
					pa[a] = 0.0;
				if (pa[b] > 1E-20 || pa[b] < -1E-20)
					pa[b] *= 0.5 * (1.0 + f);
				else
					pa[b] = 0.0;
			}
		}
	}

	float w = 0.0;
	for (int a = 0; a < num_centers; a++)
		w += pa[a];

	float res = 1.0;
	if (abs(w) > 1.0e-14)
		res = pa[center_index] / w;

	free(pa);

	return res;
}

__global__ void gpu_make_grid(
	const int center_index,
	const float* x,
	const float* y,
	const float* z,
	const double* atom_grid_x_bohr_,
	const double* atom_grid_y_bohr_,
	const double* atom_grid_z_bohr_,
	const double* atom_grid_w,
	const int* proton_charges,
	const int* asym_atom_list,
	const int* num_points,
	const int offset,
	float* grid_x_bohr,
	float* grid_y_bohr,
	float* grid_z_bohr,
	float* grid_aw,
	float* grid_mw)
{
	const int ipoint = blockIdx.x * blockDim.x + threadIdx.x;
	if (ipoint >= num_points[center_index]) return;
	if (ipoint + offset >= gpu_MaxGrid[0])
		return;
	const float x_bohr = atom_grid_x_bohr_[ipoint] + x[center_index];
	const float y_bohr = atom_grid_y_bohr_[ipoint] + y[center_index];
	const float z_bohr = atom_grid_z_bohr_[ipoint] + z[center_index];
	const float w = atom_grid_w[ipoint];

	grid_mw[ipoint + offset] = w * get_becke_w(gpu_ncen[0],
		proton_charges,
		x,
		y,
		z,
		asym_atom_list[center_index],
		x_bohr,
		y_bohr,
		z_bohr);
	grid_aw[ipoint + offset] = w;
	grid_x_bohr[ipoint + offset] = x_bohr;
	grid_y_bohr[ipoint + offset] = y_bohr;
	grid_z_bohr[ipoint + offset] = z_bohr;
}

__global__ void gpu_linear_interpolate_spherical_density(
	const int atom,
	const float* radial_dens,
	const float* spherical_dist,
	const int size,
	const bool match,
	const int offset,
	const float* gridx,
	const float* gridy,
	const float* gridz,
	const int* num_points,
	float* spherical_density,
	float* Grids,
	const float* posx,
	const float* posy,
	const float* posz
)
{
	const int indice = blockIdx.x * blockDim.x + threadIdx.x;
	if (indice >= num_points[atom]) return;
	if (indice + offset >= gpu_MaxGrid[0]) return;
	float dist[3] = { gridx[indice + offset] - posx[atom], gridy[indice + offset] - posy[atom], gridz[indice + offset] - posz[atom] };
	//storing result in dist[0] to save memory
	dist[0] = sqrt(dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2]);
	float result;
	if (dist[0] > spherical_dist[size - 1])
		return;
	else if (dist[0] < spherical_dist[0])
		result = radial_dens[0];
	else {
		float step0 = dist[0] / gpu_start_radial_dens[0];
		step0 = logf(step0);
		step0 /= gpu_log_incr[0];
		int step = floorf(step0);
		if (step > size)
			result = spherical_dist[size - 1];
		else {
			result = radial_dens[step] + (radial_dens[step + 1] - radial_dens[step]) / (spherical_dist[step] - spherical_dist[step - 1]) * (dist[0] - spherical_dist[step - 1]);
			if (result < 1E-10) return;
		}
	}
	if (match) spherical_density[indice] = result;
	Grids[indice + offset] += result;
}

__global__ void gpu_apply_weights(
	const int offset,
	const float* GridRho,
	const float* Grids,
	float* Grid_spherical_atom,
	const float* grid_aw,
	const int number_of_points)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= number_of_points) return;
	const int j = i + offset;
	if (j >= gpu_MaxGrid[0])
		return;
	if (Grid_spherical_atom[i] != 0) {
		if (Grids[j] != 0)
			Grid_spherical_atom[i] *= GridRho[j] * grid_aw[j] / Grids[j];
		else
			Grid_spherical_atom[i] = 0;
	}
}

__global__ void gpu_count_non_zero(
	const int len,
	const double* Array,
	int* result)
{
	int counter = 0;
	for (int i = 0; i < len; i++)
		if (Array[i] != 0) counter++;
	result[0] = counter;
}

__global__ void gpu_calc_charges(
	const float* GridRho,
	const float* Gridaw,
	const float* Gridmw,
	const float* Grids,
	const float* spherical_density,
	const int num_points,
	const int offset,
	const float cutoff,
	double* result)
{
	// Vector containing integrated numbers of electrons
	// dimension 0: 0=Becke grid integration 1=Summed spherical density 2=hirshfeld weighted density
	double atom_els[3] = { 0.0, 0.0, 0.0 };

	//Generate Electron sums
	for (int p = 0; p < num_points; p++) {
		int big_p = p + offset;
		if (Grids[big_p] > cutoff) {
			atom_els[0] += Gridmw[big_p] * GridRho[big_p];
			atom_els[1] += Gridmw[big_p] * Grids[big_p];
			if (Grids[big_p] != 0)
				atom_els[2] += GridRho[big_p] * Gridaw[p] * spherical_density[p] / Grids[big_p];
		}
	}
	for (int i = 0; i < 3; i++)
		result[i] = atom_els[i];
}

__global__ void gpu_calc_dens_per_MO_shared(
	float* __restrict__ GridRho,
	const float* __restrict__ Gridx,
	const float* __restrict__ Gridy,
	const float* __restrict__ Gridz,
	const float* __restrict__ x,
	const float* __restrict__ y,
	const float* __restrict__ z,
	const int* __restrict__ types,
	const int* __restrict__ centers,
	const float* __restrict__ exponents,
	const float* __restrict__ coefficients,
	const float occ)
{
	const int indice = blockIdx.x * blockDim.x + threadIdx.x;
	if (indice >= gpu_MaxGrid[0])
		return;
	const float Grid[3]{ Gridx[indice], Gridy[indice], Gridz[indice] };
	const int max_center = gpu_ncen[0];
	extern __shared__ float shared[];
	float* d0 = (float*)malloc(sizeof(float) * max_center);
	float* d1 = (float*)malloc(sizeof(float) * max_center);
	float* d2 = (float*)malloc(sizeof(float) * max_center);
	float* dsqd = (float*)malloc(sizeof(float) * max_center);
	float ex, phi = 0.0;
	int i, iat;

	if (threadIdx.x == 0)
		for (i = 0; i < gpu_nex[0]; i++)
			shared[i] = coefficients[i];
	__syncthreads();

	for (i = 0; i < max_center; i++) {
		d0[i] = Grid[0] - x[i];
		d1[i] = Grid[1] - y[i];
		d2[i] = Grid[2] - z[i];
		dsqd[i] = d0[i] * d0[i] + d1[i] * d1[i] + d2[i] * d2[i];
	}

	for (i = 0; i < gpu_nex[0]; i++) {
		iat = centers[i] - 1;
		ex = -exponents[i] * dsqd[iat];
		ex = expf(ex);
		ex *= gpu_type2mult(types[i], d0[iat], d1[iat], d2[iat]);
		phi += shared[i] * ex;      //build MO value at this point
	}

	free(d0);
	free(d1);
	free(d2);
	free(dsqd);

	GridRho[indice] += occ * phi * phi;
}

__global__ void gpu_calc_dens_per_MO_static1(
	float* __restrict__ GridRho,
	const float* __restrict__ Gridx,
	const float* __restrict__ Gridy,
	const float* __restrict__ Gridz,
	const float* __restrict__ x,
	const float* __restrict__ y,
	const float* __restrict__ z,
	const int* __restrict__ types,
	const int* __restrict__ centers,
	const float* __restrict__ exponents,
	const float* __restrict__ coefficients,
	const float occ)
{
	const int indice = blockIdx.x * blockDim.x + threadIdx.x;
	if (indice >= gpu_MaxGrid[0])
		return;
	int iat, m;
	const float Grid[3]{ Gridx[indice], Gridy[indice], Gridz[indice] };
	float d0, d1, d2, ex;
	int l0, l1, l2, i;
	float phi = 0.0;

	for (i = 0; i < gpu_nex[0]; i++) {
		iat = centers[i];
		d0 = Grid[0] - x[iat];
		d1 = Grid[1] - y[iat];
		d2 = Grid[2] - z[iat];
		ex = -exponents[i] * d0 * d0 + d1 * d1 + d2 * d2;
		ex = expf(ex);
		gpu_type2vector(types[i], l0, l1, l2);
		for (m = 0; m < l0; m++)
			ex *= d0;
		for (m = 0; m < l1; m++)
			ex *= d1;
		for (m = 0; m < l2; m++)
			ex *= d2;
		phi += coefficients[i] * ex;      //build MO value at this point
	}

	GridRho[indice] += occ * phi * phi;
}

__global__ void gpu_calc_dens_per_MO_static2(
	float* __restrict__ GridRho,
	const float* __restrict__ Gridx,
	const float* __restrict__ Gridy,
	const float* __restrict__ Gridz,
	const float* __restrict__ x,
	const float* __restrict__ y,
	const float* __restrict__ z,
	const int* __restrict__ types,
	const int* __restrict__ centers,
	const float* __restrict__ exponents,
	const float* __restrict__ coefficients,
	const float occ)
{
	const int indice = blockIdx.x * blockDim.x + threadIdx.x;
	if (indice >= gpu_MaxGrid[0])
		return;
	int iat, m;
	const float Grid[3]{ Gridx[indice], Gridy[indice], Gridz[indice] };
	const int max_center = gpu_ncen[0];
	float* d0 = (float*)malloc(sizeof(float) * max_center);
	float* d1 = (float*)malloc(sizeof(float) * max_center);
	float* d2 = (float*)malloc(sizeof(float) * max_center);
	float* dsqd = (float*)malloc(sizeof(float) * max_center);
	float ex;
	int l0, l1, l2, i;
	float phi = 0.0;

	for (i = 0; i < max_center; i++) {
		d0[i] = Grid[0] - x[i];
		d1[i] = Grid[1] - y[i];
		d2[i] = Grid[2] - z[i];
	}

	for (i = 0; i < gpu_nex[0]; i++) {
		iat = centers[i];
		ex = -exponents[i] * dsqd[iat];
		ex = expf(ex);
		gpu_type2vector(types[i], l0, l1, l2);
		for (m = 0; m < l0; m++)
			ex *= d0[iat];
		for (m = 0; m < l1; m++)
			ex *= d1[iat];
		for (m = 0; m < l2; m++)
			ex *= d2[iat];
		phi += coefficients[i] * ex;      //build MO value at this point
	}

	GridRho[indice] += occ * phi * phi;
}


__global__ void gpu_calc_dens(
	float* GridRho,
	const float* Gridx,
	const float* Gridy,
	const float* Gridz,
	const float* x,
	const float* y,
	const float* z,
	const int* types,
	const int* centers,
	const float* exponents,
	const float* coefficients,
	const float* occ)
{
	const int indice = blockIdx.x * blockDim.x + threadIdx.x;
	if (indice >= gpu_MaxGrid[0])
		return;
	const float Grid[3]{ Gridx[indice], Gridy[indice], Gridz[indice] };
	int iat, l0, l1, l2, i, m;
	const int max_center = gpu_ncen[0];
	float* d0 = (float*)malloc(sizeof(float) * max_center);
	float* d1 = (float*)malloc(sizeof(float) * max_center);
	float* d2 = (float*)malloc(sizeof(float) * max_center);
	float* dsq = (float*)malloc(sizeof(float) * max_center);
	float ex;
	float Rho = 0;
	const int max_mo = gpu_nmo[0];
	float* phi = (float*)malloc(sizeof(float) * max_mo);
	int max_moi = 0;
	for (i = 0; i < max_mo; i++)
		phi[i] = 0.0;

	for (i = 0; i < max_center; i++) {
		d0[i] = Grid[0] - x[i];
		d1[i] = Grid[1] - y[i];
		d2[i] = Grid[2] - z[i];
		dsq[i] = d0[i] * d0[i] + d1[i] * d1[i] + d2[i] * d2[i];
	}

	for (i = 0; i < gpu_nex[0]; i++) {
		iat = centers[i] - 1;
		ex = -exponents[i] * dsq[iat];
		ex = expf(ex);
		gpu_type2vector(types[i], l0, l1, l2);
		for (m = 0; m < l0; m++)
			ex *= d0[iat];
		for (m = 0; m < l1; m++)
			ex *= d1[iat];
		for (m = 0; m < l2; m++)
			ex *= d2[iat];
		max_moi += max_mo;
		for (m = 0; m < max_mo; m++)
			phi[m] += coefficients[m + max_moi] * ex;      //build MO values at this point
	}

	for (m = 0; m < max_mo; m++)
		Rho += occ[m] * phi[m] * phi[m];

	GridRho[indice] = Rho;
	free(phi);
}

void printCUDA(const cudaDeviceProp* prop, const int nDevices, ostream& file)
{
	file << endl << "Number of CUDA-capable devices : " << nDevices << endl;
	/*table with the properties of the devices*/

	for (int dev = 0; dev < nDevices; dev++) {
		file << "Device Number: " << dev << endl;
		file << "- Device name: " << (prop[dev]).name << endl;
		file << "- Max global memory   = " << prop[dev].totalGlobalMem / 1048576.0f << " MBytes" << endl;
		file << "- Max constant memory = " << prop[dev].totalConstMem / 1048576.0f << " MBytes" << endl;
		file << "- Capability : " << prop[dev].major << "." << prop[dev].minor << endl;
		file << "- Total global memory: " << (float)prop[dev].totalGlobalMem / 1048576.0f << " MByte" << endl;
		file << "_____________________________________________________" << endl;
	}
}

__host__ int cu_make_hirshfeld_grids(
    const int& pbc,
    const int& accuracy,
    cell& unit_cell,
    const WFN& wave,
    const ivec& atom_type_list,
    const ivec& cif2wfn_list,
    bvec& needs_grid,
    vec2& d1,
    vec2& d2,
    vec2& d3,
    vec2& dens,
    const vector<string>& labels,
    ostream& file,
    time_point& start,
    time_point& end_becke,
    time_point& end_prototypes,
    time_point& end_spherical,
    time_point& end_prune,
    time_point& end_aspherical,
    bool debug,
    bool no_date)
{
    using namespace std;
    int atoms_with_grids = vec_sum(needs_grid);
    err_checkf(atoms_with_grids > 0, "No atoms with grids to generate!", file);
    err_checkf(atoms_with_grids <= wave.get_ncen(), "More atoms with grids than in the wavefunction! Aborting!", file);
    err_checkf(atoms_with_grids == cif2wfn_list.size(), "Number of atoms with grids does not match the number of atoms in the CIF file!", file);
    vec3 grid(atoms_with_grids);
    const int nr_of_atoms = (wave.get_ncen() * (int)pow(pbc * 2 + 1, 3));
    vec x(nr_of_atoms), y(nr_of_atoms), z(nr_of_atoms);
    ivec atom_z(nr_of_atoms);
#pragma omp parallel for
    for (int i = 0; i < atoms_with_grids; i++)
        grid[i].resize(6);
    // GRID COORDINATES for [a][c][p] a = atom [0,ncen],
    // c = coordinate [0=x, 1=y, 2=z, 3=atomic becke weight, 4= molecular becke weight, 5=total spherical density],
    // p = point in this grid

    fill_xyzc(x, y, z, atom_z, pbc, wave, unit_cell);
    int nDevices;                  /*Number of devices available (running time)*/
    cudaGetDeviceCount(&nDevices); /*Get the number of devices*/
    int dev = 0;
    cudaDeviceProp* prop = NULL;
    cudaDeviceProp deviceProp;
    prop = (cudaDeviceProp*)malloc(sizeof(cudaDeviceProp) * nDevices);
    for (int devl = 0; devl < nDevices; devl++)
    { // Make CUDA information available in prop
        cudaGetDeviceProperties(&(deviceProp), devl);
        prop[devl] = deviceProp;
    }
    cudaSetDevice(dev);
    printCUDA(prop, nDevices, file);
    gpuErrchk(cudaSetDeviceFlags(cudaDeviceScheduleBlockingSync));

    // Make Prototype grids with only single atom weights for all elements
    vector<AtomGrid> Prototype_grids = make_Prototype_atoms(atom_type_list, cif2wfn_list, debug, file, accuracy, wave, 1);
    if (no_date)
        file << "\nMaking Becke Grids..." << flush;
    else
    {
        file << endl << "Selected accuracy: " << accuracy << "\nMaking Becke Grids..." << flush;
    }
    end_prototypes = get_time();
    if (debug)
    {
        file << "\n";
        for (int prototype = 0; prototype < Prototype_grids.size(); prototype++)
            file << "Number of gridpoints for atom type " << atom_type_list[prototype] << ": " << Prototype_grids[prototype].get_num_grid_points() << endl;

        long long int dur = get_sec(start, end_prototypes);

        if (dur < 1)
            file << "Time until prototypes are done: " << fixed << setprecision(0) << get_msec(start, end_prototypes) << " ms" << endl;
        else
            file << "Time until prototypes are done: " << fixed << setprecision(0) << dur << " s" << endl;
    }

    vector<vector<double>> radial_density;
    vector<vector<double>> radial_dist;

    radial_density.resize(atom_type_list.size());
    radial_dist.resize(atom_type_list.size());
    vector<vector<double>> spherical_density(atoms_with_grids);
    spherical_density.resize(asym_atom_list.size());
    for (int i = 0; i < asym_atom_list.size(); i++)
        spherical_density[i].resize(num_points[i]);

    const double incr = pow(1.005, max(1, accuracy - 1));
    const double lincr = log(incr);
    const double min_dist = 0.0000001;
    vector<Thakkar> sphericals;
    for (int i = 0; i < atom_type_list.size(); i++)
        sphericals.push_back(Thakkar(atom_type_list[i]));
    // Make radial grids
    if (debug)
    {
        for (int i = 0; i < atom_type_list.size(); i++)
        {
            file << "Calculating for atomic number " << atom_type_list[i] << endl;
            double current = 1;
            double dist = min_dist;
            if (accuracy > 3)
                while (current > 1E-10)
                {
                    radial_dist[i].push_back(dist);
                    current = sphericals[i].get_radial_density(dist);
                    if (current == -20)
                        return false;
                    radial_density[i].push_back(current);
                    dist *= incr;
                }
            else
                while (current > 1E-12)
                {
                    radial_dist[i].push_back(dist);
                    current = sphericals[i].get_radial_density(dist);
                    if (current == -20)
                        return false;
                    radial_density[i].push_back(current);
                    dist *= incr;
                }
            file << "Number of radial density points for atomic number " << atom_type_list[i] << ": " << radial_density[i].size() << endl;
            // for (int j = 0; j < radial_density[i].size(); j++) {
            //	if (radial_density[i][j] < 0.1)
            //		break;
            //	file << scientific << setprecision(8) << radial_density[i][j] << endl;
            // }
        }
    }
    else
    {
#pragma omp parallel for
        for (int i = 0; i < atom_type_list.size(); i++)
        {
            double current = 1;
            double dist = min_dist;
            if (accuracy > 3)
                while (current > 1E-10)
                {
                    radial_dist[i].push_back(dist);
                    current = sphericals[i].get_radial_density(dist);
                    radial_density[i].push_back(current);
                    dist *= incr;
                }
            else
                while (current > 1E-12)
                {
                    radial_dist[i].push_back(dist);
                    current = sphericals[i].get_radial_density(dist);
                    radial_density[i].push_back(current);
                    dist *= incr;
                }
        }
    }
    sphericals.clear();

    float** gpu_PosAtomsx = NULL,
        ** gpu_PosAtomsy = NULL,
        ** gpu_PosAtomsz = NULL,
        ** gpu_GridRho = NULL,
        ** gpu_Gridx = NULL,
        ** gpu_Gridy = NULL,
        ** gpu_Gridz = NULL,
        ** gpu_Gridaw = NULL,
        ** gpu_Gridmw = NULL,
        ** gpu_exponents = NULL,
        ** gpu_coefficients = NULL,
        ** gpu_occ = NULL;
    double*** gpu_atomgrid_x = NULL,
        *** gpu_atomgrid_y = NULL,
        *** gpu_atomgrid_z = NULL,
        *** gpu_atomgrid_w = NULL;
    int** gpu_types = NULL,
        ** gpu_centers = NULL,
        ** gpu_asym_atom_list = NULL,
        ** gpu_atom_type_list = NULL,
        ** gpu_numpoints = NULL,
        ** gpu_atom_z = NULL;
    vector<vector<float>> PosAtoms;
    PosAtoms.resize(3);
    for (int i = 0; i < 3; i++)
        PosAtoms[i].resize(wave.get_ncen());
    for (int a = 0; a < wave.get_ncen(); a++)
    {
        PosAtoms[0][a] = wave.atoms[a].x;
        PosAtoms[1][a] = wave.atoms[a].y;
        PosAtoms[2][a] = wave.atoms[a].z;
    }
    /*Allocation GPU Pointer*/
    gpu_PosAtomsx = (float**)malloc(sizeof(float*));
    gpu_PosAtomsy = (float**)malloc(sizeof(float*));
    gpu_PosAtomsz = (float**)malloc(sizeof(float*));
    gpu_atom_z = (int**)malloc(sizeof(int*));
    gpu_types = (int**)malloc(sizeof(int*));
    gpu_centers = (int**)malloc(sizeof(int*));
    gpu_asym_atom_list = (int**)malloc(sizeof(int*));
    gpu_atom_type_list = (int**)malloc(sizeof(int*));
    gpu_numpoints = (int**)malloc(sizeof(int*));
    gpu_exponents = (float**)malloc(sizeof(float*));
    gpu_occ = (float**)malloc(sizeof(float*));
    gpu_coefficients = (float**)malloc(sizeof(float*));
    gpu_GridRho = (float**)malloc(sizeof(float*));
    gpu_Gridx = (float**)malloc(sizeof(float*));
    gpu_Gridy = (float**)malloc(sizeof(float*));
    gpu_Gridz = (float**)malloc(sizeof(float*));
    gpu_Gridaw = (float**)malloc(sizeof(float*));
    gpu_Gridmw = (float**)malloc(sizeof(float*));
    gpu_atomgrid_x = (double***)malloc(sizeof(double**));
    gpu_atomgrid_y = (double***)malloc(sizeof(double**));
    gpu_atomgrid_z = (double***)malloc(sizeof(double**));
    gpu_atomgrid_w = (double***)malloc(sizeof(double**));
    gpu_atomgrid_x[0] = (double**)malloc(sizeof(double*) * asym_atom_list.size());
    gpu_atomgrid_y[0] = (double**)malloc(sizeof(double*) * asym_atom_list.size());
    gpu_atomgrid_z[0] = (double**)malloc(sizeof(double*) * asym_atom_list.size());
    gpu_atomgrid_w[0] = (double**)malloc(sizeof(double*) * asym_atom_list.size());

    int nex_temp = wave.get_nex();
    int nmo_temp = wave.get_nmo(true);
    int ncen_temp = wave.get_ncen();
    for (int i = 0; i < wave.get_ncen(); i++)
        atom_z[i] = wave.atoms[i].charge;
    int MaxGrid = 0;
    for (int i = 0; i < asym_atom_list.size(); i++)
    {
        int nr = asym_atom_list[i];
        int type;
        for (int j = 0; j < atom_type_list.size(); j++)
            if (atom_type_list[j] == wave.atoms[nr].charge)
                type = j;

        num_points[i] = Prototype_grids[type].get_num_grid_points();
        MaxGrid += num_points[i];
    }
    int numBlocks, blocks, gridSize;
    size_t size;
    gpuErrchk(cudaDeviceGetLimit(&size, cudaLimitMallocHeapSize));
    if (debug)
        file << "\nold Heapsize: " << size / 1024 / 1024 << " MB" << endl;
    float result;
    gpuErrchk(cudaOccupancyMaxPotentialBlockSize(
        &numBlocks,
        &blocks,
        (void*)gpu_calc_dens_per_MO_static1,
        0,
        MaxGrid));
    gridSize = (MaxGrid + blocks - 1) / blocks;
    result = ((sizeof(int) * 10 + sizeof(float) * (6 * ncen_temp + 20)) * blocks * prop[dev].multiProcessorCount) / 1024 / 1024;

    if (debug)
        file << "result: " << fixed << result << " MB for " << gridSize << " blocks with " << blocks << " threads" << endl
        << "sizeof float: " << sizeof(float) << " B" << endl;
    if (result > size / 1024 / 1024)
        gpuErrchk(cudaDeviceSetLimit(cudaLimitMallocHeapSize, result * 1024 * 1024 * 1.1));

    gpuErrchk(cudaDeviceGetLimit(&size, cudaLimitMallocHeapSize));
    if (debug)
        file << "new Heapsize: " << size / 1024 / 1024 << " MB" << endl;

    gridSize = (MaxGrid + blocks - 1) / blocks;

    // Allocate and copy vectors to device for WFN
    if (debug)
        file << "Copying WFN to devices now!" << endl;
    gpuErrchk(cudaMalloc((void**)&gpu_PosAtomsx[0], sizeof(float) * ncen_temp));
    gpuErrchk(cudaMalloc((void**)&gpu_PosAtomsy[0], sizeof(float) * ncen_temp));
    gpuErrchk(cudaMalloc((void**)&gpu_PosAtomsz[0], sizeof(float) * ncen_temp));
    gpuErrchk(cudaMalloc((void**)&gpu_atom_z[0], sizeof(int) * ncen_temp));
    gpuErrchk(cudaMalloc((void**)&gpu_Gridx[0], sizeof(float) * MaxGrid));
    gpuErrchk(cudaMalloc((void**)&gpu_Gridy[0], sizeof(float) * MaxGrid));
    gpuErrchk(cudaMalloc((void**)&gpu_Gridz[0], sizeof(float) * MaxGrid));
    gpuErrchk(cudaMalloc((void**)&gpu_Gridaw[0], sizeof(float) * MaxGrid));
    gpuErrchk(cudaMalloc((void**)&gpu_Gridmw[0], sizeof(float) * MaxGrid));
    gpuErrchk(cudaMalloc((void**)&gpu_asym_atom_list[0], sizeof(int) * asym_atom_list.size()));
    gpuErrchk(cudaMalloc((void**)&gpu_atom_type_list[0], sizeof(int) * atom_type_list.size()));
    gpuErrchk(cudaMalloc((void**)&gpu_numpoints[0], sizeof(int) * asym_atom_list.size()));
    if (debug)
        file << "Mallocs done!" << endl;
    gpuErrchk(cudaMemcpyToSymbol(gpu_nex, &nex_temp, sizeof(int)));
    gpuErrchk(cudaMemcpyToSymbol(gpu_nmo, &nmo_temp, sizeof(int)));
    gpuErrchk(cudaMemcpyToSymbol(gpu_ncen, &ncen_temp, sizeof(int)));
    gpuErrchk(cudaMemcpyToSymbol(gpu_MaxGrid, &MaxGrid, sizeof(int)));
    gpuErrchk(cudaMemcpyToSymbol(gpu_start_radial_dens, &min_dist, sizeof(float)));
    gpuErrchk(cudaMemcpyToSymbol(gpu_log_incr, &lincr, sizeof(float)));
    gpuErrchk(cudaMemcpy(gpu_PosAtomsx[0], PosAtoms[0].data(), sizeof(float) * wave.get_ncen(), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(gpu_PosAtomsy[0], PosAtoms[1].data(), sizeof(float) * wave.get_ncen(), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(gpu_PosAtomsz[0], PosAtoms[2].data(), sizeof(float) * wave.get_ncen(), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(gpu_asym_atom_list[0], asym_atom_list.data(), sizeof(int) * asym_atom_list.size(), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(gpu_atom_type_list[0], atom_type_list.data(), sizeof(int) * atom_type_list.size(), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(gpu_numpoints[0], num_points.data(), sizeof(int) * asym_atom_list.size(), cudaMemcpyHostToDevice));
    gpuErrchk(cudaPeekAtLastError());
    file << "All copying done!" << endl;

    vector<cudaStream_t> streams;
    streams.resize(asym_atom_list.size());
    int offset = 0;
    for (int i = 0; i < asym_atom_list.size(); i++)
    {
        gpuErrchk(cudaOccupancyMaxPotentialBlockSize(
            &numBlocks,
            &blocks,
            (void*)gpu_make_grid,
            0,
            num_points[i]));

        gridSize = (MaxGrid + blocks - 1) / blocks;
        if (debug)
            file << i << ": num points: " << num_points[i] << " blocks: " << gridSize << " threads: " << blocks << " grid > points? " << (blocks * gridSize > num_points[i]) << endl;
        int type;
        for (int j = 0; j < atom_type_list.size(); j++)
            if (atom_type_list[j] == wave.atoms[asym_atom_list[i]].charge)
                type = j;
        gpuErrchk(cudaStreamCreate(&streams[i]));
        gpuErrchk(cudaMalloc((void**)&gpu_atomgrid_x[0][i], sizeof(double) * num_points[i]));
        gpuErrchk(cudaMalloc((void**)&gpu_atomgrid_y[0][i], sizeof(double) * num_points[i]));
        gpuErrchk(cudaMalloc((void**)&gpu_atomgrid_z[0][i], sizeof(double) * num_points[i]));
        gpuErrchk(cudaMalloc((void**)&gpu_atomgrid_w[0][i], sizeof(double) * num_points[i]));
        gpuErrchk(cudaMemcpyAsync(gpu_atomgrid_x[0][i], Prototype_grids[type].get_gridx_ptr(), sizeof(double) * num_points[i], cudaMemcpyHostToDevice, streams[i]));
        gpuErrchk(cudaMemcpyAsync(gpu_atomgrid_y[0][i], Prototype_grids[type].get_gridy_ptr(), sizeof(double) * num_points[i], cudaMemcpyHostToDevice, streams[i]));
        gpuErrchk(cudaMemcpyAsync(gpu_atomgrid_z[0][i], Prototype_grids[type].get_gridz_ptr(), sizeof(double) * num_points[i], cudaMemcpyHostToDevice, streams[i]));
        gpuErrchk(cudaMemcpyAsync(gpu_atomgrid_w[0][i], Prototype_grids[type].get_gridw_ptr(), sizeof(double) * num_points[i], cudaMemcpyHostToDevice, streams[i]));

        gpu_make_grid << <gridSize, blocks, 0, streams[i] >> > (
            i,
            gpu_PosAtomsx[0],
            gpu_PosAtomsy[0],
            gpu_PosAtomsz[0],
            gpu_atomgrid_x[0][i],
            gpu_atomgrid_y[0][i],
            gpu_atomgrid_z[0][i],
            gpu_atomgrid_w[0][i],
            gpu_atom_z[0],
            gpu_asym_atom_list[0],
            gpu_numpoints[0],
            offset,
            gpu_Gridx[0],
            gpu_Gridy[0],
            gpu_Gridz[0],
            gpu_Gridaw[0],
            gpu_Gridmw[0]);

        offset += num_points[i];
    }
    gpuErrchk(cudaDeviceSynchronize());
    gpuErrchk(cudaPeekAtLastError());

    for (int i = 0; i < atom_type_list.size(); i++)
    {
        gpuErrchk(cudaFree(gpu_atomgrid_w[0][i]));
        gpuErrchk(cudaFree(gpu_atomgrid_x[0][i]));
        gpuErrchk(cudaFree(gpu_atomgrid_y[0][i]));
        gpuErrchk(cudaFree(gpu_atomgrid_z[0][i]));
    }

    int points = vec_sum(num_points);
    if (debug)
        file << "Becke Grid exists" << endl;
    else
        file << "                           done! Number of gridpoints: " << defaultfloat << points << endl;

    end_becke = get_time();

    file << "Calculating spherical densities..." << flush;

    vector<vector<float>> total_grid(7);
    float*** gpu_spherical_density = NULL,
        *** gpu_radial_density = NULL,
        *** gpu_radial_dist = NULL,
        ** gpu_Grids = NULL;
    gpu_radial_density = (float***)malloc(sizeof(float**));
    gpu_radial_dist = (float***)malloc(sizeof(float**));
    gpu_spherical_density = (float***)malloc(sizeof(float**));

    gpu_radial_density[0] = (float**)malloc(sizeof(float*) * atom_type_list.size());
    gpu_radial_dist[0] = (float**)malloc(sizeof(float*) * atom_type_list.size());
    gpu_spherical_density[0] = (float**)malloc(sizeof(float*) * asym_atom_list.size());
    gpu_Grids = (float**)malloc(sizeof(float*));

    for (int i = 0; i < asym_atom_list.size(); i++)
        gpuErrchk(cudaMalloc((void**)&(gpu_spherical_density[0][i]), sizeof(float) * num_points[i]));

    gpuErrchk(cudaMalloc((void**)&gpu_Grids[0], sizeof(float) * MaxGrid));

    for (int i = 0; i < atom_type_list.size(); i++)
    {
        gpuErrchk(cudaMalloc((void**)&gpu_radial_density[0][i], sizeof(float) * radial_density[i].size()));
        gpuErrchk(cudaMalloc((void**)&gpu_radial_dist[0][i], sizeof(float) * radial_dist[i].size()));
        gpuErrchk(cudaMemcpy(gpu_radial_density[0][i], radial_density[i].data(), sizeof(float) * radial_density[i].size(), cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(gpu_radial_dist[0][i], radial_dist[i].data(), sizeof(float) * radial_dist[i].size(), cudaMemcpyHostToDevice));
    }
    gpuErrchk(cudaDeviceSynchronize());
    gpuErrchk(cudaPeekAtLastError());
    for (int i = 0; i < wave.get_ncen(); i++)
    {
        int nr = all_atom_list[i];
        int type_list_number = -1;
        for (int j = 0; j < atom_type_list.size(); j++)
            if (wave.atoms[nr].charge == atom_type_list[j])
                type_list_number = j;
        offset = 0;
        for (int g = 0; g < asym_atom_list.size(); g++)
        {
            if (g == 0)
            {
                gpuErrchk(cudaOccupancyMaxPotentialBlockSize(
                    &numBlocks,
                    &blocks,
                    (void*)gpu_linear_interpolate_spherical_density,
                    0,
                    num_points[i]));

                gridSize = (MaxGrid + blocks - 1) / blocks;
                // file << i << "/" << g << ": blocks: " << gridSize << " threads: " << blocks << endl;
            }
            bool match = (all_atom_list[i] == asym_atom_list[g]);
            gpu_linear_interpolate_spherical_density << <gridSize, blocks >> > (
                i,
                gpu_radial_density[0][type_list_number],
                gpu_radial_dist[0][type_list_number],
                radial_density[type_list_number].size(),
                match,
                offset,
                gpu_Gridx[0],
                gpu_Gridy[0],
                gpu_Gridz[0],
                gpu_numpoints[0],
                gpu_spherical_density[0][g],
                gpu_Grids[0],
                gpu_PosAtomsx[0],
                gpu_PosAtomsy[0],
                gpu_PosAtomsz[0]);

            offset += num_points[i];
        }
    }
    gpuErrchk(cudaDeviceSynchronize());
    gpuErrchk(cudaPeekAtLastError());

    if (debug)
    {
        // Copy grid from GPU to print:
        //  Dimensions: [c] [p]
        //  p = the number of gridpoint
        //  c = coordinate, which is 0=x, 1=y, 2=z, 3=atomic becke weight, 4=spherical density, 5=wavefunction density, 6=molecular becke weight
        for (int i = 0; i < 7; i++)
            total_grid[i].resize(MaxGrid);
        gpuErrchk(cudaMemcpy(total_grid[0].data(), gpu_Gridx[0], sizeof(float) * MaxGrid, cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(total_grid[1].data(), gpu_Gridy[0], sizeof(float) * MaxGrid, cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(total_grid[2].data(), gpu_Gridz[0], sizeof(float) * MaxGrid, cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(total_grid[3].data(), gpu_Gridaw[0], sizeof(float) * MaxGrid, cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(total_grid[6].data(), gpu_Gridmw[0], sizeof(float) * MaxGrid, cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(total_grid[4].data(), gpu_Grids[0], sizeof(float) * MaxGrid, cudaMemcpyDeviceToHost));

        gpuErrchk(cudaDeviceSynchronize());
        gpuErrchk(cudaPeekAtLastError());

        ofstream grid0("grid0.file", ios::out);
        for (int i = 0; i < MaxGrid; i++)
        {
            for (int j = 0; j < 7; j++)
                grid0 << setw(16) << scientific << setprecision(8) << total_grid[j][i];
            grid0 << "\n";
        }
        grid0.flush();
        grid0.close();
    }

    for (int i = 0; i < atom_type_list.size(); i++)
    {
        gpuErrchk(cudaFree(gpu_radial_density[0][i]));
        gpuErrchk(cudaFree(gpu_radial_dist[0][i]));
    }

    file << "                    done!" << endl;

    shrink_vector<vec>(radial_density);
    shrink_vector<vec>(radial_dist);

    end_spherical = get_time();


    end_prune = get_time();

    file << "Calculating non-spherical densities..." << flush;
    vec2 periodic_grid;

    // Vector containing integrated numbers of electrons
    // dimension 0: 0=Becke grid integration 1=Summed spherical density 2=hirshfeld weighted density
    // dimension 1: atoms of asym_atom_list
    vector<vector<double>> atom_els;
    atom_els.resize(3);
    for (int i = 0; i < asym_atom_list.size(); i++)
        for (int n = 0; n < 3; n++)
            atom_els[n].push_back(0.0);

    vector<float> coef;
    vector<float> ex;
    for (int i = 0; i < nex_temp; i++)
        ex.push_back(wave.get_exponent(i));

    gpuErrchk(cudaMalloc((void**)&gpu_types[0], sizeof(int) * nex_temp));
    gpuErrchk(cudaMalloc((void**)&gpu_centers[0], sizeof(int) * nex_temp));
    gpuErrchk(cudaMalloc((void**)&gpu_exponents[0], sizeof(float) * nex_temp));
    gpuErrchk(cudaMalloc((void**)&gpu_GridRho[0], sizeof(float) * MaxGrid));
    gpuErrchk(cudaMemcpy(gpu_types[0], wave.get_ptr_types(), sizeof(int) * nex_temp, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(gpu_centers[0], wave.get_ptr_centers(), sizeof(int) * nex_temp, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(gpu_exponents[0], ex.data(), sizeof(float) * nex_temp, cudaMemcpyHostToDevice));

    const bool per_MO = true;
    if (per_MO)
    {
        for (int mo = 0; mo < wave.get_nmo(false); mo++)
            for (int i = 0; i < nex_temp; i++)
                if (wave.get_MO_occ(mo) != 0)
                    coef.push_back(wave.get_MO_coef(mo, i));
        if (debug)
            file << "Number of coefs: " << coef.size() << endl;
        gpuErrchk(cudaMalloc((void**)&gpu_coefficients[0], sizeof(float) * nex_temp));
        gpuErrchk(cudaMemset(gpu_GridRho[0], 0.0, sizeof(float) * MaxGrid));
        if (nex_temp < 10000)
        {
            file << "Using shared memory kernel" << endl;
            gpuErrchk(cudaOccupancyMaxPotentialBlockSize(
                &numBlocks,
                &blocks,
                (void*)gpu_calc_dens_per_MO_shared,
                0,
                MaxGrid));
        }
        else
        {
            file << "Using static memory kernel" << endl;
            gpuErrchk(cudaOccupancyMaxPotentialBlockSize(
                &numBlocks,
                &blocks,
                (void*)gpu_calc_dens_per_MO_static1,
                0,
                MaxGrid));
        }

        gridSize = (MaxGrid + blocks - 1) / blocks;
        if (debug)
            file << "running " << gridSize << " blocks with " << blocks << " threads" << endl;
        unsigned int nex_offset = 0;

        for (int mo = 0; mo < nmo_temp; mo++)
        {
            gpuErrchk(cudaPeekAtLastError());
            gpuErrchk(cudaDeviceSynchronize());

            if (debug)
                file << "MO " << mo << " starting at coef: " << nex_offset << endl;
            gpuErrchk(cudaMemcpy(gpu_coefficients[0], &coef[nex_offset], sizeof(float) * nex_temp, cudaMemcpyHostToDevice));
            if (nex_temp < 10000)
            {
                gpu_calc_dens_per_MO_shared << <gridSize, blocks, nex_temp * sizeof(float) >> > (
                    gpu_GridRho[0],
                    gpu_Gridx[0],
                    gpu_Gridy[0],
                    gpu_Gridz[0],
                    gpu_PosAtomsx[0],
                    gpu_PosAtomsy[0],
                    gpu_PosAtomsz[0],
                    gpu_types[0],
                    gpu_centers[0],
                    gpu_exponents[0],
                    gpu_coefficients[0],
                    wave.get_MO_occ(mo));
            }
            else
            {
                gpu_calc_dens_per_MO_static2 << <gridSize, blocks >> > (
                    gpu_GridRho[0],
                    gpu_Gridx[0],
                    gpu_Gridy[0],
                    gpu_Gridz[0],
                    gpu_PosAtomsx[0],
                    gpu_PosAtomsy[0],
                    gpu_PosAtomsz[0],
                    gpu_types[0],
                    gpu_centers[0],
                    gpu_exponents[0],
                    gpu_coefficients[0],
                    wave.get_MO_occ(mo));
            }
            nex_offset += nex_temp;
        }
    }
    else
    {
        for (int i = 0; i < nex_temp; i++)
            for (int mo = 0; mo < wave.get_nmo(false); mo++)
                if (wave.get_MO_occ(mo) != 0)
                    coef.push_back(wave.get_MO_coef(mo, i));
        // seems broken and is slower, especially if L1 is sufficient for coefficients with size nex_temp
        gpuErrchk(cudaMalloc((void**)&gpu_coefficients[0], sizeof(float) * nex_temp * nmo_temp));
        gpuErrchk(cudaMemcpy(gpu_coefficients[0], coef.data(), sizeof(float) * nex_temp * nmo_temp, cudaMemcpyHostToDevice));
        vector<float> occ;
        for (int i = 0; i < wave.get_nmo(false); i++)
        {
            occ.push_back(wave.get_MO_occ(i));
            if (occ[occ.size() - 1] == 0)
                occ.pop_back();
        }
        gpuErrchk(cudaMalloc((void**)&gpu_occ[0], sizeof(float) * nmo_temp));
        gpuErrchk(cudaMemcpy(gpu_occ[0], occ.data(), sizeof(float) * occ.size(), cudaMemcpyHostToDevice));
        gpuErrchk(cudaOccupancyMaxPotentialBlockSize(
            &numBlocks,
            &blocks,
            (void*)gpu_calc_dens,
            0,
            MaxGrid));

        gridSize = (MaxGrid + blocks - 1) / blocks;
        gpu_calc_dens << <gridSize, blocks >> > (
            gpu_GridRho[0],
            gpu_Gridx[0],
            gpu_Gridy[0],
            gpu_Gridz[0],
            gpu_PosAtomsx[0],
            gpu_PosAtomsy[0],
            gpu_PosAtomsz[0],
            gpu_types[0],
            gpu_centers[0],
            gpu_exponents[0],
            gpu_coefficients[0],
            gpu_occ[0]);

        occ.clear();
        gpuErrchk(cudaFree(gpu_occ[0]));
    }
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());

    coef.clear();
    ex.clear();

    if (debug)
    {
        // Copy grid from GPU to print:
        //  Dimensions: [c] [p]
        //  p = the number of gridpoint
        //  c = coordinate, which is 0=x, 1=y, 2=z, 3=atomic becke weight, 4=spherical density, 5=wavefunction density, 6=molecular becke weight
        gpuErrchk(cudaMemcpy(total_grid[5].data(), gpu_GridRho[0], sizeof(float) * MaxGrid, cudaMemcpyDeviceToHost));
        ofstream grid("grid.file", ios::out);
        for (int i = 0; i < MaxGrid; i++)
        {
            for (int j = 0; j < 7; j++)
                grid << setw(16) << scientific << setprecision(8) << total_grid[j][i];
            grid << "\n";
        }
        grid.flush();
        grid.close();
    }

    gpuErrchk(cudaFree(gpu_types[0]));
    gpuErrchk(cudaFree(gpu_centers[0]));
    gpuErrchk(cudaFree(gpu_exponents[0]));
    gpuErrchk(cudaFree(gpu_coefficients[0]));
    free(gpu_types);
    free(gpu_centers);
    free(gpu_exponents);
    free(gpu_coefficients);
    free(gpu_occ);

    offset = 0;
    double el_sum_becke = 0.0;
    double el_sum_spherical = 0.0;
    double el_sum_hirshfeld = 0.0;
    for (int i = 0; i < asym_atom_list.size(); i++)
    {
        double charges[3];
        double* result = (double*)malloc(sizeof(double));
        gpuErrchk(cudaMalloc((void**)&(result), sizeof(double) * 3));
        gpu_calc_charges << <1, 1 >> > (
            gpu_GridRho[0],
            gpu_Gridaw[0],
            gpu_Gridmw[0],
            gpu_Grids[0],
            gpu_spherical_density[0][i],
            num_points[i],
            offset,
            cutoff,
            result);
        gpuErrchk(cudaMemcpy(&charges[0], result, sizeof(double) * 3, cudaMemcpyDeviceToHost));
        gpuErrchk(cudaDeviceSynchronize());
        gpuErrchk(cudaPeekAtLastError());
        offset += num_points[i];
        el_sum_becke += charges[0];
        el_sum_spherical += charges[1];
        el_sum_hirshfeld += charges[2];
        for (int j = 0; j < 3; j++)
            atom_els[j][i] = charges[j];
    }

    offset = 0;

    file << "Applying weights..." << endl;
    for (int i = 0; i < asym_atom_list.size(); i++)
    {
        cudaOccupancyMaxPotentialBlockSize(
            &numBlocks,
            &blocks,
            (void*)gpu_apply_weights,
            0,
            num_points[i]);

        gridSize = (MaxGrid + blocks - 1) / blocks;
        // file << i << ": blocks: " << gridSize << " threads: " << blocks << endl;

        gpu_apply_weights << <gridSize, blocks, 0, streams[i] >> > (
            offset,
            gpu_GridRho[0],
            gpu_Grids[0],
            gpu_spherical_density[0][i],
            gpu_Gridaw[0],
            num_points[i]);

        offset += num_points[i];
    }
    gpuErrchk(cudaDeviceSynchronize());
    gpuErrchk(cudaPeekAtLastError());
    file << "After Applying!" << endl;

    file << endl;
    file << " done!" << endl;
    file << "Number of points evaluated: " << MaxGrid;

    file << " with " << fixed << setw(10) << setprecision(6) << el_sum_becke << " electrons in Becke Grid in total." << endl
        << endl;

    file << "Table of Charges in electrons" << endl
        << endl
        << "    Atom       Becke   Spherical Hirshfeld" << endl;

    int counter = 0;
    for (int i = 0; i < cif2wfn_list.size(); i++)
    {
        int a = cif2wfn_list[i];
        file << setw(10) << labels[i]
            << fixed << setw(10) << setprecision(3) << wave.get_atom_charge(a) - atom_els[0][counter]
            << fixed << setw(10) << setprecision(3) << wave.get_atom_charge(a) - atom_els[1][counter]
            << fixed << setw(10) << setprecision(3) << wave.get_atom_charge(a) - atom_els[2][counter];
        if (debug)
            file << " " << setw(4) << wave.get_atom_charge(a) << " " << fixed << setw(10) << setprecision(3) << wave.get_atom_charge(a) - atom_els[0][counter]
            << fixed << setw(10) << setprecision(3) << atom_els[1][counter]
            << fixed << setw(10) << setprecision(3) << atom_els[2][counter];
        counter++;
        file << endl;
    }

    file << "Total number of electrons in the wavefunction: " << el_sum_becke << endl
        << " and Hirshfeld electrons (asym unit): " << el_sum_hirshfeld << endl;

    if (debug)
    {
        file << "Taking time..." << endl;
    }
    end_aspherical = get_time();

    dens.resize(cif2wfn_list.size());
    if (debug)
        file << "resized outer dens" << endl;
    d1.resize(cif2wfn_list.size());
    d2.resize(cif2wfn_list.size());
    d3.resize(cif2wfn_list.size());
    if (debug)
        file << "resized outer d1-3" << endl;
    points = 0;
#pragma omp parallel for
    for (int i = 0; i < asym_atom_list.size(); i++)
    {
        int type;
        for (int j = 0; j < atom_type_list.size(); j++)
            if (atom_type_list[j] == wave.atoms[asym_atom_list[i]].charge)
                type = j;
        vector<float> temp_dens;
        temp_dens.resize(num_points[i]);
        // file << "At atom: " << i;
        gpuErrchk(cudaMemcpy(temp_dens.data(), gpu_spherical_density[0][i], sizeof(float) * num_points[i], cudaMemcpyDeviceToHost));
        for (int p = 0; p < 0 + num_points[i]; p++)
        {
            if (abs(temp_dens[p]) > _cutoff)
            {
                dens[i].push_back(temp_dens[p]);
                d1[i].push_back(Prototype_grids[type].get_gridx(p));
                d2[i].push_back(Prototype_grids[type].get_gridy(p));
                d3[i].push_back(Prototype_grids[type].get_gridz(p));
            }
        }
        // file << " dens size: " << dens[i].size() << " num_points[i]: " << num_points[i] << endl;
    }
    for (int i = 0; i < asym_atom_list.size(); i++)
        points += dens[i].size();

    // gpuErrchk(cudaFree(gpu_Grids[0]));
    // gpuErrchk(cudaFree(gpu_PosAtomsx[0]));
    // gpuErrchk(cudaFree(gpu_PosAtomsy[0]));
    // gpuErrchk(cudaFree(gpu_PosAtomsz[0]));
    // gpuErrchk(cudaFree(gpu_GridRho[0]));
    // gpuErrchk(cudaFree(gpu_Gridx[0]));
    // gpuErrchk(cudaFree(gpu_Gridy[0]));
    // gpuErrchk(cudaFree(gpu_Gridz[0]));
    // gpuErrchk(cudaFree(gpu_Gridaw[0]));
    // gpuErrchk(cudaFree(gpu_Gridmw[0]));
    // gpuErrchk(cudaFree(gpu_asym_atom_list[0]));
    // gpuErrchk(cudaFree(gpu_atom_type_list[0]));
    // gpuErrchk(cudaFree(gpu_numpoints[0]));
    // gpuErrchk(cudaFree(gpu_atom_z[0]));
    // for (int i = 0; i < asym_atom_list.size(); i++) {
    //	gpuErrchk(cudaFree(gpu_spherical_density[0][i]));
    //	gpuErrchk(cudaFree(gpu_atomgrid_x[0][i]));
    //	gpuErrchk(cudaFree(gpu_atomgrid_y[0][i]));
    //	gpuErrchk(cudaFree(gpu_atomgrid_z[0][i]));
    // }
    free(gpu_Grids);
    free(gpu_PosAtomsx);
    free(gpu_PosAtomsy);
    free(gpu_PosAtomsz);
    free(gpu_GridRho);
    free(gpu_Gridx);
    free(gpu_Gridy);
    free(gpu_Gridz);
    free(gpu_Gridaw);
    free(gpu_Gridmw);
    free(gpu_asym_atom_list);
    free(gpu_atom_type_list);
    free(gpu_numpoints);
    free(gpu_atom_z);
    free(gpu_radial_density);
    free(gpu_radial_dist);
    free(gpu_spherical_density);
    free(gpu_atomgrid_x);
    free(gpu_atomgrid_y);
    free(gpu_atomgrid_z);
    free(gpu_atomgrid_w);
    cudaDeviceReset();
    file << "CUDA device resetted!" << endl;
    return points;
}