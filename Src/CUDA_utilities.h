#pragma once

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

