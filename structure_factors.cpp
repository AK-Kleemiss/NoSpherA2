/*
 * structrue_factors.cpp
 *
 *  Created on: May 27, 2019
 *      Author: florian
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <list>
#include <numeric>
#include <iomanip>
#include <string>
#include <fstream>
#include <regex>
#include <omp.h>
#ifdef _WIN32
#include <windows.h>
#include <io.h>
#include <omp.h>
#else
#include <unistd.h>
#include <sys/wait.h>
#include <sys/time.h>
#endif
#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <vector>
#include <complex>
#include <algorithm>

#include "cell.h"
#include "convenience.h"
#include "structure_factors.h"
#include "wfn_class.h"
#include "spherical_density.h"
using namespace std;
# define M_PI           3.14159265358979323846  /* pi */

#include "AtomGrid.h"
#ifdef PEOJECT_NAME
#define FLO_CUDA
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#endif

bool merge_tscs(
	const string& mode,
	const vector<string>& files,
	const bool debug
) {
	//Currently only for Mode "pure merge"
	if (files.size() == 0)
		return false;
	vector <string> labels;
	vector <vector <int>> indices; //[i] is h,k or l, [j] is the number of the reflection
	vector <vector <complex<double>>> form_fact;   //[i] (len(labels)) scatterer, [j](len(indices[0])) reflection correpsonding to indices

	vector<string> header;

	size_t offset = 0;
	indices.resize(3);
	for (int f = 0; f < files.size(); f++) {
		cout << "Reading file number: " << f + 1 << endl;
		ifstream inf(files[f].c_str(), ios::in);
		string line;
		bool data = false;
		vector<bool> is_a_new_scatterer;
		int nr_scatterers = 0;
		while (!data) {
			getline(inf, line);
			if (line.find("SCATTERERS:") != string::npos) {
				string temp_labels = line.substr(12, line.size()) + " ";
				const string delimiter = " ";
				size_t pos = 0;
				string new_label;
				while ((pos = temp_labels.find(delimiter)) != std::string::npos) {
					nr_scatterers++;
					new_label = temp_labels.substr(0, pos);
					bool is_new = true;
					for (int i = 0; i < labels.size(); i++)
						if (labels[i] == new_label)	is_new = false;
					is_a_new_scatterer.push_back(is_new);
					if(is_new) labels.push_back(new_label);
					temp_labels.erase(0, pos + delimiter.length());
				}
				cout << "Read " << nr_scatterers << " atoms, " << sum_of_bools(is_a_new_scatterer) << " are new." << endl;
				form_fact.resize(labels.size());
				form_fact.resize(labels.size());
			}
			else if (line.find("DATA:") != string::npos)
				data = true;
			else if (f == 0) header.push_back(line);
		}
		cout << "Reading Data Block..." << endl;
		while (!inf.eof()) {
			getline(inf, line);
			vector<string> digest = split_string(line, " ");
			if (digest.size() == 1 && digest[0] == "")
				continue;
			if (f==0)
				for (int i = 0; i < 3; i++)
					indices[i].push_back(stoi(digest[i]));
#pragma omp parallel for
			for (int i = 0; i < nr_scatterers; i++) {
				if (!is_a_new_scatterer[i]) continue;
				size_t nr_in_new_scats = 0;
				for (int p = 0; p < i; p++)
					if(is_a_new_scatterer[p]) nr_in_new_scats++;
				vector<string> re_im = split_string(digest[i + 3], ",");
				form_fact[offset+nr_in_new_scats].push_back(complex<double>(stod(re_im[0]),stod(re_im[1])));
			}
		}
		if (form_fact[form_fact.size() - 1].size() != form_fact[0].size()) {
			cout << "ERROR! Incorrect number of reflections among files! Aborting!" << endl;
			return false;
		}
		cout << "Data for " << form_fact[0].size() << " indices read." << endl;
		offset = labels.size();
	}
	cout << "Writing combined file..." << endl;
	ofstream tsc_file("combined.tsc", ios::out);

	for (size_t h = 0; h < header.size(); h++)
		tsc_file << header[h] << endl;
	tsc_file << "    PARTS: " << files.size() << endl;
	tsc_file << "SCATTERERS:";
	for (int i = 0; i < labels.size(); i++)
		tsc_file << " " << labels[i];
	tsc_file << endl << "DATA:" << endl;

	for (int r = 0; r < indices[0].size(); r++) {
		for (int h = 0; h < 3; h++)
			tsc_file << indices[h][r] << " ";
		for (int i = 0; i < labels.size(); i++)
			tsc_file << scientific << setprecision(8) << real(form_fact[i][r]) << ","
			<< scientific << setprecision(8) << imag(form_fact[i][r]) << " ";
		tsc_file << endl;
	}
	tsc_file.flush();
	tsc_file.close();
	cout << "Done!" << endl;

	return true;
}

#ifdef FLO_CUDA
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
	const int &l) {
	if (l == 2)
		vector0 = 1, vector1 = 0, vector2 = 0;
	else if (l == 3)
		vector0 = 0, vector1 = 1, vector2 = 0;
	else if (l == 4)
		vector0 = 0, vector1 = 0, vector2 = 1;
}

__device__ float gpu_ptype(
	const float& vector0, const float& vector1, const float& vector2,
	const int& l) {
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
	const int &l) {
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
	const int& l) {
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
	const int &l) {
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
	const int& l) {
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
	const int &l) {
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
	const int& l) {
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
	int &vector0, int &vector1, int &vector2,
	const int &l) {
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
	const int& l) {
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
	const int &index,
	int &vector0, int &vector1, int &vector2) {
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
	float& vector0, float& vector1, float& vector2) {
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

__device__ int gpu_type_to_order(const int index) {
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
		else{
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
	float ex, phi=0.0;
	int i, iat;

	if(threadIdx.x == 0)
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
	float* d0 = (float*) malloc(sizeof(float) * max_center);
	float* d1 = (float*) malloc(sizeof(float) * max_center);
	float* d2 = (float*) malloc(sizeof(float) * max_center);
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
	float *d0 = (float*)malloc(sizeof(float) * max_center);
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

void printCUDA(const cudaDeviceProp* prop, const int nDevices, ofstream& file)
{
	file << endl << "Number of CUDA-capable devices : " << nDevices << endl;
	/*table with the properties of the devices*/

	for (int dev = 0; dev < nDevices; dev++) {
		file << "Device Number: " << dev << endl;
		file << "- Device name: " << (prop[dev]).name << endl;
		file << "- Max global memory   = " << prop[dev].totalGlobalMem / 1048576.0f << " MBytes" << endl;
		file << "- Max constant memory = " << prop[dev].totalConstMem / 1048576.0f << " MBytes" << endl;
		file << "- Capability : " << prop[dev].major << "." << prop[dev].minor << endl;
		file << "- Total global memory: "<< (float)prop[dev].totalGlobalMem / 1048576.0f << " MByte" << endl;
		file << "_____________________________________________________" << endl;
	}
}

#else

double compute_dens(
	WFN& wave,
	double* PosGrid,			// [3] array with current position on the grid
	double** mo_coef,
	int atom = -1
)
{
	const int nmo = wave.get_nmo(true);
	double* phi = new double[nmo];
	double Rho = 0.0;

	for (int i = 0; i < nmo; i++)
		phi[i] = 0.0;

	int iat;
	int l[3];
	double ex;
	int mo;
	double temp;

	vector<vector<float> > d;
	// x, y, z and dsqd
	d.resize(4);
	for (int i = 0; i < 4; i++)
		d[i].resize(wave.get_ncen());

	for (iat = 0; iat < wave.get_ncen(); iat++) {
		d[0][iat] = PosGrid[0] - wave.atoms[iat].x;
		d[1][iat] = PosGrid[1] - wave.atoms[iat].y;
		d[2][iat] = PosGrid[2] - wave.atoms[iat].z;
		d[3][iat] = d[0][iat] * d[0][iat] + d[1][iat] * d[1][iat] + d[2][iat] * d[2][iat];
	}

	for (int j = 0; j < wave.get_nex(); j++) {
		iat = wave.get_center(j) - 1;
//		if (iat != atom) continue;
		type2vector(wave.get_type(j), l);
		temp = -wave.get_exponent(j) * d[3][iat];
		if (temp < -46.0517) //corresponds to cutoff of ex ~< 1E-20
			continue;
		ex = exp(temp);
		for (int k = 0; k < 3; k++) {
			if (l[k] == 0)		continue;
			else if (l[k] == 1)	ex *= d[k][iat];
			else if (l[k] == 2)	ex *= d[k][iat] * d[k][iat];
			else if (l[k] == 3)	ex *= pow(d[k][iat], 3);
			else if (l[k] == 4)	ex *= pow(d[k][iat], 4);
			else if (l[k] == 5)	ex *= pow(d[k][iat], 5);
		}
		for (mo = 0; mo < nmo; mo++)
			phi[mo] += mo_coef[mo][j] * ex;      //build MO values at this point
	}

	d.clear();

	for (mo = 0; mo < nmo; mo++)
		Rho += wave.get_MO_occ(mo) * pow(phi[mo], 2);
	delete[](phi);
	
	return Rho;
}

double linear_interpolate_spherical_density(
	vector <double>& radial_dens,
	vector <double>& spherical_dist,
	const double dist,
	const double lincr,
	const double start
)
{
	double result = 0;
	if (dist > spherical_dist[spherical_dist.size() - 1])
		return 0;
	else if (dist < spherical_dist[0])
		return radial_dens[0];
	int nr = floor(log(dist / start) / lincr);
	result = radial_dens[nr] + (radial_dens[nr + 1] - radial_dens[nr]) / (spherical_dist[nr] - spherical_dist[nr - 1]) * (dist - spherical_dist[nr - 1]);
	if (result < 1E-10) result = 0;
	return result;
}
#endif

bool thakkar_sfac(
	const int atom_number,
	const double radial_acc,
	const int angular_level
) {
	Thakkar atom(atom_number);
	double alpha_max = atom.get_max_alpha();
	vector<double> alpha_min = atom.get_min_alpha();
	const int max_l = atom.get_max_l();
	int lebedev_high=146, lebedev_low=50;
	if (angular_level == 0) {
		lebedev_high = (max_l < 3) ? 0 : 10;
		lebedev_low = (max_l < 3) ? 0 : 10;
	}
	else if (angular_level == 1) {
		lebedev_high = (max_l < 3) ? 110 : 146;
		lebedev_low = (max_l < 3) ? 38 : 50;
	}
	else if (angular_level == 2) {
		lebedev_high = (max_l < 3) ? 230 : 266;
		lebedev_low = (max_l < 3) ? 110 : 146;
	}
	else if (angular_level == 3) {
		lebedev_high = (max_l < 3) ? 350 : 590;
		lebedev_low = (max_l < 3) ? 266 : 350;
	}
	else if (angular_level > 3) {
		lebedev_high = (max_l < 3) ? 590 : 5810;
		lebedev_low = (max_l < 3) ? 350 : 4802;
	}
	AtomGrid Prototype_grid(radial_acc,
		lebedev_low,
		lebedev_high,
		atom_number,
		alpha_max,
		max_l,
		alpha_min.data()
	);

	int num_points = Prototype_grid.get_num_grid_points();

	cout << "Number of points for atom: " << num_points << endl;

	double **grid = new double*[4];
	for (int n = 0; n < 4; n++)
		grid[n] = new double[num_points];

	double pos = 0.0;

	Prototype_grid.get_grid_omp(1,
		0,
		&pos,
		&pos,
		&pos,
		&atom_number,
		grid[0],
		grid[1],
		grid[2],
		grid[3]);

	return true;
}

bool calculate_structure_factors_HF(
	string& hkl_filename,
	string& cif,
	string& asym_cif,
	string& symm,
	WFN& wave,
	bool debug,
	int accuracy,
	ofstream& file,
	vector <int>& input_groups,
	vector < vector <double> >& twin_law,
	int cpus,
	bool electron_diffraction,
	int pbc
)
{
#ifdef FLO_CUDA
	if (pbc != 0) {
		file << "PERIODIC CALCULATIONS NOT IMPLEMENTED WITH CUDA YET!" << endl;
		exit(false);
	}
#endif
	if (cpus != -1) {
		omp_set_num_threads(cpus);
		omp_set_dynamic(0);
	}
	vector<double> x, y, z;
	vector<int> atom_z;
	x.resize(wave.get_ncen() * pow(pbc * 2 + 1, 3)), y.resize(wave.get_ncen() * pow(pbc * 2 + 1, 3)), z.resize(wave.get_ncen() * pow(pbc * 2 + 1, 3)), atom_z.resize(wave.get_ncen() * pow(pbc * 2 + 1, 3));
	double* alpha_max = new double[wave.get_ncen()];
	int* max_l = new int[wave.get_ncen()];
	int max_l_overall = 0;

#ifdef _WIN64
	time_t start = time(NULL);
	time_t end_becke, end_prototypes, end_spherical, end_prune, end_aspherical;
#else
	struct timeval t1, t2;

	gettimeofday(&t1, 0);
#endif

	file << "Reading: " << hkl_filename;
	file.flush();

	vector< vector <int> > hkl;
	string line;
	hkl.resize(3);
	if (!exists(hkl_filename)) {
		file << "HKL file does not exists!" << endl;
		return false;
	}
	ifstream hkl_input(hkl_filename.c_str(), ios::in);
	hkl_input.seekg(0, hkl_input.beg);
	regex r{ R"([abcdefghijklmnopqrstuvwxyz\(\)ABCDEFGHIJKLMNOPQRSTUVW])" };
	while (!hkl_input.eof()) {
		getline(hkl_input, line);
		if (hkl_input.eof())
			break;
		if (line.size() < 2)
			continue;
		cmatch result;
		if (regex_search(line.c_str(), result, r))
			continue;
		//if (debug) file << "hkl: ";
		for (int i = 0; i < 3; i++) {
			string temp = line.substr(4 * size_t(i) + 1, 3);
			temp.erase(remove_if(temp.begin(), temp.end(), ::isspace), temp.end());
			hkl[i].push_back(stoi(temp));
			//if (debug) file << setw(4) << temp;
		}
		//if (debug) file << endl;
		if (hkl[0][hkl[0].size() - 1] == 0 && hkl[1][hkl[0].size() - 1] == 0 && hkl[2][hkl[0].size() - 1] == 0) {
			if (debug) file << "popping back 0 0 0" << endl;
			for (int i = 0; i < 3; i++)
				hkl[i].pop_back();
		}
	}
	hkl_input.close();
	file << "                   ...done!" << endl;
	int refl_size = hkl[0].size();
	if (debug)
		file << "Number of reflections before twin: " << hkl[0].size() << endl;
	for (int len = 0; len < refl_size; len++)
		for (int i = 0; i < twin_law.size(); i++)
			for (int h = 0; h < 3; h++)
				hkl[h].push_back(int(twin_law[i][0 + 3 * h] * hkl[0][len] + twin_law[i][1 + 3 * h] * hkl[1][len] + twin_law[i][2 + 3 * h] * hkl[2][len]));
	if (debug)
		file << "Number of reflections after twin: " << hkl[0].size() << endl;

	file << "Remove duplicate reflections...";
	file.flush();
	for (int i = 0; i < hkl[0].size(); i++)
		for (int j = i + 1; j < hkl[0].size(); j++)
			if (hkl[0][i] == hkl[0][j] && hkl[1][i] == hkl[1][j] && hkl[2][i] == hkl[2][j])
				for (int x = 0; x < 3; x++)
					hkl[x].erase(hkl[x].begin() + j);

	file << "                   ...done!\nNr of reflections to be used: " << hkl[0].size() << endl;

	if (debug)
		file << "starting to read cif!" << endl;
	if (!exists(cif)) {
		file << "CIF does not exists!" << endl;
		return false;
	}
	if (!exists(asym_cif)) {
		file << "Asym CIF does not exists!" << endl
			<< asym_cif << endl;
		return false;
	}

	cell unit_cell(cif, file, debug);

	if (debug) {
		file << "RCM done, now labels and asym atoms!" << endl;
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j)
				file << setw(10) << fixed << unit_cell.get_rcm(i, j) / 2 / M_PI / 0.529177249 << ' ';
			file << endl;
		}
		file << "CM in 2*PI bohr:" << endl;
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j)
				file << setw(10) << fixed << unit_cell.get_cm(i, j) << ' ';
			file << endl;
		}
	}

	ifstream cif_input(cif.c_str(), std::ios::in);
	ifstream asym_cif_input(asym_cif.c_str(), std::ios::in);
	cif_input.clear();
	cif_input.seekg(0, cif_input.beg);
	vector <string> labels;
	vector <int> atom_type_list;
	vector <int> asym_atom_to_type_list;
	int count_fields = 0;
	int group_field = 0;
	int position_field[3] = { 0,0,0 };
	int label_field = 1000;
	vector <int> asym_atom_list;
	//vector <int> all_atom_list;
	vector < bool > is_asym;
	vector < vector <double > > positions;
	if (wave.get_ncen() == 0) {
		file << "No Atoms in the wavefunction, this will not work!! ABORTING!!" << endl;
		return false;
	}
	positions.resize(wave.get_ncen());
	is_asym.resize(wave.get_ncen());
	for (int i = 0; i < wave.get_ncen(); i++) {
		is_asym[i] = false;
		positions[i].resize(3);
	}
	bool atoms_read = false;
	while (!asym_cif_input.eof() && !atoms_read) {
		getline(asym_cif_input, line);
		if(debug) file << "line: "<< line << endl;
		if (line.find("loop_") != string::npos) {
			//if(debug) file << "found loop!" << endl;
			if(debug) file << trim(line).find("_") << endl;
			getline(asym_cif_input, line);
			while (trim(line).find("_") == 0) {
				if (debug) file << "line in loop field definition: " << trim(line) << endl;
				if (debug) file << trim(line).find("_") << endl;
				if (line.find("label") != string::npos)
					label_field = count_fields;
				else if (line.find("fract_x") != string::npos)
					position_field[0] = count_fields;
				else if (line.find("fract_y") != string::npos)
					position_field[1] = count_fields;
				else if (line.find("fract_z") != string::npos)
					position_field[2] = count_fields;
				else if (label_field == 1000) {
					if (debug) file << "I don't think this is the atom block.. moving on!" << endl;
					break;
				}
				getline(asym_cif_input, line);
				count_fields++;
			}
			if (label_field == 1000) continue;
			while (trim(line).find("_") > 0 && line.length() > 3) {
				//if(debug) file << "Reading atom!"<< endl;
				atoms_read = true;
				stringstream s(line);
				vector <string> fields;
				fields.resize(count_fields);
				for (int i = 0; i < count_fields; i++)
					s >> fields[i];
				if (debug) file << "label: " << fields[label_field] << " frac_position: " << stod(fields[position_field[0]]) << " " << stod(fields[position_field[1]]) << " " << stod(fields[position_field[2]]) << flush;
				positions[labels.size()] = unit_cell.get_coords_cartesian(stod(fields[position_field[0]]), stod(fields[position_field[1]]), stod(fields[position_field[2]]));
				bool found_this_one = false;
				if (debug) file << " cart. position: " << positions[labels.size()][0] << " " << positions[labels.size()][1] << " " << positions[labels.size()][2] << endl;
				for (int i = 0; i < wave.get_ncen(); i++) {
					if (is_similar(positions[labels.size()][0], wave.atoms[i].x, -0.5)
						&& is_similar(positions[labels.size()][1], wave.atoms[i].y, -0.5)
						&& is_similar(positions[labels.size()][2], wave.atoms[i].z, -0.5)) {
						if (debug) file << "WFN position: " << wave.atoms[i].x << " " << wave.atoms[i].y << " " << wave.atoms[i].z
							<< " Found an atom: " << fields[label_field] << " Corresponding to atom charge " << wave.atoms[i].charge << endl;
						wave.atoms[i].label = fields[label_field];
						//all_atom_list.push_back(i);
						found_this_one = true;
						break;
					}
				}
				if (!found_this_one && debug)
					file << "I DID NOT FIND THIS ATOM IN THE CIF?! WTF?!" << endl << positions[labels.size()][0] << " " << positions[labels.size()][1] << " " << positions[labels.size()][2] << endl;
				labels.push_back(fields[label_field]);
				getline(asym_cif_input, line);
			}
		}
	}
	//if (labels.size() != wave.get_ncen()) {
	//	file << "Number of atoms in labels: " << labels.size() << " and number of atoms in Wavefunction: " << wave.get_ncen() << "!" << endl << "This is not good, i will stop here!" << endl;
	//	return false;
	//}
	atoms_read = false;
	label_field = 1000;
	count_fields = 0;
	while (!cif_input.eof() && !atoms_read) {
		getline(cif_input, line);
		if(debug) file << "line: "<< line << endl;
		if (line.find("loop_") != string::npos) {
			//if(debug) file << "found loop!" << endl;
			getline(cif_input, line);
			if (debug) file << "line in loop field definition: " << trim(line) << endl;
			while (trim(line).find("_") == 0) {
				if (debug) file << "line in loop field definition: " << trim(line) << endl;
				if (line.find("label") != string::npos)
					label_field = count_fields;
				else if (line.find("disorder_group") != string::npos)
					group_field = count_fields;
				else if (label_field == 1000) {
					if (debug) file << "I don't think this is the atom block.. moving on!" << endl;
					break;
				}
				getline(cif_input, line);
				count_fields++;
			}
			while (trim(line).find("_") > 0 && line.length() > 3) {
				atoms_read = true;
				stringstream s(line);
				vector <string> fields;
				fields.resize(count_fields);
				for (int i = 0; i < count_fields; i++)
					s >> fields[i];
				fields[label_field].erase(remove_if(fields[label_field].begin(), fields[label_field].end(), ::isspace), fields[label_field].end());
				for (int atom = 0; atom < wave.get_ncen(); atom++) {
					if (debug) file << "Comparing atoms: " << fields[label_field] << " / " << labels[atom] << endl;
					if (fields[label_field] == labels[atom]) {
						int nr = -1;
						for (int i = 0; i < wave.get_ncen(); i++) {
							if (is_similar(positions[atom][0], wave.atoms[i].x, -0.5)
								&& is_similar(positions[atom][1], wave.atoms[i].y, -0.5)
								&& is_similar(positions[atom][2], wave.atoms[i].z, -0.5)) {
								if (debug) {
									file << "Found an asymmetric atom: " << fields[label_field] << " Corresponding to atom charge " << wave.atoms[i].charge << endl
										<< "checking disorder group: " << fields[group_field] << " vs. ";
									for (int g = 0; g < input_groups.size(); g++)
										file << input_groups[g] << "+";
									file << endl;
								}
								bool yep = false;
								for (int g = 0; g < input_groups.size(); g++)
									if (stod(fields[group_field]) == input_groups[g])
										yep = true;
								if (!yep && input_groups.size() > 0) continue;
								nr = i;
								asym_atom_list.push_back(i);
								break;
							}
						}
						is_asym[atom] = true;
						bool already_there = false;
						for (int i = 0; i < atom_type_list.size(); i++)
							if (atom_type_list[i] == wave.atoms[nr].charge) {
								already_there = true;
								asym_atom_to_type_list.push_back(i);
							}
						if (already_there == false) {
							asym_atom_to_type_list.push_back(atom_type_list.size());
							atom_type_list.push_back(wave.atoms[nr].charge);
						}
						break;
					}
				}
				getline(cif_input, line);
			}
		}
	}

	if (debug) {
		file << "There are " << atom_type_list.size() << " types of atoms" << endl;
		for (int i = 0; i < atom_type_list.size(); i++)
			file << setw(4) << atom_type_list[i];
		file << endl << "asym_atoms_to_type_list: " << endl;
		for (int i = 0; i < asym_atom_to_type_list.size(); i++)
			file << setw(4) << asym_atom_to_type_list[i];
		file << endl;
		file << "Mapping of asym atoms:" << endl;
		for (int i = 0; i < wave.get_ncen(); i++)
			file << setw(4) << wave.atoms[i].charge;
		file << endl;
		for (int i = 0; i < wave.get_ncen(); i++)
			file << setw(4) << is_asym[i];
		file << endl;
	}

	bool symm_read = false;
	//Still need to read the sym matrices
	if (symm == "") {
		if (debug) file << "No Symmetry file specified, tyring to read from the CIF!" << endl;
		symm_read = true;
	}
	else if (!exists(symm) && symm_read == true)
		return false;

	vector < vector < vector <int> > > sym;
	sym.resize(3);
	for (int i = 0; i < 3; i++)
		sym[i].resize(3);
	if (!symm_read) {
		ifstream symm_input(symm.c_str(), ios::in);
		string liny;
		int temp_int;
		while (!symm_input.eof()) {
			getline(symm_input, liny);
			stringstream streamy(liny);
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++) {
					streamy >> temp_int;
					sym[i][j].push_back(temp_int);
				}
		}
	}
	else {
		sym = unit_cell.get_sym();
	}

	if (debug) {
		file << "Read " << sym[0][0].size() << " symmetry elements!" << endl;
		for (int i = 0; i < sym[0][0].size(); i++) {
			for (int x = 0; x < 3; x++) {
				for (int y = 0; y < 3; y++)
					file << setw(3) << sym[y][x][i];
				file << endl;
			}
			file << endl;
		}
	}
	else
		file << "Number of symmetry operations: " << sym[0][0].size() << endl;

	cif_input.close();

	if (debug)
		file << "There are " << atom_type_list.size() << " Types of atoms and " << asym_atom_to_type_list.size() << " atoms in total" << endl;

	if (asym_atom_list.size() == 0) {
		file << "0 asym atoms is imposible! something is wrong with reading the CIF!" << endl;
		return false;
	}

	if (debug)
		file << "made it post CIF, now make grids!" << endl;

	double*** grid = new double** [asym_atom_list.size()];
	vector<int> num_points;
	num_points.resize(asym_atom_list.size());
	for (int i = 0; i < asym_atom_list.size(); i++)
		grid[i] = new double* [6];
	// GRID COORDINATES for [a][c][p] a = atom [0,ncen],
	// c = coordinate [0=x, 1=y, 2=z, 3=atomic becke weight, 4= molecular becke weight, 5=total spherical density],
	// p = point in this grid

#pragma omp parallel for
	for (int i = 0; i < wave.get_ncen(); i++) {
		atom_z[i] = wave.atoms[i].charge;
		x[i] = wave.atoms[i].x;
		y[i] = wave.atoms[i].y;
		z[i] = wave.atoms[i].z;
		//if(debug)
		//    file << "xyz= 000 position: " << x[i] << " " << y[i] << " " << z[i] << " Charge: " << atom_z[i] << endl;
		if (pbc != 0) {
			int j = 0;
			for (int pbc_x = -pbc; pbc_x < pbc + 1; pbc_x++)
				for (int pbc_y = -pbc; pbc_y < pbc + 1; pbc_y++)
					for (int pbc_z = -pbc; pbc_z < pbc + 1; pbc_z++) {
						if (pbc_x == 0 && pbc_y == 0 && pbc_z == 0)
							continue;
						else {
							j++;
							atom_z[i + j * wave.get_ncen()] = wave.atoms[i].charge;
							x[i + j * wave.get_ncen()] = wave.atoms[i].x + pbc_x * unit_cell.get_cm(0, 0) + pbc_y * unit_cell.get_cm(0, 1) + pbc_z * unit_cell.get_cm(0, 2);
							y[i + j * wave.get_ncen()] = wave.atoms[i].y + pbc_x * unit_cell.get_cm(1, 0) + pbc_y * unit_cell.get_cm(1, 1) + pbc_z * unit_cell.get_cm(1, 2);
							z[i + j * wave.get_ncen()] = wave.atoms[i].z + pbc_x * unit_cell.get_cm(2, 0) + pbc_y * unit_cell.get_cm(2, 1) + pbc_z * unit_cell.get_cm(2, 2);
							if (debug)
								file << "xyz= " << pbc_x << pbc_y << pbc_z << " j = " << j << " position: " << x[i + j * wave.get_ncen()] << " " << y[i + j * wave.get_ncen()] << " " << z[i + j * wave.get_ncen()] << " Charge: " << atom_z[i + j * wave.get_ncen()] << endl;
						}
					}
		}
		alpha_max[i] = 0.0;
		max_l[i] = 0;
		for (int b = 0; b < wave.get_nex(); b++) {
			if (wave.get_center(b) != i + 1)
				continue;
			if (wave.get_exponent(b) > alpha_max[i])
				alpha_max[i] = wave.get_exponent(b);
			if (wave.get_type(b) > max_l[i]) {
				int l = wave.get_type(b);
				if (l == 1)
					l = 1;
				else if (l >= 2 && l <= 4)
					l = 2;
				else if (l >= 5 && l <= 10)
					l = 3;
				else if (l >= 11 && l <= 20)
					l = 4;
				else if (l >= 21 && l <= 35)
					l = 5;
				max_l[i] = l;
				if (l > max_l_overall)
					max_l_overall = l;
			}
		}
	}

	if (debug)
		file << "Atoms are there!" << endl;

	double** alpha_min = new double* [wave.get_ncen()];
	for (int i = 0; i < wave.get_ncen(); i++)
		alpha_min[i] = new double[max_l_overall];

#pragma omp parallel for
	for (int i = 0; i < wave.get_ncen(); i++) {
		for (int b = 0; b < max_l_overall; b++)
			alpha_min[i][b] = 100000000.0;
	}

#pragma omp parallel for
	for (int i = 0; i < wave.get_ncen(); i++) {
		for (int b = 0; b < wave.get_nex(); b++) {
			if (wave.get_center(b) != i + 1)
				continue;
			int l = wave.get_type(b);
			if (l == 1)
				l = 1;
			else if (l >= 2 && l <= 4)
				l = 2;
			else if (l >= 5 && l <= 10)
				l = 3;
			else if (l >= 11 && l <= 20)
				l = 4;
			else if (l >= 21 && l <= 35)
				l = 5;
			else if (l >= 36 && l <= 56)
				l = 6;
			if (wave.get_exponent(b) < alpha_min[i][l - 1])
				alpha_min[i][l - 1] = wave.get_exponent(b);
		}
		if (debug) {
			file << "alpha_min: ";
			for (int b = 0; b < max_l_overall; b++)
				file << setw(14) << scientific << alpha_min[i][b];
			file << endl;
		}
	}

	if (debug)
		file << "alpha_min is there!" << endl
		<< "Nr of asym atoms: " << asym_atom_list.size() << " Number of atoms in wfn: " << wave.get_ncen() << endl;
	else
		file << "There are:\n" << setw(4) << wave.get_ncen() << " atoms read from the wavefunction, of which \n"
		//<< setw(4) << all_atom_list.size() << " will be used for grid setup and\n"
		<< setw(4) << asym_atom_list.size() << " are identified as asymmetric unit atoms!" << endl;
#ifdef FLO_CUDA
	int nDevices;/*Number of devices available (running time)*/
	cudaGetDeviceCount(&nDevices);/*Get the number of devices*/
	int dev = 0;
	cudaDeviceProp* prop = NULL;
	cudaDeviceProp deviceProp;
	prop = (cudaDeviceProp*)malloc(sizeof(cudaDeviceProp) * nDevices);
	for (int devl = 0; devl < nDevices; devl++) { // Make CUDA information available in prop
		cudaGetDeviceProperties(&(deviceProp), devl);
		prop[devl] = deviceProp;
	}
	cudaSetDevice(dev);
	printCUDA(prop, nDevices, file);
	gpuErrchk(cudaSetDeviceFlags(cudaDeviceScheduleBlockingSync));
#endif

	// Total grid as a sum of all atomic grids.
	// Dimensions: [c] [p]
	// p = the number of gridpoint
	// c = coordinate, which is 0=x, 1=y, 2=z, 3=atomic becke weight, 4=spherical density, 5=wavefunction density, 6=molecular becke weight
#ifdef FLO_CUDA
	vector < vector < float > > total_grid;
#else
	vector < vector < double > > total_grid;
#endif
	// density of spherical atom at each
	// Dimensions: [a] [d]
	// a = atom number in atom type list for which the weight is calcualted
	// d = distance to look at obtained from point_to_distance_map
	vector < vector < double > > spherical_density;

	if (debug)
		file << "Making Becke Grid for" << endl;
	else
		file << endl << "Making Becke Grids..." << flush;

	total_grid.resize(7);

	//Make Prototype grids with only single atom weights for all elements
	vector <AtomGrid> Prototype_grids;
	if (debug) file << "max_l_overall: " << max_l_overall << endl;
	const int lebedev_table[33] = { 6,   14,   26,   38,   50,   74,   86,  110,
								  146,  170,  194,  230,  266,  302,  350,  434,
								  590,  770,  974, 1202, 1454, 1730, 2030, 2354,
								 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810 };
	for (int i = 0; i < atom_type_list.size(); i++) {
		if (debug) file << "Atom " << i << endl;
		double alpha_max_temp;
		double max_l_temp;
		double* alpha_min_temp = new double[max_l_overall];
		for (int j = 0; j < wave.get_ncen(); j++) {
			if (wave.atoms[j].charge == atom_type_list[i]) {
				if (debug) {
					file << alpha_max[j] << " " << max_l[j] - 1 << " ";
					for (int l = 0; l < max_l_overall; l++)
						file << alpha_min[j][l] << " ";
					file << endl;
				}
				alpha_max_temp = alpha_max[j];
				max_l_temp = max_l[j] - 1;
				for (int l = 0; l <= max_l_temp; l++)
					alpha_min_temp[l] = alpha_min[j][l];
				break;
			}
		}

		if (debug) {
			file << "max_l: " << defaultfloat << max_l_temp << " alpha_max: " << scientific << alpha_max_temp << " alpha_min: ";
			for (int l = 0; l <= max_l_temp; l++)
				file << setw(14) << scientific << alpha_min_temp[l];
			file << " accuracy: " << accuracy << endl;
		}
		int lebedev_high, lebedev_low;
		double radial_acc;
		if (accuracy == 0) {
			if (atom_type_list[i] != 1) {
				lebedev_high = (max_l_temp < 3) ? lebedev_table[0] : lebedev_table[1];
				lebedev_low = (max_l_temp < 3) ? lebedev_table[0] : lebedev_table[1];
				radial_acc = 1e-2;
			}
			else {
				lebedev_high = (max_l_temp < 3) ? lebedev_table[0] : lebedev_table[1];
				lebedev_low = (max_l_temp < 3) ? lebedev_table[0] : lebedev_table[1];
				radial_acc = 1e-3;
			}
		}
		else if (accuracy == 1) {
			if (atom_type_list[i] != 1) {
				lebedev_high = (max_l_temp < 3) ? lebedev_table[7] : lebedev_table[8];
				lebedev_low = (max_l_temp < 3) ? lebedev_table[3] : lebedev_table[4];
				radial_acc = 1e-4;
			}
			else {
				lebedev_high = (max_l_temp < 3) ? lebedev_table[6] : lebedev_table[7];
				lebedev_low = (max_l_temp < 3) ? lebedev_table[2] : lebedev_table[3];
				radial_acc = 1e-5;
			}
		}
		else if (accuracy == 2) {
			if (atom_type_list[i] != 1) {
				lebedev_high = (max_l_temp < 3) ? lebedev_table[11] : lebedev_table[12];
				lebedev_low = (max_l_temp < 3) ? lebedev_table[7] : lebedev_table[8];
				radial_acc = 1e-5;
			}
			else {
				lebedev_high = (max_l_temp < 3) ? lebedev_table[10] : lebedev_table[11];
				lebedev_low = (max_l_temp < 3) ? lebedev_table[6] : lebedev_table[7];
				radial_acc = 1e-6;
			}
		}
		else if (accuracy == 3) {
			if (atom_type_list[i] != 1) {
				lebedev_high = (max_l_temp < 3) ? lebedev_table[14] : lebedev_table[16];
				lebedev_low = (max_l_temp < 3) ? lebedev_table[12] : lebedev_table[14];
				radial_acc = 1e-10;
			}
			else {
				lebedev_high = (max_l_temp < 3) ? lebedev_table[13] : lebedev_table[15];
				lebedev_low = (max_l_temp < 3) ? lebedev_table[11] : lebedev_table[13];
				radial_acc = 1e-11;
			}
		}
		else if (accuracy > 3) {
			if (atom_type_list[i] != 1) {
				lebedev_high = (max_l_temp < 3) ? lebedev_table[19] : lebedev_table[21];
				lebedev_low = (max_l_temp < 3) ? lebedev_table[14] : lebedev_table[17];
				radial_acc = 1e-20;
			}
			else {
				lebedev_high = (max_l_temp < 3) ? lebedev_table[18] : lebedev_table[20];
				lebedev_low = (max_l_temp < 3) ? lebedev_table[13] : lebedev_table[16];
				radial_acc = 1e-15;
			}
		}
		Prototype_grids.push_back(AtomGrid(radial_acc,
			lebedev_low,
			lebedev_high,
			atom_type_list[i],
			alpha_max_temp,
			max_l_temp,
			alpha_min_temp,
			debug,
			file));

	}

#ifdef _WIN64
	end_prototypes = time(NULL);
	if (debug) {

		for (int prototype = 0; prototype < Prototype_grids.size(); prototype++)
			file << "Number of gridpoints for atom type " << atom_type_list[prototype] << " :" << Prototype_grids[prototype].get_num_grid_points() << endl;

		//	int diff = end - start;
		if (end_prototypes - start < 1) file << "Time until prototypes are done: <1 s" << endl;
		else if (end_prototypes - start < 60) file << "Time until prototypes are done: " << fixed << setprecision(0) << end_prototypes - start << " s" << endl;
		else if (end_prototypes - start < 3600) file << "Time until prototypes are done: " << fixed << setprecision(0) << floor((end_prototypes - start) / 60) << " m " << (end_prototypes - start) % 60 << " s" << endl;
		else file << "Time until prototypes are done: " << fixed << setprecision(0) << floor((end_prototypes - start) / 3600) << " h " << ((end_prototypes - start) % 3600) / 60 << " m" << endl;
	}
#else
	if (debug) {
		gettimeofday(&t2, 0);

		double time3 = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000000.0;

		if (time3 < 60) printf("Time to prepare: %4.1lf s\n", time3);
		else if (time3 < 3600) printf("Time to prepare: %10.1lf m\n", time3 / 60);
		else printf("Time to prepare: %10.1lf h\n", time3 / 3600);
	}
#endif

#ifdef FLO_CUDA
	vector < vector < double > > radial_density;
	vector < vector < double > > radial_dist;

	radial_density.resize(atom_type_list.size());
	radial_dist.resize(atom_type_list.size());
	spherical_density.resize(asym_atom_list.size());
	for (int i = 0; i < asym_atom_list.size(); i++)
		spherical_density[i].resize(num_points[i]);

	const double incr = pow(1.005, max(1, accuracy - 1));
	const double lincr = log(incr);
	const double min_dist = 0.0000001;
	vector<Thakkar> sphericals;
	for (int i = 0; i < atom_type_list.size(); i++)
		sphericals.push_back(Thakkar(atom_type_list[i]));
	//Make radial grids
	if (debug) {
		for (int i = 0; i < atom_type_list.size(); i++) {
			file << "Calculating for atomic number " << atom_type_list[i] << endl;
			double current = 1;
			double dist = min_dist;
			if (accuracy > 3)
				while (current > 1E-10) {
					radial_dist[i].push_back(dist);
					current = sphericals[i].get_radial_density(dist);
					if (current == -20)
						return false;
					radial_density[i].push_back(current);
					dist *= incr;
				}
			else
				while (current > 1E-12) {
					radial_dist[i].push_back(dist);
					current = sphericals[i].get_radial_density(dist);
					if (current == -20)
						return false;
					radial_density[i].push_back(current);
					dist *= incr;
				}
			file << "Number of radial density points for atomic number " << atom_type_list[i] << ": " << radial_density[i].size() << endl;
			//for (int j = 0; j < radial_density[i].size(); j++) {
			//	if (radial_density[i][j] < 0.1)
			//		break;
			//	file << scientific << setprecision(8) << radial_density[i][j] << endl;
			//}
		}
	}
	else {
#pragma omp parallel for
		for (int i = 0; i < atom_type_list.size(); i++) {
			double current = 1;
			double dist = min_dist;
			if (accuracy > 3)
				while (current > 1E-10) {
					radial_dist[i].push_back(dist);
					current = sphericals[i].get_radial_density(dist);
					radial_density[i].push_back(current);
					dist *= incr;
				}
			else
				while (current > 1E-12) {
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
	for (int a = 0; a < wave.get_ncen(); a++) {
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
	for (int i = 0; i < asym_atom_list.size(); i++) {
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
	if (debug) file << "\nold Heapsize: " << size / 1024 / 1024 << " MB" << endl;
	float result;
	gpuErrchk(cudaOccupancyMaxPotentialBlockSize(
		&numBlocks,
		&blocks,
		(void*)gpu_calc_dens_per_MO_static1,
		0,
		MaxGrid));
	gridSize = (MaxGrid + blocks - 1) / blocks;
	result = ((sizeof(int) * 10 + sizeof(float) * (6 * ncen_temp + 20)) * blocks * prop[dev].multiProcessorCount) / 1024 / 1024;

	if (debug) file << "result: " << fixed << result << " MB for " << gridSize << " blocks with " << blocks << " threads" << endl << "sizeof float: " << sizeof(float) << " B" << endl;
	if (result > size / 1024 / 1024)
		gpuErrchk(cudaDeviceSetLimit(cudaLimitMallocHeapSize, result * 1024 * 1024 * 1.1));

	gpuErrchk(cudaDeviceGetLimit(&size, cudaLimitMallocHeapSize));
	if (debug) file << "new Heapsize: " << size / 1024 / 1024 << " MB" << endl;

	gridSize = (MaxGrid + blocks - 1) / blocks;

	//Allocate and copy vectors to device for WFN
	if (debug) file << "Copying WFN to devices now!" << endl;
	gpuErrchk(cudaMalloc((void**)&gpu_PosAtomsx[0],    sizeof(float) * ncen_temp           ));
	gpuErrchk(cudaMalloc((void**)&gpu_PosAtomsy[0],    sizeof(float) * ncen_temp           ));
	gpuErrchk(cudaMalloc((void**)&gpu_PosAtomsz[0],    sizeof(float) * ncen_temp           ));
	gpuErrchk(cudaMalloc((void**)&gpu_atom_z[0],       sizeof(int)   * ncen_temp           ));
	gpuErrchk(cudaMalloc((void**)&gpu_Gridx[0],        sizeof(float) * MaxGrid             ));
	gpuErrchk(cudaMalloc((void**)&gpu_Gridy[0],        sizeof(float) * MaxGrid             ));
	gpuErrchk(cudaMalloc((void**)&gpu_Gridz[0],        sizeof(float) * MaxGrid             ));
	gpuErrchk(cudaMalloc((void**)&gpu_Gridaw[0],       sizeof(float) * MaxGrid             ));
	gpuErrchk(cudaMalloc((void**)&gpu_Gridmw[0],       sizeof(float) * MaxGrid             ));
	gpuErrchk(cudaMalloc((void**)&gpu_asym_atom_list[0], sizeof(int) * asym_atom_list.size()));
	gpuErrchk(cudaMalloc((void**)&gpu_atom_type_list[0], sizeof(int) * atom_type_list.size()));
	gpuErrchk(cudaMalloc((void**)&gpu_numpoints[0],      sizeof(int) * asym_atom_list.size()));
	if (debug) file << "Mallocs done!" << endl;
	gpuErrchk(cudaMemcpyToSymbol(gpu_nex, &nex_temp, sizeof(int)));
	gpuErrchk(cudaMemcpyToSymbol(gpu_nmo, &nmo_temp, sizeof(int)));
	gpuErrchk(cudaMemcpyToSymbol(gpu_ncen, &ncen_temp, sizeof(int)));
	gpuErrchk(cudaMemcpyToSymbol(gpu_MaxGrid, &MaxGrid, sizeof(int)));
	gpuErrchk(cudaMemcpyToSymbol(gpu_start_radial_dens, &min_dist, sizeof(float)));
	gpuErrchk(cudaMemcpyToSymbol(gpu_log_incr, &lincr, sizeof(float)));
	gpuErrchk(cudaMemcpy(gpu_PosAtomsx[0],	    PosAtoms[0].data(),       sizeof(float) * wave.get_ncen(),     cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(gpu_PosAtomsy[0],	    PosAtoms[1].data(),       sizeof(float) * wave.get_ncen(),     cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(gpu_PosAtomsz[0],	    PosAtoms[2].data(),       sizeof(float) * wave.get_ncen(),     cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(gpu_asym_atom_list[0], asym_atom_list.data(),    sizeof(int) * asym_atom_list.size(),  cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(gpu_atom_type_list[0], atom_type_list.data(),    sizeof(int) * atom_type_list.size(),  cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(gpu_numpoints[0],      num_points.data(),        sizeof(int) * asym_atom_list.size(),  cudaMemcpyHostToDevice));
	gpuErrchk(cudaPeekAtLastError());
	file << "All copying done!" << endl;

	vector <cudaStream_t> streams;
	streams.resize(asym_atom_list.size());
	int offset = 0;
	for (int i = 0; i < asym_atom_list.size(); i++) {
		gpuErrchk(cudaOccupancyMaxPotentialBlockSize(
			&numBlocks,
			&blocks,
			(void*)gpu_make_grid,
			0,
			num_points[i]));

		gridSize = (MaxGrid + blocks - 1) / blocks;
		if(debug) file << i << ": num points: " << num_points[i] << " blocks: " << gridSize << " threads: " << blocks << " grid > points? " << (blocks*gridSize > num_points[i]) << endl;
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

		gpu_make_grid <<< gridSize, blocks, 0, streams[i] >>> (
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
	
	for (int i = 0; i < atom_type_list.size(); i++) {
		gpuErrchk(cudaFree(gpu_atomgrid_w[0][i]));
		gpuErrchk(cudaFree(gpu_atomgrid_x[0][i]));
		gpuErrchk(cudaFree(gpu_atomgrid_y[0][i]));
		gpuErrchk(cudaFree(gpu_atomgrid_z[0][i]));
	}
#else
	vector < vector < double > > radial_density;
	vector < vector < double > > radial_dist;
	
	radial_density.resize(atom_type_list.size());
	radial_dist.resize(atom_type_list.size());
	spherical_density.resize(asym_atom_list.size());
	for (int i = 0; i < asym_atom_list.size(); i++)
		spherical_density[i].resize(num_points[i]);
	
	const double incr = pow(1.005, max(1, accuracy - 1));
	const double lincr = log(incr);
	const double min_dist = 0.0000001;
	vector<Thakkar> sphericals;
	for (int i = 0; i < atom_type_list.size(); i++)
		sphericals.push_back(Thakkar(atom_type_list[i]));
	//Make radial grids
	if (debug) {
		for (int i = 0; i < atom_type_list.size(); i++) {
			file << "Calculating for atomic number " << atom_type_list[i] << endl;
			double current = 1;
			double dist = min_dist;
			if (accuracy > 3)
				while (current > 1E-10) {
					radial_dist[i].push_back(dist);
					current = sphericals[i].get_radial_density(dist);
					if (current == -20)
						return false;
					radial_density[i].push_back(current);
					dist *= incr;
				}
			else
				while (current > 1E-12) {
					radial_dist[i].push_back(dist);
					current = sphericals[i].get_radial_density(dist);
					if (current == -20)
						return false;
					radial_density[i].push_back(current);
					dist *= incr;
				}
			file << "Number of radial density points for atomic number " << atom_type_list[i] << ": " << radial_density[i].size() << endl;
			//for (int j = 0; j < radial_density[i].size(); j++) {
			//	if (radial_density[i][j] < 0.1)
			//		break;
			//	file << scientific << setprecision(8) << radial_density[i][j] << endl;
			//}
		}
	}
	else {
#pragma omp parallel for
		for (int i = 0; i < atom_type_list.size(); i++) {
			double current = 1;
			double dist = min_dist;
			if (accuracy > 3)
				while (current > 1E-10) {
					radial_dist[i].push_back(dist);
					current = sphericals[i].get_radial_density(dist);
					radial_density[i].push_back(current);
					dist *= incr;
				}
			else
				while (current > 1E-12) {
					radial_dist[i].push_back(dist);
					current = sphericals[i].get_radial_density(dist);
					radial_density[i].push_back(current);
					dist *= incr;
				}
		}
	}
	sphericals.clear();

	for (int i = 0; i < asym_atom_list.size(); i++) {
		int nr = asym_atom_list[i];
		int type;
		for (int j = 0; j < atom_type_list.size(); j++)
			if (atom_type_list[j] == wave.atoms[nr].charge)
				type = j;

		num_points[i] = Prototype_grids[type].get_num_grid_points();

		if (debug) file << "Number of points for atom " << i << ": " << num_points[i] << endl;

		for (int n = 0; n < 6; n++)
			grid[i][n] = new double[num_points[i]];
#pragma omp parallel for schedule(dynamic)
		for (int p = 0; p < num_points[i]; p++)
			grid[i][4][p] = 0.0;

		Prototype_grids[type].get_grid(int(wave.get_ncen() * pow(pbc * 2 + 1, 3)),
			nr,
			&x[0],
			&y[0],
			&z[0],
			&atom_z[0],
			grid[i][0],
			grid[i][1],
			grid[i][2],
			grid[i][3],
			grid[i][5]);
	}
	Prototype_grids.clear();
#endif

	int points = 0;
	for (int i = 0; i < asym_atom_list.size(); i++)
		points += num_points[i];
	if (debug) file << "Becke Grid exists" << endl;
	else file << "                                 done! Number of gridpoints: " << defaultfloat << points << endl;

#ifdef _WIN64
	if (debug) {
		file << "Taking time..." << endl;
	}
	end_becke = time(NULL);

#else
	if (debug) {
		gettimeofday(&t2, 0);

		double timeb = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000000.0;

		if (timeb < 60) printf("Time to prepare: %4.1lf s\n", timeb);
		else if (timeb < 3600) printf("Time to prepare: %10.1lf m\n", timeb / 60);
		else printf("Time to prepare: %10.1lf h\n", timeb / 3600);
	}
#endif

	file << "Calculating spherical densities..." << flush;

#ifdef FLO_CUDA
	float*** gpu_spherical_density = NULL,
		*** gpu_radial_density = NULL,
		*** gpu_radial_dist = NULL,
		** gpu_Grids = NULL;
	gpu_radial_density =    (float***)malloc(sizeof(float**));
	gpu_radial_dist =       (float***)malloc(sizeof(float**));
	gpu_spherical_density = (float***)malloc(sizeof(float**));

	gpu_radial_density[0] =    (float**)malloc(sizeof(float*) * atom_type_list.size());
	gpu_radial_dist[0] =       (float**)malloc(sizeof(float*) * atom_type_list.size());
	gpu_spherical_density[0] = (float**)malloc(sizeof(float*) * asym_atom_list.size());
	gpu_Grids =                (float**)malloc(sizeof(float*));

	for (int i = 0; i < asym_atom_list.size(); i++)
		gpuErrchk(cudaMalloc((void**)&(gpu_spherical_density[0][i]), sizeof(float) * num_points[i]));

	gpuErrchk(cudaMalloc((void**)&gpu_Grids[0], sizeof(float) * MaxGrid));

	for (int i = 0; i < atom_type_list.size(); i++) {
		gpuErrchk(cudaMalloc((void**)&gpu_radial_density[0][i], sizeof(float) * radial_density[i].size()));
		gpuErrchk(cudaMalloc((void**)&gpu_radial_dist[0][i],    sizeof(float) * radial_dist[i].size()));
		gpuErrchk(cudaMemcpy(gpu_radial_density[0][i], radial_density[i].data(), sizeof(float) * radial_density[i].size(), cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(gpu_radial_dist[0][i], radial_dist[i].data(),       sizeof(float) * radial_dist[i].size(), cudaMemcpyHostToDevice));
	}
	gpuErrchk(cudaDeviceSynchronize());
	gpuErrchk(cudaPeekAtLastError());
	for (int i = 0; i < wave.get_ncen(); i++) {
		int nr = all_atom_list[i];
		int type_list_number = -1;
		for (int j = 0; j < atom_type_list.size(); j++)
			if (wave.atoms[nr].charge == atom_type_list[j])
				type_list_number = j;
		offset = 0;
		for (int g = 0; g < asym_atom_list.size(); g++) {
			if (g == 0) {
				gpuErrchk(cudaOccupancyMaxPotentialBlockSize(
					&numBlocks,
					&blocks,
					(void*)gpu_linear_interpolate_spherical_density,
					0,
					num_points[i]));

				gridSize = (MaxGrid + blocks - 1) / blocks;
				//file << i << "/" << g << ": blocks: " << gridSize << " threads: " << blocks << endl;
			}
			bool match = (all_atom_list[i] == asym_atom_list[g]);
			gpu_linear_interpolate_spherical_density <<< gridSize, blocks >>> (
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

	if (debug) {
		//Copy grid from GPU to print:
		// Dimensions: [c] [p]
		// p = the number of gridpoint
		// c = coordinate, which is 0=x, 1=y, 2=z, 3=atomic becke weight, 4=spherical density, 5=wavefunction density, 6=molecular becke weight
		for (int i = 0; i < 7; i++)
			total_grid[i].resize(MaxGrid);
		gpuErrchk(cudaMemcpy(total_grid[0].data(), gpu_Gridx[0], sizeof(float) * MaxGrid, cudaMemcpyDeviceToHost));
		gpuErrchk(cudaMemcpy(total_grid[1].data(), gpu_Gridy[0], sizeof(float) * MaxGrid, cudaMemcpyDeviceToHost));
		gpuErrchk(cudaMemcpy(total_grid[2].data(), gpu_Gridz[0], sizeof(float) * MaxGrid, cudaMemcpyDeviceToHost));
		gpuErrchk(cudaMemcpy(total_grid[3].data(), gpu_Gridaw[0], sizeof(float) * MaxGrid, cudaMemcpyDeviceToHost));
		gpuErrchk(cudaMemcpy(total_grid[6].data(), gpu_Gridmw[0], sizeof(float) * MaxGrid, cudaMemcpyDeviceToHost));
		gpuErrchk(cudaMemcpy(total_grid[4].data(), gpu_Grids[0], sizeof(float)* MaxGrid, cudaMemcpyDeviceToHost));

		gpuErrchk(cudaDeviceSynchronize());
		gpuErrchk(cudaPeekAtLastError());

		ofstream grid0("grid0.file", ios::out);
		for (int i = 0; i < MaxGrid; i++) {
			for (int j = 0; j < 7; j++)
				grid0 << setw(16) << scientific << setprecision(8) << total_grid[j][i];
			grid0 << "\n";
		}
		grid0.flush();
		grid0.close();
	}

	for(int i=0; i< atom_type_list.size(); i++){
		gpuErrchk(cudaFree(gpu_radial_density[0][i]));
		gpuErrchk(cudaFree(gpu_radial_dist[0][i]));
	}

#else

	for (int i = 0; i < wave.get_ncen(); i++) {
		//int nr = all_atom_list[i];
		int type_list_number = -1;
		for (int j = 0; j < atom_type_list.size(); j++)
			if (wave.atoms[i].charge == atom_type_list[j])
				type_list_number = j;
		for (int g = 0; g < asym_atom_list.size(); g++) {
			if (i == 0) spherical_density[g].resize(num_points[g]);
#pragma omp parallel for
			for (int p = 0; p < num_points[g]; p++) {
				double temp =
					linear_interpolate_spherical_density(radial_density[type_list_number]
						, radial_dist[type_list_number]
						, sqrt(pow(grid[g][0][p] - wave.atoms[i].x, 2)
							+ pow(grid[g][1][p] - wave.atoms[i].y, 2)
							+ pow(grid[g][2][p] - wave.atoms[i].z, 2))
						, lincr
						, min_dist
					);
				if (i == asym_atom_list[g])
					spherical_density[g][p] = temp;
				grid[g][4][p] += temp;
			}
		}
	}
	//fill out with priodic information
	if (pbc != 0) {
		for (int x = -pbc; x < pbc + 1; x++)
			for (int y = -pbc; y < pbc + 1; y++)
				for (int z = -pbc; z < pbc + 1; z++) {
					if (x == 0 && y == 0 && z == 0)
						continue;
					for (int i = 0; i < wave.get_ncen(); i++) {
						int type_list_number = -1;
						//int nr = all_atom_list[i];
						for (int j = 0; j < atom_type_list.size(); j++)
							if (wave.atoms[i].charge == atom_type_list[j])
								type_list_number = j;
						for (int g = 0; g < asym_atom_list.size(); g++)
#pragma omp parallel for
							for (int p = 0; p < num_points[g]; p++)
								grid[g][4][p] += linear_interpolate_spherical_density(
									radial_density[type_list_number],
									radial_dist[type_list_number],
									sqrt(
										  pow(grid[g][0][p] - (wave.atoms[i].x + x * unit_cell.get_cm(0, 0) + y * unit_cell.get_cm(0, 1) + z * unit_cell.get_cm(0, 2)), 2)
										+ pow(grid[g][1][p] - (wave.atoms[i].y + x * unit_cell.get_cm(1, 0) + y * unit_cell.get_cm(1, 1) + z * unit_cell.get_cm(1, 2)), 2)
										+ pow(grid[g][2][p] - (wave.atoms[i].z + x * unit_cell.get_cm(2, 0) + y * unit_cell.get_cm(2, 1) + z * unit_cell.get_cm(2, 2)), 2)
									),
									lincr,
									min_dist
								);
					}
				}
	}
	if (debug)
		for (int i = 0; i < asym_atom_list.size(); i++)
			file << "number of points for atom " << i << " " << num_points[i] << " " << spherical_density[i].size() << endl;
#endif

	file << "                    done!" << endl;
	
	radial_density.clear();
	radial_dist.clear();
		
#ifdef _WIN64
	if (debug) {
		file << "Taking time..." << endl;
	}
	end_spherical = time(NULL);

#else
	if (debug) {
		gettimeofday(&t2, 0);

		double timeb = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000000.0;

		if (timeb < 60) printf("Time to prepare: %4.1lf s\n", timeb);
		else if (timeb < 3600) printf("Time to prepare: %10.1lf m\n", timeb / 60);
		else printf("Time to prepare: %10.1lf h\n", timeb / 3600);
	}
#endif

	double cutoff;
	if (accuracy < 3)
		cutoff = 1E-10;
	else if (accuracy == 3)
		cutoff = 1E-14;
	else
		cutoff = 1E-30;
#ifndef FLO_CUDA
	if (!debug) {
		file << "Pruning Grid..." << flush;
		for (int i = 0; i < asym_atom_list.size(); i++) {
			int reduction = 0;
			for (int p = 0; p < num_points[i]; p++) {
				if (grid[i][4][p] != 0.0) {
					if (abs(grid[i][3][p] * spherical_density[i][p - reduction] / grid[i][4][p]) > cutoff) {
						for (int k = 0; k < 5; k++)
							total_grid[k].push_back(grid[i][k][p]);
						total_grid[6].push_back(grid[i][5][p]);
					}
					else {
						spherical_density[i].erase(spherical_density[i].begin() + (p - reduction));
						reduction++;
					}
				}
				else {
					spherical_density[i].erase(spherical_density[i].begin() + (p - reduction));
					reduction++;
				}

			}
			num_points[i] -= reduction;
			delete[](grid[i]);
			if (debug) file << endl << "number of points for atom " << i << " " << num_points[i] << " " << spherical_density[i].size() << endl;
		}
	}
	else {
		for (int i = 0; i < asym_atom_list.size(); i++) {
			for (int p = 0; p < num_points[i]; p++) {
				for (int k = 0; k < 5; k++)
					total_grid[k].push_back(grid[i][k][p]);
				total_grid[6].push_back(grid[i][5][p]);
			}
			delete[](grid[i]);
			if (debug) file << endl << "number of points for atom " << i << " " << num_points[i] << " " << spherical_density[i].size() << endl;
		}
	}
	points = 0;
	for (int i = 0; i < asym_atom_list.size(); i++)
		points += num_points[i];

	total_grid[5].resize(total_grid[0].size());
	if (debug) file << "sphericals done!" << endl;
	else file << "                                       done! Number of gridpoints: " << defaultfloat << points << endl;
#endif

#ifdef _WIN64
	if (debug) {
		file << "Taking time..." << endl;
	}
	end_prune = time(NULL);

#else
	if (debug) {
		gettimeofday(&t2, 0);

		double timeb = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000000.0;

		if (timeb < 60) printf("Time to prepare: %4.1lf s\n", timeb);
		else if (timeb < 3600) printf("Time to prepare: %10.1lf m\n", timeb / 60);
		else printf("Time to prepare: %10.1lf h\n", timeb / 3600);
	}
#endif

	file << "Calculating non-spherical densities..." << flush;
    vector < vector < double > > periodic_grid;

#ifdef FLO_CUDA
	// Vector containing integrated numbers of electrons
	// dimension 0: 0=Becke grid integration 1=Summed spherical density 2=hirshfeld weighted density
	// dimension 1: atoms of asym_atom_list
	vector < vector <double> > atom_els;
	atom_els.resize(3);
	for (int i = 0; i < asym_atom_list.size(); i++)
		for (int n = 0; n < 3; n++)
			atom_els[n].push_back(0.0);	
	
	vector <float> coef;
	vector <float> ex;
	for (int i = 0; i < nex_temp; i++)
		ex.push_back(wave.get_exponent(i));

	gpuErrchk(cudaMalloc((void**)&gpu_types[0], sizeof(int)* nex_temp));
	gpuErrchk(cudaMalloc((void**)&gpu_centers[0], sizeof(int)* nex_temp));
	gpuErrchk(cudaMalloc((void**)&gpu_exponents[0], sizeof(float)* nex_temp));
	gpuErrchk(cudaMalloc((void**)&gpu_GridRho[0], sizeof(float)* MaxGrid));
	gpuErrchk(cudaMemcpy(gpu_types[0], wave.get_ptr_types(), sizeof(int)* nex_temp, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(gpu_centers[0], wave.get_ptr_centers(), sizeof(int)* nex_temp, cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpy(gpu_exponents[0], ex.data(), sizeof(float)* nex_temp, cudaMemcpyHostToDevice));
	
	const bool per_MO = true;
	if (per_MO) {
		for (int mo = 0; mo < wave.get_nmo(false); mo++)
			for (int i = 0; i < nex_temp; i++)
				if (wave.get_MO_occ(mo) != 0)
					coef.push_back(wave.get_MO_coef(mo, i, debug));
		if (debug) file << "Number of coefs: " << coef.size() << endl;
		gpuErrchk(cudaMalloc((void**)&gpu_coefficients[0], sizeof(float)* nex_temp));
		gpuErrchk(cudaMemset(gpu_GridRho[0], 0.0, sizeof(float) * MaxGrid));
		if (nex_temp < 10000) {
			file << "Using shared memory kernel" << endl;
			gpuErrchk(cudaOccupancyMaxPotentialBlockSize(
				&numBlocks,
				&blocks,
				(void*)gpu_calc_dens_per_MO_shared,
				0,
				MaxGrid));
		}
		else {
			file << "Using static memory kernel" << endl;
			gpuErrchk(cudaOccupancyMaxPotentialBlockSize(
				&numBlocks,
				&blocks,
				(void*)gpu_calc_dens_per_MO_static1,
				0,
				MaxGrid));
		}

		gridSize = (MaxGrid + blocks - 1) / blocks;
		if (debug) file << "running " << gridSize << " blocks with " << blocks << " threads" << endl;
		unsigned int nex_offset = 0;

		for (int mo = 0; mo < nmo_temp; mo++) {
			gpuErrchk(cudaPeekAtLastError());
			gpuErrchk(cudaDeviceSynchronize());
			
			if(debug) file << "MO " << mo << " starting at coef: " << nex_offset << endl;
			gpuErrchk(cudaMemcpy(gpu_coefficients[0], &coef[nex_offset], sizeof(float) * nex_temp, cudaMemcpyHostToDevice));
			if (nex_temp < 10000) {
				gpu_calc_dens_per_MO_shared <<< gridSize, blocks, nex_temp * sizeof(float) >>> (
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
			else {
				gpu_calc_dens_per_MO_static2 <<< gridSize, blocks >>> (
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
	else {
		for (int i = 0; i < nex_temp; i++) 
			for (int mo = 0; mo < wave.get_nmo(false); mo++) 
				if (wave.get_MO_occ(mo) != 0)
					coef.push_back(wave.get_MO_coef(mo, i, debug));
		//seems broken and is slower, especially if L1 is sufficient for coefficients with size nex_temp
		gpuErrchk(cudaMalloc((void**)&gpu_coefficients[0], sizeof(float)* nex_temp * nmo_temp));
		gpuErrchk(cudaMemcpy(gpu_coefficients[0], coef.data(), sizeof(float) * nex_temp * nmo_temp, cudaMemcpyHostToDevice));
		vector <float> occ;
		for (int i = 0; i < wave.get_nmo(false); i++) {
			occ.push_back(wave.get_MO_occ(i));
			if (occ[occ.size() - 1] == 0) occ.pop_back();
		}
		gpuErrchk(cudaMalloc((void**)&gpu_occ[0], sizeof(float)* nmo_temp));
		gpuErrchk(cudaMemcpy(gpu_occ[0], occ.data(), sizeof(float)* occ.size(), cudaMemcpyHostToDevice));
		gpuErrchk(cudaOccupancyMaxPotentialBlockSize(
			&numBlocks,
			&blocks,
			(void*)gpu_calc_dens,
			0,
			MaxGrid));

		gridSize = (MaxGrid + blocks - 1) / blocks;
		gpu_calc_dens <<<gridSize, blocks>>> (
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

	if (debug) {
		//Copy grid from GPU to print:
		// Dimensions: [c] [p]
		// p = the number of gridpoint
		// c = coordinate, which is 0=x, 1=y, 2=z, 3=atomic becke weight, 4=spherical density, 5=wavefunction density, 6=molecular becke weight
		gpuErrchk(cudaMemcpy(total_grid[5].data(), gpu_GridRho[0], sizeof(float) * MaxGrid, cudaMemcpyDeviceToHost));
		ofstream grid("grid.file", ios::out);
		for (int i = 0; i < MaxGrid; i++) {
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
	for (int i = 0; i < asym_atom_list.size(); i++) {
		double charges[3];
		double* result = (double*)malloc(sizeof(double));
		gpuErrchk(cudaMalloc((void**)&(result), sizeof(double)*3));
		gpu_calc_charges <<<1, 1 >>> (
			gpu_GridRho[0],
			gpu_Gridaw[0],
			gpu_Gridmw[0],
			gpu_Grids[0],
			gpu_spherical_density[0][i],
			num_points[i],
			offset,
			cutoff,
			result
			);
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
	for (int i = 0; i < asym_atom_list.size(); i++) {
		cudaOccupancyMaxPotentialBlockSize(
			&numBlocks,
			&blocks,
			(void*)gpu_apply_weights,
			0,
			num_points[i]);

		gridSize = (MaxGrid + blocks - 1) / blocks;
		//file << i << ": blocks: " << gridSize << " threads: " << blocks << endl;

		gpu_apply_weights <<< gridSize, blocks, 0, streams[i] >>> (
			offset,
			gpu_GridRho[0],
			gpu_Grids[0],
			gpu_spherical_density[0][i],
			gpu_Gridaw[0],
			num_points[i]
			);

		offset += num_points[i];
	}
	gpuErrchk(cudaDeviceSynchronize());
	gpuErrchk(cudaPeekAtLastError());
	file << "After Applying!" << endl;

	file << endl;
	file << " done!" << endl;
	file << "Number of points evaluated: " << MaxGrid;
#else
	double** mo_coef;
	mo_coef = new double* [wave.get_nmo()];
	for (int i = 0; i < wave.get_nmo(); i++)
		mo_coef[i] = wave.get_MO_coef_ptr(i);
#pragma omp parallel for
	for (int i = 0; i < total_grid[0].size(); i++) {
		total_grid[5][i] = compute_dens(wave, new double[3]{ 
            total_grid[0][i], 
            total_grid[1][i], 
            total_grid[2][i]},
			mo_coef);
	}
	if (debug) {
		//Copy grid from GPU to print:
		// Dimensions: [c] [p]
		// p = the number of gridpoint
		// c = coordinate, which is 0=x, 1=y, 2=z, 3=atomic becke weight, 4=spherical density, 5=wavefunction density, 6=molecular becke weight
		ofstream grid("grid.file", ios::out);
		for (int i = 0; i < total_grid[0].size(); i++) {
			for (int j = 0; j < 7; j++)
				grid << setw(16) << scientific << setprecision(8) << total_grid[j][i];
			grid << "\n";
		}
		grid.flush();
		grid.close();
	}
	if (pbc != 0) {
        periodic_grid.resize(pow(pbc*2+1,3));
        int j=0;
        for(int d=0; d<pow(pbc*2+1,3); d++)
            periodic_grid[d].resize(total_grid[5].size());
		for (int x = -pbc; x < pbc+1; x++)
			for (int y = -pbc; y < pbc+1; y++)
				for (int z = -pbc; z < pbc+1; z++) {
					if (x == 0 && y == 0 && z == 0)
						continue;
#pragma omp parallel for
					for (int i = 0; i < total_grid[0].size(); i++)
						periodic_grid[j][i] = compute_dens(wave, new double[3]{
							total_grid[0][i] + x * unit_cell.get_cm(0,0) + y * unit_cell.get_cm(0,1) + z * unit_cell.get_cm(0,2),
							total_grid[1][i] + x * unit_cell.get_cm(1,0) + y * unit_cell.get_cm(1,1) + z * unit_cell.get_cm(1,2),
							total_grid[2][i] + x * unit_cell.get_cm(2,0) + y * unit_cell.get_cm(2,1) + z * unit_cell.get_cm(2,2) },
							mo_coef);
                    j++;
				}
		if (debug){
			for (int i = 0; i < total_grid[0].size(); i++){
                if (i%1000==0)
                    file << "Old dens: " << total_grid[5][i] << " contributions of neighbour-cells:"; 
                for (int j=0; j<pow(pbc*2+1,3)-1; j++){
                    if (i%1000==0)
                        file << " " << periodic_grid[j][i];
                    total_grid[5][i] += periodic_grid[j][i];
                }
                if (i%1000==0)
				    file << endl;
            }
        }
	}

	for (int i = 0; i < wave.get_nmo(); i++)
		mo_coef[i] = nullptr;
	delete[] mo_coef;

	if (debug) file << endl << "with total number of points: " << total_grid[0].size() << endl;
	else file << "                done!" << endl;

	file << "Applying hirshfeld weights and integrating charges..." << flush;
	double el_sum_becke = 0.0;
	double el_sum_spherical = 0.0;
	double el_sum_hirshfeld = 0.0;
	// Vector containing integrated numbers of electrons
	// dimension 0: 0=Becke grid integration 1=Summed spherical density 2=hirshfeld weighted density
	// dimension 1: atoms of asym_atom_list
	vector < vector <double> > atom_els;
	atom_els.resize(3);
	for (int i = 0; i < asym_atom_list.size(); i++)
		for (int n = 0; n < 3; n++)
			atom_els[n].push_back(0.0);

	//Generate Electron sums
#pragma omp parallel for reduction(+:el_sum_becke,el_sum_spherical,el_sum_hirshfeld)
	for (int i = 0; i < asym_atom_list.size(); i++) {
		int start_p = 0;
#pragma loop(no_vector)
		for (int a = 0; a < i; a++)
			start_p += num_points[a];
		for (int p = start_p; p < start_p + num_points[i]; p++) {
			if (abs(total_grid[6][p]) > cutoff) {
				atom_els[0][i] += total_grid[6][p] * total_grid[5][p];
				atom_els[1][i] += total_grid[6][p] * total_grid[4][p];
			}
			atom_els[2][i] += total_grid[5][p] * total_grid[3][p] * spherical_density[i][p - start_p] / total_grid[4][p];
		}
		el_sum_becke += atom_els[0][i];
		el_sum_spherical += atom_els[1][i];
		el_sum_hirshfeld += atom_els[2][i];
	}

	if (debug) {
		file << "Becke grid with hirshfeld weights done!" << endl;
		file << "atom_els[2]: ";
		for (int i = 0; i < asym_atom_list.size(); i++) {
			if (isnan(atom_els[2][i]))
				file << "!!!";
			file << fixed << setw(10) << setprecision(3) << atom_els[2][i] << " ";
			if (isnan(atom_els[2][i]))
				file << "!!!";
		}
		file << endl;
	}

#pragma omp parallel for
	for (int p = 0; p < total_grid[0].size(); p++)
		total_grid[5][p] *= total_grid[3][p];
	file << " done!" << endl;
	file << "Number of points evaluated: " << total_grid[0].size();
#endif
	
	file << " with " << fixed << setw(10) << setprecision(6) << el_sum_becke << " electrons in Becke Grid in total." << endl << endl;

	file << "Table of Charges in electrons" << endl << endl << "Atom       Becke   Spherical Hirshfeld" << endl;

	int counter = 0;
	for (int i = 0; i < wave.get_ncen(); i++) {
		if (is_asym[i]) {
			file << setw(6) << labels[i]
				<< fixed << setw(10) << setprecision(3) << wave.atoms[i].charge - atom_els[0][counter]
				<< fixed << setw(10) << setprecision(3) << wave.atoms[i].charge - atom_els[1][counter]
				<< fixed << setw(10) << setprecision(3) << wave.atoms[i].charge - atom_els[2][counter];
			if (debug) file << " " << setw(4) << wave.atoms[i].charge << " " << fixed << setw(10) << setprecision(3) << wave.atoms[i].charge - atom_els[0][counter]
				<< fixed << setw(10) << setprecision(3) << atom_els[1][counter]
				<< fixed << setw(10) << setprecision(3) << atom_els[2][counter];
			counter++;
			file << endl;
		}
	}

	file << "Total number of electrons in the wavefunction: " << el_sum_becke << endl << " and Hirshfeld electrons (asym unit): " << el_sum_hirshfeld << endl;


#ifdef _WIN64
	if (debug) {
		file << "Taking time..." << endl;
	}
	end_aspherical = time(NULL);

#else
	if (debug) {
		gettimeofday(&t2, 0);

		double timeb = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000000.0;

		if (timeb < 60) printf("Time to prepare: %4.1lf s\n", timeb);
		else if (timeb < 3600) printf("Time to prepare: %10.1lf m\n", timeb / 60);
		else printf("Time to prepare: %10.1lf h\n", timeb / 3600);
	}
#endif
	vector < vector <double> > d1, d2, d3;
	vector < vector <double> > dens;
	dens.resize(asym_atom_list.size());
	if (debug)
		file << "resized outer dens" << endl;
	d1.resize(asym_atom_list.size());
	d2.resize(asym_atom_list.size());
	d3.resize(asym_atom_list.size());
	if (debug)
		file << "resized outer d1-3" << endl;
#ifdef FLO_CUDA
	points = 0;
#pragma omp parallel for
	for (int i = 0; i < asym_atom_list.size(); i++) {
		int type;
		for (int j = 0; j < atom_type_list.size(); j++)
			if (atom_type_list[j] == wave.atoms[asym_atom_list[i]].charge)
				type = j;
		vector<float> temp_dens;
		temp_dens.resize(num_points[i]);
		//file << "At atom: " << i;
		gpuErrchk(cudaMemcpy(temp_dens.data(), gpu_spherical_density[0][i], sizeof(float) * num_points[i], cudaMemcpyDeviceToHost));
		for (int p = 0; p < 0 + num_points[i]; p++) {
			if (abs(temp_dens[p]) > cutoff) {
				dens[i].push_back(temp_dens[p]);
				d1[i].push_back(Prototype_grids[type].get_gridx(p));
				d2[i].push_back(Prototype_grids[type].get_gridy(p));
				d3[i].push_back(Prototype_grids[type].get_gridz(p));
			}
		}
		//file << " dens size: " << dens[i].size() << " num_points[i]: " << num_points[i] << endl;
	}
	for (int i = 0; i < asym_atom_list.size(); i++) points += dens[i].size();

	//gpuErrchk(cudaFree(gpu_Grids[0]));
	//gpuErrchk(cudaFree(gpu_PosAtomsx[0]));
	//gpuErrchk(cudaFree(gpu_PosAtomsy[0]));
	//gpuErrchk(cudaFree(gpu_PosAtomsz[0]));
	//gpuErrchk(cudaFree(gpu_GridRho[0]));
	//gpuErrchk(cudaFree(gpu_Gridx[0]));
	//gpuErrchk(cudaFree(gpu_Gridy[0]));
	//gpuErrchk(cudaFree(gpu_Gridz[0]));
	//gpuErrchk(cudaFree(gpu_Gridaw[0]));
	//gpuErrchk(cudaFree(gpu_Gridmw[0]));
	//gpuErrchk(cudaFree(gpu_asym_atom_list[0]));
	//gpuErrchk(cudaFree(gpu_atom_type_list[0]));
	//gpuErrchk(cudaFree(gpu_numpoints[0]));
	//gpuErrchk(cudaFree(gpu_atom_z[0]));
	//for (int i = 0; i < asym_atom_list.size(); i++) {
	//	gpuErrchk(cudaFree(gpu_spherical_density[0][i]));
	//	gpuErrchk(cudaFree(gpu_atomgrid_x[0][i]));
	//	gpuErrchk(cudaFree(gpu_atomgrid_y[0][i]));
	//	gpuErrchk(cudaFree(gpu_atomgrid_z[0][i]));
	//}
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
#else
	points = 0;
#pragma omp parallel for reduction(+:points)
	for (int i = 0; i < asym_atom_list.size(); i++) {
		int start_p = 0;
		double res;
#pragma loop(no_vector)
		for (int a = 0; a < i; a++)
			start_p += num_points[a];
		for (int p = start_p; p < start_p + num_points[i]; p++) {
			res = total_grid[5][p] * spherical_density[i][p - start_p] / total_grid[4][p];
			if (abs(res) > cutoff) {
				dens[i].push_back(res);
				d1[i].push_back(total_grid[0][p] - wave.atoms[asym_atom_list[i]].x);
				d2[i].push_back(total_grid[1][p] - wave.atoms[asym_atom_list[i]].y);
				d3[i].push_back(total_grid[2][p] - wave.atoms[asym_atom_list[i]].z);
			}
		}
		points += dens[i].size();
	}

	for (int i = 0; i < asym_atom_list.size(); i++)
		spherical_density[i].clear();
	spherical_density.clear();

	for (int grid = 0; grid < total_grid.size(); grid++)
		total_grid[grid].resize(0);
	total_grid.resize(0);
#endif

	vector < vector <double> > k_pt;
	k_pt.resize(3);
#pragma omp parallel for
	for (int i = 0; i < 3; i++)
		k_pt[i].resize(sym[0][0].size()* hkl[0].size());

	if (debug)
		file << "K_point_vector is here! size: " << k_pt[0].size() << endl;

	const bool shrink = true;
	vector < vector<double> > k_pt_unique;
	vector < vector<int> > hkl_unique;
	if (shrink) {
		k_pt_unique.resize(3);
		hkl_unique.resize(3);
		for (int s = 0; s < sym[0][0].size(); s++)
			for (int ref = 0; ref < hkl[0].size(); ref++)
				for (int h = 0; h < 3; h++)
					hkl_unique[h].push_back(hkl[0][ref] * sym[h][0][s] + hkl[1][ref] * sym[h][1][s] + hkl[2][ref] * sym[h][2][s]);
	}

	for (int s = 0; s < sym[0][0].size(); s++) {
#pragma omp parallel for
		for (int ref = 0; ref < hkl[0].size(); ref++) 
			for (int x = 0; x < 3; x++)
				for (int h = 0; h < 3; h++) {
					double rcm_sym = 0.0;
					for (int j = 0; j < 3; j++)
						rcm_sym += unit_cell.get_rcm(x,j) * sym[j][h][s];
					k_pt[x][ref + s * hkl[0].size()] += rcm_sym * hkl[h][ref];
				}
	}

	
	if (shrink) {
		file << "Number of k-points from reflections: " << k_pt[0].size() << endl;
		file << "Determining unique k-points... ";
		file.flush();
		vector <bool> mask;
		mask.resize(k_pt[0].size());
		mask.assign(k_pt[0].size(), true);
		for (int i = 0; i < k_pt[0].size(); i++)
#pragma omp parallel for
			for (int j = i + 1; j < k_pt[0].size(); j++) {
				if (!mask[j])
					continue;
				if (k_pt[0][i] == k_pt[0][j] && k_pt[1][i] == k_pt[1][j] && k_pt[2][i] == k_pt[2][j])
					mask[j] = false;
				else if (k_pt[0][i] == -k_pt[0][j] && k_pt[1][i] == -k_pt[1][j] && k_pt[2][i] == -k_pt[2][j]) {
					mask[j] = false;
				}
			}
		if (debug)  file << "Mask done!" << endl;
		for (int i = k_pt[0].size() - 1; i >= 0; i--) {
			//if (debug) file << "i: " << i << " mask; " << mask[i] << endl;
			if (mask[i])
				for (int h = 0; h < 3; h++)
					k_pt_unique[h].insert(k_pt_unique[h].begin(), k_pt[h][i]);
			else
				for (int h = 0; h < 3; h++)
					hkl_unique[h].erase(hkl_unique[h].begin() + i);
		}
		file << " ... Done!";
		file << endl << "Number of k-points to evaluate: " << k_pt_unique[0].size() << " for " << points << " gridpoints." << endl;
	}
	else
		file << endl << "Number of k-points to evaluate: " << k_pt[0].size() << " for " << points << " gridpoints." << endl;

	vector< vector < complex<double> > > sf;
	sf.resize(asym_atom_list.size());
	if (shrink)
#pragma omp parallel for
		for (int i = 0; i < asym_atom_list.size(); i++)
			sf[i].resize(hkl_unique[0].size());
	else
#pragma omp parallel for
		for (int i = 0; i < asym_atom_list.size(); i++)
			sf[i].resize(sym[0][0].size() * hkl[0].size());

	if (debug)
		file << "Initialized FFs" << endl
		<< "asym atom list size: " << asym_atom_list.size() << " total grid size: " << points << endl;

#ifdef _WIN64
	time_t end1 = time(NULL);

	if (end1 - start < 60) file << "Time to prepare: " << fixed << setprecision(0) << end1 - start << " s" << endl;
	else if (end1 - start < 3600) file << "Time to prepare: " << fixed << setprecision(0) << floor((end1 - start) / 60) << " m " << (end1 - start) % 60 << " s" << endl;
	else file << "Time to prepare: " << fixed << setprecision(0) << floor((end1 - start) / 3600) << " h " << ((end1 - start) % 3600) / 60 << " m" << endl;
#else
	gettimeofday(&t2, 0);

	double timea = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000000.0;

	if (timea < 60) printf("Time to prepare: %4.1lf s\n", timea);
	else if (timea < 3600) printf("Time to prepare: %10.1lf m\n", timea / 60);
	else printf("Time to prepare: %10.1lf h\n", timea / 3600);

#endif

//#ifdef FLO_CUDA
//	double** gpu_k_pt = NULL,
//		** gpu_sf_r = NULL,
//		** gpu_sf_i = NULL;
//	vector <double> long_kpt;
//	long_kpt.resize(3 * k_pt_unique[0].size());
//	for (int i = 0; i < k_pt_unique[0].size(); i++) {
//		long_kpt[3 * i + 0] = k_pt_unique[0][i];
//		long_kpt[3 * i + 1] = k_pt_unique[1][i];
//		long_kpt[3 * i + 2] = k_pt_unique[2][i];
//	}
//	gpu_k_pt = (double**)malloc(sizeof(double*));
//	gpu_sf_r = (double**)malloc(sizeof(double*));
//	gpu_sf_i = (double**)malloc(sizeof(double*));
//	cudaMalloc((void**)&gpu_k_pt[0], sizeof(double) * k_pt_unique[0].size() * 3);
//	cudaMalloc((void**)&gpu_sf_r[0][i], sizeof(double) * asym_atom_list.size() * k_pt_unique[0].size());
//	cudaMalloc((void**)&gpu_sf_i[0][i], sizeof(double) * asym_atom_list.size() * k_pt_unique[0].size());
//	cudaMemcpy(gpu_k_pt[0], long_kpt.data(), sizeof(double) * k_pt_unique[0].size() * 3, cudaMemcpyHostToDevice);
//
//	dim3 blocks(asym_atom_list.size(), k_pt_unique[0].size());
//	gpu_make_sf <<<blocks, 1>>> (
//		gpu_sf_r[0],
//		gpu_sf_i[0],
//		gpu_k_pt[0],
//
//		);
//#else

	progress_bar * progress = new progress_bar{ file, 60u, "Calculating scattering factors" };
	const int step = max(floor(asym_atom_list.size() / 20),1.0);
	const int smax = shrink ? k_pt_unique[0].size() : k_pt[0].size();
	const int imax = asym_atom_list.size();
	int pmax;
	double* dens_local, * d1_local, * d2_local, * d3_local;
	complex<double>* sf_local;
	const double* k1_local = shrink ? k_pt_unique[0].data() : k_pt[0].data();
	const double* k2_local = shrink ? k_pt_unique[1].data() : k_pt[1].data();
	const double* k3_local = shrink ? k_pt_unique[2].data() : k_pt[2].data();
	double work, rho;
	for (int i = 0; i < imax; i++) {
		pmax = dens[i].size();
		dens_local = dens[i].data();
		d1_local = d1[i].data();
		d2_local = d2[i].data();
		d3_local = d3[i].data();
		sf_local = sf[i].data();
#pragma omp parallel for private(work,rho)
		for (int s = 0; s < smax; s++)
			for (int p = pmax-1; p >= 0; p--) {
				rho = dens_local[p];
				work = k1_local[s] * d1_local[p] + k2_local[s] * d2_local[p] + k3_local[s] * d3_local[p];
				sf_local[s] += complex<double>(rho * cos(work), rho * sin(work));
			}
		if (i != 0 && i % step == 0)
			progress->write(i / double(imax));
	}
	delete(progress);

	if (electron_diffraction) {
		const double fact = 0.023934;
		double h2;
#pragma omp parallel for private(h2)
		for (int s = 0; s < smax; s++) {
			h2 = pow(unit_cell.get_stl_of_hkl({ hkl_unique[0][s],hkl_unique[1][s],hkl_unique[2][s] }, file, debug),2);
			for (int i = 0; i < imax; i++)
				sf[i][s] = complex<double>(fact * (wave.atoms[asym_atom_list[i]].charge - sf[i][s].real()) / h2, -fact * sf[i][s].imag() / h2);
		}
	}

//#endif

#ifdef _WIN64

#else
	gettimeofday(&t2, 0);

	double time4 = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000000.0;

	if (time4 < 60) printf("Time to calculate: %4.1lf s\n", time4);
	else if (time4 < 3600) printf("Time to calculate: %10.1lf m\n", time4 / 60);
	else printf("Time to calculate: %10.1lf h\n", time4 / 3600);

#endif

	if (debug)
		file << endl << "SFs are made, now just write them!" << endl;
	else
		file << endl << "Writing tsc file..." << endl;

	ofstream tsc_file("experimental.tsc", ios::out);

	tsc_file << "TITLE: " << get_filename_from_path(cif).substr(0,cif.find(".cif")) << endl << "SYMM: ";

	if (shrink)
		tsc_file << "expanded";
	else
		for (int s = 0; s < sym[0][0].size(); s++) {
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++) {
					if (i != 0 || j != 0) tsc_file << " ";
					tsc_file << sym[j][i][s];
				}
			if (s != sym[0][0].size() - 1)
				tsc_file << ";";
		}

	tsc_file << endl << "SCATTERERS:";
	for (int i = 0; i < asym_atom_list.size(); i++)
		tsc_file << " " << labels[i];
	tsc_file << endl << "DATA:" << endl;

	if (shrink) {
		for (int r = 0; r < hkl_unique[0].size(); r++) {
			for (int h = 0; h < 3; h++)
				tsc_file << hkl_unique[h][r] << " ";
			for (int i = 0; i < asym_atom_list.size(); i++)
				tsc_file << scientific << setprecision(8) << real(sf[i][r]) << ","
				<< scientific << setprecision(8) << imag(sf[i][r]) << " ";
			tsc_file << endl;
		}
		//for (int r = 0; r < hkl_unique[0].size(); r++) {
		//	for (int h = 0; h < 3; h++)
		//		tsc_file << -hkl_unique[h][r] << " ";
		//	for (int i = 0; i < asym_atom_list.size(); i++)
		//		tsc_file << scientific << setprecision(8) << real(sf[i][r]) << ","
		//		<< scientific << setprecision(8) << -imag(sf[i][r]) << " ";
		//	tsc_file << endl;
		//}
	}

	else {
		if (debug) file << "Writing a total of " << sym[0][0].size() << " symmetry generated sets of hkls!" << endl;
		for (int s = 0; s < sym[0][0].size(); s++) {
			//if (debug) file << "writing symmetry: " << s << endl;
			for (int p = 0; p < hkl[0].size(); p++) {
				tsc_file
					<< hkl[0][p] * sym[0][0][s] + hkl[1][p] * sym[0][1][s] + hkl[2][p] * sym[0][2][s] << " "
					<< hkl[0][p] * sym[1][0][s] + hkl[1][p] * sym[1][1][s] + hkl[2][p] * sym[1][2][s] << " "
					<< hkl[0][p] * sym[2][0][s] + hkl[1][p] * sym[2][1][s] + hkl[2][p] * sym[2][2][s] << " ";
				for (int i = 0; i < asym_atom_list.size(); i++)
					tsc_file << scientific << setprecision(8) << real(sf[i][p + hkl[0].size() * s]) << ","
					<< scientific << setprecision(8) << imag(sf[i][p + hkl[0].size() * s]) << " ";
				tsc_file << endl;
			}
		}
	}
	tsc_file.close();

#ifdef _WIN64
	time_t end = time(NULL);

	if (end - start < 60) file << "Total Time: "<< fixed << setprecision(0) << end - start << " s\n";
	else if (end - start < 3600) file << "Total Time: " << fixed << setprecision(0) << floor((end - start)/60) << " m " << (end - start) % 60 << " s\n";
	else file << "Total Time: " << fixed << setprecision(0) << floor((end - start)/3600) << " h " << ((end - start) % 3600)/60 << " m\n";
	file << endl;
	file << "Time Breakdown:" << endl;
	file << " ... for Prototype Grid setup:" << setw(6) << end_prototypes - start << " s" << endl;
	file << " ... for Becke Grid setup:    " << setw(6) << end_becke - end_prototypes << " s" << endl;
	file << " ... for spherical density:   " << setw(6) << end_spherical - end_becke << " s" << endl;
	file << " ... for Grid Pruning:        " << setw(6) << end_prune - end_spherical << " s" << endl;
	file << " ... for aspherical density:  " << setw(6) << end_aspherical - end_prune << " s" << endl;
	file << " ... for final preparation:   " << setw(6) << end1 - end_aspherical << " s" << endl;
	file << " ... for tsc calculation:     " << setw(6) << end - end1 << " s" << endl;
#else
	gettimeofday(&t2, 0);

	double time2 = (1000000.0*(t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000000.0;

	if (time2<60) printf("Total Time: %4.1lf s\n", time2);
	else if (time2 < 3600) printf("Total Time: %10.1lf m\n", time2 / 60);
	else printf("Total Time: %10.1lf h\n", time2 / 3600);

#endif

#ifdef PEOJECT_NAME
#undef FLO_CUDA
#endif

	return true;
}

//bool calculate_structure_factors_RI(
//	string& hkl_filename,
//	string& cif,
//	string& asym_cif,
//	string& symm,
//	WFN& wave,
//	bool debug,
//	int accuracy,
//	ofstream& file,
//	int cpus,
//	bool electron_diffraction,
//	int pbc
//)
//{
//	if (cpus != -1) {
//		omp_set_num_threads(cpus);
//		omp_set_dynamic(0);
//	}
//	int* atom_z = new int[int(wave.get_ncen() * pow(pbc * 2 + 1, 3))];
//	double* x, * y, * z;
//	x = new double[int(wave.get_ncen() * pow(pbc * 2 + 1, 3))];
//	y = new double[int(wave.get_ncen() * pow(pbc * 2 + 1, 3))];
//	z = new double[int(wave.get_ncen() * pow(pbc * 2 + 1, 3))];
//	double* alpha_max = new double[wave.get_ncen()];
//	int* max_l = new int[wave.get_ncen()];
//	int max_l_overall = 0;
//
//#ifdef _WIN64
//	time_t start = time(NULL);
//	time_t end_becke, end_prototypes, end_spherical = 0, end_prune = 0, end_aspherical, end_tsc;
//#else
//	struct timeval t1, t2;
//
//	gettimeofday(&t1, 0);
//#endif
//
//	if (debug)
//		file << "Reading hkl now" << endl;
//
//	vector< vector <int> > hkl;
//	string line;
//	hkl.resize(3);
//	if (!exists(hkl_filename)) {
//		file << "HKL file does not exists!" << endl;
//		return false;
//	}
//	ifstream hkl_input(hkl_filename.c_str(), ios::in);
//	hkl_input.seekg(0, hkl_input.beg);
//	regex r{ R"([abcdefghijklmnopqrstuvwxyz\(\)ABCDEFGHIJKLMNOPQRSTUVW])" };
//	while (!hkl_input.eof()) {
//		getline(hkl_input, line);
//		if (hkl_input.eof())
//			break;
//		if (line.size() < 2)
//			continue;
//		cmatch result;
//		if (regex_search(line.c_str(), result, r))
//			continue;
//		//if (debug) file << "hkl: ";
//		for (int i = 0; i < 3; i++) {
//			string temp = line.substr(4 * size_t(i) + 1, 3);
//			temp.erase(remove_if(temp.begin(), temp.end(), ::isspace), temp.end());
//			hkl[i].push_back(stoi(temp));
//			//if (debug) file << setw(4) << temp;
//		}
//		//if (debug) file << endl;
//		if (hkl[0][hkl[0].size() - 1] == 0 && hkl[1][hkl[0].size() - 1] == 0 && hkl[2][hkl[0].size() - 1] == 0) {
//			if (debug) file << "popping back 0 0 0" << endl;
//			for (int i = 0; i < 3; i++)
//				hkl[i].pop_back();
//		}
//	}
//	hkl_input.close();
//	// Remove duplicate reflections
//	for (int i = 0; i < hkl[0].size(); i++)
//		for (int j = i + 1; j < hkl[0].size(); j++)
//			if (hkl[0][i] == hkl[0][j] && hkl[1][i] == hkl[1][j] && hkl[2][i] == hkl[2][j])
//				for (int x = 0; x < 3; x++)
//					hkl[x].erase(hkl[x].begin() + j);
//
//	if (debug)
//		file << "Reflections read! Nr of reflections: " << hkl[0].size() << endl;
//
//	double rcm[3][3];
//	double cm[3][3];
//	if (debug)
//		file << "starting to read cif!" << endl;
//	if (!exists(cif)) {
//		file << "CIF does not exists!" << endl;
//		return false;
//	}
//	if (!exists(asym_cif)) {
//		file << "Asym CIF does not exists!" << endl
//			<< asym_cif << endl;
//		return false;
//	}
//	ifstream cif_input(cif.c_str(), ios::in);
//	ifstream asym_cif_input(asym_cif.c_str(), ios::in);
//	vector<bool> found;
//	found.resize(7);
//	for (int k = 0; k < 7; k++)
//		found[k] = false;
//	double a = 0.0, b = 0.0, c = 0.0, v = 0.0;
//	double alpha = 0.0, beta = 0.0, gamma = 0.0;
//	vector <string> cell_keywords;
//	cell_keywords.push_back("_cell_length_a");
//	cell_keywords.push_back("_cell_length_b");
//	cell_keywords.push_back("_cell_length_c");
//	cell_keywords.push_back("_cell_angle_alpha");
//	cell_keywords.push_back("_cell_angle_beta");
//	cell_keywords.push_back("_cell_angle_gamma");
//	cell_keywords.push_back("_cell_volume");
//	if (debug)
//		file << "Starting while !.eof()" << endl;
//	while (!cif_input.eof()) {
//		if (debug)
//			file << "While line! " << setw(80) << line
//			<< setw(10) << a << found[0]
//			<< setw(10) << b << found[1]
//			<< setw(10) << c << found[2]
//			<< setw(10) << alpha << found[3]
//			<< setw(10) << beta << found[4]
//			<< setw(10) << gamma << found[5]
//			<< setw(10) << v << found[6] << endl;
//		getline(cif_input, line);
//		for (int k = 0; k < cell_keywords.size(); k++) {
//			if (line.find(cell_keywords[k]) != string::npos) {
//				switch (k) {
//				case 0:
//					a = stod(line.substr(cell_keywords[k].length(), line.find("(")));
//					break;
//				case 1:
//					b = stod(line.substr(cell_keywords[k].length(), line.find("(")));
//					break;
//				case 2:
//					c = stod(line.substr(cell_keywords[k].length(), line.find("(")));
//					break;
//				case 3:
//					alpha = stod(line.substr(cell_keywords[k].length(), line.find("(")));
//					break;
//				case 4:
//					beta = stod(line.substr(cell_keywords[k].length(), line.find("(")));
//					break;
//				case 5:
//					gamma = stod(line.substr(cell_keywords[k].length(), line.find("(")));
//					break;
//				case 6:
//					v = stod(line.substr(cell_keywords[k].length(), line.find("(")));
//					break;
//				default:
//					file << "This is weird... should never get here... aborting!" << endl;
//					return false;
//				}
//				found[k] = true;
//			}
//		}
//		if (found[0] == true && found[1] == true && found[2] == true && found[3] == true && found[4] == true && found[5] == true && found[6] == true)
//			break;
//	}
//	double ca = cos(0.0174532925199432944444444444444 * alpha);
//	double cb = cos(0.0174532925199432944444444444444 * beta);
//	double cg = cos(0.0174532925199432944444444444444 * gamma);
//	double sa = sin(0.0174532925199432944444444444444 * alpha);
//	double sb = sin(0.0174532925199432944444444444444 * beta);
//	double sg = sin(0.0174532925199432944444444444444 * gamma);
//	double V = a * b * c * sqrt(1 + 2 * ca * cb * cg - ca * ca - cb * cb - cg * cg);
//	if (V / v > 1.1 || V / v < 0.9) {
//		file << "Volume computed is more than 10% off, please check!" << endl;
//		return false;
//	}
//	double a_star = b * c * sa / V;
//	double b_star = a * c * sb / V;
//	double c_star = a * b * sg / V;
//
//	if (debug)
//		file << "Making cm and rcm" << endl
//		<< ca << " " << cb << " " << cg << " " << sa << " " << sb << " " << sg << " " << V << endl
//		<< a_star << " " << b_star << " " << c_star << endl;
//
//	cm[0][0] = a / 0.529177249;
//	cm[0][1] = sqrt(abs(a * b * cg)) / 0.529177249 * pow(-1, 1 + (cg > 0));
//	cm[0][2] = sqrt(abs(a * c * cb)) / 0.529177249 * pow(-1, 1 + (cb > 0));
//
//	cm[1][0] = sqrt(abs(a * b * cg)) / 0.529177249 * pow(-1, 1 + (cg > 0));
//	cm[1][1] = b / 0.529177249;
//	cm[1][2] = sqrt(abs(b * c * cb)) / 0.529177249 * pow(-1, 1 + (cb > 0));
//
//	cm[2][0] = sqrt(abs(a * c * cg)) / 0.529177249 * pow(-1, 1 + (cg > 0));
//	cm[2][1] = sqrt(abs(b * c * cb)) / 0.529177249 * pow(-1, 1 + (cb > 0));
//	cm[2][2] = c / 0.529177249;
//
//	rcm[0][0] = 2 * M_PI / a;
//	rcm[0][1] = 0;
//	rcm[0][2] = 0;
//
//	rcm[1][0] = 2 * M_PI * -cg / (a * sg);
//	rcm[1][1] = 2 * M_PI * 1 / (b * sg);
//	rcm[1][2] = 0;
//
//	rcm[2][0] = 2 * M_PI * b * c * (ca * cg - cb) / V / sg;
//	rcm[2][1] = 2 * M_PI * a * c * (cb * cg - ca) / V / sg;
//	rcm[2][2] = 2 * M_PI * a * b * sg / V;
//
//	for (int i = 0; i < 3; i++)
//		for (int j = 0; j < 3; j++)
//			if (abs(rcm[i][j]) < 10e-10) {
//				rcm[i][j] = 0.0;
//				//cm[i][j] = 0.0;
//			}
//			else {
//				rcm[i][j] *= 0.529177249;
//				//cm[i][j] *= 0.529177249;
//			}
//
//	if (debug) {
//		file << "RCM done, now labels and asym atoms!" << endl;
//		for (int i = 0; i < 3; ++i) {
//			for (int j = 0; j < 3; ++j)
//				file << setw(10) << fixed << rcm[i][j] / 2 / M_PI / 0.529177249 << ' ';
//			file << endl;
//		}
//		file << "CM in bohr:" << endl;
//		for (int i = 0; i < 3; ++i) {
//			for (int j = 0; j < 3; ++j)
//				file << setw(10) << fixed << cm[i][j] << ' ';
//			file << endl;
//		}
//	}
//	cif_input.clear();
//	cif_input.seekg(0, cif_input.beg);
//	vector <string> labels;
//	vector <int> atom_type_list;
//	vector <int> asym_atom_to_type_list;
//	int count_fields = 0;
//	int position_field[3] = { 0,0,0 };
//	int label_field = 1000;
//	vector <int> asym_atom_list;
//	vector <int> all_atom_list;
//	vector < bool > is_asym;
//	vector < vector <double > > positions;
//	positions.resize(wave.get_ncen());
//	is_asym.resize(wave.get_ncen());
//	for (int i = 0; i < wave.get_ncen(); i++) {
//		is_asym[i] = false;
//		positions[i].resize(3);
//	}
//	bool atoms_read = false;
//	while (!asym_cif_input.eof() && !atoms_read) {
//		getline(asym_cif_input, line);
//		//if(debug) file << "line: "<< line << endl;
//		if (line.find("loop_") != string::npos) {
//			//if(debug) file << "found loop!" << endl;
//			while (line.find("_") != string::npos) {
//				getline(asym_cif_input, line);
//				if (debug) file << "line in loop field definition: " << line << endl;
//				if (line.find("label") != string::npos)
//					label_field = count_fields;
//				else if (line.find("fract_x") != string::npos)
//					position_field[0] = count_fields;
//				else if (line.find("fract_y") != string::npos)
//					position_field[1] = count_fields;
//				else if (line.find("fract_z") != string::npos)
//					position_field[2] = count_fields;
//				else if (label_field == 1000) {
//					if (debug) file << "I don't think this is the atom block.. moving on!" << endl;
//					break;
//				}
//				count_fields++;
//			}
//			while (line.find("_") == string::npos && line.length() > 3) {
//				//if(debug) file << "Reading atom!"<< endl;
//				atoms_read = true;
//				stringstream s(line);
//				vector <string> fields;
//				fields.resize(count_fields);
//				for (int i = 0; i < count_fields; i++)
//					s >> fields[i];
//				positions[labels.size()][0] = (a * stod(fields[position_field[0]]) + b * cg * stod(fields[position_field[1]]) + c * cb * stod(fields[position_field[2]])) / 0.529177249;
//				positions[labels.size()][1] = (b * sg * stod(fields[position_field[1]]) + c * (ca - cb * cg) / sg * stod(fields[position_field[2]])) / 0.529177249;
//				positions[labels.size()][2] = V / (a * b * sg) * stod(fields[position_field[2]]) / 0.529177249;
//				bool found_this_one = false;
//				if (debug) file << "label: " << fields[label_field] << " position: " << positions[labels.size()][0] << " " << positions[labels.size()][1] << " " << positions[labels.size()][2] << endl;
//				for (int i = 0; i < wave.get_ncen(); i++) {
//					if (is_similar(positions[labels.size()][0], wave.atoms[i].x, -1)
//						&& is_similar(positions[labels.size()][1], wave.atoms[i].y, -1)
//						&& is_similar(positions[labels.size()][2], wave.atoms[i].z, -1)) {
//						if (debug) file << "WFN position: " << wave.atoms[i].x << " " << wave.atoms[i].y << " " << wave.atoms[i].z << endl
//							<< "Found an atom: " << fields[label_field] << " Corresponding to atom charge " << wave.atoms[i].charge << endl;
//						all_atom_list.push_back(i);
//						found_this_one = true;
//						break;
//					}
//				}
//				if (!found_this_one && debug)
//					file << "I DID NOT FIND THIS ATOM IN THE CIF?! WTF?!" << endl << positions[labels.size()][0] << " " << positions[labels.size()][1] << " " << positions[labels.size()][2] << endl;
//				labels.push_back(fields[label_field]);
//				getline(asym_cif_input, line);
//			}
//		}
//	}
//	if (labels.size() != wave.get_ncen()) {
//		file << "Number of atoms in labels: " << labels.size() << " and number of atoms in Wavefunction: " << wave.get_ncen() << "!" << endl << "This is not good, i will stop here!" << endl;
//		return false;
//	}
//	atoms_read = false;
//	label_field = 1000;
//	count_fields = 0;
//	while (!cif_input.eof() && !atoms_read) {
//		getline(cif_input, line);
//		//if(debug) file << "line: "<< line << endl;
//		if (line.find("loop_") != string::npos) {
//			//if(debug) file << "found loop!" << endl;
//			while (line.find("_") != string::npos) {
//				getline(cif_input, line);
//				if (debug) file << "line in loop field definition: " << line << endl;
//				if (line.find("label") != string::npos)
//					label_field = count_fields;
//				else if (label_field == 1000) {
//					if (debug) file << "I don't think this is the atom block.. moving on!" << endl;
//					break;
//				}
//				count_fields++;
//			}
//			while (line.find("_") == string::npos && line.length() > 3) {
//				atoms_read = true;
//				stringstream s(line);
//				vector <string> fields;
//				fields.resize(count_fields);
//				for (int i = 0; i < count_fields; i++)
//					s >> fields[i];
//				fields[label_field].erase(remove_if(fields[label_field].begin(), fields[label_field].end(), ::isspace), fields[label_field].end());
//				for (int atom = 0; atom < wave.get_ncen(); atom++) {
//					if (debug) file << "Comparing atoms: " << fields[label_field] << " / " << labels[atom] << endl;
//					if (fields[label_field] == labels[atom]) {
//						int nr = -1;
//						for (int i = 0; i < wave.get_ncen(); i++) {
//							if (is_similar(positions[atom][0], wave.atoms[i].x, -1)
//								&& is_similar(positions[atom][1], wave.atoms[i].y, -1)
//								&& is_similar(positions[atom][2], wave.atoms[i].z, -1)) {
//								if (debug) file << "Found an asymmetric atom: " << fields[label_field] << " Corresponding to atom charge " << wave.atoms[i].charge << endl;
//								nr = i;
//								asym_atom_list.push_back(i);
//								break;
//							}
//						}
//						is_asym[atom] = true;
//						bool already_there = false;
//						for (int i = 0; i < atom_type_list.size(); i++)
//							if (atom_type_list[i] == wave.atoms[nr].charge) {
//								already_there = true;
//								asym_atom_to_type_list.push_back(i);
//							}
//						if (already_there == false) {
//							asym_atom_to_type_list.push_back(atom_type_list.size());
//							atom_type_list.push_back(wave.atoms[nr].charge);
//						}
//						break;
//					}
//				}
//				getline(cif_input, line);
//			}
//		}
//	}
//
//	if (debug) {
//		file << "There are " << atom_type_list.size() << " types of atoms" << endl;
//		for (int i = 0; i < atom_type_list.size(); i++)
//			file << setw(4) << atom_type_list[i];
//		file << endl << "asym_atoms_to_type_list: " << endl;
//		for (int i = 0; i < asym_atom_to_type_list.size(); i++)
//			file << setw(4) << asym_atom_to_type_list[i];
//		file << endl;
//		file << "Mapping of asym atoms:" << endl;
//		for (int i = 0; i < wave.get_ncen(); i++)
//			file << setw(4) << wave.atoms[all_atom_list[i]].charge;
//		file << endl;
//		for (int i = 0; i < wave.get_ncen(); i++)
//			file << setw(4) << is_asym[i];
//		file << endl;
//	}
//
//	bool symm_read = false;
//	//Still need to read the sym matrices
//	if (symm == "") {
//		if (debug) file << "No Symmetry file specified, tyring to read from the CIF!" << endl;
//		symm_read = true;
//	}
//	else if (!exists(symm) && symm_read == true)
//		return false;
//
//	vector < vector < vector <int> > > sym;
//	sym.resize(3);
//	for (int i = 0; i < 3; i++)
//		sym[i].resize(3);
//	if (!symm_read) {
//		ifstream symm_input(symm.c_str(), ios::in);
//		string liny;
//		int temp_int;
//		while (!symm_input.eof()) {
//			getline(symm_input, liny);
//			stringstream streamy(liny);
//			for (int i = 0; i < 3; i++)
//				for (int j = 0; j < 3; j++) {
//					streamy >> temp_int;
//					sym[i][j].push_back(temp_int);
//				}
//		}
//	}
//	else {
//		cif_input.clear();
//		cif_input.seekg(0, cif_input.beg);
//		bool symm_found = false;
//		int operation_field = 200;
//		count_fields = 0;
//		while (!cif_input.eof() && !symm_found) {
//			getline(cif_input, line);
//			if (line.find("loop_") != string::npos) {
//				//if(debug) file << "found loop!" << endl;
//				while (line.find("_") != string::npos) {
//					getline(cif_input, line);
//					if (debug) file << "line in loop field definition: " << line << endl;
//					if (line.find("space_group_symop_operation_xyz") != string::npos)
//						operation_field = count_fields;
//					else if (count_fields > 2 || (operation_field == 200 && count_fields != 0)) {
//						if (debug) file << "I don't think this is the symmetry block.. moving on!" << endl;
//						count_fields = 0;
//						break;
//					}
//					count_fields++;
//				}
//				while (line.find("_") == string::npos && line.length() > 3 && count_fields != 0) {
//					if (debug) file << "Reading operation!" << line << endl;
//					symm_found = true;
//					stringstream s(line);
//					vector <string> fields;
//					fields.resize(count_fields);
//					int sym_from_cif[3][3];
//					for (int x = 0; x < 3; x++)
//						for (int y = 0; y < 3; y++)
//							sym_from_cif[x][y] = 0;
//					for (int i = 0; i < count_fields; i++)
//						s >> fields[i];
//					vector<string> vectors;
//					vectors.resize(3);
//					int column = 0;
//					for (int c = 0; c < fields[operation_field].length(); c++) {
//						if (fields[operation_field][c] != ',')
//							vectors[column].push_back(fields[operation_field][c]);
//						else column++;
//					}
//
//					for (int x = 0; x < 3; x++) {
//						if (vectors[x].find("X") != string::npos || vectors[x].find("x") != string::npos) {
//							char sign = ' ';
//							if (vectors[x].find("X") != string::npos && vectors[x].find("X") != 0)
//								sign = vectors[x].at(vectors[x].find("X") - 1);
//							else if (vectors[x].find("X") == 0)
//								sign = '+';
//							if (vectors[x].find("x") != string::npos && vectors[x].find("x") != 0)
//								sign = vectors[x].at(vectors[x].find("x") - 1);
//							else if (vectors[x].find("x") == 0)
//								sign = '+';
//							if (sign == '-')
//								sym_from_cif[x][0] = -1;
//							if (sign == '+')
//								sym_from_cif[x][0] = 1;
//						}
//						if (vectors[x].find("Y") != string::npos || vectors[x].find("y") != string::npos) {
//							char sign = ' ';
//							if (vectors[x].find("Y") != string::npos && vectors[x].find("Y") != 0)
//								sign = vectors[x].at(vectors[x].find("Y") - 1);
//							else if (vectors[x].find("Y") == 0)
//								sign = '+';
//							if (vectors[x].find("y") != string::npos && vectors[x].find("y") != 0)
//								sign = vectors[x].at(vectors[x].find("y") - 1);
//							else if (vectors[x].find("y") == 0)
//								sign = '+';
//							if (sign == '-')
//								sym_from_cif[x][1] = -1;
//							if (sign == '+')
//								sym_from_cif[x][1] = 1;
//						}
//						if (vectors[x].find("Z") != string::npos || vectors[x].find("z") != string::npos) {
//							char sign = ' ';
//							if (vectors[x].find("Z") != string::npos && vectors[x].find("Z") != 0)
//								sign = vectors[x].at(vectors[x].find("Z") - 1);
//							else if (vectors[x].find("Z") == 0)
//								sign = '+';
//							if (vectors[x].find("z") != string::npos && vectors[x].find("z") != 0)
//								sign = vectors[x].at(vectors[x].find("z") - 1);
//							else if (vectors[x].find("z") == 0)
//								sign = '+';
//							if (sign == '-')
//								sym_from_cif[x][2] = -1;
//							if (sign == '+')
//								sym_from_cif[x][2] = 1;
//						}
//					}
//					if (debug) {
//						file << "Comparing ";
//						for (int x = 0; x < 3; x++)
//							for (int y = 0; y < 3; y++)
//								file << sym_from_cif[x][y] << " ";
//						file << endl;
//					}
//					bool already_known = false;
//					for (int s = 0; s < sym[0][0].size(); s++) {
//						bool identical = true;
//						bool inverse = true;
//						for (int x = 0; x < 3; x++)
//							for (int y = 0; y < 3; y++)
//								if (sym[y][x][s] != sym_from_cif[x][y])
//									identical = false;
//						if (!identical)
//							for (int x = 0; x < 3; x++)
//								for (int y = 0; y < 3; y++)
//									if (sym[y][x][s] != sym_from_cif[x][y] * -1)
//										inverse = false;
//						/*if (debug) {
//							file << "Comparison with ";
//							for (int x = 0; x < 3; x++)
//								for (int y = 0; y < 3; y++)
//									file << sym[x][y][s] << " ";
//							file << "resulted in ";
//							file << identical << " " << inverse << endl;
//						}*/
//						if (identical || inverse) {
//							already_known = true;
//							break;
//						}
//					}
//					if (!already_known) {
//						if (debug) file << "This is a new symmetry operation!" << endl;
//						for (int x = 0; x < 3; x++)
//							for (int y = 0; y < 3; y++)
//								sym[y][x].push_back(sym_from_cif[x][y]);
//					}
//					getline(cif_input, line);
//				}
//			}
//		}
//	}
//
//	if (debug) {
//		file << "Read " << sym[0][0].size() << " symmetry elements! Size of sym: " << sym[0][0].size() << endl;
//		for (int i = 0; i < sym[0][0].size(); i++) {
//			for (int x = 0; x < 3; x++) {
//				for (int y = 0; y < 3; y++)
//					file << setw(3) << sym[y][x][i];
//				file << endl;
//			}
//			file << endl;
//		}
//	}
//
//	cif_input.close();
//
//	if (debug)
//		file << "There are " << atom_type_list.size() << " Types of atoms and " << asym_atom_to_type_list.size() << " atoms in total" << endl;
//
//	if (asym_atom_list.size() == 0) {
//		file << "0 asym atoms is imposible! something is wrong with reading the CIF!" << endl;
//		return false;
//	}
//
//	if (debug)
//		file << "made it post CIF, now make grids!" << endl;
//
//	double*** grid = new double** [asym_atom_list.size()];
//	int* num_points = new int[asym_atom_list.size()];
//	for (int i = 0; i < asym_atom_list.size(); i++)
//		grid[i] = new double* [6];
//	// GRID COORDINATES for [a][c][p] a = atom [0,ncen],
//	// c = coordinate [0=x, 1=y, 2=z, 3=atomic becke weight and molecular becke weight, 5=total spherical density],
//	// p = point in this grid
//
//#pragma omp parallel for
//	for (int i = 0; i < wave.get_ncen(); i++) {
//		atom_z[i] = wave.atoms[i].charge;
//		x[i] = wave.atoms[i].x;
//		y[i] = wave.atoms[i].y;
//		z[i] = wave.atoms[i].z;
//		//if(debug)
//		//    file << "xyz= 000 position: " << x[i] << " " << y[i] << " " << z[i] << " Charge: " << atom_z[i] << endl;
//		if (pbc != 0) {
//			int j = 0;
//			for (int pbc_x = -pbc; pbc_x < pbc + 1; pbc_x++)
//				for (int pbc_y = -pbc; pbc_y < pbc + 1; pbc_y++)
//					for (int pbc_z = -pbc; pbc_z < pbc + 1; pbc_z++) {
//						if (pbc_x == 0 && pbc_y == 0 && pbc_z == 0)
//							continue;
//						else {
//							j++;
//							atom_z[i + j * wave.get_ncen()] = wave.atoms[i].charge;
//							x[i + j * wave.get_ncen()] = wave.atoms[i].x + pbc_x * cm[0][0] + pbc_y * cm[0][1] + pbc_z * cm[0][2];
//							y[i + j * wave.get_ncen()] = wave.atoms[i].y + pbc_x * cm[1][0] + pbc_y * cm[1][1] + pbc_z * cm[1][2];
//							z[i + j * wave.get_ncen()] = wave.atoms[i].z + pbc_x * cm[2][0] + pbc_y * cm[2][1] + pbc_z * cm[2][2];
//							if (debug)
//								file << "xyz= " << pbc_x << pbc_y << pbc_z << " j = " << j << " position: " << x[i + j * wave.get_ncen()] << " " << y[i + j * wave.get_ncen()] << " " << z[i + j * wave.get_ncen()] << " Charge: " << atom_z[i + j * wave.get_ncen()] << endl;
//						}
//					}
//		}
//		alpha_max[i] = 0.0;
//		max_l[i] = 0;
//		for (int b = 0; b < wave.get_nex(); b++) {
//			if (wave.get_center(b) != i + 1)
//				continue;
//			if (wave.get_exponent(b) > alpha_max[i])
//				alpha_max[i] = wave.get_exponent(b);
//			if (wave.get_type(b) > max_l[i]) {
//				int l = wave.get_type(b);
//				if (l == 1)
//					l = 1;
//				else if (l >= 2 && l <= 4)
//					l = 2;
//				else if (l >= 5 && l <= 10)
//					l = 3;
//				else if (l >= 11 && l <= 20)
//					l = 4;
//				else if (l >= 21 && l <= 35)
//					l = 5;
//				max_l[i] = l;
//				if (l > max_l_overall)
//					max_l_overall = l;
//			}
//		}
//	}
//
//	if (debug)
//		file << "Atoms are there!" << endl;
//
//	double** alpha_min = new double* [wave.get_ncen()];
//	for (int i = 0; i < wave.get_ncen(); i++)
//		alpha_min[i] = new double[max_l_overall];
//
//#pragma omp parallel for
//	for (int i = 0; i < wave.get_ncen(); i++) {
//		for (int b = 0; b < max_l_overall; b++)
//			alpha_min[i][b] = 100000000.0;
//	}
//
//#pragma omp parallel for
//	for (int i = 0; i < wave.get_ncen(); i++) {
//		for (int b = 0; b < wave.get_nex(); b++) {
//			if (wave.get_center(b) != i + 1)
//				continue;
//			int l = wave.get_type(b);
//			if (l == 1)
//				l = 1;
//			else if (l >= 2 && l <= 4)
//				l = 2;
//			else if (l >= 5 && l <= 10)
//				l = 3;
//			else if (l >= 11 && l <= 20)
//				l = 4;
//			else if (l >= 21 && l <= 35)
//				l = 5;
//			else if (l >= 36 && l <= 56)
//				l = 6;
//			if (wave.get_exponent(b) < alpha_min[i][l - 1])
//				alpha_min[i][l - 1] = wave.get_exponent(b);
//		}
//		if (debug) {
//			file << "alpha_min: ";
//			for (int b = 0; b < max_l_overall; b++)
//				file << setw(14) << scientific << alpha_min[i][b];
//			file << endl;
//		}
//	}
//
//	if (debug)
//		file << "alpha_min is there!" << endl
//		<< "Nr of asym atoms: " << asym_atom_list.size() << " Number of all atoms: " << all_atom_list.size() << " Number of atoms in wfn: " << wave.get_ncen() << endl;
//	else
//		file << "There are:\n" << setw(4) << wave.get_ncen() << " atoms read from the wavefunction, of which \n"
//		//<< setw(4) << all_atom_list.size() << " will be used for grid setup and\n"
//		<< setw(4) << asym_atom_list.size() << " are identified as asymmetric unit atoms!" << endl;
//	// Total grid as a sum of all atomic grids.
//	// Dimensions: [c] [p]
//	// p = the number of gridpoint
//	// c = coordinate, which is 0=x, 1=y, 2=z, 3=atomic becke weight, 4=spherical density, 5=wavefunction density, 6=molecular becke weight
//	vector < vector <  double > > total_grid;
//
//	if (debug)
//		file << "Making Becke Grid for" << endl;
//	else
//		file << endl << "Making Becke Grids..." << flush;
//
//	total_grid.resize(6);
//
//	//Make Prototype grids with only single atom weights for all elements
//	vector <AtomGrid> Prototype_grids;
//	if (debug) file << "max_l_overall: " << max_l_overall << endl;
//	for (int i = 0; i < atom_type_list.size(); i++) {
//		if (debug) file << "Atom " << i << endl;
//		double alpha_max_temp;
//		double max_l_temp;
//		double* alpha_min_temp = new double[max_l_overall];
//		for (int j = 0; j < wave.get_ncen(); j++) {
//			if (wave.atoms[j].charge == atom_type_list[i]) {
//				if (debug) {
//					file << alpha_max[j] << " " << max_l[j] - 1 << " ";
//					for (int l = 0; l < max_l_overall; l++)
//						file << alpha_min[j][l] << " ";
//					file << endl;
//				}
//				alpha_max_temp = alpha_max[j];
//				max_l_temp = max_l[j] - 1;
//				for (int l = 0; l <= max_l_temp; l++)
//					alpha_min_temp[l] = alpha_min[j][l];
//				break;
//			}
//		}
//
//		if (debug) {
//			file << "max_l: " << defaultfloat << max_l_temp << " alpha_max: " << scientific << alpha_max_temp << " alpha_min: ";
//			for (int l = 0; l <= max_l_temp; l++)
//				file << setw(14) << scientific << alpha_min_temp[l];
//			file << " accuracy: " << accuracy << endl;
//		}
//		int lebedev_high, lebedev_low;
//		double radial_acc;
//		if (accuracy == 0) {
//			lebedev_high = (max_l_temp < 3) ? 0 : 10;
//			lebedev_low = (max_l_temp < 3) ? 0 : 10;
//			radial_acc = 1e-5;
//		}
//		else if (accuracy == 1) {
//			lebedev_high = (max_l_temp < 3) ? 110 : 146;
//			lebedev_low = (max_l_temp < 3) ? 38 : 50;
//			radial_acc = 1e-8;
//		}
//		else if (accuracy == 2) {
//			lebedev_high = (max_l_temp < 3) ? 230 : 266;
//			lebedev_low = (max_l_temp < 3) ? 110 : 146;
//			radial_acc = 1e-10;
//		}
//		else if (accuracy == 3) {
//			lebedev_high = (max_l_temp < 3) ? 350 : 590;
//			lebedev_low = (max_l_temp < 3) ? 266 : 350;
//			radial_acc = 1e-15;
//		}
//		else if (accuracy > 3) {
//			lebedev_high = (max_l_temp < 3) ? 590 : 5810;
//			lebedev_low = (max_l_temp < 3) ? 350 : 4802;
//			radial_acc = 1e-20;
//		}
//		Prototype_grids.push_back(AtomGrid(radial_acc,
//			lebedev_low,
//			lebedev_high,
//			atom_type_list[i],
//			alpha_max_temp,
//			max_l_temp,
//			alpha_min_temp,
//			debug,
//			file));
//
//	}
//
//#ifdef _WIN64
//	end_prototypes = time(NULL);
//	if (debug) {
//
//		for (int prototype = 0; prototype < Prototype_grids.size(); prototype++)
//			file << "Number of gridpoints for atom type " << atom_type_list[prototype] << " :" << Prototype_grids[prototype].get_num_grid_points() << endl;
//
//		//	int diff = end - start;
//		if (end_prototypes - start < 1) file << "Time until prototypes are done: <1 s" << endl;
//		else if (end_prototypes - start < 60) file << "Time until prototypes are done: " << fixed << setprecision(0) << end_prototypes - start << " s" << endl;
//		else if (end_prototypes - start < 3600) file << "Time until prototypes are done: " << fixed << setprecision(0) << floor((end_prototypes - start) / 60) << " m " << (end_prototypes - start) % 60 << " s" << endl;
//		else file << "Time until prototypes are done: " << fixed << setprecision(0) << floor((end_prototypes - start) / 3600) << " h " << ((end_prototypes - start) % 3600) / 60 << " m" << endl;
//	}
//#else
//	if (debug) {
//		gettimeofday(&t2, 0);
//
//		double time3 = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000000.0;
//
//		if (time3 < 60) printf("Time to prepare: %4.1lf s\n", time3);
//		else if (time3 < 3600) printf("Time to prepare: %10.1lf m\n", time3 / 60);
//		else printf("Time to prepare: %10.1lf h\n", time3 / 3600);
//	}
//#endif
//
//#pragma omp parallel for schedule(dynamic)
//	for (int i = 0; i < asym_atom_list.size(); i++) {
//		int nr = asym_atom_list[i];
//		int type;
//		for (int i = 0; i < atom_type_list.size(); i++)
//			if (atom_type_list[i] == wave.atoms[nr].charge)
//				type = i;
//
//		num_points[i] = Prototype_grids[type].get_num_grid_points();
//
//		if (debug) file << "Number of points for atom " << i << ": " << num_points[i] << endl;
//
//		for (int n = 0; n < 6; n++)
//			grid[i][n] = new double[num_points[i]];
//		for (int p = 0; p < num_points[i]; p++)
//			grid[i][4][p] = 0.0;
//
//		Prototype_grids[type].get_grid(wave.get_ncen() * pow(pbc * 2 + 1, 3),
//			nr,
//			x,
//			y,
//			z,
//			atom_z,
//			grid[i][0],
//			grid[i][1],
//			grid[i][2],
//			grid[i][3],
//			grid[i][5]);
//	}
//
//	Prototype_grids.clear();
//	int points = 0;
//	for (int i = 0; i < asym_atom_list.size(); i++)
//		points += num_points[i];
//	if (debug) file << "Becke Grid exists" << endl;
//	else file << "                                 done! Number of gridpoints: " << defaultfloat << points << endl;
//
//#ifdef _WIN64
//	if (debug) {
//		file << "Taking time..." << endl;
//	}
//	end_becke = time(NULL);
//
//#else
//	if (debug) {
//		gettimeofday(&t2, 0);
//
//		double timeb = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000000.0;
//
//		if (timeb < 60) printf("Time to prepare: %4.1lf s\n", timeb);
//		else if (timeb < 3600) printf("Time to prepare: %10.1lf m\n", timeb / 60);
//		else printf("Time to prepare: %10.1lf h\n", timeb / 3600);
//	}
//#endif
//
//	if (debug) {
//		file << "Asym atom list: ";
//		for (int g = 0; g < asym_atom_list.size(); g++)
//			file << setw(4) << asym_atom_list[g];
//		file << endl;
//		file << "All atom list: ";
//		for (int g = 0; g < all_atom_list.size(); g++)
//			file << setw(4) << all_atom_list[g];
//		file << endl;
//	}
//	for (int i = 0; i < asym_atom_list.size(); i++) {
//		int reduction = 0;
//		for (int p = 0; p < num_points[i]; p++) {
//			for (int k = 0; k < 4; k++)
//				total_grid[k].push_back(grid[i][k][p]);
//			total_grid[5].push_back(grid[i][5][p]);
//
//		}
//		num_points[i] -= reduction;
//		delete[](grid[i]);
//		if (debug) file << endl << "number of points for atom " << i << " " << num_points[i]<< endl;
//	}
//
//
//	file << "Calculating aspherical densities..." << flush;
//	vector < vector < double > > periodic_grid;
//
//	total_grid[4].resize(total_grid[0].size());
//#pragma omp parallel for
//	for (int i = 0; i < total_grid[0].size(); i++) {
//		total_grid[4][i] = compute_dens(wave, new double[3]{
//			total_grid[0][i],
//			total_grid[1][i],
//			total_grid[2][i] });
//	}
//	if (pbc != 0) {
//		periodic_grid.resize(pow(pbc * 2 + 1, 3));
//		int j = 0;
//		for (int d = 0; d < pow(pbc * 2 + 1, 3); d++)
//			periodic_grid[d].resize(total_grid[5].size());
//		for (int x = -pbc; x < pbc + 1; x++)
//			for (int y = -pbc; y < pbc + 1; y++)
//				for (int z = -pbc; z < pbc + 1; z++) {
//					if (x == 0 && y == 0 && z == 0)
//						continue;
//#pragma omp parallel for
//					for (int i = 0; i < total_grid[0].size(); i++)
//						periodic_grid[j][i] = compute_dens(wave, new double[3]{
//							total_grid[0][i] + x * cm[0][0] + y * cm[0][1] + z * cm[0][2],
//							total_grid[1][i] + x * cm[1][0] + y * cm[1][1] + z * cm[1][2],
//							total_grid[2][i] + x * cm[2][0] + y * cm[2][1] + z * cm[2][2] });
//					j++;
//				}
//		if (debug) {
//			for (int i = 0; i < total_grid[0].size(); i++) {
//				if (i % 1000 == 0)
//					file << "Old dens: " << total_grid[5][i] << " contributions of neighbour-cells:";
//				for (int j = 0; j < pow(pbc * 2 + 1, 3) - 1; j++) {
//					if (i % 1000 == 0)
//						file << " " << periodic_grid[j][i];
//					total_grid[4][i] += periodic_grid[j][i];
//				}
//				if (i % 1000 == 0)
//					file << endl;
//			}
//		}
//	}
//
//	if (debug) file << endl << "with total number of points: " << total_grid[0].size() << endl;
//	else file << "                   done!" << endl;
//
//	file << "Applying hirshfeld weights and integrating charges..." << flush;
//	double el_sum_becke = 0.0;
//	// Vector containing integrated numbers of electrons
//	// dimension 0: 0=Becek grid integration 1=Summed spherical density 2= hirshfeld weighted density
//	// dimension 1: atoms of asym_atom_list
//	vector < vector <double> > atom_els;
//	atom_els.resize(3);
//	for (int i = 0; i < asym_atom_list.size(); i++)
//		for (int n = 0; n < 3; n++)
//			atom_els[n].push_back(0.0);
//
//	//Generate Electron sums
//	for (int i = 0; i < asym_atom_list.size(); i++) {
//		int start_p = 0;
//#pragma loop(no_vector)
//		for (int a = 0; a < i; a++)
//			start_p += num_points[a];
//		for (int p = start_p; p < start_p + num_points[i]; p++) {
//			if (abs(total_grid[6][p]) > 1E-20) {
//				atom_els[0][i] += total_grid[5][p] * total_grid[4][p];
//			}
//		}
//		el_sum_becke += atom_els[0][i];
//	}
//
//	if (debug)
//		file << "Becke grid done!" << endl;
//
//#ifdef _WIN64
//	if (debug) {
//		file << "Taking time..." << endl;
//	}
//	end_aspherical = time(NULL);
//
//#else
//	if (debug) {
//		gettimeofday(&t2, 0);
//
//		double timeb = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000000.0;
//
//		if (timeb < 60) printf("Time to prepare: %4.1lf s\n", timeb);
//		else if (timeb < 3600) printf("Time to prepare: %10.1lf m\n", timeb / 60);
//		else printf("Time to prepare: %10.1lf h\n", timeb / 3600);
//	}
//#endif
//
//#pragma omp parallel for
//	for (int p = 0; p < total_grid[0].size(); p++)
//		total_grid[5][p] *= total_grid[3][p];
//
//	file << " done!" << endl;
//	file << "Number of points evaluated: " << total_grid[0].size() << " with " << el_sum_becke << " electrons in Becke Grid in total." << endl << endl;
//
//	file << "Table of Charges in electrons" << endl << endl << "Atom       Becke   Spherical Hirshfeld" << endl;
//
//	int counter = 0;
//	for (int i = 0; i < wave.get_ncen(); i++) {
//		if (is_asym[i]) {
//			file << setw(6) << labels[i]
//				<< fixed << setw(10) << setprecision(3) << wave.atoms[all_atom_list[i]].charge - atom_els[0][counter];
//			if (debug) file << " " << wave.atoms[all_atom_list[i]].charge << " " << fixed << setw(10) << setprecision(3) << wave.atoms[all_atom_list[i]].charge - atom_els[0][counter];
//			counter++;
//			file << endl;
//		}
//	}
//
//	file << "Total number of electrons in the wavefunction: " << el_sum_becke << endl;
//
//	vector < vector <double> > d1, d2, d3;
//	vector < vector <double> > dens;
//	points = 0;
//	dens.resize(asym_atom_list.size());
//	d1.resize(asym_atom_list.size());
//	d2.resize(asym_atom_list.size());
//	d3.resize(asym_atom_list.size());
//	if (debug)
//		file << "resized outer d1-3" << endl;
//
//#pragma omp parallel for reduction(+:points)
//	for (int i = 0; i < asym_atom_list.size(); i++) {
//		int start_p = 0;
//#pragma loop(no_vector)
//		for (int a = 0; a < i; a++)
//			start_p += num_points[a];
//		for (int p = start_p; p < start_p + num_points[i]; p++) {
//
//			dens[i].push_back(total_grid[5][p]);
//			d1[i].push_back(total_grid[0][p] - wave.atoms[asym_atom_list[i]].x);
//			d2[i].push_back(total_grid[1][p] - wave.atoms[asym_atom_list[i]].y);
//			d3[i].push_back(total_grid[2][p] - wave.atoms[asym_atom_list[i]].z);
//			
//		}
//		points += dens[i].size();
//	}
//
//	for (int grid = 0; grid < total_grid.size(); grid++)
//		total_grid[grid].clear();
//	total_grid.clear();
//
//	vector < vector <double> > k_pt;
//	k_pt.resize(3);
//#pragma omp parallel for
//	for (int i = 0; i < 3; i++)
//		k_pt[i].resize(sym[0][0].size() * hkl[0].size());
//
//	if (debug)
//		file << "K_point_vector is here" << endl;
//
//	for (int s = 0; s < sym[0][0].size(); s++) {
//#pragma omp parallel for
//		for (int ref = 0; ref < hkl[0].size(); ref++)
//			for (int x = 0; x < 3; x++)
//				for (int h = 0; h < 3; h++) {
//					double rcm_sym = 0.0;
//					for (int j = 0; j < 3; j++)
//						rcm_sym += rcm[x][j] * sym[j][h][s];
//					k_pt[x][ref + s * hkl[0].size()] += rcm_sym * hkl[h][ref];
//				}
//	}
//
//	file << endl << "Number of k points to evaluate: " << k_pt[0].size() << " for " << points << " gridpoints." << endl;
//
//	vector< vector < complex<double> > > sf;
//	sf.resize(asym_atom_list.size());
//
//#pragma omp parallel for
//	for (int i = 0; i < asym_atom_list.size(); i++)
//		sf[i].resize(sym[0][0].size() * hkl[0].size());
//
//	if (debug)
//		file << "Initialized FFs" << endl
//		<< "asym atom list size: " << asym_atom_list.size() << " total grid size: " << points << endl;
//
//#ifdef _WIN64
//	time_t end1 = time(NULL);
//
//	if (end1 - start < 60) file << "Time to prepare: " << fixed << setprecision(0) << end1 - start << " s" << endl;
//	else if (end1 - start < 3600) file << "Time to prepare: " << fixed << setprecision(0) << floor((end1 - start) / 60) << " m " << (end1 - start) % 60 << " s" << endl;
//	else file << "Time to prepare: " << fixed << setprecision(0) << floor((end1 - start) / 3600) << " h " << ((end1 - start) % 3600) / 60 << " m" << endl;
//#else
//	gettimeofday(&t2, 0);
//
//	double timea = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000000.0;
//
//	if (timea < 60) printf("Time to prepare: %4.1lf s\n", timea);
//	else if (timea < 3600) printf("Time to prepare: %10.1lf m\n", timea / 60);
//	else printf("Time to prepare: %10.1lf h\n", timea / 3600);
//
//#endif
//
//	progress_bar* progress = new progress_bar{ file, 60u, "Calculating scattering factors" };
//	const int step = max(floor(asym_atom_list.size() / 20), 1.0);
//	int smax = k_pt[0].size();
//	int imax = (int)asym_atom_list.size();
//	for (int i = 0; i < imax; i++) {
//		int pmax = (int)dens[i].size();
//#pragma omp parallel for
//		for (int s = 0; s < smax; s++) {
//			double work;
//			double rho;
//			for (int p = pmax - 1; p >= 0; p--) {
//				rho = dens[i][p];
//				work = k_pt[0][s] * d1[i][p] + k_pt[1][s] * d2[i][p] + k_pt[2][s] * d3[i][p];
//				sf[i][s] += complex<double>(rho * cos(work), rho * sin(work));
//			}
//		}
//		if (i != 0 && i % step == 0)
//			progress->write(i / double(imax));
//	}
//	delete(progress);
//
//	if (electron_diffraction) {
//		double fact = 2 * M_PI * 9.1093837015E-31 * pow(1.602176634E-19, 2) / (pow(6.62607015E-34, 2) * 8.8541878128E-12);
//		for (unsigned int s = 0; s < sym[0][0].size(); s++) {
//#pragma omp parallel for
//			for (int p = 0; p < (int)hkl[0].size(); p++) {
//				double h2 = pow(hkl[0][p] * sym[0][0][s] + hkl[1][p] * sym[0][1][s] + hkl[2][p] * sym[0][2][s], 2)
//					+ pow(hkl[0][p] * sym[1][0][s] + hkl[1][p] * sym[1][1][s] + hkl[2][p] * sym[1][2][s], 2)
//					+ pow(hkl[0][p] * sym[2][0][s] + hkl[1][p] * sym[2][1][s] + hkl[2][p] * sym[2][2][s], 2);
//				for (unsigned int i = 0, imax = asym_atom_list.size(); i < imax; i++)
//					sf[i][p + hkl[0].size() * s] = fact * ((double)wave.atoms[asym_atom_list[i]].charge - sf[i][p + hkl[0].size() * s]) / h2;
//			}
//		}
//	}
//
//
//#ifdef _WIN64
//	end_tsc = time(NULL);
//
//#else
//	gettimeofday(&t2, 0);
//
//	double time4 = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000000.0;
//
//	if (time4 < 60) printf("Time to calculate: %4.1lf s\n", time4);
//	else if (time4 < 3600) printf("Time to calculate: %10.1lf m\n", time4 / 60);
//	else printf("Time to calculate: %10.1lf h\n", time4 / 3600);
//
//#endif
//
//	if (debug)
//		file << endl << "SFs are made, now just write them!" << endl;
//	else
//		file << endl << "Writing tsc file..." << endl;
//
//	ofstream tsc_file("experimental.tsc", ios::out);
//
//	tsc_file << "TITLE: " << get_filename_from_path(cif).substr(0, cif.find(".cif")) << endl << "SYMM: ";
//
//	for (unsigned int s = 0; s < sym[0][0].size(); s++) {
//		for (int i = 0; i < 3; i++)
//			for (int j = 0; j < 3; j++) {
//				if (i != 0 || j != 0) tsc_file << " ";
//				tsc_file << sym[j][i][s];
//			}
//		if (s != sym[0][0].size() - 1)
//			tsc_file << ";";
//	}
//	tsc_file << endl << "SCATTERERS:";
//	for (int i = 0; i < asym_atom_list.size(); i++)
//		tsc_file << " " << labels[i];
//	tsc_file << endl << "DATA:" << endl;
//
//	if (debug) file << "Writing a total of " << sym[0][0].size() << " symmetry generated sets of hkls!" << endl;
//	for (unsigned int s = 0; s < sym[0][0].size(); s++) {
//		//if (debug) file << "writing symmetry: " << s << endl;
//		for (unsigned int p = 0; p < hkl[0].size(); p++) {
//			tsc_file
//				<< hkl[0][p] * sym[0][0][s] + hkl[1][p] * sym[0][1][s] + hkl[2][p] * sym[0][2][s] << " "
//				<< hkl[0][p] * sym[1][0][s] + hkl[1][p] * sym[1][1][s] + hkl[2][p] * sym[1][2][s] << " "
//				<< hkl[0][p] * sym[2][0][s] + hkl[1][p] * sym[2][1][s] + hkl[2][p] * sym[2][2][s] << " ";
//			for (unsigned int i = 0; i < asym_atom_list.size(); i++)
//				tsc_file << scientific << setprecision(8) << real(sf[i][p + hkl[0].size() * s]) << ","
//				<< scientific << setprecision(8) << imag(sf[i][p + hkl[0].size() * s]) << " ";
//			tsc_file << endl;
//		}
//	}
//	tsc_file.close();
//
//#ifdef _WIN64
//	time_t end = time(NULL);
//
//	if (end - start < 60) file << "Total Time: " << fixed << setprecision(0) << end - start << " s\n";
//	else if (end - start < 3600) file << "Total Time: " << fixed << setprecision(0) << floor((end - start) / 60) << " m " << (end - start) % 60 << " s\n";
//	else file << "Total Time: " << fixed << setprecision(0) << floor((end - start) / 3600) << " h " << ((end - start) % 3600) / 60 << " m\n";
//	file << endl;
//	file << "Time Breakdown:" << endl;
//	file << " ... for Prototype Grid setup:" << setw(6) << end_prototypes - start << " s" << endl;
//	file << " ... for Becke Grid setup:    " << setw(6) << end_becke - end_prototypes << " s" << endl;
//	file << " ... for spherical density:   " << setw(6) << end_spherical - end_becke << " s" << endl;
//	file << " ... for Grid Pruning:        " << setw(6) << end_prune - end_spherical << " s" << endl;
//	file << " ... for aspherical density:  " << setw(6) << end_aspherical - end_prune << " s" << endl;
//	file << " ... for final preparation:   " << setw(6) << end1 - end_aspherical << " s" << endl;
//	file << " ... for tsc calculation:     " << setw(6) << end - end1 << " s" << endl;
//#else
//	gettimeofday(&t2, 0);
//
//	double time2 = (1000000.0 * (t2.tv_sec - t1.tv_sec) + t2.tv_usec - t1.tv_usec) / 1000000.0;
//
//	if (time2 < 60) printf("Total Time: %4.1lf s\n", time2);
//	else if (time2 < 3600) printf("Total Time: %10.1lf m\n", time2 / 60);
//	else printf("Total Time: %10.1lf h\n", time2 / 3600);
//
//#endif
//
//	return true;
//}
