#pragma once

#include "gpu_api.h"

//Becke and TFVC partitioning weights per grid point, O(centers^2) each; the device code is a
//verbatim transcription of get_integration_weights, keep the two in step.
//One call covers every atom: prototype points of all grids concatenated, pcen the owning centre.
//chi is null without TFVC; otherwise exactly num_centers^2, while make_chi strides by
//wfn.get_ncen(), and a wrongly sized chi comes back wrong without complaint.
//False without a device, room, or per-thread array space for num_centers; the caller keeps the CPU loop.

NOSPHERA2_GPU_API_BEGIN

bool grid_gpu_available();
//"CUDA" or "HIP"
const char* grid_gpu_backend();

//On by default, off with -no_gpu_grid
void grid_gpu_set_enabled(bool on);
bool grid_gpu_enabled();

bool grid_gpu_becke_weights(
	int np, int num_centers,
	const int* pcen,            //[np] owning centre of each point
	const double* proto_x, const double* proto_y, const double* proto_z, const double* proto_w,
	const double* cx, const double* cy, const double* cz,
	const double* R_v,          //[num_centers] Bragg radii, already looked up by Z
	const double* chi,          //[num_centers * num_centers] or null
	double far_away, double cutoff,
	double* out_x, double* out_y, double* out_z,
	double* out_aw, double* out_becke, double* out_tfvc);

NOSPHERA2_GPU_API_END
