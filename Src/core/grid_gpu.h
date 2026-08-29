#pragma once

//GPU path for the atomic integration grid weights. generateIntegrationGrids is 46% of a
//production tsc run and 98.8% of it is AtomGrid::get_grid, which evaluates the Becke and
//TFVC partitioning weights at every grid point. Each point is independent and the work is
//O(centers^2) of plain arithmetic per point, which is the shape a device wants.
//
//This is partitioning physics, not a contraction: a slip changes every partitioned charge.
//The device code is a verbatim transcription of the chi-present branch of
//get_integration_weights.
//
//That branch is not the only one get_grid can reach: make_chi returns empty when the
//wavefunction carries no MOs, and the CPU has a chi-absent branch this kernel does not
//implement. The caller checks chi is exactly num_centers^2 first - make_chi lays its rows
//out with a stride of wfn.get_ncen(), which is not num_centers in general, and a wrongly
//sized chi copies without complaint and comes back wrong.
//
//Returns false if there is no device, if the scratch will not fit, or if num_centers is
//too large for the per-thread arrays, and the caller keeps the CPU loop.

bool grid_gpu_available();
//"CUDA" or "HIP", for the log line
const char* grid_gpu_backend();

//-gpu_grid turns this on; off unless asked, like the other GPU paths
void grid_gpu_set_enabled(bool on);
bool grid_gpu_enabled();

bool grid_gpu_becke_weights(
	int np, int num_centers, int center_index,
	const double* proto_x, const double* proto_y, const double* proto_z, const double* proto_w,
	const double* cx, const double* cy, const double* cz,
	const double* R_v,          //[num_centers] Bragg radii, already looked up by Z
	const double* chi,          //[num_centers * num_centers]
	double far_away, double cutoff,
	double* out_x, double* out_y, double* out_z,
	double* out_aw, double* out_becke, double* out_tfvc);
