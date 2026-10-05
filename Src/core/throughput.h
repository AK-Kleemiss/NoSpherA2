#pragma once

#include <iosfwd>

//Achieved flop rate per stage, CPU and GPU. Record from serial code around the whole parallel region (a per-thread
//timer sums concurrent time). The flop counts are conventions shared by both paths, so the ratio holds where the rate may not
namespace throughput {

//-gflops; one guarded add per launch, never per element
void set_enabled(bool on);
bool enabled();

void record(const char* stage, bool on_device, double flops, double ms);

//No honest flop count: the rate columns read "-"
void record_time(const char* stage, bool on_device, double ms);

//Rows in first-seen order, with a CPU/GPU ratio where a stage has both
void report(std::ostream& out);

void reset();

//Non-uniform DFT per (atom, k-point, grid point): 3-term dot product, sincos counted as its two results, complex accumulate
inline double flops_ndft(double atoms_times_points, double k_points)
{
	return 11.0 * atoms_times_points * k_points;
}

inline double flops_gemm(double m, double n, double k)
{
	return 2.0 * m * n * k;
}

//equicomb: Wigner-3j contraction over (atom, nrad1, nrad2, ll, mu), complex MAC = 8; dense extent of a sparse loop, a lower bound
inline double flops_equicomb(double natoms, double nrad1, double nrad2, double llmax, double l21)
{
	return 8.0 * natoms * nrad1 * nrad2 * llmax * l21;
}

} //namespace throughput
