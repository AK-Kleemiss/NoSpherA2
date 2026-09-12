#pragma once

#include <iosfwd>

//Achieved arithmetic throughput per stage, for the CPU and GPU paths.
//
//Record once per stage from serial code, timing the whole parallel region; a per-thread
//timer inside an omp region sums concurrent time and exceeds the wall clock. The flop
//counts below are conventions applied identically to both paths, so the ratio is
//meaningful even where the absolute number is not.
namespace throughput {

//-gflops turns the report on. Recording is one guarded add per launch, never per element.
void set_enabled(bool on);
bool enabled();

//stage labels the row; device separates the two rows for a stage that ran in both places.
void record(const char* stage, bool on_device, double flops, double ms);

//For stages with no honest flop count. The rate columns read "-" rather than inventing one.
void record_time(const char* stage, bool on_device, double ms);

//Rows in first-seen order, with a CPU/GPU ratio wherever a stage has both.
void report(std::ostream& out);

void reset();

//The counting conventions, in one place so the CPU and GPU sides cannot drift apart.

//Non-uniform DFT: per (atom, k-point, grid point) a 3-term dot product, a sincos counted
//as its two results, and a complex accumulate. Eleven.
inline double flops_ndft(double atoms_times_points, double k_points)
{
    return 11.0 * atoms_times_points * k_points;
}

//Dense matrix product, the one unambiguous count here
inline double flops_gemm(double m, double n, double k)
{
    return 2.0 * m * n * k;
}

//equicomb: Wigner-3j weighted contraction over (atom, nrad1, nrad2, ll, mu), a complex
//multiply-accumulate at 8 flops. Dense extent; the real loop is sparse, so this is a
//lower bound on the rate and both paths are counted the same way.
inline double flops_equicomb(double natoms, double nrad1, double nrad2, double llmax, double l21)
{
    return 8.0 * natoms * nrad1 * nrad2 * llmax * l21;
}

} //namespace throughput
