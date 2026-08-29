#pragma once

#include <iosfwd>

//Achieved arithmetic throughput per stage, for both the CPU and the GPU paths.
//
//The point is not curiosity. Every GPU decision in this codebase has a threshold in it -
//which precision the transform uses, how large a GEMM has to be before it is worth shipping
//to the device - and each was calibrated on one machine. A consumer card with fp64 at 1/64
//rate and a V100 at 1/2 disagree about all of them. Numbers a run reports about itself are
//what let those thresholds be re-derived on the machine that matters rather than guessed.
//
//Two rules about using it, both learned the hard way:
//
//Record once per stage, from serial code, timing the whole parallel region. Recording per
//thread inside an OpenMP region sums concurrent time and reports a total larger than the
//wall clock, which makes the profile impossible rather than merely wrong.
//
//Compare like with like. The FLOP counts below are conventions, not physics - the count for
//the transform decides what a sincos is worth, and any choice is arguable. They are applied
//identically to the CPU and GPU paths, so the ratio between them is meaningful even where
//the absolute number is a convention.
namespace throughput {

//-gflops turns the report on. Recording happens regardless: it is one mutex-guarded add per
//kernel launch or per BLAS call, never per element, so it cannot show up in a measurement.
void set_enabled(bool on);
bool enabled();

//stage is a stable literal used as the row label; device separates the two rows for a stage
//that ran in both places. flops is the convention-counted work, ms the wall time of the
//whole stage.
void record(const char* stage, bool on_device, double flops, double ms);

//Rows in first-seen order, with a CPU/GPU speed ratio wherever a stage has both.
void report(std::ostream& out);

void reset();

//The counting conventions, in one place so the CPU and GPU sides cannot drift apart.

//Non-uniform DFT: per (atom, k-point, grid point) a 3-term dot product (3 mul + 2 add), a
//sincos counted as its two results, and a complex accumulate (2 mul + 2 add). Eleven. The
//sincos is the arguable part - it is one instruction on the device and a call on the host -
//so it is counted the same on both and the ratio carries the meaning.
inline double flops_ndft(double atoms_times_points, double k_points)
{
    return 11.0 * atoms_times_points * k_points;
}

//Dense matrix product, the one unambiguous count here
inline double flops_gemm(double m, double n, double k)
{
    return 2.0 * m * n * k;
}

//SALTED's equicomb: a Wigner-3j weighted contraction over
//(atom, nrad1, nrad2, ll, mu) with a complex multiply-accumulate at the centre, counted
//at 8 flops. This is the *dense* extent - the real loop is sparse, skipping runs of zero
//w3j and screened environments - so the absolute number overstates the work, by a factor
//that depends on the model. Both paths walk the same sparsity and are counted the same
//way, so the CPU/GPU ratio is unaffected; treat the absolute rate as a lower bound.
inline double flops_equicomb(double natoms, double nrad1, double nrad2, double llmax, double l21)
{
    return 8.0 * natoms * nrad1 * nrad2 * llmax * l21;
}

} //namespace throughput
