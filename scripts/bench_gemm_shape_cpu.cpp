// The CPU side of the same question: how much of eval_I's 7.2 GFLOP/s per core is the
// SHAPE rather than the hardware.
//
// eval_I calls cblas_dgemm with m = n = 64 and k ~ 4475, once per tile per reflection,
// with MKL pinned to one thread inside each OpenMP worker. That is a rank-k update: the
// operands are 4.6 MB while the work is only 36.6 MFLOP, so MKL spends much of its time
// packing rather than multiplying. The AO operand is the same for every reflection, so
// batching reflections would grow n and amortise that. This prices it.
//
// Build (from a vcvars shell, adjust the env prefix):
//   cl /O2 /EHsc /I"%MKL%\include" bench_gemm_shape_cpu.cpp /link /LIBPATH:"%MKL%\lib" ^
//      mkl_intel_lp64.lib mkl_sequential.lib mkl_core.lib

#include <mkl.h>
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <vector>

// Reusing one operand pair measures the cache-hot rate. eval_I walks a different block
// and a different reflection every call, so its operands arrive cold - cycling through a
// pool larger than L3 is the representative case.
static double bench(int m, int n, int k, int pool, int threads)
{
    const size_t na = (size_t)m * k, nb = (size_t)n * k, nc = (size_t)m * n;
    std::vector<std::vector<double>> A(pool), B(pool);
    for (int p = 0; p < pool; p++) { A[p].assign(na, 1.0); B[p].assign(nb, 1.0); }
    std::vector<double> C(nc, 0.0);
    mkl_set_num_threads_local(threads);
    const int iters = std::max(pool, (int)(4e9 / (2.0 * m * n * k)));
    for (int w = 0; w < 2; w++)
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, n, k, 1.0,
            A[0].data(), k, B[0].data(), k, 0.0, C.data(), n);
    const auto t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iters; i++) {
        const int p = i % pool;
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, n, k, 1.0,
            A[p].data(), k, B[p].data(), k, 0.0, C.data(), n);
    }
    const auto t1 = std::chrono::high_resolution_clock::now();
    const double s = std::chrono::duration<double>(t1 - t0).count();
    return 2.0 * (double)m * n * k * iters / s / 1e9;
}

// eval_I runs 16 of these concurrently, one per OpenMP worker with MKL pinned to one
// thread. That is a different machine from a lone core: the operands compete for L3 and
// DRAM. Growing n raises arithmetic intensity (A is read once for more work), which is
// worth nothing when bandwidth is free and should be worth a lot when it is not.
static double bench_threaded(int m, int n, int k, int threads)
{
    const size_t na = (size_t)m * k, nb = (size_t)n * k, nc = (size_t)m * n;
    const int iters = std::max(2, (int)(2e9 / (2.0 * m * n * k)));
    double elapsed = 0.0;
    const auto t0 = std::chrono::high_resolution_clock::now();
#pragma omp parallel num_threads(threads)
    {
        std::vector<double> A(na, 1.0), B(nb, 1.0), C(nc, 0.0);
        mkl_set_num_threads_local(1);
        for (int i = 0; i < iters; i++)
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, n, k, 1.0,
                A.data(), k, B.data(), k, 0.0, C.data(), n);
    }
    const auto t1 = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration<double>(t1 - t0).count();
    return 2.0 * (double)m * n * k * iters * threads / elapsed / 1e9;
}

int main()
{
    const int m = 64, k = 4475;
    std::printf("m=%d k=%d, single-threaded MKL (as eval_I runs it)\n\n", m, k);

    std::printf("shape sweep, operands cache-hot (one buffer pair reused):\n");
    std::printf("%8s %14s %10s\n", "n", "GFLOP/s", "vs n=64");
    double base = 0.0;
    for (int n : { 64, 256, 1024, 4096 }) {
        const double g = bench(m, n, k, 1, 1);
        if (!base) base = g;
        std::printf("%8d %14.1f %9.2fx\n", n, g, g / base);
    }

    // 40 pairs at n=64 is 40 * 2 * 2.3 MB = ~183 MB, well past any L3 here
    std::printf("\nn=64, cold operands (pool of distinct buffers):\n");
    std::printf("%8s %14s\n", "pool", "GFLOP/s");
    for (int pool : { 1, 4, 16, 40 })
        std::printf("%8d %14.1f\n", pool, bench(m, 64, k, pool, 1));

    std::printf("\n16 concurrent GEMMs, as eval_I actually runs (aggregate GFLOP/s):\n");
    std::printf("%8s %14s %12s\n", "n", "aggregate", "per core");
    for (int n : { 64, 256, 1024, 4096 }) {
        const double g = bench_threaded(m, n, k, 16);
        std::printf("%8d %14.1f %12.1f\n", n, g, g / 16.0);
    }
    return 0;
}
