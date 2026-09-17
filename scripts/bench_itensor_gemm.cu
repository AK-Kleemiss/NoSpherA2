// Throughput of the I tensor's GEMM shape on this device, in fp64 and fp32.
//
// eval_I spends 99.4% of its time in cblas_dgemm with m = n = 64 (the screened tile
// size) and k = the grid block's point count, measured at a mean of 4475. That is a
// rank-k update, not a square GEMM, so peak-flops arithmetic is not a safe guide -
// this measures the shape that is actually issued.
//
// Build: nvcc -O3 -arch=sm_89 bench_itensor_gemm.cu -lcublas -o bench_itensor_gemm

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cstdio>
#include <vector>

#define M 64
#define N 64
#define K 4475
#define BATCH 64
#define ITERS 20

#define CK(x) do { if ((x) != cudaSuccess) { std::printf("cuda error line %d\n", __LINE__); return 1; } } while (0)

template <typename T>
static double run(cublasHandle_t h, bool f32, T* A, T* B, T* C)
{
	cudaEvent_t t0, t1;
	cudaEventCreate(&t0);
	cudaEventCreate(&t1);
	const T alpha = (T)1, beta = (T)0;
	//Warm up so the timing is not the first-call cuBLAS setup
	for (int i = 0; i < 3; i++) {
		if (f32)
			cublasSgemmStridedBatched(h, CUBLAS_OP_T, CUBLAS_OP_N, M, N, K,
				(const float*)&alpha, (const float*)A, K, (long long)M * K,
				(const float*)B, K, (long long)N * K, (const float*)&beta,
				(float*)C, M, (long long)M * N, BATCH);
		else
			cublasDgemmStridedBatched(h, CUBLAS_OP_T, CUBLAS_OP_N, M, N, K,
				(const double*)&alpha, (const double*)A, K, (long long)M * K,
				(const double*)B, K, (long long)N * K, (const double*)&beta,
				(double*)C, M, (long long)M * N, BATCH);
	}
	cudaDeviceSynchronize();
	cudaEventRecord(t0);
	for (int i = 0; i < ITERS; i++) {
		if (f32)
			cublasSgemmStridedBatched(h, CUBLAS_OP_T, CUBLAS_OP_N, M, N, K,
				(const float*)&alpha, (const float*)A, K, (long long)M * K,
				(const float*)B, K, (long long)N * K, (const float*)&beta,
				(float*)C, M, (long long)M * N, BATCH);
		else
			cublasDgemmStridedBatched(h, CUBLAS_OP_T, CUBLAS_OP_N, M, N, K,
				(const double*)&alpha, (const double*)A, K, (long long)M * K,
				(const double*)B, K, (long long)N * K, (const double*)&beta,
				(double*)C, M, (long long)M * N, BATCH);
	}
	cudaEventRecord(t1);
	cudaEventSynchronize(t1);
	float ms = 0;
	cudaEventElapsedTime(&ms, t0, t1);
	const double flop = 2.0 * M * N * K * BATCH * ITERS;
	return flop / (ms * 1e-3) / 1e9;
}

int main()
{
	cublasHandle_t h;
	cublasCreate(&h);
	const size_t na = (size_t)M * K * BATCH, nb = (size_t)N * K * BATCH, nc = (size_t)M * N * BATCH;

	double *dA, *dB, *dC;
	CK(cudaMalloc(&dA, na * sizeof(double)));
	CK(cudaMalloc(&dB, nb * sizeof(double)));
	CK(cudaMalloc(&dC, nc * sizeof(double)));
	CK(cudaMemset(dA, 0x3f, na * sizeof(double)));
	CK(cudaMemset(dB, 0x3f, nb * sizeof(double)));
	const double g64 = run<double>(h, false, dA, dB, dC);
	cudaFree(dA); cudaFree(dB); cudaFree(dC);

	float *fA, *fB, *fC;
	CK(cudaMalloc(&fA, na * sizeof(float)));
	CK(cudaMalloc(&fB, nb * sizeof(float)));
	CK(cudaMalloc(&fC, nc * sizeof(float)));
	CK(cudaMemset(fA, 0x3f, na * sizeof(float)));
	CK(cudaMemset(fB, 0x3f, nb * sizeof(float)));
	const double g32 = run<float>(h, true, fA, fB, fC);
	cudaFree(fA); cudaFree(fB); cudaFree(fC);

	std::printf("m=%d n=%d k=%d batch=%d\n", M, N, K, BATCH);
	std::printf("cublas DgemmStridedBatched : %8.1f GFLOP/s\n", g64);
	std::printf("cublas SgemmStridedBatched : %8.1f GFLOP/s\n", g32);
	std::printf("fp32:fp64 measured ratio   : %8.1f\n", g32 / g64);
	cublasDestroy(h);
	return 0;
}
