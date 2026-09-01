// How much of the I tensor's GEMM cost is the SHAPE rather than the hardware.
//
// eval_I issues m = n = 64 with k ~ 4475: a rank-k update. Both MKL and cuBLAS are far
// from peak on it. The AO operand A is identical for every reflection - only the
// phase-weighted operand changes - so several reflections could share one call, which
// grows n while leaving m and k alone. This sweeps n to price that.
//
// Build: nvcc -O3 -arch=sm_89 bench_gemm_shape.cu -lcublas -o bench_gemm_shape

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cstdio>

#define M 64
#define K 4475
#define ITERS 10

static double run(cublasHandle_t h, int n, float* A, float* B, float* C)
{
	cudaEvent_t t0, t1;
	cudaEventCreate(&t0);
	cudaEventCreate(&t1);
	const float alpha = 1.0f, beta = 0.0f;
	for (int i = 0; i < 3; i++)
		cublasSgemm(h, CUBLAS_OP_T, CUBLAS_OP_N, M, n, K, &alpha, A, K, B, K, &beta, C, M);
	cudaDeviceSynchronize();
	cudaEventRecord(t0);
	for (int i = 0; i < ITERS; i++)
		cublasSgemm(h, CUBLAS_OP_T, CUBLAS_OP_N, M, n, K, &alpha, A, K, B, K, &beta, C, M);
	cudaEventRecord(t1);
	cudaEventSynchronize(t1);
	float ms = 0;
	cudaEventElapsedTime(&ms, t0, t1);
	return 2.0 * M * n * K * ITERS / (ms * 1e-3) / 1e9;
}

int main()
{
	cublasHandle_t h;
	cublasCreate(&h);
	const int nmax = 8192;
	float *A, *B, *C;
	if (cudaMalloc(&A, (size_t)M * K * sizeof(float)) != cudaSuccess) return 1;
	if (cudaMalloc(&B, (size_t)nmax * K * sizeof(float)) != cudaSuccess) return 1;
	if (cudaMalloc(&C, (size_t)M * nmax * sizeof(float)) != cudaSuccess) return 1;
	cudaMemset(A, 0x3f, (size_t)M * K * sizeof(float));
	cudaMemset(B, 0x3f, (size_t)nmax * K * sizeof(float));

	std::printf("m=%d k=%d, sweeping n (n=64 is what eval_I issues today)\n", M, K);
	std::printf("%8s %14s %10s\n", "n", "GFLOP/s", "vs n=64");
	double base = 0;
	for (int n = 64; n <= nmax; n *= 4) {
		const double g = run(h, n, A, B, C);
		if (n == 64) base = g;
		std::printf("%8d %14.1f %9.2fx\n", n, g, g / base);
	}
	cudaFree(A); cudaFree(B); cudaFree(C);
	cublasDestroy(h);
	return 0;
}
