//The global GPU entry points of a build that carries both CUDA and HIP kernels. Every kernel
//source is compiled once per backend with NOSPHERA2_GPU_BACKEND_NS naming its namespace
//(gpu_api.h); the global names the six headers declare are defined here and forward to the
//backend that has a device. Single-backend builds do not compile this file. The backend is
//chosen once and sticks, since the two sides keep separate state: CUDA first, then HIP, and
//neither gives the CUDA side, which reports no device as a CUDA-only build would.
//NOSPHERA2_GPU_BACKEND=cuda|hip in the environment overrides the probe.
#include "sf_gpu.h"
#include "itensor_gpu.h"
#include "grid_gpu.h"
#include "salted_gpu.h"
#include "aux_density_gpu.h"
#include "spherical_density_gpu.h"
#include "blas_gpu.h"

#include <cctype>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <string>

//ret, name, parameter list, argument list: every function between the headers' NOSPHERA2_GPU_API
//markers, defaults omitted. One missing here is an unresolved global at link time of a fat build.
#define NOSPHERA2_GPU_ENTRIES(F) \
	F(bool, sf_gpu_available, (), ()) \
	F(void, sf_gpu_warmup_start, (), ()) \
	F(void, sf_gpu_warmup_wait, (), ()) \
	F(const char*, sf_gpu_backend, (), ()) \
	F(int, sf_gpu_fp64_ratio, (), ()) \
	F(bool, sf_gpu_uses_fp32, (const sf_precision prec), (prec)) \
	F(bool, sf_gpu_run, (const int imax, const long long smax, const double* k1, const double* k2, const double* k3, const double* d1, const double* d2, const double* d3, const double* dens, const int* offs, const long long total_points, double* const* sf_rows, const sf_precision prec), (imax, smax, k1, k2, k3, d1, d2, d3, dens, offs, total_points, sf_rows, prec)) \
	F(bool, itensor_gpu_available, (), ()) \
	F(const char*, itensor_gpu_gemm_name, (), ()) \
	F(double, itensor_gpu_issued_flops, (), ()) \
	F(int, itensor_gpu_batch, (int num_syms), (num_syms)) \
	F(bool, itensor_gpu_init, (const itensor_gpu_layout& L, sf_precision prec, bool tensor), (L, prec, tensor)) \
	F(bool, itensor_gpu_submit, (int slot, int n_refl, int num_syms, const double* kx, const double* ky, const double* kz, const std::complex<double>* factors), (slot, n_refl, num_syms, kx, ky, kz, factors)) \
	F(bool, itensor_gpu_collect, (int slot, std::complex<double>* I_r, long long row_stride), (slot, I_r, row_stride)) \
	F(void, itensor_gpu_free, (), ()) \
	F(bool, itensor_gpu_hold, (const std::complex<float>* I, int nr, int packed), (I, nr, packed)) \
	F(bool, itensor_gpu_hold, (const std::complex<double>* I, int nr, int packed), (I, nr, packed)) \
	F(bool, itensor_gpu_held, (), ()) \
	F(bool, itensor_gpu_rows, (const double* w, const std::complex<double>* F0, std::complex<double>* Fv), (w, F0, Fv)) \
	F(bool, itensor_gpu_cols, (const std::complex<double>* pre, double* out), (pre, out)) \
	F(void, itensor_gpu_release, (), ()) \
	F(bool, eri_gpu_hold, (const double* eri, int n, int npair, const int* pa, const int* pb, const int* first, const int* idx), (eri, n, npair, pa, pb, first, idx)) \
	F(bool, eri_gpu_JK, (const double* D, double* J, double* K), (D, J, K)) \
	F(void, eri_gpu_release, (), ()) \
	F(bool, grid_gpu_available, (), ()) \
	F(const char*, grid_gpu_backend, (), ()) \
	F(void, grid_gpu_set_enabled, (bool on), (on)) \
	F(bool, grid_gpu_enabled, (), ()) \
	F(bool, grid_gpu_becke_weights, (int np, int num_centers, const int* pcen, const double* proto_x, const double* proto_y, const double* proto_z, const double* proto_w, const double* cx, const double* cy, const double* cz, const double* R_v, const double* chi, double far_away, double cutoff, double* out_x, double* out_y, double* out_z, double* out_aw, double* out_becke, double* out_tfvc), (np, num_centers, pcen, proto_x, proto_y, proto_z, proto_w, cx, cy, cz, R_v, chi, far_away, cutoff, out_x, out_y, out_z, out_aw, out_becke, out_tfvc)) \
	F(bool, salted_gpu_available, (), ()) \
	F(void, salted_gpu_clear_cache, (), ()) \
	F(bool, salted_gpu_equicomb, (const salted_gpu_problem& prob, int* empty_environments), (prob, empty_environments)) \
	F(bool, aux_density_gpu_available, (), ()) \
	F(void, aux_density_gpu_set_enabled, (bool on), (on)) \
	F(bool, aux_density_gpu_enabled, (), ()) \
	F(bool, aux_density_gpu_eval, (int n_at, const double* cx, const double* cy, const double* cz, const double* r2_max, int n_sh, const int* sh_start, const int* sh_l, const int* pr_start, const int* coef_off, int n_pr, const double* pr_exp, const double* pr_norm, int n_coef, const double* coefs, int np, const double* x, const double* y, const double* z, double* rho, double* gx, double* gy, double* gz, double* lap, double* hess), (n_at, cx, cy, cz, r2_max, n_sh, sh_start, sh_l, pr_start, coef_off, n_pr, pr_exp, pr_norm, n_coef, coefs, np, x, y, z, rho, gx, gy, gz, lap, hess)) \
	F(bool, spherical_density_gpu_eval, (int nx, int ny, int nz, const double* origin, const double* vectors, int n_at, const double* ax, const double* ay, const double* az, const int* at_tab, int n_tab, const int* tab_off, const double* r_tab, const double* rho_tab, double lincr, double start, double radius_bohr, double* out), (nx, ny, nz, origin, vectors, n_at, ax, ay, az, at_tab, n_tab, tab_off, r_tab, rho_tab, lincr, start, radius_bohr, out)) \
	F(bool, blas_gpu_available, (), ()) \
	F(void, blas_gpu_set_enabled, (bool on), (on)) \
	F(bool, blas_gpu_enabled, (), ()) \
	F(bool, blas_gpu_dgemm, (bool transA, bool transB, int m, int n, int k, double alpha, const double* A, int lda, const double* B, int ldb, double beta, double* C, int ldc), (transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc))

#define NOSPHERA2_GPU_DECLARE(ret, name, params, args) ret name params;
namespace nosphera2_cuda {
NOSPHERA2_GPU_ENTRIES(NOSPHERA2_GPU_DECLARE)
}
namespace nosphera2_hip {
NOSPHERA2_GPU_ENTRIES(NOSPHERA2_GPU_DECLARE)
}
#undef NOSPHERA2_GPU_DECLARE

namespace {

enum class gpu_backend { cuda, hip };

gpu_backend chosen_backend()
{
	static const gpu_backend chosen = [] {
		if (const char* env = std::getenv("NOSPHERA2_GPU_BACKEND")) {
			std::string v(env);
			for (char& c : v) c = (char)std::tolower((unsigned char)c);
			if (v == "cuda") return gpu_backend::cuda;
			if (v == "hip") return gpu_backend::hip;
			if (!v.empty())
				std::fprintf(stderr, "NoSpherA2: NOSPHERA2_GPU_BACKEND=%s is neither cuda nor hip; probing instead\n", env);
		}
		if (nosphera2_cuda::sf_gpu_available()) return gpu_backend::cuda;
		if (nosphera2_hip::sf_gpu_available()) return gpu_backend::hip;
		return gpu_backend::cuda;
	}();
	return chosen;
}

} //namespace

#define NOSPHERA2_GPU_FORWARD(ret, name, params, args) \
	ret name params \
	{ \
		return chosen_backend() == gpu_backend::hip ? nosphera2_hip::name args : nosphera2_cuda::name args; \
	}
NOSPHERA2_GPU_ENTRIES(NOSPHERA2_GPU_FORWARD)
#undef NOSPHERA2_GPU_FORWARD
