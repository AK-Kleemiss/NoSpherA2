#ifndef NOSPHERA2_PCH_H
#define NOSPHERA2_PCH_H

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#ifdef __cplusplus__
#include <cstdlib>
#else
#include <stdlib.h>
#endif
#include <functional>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif

#define HAVE_ECPINT 1
#include <occ/main/occ_scf.h>
#include <occ/gto/gto.h>
#include <occ/core/parallel.h>
#include <occ/qm/io/conversion.h>
#include <occ/io/occ_input.h>
#include <occ/qm/wavefunction.h>
#include <occ/io/xyz.h>
#include <occ/core/molecule.h>
#include <occ/qm/scf.h>
#include <occ/qm/scf_impl.h>
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
#if defined(__APPLE__)
// On macOS we are using Accelerate for BLAS/LAPACK
#include <Accelerate/Accelerate.h>
#define lapack_int int
#define MKL_Set_Num_Threads(num) omp_set_num_threads(num)
#elif defined(NSA2_OPENBLAS)
// Windows on ARM: oneMKL has no ARM64 build, OpenBLAS supplies CBLAS and LAPACKE
#include <cblas.h>
// MSVC C++ cannot parse the C99 _Complex default
#include <complex>
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include <lapacke.h>
#if defined(NSA2_ARMPL) && defined(_WIN32)
//serial armpl_lp64, nothing to set
#define MKL_Set_Num_Threads(num) ((void)(num))
#elif defined(NSA2_ARMPL)
//armpl_lp64_mp takes its thread count from OpenMP, as Accelerate's macro above does
#define MKL_Set_Num_Threads(num) omp_set_num_threads(num)
#elif defined(__ANDROID__)
#include <cstdio>
//the NDK OpenBLAS is built USE_OPENMP=1 and runs serial inside OpenMP regions on its own, so
//the serial SALTED GEMMs get -cpus threads, but only within cpu0's core cluster: Exynos 7870
//sucrose kernels 8.5 s on 1 thread, 3.0 s on 4, 14-16 s on 5-8 (spilling into the 2nd cluster)
inline int nsa2_cluster_threads(int num) {
  int a = 0, b = 3;
  if (FILE* f = fopen("/sys/devices/system/cpu/cpu0/topology/core_siblings_list", "r")) {
    if (fscanf(f, "%d-%d", &a, &b) != 2) b = a + 3;
    fclose(f);
  }
  return std::max(1, std::min(num, b - a + 1));
}
#define MKL_Set_Num_Threads(num) openblas_set_num_threads(nsa2_cluster_threads(num))
#else
//pthreads OpenBLAS cannot tell it runs inside an OpenMP region (MKL can), so N threads
//here would mean N BLAS threads per OpenMP thread: keep it serial, -cpus goes to OpenMP
#define MKL_Set_Num_Threads(num) ((void)(num), openblas_set_num_threads(1))
#endif
#else
// Linux/Windows with oneMKL
#include <mkl.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif
#include <regex>
#include <set>
#include <map>
#include <string>
#include <stdexcept>
#include <utility>
#include <sstream>
#include <typeinfo>
#include <vector>
#include <array>
#include <cassert>
#include <float.h>
#include <atomic>
#include <deque>
#include <filesystem>
#include <source_location>
#include <memory>
#include <cstddef>
#include <limits>
#include <cstdio>
#include <ranges>
#ifdef __SSE2__
#include <emmintrin.h>
#endif
#ifdef __AVX__
#include <immintrin.h>
#endif
#if defined(__aarch64__) || defined(_M_ARM64)
#include <arm_neon.h>
#endif
#define MDSPAN_USE_BRACKET_OPERATOR 0
#define MDSPAN_USE_PAREN_OPERATOR 1
#ifdef __CMAKE_BUILD__
#include <mdspan/mdarray.hpp>
#else
#include  "../../mdspan/include/mdspan/mdarray.hpp"
#endif


// Here are the system specific libaries
#ifdef _WIN32
//#define WIN32_LEAN_AND_MEAN
#include <direct.h>
#define GetCurrentDir _getcwd(NULL, 0)
#include <io.h>
#define NOMINMAX
#include <windows.h>
#include <shobjidl.h>
#include <algorithm>
#else
#define GetCurrentDir std::filesystem::current_path().string() // bare getcwd streamed a function pointer ("1" on glibc, an error on bionic)
#include <optional>
#include <unistd.h>
#include <cfloat>
#include <sys/wait.h>
#include <termios.h>
#include <cstring>
#if defined(__APPLE__)
#include <mach-o/dyld.h>
#endif
#endif


// Must be last: redefines exit() as a throw when NOSPHERA2_IN_PROCESS is set.
#include "early_exit.h"

#endif // NOSPHERA2_PCH_H

