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
#define MDSPAN_USE_BRACKET_OPERATOR 0
#define MDSPAN_USE_PAREN_OPERATOR 1
#ifdef __CMAKE_BUILD__
#include <mdspan/mdarray.hpp>
#else
#include  "../mdspan/include/mdspan/mdarray.hpp"
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
#define GetCurrentDir getcwd
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

// Route malloc/free/new/delete through TBB's scalable allocator.
// This is needed in BOTH the exe and the DLL:
//   - DLL: ensures all in-process allocations share one heap (avoids
//     cross-heap free crashes when OCC and NoSpherA2 exchange objects).
//   - EXE: OCC uses Eigen's handmade_aligned_free() on clang-cl (because
//     __clang__ being defined prevents EIGEN_HAS_WIN32_MALLOC from being set).
//     handmade_aligned_free calls std::free(*(ptr-1)) on the original malloc
//     pointer.  OCC's internal buffer overflow can corrupt that stored pointer.
//     scalable_free() on a non-TBB pointer falls back to HeapFree via a
//     secondary path that does NOT raise STATUS_HEAP_CORRUPTION for a
//     garbage-value pointer, whereas ucrtbase's direct HeapFree does.
//     Linking tbbmalloc_proxy therefore lets run_method<DFT> return normally
//     instead of crashing before the computation result can be used.
#include <tbb/tbbmalloc_proxy.h>

// Must be last: redefines exit() as a throw when NOSPHERA2_IN_PROCESS is set.
#include "early_exit.h"

#endif // NOSPHERA2_PCH_H

