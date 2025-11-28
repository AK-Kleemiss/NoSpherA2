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
#else
#define GetCurrentDir getcwd
#include <optional>
#include <unistd.h>
#include <cfloat>
#include <sys/wait.h>
#include <termios.h>
#include <cstring>
#endif
