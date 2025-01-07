#pragma once
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
#include <numeric>
#include <cassert>
#include <float.h>
#include <algorithm>
#include <atomic>
#include <deque>
#include <filesystem>
#ifdef __SSE2__
#include <emmintrin.h>
#endif
#ifdef __AVX__
#include <immintrin.h>
#endif

// Here are the system specific libaries
#ifdef _WIN32
#include <direct.h>
#define GetCurrentDir _getcwd(NULL, 0)
#include <io.h>
#include "OpenBLAS.h"
#else
#define GetCurrentDir getcwd
#include <optional>
#include <unistd.h>
#include <cfloat>
#include <sys/wait.h>
#include <termios.h>
#include <cstring>
#endif

extern void* BLAS_pointer;
extern bool has_BLAS;