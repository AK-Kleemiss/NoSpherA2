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
#include <array>
#include <cassert>
#include <float.h>
#include <atomic>
#include <deque>
#include <filesystem>
#include <source_location>
#ifdef __SSE2__
#include <emmintrin.h>
#endif
#ifdef __AVX__
#include <immintrin.h>
#endif

// Here are the system specific libaries
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
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

extern bool has_BLAS;