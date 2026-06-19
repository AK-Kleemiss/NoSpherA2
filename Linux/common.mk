# Shared Linux make configuration for NoSpherA2
#
# This file is intended to be included by thin makefiles for:
# - building the core static library (libnosphera2.a)
# - linking the CLI executable (NoSpherA2)
# - (later) linking a shared library / DLL and tests

.DEFAULT_GOAL := all
RM := rm -f

# Absolute/relative path to the repository root from the including makefile.
# Linux/makefile includes this from Linux/, so the default is `..`.
REPO_ROOT ?= ..

# Build configuration (Release or Debug)
CONFIG ?= Release

BUILD_TESTS ?= 0

# Central build output directory (repo-root). All artifacts go under build/<CONFIG>/.
BUILD_DIR ?= $(REPO_ROOT)/build/$(CONFIG)

COMP   := g++
C_COMP := gcc

# oneAPI env (source setvars.sh before each compile/link)
SETVARS ?= $(REPO_ROOT)/Lib/MKL/setvars.sh
WITH_ONEAPI = . "$(SETVARS)" >/dev/null 2>&1 &&

# ---- Optional MKL knobs ----
MKL_THREADING ?= INTEL     # INTEL | SEQ  (INTEL = mkl_intel_thread + libiomp5)

# ---- Source & include roots ----
# Note: Linux paths are case-sensitive; the directory is `Src/core`.
SRC_CORE_DIR ?= $(REPO_ROOT)/Src/core

# Target-specific source directory (override in exe/dll/tests makefiles).
# Default keeps compatibility with the existing core build.
SRC_DIR ?= $(SRC_CORE_DIR)

# Precompiled header to use for this build unit.
# Core uses Src/core/pch.h; DLL/tests have their own pch.h.
PCH_HEADER ?= $(SRC_DIR)/pch.h
PCH_BASENAME := $(notdir $(PCH_HEADER))

# Add SRC_DIR so `-include <pch>` can find the header.
USER_INC_DIRS ?= $(SRC_DIR) $(REPO_ROOT)/Src $(REPO_ROOT)/Lib $(SRC_CORE_DIR)/occ $(REPO_ROOT)/Src/occ
SYSTEM_INC_DIRS ?= $(REPO_ROOT)/Lib/featomic_install/include $(REPO_ROOT)/Lib/LibCint/include $(REPO_ROOT)/Lib/occ/include

INCLUDES := $(addprefix -I,$(USER_INC_DIRS)) $(foreach d,$(SYSTEM_INC_DIRS),-isystem $(d))

# ---- Local libs (built/installed by top-level makefile) ----
LIBCINT ?= $(REPO_ROOT)/Lib/LibCint/lib/libcint.a
LIBRASCALINE ?= $(REPO_ROOT)/Lib/featomic_install/lib/libfeatomic.a \
			   $(REPO_ROOT)/Lib/featomic_install/lib/libmetatensor.a

# ---- Vector instructions ----
# NOS_AVX=ON/OFF overrides (CI sets OFF for max compatibility); unset = detect host.
# Must match the setting OCC was built with (Eigen ABI differs with AVX).
NOS_AVX ?= $(shell grep -q avx /proc/cpuinfo 2>/dev/null && echo ON || echo OFF)
ifeq ($(NOS_AVX),OFF)
  VEC_OPTS := -msse2 -msse3 -msse4.1 -msse4.2
else
  VEC_OPTS := -msse2 -msse3 -msse4.1 -msse4.2 -mavx -mfma
endif

# ---- Compiler flags ----
# no -fopenmp at link time; only at compile time
GCC_OPTS_RELEASE ?= -std=c++2b -O3 -c -fmessage-length=0 -fopenmp -static -MMD -MP \
					$(VEC_OPTS) -ffast-math \
          -include $(PCH_BASENAME)
GCC_OPTS_DEBUG   ?= -std=c++2b -Og -g -c -fmessage-length=0 -static -MMD -MP \
					$(VEC_OPTS) -ffast-math \
          -include $(PCH_BASENAME) -mfma

ifeq ($(CONFIG),Debug)
  GCC_OPTS := $(GCC_OPTS_DEBUG)
else
  GCC_OPTS := $(GCC_OPTS_RELEASE)
endif

# ---- MKL link flags ----
ifeq ($(MKL_THREADING),SEQ)
  MKL_STATIC_LIBS = \
    -Wl,--start-group \
    "$$MKLROOT/lib/intel64/libmkl_intel_lp64.a" \
    "$$MKLROOT/lib/intel64/libmkl_sequential.a" \
    "$$MKLROOT/lib/intel64/libmkl_core.a" \
    -Wl,--end-group
  IOMP_LIB =
else
  MKL_STATIC_LIBS = \
    -Wl,--start-group \
    "$$MKLROOT/lib/intel64/libmkl_intel_lp64.a" \
    "$$MKLROOT/lib/intel64/libmkl_intel_thread.a" \
    "$$MKLROOT/lib/intel64/libmkl_core.a" \
    -Wl,--end-group
  IOMP_LIB = -L"$$ONEAPI_ROOT/compiler/latest/linux/compiler/lib/intel64_lin" \
             -Wl,-rpath,"$$ONEAPI_ROOT/compiler/latest/linux/compiler/lib/intel64_lin" \
             -Wl,-rpath,'$$ORIGIN' \
             -liomp5
endif

MKL_LINK_FLAGS = $(MKL_STATIC_LIBS) $(IOMP_LIB) -lpthread -lm -ldl

# This file intentionally does not define SRCS/OBJS or any pattern rules.
# Each target makefile (static/exe/shared/tests) owns:
# - SRCS, OBJS (and any exclusions)
# - pattern rules to compile sources into OBJS
# - dependency includes for incremental rebuilds