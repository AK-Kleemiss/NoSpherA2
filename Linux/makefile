# === LINUX MAKEFILE ===
RM := rm -f

COMP   := g++
C_COMP := gcc

LIBS := -Wl,-Bstatic -lstdc++ -lpthread -Wl,-Bdynamic -lm -lc

# oneAPI env (source setvars.sh before each compile/link)
SETVARS ?= ../Lib/MKL/setvars.sh
WITH_ONEAPI = . "$(SETVARS)" >/dev/null 2>&1 &&

# ---- Optional MKL knobs ----
MKL_THREADING ?= INTEL     # INTEL | SEQ  (INTEL = mkl_intel_thread + libiomp5)

# Source & include roots
SRC_DIR := ../Src
INC_DIRS := ../Src ../Lib ../Lib/featomic_install/include

# Featomic local libs
LIBRASCALINE   := ../Lib/featomic_install/lib/libfeatomic.a
LIBRASCALINE   += ../Lib/featomic_install/lib/libmetatensor.a
LIBRASCALINE_I := ../Lib/featomic_install/include

# ---- Compiler flags (no -fopenmp at link time; only at compile time) ----
GCC_OPTS       := -std=c++2b -O3 -c -fmessage-length=0 -fopenmp -static -MMD -MP \
                  -msse2 -msse3 -msse4.1 -msse4.2 -mavx -ffast-math
GCC_OPTS_DEBUG := -std=c++2b -Og -g -c -fmessage-length=0 -static -MMD -MP \
                  -msse2 -msse3 -msse4.1 -msse4.2 -mavx -ffast-math

# Include directories
INCLUDES := $(addprefix -I,$(INC_DIRS))

# Threading backend selection - STATIC MKL but DYNAMIC libomp
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
  # Use dynamic linking for libiomp5 with rpath for runtime location
  IOMP_LIB = -L"$$ONEAPI_ROOT/compiler/latest/linux/compiler/lib/intel64_lin" \
             -Wl,-rpath,"$$ONEAPI_ROOT/compiler/latest/linux/compiler/lib/intel64_lin" \
             -Wl,-rpath,'$$ORIGIN' \
             -liomp5
endif

# MKL link flags - static MKL libraries + dynamic libiomp5
MKL_LINK_FLAGS = $(MKL_STATIC_LIBS) $(IOMP_LIB) -lpthread -lm -ldl

# ---- Objects ----
# Sources - automatically discover all .cpp files in SRC_DIR
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
ifneq ($(filter $(SRC_DIR)/auxiliary_basis.cpp,$(SRCS)), $(SRC_DIR)/auxiliary_basis.cpp)
  SRCS += $(SRC_DIR)/auxiliary_basis.cpp
endif

# Convert sources to object names (extract basename without path/extension)
OBJ_NAMES := $(basename $(notdir $(SRCS)))
OBJS   := $(OBJ_NAMES:=.o)
OBJS_d := $(OBJ_NAMES:=.o_d)

all: NoSpherA2
