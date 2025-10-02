# === macOS Makefile (Accelerate + libomp, universal) ===
SHELL := /bin/bash

RM   := rm -f
COMP ?= clang++
APP  ?= NoSpherA2

# Build both architectures for universal binary
ARCHS := arm64 x86_64

# Source & include roots
SRC_DIR := ../Src
INC_DIRS := ../Src ../Lib ../Lib/featomic_install/include

# Local static libs (adjust as needed)
EXTRA_LIBS := ../Lib/featomic_install/lib/libfeatomic.a \
			  ../Lib/featomic_install/lib/libmetatensor.a

# Sources - automatically discover all .cpp files in SRC_DIR
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
ifneq ($(filter $(SRC_DIR)/auxiliary_basis.cpp,$(SRCS)), $(SRC_DIR)/auxiliary_basis.cpp)
  SRCS += $(SRC_DIR)/auxiliary_basis.cpp
endif

# Per-arch objects
OBJS_arm64  := $(patsubst $(SRC_DIR)/%.cpp, build/arm64/%.o,  $(SRCS))
OBJS_x86_64 := $(patsubst $(SRC_DIR)/%.cpp, build/x86_64/%.o, $(SRCS))

# Common compile settings
CXXSTD    := -std=c++20 -DACCELERATE_NEW_LAPACK
OPTFLAGS  := -O3 -fvisibility=hidden
WARNFLAGS := -fmessage-length=0
INCLUDES  := $(addprefix -I,$(INC_DIRS))

# macOS SDK and deployment target
SDKROOT := $(shell xcrun --sdk macosx --show-sdk-path 2>/dev/null)
SYSROOT := $(if $(SDKROOT),-isysroot $(SDKROOT),)
# Set deployment target to 13.3 to match LAPACK requirements
DEPLOYMENT_TARGET := -mmacosx-version-min=13.3

# ---- OpenMP via libomp (Homebrew) ----
# Override these if your Homebrew lives elsewhere.
LIBOMP_PREFIX_arm64  ?= /opt/homebrew/opt/libomp
LIBOMP_PREFIX_x86_64 ?= /usr/local/opt/libomp

# Compile flags to enable OpenMP with Apple clang
OMPFLAGS_arm64  := -Xpreprocessor -fopenmp -I$(LIBOMP_PREFIX_arm64)/include
OMPFLAGS_x86_64 := -Xpreprocessor -fopenmp -I$(LIBOMP_PREFIX_x86_64)/include

# Linker flags for libomp (+rpath so users don't need DYLD_LIBRARY_PATH)
OMPLIBS_arm64   := -L$(LIBOMP_PREFIX_arm64)/lib   -Wl,-rpath,$(LIBOMP_PREFIX_arm64)/lib   -lomp
OMPLIBS_x86_64  := -L$(LIBOMP_PREFIX_x86_64)/lib  -Wl,-rpath,$(LIBOMP_PREFIX_x86_64)/lib  -lomp

# Common flags
CXXFLAGS_COMMON := $(CXXSTD) $(OPTFLAGS) $(WARNFLAGS) $(INCLUDES) $(SYSROOT) $(DEPLOYMENT_TARGET)
LDFLAGS_COMMON  := $(SYSROOT) $(DEPLOYMENT_TARGET) -framework Accelerate

# Per-arch final flags
CXXFLAGS_arm64  := $(CXXFLAGS_COMMON)  -arch arm64  $(OMPFLAGS_arm64)
CXXFLAGS_x86_64 := $(CXXFLAGS_COMMON)  -arch x86_64 $(OMPFLAGS_x86_64)
LDFLAGS_arm64   := $(LDFLAGS_COMMON)   -arch arm64  $(OMPLIBS_arm64)
LDFLAGS_x86_64  := $(LDFLAGS_COMMON)   -arch x86_64 $(OMPLIBS_x86_64)

# --------------------------------------------------------------------
# Targets
# --------------------------------------------------------------------

# Default target is universal binary
all: NoSpherA2

# Universal binary (both architectures)
NoSpherA2: NoSpherA2_arm64 NoSpherA2_x86_64
	@echo "Building MAC executable"
	@echo "Start making Mac executable"
	@echo "Lipo → $@"
	@lipo -create -output $@ NoSpherA2_arm64 NoSpherA2_x86_64
	@cp $@ ../
	@echo "Done."

# Single architecture builds
NoSpherA2_arm64: $(OBJS_arm64)
	@echo "Building MAC executable for arm64"
	@echo "Link (arm64)…"
	@$(COMP) $^ $(EXTRA_LIBS) $(LDFLAGS_arm64) -o $@

NoSpherA2_x86_64: $(OBJS_x86_64)
	@echo "Building MAC executable for x86_64"
	@echo "Link (x86_64)…"
	@$(COMP) $^ $(EXTRA_LIBS) $(LDFLAGS_x86_64) -o $@

# Architecture-specific basis set converter
basis_set_converter:
	@echo "Building basis set converter for current architecture: $(shell uname -m)"
	@$(MAKE) -f makefile_basis_set_converter $(shell uname -m)

# Generate auxiliary_basis.cpp when needed
$(SRC_DIR)/auxiliary_basis.cpp: basis_set_converter
	@echo "auxiliary_basis.cpp generated"

# Special rule for auxiliary_basis.o - fix the source file reference
build/arm64/auxiliary_basis.o: $(SRC_DIR)/auxiliary_basis.cpp | build/arm64
	@echo "Compiling auxiliary_basis.o (arm64) after BasisSetConverter finishes"
	@$(COMP) $(CXXFLAGS_arm64) -c $< -o $@
	@echo "Finished building: auxiliary_basis.cpp (arm64)"

build/x86_64/auxiliary_basis.o: $(SRC_DIR)/auxiliary_basis.cpp | build/x86_64
	@echo "Compiling auxiliary_basis.o (x86_64) after BasisSetConverter finishes"
	@$(COMP) $(CXXFLAGS_x86_64) -c $< -o $@
	@echo "Finished building: auxiliary_basis.cpp (x86_64)"

# Pattern rules
build/arm64/%.o: $(SRC_DIR)/%.cpp | build/arm64
	@echo "Compiling (arm64): $<"
	@$(COMP) $(CXXFLAGS_arm64) -c $< -o $@

build/x86_64/%.o: $(SRC_DIR)/%.cpp | build/x86_64
	@echo "Compiling (x86_64): $<"
	@$(COMP) $(CXXFLAGS_x86_64) -c $< -o $@

build/arm64  build/x86_64:
	@mkdir -p $@

# Debug variant
debug: CXXFLAGS_COMMON := -std=c++20 -O0 -g $(WARNFLAGS) $(INCLUDES) $(SYSROOT) $(DEPLOYMENT_TARGET)
debug: LDFLAGS_COMMON  := $(SYSROOT) $(DEPLOYMENT_TARGET) -framework Accelerate
debug: clean NoSpherA2_arm64

clean:
	-$(RM) -r build NoSpherA2 NoSpherA2_arm64 NoSpherA2_x86_64 ../NoSpherA2 BasisSetConverter
	@echo "Cleaned."

.PHONY: all clean debug basis_set_converter NoSpherA2 NoSpherA2_arm64 NoSpherA2_x86_64
