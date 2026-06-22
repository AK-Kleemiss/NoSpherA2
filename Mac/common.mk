# Shared macOS make configuration for NoSpherA2
#
# This mirrors Linux/common.mk: thin makefiles for static/exe/tests include this
# to share flags, include paths, and output directory conventions.

.DEFAULT_GOAL := all
RM := rm -f

# Absolute/relative path to repository root from the including makefile.
REPO_ROOT ?= ..

# Build configuration (Release or Debug)
CONFIG ?= Release

# Central build output directory (repo-root). All artifacts go under build/<CONFIG>/.
BUILD_DIR ?= $(REPO_ROOT)/build/$(CONFIG)

COMP ?= clang++
C_COMP ?= clang

# ---- macOS SDK and deployment target ----
SDKROOT := $(shell xcrun --sdk macosx --show-sdk-path 2>/dev/null)
SYSROOT := $(if $(SDKROOT),-isysroot $(SDKROOT),)

# Use the same deployment target as the top-level makefile env exports.
# (Mac/makefile can override DEPLOYMENT_TARGET if needed.)
DEPLOYMENT_TARGET ?= -mmacosx-version-min=15.0

# Target architecture: default to host arch
NATIVE_ARCH ?= $(shell uname -m)
ifeq ($(NATIVE_ARCH),arm64)
  ARCH := arm64
else
  ARCH := x86_64
endif

ARCHFLAGS := -arch $(ARCH)

# Make obj/PCH dirs stable across config changes but distinct across arch + deployment target.
DEPLOYMENT_TAG := $(subst -mmacosx-version-min=,macos,$(DEPLOYMENT_TARGET))

# ---- OpenMP via libomp (Homebrew) ----
HOMEBREW_PREFIX ?= $(shell brew --prefix 2>/dev/null)
LIBOMP_PREFIX ?= $(HOMEBREW_PREFIX)/opt/libomp

OMPFLAGS := -Xpreprocessor -fopenmp -I$(LIBOMP_PREFIX)/include
OMPLIBS  := -L$(LIBOMP_PREFIX)/lib -Wl,-rpath,$(LIBOMP_PREFIX)/lib -lomp

# ---- Includes ----
# Note: we keep separate per-arch installs for featomic/libcint/occ.
USER_INC_DIRS ?= $(REPO_ROOT)/Src/core $(REPO_ROOT)/Src $(REPO_ROOT)/Lib $(REPO_ROOT)/Src/occ
SYSTEM_INC_DIRS ?= \
	$(REPO_ROOT)/Lib/featomic_install_$(ARCH)/include \
	$(REPO_ROOT)/Lib/LibCint_$(ARCH)/include \
	$(REPO_ROOT)/Lib/occ_$(ARCH)/include

INCLUDES := $(addprefix -I,$(USER_INC_DIRS)) $(foreach d,$(SYSTEM_INC_DIRS),-isystem $(d))

# ---- Local libs (built/installed by top-level makefile) ----
LIBCINT ?= $(REPO_ROOT)/Lib/LibCint_$(ARCH)/lib/libcint.a
LIBRASCALINE ?= \
	$(REPO_ROOT)/Lib/featomic_install_$(ARCH)/lib/libfeatomic.a \
	$(REPO_ROOT)/Lib/featomic_install_$(ARCH)/lib/libmetatensor.a

LIBOCC ?= $(REPO_ROOT)/Lib/occ_$(ARCH)/lib/libocc.a

# ---- Compiler flags ----
CXXSTD := -std=c++20 -DACCELERATE_NEW_LAPACK
WARNFLAGS := -fmessage-length=0
DEPFLAGS := -MMD -MP

ifeq ($(CONFIG),Debug)
  OPTFLAGS := -O0 -g
else
  OPTFLAGS := -O3 -fvisibility=hidden
endif

CXXFLAGS_COMMON := $(CXXSTD) $(OPTFLAGS) $(WARNFLAGS) $(DEPFLAGS) $(SYSROOT) $(DEPLOYMENT_TARGET) $(ARCHFLAGS) $(OMPFLAGS) $(INCLUDES)

# Linker flags
LDFLAGS_COMMON := $(SYSROOT) $(DEPLOYMENT_TARGET) $(ARCHFLAGS) -framework Accelerate $(OMPLIBS) \
	-L$(REPO_ROOT)/Lib/occ_$(ARCH)/lib -ltbb -Wl,-rpath,@executable_path

# common.mk intentionally does not define SRCS/OBJS rules; wrappers own them.
