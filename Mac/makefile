# === macOS Makefile (Accelerate + libomp, universal) ===
SHELL := /bin/bash

RM   := rm -f
COMP ?= clang++
APP  ?= NoSpherA2

# Build both architectures; set ARCHS := arm64 if you don't have x86_64 libomp installed
#ARCHS := arm64 x86_64
ARCHS := arm64

# Source & include roots
SRC_DIR := ../Src
INC_DIRS := ../Src ../Lib ../Lib/featomic_install/include

# Local static libs (adjust as needed)
EXTRA_LIBS := ../Lib/featomic_install/lib/libfeatomic.a \
              ../Lib/featomic_install/lib/libmetatensor.a

# Sources
SRCS := $(wildcard $(SRC_DIR)/*.cpp)

# Per-arch objects
OBJS_arm64  := $(patsubst $(SRC_DIR)/%.cpp, build/arm64/%.o,  $(SRCS))
OBJS_x86_64 := $(patsubst $(SRC_DIR)/%.cpp, build/x86_64/%.o, $(SRCS))

# Common compile settings
CXXSTD    := -std=c++20
OPTFLAGS  := -O3 -fvisibility=hidden
WARNFLAGS := -fmessage-length=0
INCLUDES  := $(addprefix -I,$(INC_DIRS))

# macOS SDK
SDKROOT := $(shell xcrun --sdk macosx --show-sdk-path 2>/dev/null)
SYSROOT := $(if $(SDKROOT),-isysroot $(SDKROOT),)

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
CXXFLAGS_COMMON := $(CXXSTD) $(OPTFLAGS) $(WARNFLAGS) $(INCLUDES) $(SYSROOT)
LDFLAGS_COMMON  := $(SYSROOT) -framework Accelerate

# Per-arch final flags
CXXFLAGS_arm64  := $(CXXFLAGS_COMMON)  -arch arm64  $(OMPFLAGS_arm64)
CXXFLAGS_x86_64 := $(CXXFLAGS_COMMON)  -arch x86_64 $(OMPFLAGS_x86_64)
LDFLAGS_arm64   := $(LDFLAGS_COMMON)   -arch arm64  $(OMPLIBS_arm64)
LDFLAGS_x86_64  := $(LDFLAGS_COMMON)   -arch x86_64 $(OMPLIBS_x86_64)

# --------------------------------------------------------------------
# Targets
# --------------------------------------------------------------------
all: $(APP)

# Keep your original converter step
basis_set_converter:
	@$(MAKE) -f makefile_basis_set_converter

# Ensure aux file builds after converter
build/%/auxiliary_basis.o: basis_set_converter

# Pattern rules
build/arm64/%.o: $(SRC_DIR)/%.cpp | build/arm64
	@echo "Compiling (arm64): $<"
	@$(COMP) $(CXXFLAGS_arm64) -c $< -o $@

build/x86_64/%.o: $(SRC_DIR)/%.cpp | build/x86_64
	@echo "Compiling (x86_64): $<"
	@$(COMP) $(CXXFLAGS_x86_64) -c $< -o $@

build/arm64  build/x86_64:
	@mkdir -p $@

# Link per-arch
NoSpherA2_arm64: $(OBJS_arm64)
	@echo "Link (arm64)…"
	@$(COMP) $^ $(EXTRA_LIBS) $(LDFLAGS_arm64) -o $@

NoSpherA2_x86_64: $(OBJS_x86_64)
	@echo "Link (x86_64)…"
	@$(COMP) $^ $(EXTRA_LIBS) $(LDFLAGS_x86_64) -o $@

# Create universal
$(APP): $(filter NoSpherA2_%,$(foreach A,$(ARCHS),NoSpherA2_$(A)))
	@echo "Lipo → $@"
	@lipo -create -output $@ $^
	@cp $@ ../
	@echo "Done."

# Debug variant
debug: CXXFLAGS_COMMON := -std=c++20 -O0 -g $(WARNFLAGS) $(INCLUDES) $(SYSROOT)
debug: LDFLAGS_COMMON  := $(SYSROOT) -framework Accelerate
debug: clean $(APP)

clean:
	-$(RM) -r build NoSpherA2 NoSpherA2_arm64 NoSpherA2_x86_64 ../NoSpherA2
	@echo "Cleaned."

.PHONY: all clean debug basis_set_converter
