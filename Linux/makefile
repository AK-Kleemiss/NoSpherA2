# === LINUX MAKEFILE ===
RM := rm -f

COMP   := g++
C_COMP := gcc

LIBS := -Wl,-Bstatic -lstdc++ -lpthread -Wl,-Bdynamic -lm -lc -lquadmath

# oneAPI env (source setvars.sh before each compile/link)
SETVARS ?= ../Lib/MKL/setvars.sh
WITH_ONEAPI = . "$(SETVARS)" >/dev/null 2>&1 &&

# ---- Optional MKL knobs ----
MKL_THREADING ?= INTEL     # INTEL | SEQ  (INTEL = mkl_intel_thread + libiomp5)

# Source & include roots
SRC_DIR := ../Src
INC_DIRS := ../Src ../Lib ../Lib/featomic_install/include ../Lib/LibCint/include ../Lib/occ/include
LIB_OCC := $(wildcard ../Lib/occ/lib/*.a)
#LibCint local libs
LIBCINT		   := ../Lib/LibCint/lib/libcint.a

# Featomic local libs
LIBRASCALINE   := ../Lib/featomic_install/lib/libfeatomic.a
LIBRASCALINE   += ../Lib/featomic_install/lib/libmetatensor.a
LIBRASCALINE_I := ../Lib/featomic_install/include

# ---- Compiler flags (no -fopenmp at link time; only at compile time) ----
GCC_OPTS       := -std=c++2b -O3 -c -fmessage-length=0 -fopenmp -static -MMD -MP \
                  -msse2 -msse3 -msse4.1 -msse4.2 -mavx -ffast-math \
                  -include pch.h
GCC_OPTS_DEBUG := -std=c++2b -Og -g -c -fmessage-length=0 -static -MMD -MP \
                  -msse2 -msse3 -msse4.1 -msse4.2 -mavx -ffast-math \
                  -include pch.h

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

#remove mat_nos_math.cpp and vec_nos_math.cpp from SRCS
SRCS := $(filter-out $(SRC_DIR)/mat_nos_math.cpp $(SRC_DIR)/vec_nos_math.cpp, $(SRCS))

# Convert sources to object names (extract basename without path/extension)
OBJ_NAMES := $(basename $(notdir $(SRCS)))
OBJS   := $(OBJ_NAMES:=.o)
OBJS_d := $(OBJ_NAMES:=.o_d)

all: NoSpherA2

# ---- Basis set converter (kept) ----
../BasisSetConverter: makefile_basis_set_converter
	@echo "Building basis set converter"
	@$(MAKE) -f makefile_basis_set_converter

# Actually run the converter to generate the file
$(SRC_DIR)/auxiliary_basis.cpp: ../BasisSetConverter
	@echo "Running BasisSetConverter..."
	@./../BasisSetConverter ../Src/basis_set_helper ../Src
	@echo "auxiliary_basis.cpp generated"

# Flags without -include pch.h for compiling the PCH itself
GCC_OPTS_NO_PCH       := $(filter-out -include pch.h,$(GCC_OPTS))
GCC_OPTS_DEBUG_NO_PCH := $(filter-out -include pch.h,$(GCC_OPTS_DEBUG))

# ---- Precompiled headers ----
# GCC with -include pch.h searches include paths for pch.h, then checks for pch.h.gch alongside it.
# If pch.h.gch is a directory, GCC searches inside it for a valid PCH matching the current flags.
# We output the .gch next to the source header so GCC auto-discovers it.
PCH_DIR   := $(SRC_DIR)/pch.h.gch
PCH_GCH   := $(PCH_DIR)/release.gch
PCH_GCH_D := $(PCH_DIR)/debug.gch

$(PCH_GCH): ../Src/pch.h
	@mkdir -p $(PCH_DIR)
	@echo 'Building precompiled header (release)...'
	@$(WITH_ONEAPI) ${COMP} $(GCC_OPTS_NO_PCH) -I"$$MKLROOT/include" $(INCLUDES) -x c++-header -o "$@" "$<"

$(PCH_GCH_D): ../Src/pch.h
	@mkdir -p $(PCH_DIR)
	@echo 'Building precompiled header (debug)...'
	@$(WITH_ONEAPI) ${COMP} $(GCC_OPTS_DEBUG_NO_PCH) -I"$$MKLROOT/include" $(INCLUDES) -x c++-header -o "$@" "$<"

# ---- auxiliary_basis depends on converter ----
auxiliary_basis.o: $(SRC_DIR)/auxiliary_basis.cpp $(PCH_GCH)
	@echo 'Compiling auxiliary_basis.o after BasisSetConverter finishes'
	@$(WITH_ONEAPI) ${COMP} ${GCC_OPTS} -I"$$MKLROOT/include" $(INCLUDES) \
	  -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -O0 -o "$@" ../Src/auxiliary_basis.cpp
	@echo 'Finished building: ../Src/auxiliary_basis.cpp'

auxiliary_basis.o_d: $(SRC_DIR)/auxiliary_basis.cpp $(PCH_GCH_D)
	@echo 'Compiling auxiliary_basis.o_d after BasisSetConverter finishes'
	@$(WITH_ONEAPI) ${COMP} ${GCC_OPTS_DEBUG} -I"$$MKLROOT/include" $(INCLUDES) \
	  -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -O0 -o "$@" ../Src/auxiliary_basis.cpp
	@echo 'Finished building: ../Src/auxiliary_basis.cpp'

# ---- Pattern rules ----
%.o: ../Src/%.cpp $(PCH_GCH)
	@echo 'Building file: $<'
	@$(WITH_ONEAPI) ${COMP} ${GCC_OPTS} -I"$$MKLROOT/include" $(INCLUDES) \
	  -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

%.o_d: ../Src/%.cpp $(PCH_GCH_D)
	@echo 'Building file (debug): $<'
	@$(WITH_ONEAPI) ${COMP} ${GCC_OPTS_DEBUG} -I"$$MKLROOT/include" $(INCLUDES) \
	  -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

# ---- Link ----
NoSpherA2: ${OBJS}
	@echo 'Building target: $@'
	@$(WITH_ONEAPI) ${COMP} -o "NoSpherA2" ${OBJS} ${LIBRASCALINE} $(MKL_LINK_FLAGS) $(LIBCINT) ../Lib/occ/lib/libocc.a ${LIBS}
	cp NoSpherA2 ../
	cp ../Lib/MKL/compiler/latest/lib/libiomp5.so ../
	@echo 'Done building NoSpherA2'

NoSpherA2_Debug: ${OBJS_d}
	@echo 'Building target: $@'
	@$(WITH_ONEAPI) ${COMP} -o "NoSpherA2_Debug" ${OBJS_d} ${LIBRASCALINE} $(MKL_LINK_FLAGS) $(LIBCINT) ../Lib/occ/lib/libocc.a ${LIBS}
	cp NoSpherA2_Debug ../
	cp ../Lib/MKL/compiler/latest/lib/libiomp5.so ../
	@echo 'Done building Debug version of NoSpherA2'

clean:
	-${RM} NoSpherA2 ../NoSpherA2 NoSpherA2_Debug ../NoSpherA2_Debug ../BasisSetConverter ${OBJS} ${OBJS_d} $(OBJ_NAMES:=.d)
	-${RM} -r $(PCH_DIR)
	-${RM} pch.h.gch pch.h.gch_d pch.h.d pch.d
	-@echo ' '

.PHONY: all clean
.SECONDARY:

# Include generated dependency files for correct incremental rebuilds.
# Only include .d files for our own objects — avoid stale/corrupt .d files
# from sub-builds (e.g. BasisSetConverter.d, basis_set_helper.d).
-include $(wildcard $(addsuffix .d,$(OBJ_NAMES)))

