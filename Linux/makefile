# === LINUX MAKEFILE ===
RM := rm -f

COMP   := g++
C_COMP := gcc

LIBS := -Wl,-Bstatic -lstdc++ -lpthread -Wl,-Bdynamic -lm -lc

# oneAPI env (source setvars.sh before each compile/link)
SETVARS ?= ../Lib/MKL/setvars.sh
WITH_ONEAPI = . "$(SETVARS)" >/dev/null 2>&1 &&

# ---- Optional MKL knobs ----
# LP64 vs ILP64 (add -DMKL_ILP64 automatically for ILP64)
MKL_INTERFACE ?= LP64      # LP64 | ILP64
MKL_THREADING ?= INTEL     # INTEL | SEQ  (INTEL = mkl_intel_thread + libiomp5)

# Featomic local libs
LIBRASCALINE   := ../Lib/featomic_install/lib/libfeatomic.a
LIBRASCALINE   += ../Lib/featomic_install/lib/libmetatensor.a
LIBRASCALINE_I := ../Lib/featomic_install/include

# ---- Compiler flags (no -fopenmp at link time; only at compile time) ----
GCC_OPTS       := -std=c++2b -O3 -c -fmessage-length=0 -fopenmp -static -MMD -MP \
                  -msse2 -msse3 -msse4.1 -msse4.2 -mavx -ffast-math
GCC_OPTS_DEBUG := -std=c++2b -Og -g -c -fmessage-length=0 -static -MMD -MP \
                  -msse2 -msse3 -msse4.1 -msse4.2 -mavx -ffast-math

# Threading backend selection
ifeq ($(MKL_THREADING),SEQ)
  MKL_THR_LIB = -lmkl_sequential
  IOMP_LIB    =
else
  MKL_THR_LIB = -lmkl_intel_thread
  IOMP_LIB    = -liomp5
endif

# MKL link set per Intel guidance (+ rpaths so setvars isn't required at runtime)
MKL_LINK_FLAGS = \
  -L"$$MKLROOT/lib/intel64" -Wl,-rpath,"$$MKLROOT/lib/intel64" \
  -Wl,--start-group -lmkl_intel_lp64 $(MKL_THR_LIB) -lmkl_core -Wl,--end-group \
  $(IOMP_LIB) -lpthread -lm -ldl \
  -L"$$ONEAPI_ROOT/compiler/latest/linux/lib" -Wl,-rpath,"$$ONEAPI_ROOT/compiler/latest/linux/lib" \
  -L"$$ONEAPI_ROOT/compiler/latest/linux/compiler/lib/intel64_lin" -Wl,-rpath,"$$ONEAPI_ROOT/compiler/latest/linux/compiler/lib/intel64_lin"

# ---- Objects ----
OBJ_NAMES := wfn_class atoms properties spherical_density AtomGrid basis_set convenience sphere_lebedev_rule scattering_factors cube fchk JKFit NoSpherA2 SALTED_utilities constants ECPs_corrections libCintMain nos_math SALTED_io SALTED_predictor SALTED_equicomb isosurface cart2sph int_g2e int_optimizer integration_params integrator rys_roots auxiliary_basis qct b2c bondwise_analysis
OBJS   := $(OBJ_NAMES:=.o)
OBJS_d := $(OBJ_NAMES:=.o_d)

all: NoSpherA2

# ---- Basis set converter (kept) ----
basis_set_converter:
	@$(MAKE) -f makefile_basis_set_converter

# ---- Precompiled headers ----
pch.h.gch: ../Src/pch.h
	@$(WITH_ONEAPI) ${COMP} ${GCC_OPTS} -x c++-header -o "$@" "$<"

pch.h.gch_d: ../Src/pch.h
	@$(WITH_ONEAPI) ${COMP} ${GCC_OPTS_DEBUG} -x c++-header -o "$@" "$<"

# ---- auxiliary_basis depends on converter ----
auxiliary_basis.o: basis_set_converter
	@echo 'Compiling auxiliary_basis.o after BasisSetConverter finishes'
	@$(WITH_ONEAPI) ${COMP} ${GCC_OPTS} -I"$$MKLROOT/include" -I${LIBRASCALINE_I} \
	  -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" ../Src/auxiliary_basis.cpp
	@echo 'Finished building: auxiliary_basis.o'

auxiliary_basis.o_d: basis_set_converter
	@echo 'Compiling auxiliary_basis.o_d after BasisSetConverter finishes'
	@$(WITH_ONEAPI) ${COMP} ${GCC_OPTS_DEBUG} -I"$$MKLROOT/include" -I${LIBRASCALINE_I} \
	  -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" ../Src/auxiliary_basis.cpp
	@echo 'Finished building: auxiliary_basis.o_d'

# ---- Pattern rules ----
%.o: ../Src/%.cpp pch.h.gch
	@echo 'Building file: $<'
	@$(WITH_ONEAPI) ${COMP} ${GCC_OPTS} -I"$$MKLROOT/include" -I${LIBRASCALINE_I} \
	  -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

%.o_d: ../Src/%.cpp pch.h.gch_d
	@echo 'Building file (debug): $<'
	@$(WITH_ONEAPI) ${COMP} ${GCC_OPTS_DEBUG} -I"$$MKLROOT/include" -I${LIBRASCALINE_I} \
	  -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

# ---- Link ----
NoSpherA2: ${OBJS}
	@echo 'Building target: $@'
	@$(WITH_ONEAPI) ${COMP} -o "NoSpherA2" ${OBJS} ${LIBRASCALINE} $(MKL_LINK_FLAGS) ${LIBS}
	cp NoSpherA2 ../
	@echo 'Done building NoSpherA2'

NoSpherA2_Debug: ${OBJS_d}
	@echo 'Building target: $@'
	@$(WITH_ONEAPI) ${COMP} -o "NoSpherA2_Debug" ${OBJS_d} ${LIBRASCALINE} $(MKL_LINK_FLAGS) ${LIBS}
	cp NoSpherA2_Debug ../
	@echo 'Done building Debug version of NoSpherA2'

clean:
	-${RM} NoSpherA2 ../NoSpherA2 NoSpherA2_Debug ../NoSpherA2_Debug ../BasisSetConverter ${OBJS} ${OBJS_d}
	-@echo ' '

.PHONY: all clean basis_set_converter
.SECONDARY:
