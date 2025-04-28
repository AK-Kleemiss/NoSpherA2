#THIS IS THE MACOS MAKEFILE!#
RM := rm

# Detect architecture
ARCH := $(shell uname -m)

all: NoSpherA2

COMP := g++

OBJ_NAMES := wfn_class atoms properties spherical_density AtomGrid basis_set convenience sphere_lebedev_rule scattering_factors cube fchk JKFit NoSpherA2 SALTED_utilities constants ECPs_corrections libCintMain nos_math SALTED_io SALTED_predictor SALTED_equicomb isosurface cart2sph int_g2e int_optimizer integration_params integrator rys_roots b2c bondwise_analysis
OBJS := $(OBJ_NAMES:=.o_$(ARCH))
OBJS_DEBUG := $(OBJ_NAMES:=.o_$(ARCH)_d)
OMP_BUILD := ../llvm-project/openmp/build/runtime/src
blas := ../OpenBLAS/installation/

GCC_OPTS := -std=c++2a -O3 -c -fmessage-length=0 -Xpreprocessor -fopenmp -static -MMD -MP
GCC_OPTS_DEBUG := -std=c++2a -Og -g -c -fmessage-length=0 -Xpreprocessor -MP

# Add architecture-specific flags
ifeq ($(ARCH), arm64)
    GCC_OPTS += -mcpu=apple-m1
    GCC_OPTS_DEBUG += -mcpu=apple-m1
    TARGET_FLAG := -target aarch64-apple-macos14.0
	LIBRASCALINE := ../Lib/librascaline_Mac.a
	LIBRASCALINE += ../Lib/libmetatensor_Mac.a
else ifeq ($(ARCH), x86_64)
    GCC_OPTS += -march=native
    GCC_OPTS_DEBUG += -march=native
    TARGET_FLAG := -target x86_64-apple-macos14.0
	LIBRASCALINE := ../Lib/librascaline_Mac_x86_64.a
	LIBRASCALINE += ../Lib/libmetatensor_Mac_x86_64.a
else
    $(error Unsupported architecture: $(ARCH))
endif

%.o_arm64: ../Src/%.cpp
	@echo 'Building file: $<'
	${COMP} ${GCC_OPTS} -I${OMP_BUILD} -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -I../Lib -I${blas}/include -o "$@" "$<" ${TARGET_FLAG}
	@echo 'Finished building: $<'

NoSpherA2: ${OBJS}
	@echo 'Building target: $@'
	${COMP} -o "$@" ${OMP_BUILD}/libomp.a ${blas}/lib/libopenblas.a ${OBJS} ${LIBRASCALINE} -I../ -I${blas}/include -L${OMP_BUILD} -L${blas}/lib -lcblas ${TARGET_FLAG}
	${RM} ${OBJS}
	cp NoSpherA2 ../
	@echo 'Finished building target: $@'

%.o_arm64_d: ../Src/%.cpp
	@echo 'Building file: $<'
	${COMP} ${GCC_OPTS_DEBUG} -I${OMP_BUILD} -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -I../Lib -I${blas}/include -o "$@" "$<" ${TARGET_FLAG}
	@echo 'Finished building: $<'

NoSpherA2_Debug: ${OBJS_DEBUG}
	@echo 'Building target: $@'
	${COMP} -o "$@" ${OMP_BUILD}/libomp.a ${blas}/lib/libopenblas.a ${OBJS_DEBUG} ${LIBRASCALINE} -I../ -I${blas}/include -L${OMP_BUILD} -L${blas}/lib -lcblas ${TARGET_FLAG}
	${RM} ${OBJS_DEBUG}
	cp NoSpherA2_Debug ../
	@echo 'Finished building target: $@'

clean:
	-${RM} NoSpherA2
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:


