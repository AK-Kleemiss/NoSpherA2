#THIS IS THE MACOS MAKEFILE!#
RM := rm

# Detect architecture
ARCH := $(shell uname -m)

all: NoSpherA2

COMP := g++

OBJ_NAMES := wfn_class atoms properties spherical_density AtomGrid basis_set convenience sphere_lebedev_rule scattering_factors cube fchk JKFit NoSpherA2 SALTED_utilities constants ECPs_corrections libCintMain nos_math SALTED_io SALTED_predictor SALTED_equicomb isosurface cart2sph int_g2e int_optimizer integration_params integrator rys_roots b2c bondwise_analysis qct
OBJS := $(OBJ_NAMES:=.o_$(ARCH))
OBJS_DEBUG := $(OBJ_NAMES:=.o_$(ARCH)_d)

GCC_OPTS := -std=c++2a -O3 -c -fmessage-length=0 -Xpreprocessor -fopenmp -static -MMD -MP
GCC_OPTS_DEBUG := -std=c++2a -Og -g -c -fmessage-length=0 -Xpreprocessor -MP

# Add architecture-specific flags
OMP_BUILD_a := ../llvm-project/openmp/build/runtime/src
blas_a := ../OpenBLAS/installation/
GCC_OPTS_a := ${GCC_OPTS} -mcpu=apple-m1
GCC_OPTS_DEBUG_a := ${GCC_OPTS_DEBUG} -mcpu=apple-m1
TARGET_FLAG_a := -target aarch64-apple-macos14.0
LIBFEATOMIC_a := ../featomic/featomic_install_arm/lib/libfeatomic.a
LIBFEATOMIC_a += ../featomic/featomic_install_arm/lib/libmetatensor.a
LIBFEATOMIC_a_inc += ../featomic/featomic_install_arm/include

OMP_BUILD_x := ../llvm-project/openmp/build_x86_64/runtime/src
blas_x := ../OpenBLAS/installationx86_64/
GCC_OPTS_x := ${GCC_OPTS} -march=x86_64
GCC_OPTS_DEBUG_x := ${GCC_OPTS_DEBUG} -march=x86_64
TARGET_FLAG_x := -target x86_64-apple-macos14.0
LIBFEATOMIC_x := ../featomic/featomic_install_x86/lib/libfeatomic.a
LIBFEATOMIC_x += ../featomic/featomic_install_x86/lib/libmetatensor.a
LIBFEATOMIC_x_inc += ../featomic/featomic_install_x86/include

%.o_arm64: ../Src/%.cpp
	@echo 'Building file: $<'
	@${COMP} ${GCC_OPTS} -I${OMP_BUILD} -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -I../featomic/featomic_install/include -I${blas}/include -o "$@" "$<" ${TARGET_FLAG}
	@echo 'Finished building: $<'

NoSpherA2: ${OBJS}
	@echo 'Building target: $@'
	@${COMP} -o "$@" ${OMP_BUILD}/libomp.a ${blas}/lib/libopenblas.a ${OBJS} ${LIBFEATOMIC} -I../ -I${blas}/include -L${OMP_BUILD} -L${blas}/lib -lcblas ${TARGET_FLAG}
	${RM} ${OBJS}
	@cp NoSpherA2 ../
	@echo 'Finished building target: $@'

%.o_arm64_d: ../Src/%.cpp
	@echo 'Building file: $<'
	@${COMP} ${GCC_OPTS_DEBUG} -I${OMP_BUILD} -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -I../featomic/featomic_install/include -I${blas}/include -o "$@" "$<" ${TARGET_FLAG}
	@echo 'Finished building: $<'

NoSpherA2_Debug: ${OBJS_DEBUG}
	@echo 'Building target: $@'
	@${COMP} -o "$@" ${OMP_BUILD}/libomp.a ${blas}/lib/libopenblas.a ${OBJS_DEBUG} ${LIBFEATOMIC} -I../ -I${blas}/include -L${OMP_BUILD} -L${blas}/lib -lcblas ${TARGET_FLAG}
	${RM} ${OBJS_DEBUG}
	@cp NoSpherA2_Debug ../
	@echo 'Finished building target: $@'

clean:
	-${RM} NoSpherA2
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:


