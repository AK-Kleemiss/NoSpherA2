#THIS IS THE MACOS MAKEFILE!#
RM := rm

all: NoSpherA2

COMP := g++

OBJ_NAMES := wfn_class atoms properties spherical_density AtomGrid basis_set convenience sphere_lebedev_rule scattering_factors cube fchk JKFit NoSpherA2 SALTED_utilities constants ECPs_corrections libCintMain nos_math SALTED_io SALTED_predictor SALTED_equicomb isosurface cart2sph int_g2e int_optimizer integration_params integrator rys_roots b2c bondwise_analysis qct auxiliary_basis
OBJS_a := $(OBJ_NAMES:=.o_arm64)
OBJS_DEBUG_a := $(OBJ_NAMES:=.o_arm64_d)

OBJS_x := $(OBJ_NAMES:=.o_x86_64)
OBJS_DEBUG_x := $(OBJ_NAMES:=.o_x86_64_d)

GCC_OPTS := -std=c++2a -O3 -c -fmessage-length=0 -Xpreprocessor -fopenmp -static -MMD -MP
GCC_OPTS_DEBUG := -std=c++2a -Og -g -c -fmessage-length=0 -Xpreprocessor -MP

# Add architecture-specific flags
OMP_BUILD_a := ../llvm-project/openmp/build/runtime/src
blas_a := ../OpenBLAS/installation/
GCC_OPTS_a := ${GCC_OPTS} -mcpu=apple-m1
GCC_OPTS_DEBUG_a := ${GCC_OPTS_DEBUG} -mcpu=apple-m1
TARGET_FLAG_a := -target aarch64-apple-macos13.0
LIBFEATOMIC_a := ../featomic/featomic_install_arm/lib/libfeatomic.a
LIBFEATOMIC_a += ../featomic/featomic_install_arm/lib/libmetatensor.a
LIBFEATOMIC_a_inc += ../featomic/featomic_install_arm/include

OMP_BUILD_x := ../llvm-project/openmp/build_x86_64/runtime/src
blas_x := ../OpenBLAS/installationx86_64/
GCC_OPTS_x := ${GCC_OPTS} -march=x86_64
GCC_OPTS_DEBUG_x := ${GCC_OPTS_DEBUG} -march=x86_64
TARGET_FLAG_x := -target x86_64-apple-macos13.0
LIBFEATOMIC_x := ../featomic/featomic_install_x86/lib/libfeatomic.a
LIBFEATOMIC_x += ../featomic/featomic_install_x86/lib/libmetatensor.a
LIBFEATOMIC_x_inc += ../featomic/featomic_install_x86/include

#Basis set converter
basis_set_converter:
	@make -f makefile_basis_set_converter

# Ensure auxiliary_basis.o depends on basis_set_converter to ensure proper sequence
auxiliary_basis.o_arm64: basis_set_converter
	@echo 'Compiling auxiliary_basis.o after BasisSetConverter finishes'
	@${COMP} ${GCC_OPTS} -I${OMP_BUILD_arm64} -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -I../Lib -I${blas_arm64}/include -o "$@" ../Src/auxiliary_basis.cpp -target aarch64-apple-macos14.0
	@echo 'Finished building: auxiliary_basis.o'

%.o_arm64: ../Src/%.cpp
	@echo 'Building file: $<'
	@${COMP} ${GCC_OPTS_a} -I${OMP_BUILD_a} -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -I../featomic/featomic_install/include -I${blas_a}/include -o "$@" "$<" ${TARGET_FLAG_a}
	@echo 'Finished building: $<'

NoSpherA2_arm: ${OBJS_a}
	@echo 'Building target: $@'
	@${COMP} -o "$@" ${OMP_BUILD_a}/libomp.a ${blas_a}/lib/libopenblas.a ${OBJS_a} ${LIBFEATOMIC_a} -I../ -I${blas_a}/include -L${OMP_BUILD_a} -L${blas_a}/lib -lcblas ${TARGET_FLAG_a}
	${RM} ${OBJS_a}
	@echo 'Finished building target: $@'

# Ensure auxiliary_basis.o_d depends on basis_set_converter to ensure proper sequence
auxiliary_basis.o_arm64_d: basis_set_converter
	@echo 'Compiling auxiliary_basis.o after BasisSetConverter finishes'
	@${COMP} ${GCC_OPTS_DEBUG} -I${OMP_BUILD_arm64} -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -I../Lib -I${blas_arm64}/include -o "$@" ../Src/auxiliary_basis.cpp -target aarch64-apple-macos14.0
	@echo 'Finished building: auxiliary_basis.o'

%.o_arm64_d: ../Src/%.cpp
	@echo 'Building file: $<'
	@${COMP} ${GCC_OPTS_DEBUG_a} -I${OMP_BUILD_a} -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -I../featomic/featomic_install/include -I${blas_a}/include -o "$@" "$<" ${TARGET_FLAG_a}
	@echo 'Finished building: $<'

NoSpherA2_Debug_arm: ${OBJS_DEBUG_a}
	@echo 'Building target: $@'
	@${COMP} -o "$@" ${OMP_BUILD_a}/libomp.a ${blas_a}/lib/libopenblas.a ${OBJS_DEBUG_a} ${LIBFEATOMIC_a} -I../ -I${blas_a}/include -L${OMP_BUILD_a} -L${blas_a}/lib -lcblas ${TARGET_FLAG_a}
	@${RM} ${OBJS_DEBUG_a}
	@echo 'Finished building target: $@'

%.o_x86: ../Src/%.cpp
	@echo 'Building file: $<'
	@${COMP} ${GCC_OPTS_x} -I${OMP_BUILD_x} -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -I../featomic/featomic_install/include -I${blas_x}/include -o "$@" "$<" ${TARGET_FLAG_x}
	@echo 'Finished building: $<'

NoSpherA2_x86: ${OBJS_x}
	@echo 'Building target: $@'
	@${COMP} -o "$@" ${OMP_BUILD_x}/libomp.a ${blas_x}/lib/libopenblas.a ${OBJS_x} ${LIBFEATOMIC_x} -I../ -I${blas_x}/include -L${OMP_BUILD_x} -L${blas_x}/lib -lcblas ${TARGET_FLAG_x}
	${RM} ${OBJS_x}
	@echo 'Finished building target: $@'

%.o_x86_d: ../Src/%.cpp
	@echo 'Building file: $<'
	@${COMP} ${GCC_OPTS_DEBUG_x} -I${OMP_BUILD_x} -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -I../featomic/featomic_install/include -I${blas_x}/include -o "$@" "$<" ${TARGET_FLAG_x}
	@echo 'Finished building: $<'

NoSpherA2_Debug_x86: ${OBJS_DEBUG_x}
	@echo 'Building target: $@'
	@${COMP} -o "$@" ${OMP_BUILD_x}/libomp.a ${blas_x}/lib/libopenblas.a ${OBJS_DEBUG_x} ${LIBFEATOMIC_x} -I../ -I${blas_x}/include -L${OMP_BUILD_x} -L${blas_x}/lib -lcblas ${TARGET_FLAG_x}
	@${RM} ${OBJS_DEBUG_x}
	@echo 'Finished building target: $@'

NoSpherA2: NoSpherA2_arm NoSpherA2_x86
	@echo 'Building target: $@'
	@lipo -create -output NoSpherA2 NoSpherA2_arm NoSpherA2_x86
	@${RM} NoSpherA2_arm NoSpherA2_x86
	@echo 'Finished building target: $@'

NoSpherA2_Debug: NoSpherA2_Debug_arm NoSpherA2_Debug_x86
	@echo 'Building target: $@'
	@lipo -create -output NoSpherA2_Debug NoSpherA2_Debug_arm NoSpherA2_Debug_x86
	@${RM} NoSpherA2_Debug_arm NoSpherA2_Debug_x86
	@echo 'Finished building target: $@'

clean:
	-${RM} NoSpherA2
	-@echo ' '

.PHONY: all clean dependents basis_set_converter
.SECONDARY:


