#THIS IS THE MACOS MAKEFILE!#
RM := rm

all: NoSpherA2

OBJ_NAMES := wfn_class atoms properties spherical_density AtomGrid basis_set convenience sphere_lebedev_rule scattering_factors cube fchk JKFit NoSpherA2 SALTED_utilities constants ECPs_corrections libCintMain nos_math SALTED_io SALTED_predictor SALTED_equicomb isosurface cart2sph int_g2e int_optimizer integration_params integrator rys_roots
OBJS_arm64 := $(OBJ_NAMES:=.o_arm64)
OBJS_arm64_d := $(OBJ_NAMES:=.o_arm64_d)
OMP_BUILD_arm64 := ../llvm-project/openmp/build/runtime/src
blas_arm64 := ../OpenBLAS/installation/
LIBRASCALINE := ../Lib/librascaline_Mac.a
LIBRASCALINE += ../Lib/libmetatensor_Mac.a

GCC_OPTS := -std=c++2a -O3 -c -fmessage-length=0 -Xpreprocessor -fopenmp -static -MMD -MP
GCC_OPTS_DEBUG := -std=c++20 -Og -g -c -fmessage-length=0 -Xpreprocessor -MP

%.o_arm64: ../Src/%.cpp
	@echo 'Building file: $<'
	g++ ${GCC_OPTS} -I${OMP_BUILD_arm64} -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -I../Lib -I${blas_arm64}/include -o "$@" "$<" -target aarch64-apple-macos14.0
	@echo 'Finished building: $<'

NoSpherA2: ${OBJS_arm64}
	@echo 'Building target: $@'
	g++ -o "$@" ${OMP_BUILD_arm64}/libomp.a ${blas_arm64}/lib/libopenblas.a ${OBJS_arm64} ${LIBRASCALINE} -I../ -I${blas_arm64}/include -L${OMP_BUILD_arm64} -L${blas_arm64}/lib -lcblas -target aarch64-apple-macos14.0
	${RM} ${OBJS_arm64}
	cp NoSpherA2 ../
	@echo 'Finished building target: $@'

%.o_arm64_d: ../Src/%.cpp
	@echo 'Building file: $<'
	g++ ${GCC_OPTS_DEBUG} -I${OMP_BUILD_arm64} -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -I../Lib -I${blas_arm64}/include -o "$@" "$<" -target aarch64-apple-macos14.0
	@echo 'Finished building: $<'

NoSpherA2_Debug: ${OBJS_arm64_d}
	@echo 'Building target: $@'
	g++ -o "$@" ${OMP_BUILD_arm64}/libomp.a ${blas_arm64}/lib/libopenblas.a ${OBJS_arm64_d} ${LIBRASCALINE} -I../ -I${blas_arm64}/include -L${OMP_BUILD_arm64} -L${blas_arm64}/lib -lcblas -target aarch64-apple-macos14.0
	${RM} ${OBJS_arm64_d}
	cp NoSpherA2_Debug ../
	@echo 'Finished building target: $@'

clean:
	-${RM} NoSpherA2
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:


