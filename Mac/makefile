#THIS IS THE MACOS MAKEFILE!#
RM := rm

all: lipo

OBJ_NAMES := wfn_class atoms properties spherical_density AtomGrid basis_set convenience sphere_lebedev_rule scattering_factors cube fchk JKFit NoSpherA2 SALTED_utilities constants ECPs_corrections integrator math isosurface
OBJS_arm64 := $(OBJ_NAMES:=.o_arm64)
OBJS_x86_64 := $(OBJ_NAMES:=.o_x86_64)

GCC_OPTS := -std=c++2a -O3 -c -fmessage-length=0 -Xpreprocessor -fopenmp -static -MMD -MP

%.o_x86_64: ../Src/%.cpp
	@echo 'Building file: $<'
	g++ ${GCC_OPTS} -I../llvm-project/openmp/build_x64/runtime/src -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<" -target x86_64-apple-macos14.0
	@echo 'Finished building: $<'

NoSpherA2_x86_64: ${OBJS_x86_64}
	@echo 'Building target: $@'
	g++ -o "NoSpherA2_64" ../llvm-project/openmp/build_x64/runtime/src/libomp.a ${OBJS_x86_64} -I../ -I../llvm-project/openmp/build_x64/runtime/src -L../llvm-project/openmp/build_x64/runtime/src -target x86_64-apple-macos14.0
	${RM} ${OBJS_x86_64}
	@echo 'Finished building target: $@'

%.o_arm64: ../Src/%.cpp
	@echo 'Building file: $<'
	g++ ${GCC_OPTS} -I../llvm-project/openmp/build/runtime/src -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<" -target aarch64-apple-macos14.0
	@echo 'Finished building: $<'

NoSpherA2_arm64: ${OBJS_arm64}
	@echo 'Building target: $@'
	g++ -o "$@" ../llvm-project/openmp/build/runtime/src/libomp.a ${OBJS_arm64} -I../ -L../llvm-project/openmp/build/runtime/src -target aarch64-apple-macos14.0
	${RM} ${OBJS_arm64}
	@echo 'Finished building target: $@'

lipo: NoSpherA2_x86_64 NoSpherA2_arm64
	lipo -create -output universal_NoSpherA2 NoSpherA2_arm64 NoSpherA2_64


clean:
	-${RM} NoSpherA2
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:


