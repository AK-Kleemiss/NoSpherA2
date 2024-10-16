#THIS IS THE MACOS MAKEFILE!#
RM := rm

all: lipo

OBJS_x86_64 := wfn_class.o_x86_64 atoms.o_x86_64 properties.o_x86_64 spherical_density.o_x86_64 AtomGrid.o_x86_64 basis_set.o_x86_64 convenience.o_x86_64 sphere_lebedev_rule.o_x86_64 scattering_factors.o_x86_64 cube.o_x86_64 fchk.o_x86_64 JKFit.o_x86_64 NoSpherA2.o_x86_64 SALTED_utilities.o_x86_64 constants.o_x86_64 ECPs_corrections.o_x86_64
OBJS_arm64 := wfn_class.o_arm64 atoms.o_arm64 properties.o_arm64 spherical_density.o_arm64 AtomGrid.o_arm64 basis_set.o_arm64 convenience.o_arm64 sphere_lebedev_rule.o_arm64 scattering_factors.o_arm64 cube.o_arm64 fchk.o_arm64 JKFit.o_arm64 NoSpherA2.o_arm64 SALTED_utilities.o_arm64 constants.o_arm64 ECPs_corrections.o_arm64

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


