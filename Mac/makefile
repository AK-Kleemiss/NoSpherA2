#THIS IS THE MACOS MAKEFILE!#
RM := rm

all: lipo

GCC_OPTS := -std=c++2a -O3 -c -fmessage-length=0 -Xpreprocessor -fopenmp -static -MMD -MP

cerf.o_x86_64: ../Src/cerf.cpp
	@echo 'Building file: $<'
	g++ ${GCC_OPTS} -I../llvm-project/openmp/build_x64/runtime/src -I../Lib/cerf -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<" -target x86_64-apple-macos14.0
	@echo 'Finished building: $<'

NoSpherA2_x86_64: ../Src/unity.cpp cerf.o_x86_64
	@echo 'Building target: $@'
	g++ ${GCC_OPTS} -I../llvm-project/openmp/build_x64/runtime/src -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o unity.o_x64 "$<" -target x86_64-apple-macos14.0
	g++ -o "NoSpherA2_64" ../llvm-project/openmp/build_x64/runtime/src/libomp.a unity.o_x64 cerf.o_x86_64 -I../ -I../llvm-project/openmp/build_x64/runtime/src -L../llvm-project/openmp/build_x64/runtime/src -target x86_64-apple-macos14.0
	${RM} unity.o_x64 cerf.o_x86_64
	@echo 'Finished building target: $@'

cerf.o_arm64: ../Src/cerf.cpp
	@echo 'Building file: $<'
	g++ ${GCC_OPTS} -I../llvm-project/openmp/build/runtime/src -I../Lib/cerf -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<" -target aarch64-apple-macos14.0
	@echo 'Finished building: $<'

NoSpherA2_arm64: ../Src/unity.cpp
	@echo 'Building target: $@'
	g++ ${GCC_OPTS} -I../llvm-project/openmp/build/runtime/src -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o unity.o_arm "$<" -target aarch64-apple-macos14.0
	g++ -o "$@" ../llvm-project/openmp/build/runtime/src/libomp.a unity.o_arm cerf.o_arm64 -I../ -L../llvm-project/openmp/build/runtime/src -target aarch64-apple-macos14.0
	${RM} unity.o_arm
	@echo 'Finished building target: $@'

lipo: NoSpherA2_x86_64 NoSpherA2_arm64
	lipo -create -output universal_NoSpherA2 NoSpherA2_arm64 NoSpherA2_64


clean:
	-${RM} NoSpherA2
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:


