#THIS IS THE MACOS MAKEFILE!#
RM := rm

all: lipo

GCC_OPTS := -std=c++2a -O3 -c -fmessage-length=0 -Xpreprocessor -fopenmp -static -MMD -MP

NoSpherA2_x86_64: ../Src/unity.cpp
	@echo 'Building target: $@'
	g++ ${GCC_OPTS} -I../llvm-project/openmp/build_x64/runtime/src -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o unity.o_x64 "$<" -target x86_64-apple-macos12.6
	g++ -o "NoSpherA2_64" ../llvm-project/openmp/build_x64/runtime/src/libomp.a unity.o_x64 -I../ -I../llvm-project/openmp/build_x64/runtime/src -L../llvm-project/openmp/build_x64/runtime/src -target x86_64-apple-macos12.6
	${RM} unity.o_x64
	@echo 'Finished building target: $@'

NoSpherA2_arm64: ../Src/unity.cpp
	@echo 'Building target: $@'
	g++ ${GCC_OPTS} -I../llvm-project/openmp/build/runtime/src -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o unity.o_arm "$<" -target aarch64-apple-macos12.6
	g++ -o "$@" ../llvm-project/openmp/build/runtime/src/libomp.a unity.o_arm -I../ -L../llvm-project/openmp/build/runtime/src -target aarch64-apple-macos12.6
	${RM} unity.o_arm
	@echo 'Finished building target: $@'

lipo: NoSpherA2_x86_64 NoSpherA2_arm64
	lipo -create -output universal_NoSpherA2 NoSpherA2_arm64 NoSpherA2_64


clean:
	-${RM} NoSpherA2
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:


