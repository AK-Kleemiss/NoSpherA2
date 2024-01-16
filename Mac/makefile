#THIS IS THE MACOS MAKEFILE!#
RM := rm

LIBOMP := /usr/local/opt/libomp/lib/libomp.a

all: NoSpherA2 lipo

CPP_SRCS += \
../Src/unity.cpp 

OBJS += \
./unity.o 

CPP_DEPS += \
./unity.d 

-include ${CPP_DEPS}

GCC_OPTS := -std=c++2a -O3 -c -fmessage-length=0 -Xpreprocessor -fopenmp -static -MMD -MP

NoSpherA2: ${OBJS} ${USER_OBJS}
	@echo 'Building target: $@'
	g++ -o "NoSpherA2" ${LIBOMP} ${OBJS} -I../ -L/usr/local/opt/libomp/lib
	@echo 'Finished building target: $@'

./unity.o: ../Src/unity.cpp
	@echo 'Building file: $<'
	g++ ${GCC_OPTS} -I/usr/local/opt/libomp/include -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

NoSpherA2_x86_64: ${OBJS} ${USER_OBJS}
	@echo 'Building target: $@'
	g++ -o "NoSpherA2_64" ${LIBOMP} ${OBJS} -I../ -L/usr/local/opt/libomp/lib -target x86_64-apple-macos12.6
	${RM} ${CPP_DEPS} ${OBJS}
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


