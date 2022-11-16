#THIS IS THE MACOS MAKEFILE!#
RM := rm

LIBOMP := /usr/local/opt/libomp/lib/libomp.a

all: NoSpherA2

CPP_SRCS += \
../unity.cpp 

OBJS += \
./unity.o 

CPP_DEPS += \
./unity.d 

-include ${CPP_DEPS}

GCC_OPTS := -std=c++2a -O3 -c -fmessage-length=0  -Xpreprocessor -fopenmp -static -MMD -MP

NoSpherA2: ${OBJS} ${USER_OBJS}
	@echo 'Building target: $@'
	g++ -o "NoSpherA2" ${LIBOMP} ${OBJS} -I../ -L/usr/local/opt/libomp/lib
	@echo 'Finished building target: $@'

./unity.o: ../unity.cpp
	@echo 'Building file: $<'
	g++ ${GGC_OPTS} -I/usr/local/opt/libomp/include -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

clean:
	-${RM} ${CPP_DEPS} ${OBJS} NoSpherA2
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:


