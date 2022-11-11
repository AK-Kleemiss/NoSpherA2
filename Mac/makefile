#THIS IS THE MACOS MAKEFILE!#
RM := rm

#LIBOMP := /usr/local/lib/libomp.a

all: NoSpherA2

CPP_SRCS += \
../unity.cpp 

OBJS += \
./unity.o 

CPP_DEPS += \
./unity.d 

-include ${CPP_DEPS}

NoSpherA2: ${OBJS} ${USER_OBJS}
	@echo 'Building target: $@'
	g++ -o "NoSpherA2" ${OBJS} -I../ -L../ -L/usr/local/opt/libomp/lib -L/usr/local/libomp/lib
	@echo 'Finished building target: $@'

./unity.o: ../unity.cpp
	@echo 'Building file: $<'
	g++ -std=c++17 -O3 -c -fmessage-length=0  -Xpreprocessor -fopenmp -lomp -static -I/usr/local/opt/libomp/include -I/usr/local/libomp/include -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

clean:
	-${RM} ${CPP_DEPS} ${OBJS} NoSpherA2
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:


