#THIS IS THE MACOS MAKEFILE!#
RM := rm

#LIBOMP := /usr/local/opt/libomp/lib/libomp.a

all: NoSpherA2

CPP_SRCS += \
../unity.cpp 

OBJS += \
./unity.o 

CPP_DEPS += \
./unity.d 

-include ${CPP_DEPS}

GCC_OPTS := -std=c++11 -O3 -fopenmp -static

NoSpherA2: ${OBJS} ${USER_OBJS}
	@echo 'Building target: $@'
	clang -o "NoSpherA2" ${LIBOMP} ${OBJS} -I../ -I/usr/local/include
	${RM} ${CPP_DEPS} ${OBJS}
	@echo 'Finished building target: $@'

./unity.o: ../unity.cpp
	@echo 'Building file: $<'
	clang ${GCC_OPTS} -I/usr/local/include -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

clean:
	-${RM} NoSpherA2
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:


