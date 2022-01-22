RM := rm

#LIBS := -lgomp -static-libgcc -static-libstdc++
#LIBS := -lgomp

LIBOMP := /usr/local/lib/libomp.a

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
	g++ -o "NoSpherA2" ${LIBOMP} ${OBJS} ${LIBS} -I../ -L../
	@echo 'Finished building target: $@'

./unity.o: ../unity1.cpp
	@echo 'Building file: $<'
	g++ -std=c++17 -O3   -c -fmessage-length=0  -Xpreprocessor -fopenmp -static -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

clean:
	-${RM} ${CPP_DEPS} ${OBJS} NoSpherA2
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:


