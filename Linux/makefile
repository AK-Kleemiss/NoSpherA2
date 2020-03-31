RM := rm

#LIBS := -lgomp -static-libgcc -static-libstdc++
LIBS := -lgomp

LIBOMP := /usr/lib/gcc/x86_64-linux-gnu/7/libgomp.a

all: NoSpherA2

CPP_SRCS += \
./unity1.cpp \
./unity2.cpp \
./unity3.cpp \
./unity4.cpp \
./unity5.cpp

OBJS += \
./unity1.o \
./unity2.o \
./unity3.o \
./unity4.o \
./unity5.o

CPP_DEPS += \
./unity1.d \
./unity2.d \
./unity3.d \
./unity4.d \
./unity5.d

-include ${CPP_DEPS}

NoSpherA2: ${OBJS} ${USER_OBJS}
	@echo 'Building target: $@'
	g++ -o "NoSpherA2" ${LIBOMP} ${OBJS} ${LIBS} -I../ -L../
	@echo 'Finished building target: $@'

./unity1.o: ./unity1.cpp
	@echo 'Building file: $<'
	g++ -std=c++0x -O3   -c -fmessage-length=0  -fopenmp -static -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

./unity2.o: ./unity2.cpp
	@echo 'Building file: $<'
	g++ -std=c++0x -O3   -c -fmessage-length=0  -fopenmp -static -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

./unity3.o: ./unity3.cpp
	@echo 'Building file: $<'
	g++ -std=c++0x -O3   -c -fmessage-length=0  -fopenmp -static -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

./unity4.o: ./unity4.cpp
	@echo 'Building file: $<'
	g++ -std=c++0x -O3   -c -fmessage-length=0  -fopenmp -static -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

./unity5.o: ./unity5.cpp
	@echo 'Building file: $<'
	g++ -std=c++0x -O3   -c -fmessage-length=0  -fopenmp -static -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

clean:
	-${RM} ${CC_DEPS}${C++_DEPS}${EXECUTABLES} NoSpherA2
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:


