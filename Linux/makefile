#THIS IS THE LINUX MAKEFILE#
RM := rm -f

LIBS := -lgomp -static-libgcc -static-libstdc++
LIBOMP := /usr/lib/gcc/x86_64-linux-gnu/$(shell expr `g++ -dumpversion`)/libgomp.a
LIBRASCALINE := ../Lib/librascaline.a
LIBMETATENSOR := ../Lib/libmetatensor.a
LIBHDF5 := ../Lib/libhdf5_serial.a

all: clean NoSpherA2

CPP_SRCS += ../unity.cpp 

OBJS += ./unity.o 

CPP_DEPS += ./unity.d 

-include ${CPP_DEPS}

GCC_OPTS := -std=c++17 -O3 -c -fmessage-length=0 -fopenmp -static -MMD -MP

NoSpherA2: ${OBJS} ${USER_OBJS}
	@echo 'Building target: $@'
	g++ -o "NoSpherA2" ${LIBOMP} ${OBJS} ${LIBS} -I../
	${RM} ${CPP_DEPS} ${OBJS}
	cp NoSpherA2 ../
	@echo 'Finished building target: $@'

./unity.o: ../Src/unity.cpp
	@echo 'Building file: $<'
	g++ ${GCC_OPTS} -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

NoSpherA2_rascaline: ${OBJS}
	@echo 'Building target: $@'
	g++ -o "NoSpherA2_rascaline" ${LIBOMP} ${LIBRASCALINE} ${LIBMETATENSOR} ${LIBHDF5} ${OBJS} ${LIBS} -lrascaline -lmetatensor -lhdf5_serial -I../ -I../Lib/ -L../Lib/
	cp NoSpherA2_rascaline ../
	@echo 'Done building rascaline version of NoSpherA2'

clean:
	-${RM} NoSpherA2 ../NoSpherA2
	-@echo ' '

.PHONY: all
.SECONDARY:


