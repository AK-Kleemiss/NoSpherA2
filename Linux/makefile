#THIS IS THE LINUX MAKEFILE#
RM := rm -f

LIBS := -lgomp -static-libgcc -static-libstdc++
LIBOMP := /usr/lib/gcc/x86_64-linux-gnu/$(shell expr `g++ -dumpversion`)/libgomp.a
LIBRASCALINE := ../Lib/librascaline.a
LIBRASCALINE += ../Lib/libmetatensor.a
LIBHDF5 := ../Lib/libhdf5.a

all: clean NoSpherA2 

OBJS := ./unity.o
CPP_DEPS := ./unity.d
GCC_OPTS := -std=c++17 -O3 -c -fmessage-length=0 -fopenmp -static -MMD -MP

NoSpherA2: ${OBJS} ${USER_OBJS}
	@echo 'Building target: $@'
	g++ -o "NoSpherA2" ${LIBOMP} ${OBJS} ${LIBS}
	${RM} ${CPP_DEPS} ${OBJS}
	cp NoSpherA2 ../
	@echo 'Finished building target: $@'

./unity.o: ../Src/unity.cpp
	@echo 'Building file: $<'
	g++ ${GCC_OPTS} -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

./unity_rasc.o: ../Src/unity.cpp
	@echo 'Building file: $< with rascaline support'
	g++ ${GCC_OPTS} -D__RASCALINE__ -I../Lib/ -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

./hdf5.o: ../Src/hdf5.cpp
	@echo 'Building file: $<'
	g++ ${GCC_OPTS} -I../Lib/ -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

NoSpherA2_rascaline: ./unity_rasc.o ./hdf5.o
	@echo 'Building target: $@'
	g++ -o "NoSpherA2_rascaline" ./unity_rasc.o ./hdf5.o ${LIBOMP} ${LIBRASCALINE} ${LIBHDF5} ${LIBS}
	${RM} unity_rasc.o hdf5.o unity_rasc.d hdf5.d
	cp NoSpherA2_rascaline ../
	@echo 'Done building rascaline version of NoSpherA2'

clean:
	-${RM} NoSpherA2 ../NoSpherA2 NoSpherA2_rascaline ../NoSpherA2_rascaline
	-@echo ' '

.PHONY: all
.SECONDARY:


