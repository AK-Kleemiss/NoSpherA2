RM := rm

#LIBS := -lgomp -static-libgcc -static-libstdc++
LIBS := -lgomp

LIBOMP := $(shell locate 'libgomp.a')

all: Wfn_2_Fchk

CPP_SRCS += \
../AtomGrid.cpp \
../basis_set.cpp \
../becke_partitioning.cpp \
../bragg.cpp \
../convenience.cpp \
../cube.cpp \
../fchk.cpp \
../grid_radial.cpp \
../mo_class.cpp \
../sphere_lebedev_rule.cpp \
../structure_factors.cpp \
../version.cpp \
../wfn_2_fchk.cpp \
../wfn_class.cpp

OBJS += \
./AtomGrid.o \
./basis_set.o \
./becke_partitioning.o \
./bragg.o \
./convenience.o \
./cube.o \
./fchk.o \
./grid_radial.o \
./mo_class.o \
./sphere_lebedev_rule.o \
./structure_factors.o \
./version.o \
./wfn_2_fchk.o \
./wfn_class.o

CPP_DEPS += \
./AtomGrid.d \
./basis_set.d \
./becke_partitioning.d \
./bragg.d \
./convenience.d \
./cube.d \
./fchk.d \
./grid_radial.d \
./mo_class.d \
./sphere_lebedev_rule.d \
./structure_factors.d \
./version.d \
./wfn_2_fchk.d \
./wfn_class.d

-include ${CPP_DEPS}

Wfn_2_Fchk: ${OBJS} ${USER_OBJS}
	@echo 'Building target: $@'
	g++ -o "Wfn_2_Fchk" ${LIBOMP} ${OBJS} ${LIBS} -I../ -L../
	@echo 'Finished building target: $@'

./AtomGrid.o: ../AtomGrid.cpp
	@echo 'Building file: $<'
	g++ -std=c++0x -O3 -Wall -c -fmessage-length=0  -fopenmp -static -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

./basis_set.o: ../basis_set.cpp
	@echo 'Building file: $<'
	g++ -std=c++0x -O3 -Wall -c -fmessage-length=0  -fopenmp -static -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

./becke_partitioning.o: ../becke_partitioning.cpp
	@echo 'Building file: $<'
	g++ -std=c++0x -O3 -Wall -c -fmessage-length=0  -fopenmp -static -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

./bragg.o: ../bragg.cpp
	@echo 'Building file: $<'
	g++ -std=c++0x -O3 -Wall -c -fmessage-length=0  -fopenmp -static -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

./convenience.o: ../convenience.cpp
	@echo 'Building file: $<'
	g++ -std=c++0x -O3 -Wall -c -fmessage-length=0  -fopenmp -static -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

./cube.o: ../cube.cpp
	@echo 'Building file: $<'
	g++ -std=c++0x -O3 -Wall -c -fmessage-length=0  -fopenmp -static -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

./fchk.o: ../fchk.cpp
	@echo 'Building file: $<'
	g++ -std=c++0x -O3 -Wall -c -fmessage-length=0  -fopenmp -static -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

./grid_radial.o: ../grid_radial.cpp
	@echo 'Building file: $<'
	g++ -std=c++0x -O3 -Wall -c -fmessage-length=0  -fopenmp -static -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

./mo_class.o: ../mo_class.cpp
	@echo 'Building file: $<'
	g++ -std=c++0x -O3 -Wall -c -fmessage-length=0  -fopenmp -static -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

./sphere_lebedev_rule.o: ../sphere_lebedev_rule.cpp
	@echo 'Building file: $<'
	g++ -std=c++0x -O3 -Wall -c -fmessage-length=0  -fopenmp -static -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

./structure_factors.o: ../structure_factors.cpp
	@echo 'Building file: $<'
	g++ -std=c++0x -O3 -Wall -c -fmessage-length=0  -fopenmp -static -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

./version.o: ../version.cpp
	@echo 'Building file: $<'
	g++ -std=c++0x -O3 -Wall -c -fmessage-length=0  -fopenmp -static -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

./wfn_2_fchk.o: ../wfn_2_fchk.cpp
	@echo 'Building file: $<'
	g++ -std=c++0x -O3 -Wall -c -fmessage-length=0  -fopenmp -static -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

./wfn_class.o: ../wfn_class.cpp
	@echo 'Building file: $<'
	g++ -std=c++0x -O3 -Wall -c -fmessage-length=0  -fopenmp -static -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'



clean:
	-${RM} ${CC_DEPS}${C++_DEPS}${EXECUTABLES} Wfn_2_Fchk
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:


