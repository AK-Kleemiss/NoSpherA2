#THIS IS THE LINUX MAKEFILE#
RM := rm -f

LIBS := -lgomp -static-libgcc -static-libstdc++ -ltbb
LIBOMP := /usr/lib/gcc/x86_64-linux-gnu/$(shell expr `g++ -dumpversion`)/libgomp.a
LIBRASCALINE := ../Lib/librascaline.a
LIBRASCALINE += ../Lib/libmetatensor.a
LIBHDF5 := ../Lib/libhdf5.a
LIBOPENBLAS := ../Lib/OpenBlas/lib/libopenblas_openmp.a

all: NoSpherA2 

OBJS := cerf.o wfn_class.o atoms.o properties.o spherical_density.o AtomGrid.o basis_set.o convenience.o sphere_lebedev_rule.o scattering_factors.o cube.o fchk.o JKFit.o NoSpherA2.o SALTED_utilities.o constants.o
OBJS_RASC := cerf.o_RASC wfn_class.o_RASC atoms.o_RASC properties.o_RASC spherical_density.o_RASC AtomGrid.o_RASC basis_set.o_RASC convenience.o_RASC sphere_lebedev_rule.o_RASC scattering_factors.o_RASC cube.o_RASC fchk.o_RASC JKFit.o_RASC NoSpherA2.o_RASC SALTED_utilities.o_RASC constants.o_RASC
CPP_DEPS := unity.d
GCC_OPTS := -std=c++17 -O3 -c -fmessage-length=0 -fopenmp -static -MMD -MP

NoSpherA2: ${OBJS} ${USER_OBJS}
	@echo 'Building target: $@'
	g++ -o "NoSpherA2" ${LIBOMP} ${OBJS} ${LIBS}
	${RM} ${CPP_DEPS} ${OBJS}
	cp NoSpherA2 ../
	@echo 'Finished building target: $@'

# Pattern rule to compile all .cpp files
%.o: ../Src/%.cpp
	@echo 'Building file: $<'
	g++ ${GCC_OPTS} -I../Lib/cerf -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

# Pattern rule to compile all .cpp files
%.o_RASC: ../Src/%.cpp
	@echo 'Building file: $<'
	g++ ${GCC_OPTS} -D__RASCALINE__ -I../Lib/ -I../Lib/cerf -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

SALTED.o: ../Src/SALTED.cpp
	@echo 'Building file: $< with rascaline support'
	g++ ${GCC_OPTS} -D__RASCALINE__ -I../Lib/ -I../Lib/OpenBlas/include_Linux -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

hdf5.o: ../Src/hdf5.cpp
	@echo 'Building file: $<'
	g++ ${GCC_OPTS} -I../Lib/ -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

NoSpherA2_Rascaline: ${OBJS_RASC} hdf5.o SALTED.o cerf.o
	@echo 'Building target: $@'
	g++ -o "NoSpherA2_rascaline" ${OBJS_RASC} hdf5.o SALTED.o ${LIBOMP} ${LIBRASCALINE} ${LIBHDF5} ${LIBOPENBLAS} ${LIBS}
	${RM} ${OBJS_RASC} hdf5.o unity_rasc.d hdf5.d SALTED.o SALTED.d cerf.o cerf.d
	cp NoSpherA2_rascaline ../
	@echo 'Done building rascaline version of NoSpherA2'

clean:
	-${RM} NoSpherA2 ../NoSpherA2 NoSpherA2_rascaline ../NoSpherA2_rascaline
	-@echo ' '

.PHONY: all
.SECONDARY:


