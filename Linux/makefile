#THIS IS THE LINUX MAKEFILE#
RM := rm -f

LIBS := -lgomp -static-libgcc -static-libstdc++
LIBOMP := /usr/lib/gcc/x86_64-linux-gnu/$(shell expr `g++ -dumpversion`)/libgomp.a
LIBRASCALINE := ../Lib/librascaline.a
LIBRASCALINE += ../Lib/libmetatensor.a
LIBHDF5 := ../Lib/libhdf5.a
LIBOPENBLAS := ../Lib/OpenBlas/lib/libopenblas_openmp.a /usr/lib/x86_64-linux-gnu/liblapacke.a /usr/lib/x86_64-linux-gnu/liblapack.a

all: NoSpherA2 

OBJS := wfn_class.o atoms.o properties.o spherical_density.o AtomGrid.o basis_set.o convenience.o sphere_lebedev_rule.o scattering_factors.o cube.o fchk.o JKFit.o NoSpherA2.o SALTED_utilities.o constants.o ECPs_corrections.o integrator.o SALTED_predictor.o
OBJS_RASC := wfn_class.o_RASC atoms.o_RASC properties.o_RASC spherical_density.o_RASC AtomGrid.o_RASC basis_set.o_RASC convenience.o_RASC sphere_lebedev_rule.o_RASC scattering_factors.o_RASC cube.o_RASC fchk.o_RASC JKFit.o_RASC NoSpherA2.o_RASC SALTED_utilities.o_RASC constants.o_RASC ECPs_corrections.o_RASC integrator.o_RASC hdf5.o_RASC SALTED.o_RASC
OBJS_RASC_d := wfn_class.o_RASC_d atoms.o_RASC_d properties.o_RASC_d spherical_density.o_RASC_d AtomGrid.o_RASC_d basis_set.o_RASC_d convenience.o_RASC_d sphere_lebedev_rule.o_RASC_d scattering_factors.o_RASC_d cube.o_RASC_d fchk.o_RASC_d JKFit.o_RASC_d NoSpherA2.o_RASC_d SALTED_utilities.o_RASC_d constants.o_RASC_d SALTED_math.o_RASC_d SALTED_equicomb.o_RASC_d SALTED_io.o_RASC_d SALTED_predictor.o_RASC_d hdf5.o_RASC_d ECPs_corrections.o_RASC_d integrator.o_RASC_d
GCC_OPTS := -std=c++17 -O3 -c -fmessage-length=0 -fopenmp -static -MMD -MP
GCC_OPTS_DEBUG := -std=c++17 -Og -g -c -fmessage-length=0 -static -MMD -MP

NoSpherA2_old: ${OBJS} ${USER_OBJS}
	@echo 'Building target: $@'
	g++ -o "NoSpherA2_no_RASC" ${LIBOMP} ${OBJS} ${LIBS}
	${RM} ${OBJS}
	cp NoSpherA2 ../
	@echo 'Finished building target: $@'

# Pattern rule to compile all .cpp files
%.o: ../Src/%.cpp
	@echo 'Building file: $<'
	g++ ${GCC_OPTS} -I../Lib/OpenBlas/include_Linux -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

# Pattern rule to compile all .cpp files
%.o_RASC: ../Src/%.cpp
	@echo 'Building file: $<'
	g++ ${GCC_OPTS} -D__RASCALINE__ -I../Lib/OpenBlas/include_Linux -I../Lib/ -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

# Pattern rule to compile all .cpp files
%.o_RASC_d: ../Src/%.cpp
	@echo 'Building file: $<'
	g++ ${GCC_OPTS_DEBUG} -D__RASCALINE__ -I../Lib/ -I../Lib/OpenBlas/include_Linux -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

NoSpherA2: ${OBJS_RASC}
	@echo 'Building target: $@'
	g++ -o "NoSpherA2" ${OBJS_RASC} ${LIBOMP} ${LIBRASCALINE} ${LIBHDF5} ${LIBOPENBLAS} ${LIBS}
	${RM} ${OBJS_RASC}
	cp NoSpherA2 ../
	@echo 'Done building rascaline version of NoSpherA2'

NoSpherA2_Debug: ${OBJS_RASC_d}
	@echo 'Building target: $@'
	g++ -o "NoSpherA2_Debug" ${OBJS_RASC_d} ${LIBOMP} ${LIBRASCALINE} ${LIBHDF5} ${LIBOPENBLAS} ${LIBS}
	${RM} ${OBJS_RASC_d}
	cp NoSpherA2_Rascaline_Debug ../
	@echo 'Done building rascaline version of NoSpherA2'

clean:
	-${RM} NoSpherA2 ../NoSpherA2 NoSpherA2_Debug ../NoSpherA2_Debug ${OBJS} ${OBJS_RASC} ${OBJS_RASC_d}
	-@echo ' '

.PHONY: all
.SECONDARY:


