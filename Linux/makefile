#THIS IS THE LINUX MAKEFILE#
RM := rm -f

LIBS := -Wl,-Bstatic -lgomp -static-libgcc -static-libstdc++ -Wl,-Bstatic -llapacke -Wl,-Bstatic -llapack -Wl,-Bstatic -lgfortran -Wl,-Bstatic -lquadmath -Wl,-Bstatic -lm
LIBOMP := /usr/lib/gcc/x86_64-linux-gnu/$(shell expr `g++ -dumpversion`)/libgomp.a /usr/lib/gcc/x86_64-linux-gnu/$(shell expr `g++ -dumpversion`)/libgfortran.a /usr/lib/gcc/x86_64-linux-gnu/$(shell expr `g++ -dumpversion`)/libquadmath.a
LIBRASCALINE := ../Lib/librascaline.a
LIBRASCALINE += ../Lib/libmetatensor.a
LIBHDF5 := ../Lib/libhdf5.a
LIBOPENBLAS := ../Lib/OpenBlas/lib/libopenblas_openmp.a

all: NoSpherA2 

OBJS := wfn_class.o atoms.o properties.o spherical_density.o AtomGrid.o basis_set.o convenience.o sphere_lebedev_rule.o scattering_factors.o cube.o fchk.o JKFit.o NoSpherA2.o SALTED_utilities.o constants.o ECPs_corrections.o integrator.o hdf5.o SALTED.o math.o
OBJS_RASC_d := wfn_class.o_RASC_d atoms.o_RASC_d properties.o_RASC_d spherical_density.o_RASC_d AtomGrid.o_RASC_d basis_set.o_RASC_d convenience.o_RASC_d sphere_lebedev_rule.o_RASC_d scattering_factors.o_RASC_d cube.o_RASC_d fchk.o_RASC_d JKFit.o_RASC_d NoSpherA2.o_RASC_d SALTED_utilities.o_RASC_d constants.o_RASC_d math.o_RASC_d SALTED_equicomb.o_RASC_d SALTED_io.o_RASC_d SALTED_predictor.o_RASC_d hdf5.o_RASC_d ECPs_corrections.o_RASC_d integrator.o_RASC_d
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
	g++ ${GCC_OPTS} -D__RASCALINE__ -I../Lib/OpenBlas/include_Linux -I../Lib/ -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

# Pattern rule to compile all .cpp files
%.o_RASC_d: ../Src/%.cpp
	@echo 'Building file: $<'
	g++ ${GCC_OPTS_DEBUG} -D__RASCALINE__ -I../Lib/ -I../Lib/OpenBlas/include_Linux -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

NoSpherA2: ${OBJS}
	@echo 'Building target: $@'
	g++ -o "NoSpherA2" ${OBJS} ${LIBOMP} ${LIBRASCALINE} ${LIBHDF5} ${LIBS} ${LIBOPENBLAS}
	cp NoSpherA2 ../
	@echo 'Done building rascaline version of NoSpherA2'

NoSpherA2_Debug: ${OBJS_RASC_d}
	@echo 'Building target: $@'
	g++ -o "NoSpherA2_Debug" ${OBJS_RASC_d} ${LIBOMP} ${LIBRASCALINE} ${LIBHDF5} ${LIBS} ${LIBOPENBLAS}
	cp NoSpherA2_Rascaline_Debug ../
	@echo 'Done building rascaline version of NoSpherA2'

clean:
	-${RM} NoSpherA2 ../NoSpherA2 NoSpherA2_Debug ../NoSpherA2_Debug ${OBJS} ${OBJS_RASC_d}
	-@echo ' '

.PHONY: all
.SECONDARY:


