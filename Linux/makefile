#THIS IS THE LINUX MAKEFILE#
RM := rm -f

LIBS := -Wl,-Bstatic -lgomp -static-libgcc -static-libstdc++ -Wl,-Bstatic -llapacke -Wl,-Bstatic -llapack -Wl,-Bstatic -lgfortran -Wl,-Bstatic -lquadmath -Wl,-Bstatic -lm
LIBOMP := /usr/lib/gcc/x86_64-linux-gnu/$(shell expr `g++ -dumpversion`)/libgomp.a /usr/lib/gcc/x86_64-linux-gnu/$(shell expr `g++ -dumpversion`)/libgfortran.a /usr/lib/gcc/x86_64-linux-gnu/$(shell expr `g++ -dumpversion`)/libquadmath.a
LIBRASCALINE := ../Lib/librascaline.a
LIBRASCALINE += ../Lib/libmetatensor.a
LIBHDF5 := ../Lib/libhdf5.a
LIBOPENBLAS := ../Lib/OpenBlas/lib/libopenblas_openmp.a

all: NoSpherA2 

# List of object file names without extensions
OBJ_NAMES := wfn_class atoms properties spherical_density AtomGrid basis_set convenience sphere_lebedev_rule scattering_factors cube fchk JKFit NoSpherA2 SALTED_utilities constants ECPs_corrections libCintMain hdf5 math SALTED_io SALTED_predictor SALTED_equicomb isosurface cart2sph int_g2e int_optimizer integration_params integrator rys_roots

# Define the two separate lists with appropriate endings
OBJS := $(OBJ_NAMES:=.o)
OBJS_d := $(OBJ_NAMES:=.o_RASC_d)
GCC_OPTS := -std=c++20 -O3 -c -fmessage-length=0 -fopenmp -static -MMD -MP
GCC_OPTS_DEBUG := -std=c++20 -Og -g -c -fmessage-length=0 -static -MMD -MP

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

NoSpherA2_Debug: ${OBJS_d}
	@echo 'Building target: $@'
	g++ -o "NoSpherA2_Debug" ${OBJS_d} ${LIBOMP} ${LIBRASCALINE} ${LIBHDF5} ${LIBS} ${LIBOPENBLAS}
	cp NoSpherA2_Rascaline_Debug ../
	@echo 'Done building rascaline version of NoSpherA2'

clean:
	-${RM} NoSpherA2 ../NoSpherA2 NoSpherA2_Debug ../NoSpherA2_Debug ${OBJS} ${OBJS_d}
	-@echo ' '

.PHONY: all
.SECONDARY:


