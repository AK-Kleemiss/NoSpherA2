#THIS IS THE LINUX MAKEFILE#
RM := rm -f

COMP := g++
C_COMP := gcc

LIBS := -static-libgcc -Wl,-Bstatic -static-libstdc++ -Wl,-Bstatic
ver := $(shell expr `${COMP} -dumpversion`)
LIBOMP := ../llvm-project/openmp/build/runtime/src/libgomp.a
LIBRASCALINE := ../Lib/librascaline.a
LIBRASCALINE += ../Lib/libmetatensor.a
LIBHDF5 := ../Lib/libhdf5.a
LIBOPENBLAS := ../OpenBLAS/installation/lib/libopenblas.a
LIBOPENBLAS_I := -I../OpenBLAS/installation/include

all: NoSpherA2 

# List of object file names without extensions
OBJ_NAMES := wfn_class atoms properties spherical_density AtomGrid basis_set convenience sphere_lebedev_rule scattering_factors cube fchk JKFit NoSpherA2 SALTED_utilities constants ECPs_corrections libCintMain hdf5 nos_math SALTED_io SALTED_predictor SALTED_equicomb isosurface cart2sph int_g2e int_optimizer integration_params integrator rys_roots

# Define the two separate lists with appropriate endings
OBJS := $(OBJ_NAMES:=.o)
OBJS_d := $(OBJ_NAMES:=.o_d)
GCC_OPTS := -std=c++2b -O3 -c -fmessage-length=0 -fopenmp -static -MMD -MP
GCC_OPTS_DEBUG := -std=c++2b -Og -g -c -fmessage-length=0 -static -MMD -MP

omp: 
	@if [ ! -f $(LIBOMP) ]; then \
		echo 'Building OpenMP, since $(LIBOMP) doesnt exist'; \
		cd ../llvm-project/openmp && mkdir -p build && cd build && cmake .. -DCMAKE_C_COMPILER=${C_COMP} -DCMAKE_CXX_COMPILER=${COMP} -DCMAKE_BUILD_TYPE=Release -DLIBOMP_ENABLE_SHARED=OFF && make -j; \
	else \
		echo 'Skipping OpenMP build, $(LIBOMP) already exists'; \
	fi

OpenBLAS: 
	@if [ ! -f $(LIBOPENBLAS) ]; then \
		echo 'Building OpenBLAS, since $(LIBOPENBLAS) doesnt exist'; \
		cd ../OpenBLAS && make -j && ${RM} installation && make install PREFIX=installation; \
	else \
		echo 'Skipping OpenBLAS build, $(LIBOPENBLAS) already exists'; \
	fi

# Precompiled header
pch.h.gch: ../Src/pch.h omp OpenBLAS
	${COMP} -std=c++23 -x c++-header -o "$@" "$<"

# Pattern rule to compile all .cpp files
%.o: ../Src/%.cpp pch.h.gch omp OpenBLAS
	@echo 'Building file: $<'
	${COMP} ${GCC_OPTS} ${LIBOPENBLAS_I} -I../Lib/ -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

# Pattern rule to compile all .cpp files
%.o_d: ../Src/%.cpp pch.h.gch omp OpenBLAS
	@echo 'Building file: $<'
	${COMP} ${GCC_OPTS_DEBUG} ${LIBOPENBLAS_I} -I../Lib/ -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

NoSpherA2: ${OBJS} omp OpenBLAS
	@echo 'Building target: $@'
	${COMP} -o "NoSpherA2" ${OBJS} ${LIBOMP} ${LIBRASCALINE} ${LIBHDF5} ${LIBS} ${LIBOPENBLAS}
	cp NoSpherA2 ../
	@echo 'Done building NoSpherA2'

NoSpherA2_Debug: ${OBJS_d} omp OpenBLAS
	@echo 'Building target: $@'
	${COMP} -o "NoSpherA2_Debug" ${OBJS_d} ${LIBOMP} ${LIBRASCALINE} ${LIBHDF5} ${LIBS} ${LIBOPENBLAS}
	cp NoSpherA2_Debug ../
	@echo 'Done building Debug version of NoSpherA2'

clean:
	-${RM} NoSpherA2 ../NoSpherA2 NoSpherA2_Debug ../NoSpherA2_Debug ${OBJS} ${OBJS_d}
	-@echo ' '

.PHONY: all clean
.SECONDARY:


