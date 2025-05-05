#THIS IS THE LINUX MAKEFILE#
RM := rm -f

COMP := g++
C_COMP := gcc

LIBS := -static-libgcc -Wl,-Bstatic -static-libstdc++ -Wl,-Bstatic -ldl -lpthread
ver := $(shell expr `${COMP} -dumpversion`)
LIBOMP := ../llvm-project/openmp/build/runtime/src/libgomp.a
LIBRASCALINE := ../Lib/librascaline.a
LIBRASCALINE += ../Lib/libmetatensor.a
LIBOPENBLAS := ../OpenBLAS/installation/lib/libopenblas.a
LIBOPENBLAS_I := -I../OpenBLAS/installation/include

all: NoSpherA2 

# List of object file names without extensions
OBJ_NAMES := wfn_class atoms properties spherical_density AtomGrid basis_set convenience sphere_lebedev_rule scattering_factors cube fchk JKFit NoSpherA2 SALTED_utilities constants ECPs_corrections libCintMain nos_math SALTED_io SALTED_predictor SALTED_equicomb isosurface cart2sph int_g2e int_optimizer integration_params integrator rys_roots auxiliary_basis qct b2c bondwise_analysis

# Define the two separate lists with appropriate endings
OBJS := $(OBJ_NAMES:=.o)
OBJS_d := $(OBJ_NAMES:=.o_d)
GCC_OPTS := -std=c++2b -O3 -c -fmessage-length=0 -fopenmp -static -MMD -MP -msse2 -msse3 -msse4.1 -msse4.2 -mavx -ffast-math 
GCC_OPTS_DEBUG := -std=c++2b -Og -g -c -fmessage-length=0 -static -MMD -MP -msse2 -msse3 -msse4.1 -msse4.2 -mavx -ffast-math

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


#Basis set converter
basis_set_converter:
	@make -f makefile_basis_set_converter
	
# Precompiled header
pch.h.gch: ../Src/pch.h omp OpenBLAS
	${COMP} -std=c++23 -x c++-header -o "$@" "$<"

# Ensure auxiliary_basis.o depends on basis_set_converter to ensure proper sequence
auxiliary_basis.o: basis_set_converter
	@echo 'Compiling auxiliary_basis.o after BasisSetConverter finishes'
	@${COMP} ${GCC_OPTS} ${LIBOPENBLAS_I} -I../Lib/ -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" ../Src/auxiliary_basis.cpp
	@echo 'Finished building: auxiliary_basis.o'

# Pattern rule to compile all .cpp files
%.o: ../Src/%.cpp pch.h.gch omp OpenBLAS
	@echo 'Building file: $<'
	@${COMP} ${GCC_OPTS} ${LIBOPENBLAS_I} -I../Lib/ -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

# Ensure auxiliary_basis.o depends on basis_set_converter to ensure proper sequence
auxiliary_basis.o_d: basis_set_converter
	@echo 'Compiling auxiliary_basis.o after BasisSetConverter finishes'
	@${COMP} ${GCC_OPTS_DEBUG} ${LIBOPENBLAS_I} -I../Lib/ -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" ../Src/auxiliary_basis.cpp
	@echo 'Finished building: auxiliary_basis.o_d'

# Pattern rule to compile all .cpp files
%.o_d: ../Src/%.cpp pch.h.gch omp OpenBLAS
	@echo 'Building file: $<'
	@${COMP} ${GCC_OPTS_DEBUG} ${LIBOPENBLAS_I} -I../Lib/ -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

NoSpherA2: OpenBLAS omp ${OBJS}
	@echo 'Building target: $@'
	@${COMP} -o "NoSpherA2" ${OBJS} ${LIBOMP} ${LIBRASCALINE} ${LIBS} ${LIBOPENBLAS}
	cp NoSpherA2 ../
	@echo 'Done building NoSpherA2'

NoSpherA2_Debug: OpenBLAS omp ${OBJS_d} 
	@echo 'Building target: $@'
	@${COMP} -o "NoSpherA2_Debug" ${OBJS_d} ${LIBOMP} ${LIBRASCALINE} ${LIBS} ${LIBOPENBLAS}
	cp NoSpherA2_Debug ../
	@echo 'Done building Debug version of NoSpherA2'

clean:
	-${RM} NoSpherA2 ../NoSpherA2 NoSpherA2_Debug ../NoSpherA2_Debug ../BasisSetConverter ${OBJS} ${OBJS_d} ${BASIS_OBJS}
	-@echo ' '

.PHONY: all clean basis_set_converter
.SECONDARY:


