#THIS IS THE LINUX MAKEFILE#
RM := rm -f

COMP := g++
C_COMP := gcc

LIBS := -Wl,-Bstatic -lstdc++ -lpthread -Wl,-Bdynamic -lm -lc
ver := $(shell expr `${COMP} -dumpversion`)
LIBOMP := ../Lib/OpenMP_build/lib/libomp.a
LIBRASCALINE := ../Lib/featomic_install/lib/libfeatomic.a
LIBRASCALINE += ../Lib/featomic_install/lib/libmetatensor.a
LIBRASCALINE_I := ../Lib/featomic_install/include
LIBOPENBLAS := ../Lib/OpenBLAS_build/lib/libopenblas.a
LIBOPENBLAS_I := -I../Lib/OpenBLAS_build/include

all: NoSpherA2 

# List of object file names without extensions
OBJ_NAMES := wfn_class atoms properties spherical_density AtomGrid basis_set convenience sphere_lebedev_rule scattering_factors cube fchk JKFit NoSpherA2 SALTED_utilities constants ECPs_corrections libCintMain nos_math SALTED_io SALTED_predictor SALTED_equicomb isosurface cart2sph int_g2e int_optimizer integration_params integrator rys_roots auxiliary_basis qct b2c bondwise_analysis

# Define the two separate lists with appropriate endings
OBJS := $(OBJ_NAMES:=.o)
OBJS_d := $(OBJ_NAMES:=.o_d)
GCC_OPTS := -std=c++2b -O3 -c -fmessage-length=0 -fopenmp -static -MMD -MP -msse2 -msse3 -msse4.1 -msse4.2 -mavx -ffast-math 
GCC_OPTS_DEBUG := -std=c++2b -Og -g -c -fmessage-length=0 -static -MMD -MP -msse2 -msse3 -msse4.1 -msse4.2 -mavx -ffast-math

#Basis set converter
basis_set_converter:
	@make -f makefile_basis_set_converter
	
# Precompiled header
pch.h.gch: ../Src/pch.h
	@${COMP} ${GCC_OPTS} -x c++-header -o "$@" "$<"

# Precompiled header for debug
pch.h.gch_d: ../Src/pch.h
	@${COMP} ${GCC_OPTS_DEBUG} -x c++-header -o "$@" "$<"

# Ensure auxiliary_basis.o depends on basis_set_converter to ensure proper sequence
auxiliary_basis.o: basis_set_converter
	@echo 'Compiling auxiliary_basis.o after BasisSetConverter finishes'
	@${COMP} ${GCC_OPTS} ${LIBOPENBLAS_I} -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" ../Src/auxiliary_basis.cpp
	@echo 'Finished building: auxiliary_basis.o'

# Pattern rule to compile all .cpp files
%.o: ../Src/%.cpp pch.h.gch
	@echo 'Building file: $<'
	@${COMP} ${GCC_OPTS} ${LIBOPENBLAS_I} -I${LIBRASCALINE_I} -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

# Ensure auxiliary_basis.o depends on basis_set_converter to ensure proper sequence
auxiliary_basis.o_d: basis_set_converter
	@echo 'Compiling auxiliary_basis.o after BasisSetConverter finishes'
	@${COMP} ${GCC_OPTS_DEBUG} ${LIBOPENBLAS_I} -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" ../Src/auxiliary_basis.cpp
	@echo 'Finished building: auxiliary_basis.o_d'

# Pattern rule to compile all .cpp files
%.o_d: ../Src/%.cpp pch.h.gch_d
	@echo 'Building file: $<'
	@${COMP} ${GCC_OPTS_DEBUG} ${LIBOPENBLAS_I} -I${LIBRASCALINE_I} -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'

NoSpherA2: ${OBJS}
	@echo 'Building target: $@'
	@${COMP} -o "NoSpherA2" ${OBJS} ${LIBOMP} ${LIBRASCALINE} ${LIBOPENBLAS} ${LIBS}
	cp NoSpherA2 ../
	@echo 'Done building NoSpherA2'

NoSpherA2_Debug: ${OBJS_d} 
	@echo 'Building target: $@'
	@${COMP} -o "NoSpherA2_Debug" ${OBJS_d} ${LIBOMP} ${LIBRASCALINE} ${LIBOPENBLAS} ${LIBS}
	cp NoSpherA2_Debug ../
	@echo 'Done building Debug version of NoSpherA2'

clean:
	-${RM} NoSpherA2 ../NoSpherA2 NoSpherA2_Debug ../NoSpherA2_Debug ../BasisSetConverter ${OBJS} ${OBJS_d} ${BASIS_OBJS}
	-@echo ' '

.PHONY: all clean basis_set_converter
.SECONDARY:


