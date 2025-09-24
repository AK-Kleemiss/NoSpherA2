add_subdirectory(basis_set_helper)

set(GENERATED_AUX_CPP "${CMAKE_CURRENT_BINARY_DIR}/auxiliary_basis.cpp")
set_source_files_properties(${GENERATED_AUX_CPP} PROPERTIES GENERATED TRUE)
add_executable(NoSpherA2
    "${GENERATED_AUX_CPP}"
    "${CMAKE_CURRENT_SOURCE_DIR}/wfn_class.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/atoms.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/cif.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/properties.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/spherical_density.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/AtomGrid.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/basis_set.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/convenience.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/sphere_lebedev_rule.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/scattering_factors.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/cube.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/fchk.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/JKFit.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/NoSpherA2.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/SALTED_utilities.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/constants.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/ECPs_corrections.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/libCintMain.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/nos_math.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/SALTED_io.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/SALTED_predictor.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/SALTED_equicomb.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/isosurface.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/cart2sph.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/int_g2e.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/int_optimizer.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/integration_params.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/integrator.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/rys_roots.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/qct.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/b2c.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/bondwise_analysis.cpp"
)
add_dependencies(NoSpherA2 generate_auxiliary_basis)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

target_precompile_headers(NoSpherA2 PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/pch.h
)

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(OPENMP_FLAG "-fopenmp=libomp")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(OPENMP_FLAG "-fopenmp")
endif ()

target_compile_options(NoSpherA2 PRIVATE
    $<$<CONFIG:Debug>:-g -Og>
    $<$<CONFIG:Release>:-O3>
    ${OPENMP_FLAG}
    -static
    -MMD
    -MP
    -fmessage-length=0
)
# This is done so Makefile can keep working the way it does. I have no Idea how to make windows work if I change the location of mdspan
add_compile_definitions(NoSpherA2 "-D__CMAKE_BUILD__")

if (${USE_VECTOR_INSTRUCTIONS})
    target_compile_options(NoSpherA2 PRIVATE
        -msse2
        -msse3
        -msse4.1
        -msse4.2
        -mavx
        -ffast-math
    )
endif ()
target_link_options(NoSpherA2 PRIVATE
    -static-libgcc
    -static-libstdc++
)

target_link_libraries(NoSpherA2 PUBLIC
    MKL::MKL
    basis_set_helper
    mdspan
    metatensor
    featomic_static
    -ldl
    -lpthread
    -latomic
)

