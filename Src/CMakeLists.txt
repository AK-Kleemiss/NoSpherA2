add_subdirectory(basis_set_helper)

set(GENERATED_AUX_CPP "${CMAKE_CURRENT_BINARY_DIR}/auxiliary_basis.cpp"
        debug_utils.h)
set_source_files_properties(${GENERATED_AUX_CPP} PROPERTIES GENERATED TRUE)
add_library(libnosphera2 STATIC
    "${GENERATED_AUX_CPP}"
    "${CMAKE_CURRENT_SOURCE_DIR}/wfn_class.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/atoms.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/cif.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/properties.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/spherical_density.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/AtomGrid.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/GridManager.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/mat_nos_math.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/molecule.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/pch.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/vec_nos_math.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/properties.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/basis_set.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/convenience.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/sphere_lebedev_rule.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/scattering_factors.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/cube.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/fchk.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/JKFit.cpp"
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
add_dependencies(libnosphera2 generate_auxiliary_basis)

if (NOT (CMAKE_SYSTEM_NAME STREQUAL "Darwin"))
    if (MKL_INTERFACE STREQUAL "ilp64")
        target_compile_definitions(libnosphera2 PUBLIC USE_64BITMKL)
    endif ()
endif ()
#target_link_libraries(libnosphera2 PRIVATE "-fsanitize=address")
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

target_precompile_headers(libnosphera2 PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/pch.h
)

if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    find_package(OpenMP REQUIRED)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(OPENMP_FLAG "-fopenmp=libomp")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(OPENMP_FLAG "-fopenmp")
endif()
if (MSVC)
    target_compile_options(libnosphera2 PRIVATE $<$<CONFIG:Release>:/MT> $<$<CONFIG:Debug>:/MTd>)
endif()
target_compile_options(libnosphera2 PUBLIC
    $<$<AND:$<CXX_COMPILER_ID:GNU,Clang>,$<CONFIG:Debug>>:-g -Og>
    $<$<AND:$<CXX_COMPILER_ID:GNU,Clang>,$<CONFIG:RelWithDebInfo>>:-O2 -g>
    $<$<AND:$<CXX_COMPILER_ID:GNU,Clang>,$<CONFIG:Release>>:-O2>
    $<$<AND:$<CXX_COMPILER_ID:GNU,Clang>>:-MMD -MP -fmessage-length=0>
    $<$<AND:$<CXX_COMPILER_ID:GNU,Clang>,$<NOT:$<PLATFORM_ID:Darwin>>>:${OPENMP_FLAG}>
    $<$<PLATFORM_ID:Darwin>:-DACCELERATE_NEW_LAPACK>
)
if (MSVC)
    target_compile_options(libnosphera2 PUBLIC /openmp)
    target_compile_options(libnosphera2 PUBLIC
        $<$<CXX_COMPILER_ID:Clang>:-Xclang -fopenmp>
    )
#    target_compile_definitions(libnosphera2 PUBLIC MKL_ILP64 MKL_Direct_Call)
endif()
# This is done so Makefile can keep working the way it does. I have no Idea how to make windows work if I change the location of mdspan
#add_compile_definitions(libnosphera2 "-D__CMAKE_BUILD__")
if (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    target_link_libraries(libnosphera2 PUBLIC
        $<LINK_LIBRARY:FRAMEWORK,Accelerate>
    )
else ()
    target_compile_definitions(libnosphera2 PRIVATE __CMAKE_BUILD__  MKL_LP64 MKL_Direct_Call)
    target_include_directories(libnosphera2 PRIVATE ${MKL_INCLUDE})
endif ()
if (${USE_VECTOR_INSTRUCTIONS})
    target_compile_options(libnosphera2 PRIVATE
        $<$<AND:$<CXX_COMPILER_ID:GNU,Clang>,$<NOT:$<PLATFORM_ID:Darwin>>>:-msse2 -msse3 -msse4.1 -msse4.2 -mavx -ffast-math>
        $<$<PLATFORM_ID:Windows>:/arch:AVX>
    )
endif ()

target_link_options(libnosphera2 PRIVATE
    $<$<AND:$<NOT:$<CXX_COMPILER_ID:MSVC>>,$<LINK_LANGUAGE:CXX>>:-static>
    $<$<CXX_COMPILER_ID:GNU>:-static-libstdc++ -static-libgcc -MMD -MP -fopenmp -ffast-math>
)
target_link_libraries(libnosphera2 PUBLIC
    $<$<OR:$<PLATFORM_ID:Linux>,$<PLATFORM_ID:Windows>>:MKL::MKL>
    $<$<PLATFORM_ID:Darwin>:OpenMP::OpenMP_C OpenMP::OpenMP_CXX>
    occ::qm
    occ::main
    basis_set_helper
    mdspan
    featomic::static
    metatensor::static
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-ldl>
    $<$<CXX_COMPILER_ID:MSVC>:userenv>
    $<$<CXX_COMPILER_ID:MSVC>:ntdll>
    $<$<CXX_COMPILER_ID:MSVC>:ws2_32>
)
set(CMAKE_VS_NUGET_PACKAGE_RESTORE ON)
set_property(TARGET libnosphera2
    PROPERTY VS_PACKAGE_REFERENCES
        "intelmkl.static.win-x64_2025.3.0.453;inteltbb.redist.win_2022.2.0.503;inteltbb.devel.win_2022.2.0.503;intelopenmp.redist.win_2025.2.1.1000;intelopenmp.devel.win_2025.2.1.1000"
)
set_target_properties(libnosphera2 PROPERTIES PREFIX "")
target_include_directories(libnosphera2 PUBLIC ${MKL_INCLUDE})
add_executable(NoSpherA2 "${CMAKE_CURRENT_SOURCE_DIR}/NoSpherA2.cpp")
target_link_libraries(NoSpherA2 PRIVATE libnosphera2)
target_include_directories(NoSpherA2 PRIVATE ${MKL_INCLUDE})
if(WIN32)
    set_target_properties(libnosphera2 PROPERTIES
            VS_DEBUGGER_ENVIRONMENT "PATH=${MKL_DLL_PATH};$ENV{PATH}"
    )
    set_target_properties(NoSpherA2 PROPERTIES
            VS_DEBUGGER_ENVIRONMENT "PATH=${MKL_DLL_PATH};$ENV{PATH}"
    )
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        set(MKL_DLLS "${MKL_DLL_PATH}/libiomp5md.dll;${MKL_DLL_PATH}/tbb12_debug.dll")
    else()
        set(MKL_DLLS "${MKL_DLL_PATH}/libiomp5md.dll;${MKL_DLL_PATH}/tbb12.dll")
    endif ()
    file(COPY ${MKL_DLLS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()
if (MSVC)
    target_compile_options(NoSpherA2 PRIVATE $<$<CONFIG:Release>:/MT> $<$<CONFIG:Debug>:/MTd>)
endif()
if (LINUX)
    get_target_property(MKL_LIBRARIES MKL::MKL IMPORTED_LOCATION)
endif ()

install(TARGETS NoSpherA2 DESTINATION ${CMAKE_SOURCE_DIR} OPTIONAL)


