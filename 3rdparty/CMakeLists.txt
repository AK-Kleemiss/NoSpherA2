include(ExternalProject)

set(CMAKE_TRY_COMPILE_CONFIGURATION Release)
set(ENV{CMAKE_TRY_COMPILE_CONFIGURATION} Release)
if (MSVC)
    add_compile_options(/EHsc)
endif ()
# Ccache.cmake - include directly from submodule
if (NOT WIN32 AND NOT IS_ACTIONS)
    CPMAddPackage("gh:TheLartians/Ccache.cmake")
endif ()

if (BUILD_RELEASE_OCC)
    set(OCC_BUILD_MODE "Release")
else ()
    set(OCC_BUILD_MODE ${CMAKE_BUILD_TYPE})
endif ()

include(CMakePrintHelpers)
include(${CMAKE_SOURCE_DIR}/cmake/get_all_properties.cmake)

# featomic - from submodule
set(featomic_SOURCE_DIR ${CMAKE_SOURCE_DIR}/featomic)
set(featomic_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/featomic-build)
set(featomic_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/featomic-install)

if (NOT PATCHED_FEATOMIC)
    set(PATCHED_FEATOMIC ON CACHE INTERNAL "Whether to apply patches to featomic")
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${featomic_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/featomic-patched
    )
    execute_process(
        COMMAND git apply ${CMAKE_SOURCE_DIR}/cmake/patches/featomic-norebuild.patch
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/featomic-patched
        RESULT_VARIABLE PATCH_RESULT
    )
    set(featomic_PATCHED_SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/featomic-patched)
else()
    set(featomic_PATCHED_SOURCE_DIR ${featomic_SOURCE_DIR})
endif ()

ExternalProject_Add(featomic_external
    SOURCE_DIR ${featomic_PATCHED_SOURCE_DIR}/featomic
    BINARY_DIR ${featomic_BINARY_DIR}
    INSTALL_DIR ${featomic_INSTALL_DIR}
    CMAKE_GENERATOR Ninja
    CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=Release
        -DFEATOMIC_FETCH_METATENSOR=ON
        -DBUILD_SHARED_LIBS=OFF
        -DFEATOMIC_INSTALL_BOTH_STATIC_SHARED=OFF
        -DFEATOMIC_USE_STATIC_METATENSOR=ON
        -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
    BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --config Release
    INSTALL_COMMAND ${CMAKE_COMMAND} --install <BINARY_DIR> --config Release
)

if (NOT PATCHED_FEATOMIC)
    execute_process(
        COMMAND patch -p1 < ${CMAKE_SOURCE_DIR}/cmake/patches/metatensor-norebuild.patch
        WORKING_DIRECTORY ${featomic_BINARY_DIR}/_deps/metatensor-src
    )
endif ()

list(APPEND CMAKE_PREFIX_PATH "${featomic_BINARY_DIR}/_deps/metatensor-build")

if (WIN32)
    set(featomic_st_name "featomic.lib")
else ()
    set(featomic_st_name "libfeatomic.a")
endif ()

if (NOT (CMAKE_SYSTEM_NAME STREQUAL "Darwin"))
    set(MKL_INSTALL_DIR ${CMAKE_SOURCE_DIR}/Lib/MKL)
    set(MKL_LINK "static")
    set(MKL_INTERFACE lp64)
    include(${CMAKE_SOURCE_DIR}/cmake/micromamba.cmake)
    file(SHA1  ${CMAKE_SOURCE_DIR}/environment.yml _ENV_HASH)
    set(_SAVED_ENV_HASH "")
    if (EXISTS ${CMAKE_BINARY_DIR}/env.sha1)
        file(READ ${CMAKE_BINARY_DIR}/env.sha1 _SAVED_ENV_HASH OPTIONAL)
    endif ()
    set(_ONEAPI_MKL_ROOT "${MICROMAMBA_ENV_PATH}/Library")
    if (NOT "${_SAVED_ENV_HASH}" STREQUAL "${_ENV_HASH}")
        set(ENV_PATH "${MICROMAMBA_ENV_PATH}")
        micromamba_environment(
            NAME mambaenv
            SPEC_FILE ${CMAKE_SOURCE_DIR}/environment.yml
        )
        file(SHA1 ${CMAKE_SOURCE_DIR}/environment.yml _ENV_HASH)
        file(WRITE ${CMAKE_BINARY_DIR}/env.sha1 ${_ENV_HASH})
        # Update the path after environment creation
        if (WIN32)
            set(_ONEAPI_MKL_ROOT "${MICROMAMBA_ENV_PATH}/Library")
        else ()
            set(_ONEAPI_MKL_ROOT "${MICROMAMBA_ENV_PATH}")
        endif ()
    else()
        if (WIN32)
            set(_ONEAPI_MKL_ROOT "${MICROMAMBA_ENV_PATH}/Library")
        else ()
            set(_ONEAPI_MKL_ROOT "${MICROMAMBA_ENV_PATH}")
        endif ()
    endif ()
    set(MKL_THREADING intel_thread)
    set(MKL_MPI intelmpi)
    find_package(MKL CONFIG REQUIRED
        HINTS "${_ONEAPI_MKL_ROOT}/lib/cmake"
        PATH_SUFFIXES "mkl"
    )
    if (WIN32)
        file(COPY "${_ONEAPI_MKL_ROOT}/bin/libiomp5md.dll" DESTINATION "${CMAKE_BINARY_DIR}")
    else ()
        file(COPY "${_ONEAPI_MKL_ROOT}/lib/libiomp5.so" DESTINATION "${CMAKE_BINARY_DIR}")
    endif ()
    message(STATUS "Imported oneMKL targets: ${MKL_IMPORTED_TARGETS}")
    set(MKL_INCLUDE "${_ONEAPI_MKL_ROOT}/include" PARENT_SCOPE)
    set(MKL_DLL_PATH "${_ONEAPI_MKL_ROOT}/bin" PARENT_SCOPE)
    set(_ONEAPI_MKL_ROOT ${_ONEAPI_MKL_ROOT} PARENT_SCOPE)
    list(APPEND CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/cmake")
endif ()

if (WIN32)
    set(WIN_FLAGS
        -DGG_NO_PRAGMA=ON
        -DCPACK_SYSTEM_NAME=windows
        "-DCMAKE_CXX_FLAGS=-D__TBB_DYNAMIC_LOAD_ENABLED=0"
    )
endif ()

if (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(Darwin_FLAGS
        -DUSE_SYSTEM_BLAS=OFF
        -DUSE_OPENMP=OFF
        "-DCMAKE_CXX_FLAGS=-O2"
        "-DCMAKE_C_FLAGS=-O2"
    )
else ()
    set(LINUX_WINDOWS_COMMON_FLAGS
        -DBLA_VENDOR=Intel10_64lp
        -DUSE_SYSTEM_TBB=OFF
        -DTBB_PREFIX_PATH=${_ONEAPI_MKL_ROOT}
        -DBLAS_ROOT=${_ONEAPI_MKL_ROOT}/lib/cmake/mkl
        -DUSE_QCINT=OFF
        -DENABLE_HOST_OPT=OFF
        -DUSE_FORTRAN=OFF
        -DBLA_STATIC=ON
    )
endif ()
if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(LINUX_FLAGS
        "-DCMAKE_CXX_FLAGS=-O3 -ffast-math"
        "-DCMAKE_C_FLAGS=-O3 -ffast-math"
    )
endif ()

# occ - from submodule
set(occ_SOURCE_DIR ${CMAKE_SOURCE_DIR}/occ)
set(occ_BINARY_DIR ${CMAKE_BINARY_DIR}/_deps/occ-build)
set(occ_INSTALL_DIR ${CMAKE_BINARY_DIR}/_deps/occ-install)

ExternalProject_Add(occ_external
    SOURCE_DIR ${occ_SOURCE_DIR}
    BINARY_DIR ${occ_BINARY_DIR}
    INSTALL_DIR ${occ_INSTALL_DIR}
    CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=${OCC_BUILD_MODE}
        -DWITH_PYTHON_BINDINGS=OFF
        -DCMAKE_INSTALL_PREFIX=${occ_INSTALL_DIR}
        ${Darwin_FLAGS}
        ${WIN_FLAGS}
        ${LINUX_FLAGS}
        ${LINUX_WINDOWS_COMMON_FLAGS}
    BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --config ${OCC_BUILD_MODE}
    INSTALL_COMMAND ${CMAKE_COMMAND} --install <INSTALL_DIR> --config ${OCC_BUILD_MODE}
)
FetchContent_Declare(
    occ
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/occ
)
file(COPY ${occ_SOURCE_DIR}/share/ DESTINATION ${CMAKE_BINARY_DIR}/occ_data)
set(OCC_DATA_PATH "${CMAKE_BINARY_DIR}/occ_data" PARENT_SCOPE)

# Create imported target for occ
add_library(occ INTERFACE IMPORTED GLOBAL)
add_dependencies(occ occ_external)
target_include_directories(occ INTERFACE ${occ_SOURCE_DIR}/include)
target_include_directories(occ INTERFACE ${occ_INCLUDE_DIR})
target_link_libraries(occ INTERFACE
    occ::main
    occ::qm
    occ::core
    occ::crystal
    occ::dft
    occ::io
    occ::interaction
    occ::isosurface
    occ::slater
    occ::solvent
)
# Fix libcint quadmath Fortran bindings issue in debug builds with GCC
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND TARGET cint)
    target_link_libraries(cint gfortran)
endif()

print_target_properties(occ)
