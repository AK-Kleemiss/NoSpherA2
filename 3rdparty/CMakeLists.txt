include(ExternalProject)

set(CMAKE_TRY_COMPILE_CONFIGURATION Release)
set(ENV{CMAKE_TRY_COMPILE_CONFIGURATION} Release)
if (NOT WIN32)
    CPMAddPackage(
        NAME Ccache.cmake
        GITHUB_REPOSITORY TheLartians/Ccache.cmake
        VERSION 1.2.5
    )
endif ()
if (BUILD_RELEASE_OCC)
    set(OCC_BUILD_MODE "Release")
else ()
    set(OCC_BUILD_MODE ${CMAKE_BUILD_TYPE})
endif ()

# Fix libcint quadmath Fortran bindings issue in debug builds with GCC
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND TARGET cint)
    target_link_libraries(cint gfortran)
endif()

include(CMakePrintHelpers)
include(${CMAKE_SOURCE_DIR}/cmake/get_all_properties.cmake)
print_target_properties(occ)

CPMAddPackage(
    NAME mdspan
    GITHUB_REPOSITORY "kokkos/mdspan"
    GIT_TAG "stable"
)

list(APPEND CMAKE_PREFIX_PATH "${metatensor_BINARY_DIR}")
CPMAddPackage(
    NAME featomic
    URL "https://github.com/metatensor/featomic/releases/download/featomic-v0.6.4/featomic-cxx-0.6.4.tar.gz"
    OPTIONS
        "CMAKE_BUILD_TYPE Release"
        "FEATOMIC_FETCH_METATENSOR ON"
        "BUILD_SHARED_LIBS OFF"
        "FEATOMIC_INSTALL_BOTH_STATIC_SHARED OFF"
        "FEATOMIC_USE_STATIC_METATENSOR ON"
#        "CMAKE_TRY_COMPILE_CONFIGURATION Release"
)
if (WIN32)
    set(featomic_st_name "featomic.lib")
else ()
    set(featomic_st_name "libfeatomic.a")
endif ()

set(MKL_INSTALL_DIR ${CMAKE_SOURCE_DIR}/Lib/MKL)
if (${FAKE_HOME_TEST})
    set(HOME "/tmp/home")
else ()
    set(HOME $ENV{HOME})
endif ()

set(MKL_LINK "static")
set(MKL_INTERFACE lp64)
if (NOT CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    include(${CMAKE_SOURCE_DIR}/cmake/micromamba.cmake)
    set(_ONEAPI_MKL_ROOT "${CMAKE_BINARY_DIR}/environments/mambaenv/Library")
    if (NOT EXISTS "${_ONEAPI_MKL_ROOT}/lib/cmake/mkl/MKLConfig.cmake")
        micromamba_environment(
            NAME mambaenv
            SPEC_FILE ${CMAKE_SOURCE_DIR}/environment.yml
        )
    endif ()
    set(MKL_THREADING intel_thread)
    set(MKL_MPI intelmpi)
    find_package(MKL CONFIG REQUIRED
        HINTS "${_ONEAPI_MKL_ROOT}/lib/cmake"
        PATH_SUFFIXES "mkl"
    )
    message(STATUS "Imported oneMKL targets: ${MKL_IMPORTED_TARGETS}")
    set(MKL_INCLUDE "${_ONEAPI_MKL_ROOT}/include" PARENT_SCOPE)
    set(MKL_DLL_PATH "${_ONEAPI_MKL_ROOT}/bin" PARENT_SCOPE)
endif ()
if (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(USE_SYSTEM_BLAS ON)
else ()
    set(USE_SYSTEM_BLAS OFF)
endif ()
if (WIN32)
    set(WIN_FLAGS "GG_NO_PRAGMA ON;CPACK_SYSTEM_NAME windows")
endif ()
if (NOT CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(MKL_FLAGS
        "BLA_VENDOR Intel10_64lp_seq"
        "BLAS_ROOT ${_ONEAPI_MKL_ROOT}/lib/cmake/mkl"
        "BLA_STATIC ON"
    )
endif ()

CPMAddPackage(
    NAME occ
    GITHUB_REPOSITORY MilitaoLucas/occ
    GIT_TAG "main_with_libs"
    OPTIONS
        "USE_QCINT OFF"
        "ENABLE_HOST_OPT OFF"
        "WITH_PYTHON_BINDINGS OFF"
        "USE_FORTRAN OFF"
        "USE_SYSTEM_BLAS ${USE_SYSTEM_BLAS}"
        "CMAKE_BUILD_TYPE ${OCC_BUILD_MODE}"
        ${MKL_FLAGS}
        ${WIN_FLAGS}
)
file(COPY ${occ_SOURCE_DIR}/share/ DESTINATION ${CMAKE_BINARY_DIR}/occ_data)
set(OCC_DATA_PATH "${CMAKE_BINARY_DIR}/occ_data" PARENT_SCOPE)

# Fix libcint quadmath Fortran bindings issue in debug builds with GCC
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND TARGET cint)
    target_link_libraries(cint)
endif()
