include(ExternalProject)
set(CMAKE_TRY_COMPILE_CONFIGURATION Release)
set(ENV{CMAKE_TRY_COMPILE_CONFIGURATION} Release)
add_link_options(-static)
CPMAddPackage(
    NAME Ccache.cmake
    GITHUB_REPOSITORY TheLartians/Ccache.cmake
    VERSION 1.2.5
)
if (BUILD_RELEASE_OCC)
    set(OCC_BUILD_MODE "Release")
else ()
    set(OCC_BUILD_MODE ${CMAKE_BUILD_TYPE})
endif ()


#ExternalProject_Add(occ
#    SOURCE_DIR /home/lucas/CLionProjects/occ
##    CMAKE_ARGS -DUSE_QCINT=ON -DENABLE_HOST_OPT=ON
#)
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
        "CMAKE_TRY_COMPILE_CONFIGURATION Release"
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

if (LINUX)
    set(MKL_LINK "static")
    if (NOT EXISTS ${ONEAPI_DIRECTORY})
        set(_MKL_ROOT "${MKL_INSTALL_DIR}/mkl/latest")
        set(_MKL_CONFIG "${_MKL_ROOT}/lib/cmake/mkl/MKLConfig.cmake")

        if(NOT EXISTS "${_MKL_CONFIG}")
            file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/Lib/MKL")
            set(_MKL_SH "${CMAKE_SOURCE_DIR}/Lib/intel-mkl-installer.sh")
            if (NOT EXISTS ${_MKL_SH})
                file(DOWNLOAD
                  "https://registrationcenter-download.intel.com/akdlm/IRC_NAS/47c7d946-fca1-441a-b0df-b094e3f045ea/intel-onemkl-2025.2.0.629_offline.sh"
                  "${_MKL_SH}"
                  EXPECTED_HASH SHA256=7bfdde1379a9cd3c2784af4352931a0527a5c483c9da71e851a3f06e055af7c7
                  SHOW_PROGRESS
                  STATUS _DL_STATUS
                )
                list(GET _DL_STATUS 0 _DL_CODE)
                if(NOT _DL_CODE EQUAL 0)
                  message(FATAL_ERROR "Failed to download oneMKL: ${_DL_STATUS}")
                endif()
            endif ()
            set(OLD_HOME $ENV{HOME})
            set(ENV{HOME} "/tmp/home")
            message(STATUS "INSTALLING MKL IN: ${MKL_INSTALL_DIR}")
            file(MAKE_DIRECTORY $ENV{HOME})
            execute_process(COMMAND chmod +x ${_MKL_SH})
            execute_process(
                COMMAND ${_MKL_SH} -a --silent --eula accept --install-dir "${MKL_INSTALL_DIR}"

              RESULT_VARIABLE _MKL_RC
            )
            FILE(REMOVE $ENV{HOME})
            set(ENV{HOME} ${OLD_HOME})
            if(NOT _MKL_RC EQUAL 0)
              message(FATAL_ERROR "oneMKL installer failed with code ${_MKL_RC}")
            endif()
        endif()
        list(APPEND CMAKE_PREFIX_PATH "${_MKL_ROOT}")
        find_package(MKL CONFIG REQUIRED
            HINTS "${_MKL_ROOT}"
            PATH_SUFFIXES "lib/cmake/mkl"
        )
        message(STATUS "Imported oneMKL targets: ${MKL_IMPORTED_TARGETS}")

    else ()
        set(_ONEAPI_MKL_ROOT "${ONEAPI_DIRECTORY}/mkl/latest")
        set(MKL_INTERFACE "lp64")
#        set(MKL_THREADING "tbb_thread")

        find_package(MKL CONFIG REQUIRED
            HINTS "${_ONEAPI_MKL_ROOT}" "${ONEAPI_DIRECTORY}"
            PATH_SUFFIXES "lib/cmake/mkl"
        )
        print_target_properties(MKL::MKL)
        get_target_property(MKL_LINK_LIBRARIES MKL::MKL INTERFACE_LINK_LIBRARIES)
        set(MKL_LINK_LIBRARIES ${MKL_LINK_LIBRARIES} PARENT_SCOPE)
        message(STATUS "Imported oneMKL targets: ${MKL_IMPORTED_TARGETS}")
    endif ()
elseif (WIN32)

    include("${CMAKE_SOURCE_DIR}/cmake/NuGetTools.cmake")
    nuget_initialize()
    nuget_add_dependencies(
        PACKAGE intelmkl.devel.win-x64 VERSION 2025.3.0.453 CMAKE_PREFIX_PATHS installed/x64-windows
    )
    find_package(MKL CONFIG REQUIRED)
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
        "USE_SYSTEM_BLAS ON"
        "USE_SYSTEM_TBB OFF"
        "CMAKE_BUILD_TYPE ${OCC_BUILD_MODE}"
        "BLA_VENDOR Intel10_64lp_seq"
        "BLA_STATIC ON"
        "BLAS_ROOT ${_ONEAPI_MKL_ROOT}" #${CMAKE_SOURCE_DIR}/Lib/MKL/mkl/latest"
)

# Fix libcint quadmath Fortran bindings issue in debug builds with GCC
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND TARGET cint)
    target_link_libraries(cint)
endif()
