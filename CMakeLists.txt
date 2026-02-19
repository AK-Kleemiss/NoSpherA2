cmake_minimum_required(VERSION 3.15)
project(NoSpherA2 CXX C)
message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
set(CMAKE_CXX_STANDARD 23)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
if (NO_TRIVIAL_INIT)
    add_compile_options(-ftrivial-auto-var-init=pattern)
endif()

if (LINUX)
    add_compile_options(-msse2 -msse3 -msse4.1 -msse4.2 -mavx)
endif()

if (WIN32)
    set(CMAKE_TRY_COMPILE_CONFIGURATION Release)
    add_compile_options(/arch:AVX)
endif()
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(BUILD_DIR ${CMAKE_SOURCE_DIR}/build)
set(DEBUG_DIR ${CMAKE_SOURCE_DIR}/build_debug)
set(CMAKE_INSTALL_PREFIX $ENV{HOME} CACHE PATH "Install path prefix, prepended onto install directories.")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    option(BUILD_TESTING "Run tests" ON)
    option(BUILD_RELEASE_OCC "Build OCC in release mode even if NoSpherA2 is built in Debug" OFF)
    option(FROM_SOURCE_DEPENDENCIES "BUILD EVERY DEPENDENCY FROM SOURCE" OFF)
    option(WINDOWS_CROSS_COMPILE "CROSS COMPILING TO WINDOWS" OFF)
    option(BUILD_OPENMP "Necessary if cross compiling from windows to Linux/Darwin" OFF)
    option(USE_VECTOR_INSTRUCTIONS "If you want to use vector instructions. Defaults to ON" ON)
    option(WITH_CLANG_TIDY "Enable clang tidy" OFF)
    if (NOT (CMAKE_SYSTEM_NAME STREQUAL "Darwin"))
        set(MKL_INTERFACE "lp64" CACHE STRING "ilp64 uses double precision for indices in MKL, the other option is lp64." )
        set(ONEAPI_DIRECTORY $ENV{MKLROOT} CACHE STRING
            "Location of OneAPI. If empty and you didn't loaded modules or used a setvars.sh, OneAPI will be downloaded.")
        set(FAKE_HOME_TEST "ON" CACHE STRING
            "This is for testing the automatic install of MKL. It pretends the home directory is at /tmp/home.")
    endif ()
    set(CONDA_ENV "base" CACHE STRING "Conda env to be used.")
    set(CONDA_COMMAND "conda" CACHE STRING "CONDA command (used to download dependencies.)")
    set(OCC_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING "Build type for OCC")
    option(IS_ACTIONS "Building in GitHub Actions" OFF)
    option(USE_CUSTOM_ENV OFF)
    if (NOT CONDA_ENV STREQUAL "base")
        option(USE_CUSTOM_ENV "Use a user provided environment for conda." ON)
    endif ()
endif()
if (MSVC AND NOT ("${OCC_BUILD_TYPE}" STREQUAL "${CMAKE_BUILD_TYPE}"))
    message(FATAL_ERROR "Building OCC in a different mode than NoSpherA2 is not supported because MSVC containers have
            different sizes between Release and Debug builds.")
endif ()
#find_package(Git QUIET)
#    if(GIT_FOUND)
#      message(STATUS "Git found: initializing/updating submodules...")
#      execute_process(
#        COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
#        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
#        RESULT_VARIABLE _git_submodule_result
#        OUTPUT_QUIET ERROR_QUIET
#      )
#      if(NOT _git_submodule_result EQUAL 0)
#        message(WARNING "git submodule update --init --recursive failed with exit code ${_git_submodule_result}")
#      endif()
#    else()
#      message(STATUS "Git not found; skipping automatic git submodule initialization")
#    endif()

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
#include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/CPM.cmake")
include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/NOSTEST.cmake")
if (NOT ONLY_OCC)
    add_subdirectory(3rdparty EXCLUDE_FROM_ALL)
endif ()

# external dependencieso
set(WITH_FORTRAN "OFF")
set(WITH_CINT2_INTERFACE "OFF")
set(ENABLE_STATIC "ON")
set(BUILD_SHARED_LIBS "OFF")
set(PYPZPX "ON")
set(WITH_RANGE_COULOMB "ON")
 if (WIN32 AND ONLY_OCC)
     set(libcint_IMPORTED_LOCATION "${CMAKE_SOURCE_DIR}/Lib/LibCint/lib/cint.lib")
     if (EXISTS ${libcint_IMPORTED_LOCATION})
         message(STATUS "Using precompiled libcint from ${libcint_IMPORTED_LOCATION}")
     else ()
       message(FATAL_ERROR "Need to run '.\Windows\Build_Dependencies\build_deps.ps1 .'")
     endif ()
     add_library(cint STATIC IMPORTED)
     set_target_properties(cint PROPERTIES
         IMPORTED_LOCATION "${libcint_IMPORTED_LOCATION}"
         INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_SOURCE_DIR}/Lib/LibCint/include"
     )
 else()
     add_subdirectory(libcint)

 endif ()

#add_subdirectory(libcint)
add_library(libcint::libcint INTERFACE IMPORTED GLOBAL)
target_link_libraries(libcint::libcint INTERFACE cint)
target_include_directories(
    libcint::libcint INTERFACE
    "${CMAKE_SOURCE_DIR}/libcint/include"
    "${CMAKE_SOURCE_DIR}/libcint/src"
)
set(USE_SYSTEM_LIBCINT "ON")
set(BUILD_FPIC "ON")
set(WITH_PYTHON_BINDINGS "OFF")
set(BUILD_SHARED_LIBS "FALSE")
if (NOT (CMAKE_SYSTEM_NAME STREQUAL "Darwin"))
    if (WIN32)
        set(Python3_EXECUTABLE ${MICROMAMBA_ENV_PATH}/python.exe)
        set(Python_EXECUTABLE ${MICROMAMBA_ENV_PATH}/python.exe)
    else ()
        set(Python3_EXECUTABLE ${MICROMAMBA_ENV_PATH}/bin/python)
    endif ()
endif ()
add_subdirectory(occ)
if (ONLY_OCC)
    return()
endif ()
add_subdirectory(mdspan)
set(FEATOMIC_FETCH_METATENSOR "ON")
set(BUILD_SHARED_LIBS "OFF")
set(FEATOMIC_INSTALL_BOTH_STATIC_SHARED "OFF")
set(FEATOMIC_USE_STATIC_METATENSOR "ON")
add_subdirectory(featomic/featomic)



add_subdirectory(Src)
enable_testing()
add_subdirectory(tests)
add_subdirectory(occ_integration_tests)
