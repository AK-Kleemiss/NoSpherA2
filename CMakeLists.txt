cmake_minimum_required(VERSION 3.25)

project(NoSpherA2 LANGUAGES C CXX)

message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

option(NOSPHERA2_BUILD_DLL "Build the NoSpherA2 DLL" OFF)
option(NOSPHERA2_BUILD_TESTS "Build NoSpherA2 tests" OFF)
option(NOSPHERA2_DEPENDENCIES_ONLY "Build dependencies without building NoSpherA2" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NO_TRIVIAL_INIT)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        add_compile_options(-ftrivial-auto-var-init=pattern)
    else()
        message(WARNING
            "NO_TRIVIAL_INIT is enabled, but "
            "${CMAKE_CXX_COMPILER_ID} does not support "
            "-ftrivial-auto-var-init=pattern"
        )
    endif()
endif()

if(WIN32)
    set(CMAKE_TRY_COMPILE_CONFIGURATION Release)
endif()

list(PREPEND CMAKE_MODULE_PATH
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake"
)

include(NoSpherA2Avx)
include(NoSpherA2Sanitizers)
include(ValidateMambaEnvironemt)
include(CopyRuntimeLibraries)

validate_micromamba_environment()


# ---------------------------------------------------------------------------
# libcint
# ---------------------------------------------------------------------------

set(WITH_FORTRAN OFF)
set(WITH_CINT2_INTERFACE OFF)
set(ENABLE_STATIC ON)
set(PYPZPX ON)
set(WITH_RANGE_COULOMB ON)
set(WITH_PYTHON_BINDINGS OFF)
set(BUILD_SHARED_LIBS OFF)

add_subdirectory(
    libcint
    "${CMAKE_BINARY_DIR}/_deps/libcint"
    EXCLUDE_FROM_ALL
)

if(NOT TARGET cint)
    message(FATAL_ERROR
        "The bundled libcint project did not define the target 'cint'"
    )
endif()

# libcint's current CMake target does not expose all source-tree public
# headers when consumed directly with add_subdirectory().
target_include_directories(
    cint
    PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/libcint/include>"
)

if(NOT TARGET libcint::libcint)
    add_library(libcint::libcint ALIAS cint)
endif()

# ---------------------------------------------------------------------------
# OCC
# ---------------------------------------------------------------------------

set(USE_SYSTEM_LIBCINT ON)
set(BUILD_FPIC ON)
set(WITH_PYTHON_BINDINGS OFF)
set(BUILD_SHARED_LIBS OFF)

add_subdirectory(
    occ
    "${CMAKE_BINARY_DIR}/_deps/occ"
    EXCLUDE_FROM_ALL
)

# ---------------------------------------------------------------------------
# mdspan
# ---------------------------------------------------------------------------

add_subdirectory(
    mdspan
    "${CMAKE_BINARY_DIR}/_deps/mdspan"
    EXCLUDE_FROM_ALL
)

# ---------------------------------------------------------------------------
# featomic
# ---------------------------------------------------------------------------

set(FEATOMIC_FETCH_METATENSOR ON CACHE BOOL
    "Fetch the matching metatensor-core source" FORCE
)
set(FEATOMIC_INSTALL_BOTH_STATIC_SHARED OFF CACHE BOOL
    "Build only the selected featomic library type" FORCE
)
set(FEATOMIC_USE_STATIC_METATENSOR ON CACHE BOOL
    "Link featomic against static metatensor-core" FORCE
)
set(BUILD_SHARED_LIBS OFF)

add_subdirectory(
    featomic/featomic
    "${CMAKE_BINARY_DIR}/_deps/featomic"
    EXCLUDE_FROM_ALL
)

# ---------------------------------------------------------------------------
# Intel MKL
# ---------------------------------------------------------------------------

set(MKL_LINK "static")
set(MKL_INTERFACE lp64)
set(MKL_THREADING intel_thread)
set(MKL_MPI intelmpi)
if(WIN32)
    set(
        MKL_DIR
        "${MICROMAMBA_ENV_PREFIX}/Library/lib/cmake/mkl"
        CACHE PATH
        "MKL CMake package directory"
    )
else()
    set(
        MKL_DIR
        "${MICROMAMBA_ENV_PREFIX}/lib/cmake/mkl"
        CACHE PATH
        "MKL CMake package directory"
    )
endif()

if(APPLE)
    find_library(
        NOSPHERA2_ACCELERATE_FRAMEWORK
        Accelerate
        REQUIRED
    )
    find_package(OpenMP REQUIRED COMPONENTS CXX)
else()
    find_package(MKL CONFIG REQUIRED)
endif()


# ---------------------------------------------------------------------------
# Basis Set Generation
# ---------------------------------------------------------------------------
add_subdirectory(
    BasisSetGenerator
    "${CMAKE_BINARY_DIR}/tools/BasisSetConverter"
)

# ---------------------------------------------------------------------------
# NoSpherA2 targets
# ---------------------------------------------------------------------------
if(NOSPHERA2_DEPENDENCIES_ONLY)
    include(cmake/InstallDependenciesOnly.cmake)
    return()
endif()

add_subdirectory(Src/core)
add_subdirectory(app)

if(NOSPHERA2_BUILD_DLL)
    add_subdirectory(Src/dll)
endif()

include(CTest)

if(NOSPHERA2_BUILD_TESTS)
    add_subdirectory(tests/src)
endif()

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    add_custom_target(
        nos-build
        DEPENDS
            NoSpherA2
    )

    if(TARGET NoSpherA2_Tests)
        add_custom_target(
            nos-test
            COMMAND
                "${CMAKE_CTEST_COMMAND}"
                --output-on-failure
            DEPENDS
                NoSpherA2_Tests
            USES_TERMINAL
        )
    endif()
endif()

if(APPLE)
    add_custom_target(
        NoSpherA2Universal

        COMMAND
            "${CMAKE_COMMAND}"
            "-DARM64_DIRECTORY=${PROJECT_SOURCE_DIR}/build/release-macos-arm64/bin"
            "-DX86_64_DIRECTORY=${PROJECT_SOURCE_DIR}/build/release-macos-x86_64/bin"
            "-DOUTPUT_DIRECTORY=${PROJECT_SOURCE_DIR}/build/release-macos-universal/bin"
            -P
            "${PROJECT_SOURCE_DIR}/cmake/CreateMacOSUniversal.cmake"

        COMMENT
            "Creating universal macOS binaries"

        VERBATIM
    )
endif()