cmake_minimum_required(VERSION 3.25)

project(NoSpherA2 LANGUAGES C CXX)

message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

option(NOSPHERA2_BUILD_DLL "Build the NoSpherA2 DLL" OFF)
option(NOSPHERA2_BUILD_TESTS "Build NoSpherA2 tests" OFF)

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
include(NoSpherA2Deps)

nosphera2_setup_dependencies()


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

find_package(MKL CONFIG REQUIRED)


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

function(nosphera2_copy_runtime_libraries target)
    if(NOT TARGET "${target}")
        message(FATAL_ERROR
            "nosphera2_copy_runtime_libraries(): target '${target}' does not exist"
        )
    endif()

    if(WIN32)
        add_custom_command(
            TARGET "${target}"
            POST_BUILD
            COMMAND
                "${CMAKE_COMMAND}" -E copy_if_different
                $<TARGET_RUNTIME_DLLS:${target}>
                "$<TARGET_FILE_DIR:${target}>"
            COMMAND_EXPAND_LISTS
            VERBATIM
        )

        find_file(
            NOSPHERA2_INTEL_OPENMP_DLL
            NAMES
                libiomp5md.dll
                libiomp5.dll
            PATHS
                "${MICROMAMBA_ENV_PREFIX}/Library/bin"
                "${MICROMAMBA_ENV_PREFIX}/bin"
            NO_DEFAULT_PATH
            REQUIRED
        )

        add_custom_command(
            TARGET "${target}"
            POST_BUILD
            COMMAND
                "${CMAKE_COMMAND}" -E copy_if_different
                "${NOSPHERA2_INTEL_OPENMP_DLL}"
                "$<TARGET_FILE_DIR:${target}>"
            COMMENT
                "Copying Intel OpenMP runtime for ${target}"
            VERBATIM
        )

    elseif(APPLE)
        set_target_properties(
            "${target}"
            PROPERTIES
                BUILD_RPATH "${MICROMAMBA_ENV_PREFIX}/lib"
                INSTALL_RPATH "@loader_path"
                INSTALL_RPATH_USE_LINK_PATH TRUE
        )

    elseif(UNIX)
        set_target_properties(
            "${target}"
            PROPERTIES
                BUILD_RPATH "${MICROMAMBA_ENV_PREFIX}/lib"
                INSTALL_RPATH "$ORIGIN"
                INSTALL_RPATH_USE_LINK_PATH TRUE
        )
    endif()
endfunction()

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

