# FindHomebrewOpenMP.cmake
# Find libomp installed via Homebrew on macOS (supports both Intel and Apple Silicon)
# For universal binaries, this creates a fat library combining both architectures.
#
# This module sets the following variables:
#   HOMEBREW_OPENMP_FOUND        - True if libomp was found
#   HOMEBREW_OPENMP_INCLUDE_DIRS - Include directories for libomp
#   HOMEBREW_OPENMP_LIBRARIES    - Libraries to link against
#   HOMEBREW_OPENMP_FLAGS        - Compiler flags for OpenMP
#
# This module creates the following imported targets:
#   Homebrew::OpenMP             - Imported target for libomp

include(FindPackageHandleStandardArgs)

# Homebrew paths for each architecture
set(_HOMEBREW_ARM64_PREFIX "/opt/homebrew")
set(_HOMEBREW_X86_64_PREFIX "/usr/local")

# Check if we're building a universal binary
set(_BUILDING_UNIVERSAL OFF)
if(CMAKE_OSX_ARCHITECTURES)
    list(LENGTH CMAKE_OSX_ARCHITECTURES _ARCH_COUNT)
    if(_ARCH_COUNT GREATER 1)
        set(_BUILDING_UNIVERSAL ON)
    endif()
endif()

message(STATUS "Building universal binary: ${_BUILDING_UNIVERSAL}")

if(_BUILDING_UNIVERSAL)
    # Universal binary: need both ARM64 and x86_64 versions
    message(STATUS "Looking for Homebrew libomp for universal binary...")

    # Find ARM64 libomp
    find_path(HOMEBREW_OPENMP_INCLUDE_DIR_ARM64
        NAMES omp.h
        PATHS "${_HOMEBREW_ARM64_PREFIX}/opt/libomp/include"
        NO_DEFAULT_PATH
    )
    find_library(HOMEBREW_OPENMP_LIBRARY_ARM64
        NAMES omp
        PATHS "${_HOMEBREW_ARM64_PREFIX}/opt/libomp/lib"
        NO_DEFAULT_PATH
    )

    # Find x86_64 libomp
    find_path(HOMEBREW_OPENMP_INCLUDE_DIR_X86_64
        NAMES omp.h
        PATHS "${_HOMEBREW_X86_64_PREFIX}/opt/libomp/include"
        NO_DEFAULT_PATH
    )
    find_library(HOMEBREW_OPENMP_LIBRARY_X86_64
        NAMES omp
        PATHS "${_HOMEBREW_X86_64_PREFIX}/opt/libomp/lib"
        NO_DEFAULT_PATH
    )

    # Use ARM64 include dir (headers are architecture-independent)
    set(HOMEBREW_OPENMP_INCLUDE_DIR "${HOMEBREW_OPENMP_INCLUDE_DIR_ARM64}")

    # Create a universal (fat) library if both are found
    if(HOMEBREW_OPENMP_LIBRARY_ARM64 AND HOMEBREW_OPENMP_LIBRARY_X86_64)
        set(_FAT_LIBOMP "${CMAKE_BINARY_DIR}/libomp_universal.dylib")

        # Check if fat library already exists and is up to date
        set(_NEED_LIPO ON)
        if(EXISTS "${_FAT_LIBOMP}")
            # Check if source libraries are newer
            file(TIMESTAMP "${_FAT_LIBOMP}" _FAT_TIME)
            file(TIMESTAMP "${HOMEBREW_OPENMP_LIBRARY_ARM64}" _ARM64_TIME)
            file(TIMESTAMP "${HOMEBREW_OPENMP_LIBRARY_X86_64}" _X86_64_TIME)
            if(_FAT_TIME GREATER _ARM64_TIME AND _FAT_TIME GREATER _X86_64_TIME)
                set(_NEED_LIPO OFF)
            endif()
        endif()

        if(_NEED_LIPO)
            message(STATUS "Creating universal libomp with lipo...")
            execute_process(
                COMMAND lipo -create
                    "${HOMEBREW_OPENMP_LIBRARY_ARM64}"
                    "${HOMEBREW_OPENMP_LIBRARY_X86_64}"
                    -output "${_FAT_LIBOMP}"
                RESULT_VARIABLE _LIPO_RESULT
                ERROR_VARIABLE _LIPO_ERROR
            )
            if(NOT _LIPO_RESULT EQUAL 0)
                message(FATAL_ERROR "Failed to create universal libomp: ${_LIPO_ERROR}")
            endif()
        endif()

        set(HOMEBREW_OPENMP_LIBRARY "${_FAT_LIBOMP}")
        message(STATUS "Created universal libomp at: ${_FAT_LIBOMP}")
    else()
        message(WARNING "Could not find both ARM64 and x86_64 libomp for universal binary")
        message(STATUS "  ARM64: ${HOMEBREW_OPENMP_LIBRARY_ARM64}")
        message(STATUS "  x86_64: ${HOMEBREW_OPENMP_LIBRARY_X86_64}")
    endif()

    mark_as_advanced(
        HOMEBREW_OPENMP_INCLUDE_DIR_ARM64
        HOMEBREW_OPENMP_INCLUDE_DIR_X86_64
        HOMEBREW_OPENMP_LIBRARY_ARM64
        HOMEBREW_OPENMP_LIBRARY_X86_64
    )
else()
    # Single architecture build: detect based on current architecture
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
        set(_HOMEBREW_PREFIX "${_HOMEBREW_ARM64_PREFIX}")
    else()
        set(_HOMEBREW_PREFIX "${_HOMEBREW_X86_64_PREFIX}")
    endif()

    # Allow user override via environment variable or CMake variable
    if(DEFINED ENV{HOMEBREW_PREFIX})
        set(_HOMEBREW_PREFIX "$ENV{HOMEBREW_PREFIX}")
    endif()
    if(DEFINED HOMEBREW_PREFIX)
        set(_HOMEBREW_PREFIX "${HOMEBREW_PREFIX}")
    endif()

    message(STATUS "Looking for Homebrew libomp in: ${_HOMEBREW_PREFIX}")

    # Find libomp include directory
    find_path(HOMEBREW_OPENMP_INCLUDE_DIR
        NAMES omp.h
        PATHS
            "${_HOMEBREW_PREFIX}/opt/libomp/include"
            "${_HOMEBREW_PREFIX}/include"
        NO_DEFAULT_PATH
    )

    # Find libomp library
    find_library(HOMEBREW_OPENMP_LIBRARY
        NAMES omp
        PATHS
            "${_HOMEBREW_PREFIX}/opt/libomp/lib"
            "${_HOMEBREW_PREFIX}/lib"
        NO_DEFAULT_PATH
    )
endif()

# Handle standard find_package arguments
find_package_handle_standard_args(HomebrewOpenMP
    REQUIRED_VARS
        HOMEBREW_OPENMP_INCLUDE_DIR
        HOMEBREW_OPENMP_LIBRARY
    FAIL_MESSAGE
        "Could not find Homebrew libomp. Please install it via 'brew install libomp'"
)

if(HOMEBREW_OPENMP_FOUND)
    set(HOMEBREW_OPENMP_INCLUDE_DIRS "${HOMEBREW_OPENMP_INCLUDE_DIR}")
    set(HOMEBREW_OPENMP_LIBRARIES "${HOMEBREW_OPENMP_LIBRARY}")
    set(HOMEBREW_OPENMP_FLAGS "-Xpreprocessor -fopenmp")

    # Create imported target if it doesn't exist
    if(NOT TARGET Homebrew::OpenMP)
        add_library(Homebrew::OpenMP SHARED IMPORTED)
        set_target_properties(Homebrew::OpenMP PROPERTIES
            IMPORTED_LOCATION "${HOMEBREW_OPENMP_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${HOMEBREW_OPENMP_INCLUDE_DIR}"
            INTERFACE_COMPILE_OPTIONS "-Xpreprocessor;-fopenmp"
        )
    endif()

    message(STATUS "Found Homebrew libomp:")
    message(STATUS "  Include dir: ${HOMEBREW_OPENMP_INCLUDE_DIR}")
    message(STATUS "  Library: ${HOMEBREW_OPENMP_LIBRARY}")
endif()

mark_as_advanced(
    HOMEBREW_OPENMP_INCLUDE_DIR
    HOMEBREW_OPENMP_LIBRARY
)
