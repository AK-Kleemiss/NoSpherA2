cmake_minimum_required(VERSION 3.25)

if(NOT APPLE)
    message(FATAL_ERROR
        "BuildMacOSUniversal.cmake can only be run on macOS."
    )
endif()

get_filename_component(
    NOSPHERA2_SOURCE_DIR
    "${CMAKE_CURRENT_LIST_DIR}/.."
    ABSOLUTE
)

set(
    NOSPHERA2_ARM64_ENV
    "${NOSPHERA2_SOURCE_DIR}/.mambaenv/env-arm64"
)

set(
    NOSPHERA2_X86_64_ENV
    "${NOSPHERA2_SOURCE_DIR}/.mambaenv/env-x86_64"
)

function(require_micromamba_environment prefix architecture)
    if(NOT EXISTS "${prefix}/conda-meta/history")
        message(FATAL_ERROR
            "The ${architecture} micromamba environment does not exist:\n"
            "  ${prefix}\n\n"
            "Create the macOS environments first with:\n"
            "  cmake -P scripts/BootstrapMicromamba.cmake"
        )
    endif()
endfunction()

function(run_checked)
    execute_process(
        COMMAND ${ARGN}
        WORKING_DIRECTORY
            "${NOSPHERA2_SOURCE_DIR}"
        RESULT_VARIABLE
            _result
        COMMAND_ECHO
            STDOUT
    )

    if(NOT _result EQUAL 0)
        string(JOIN " " _command ${ARGN})

        message(FATAL_ERROR
            "Command failed with exit code ${_result}:\n"
            "  ${_command}"
        )
    endif()
endfunction()

require_micromamba_environment(
    "${NOSPHERA2_ARM64_ENV}"
    "arm64"
)

require_micromamba_environment(
    "${NOSPHERA2_X86_64_ENV}"
    "x86_64"
)

message(STATUS "Configuring macOS arm64 build")

run_checked(
    "${CMAKE_COMMAND}"
    --preset
    release-macos-arm64
)

message(STATUS "Building macOS arm64 target")

run_checked(
    "${CMAKE_COMMAND}"
    --build
    --preset
    release-macos-arm64
)

message(STATUS "Configuring macOS x86_64 build")

run_checked(
    "${CMAKE_COMMAND}"
    --preset
    release-macos-x86_64
)

message(STATUS "Building macOS x86_64 target")

run_checked(
    "${CMAKE_COMMAND}"
    --build
    --preset
    release-macos-x86_64
)

message(STATUS "Creating universal macOS binary")

run_checked(
    "${CMAKE_COMMAND}"
    --build
    --preset
    release-macos-universal
)

set(
    _universal_executable
    "${NOSPHERA2_SOURCE_DIR}/build/release-macos-universal/bin/NoSpherA2"
)

if(NOT EXISTS "${_universal_executable}")
    message(FATAL_ERROR
        "Universal executable was not created:\n"
        "  ${_universal_executable}"
    )
endif()

execute_process(
    COMMAND
        /usr/bin/lipo
        -archs
        "${_universal_executable}"
    OUTPUT_VARIABLE
        _architectures
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE
        _lipo_result
)

if(NOT _lipo_result EQUAL 0)
    message(FATAL_ERROR
        "Could not inspect the universal executable:\n"
        "  ${_universal_executable}"
    )
endif()

message(STATUS "Universal macOS build complete")
message(STATUS "Architectures: ${_architectures}")
message(STATUS "Executable: ${_universal_executable}")