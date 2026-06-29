cmake_minimum_required(VERSION 3.25)

# Usage:
#
#   cmake -P scripts/SetupVSEnvironment.cmake
#
# Optional overrides:
#
#   cmake
#       -DSOURCE_DIR=C:/path/to/NoSpherA2
#       -DRELEASE_PRESET=release-windows
#       -DDEBUG_PRESET=debug-windows
#       -P scripts/SetupVSEnvironment.cmake

if(NOT WIN32)
    message(FATAL_ERROR
        "SetupVSEnvironment.cmake is intended for Windows."
    )
endif()

# Assume this script is located in <repository>/cmake/.
if(NOT DEFINED SOURCE_DIR)
    get_filename_component(
        SOURCE_DIR
        "${CMAKE_CURRENT_LIST_DIR}/.."
        ABSOLUTE
    )
endif()

if(NOT DEFINED RELEASE_PRESET)
    set(RELEASE_PRESET "release-windows")
endif()

if(NOT DEFINED DEBUG_PRESET)
    set(DEBUG_PRESET "debug-windows")
endif()

if(NOT DEFINED RELEASE_BUILD_DIR)
    set(
        RELEASE_BUILD_DIR
        "${SOURCE_DIR}/build/release-windows"
    )
endif()

if(NOT DEFINED DEBUG_BUILD_DIR)
    set(
        DEBUG_BUILD_DIR
        "${SOURCE_DIR}/build/debug-windows"
    )
endif()

if(NOT DEFINED RELEASE_INSTALL_DIR)
    set(
        RELEASE_INSTALL_DIR
        "${SOURCE_DIR}/deps-install-release"
    )
endif()

if(NOT DEFINED DEBUG_INSTALL_DIR)
    set(
        DEBUG_INSTALL_DIR
        "${SOURCE_DIR}/deps-install-debug"
    )
endif()

# Number of parallel build jobs. An empty value lets CMake choose.
if(NOT DEFINED PARALLEL_JOBS)
    set(PARALLEL_JOBS "")
endif()


function(run_checked)
    execute_process(
        COMMAND
            ${ARGN}
        WORKING_DIRECTORY
            "${SOURCE_DIR}"
        COMMAND_ECHO
            STDOUT
        RESULT_VARIABLE
            result
    )

    if(NOT result EQUAL 0)
        string(JOIN " " command_text ${ARGN})

        message(FATAL_ERROR
            "Command failed with exit code ${result}:\n"
            "${command_text}"
        )
    endif()
endfunction()


function(build_dependency_configuration)
    set(options)
    set(one_value_arguments
        NAME
        PRESET
        CONFIG
        BUILD_DIR
        INSTALL_DIR
    )

    cmake_parse_arguments(
        DEP
        "${options}"
        "${one_value_arguments}"
        ""
        ${ARGN}
    )

    foreach(required_argument
        NAME
        PRESET
        CONFIG
        BUILD_DIR
        INSTALL_DIR
    )
        if(NOT DEP_${required_argument})
            message(FATAL_ERROR
                "build_dependency_configuration() is missing "
                "${required_argument}"
            )
        endif()
    endforeach()

    message(STATUS "")
    message(STATUS "========================================")
    message(STATUS "Building ${DEP_NAME} dependencies")
    message(STATUS "========================================")
    message(STATUS "Preset:       ${DEP_PRESET}")
    message(STATUS "Configuration:${DEP_CONFIG}")
    message(STATUS "Build tree:   ${DEP_BUILD_DIR}")
    message(STATUS "Install tree: ${DEP_INSTALL_DIR}")
    message(STATUS "")

    # Configure the dependency-only build.
    run_checked(
        "${CMAKE_COMMAND}"
        "--preset"
        "${DEP_PRESET}"
        "-DNOSPHERA2_DEPENDENCIES_ONLY=ON"
    )

    # Build.
    if(PARALLEL_JOBS)
        run_checked(
            "${CMAKE_COMMAND}"
            "--build"
            "--preset"
            "${DEP_PRESET}"
            "--parallel"
            "${PARALLEL_JOBS}"
        )
    else()
        run_checked(
            "${CMAKE_COMMAND}"
            "--build"
            "--preset"
            "${DEP_PRESET}"
            "--parallel"
        )
    endif()

    # Install the selected configuration.
    run_checked(
        "${CMAKE_COMMAND}"
        "--install"
        "${DEP_BUILD_DIR}"
        "--config"
        "${DEP_CONFIG}"
        "--prefix"
        "${DEP_INSTALL_DIR}"
    )

    message(STATUS "")
    message(STATUS
        "${DEP_NAME} dependencies installed to:\n"
        "  ${DEP_INSTALL_DIR}"
    )
endfunction()


build_dependency_configuration(
    NAME
        "Release"
    PRESET
        "${RELEASE_PRESET}"
    CONFIG
        "Release"
    BUILD_DIR
        "${RELEASE_BUILD_DIR}"
    INSTALL_DIR
        "${RELEASE_INSTALL_DIR}"
)

build_dependency_configuration(
    NAME
        "Debug"
    PRESET
        "${DEBUG_PRESET}"
    CONFIG
        "Debug"
    BUILD_DIR
        "${DEBUG_BUILD_DIR}"
    INSTALL_DIR
        "${DEBUG_INSTALL_DIR}"
)

message(STATUS "")
message(STATUS "========================================")
message(STATUS "All Windows dependencies were installed")
message(STATUS "========================================")
message(STATUS "Release: ${RELEASE_INSTALL_DIR}")
message(STATUS "Debug:   ${DEBUG_INSTALL_DIR}")