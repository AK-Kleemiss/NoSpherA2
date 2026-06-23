include_guard(GLOBAL)

# ---------------------------------------------------------------------------
# setup_micromamba_environment()
#
# Creates a shared micromamba environment, normally under:
#
#   <project>/.mambaenv/
#
# Usage:
#
#   setup_micromamba_environment(
#       ENVIRONMENT_FILE "${CMAKE_SOURCE_DIR}/environment.yaml"
#   )
#
# Optional arguments:
#
#   PREFIX             Environment installation directory
#   ROOT_PREFIX        Micromamba package cache/root directory
#   DOWNLOAD_DIRECTORY Micromamba executable directory
#   MICROMAMBA_VERSION Version such as "2.3.0"; defaults to "latest"
#   FORCE_UPDATE       Always run micromamba install
#
# Exported cache variables:
#
#   MICROMAMBA_EXECUTABLE
#   MICROMAMBA_ROOT_PREFIX
#   MICROMAMBA_ENV_PREFIX
#
# This module intentionally does not modify PATH, CMAKE_PREFIX_PATH,
# Python3_EXECUTABLE, CARGO_EXE, RUSTC, or other build settings. These should
# be supplied through CMakePresets.json.
# ---------------------------------------------------------------------------

function(setup_micromamba_environment)
    set(options
        FORCE_UPDATE
    )

    set(one_value_arguments
        ENVIRONMENT_FILE
        PREFIX
        ROOT_PREFIX
        DOWNLOAD_DIRECTORY
        MICROMAMBA_VERSION
    )

    cmake_parse_arguments(
        MAMBA
        "${options}"
        "${one_value_arguments}"
        ""
        ${ARGN}
    )

    if(MAMBA_UNPARSED_ARGUMENTS)
        message(FATAL_ERROR
            "Unknown arguments passed to setup_micromamba_environment(): "
            "${MAMBA_UNPARSED_ARGUMENTS}"
        )
    endif()

    if(NOT MAMBA_ENVIRONMENT_FILE)
        message(FATAL_ERROR
            "setup_micromamba_environment() requires ENVIRONMENT_FILE"
        )
    endif()

    # -----------------------------------------------------------------------
    # Resolve paths
    # -----------------------------------------------------------------------

    get_filename_component(
        environment_file
        "${MAMBA_ENVIRONMENT_FILE}"
        ABSOLUTE
        BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    )

    if(NOT EXISTS "${environment_file}")
        message(FATAL_ERROR
            "Micromamba environment file does not exist:\n"
            "  ${environment_file}"
        )
    endif()

    set(default_mamba_directory
        "${CMAKE_SOURCE_DIR}/.mambaenv"
    )

    if(MAMBA_PREFIX)
        get_filename_component(
            environment_prefix
            "${MAMBA_PREFIX}"
            ABSOLUTE
            BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
        )
    else()
        set(environment_prefix
            "${default_mamba_directory}/env"
        )
    endif()

    if(MAMBA_ROOT_PREFIX)
        get_filename_component(
            micromamba_root_prefix
            "${MAMBA_ROOT_PREFIX}"
            ABSOLUTE
            BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
        )
    else()
        set(micromamba_root_prefix
            "${default_mamba_directory}/root"
        )
    endif()

    if(MAMBA_DOWNLOAD_DIRECTORY)
        get_filename_component(
            micromamba_download_directory
            "${MAMBA_DOWNLOAD_DIRECTORY}"
            ABSOLUTE
            BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
        )
    else()
        set(micromamba_download_directory
            "${default_mamba_directory}/bootstrap"
        )
    endif()

    file(MAKE_DIRECTORY
        "${default_mamba_directory}"
        "${micromamba_root_prefix}"
        "${micromamba_download_directory}"
    )

    # Prevent multiple build directories from creating/updating the same
    # environment simultaneously.
    file(
        LOCK "${default_mamba_directory}/setup.lock"
        GUARD FUNCTION
        TIMEOUT 600
        RESULT_VARIABLE lock_result
    )

    if(NOT lock_result STREQUAL "0")
        message(FATAL_ERROR
            "Could not lock the shared micromamba environment:\n"
            "  ${default_mamba_directory}/setup.lock\n"
            "Result: ${lock_result}"
        )
    endif()

    # -----------------------------------------------------------------------
    # Determine the micromamba platform identifier
    # -----------------------------------------------------------------------

    string(TOLOWER
        "${CMAKE_HOST_SYSTEM_PROCESSOR}"
        host_processor
    )

    if(WIN32)
        if(host_processor MATCHES "^(arm64|aarch64)$")
            set(micromamba_platform "win-arm64")
        else()
            set(micromamba_platform "win-64")
        endif()

    elseif(APPLE)
        if(host_processor MATCHES "^(arm64|aarch64)$")
            set(micromamba_platform "osx-arm64")
        else()
            set(micromamba_platform "osx-64")
        endif()

    elseif(CMAKE_HOST_SYSTEM_NAME STREQUAL "Linux")
        if(host_processor MATCHES "^(arm64|aarch64)$")
            set(micromamba_platform "linux-aarch64")
        elseif(host_processor MATCHES "^(ppc64le|powerpc64le)$")
            set(micromamba_platform "linux-ppc64le")
        else()
            set(micromamba_platform "linux-64")
        endif()

    else()
        message(FATAL_ERROR
            "Unsupported host platform for micromamba:\n"
            "  System: ${CMAKE_HOST_SYSTEM_NAME}\n"
            "  Processor: ${CMAKE_HOST_SYSTEM_PROCESSOR}"
        )
    endif()

    # -----------------------------------------------------------------------
    # Determine executable and archive paths
    # -----------------------------------------------------------------------

    if(WIN32)
        set(micromamba_executable
            "${micromamba_download_directory}/Library/bin/micromamba.exe"
        )
    else()
        set(micromamba_executable
            "${micromamba_download_directory}/bin/micromamba"
        )
    endif()

    set(micromamba_archive
        "${micromamba_download_directory}/micromamba.tar.bz2"
    )

    if(MAMBA_MICROMAMBA_VERSION)
        string(CONCAT micromamba_url
            "https://micro.mamba.pm/api/micromamba/"
            "${micromamba_platform}/"
            "${MAMBA_MICROMAMBA_VERSION}"
        )
    else()
        string(CONCAT micromamba_url
            "https://micro.mamba.pm/api/micromamba/"
            "${micromamba_platform}/latest"
        )
    endif()

    # -----------------------------------------------------------------------
    # Download and extract micromamba
    # -----------------------------------------------------------------------

    if(NOT EXISTS "${micromamba_executable}")
        message(STATUS
            "Downloading micromamba for ${micromamba_platform}"
        )
        message(STATUS
            "Micromamba URL: ${micromamba_url}"
        )

        file(
            DOWNLOAD
            "${micromamba_url}"
            "${micromamba_archive}"
            STATUS download_status
            LOG download_log
            SHOW_PROGRESS
            TLS_VERIFY ON
        )

        list(GET download_status 0 download_result)
        list(GET download_status 1 download_message)

        if(NOT download_result EQUAL 0)
            file(REMOVE "${micromamba_archive}")

            message(FATAL_ERROR
                "Could not download micromamba.\n"
                "URL: ${micromamba_url}\n"
                "Status: ${download_message}\n"
                "Log:\n${download_log}"
            )
        endif()

        execute_process(
            COMMAND
                "${CMAKE_COMMAND}" -E tar xjf
                "${micromamba_archive}"
            WORKING_DIRECTORY
                "${micromamba_download_directory}"
            RESULT_VARIABLE extraction_result
            OUTPUT_VARIABLE extraction_output
            ERROR_VARIABLE extraction_error
        )

        file(REMOVE "${micromamba_archive}")

        if(NOT extraction_result EQUAL 0)
            message(FATAL_ERROR
                "Could not extract micromamba.\n"
                "Output:\n${extraction_output}\n"
                "Error:\n${extraction_error}"
            )
        endif()

        if(NOT EXISTS "${micromamba_executable}")
            message(FATAL_ERROR
                "Micromamba was extracted, but the executable was not found.\n"
                "Expected:\n"
                "  ${micromamba_executable}"
            )
        endif()

        if(UNIX)
            file(
                CHMOD "${micromamba_executable}"
                PERMISSIONS
                    OWNER_READ
                    OWNER_WRITE
                    OWNER_EXECUTE
                    GROUP_READ
                    GROUP_EXECUTE
                    WORLD_READ
                    WORLD_EXECUTE
            )
        endif()
    endif()

    # -----------------------------------------------------------------------
    # Check whether environment.yaml changed
    # -----------------------------------------------------------------------

    file(
        SHA256
        "${environment_file}"
        current_environment_hash
    )

    set(environment_hash_file
        "${environment_prefix}/.environment-yaml.sha256"
    )

    set(environment_needs_update FALSE)

    if(MAMBA_FORCE_UPDATE)
        set(environment_needs_update TRUE)

    elseif(NOT EXISTS
        "${environment_prefix}/conda-meta/history"
    )
        set(environment_needs_update TRUE)

    elseif(NOT EXISTS "${environment_hash_file}")
        set(environment_needs_update TRUE)

    else()
        file(
            READ
            "${environment_hash_file}"
            installed_environment_hash
        )

        string(
            STRIP
            "${installed_environment_hash}"
            installed_environment_hash
        )

        if(NOT installed_environment_hash
           STREQUAL current_environment_hash)
            set(environment_needs_update TRUE)
        endif()
    endif()

    # -----------------------------------------------------------------------
    # Create or update the environment
    # -----------------------------------------------------------------------

    if(environment_needs_update)
        if(EXISTS "${environment_prefix}/conda-meta/history")
            set(micromamba_operation "install")
            message(STATUS
                "Updating micromamba environment: "
                "${environment_prefix}"
            )
        else()
            set(micromamba_operation "create")
            message(STATUS
                "Creating micromamba environment: "
                "${environment_prefix}"
            )
        endif()

        execute_process(
            COMMAND
                "${CMAKE_COMMAND}" -E env
                "MAMBA_ROOT_PREFIX=${micromamba_root_prefix}"
                "${micromamba_executable}"
                "${micromamba_operation}"
                --yes
                --prefix "${environment_prefix}"
                --file "${environment_file}"
            RESULT_VARIABLE environment_result
            OUTPUT_VARIABLE environment_output
            ERROR_VARIABLE environment_error
            COMMAND_ECHO STDOUT
        )

        if(NOT environment_result EQUAL 0)
            message(FATAL_ERROR
                "Micromamba failed to create/update the environment.\n"
                "Environment file:\n"
                "  ${environment_file}\n"
                "Environment prefix:\n"
                "  ${environment_prefix}\n"
                "Output:\n${environment_output}\n"
                "Error:\n${environment_error}"
            )
        endif()

        file(
            WRITE
            "${environment_hash_file}"
            "${current_environment_hash}\n"
        )
    else()
        message(STATUS
            "Micromamba environment is up to date: "
            "${environment_prefix}"
        )
    endif()

    # -----------------------------------------------------------------------
    # Export useful paths
    # -----------------------------------------------------------------------

    set(
        MICROMAMBA_EXECUTABLE
        "${micromamba_executable}"
        CACHE FILEPATH
        "Micromamba executable"
        FORCE
    )

    set(
        MICROMAMBA_ROOT_PREFIX
        "${micromamba_root_prefix}"
        CACHE PATH
        "Micromamba root prefix and package cache"
        FORCE
    )

    set(
        MICROMAMBA_ENV_PREFIX
        "${environment_prefix}"
        CACHE PATH
        "Micromamba environment prefix"
        FORCE
    )

    mark_as_advanced(
        MICROMAMBA_EXECUTABLE
        MICROMAMBA_ROOT_PREFIX
        MICROMAMBA_ENV_PREFIX
    )

    message(STATUS
        "Micromamba environment: ${MICROMAMBA_ENV_PREFIX}"
    )
endfunction()
