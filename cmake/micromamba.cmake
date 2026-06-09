# Returns the micromamba download URL for the current host platform/arch.
# Arguments:
#   OUT_VAR: Name of the variable to store the URL in
function(_mamba_get_url OUT_VAR)
    if (CMAKE_HOST_APPLE)
        if (CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "arm|arm64")
            set(_url "https://micro.mamba.pm/api/micromamba/osx-arm64/latest")
        else ()
            set(_url "https://micro.mamba.pm/api/micromamba/osx-64/latest")
        endif ()
    elseif (LINUX)
        if (CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "x86_64|x64")
            set(_url "https://micro.mamba.pm/api/micromamba/linux-64/latest")
        elseif (CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "aarch64|arm64")
            set(_url "https://micro.mamba.pm/api/micromamba/linux-aarch64/latest")
        else ()
            message(FATAL_ERROR "Unsupported Linux architecture: ${CMAKE_HOST_SYSTEM_PROCESSOR}")
        endif ()
    elseif (CMAKE_HOST_WIN32)
        set(_url "https://micro.mamba.pm/api/micromamba/win-64/latest")
    else ()
        message(FATAL_ERROR "Unsupported platform for micromamba download")
    endif ()
    set(${OUT_VAR} "${_url}" PARENT_SCOPE)
endfunction()


# Ensures the micromamba binary is available, first searching the system and
# then falling back to downloading it into CMAKE_BINARY_DIR.
#
# Search order:
#   1. System PATH (any micromamba already installed on the machine)
#   2. ${CMAKE_SOURCE_DIR}/micromamba  (vendored copy in-tree)
#   3. ${CMAKE_BINARY_DIR}/micromamba_bin  (previously downloaded by this script)
#   4. Download from micro.mamba.pm
#
# Arguments:
#   OUT_VAR: Name of the variable to store the path to the micromamba binary
function(ensure_micromamba OUT_VAR)
    set(_download_dir "${CMAKE_BINARY_DIR}/micromamba_bin")

    # Platform-specific sub-directory layout inside the downloaded tarball
    if (CMAKE_HOST_WIN32)
        set(_bin_subdir "Library/bin")
        set(_archive_entry "Library/bin/micromamba.exe")
    else ()
        set(_bin_subdir "bin")
        set(_archive_entry "bin/micromamba")
    endif ()

    # Steps 1-3: search system PATH, vendored copy, and previously downloaded binary.
    # NO_CACHE prevents CMake from storing a NOTFOUND result, which would cause the
    # post-download search to return the cached miss instead of the real binary.
    find_program(_micromamba_bin micromamba
        PATHS
            "${CMAKE_SOURCE_DIR}/micromamba"
            "${_download_dir}/${_bin_subdir}"
        NO_CACHE
    )

    if (_micromamba_bin)
        message(STATUS "Found micromamba: ${_micromamba_bin}")
    else ()
        # Step 4: download
        _mamba_get_url(_mamba_url)
        message(STATUS "micromamba not found on system – downloading from ${_mamba_url}")

        file(MAKE_DIRECTORY "${_download_dir}/${_bin_subdir}")

        set(_archive "${_download_dir}/micromamba.tar.bz2")
        file(DOWNLOAD "${_mamba_url}" "${_archive}"
            SHOW_PROGRESS
            STATUS _dl_status
        )
        list(GET _dl_status 0 _dl_error_code)
        if (_dl_error_code)
            list(GET _dl_status 1 _dl_error_msg)
            message(FATAL_ERROR "Failed to download micromamba: ${_dl_error_msg}")
        endif ()

        # Use 'xf' without an explicit compression flag: cmake -E tar auto-detects
        # the format from the file header, which is more portable (especially on Windows).
        execute_process(
            COMMAND ${CMAKE_COMMAND} -E tar xf "${_archive}" -- "${_archive_entry}"
            WORKING_DIRECTORY "${_download_dir}"
            RESULT_VARIABLE _extract_result
        )
        if (_extract_result)
            message(FATAL_ERROR "Failed to extract micromamba archive (exit code: ${_extract_result})")
        endif ()

        find_program(_micromamba_bin micromamba
            PATHS "${_download_dir}/${_bin_subdir}"
            NO_DEFAULT_PATH
            NO_CACHE
        )
        if (NOT _micromamba_bin)
            message(FATAL_ERROR "micromamba binary not found after download/extract – check the archive layout")
        endif ()

        message(STATUS "micromamba downloaded to: ${_micromamba_bin}")
    endif ()

    set(${OUT_VAR} "${_micromamba_bin}" PARENT_SCOPE)
endfunction()


# Checks whether a micromamba environment at ENV_PATH has been created.
# Arguments:
#   ENV_PATH: Filesystem path to the environment
#   OUT_VAR:  Name of the boolean variable to set (TRUE if the env exists)
function(_check_micromamba_environment ENV_PATH OUT_VAR)
    if (EXISTS "${ENV_PATH}/conda-meta/history")
        set(${OUT_VAR} TRUE PARENT_SCOPE)
    else ()
        set(${OUT_VAR} FALSE PARENT_SCOPE)
    endif ()
endfunction()


# Create (or update) a micromamba environment.
#
# Named arguments:
#   NAME           Name of the environment (default: "environment")
#   VERBOSE        Verbosity level: 0=quiet, 1=normal, 2=verbose, 3=very verbose
#   RUN_CMD        Variable to receive the command prefix for running inside the env
#   ENV_PATH       Path for the environment. If given, used directly as the
#                  environment location (input). If omitted, defaults to
#                  ${CMAKE_BINARY_DIR}/environments/${NAME}.
#   CHANNELS       List of conda channels (-c flags)
#   DEPENDENCIES   Package names to install
#   SPEC_FILE      Conda spec / environment YAML files (-f flags)
#
# Boolean flags:
#   NO_PREFIX_PATH  Do not prepend the env path to CMAKE_PREFIX_PATH
#   NO_DEPENDS      Do not register spec files as CMake configure dependencies
#
function(micromamba_environment)
    cmake_parse_arguments(PARSED
        "NO_PREFIX_PATH;NO_DEPENDS"
        "NAME;VERBOSE;RUN_CMD;ENV_PATH"
        "CHANNELS;DEPENDENCIES;SPEC_FILE"
        ${ARGN}
    )

    # ---- Resolve environment name & path --------------------------------
    if (DEFINED PARSED_NAME)
        set(_env_name "${PARSED_NAME}")
    else ()
        set(_env_name "environment")
    endif ()

    # ENV_PATH can be provided directly as an input to override the default location.
    # If not given, fall back to a path inside the build directory.
    if (DEFINED PARSED_ENV_PATH)
        set(_env_path "${PARSED_ENV_PATH}")
    else ()
        set(_env_path "${CMAKE_BINARY_DIR}/environments/${_env_name}")
    endif ()

    # ---- Verbosity flag -------------------------------------------------
    set(_verbose_flag "--quiet")
    if (DEFINED PARSED_VERBOSE)
        if (PARSED_VERBOSE EQUAL 1)
            set(_verbose_flag "")
        elseif (PARSED_VERBOSE EQUAL 2)
            set(_verbose_flag "--verbose")
        elseif (PARSED_VERBOSE EQUAL 3)
            set(_verbose_flag "-vv")
        endif ()
    endif ()

    # ---- Channel & spec-file argument lists -----------------------------
    set(_channel_args "")
    foreach (_ch ${PARSED_CHANNELS})
        list(APPEND _channel_args "-c" "${_ch}")
    endforeach ()

    set(_file_args "")
    foreach (_f ${PARSED_SPEC_FILE})
        list(APPEND _file_args "-f" "${_f}")
    endforeach ()

    # ---- Resolve micromamba binary --------------------------------------
    if (OVERWRITE_MICROMAMBA)
        set(_mamba_bin "${OVERWRITE_MICROMAMBA}")
    else ()
        ensure_micromamba(_mamba_bin)
    endif ()

    # Flags applied to every micromamba invocation:
    #   --no-rc          ignore all system/user .mambarc and condarc files
    #   --root-prefix    store micromamba's own state (pkgs cache, etc.) inside
    #                    the build tree instead of the default location, which may
    #                    be read-only (e.g. a Nix store path)
    set(_mamba_base_args
        --no-rc
        --root-prefix "${CMAKE_BINARY_DIR}/micromamba_root"
    )

    # ---- Create or update the environment -------------------------------
    _check_micromamba_environment("${_env_path}" _env_exists)

    if (_env_exists)
        message(STATUS "Updating micromamba environment '${_env_name}' at ${_env_path}")
        execute_process(
            COMMAND ${_mamba_bin} update
                ${_mamba_base_args}
                -p "${_env_path}"
                ${_channel_args}
                ${_file_args}
                ${PARSED_DEPENDENCIES}
                --yes ${_verbose_flag}
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
            RESULT_VARIABLE _mamba_result
        )
    else ()
        message(STATUS "Creating micromamba environment '${_env_name}' at ${_env_path} – this may take a while")
        execute_process(
            COMMAND ${_mamba_bin} create
                ${_mamba_base_args}
                -p "${_env_path}"
                ${_channel_args}
                ${_file_args}
                ${PARSED_DEPENDENCIES}
                --yes ${_verbose_flag}
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
            RESULT_VARIABLE _mamba_result
        )
    endif ()

    if (_mamba_result)
        message(FATAL_ERROR "micromamba exited with error code ${_mamba_result}")
    endif ()

    # ---- Register spec files as configure-time dependencies -------------
    if (NOT PARSED_NO_DEPENDS AND PARSED_SPEC_FILE)
        set_property(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
            APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${PARSED_SPEC_FILE}
        )
    endif ()

    # ---- Propagate prefix path ------------------------------------------
    if (NOT PARSED_NO_PREFIX_PATH)
        set(CMAKE_INSTALL_PREFIX "${_env_path}" CACHE INTERNAL "CMAKE_INSTALL_PREFIX")
        list(PREPEND CMAKE_PREFIX_PATH "${_env_path}")
        set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}" PARENT_SCOPE)
    endif ()

    # ---- Propagate output variables to the caller -----------------------
    if (DEFINED PARSED_RUN_CMD)
        set(${PARSED_RUN_CMD} "${_mamba_bin}" ${_mamba_base_args} -p "${_env_path}" run PARENT_SCOPE)
    endif ()

    set(MICROMAMBA_ENV_PATH "${_env_path}" PARENT_SCOPE)
endfunction()
