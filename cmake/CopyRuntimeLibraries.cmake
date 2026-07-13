include_guard(GLOBAL)

function(nosphera2_copy_runtime_libraries target)
    if(NOT TARGET "${target}")
        message(FATAL_ERROR
            "nosphera2_copy_runtime_libraries(): "
            "target '${target}' does not exist"
        )
    endif()

    if(NOT TARGET TBB::tbb)
        message(FATAL_ERROR
            "nosphera2_copy_runtime_libraries(): "
            "TBB::tbb does not exist"
        )
    endif()

    if(NOT DEFINED MICROMAMBA_ENV_PREFIX)
        message(FATAL_ERROR
            "MICROMAMBA_ENV_PREFIX is not defined"
        )
    endif()

    # Select the OpenMP runtime provided by the Micromamba environment.
    if(WIN32)
        set(_openmp_source
            "${MICROMAMBA_ENV_PREFIX}/Library/bin/libiomp5md.dll"
        )
        set(_tbb_destination_name
            "$<TARGET_FILE_NAME:TBB::tbb>"
        )
    elseif(APPLE)
        set(_openmp_source
            "${MICROMAMBA_ENV_PREFIX}/lib/libomp.dylib"
        )
        set(_tbb_destination_name
            "$<TARGET_SONAME_FILE_NAME:TBB::tbb>"
        )
        set(_runtime_rpath
            "@loader_path"
        )
    elseif(UNIX)
        set(_openmp_source
            "${MICROMAMBA_ENV_PREFIX}/lib/libiomp5.so"
        )
        set(_tbb_destination_name
            "$<TARGET_SONAME_FILE_NAME:TBB::tbb>"
        )
        set(_runtime_rpath
            "$ORIGIN"
        )
    else()
        message(FATAL_ERROR
            "Unsupported platform: ${CMAKE_SYSTEM_NAME}"
        )
    endif()

    if(NOT EXISTS "${_openmp_source}")
        message(FATAL_ERROR
            "OpenMP runtime does not exist: ${_openmp_source}"
        )
    endif()

    # The source may be a symlink. Copy the actual file while retaining the
    # public runtime filename, such as libiomp5.so.
    get_filename_component(
        _openmp_destination_name
        "${_openmp_source}"
        NAME
    )

    file(
        REAL_PATH
        "${_openmp_source}"
        _openmp_real_source
    )

    # Windows searches the executable directory automatically.
    if(DEFINED _runtime_rpath)
        set_target_properties(
            "${target}"
            PROPERTIES
                BUILD_WITH_INSTALL_RPATH TRUE
                INSTALL_RPATH "${_runtime_rpath}"
                INSTALL_RPATH_USE_LINK_PATH FALSE
        )
    endif()

    add_custom_command(
        TARGET "${target}"
        POST_BUILD

        # TBB: copy the real file under its runtime SONAME.
        #
        # Example:
        #   libtbb.so.12.16 -> libtbb.so.12
        COMMAND
            "${CMAKE_COMMAND}" -E copy_if_different
            "$<TARGET_FILE:TBB::tbb>"
            "$<TARGET_FILE_DIR:${target}>/${_tbb_destination_name}"

        # OpenMP: copy the real file under its public runtime name.
        COMMAND
            "${CMAKE_COMMAND}" -E copy_if_different
            "${_openmp_real_source}"
            "$<TARGET_FILE_DIR:${target}>/${_openmp_destination_name}"

        COMMENT
            "Copying TBB and OpenMP runtimes for ${target}"

        VERBATIM
    )
endfunction()