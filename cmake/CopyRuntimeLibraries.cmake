include_guard(GLOBAL)

function(nosphera2_copy_runtime_libraries target)
    if(NOT TARGET "${target}")
        message(FATAL_ERROR
            "nosphera2_copy_runtime_libraries(): "
            "target '${target}' does not exist"
        )
    endif()

    # Android: TBB, OpenMP and OpenBLAS are linked statically; there is nothing to copy
    if(ANDROID)
        return()
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
    if(WIN32 AND NOSPHERA2_OPENBLAS)
        # ARM64: vcomp is OpenMP and comes with the VC redistributable; the DLL to ship is OpenBLAS
        set(_openmp_source
            "${NOSPHERA2_OPENBLAS_DLL}"
        )
        set(_tbb_destination_name
            "$<TARGET_FILE_NAME:TBB::tbb>"
        )
    elseif(WIN32)
        set(_openmp_source
            "${MICROMAMBA_ENV_PREFIX}/Library/bin/libiomp5md.dll"
        )
        set(_tbb_destination_name
            "$<TARGET_FILE_NAME:TBB::tbb>"
        )
    elseif(UNIX AND NOSPHERA2_OPENBLAS)
        # aarch64: libgomp is the system's; ship OpenBLAS and the gfortran runtime its LAPACK needs
        set(_openmp_source
            "${MICROMAMBA_ENV_PREFIX}/lib/libopenblas.so.0"
        )
        set(_extra_runtime_source
            "${MICROMAMBA_ENV_PREFIX}/lib/libgfortran.so.5"
        )
        set(_tbb_destination_name
            "$<TARGET_SONAME_FILE_NAME:TBB::tbb>"
        )
        set(_runtime_rpath
            "$ORIGIN"
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

    # ArmPL is linked statically and OpenMP is the system's (libgomp, vcomp): only TBB to ship, so
    # the OpenMP copy below repeats the TBB one, which copy_if_different skips
    if(NOSPHERA2_ARMPL)
        set(_openmp_source "$<TARGET_FILE:TBB::tbb>")
        set(_openmp_real_source "${_openmp_source}")
        set(_openmp_destination_name "${_tbb_destination_name}")
        unset(_extra_runtime_source)
    elseif(NOT EXISTS "${_openmp_source}")
        message(FATAL_ERROR
            "OpenMP runtime does not exist: ${_openmp_source}"
        )
    else()
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
    endif()

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

    if(DEFINED _extra_runtime_source)
        get_filename_component(_extra_destination_name "${_extra_runtime_source}" NAME)
        file(REAL_PATH "${_extra_runtime_source}" _extra_real_source)
        add_custom_command(
            TARGET "${target}"
            POST_BUILD
            COMMAND
                "${CMAKE_COMMAND}" -E copy_if_different
                "${_extra_real_source}"
                "$<TARGET_FILE_DIR:${target}>/${_extra_destination_name}"
            VERBATIM
        )
    endif()

    # No CUDA or ROCm runtime is copied. cuBLAS was the only one that ever needed to be, and
    # shipping it cost half a gigabyte for two GEMM calls; gemm_gpu.cuh replaced it.
endfunction()