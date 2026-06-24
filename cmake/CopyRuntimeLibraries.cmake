

function(nosphera2_copy_tbb_runtime target)
    if(TARGET TBB::tbb)
        set(_tbb_target TBB::tbb)
    elseif(TARGET tbb)
        set(_tbb_target tbb)
    else()
        message(FATAL_ERROR
            "No TBB target found. Expected TBB::tbb or tbb."
        )
    endif()

    if(WIN32)
        set(_tbb_runtime_name "tbb12.dll")
    elseif(APPLE)
        set(_tbb_runtime_name "libtbb.12.dylib")
    elseif(UNIX)
        set(_tbb_runtime_name "libtbb.so.12")
    else()
        message(FATAL_ERROR "Unsupported platform")
    endif()

    add_custom_command(
        TARGET "${target}"
        POST_BUILD
        COMMAND
            "${CMAKE_COMMAND}" -E copy_if_different
            "$<TARGET_FILE:${_tbb_target}>"
            "$<TARGET_FILE_DIR:${target}>/${_tbb_runtime_name}"
        COMMENT
            "Copying TBB runtime for ${target}"
        VERBATIM
    )
endfunction()

function(nosphera2_copy_openmp_runtime target)
    if(WIN32)
        set(_openmp_runtime
            "${MICROMAMBA_ENV_PREFIX}/Library/bin/libiomp5md.dll"
        )
    elseif(APPLE)
        set(_openmp_runtime
            "${MICROMAMBA_ENV_PREFIX}/lib/libomp.dylib"
        )
    else()
        set(_openmp_runtime
            "${MICROMAMBA_ENV_PREFIX}/lib/libiomp5.so"
        )
    endif()

    if(NOT EXISTS "${_openmp_runtime}")
        message(FATAL_ERROR
            "OpenMP runtime not found:\n"
            "  ${_openmp_runtime}"
        )
    endif()

    add_custom_command(
        TARGET "${target}"
        POST_BUILD
        COMMAND
            "${CMAKE_COMMAND}" -E copy_if_different
            "${_openmp_runtime}"
            "$<TARGET_FILE_DIR:${target}>"
        VERBATIM
    )
endfunction()

function(nosphera2_copy_runtime_libraries target)
    if(NOT TARGET "${target}")
        message(FATAL_ERROR
            "nosphera2_copy_runtime_libraries(): "
            "target '${target}' does not exist"
        )
    endif()

    nosphera2_copy_tbb_runtime("${target}")
    nosphera2_copy_openmp_runtime("${target}")
endfunction()