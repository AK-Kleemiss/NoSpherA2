set(_NOSPHERA2_RUNTIME_COPY_SCRIPT
    "${CMAKE_CURRENT_LIST_DIR}/FallbackFindAndCopyRuntime.cmake"
)

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
        set(_tbb_runtime_names "tbb12.dll")
    elseif(APPLE)
        # Different TBB packages use either the ABI soname or the full version.
        set(_tbb_runtime_names
            "libtbb.12.dylib;libtbb.dylib;libtbb.*.dylib"
        )
    elseif(UNIX)
        set(_tbb_runtime_names
            "libtbb.so.12;libtbb.so;libtbb.so.*"
        )
    else()
        message(FATAL_ERROR "Unsupported platform")
    endif()

    set(_search_roots
        "${CMAKE_BINARY_DIR}"
        "${CMAKE_SOURCE_DIR}"
        "${MICROMAMBA_ENV_PREFIX}"
    )

    foreach(_prefix IN LISTS CMAKE_PREFIX_PATH)
        if(IS_ABSOLUTE "${_prefix}")
            list(APPEND _search_roots "${_prefix}")
        endif()
    endforeach()

    list(REMOVE_DUPLICATES _search_roots)
    list(JOIN _tbb_runtime_names "|" _tbb_runtime_names_arg)
    list(JOIN _search_roots "|" _search_roots_arg)

    add_custom_command(
        TARGET "${target}"
        POST_BUILD
        COMMAND
            "${CMAKE_COMMAND}"
            "-DRUNTIME_CANDIDATE=$<TARGET_FILE:${_tbb_target}>"
            "-DRUNTIME_DESTINATION=$<TARGET_FILE_DIR:${target}>"
            "-DRUNTIME_NAMES=${_tbb_runtime_names_arg}"
            "-DSEARCH_ROOTS=${_search_roots_arg}"
            -P "${_NOSPHERA2_RUNTIME_COPY_SCRIPT}"
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
        set(_openmp_runtime_names "libiomp5md.dll")
    elseif(APPLE)
        set(_openmp_runtime
            "${MICROMAMBA_ENV_PREFIX}/lib/libomp.dylib"
        )
        set(_openmp_runtime_names "libomp.dylib;libomp.*.dylib")
    else()
        set(_openmp_runtime
            "${MICROMAMBA_ENV_PREFIX}/lib/libiomp5.so"
        )
        set(_openmp_runtime_names "libiomp5.so;libiomp5.so.*;libomp.so;libomp.so.*")
    endif()

    set(_search_roots
        "${MICROMAMBA_ENV_PREFIX}"
        "${CMAKE_BINARY_DIR}"
        "${CMAKE_SOURCE_DIR}"
    )

    foreach(_prefix IN LISTS CMAKE_PREFIX_PATH)
        if(IS_ABSOLUTE "${_prefix}")
            list(APPEND _search_roots "${_prefix}")
        endif()
    endforeach()

    list(REMOVE_DUPLICATES _search_roots)
    list(JOIN _openmp_runtime_names "|" _openmp_runtime_names_arg)
    list(JOIN _search_roots "|" _search_roots_arg)

    add_custom_command(
        TARGET "${target}"
        POST_BUILD
        COMMAND
            "${CMAKE_COMMAND}"
            "-DRUNTIME_CANDIDATE=${_openmp_runtime}"
            "-DRUNTIME_DESTINATION=$<TARGET_FILE_DIR:${target}>"
            "-DRUNTIME_NAMES=${_openmp_runtime_names_arg}"
            "-DSEARCH_ROOTS=${_search_roots_arg}"
            -P "${_NOSPHERA2_RUNTIME_COPY_SCRIPT}"
        COMMENT
            "Copying OpenMP runtime for ${target}"
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
