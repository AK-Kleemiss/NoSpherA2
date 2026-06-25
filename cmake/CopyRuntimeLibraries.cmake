set(_NOSPHERA2_RUNTIME_COPY_SCRIPT
    "${CMAKE_CURRENT_LIST_DIR}/FindAndCopyRuntime.cmake"
)

function(nosphera2_copy_tbb_runtime target)
    if(NOT TARGET "${target}")
        message(FATAL_ERROR "Unknown target: ${target}")
    endif()

    set(search_roots
        "${CMAKE_BINARY_DIR}"
        "${CMAKE_SOURCE_DIR}"
    )

    if(DEFINED MICROMAMBA_ENV_PREFIX)
        list(APPEND search_roots "${MICROMAMBA_ENV_PREFIX}")
    endif()

    list(JOIN search_roots "@@" search_roots_argument)

    if(WIN32)
        set(runtime_names
            "tbb12_debug.dll"
            "tbb12.dll"
            "tbb.dll"
        )
    elseif(APPLE)
        set(runtime_names
            "libtbb.12.dylib"
            "libtbb.dylib"
            "libtbb.*.dylib"
        )
    else()
        set(runtime_names
            "libtbb.so.12"
            "libtbb.so"
            "libtbb.so.*"
        )
    endif()

    list(JOIN runtime_names "@@" runtime_names_argument)

    add_custom_command(
        TARGET "${target}"
        POST_BUILD
        COMMAND "${CMAKE_COMMAND}"
            "-DRUNTIME_CANDIDATE=$<TARGET_FILE:TBB::tbb>"
            "-DRUNTIME_DESTINATION=$<TARGET_FILE_DIR:${target}>"
            "-DRUNTIME_NAMES=${runtime_names_argument}"
            "-DSEARCH_ROOTS=${search_roots_argument}"
            -P "${_NOSPHERA2_RUNTIME_COPY_SCRIPT}"
        COMMENT "Copying TBB runtime for ${target}"
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
    list(JOIN _openmp_runtime_names "@@" _openmp_runtime_names_arg)
    list(JOIN _search_roots "@@" _search_roots_arg)

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
