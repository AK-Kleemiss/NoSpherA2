if(NOT DEFINED RUNTIME_CANDIDATE OR RUNTIME_CANDIDATE STREQUAL "")
    message(FATAL_ERROR "RUNTIME_CANDIDATE was not provided")
endif()

if(NOT DEFINED RUNTIME_DESTINATION OR RUNTIME_DESTINATION STREQUAL "")
    message(FATAL_ERROR "RUNTIME_DESTINATION was not provided")
endif()

set(_runtime "${RUNTIME_CANDIDATE}")

if(NOT EXISTS "${_runtime}")
    get_filename_component(_candidate_name "${RUNTIME_CANDIDATE}" NAME)

    set(_names "${_candidate_name}")
    if(DEFINED RUNTIME_NAMES AND NOT RUNTIME_NAMES STREQUAL "")
        string(REPLACE "|" ";" RUNTIME_NAMES "${RUNTIME_NAMES}")
        list(APPEND _names ${RUNTIME_NAMES})
    endif()
    list(REMOVE_DUPLICATES _names)

    set(_roots)
    if(DEFINED SEARCH_ROOTS AND NOT SEARCH_ROOTS STREQUAL "")
        string(REPLACE "|" ";" SEARCH_ROOTS "${SEARCH_ROOTS}")
        list(APPEND _roots ${SEARCH_ROOTS})
    endif()
    list(REMOVE_DUPLICATES _roots)

    message(STATUS
        "Runtime library not found at configured path:\n"
        "  ${RUNTIME_CANDIDATE}\n"
        "Searching under:\n  ${_roots}"
    )

    set(_matches)
    foreach(_root IN LISTS _roots)
        if(NOT IS_DIRECTORY "${_root}")
            continue()
        endif()

        foreach(_name IN LISTS _names)
            file(GLOB_RECURSE _found
                LIST_DIRECTORIES FALSE
                "${_root}/${_name}"
            )
            list(APPEND _matches ${_found})
        endforeach()
    endforeach()

    list(REMOVE_DUPLICATES _matches)

    if(NOT _matches)
        message(FATAL_ERROR
            "Could not find the runtime library.\n"
            "Configured path:\n  ${RUNTIME_CANDIDATE}\n"
            "Accepted names:\n  ${_names}\n"
            "Search roots:\n  ${_roots}"
        )
    endif()

    # Prefer a real file over a dangling symlink, then use the first match.
    foreach(_match IN LISTS _matches)
        if(EXISTS "${_match}")
            set(_runtime "${_match}")
            break()
        endif()
    endforeach()

    message(STATUS "Found runtime library: ${_runtime}")
endif()

file(MAKE_DIRECTORY "${RUNTIME_DESTINATION}")

execute_process(
    COMMAND "${CMAKE_COMMAND}" -E copy_if_different
        "${_runtime}"
        "${RUNTIME_DESTINATION}"
    RESULT_VARIABLE _copy_result
    ERROR_VARIABLE _copy_error
)

if(NOT _copy_result EQUAL 0)
    message(FATAL_ERROR
        "Failed to copy runtime library:\n"
        "  from: ${_runtime}\n"
        "  to:   ${RUNTIME_DESTINATION}\n"
        "${_copy_error}"
    )
endif()
