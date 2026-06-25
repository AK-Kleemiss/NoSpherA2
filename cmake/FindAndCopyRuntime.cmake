cmake_minimum_required(VERSION 3.20)

if(NOT DEFINED RUNTIME_CANDIDATE)
    message(FATAL_ERROR "RUNTIME_CANDIDATE was not provided")
endif()

if(NOT DEFINED RUNTIME_DESTINATION)
    message(FATAL_ERROR "RUNTIME_DESTINATION was not provided")
endif()

if(DEFINED RUNTIME_NAMES)
    string(REPLACE "@@" ";" RUNTIME_NAMES "${RUNTIME_NAMES}")
else()
    set(RUNTIME_NAMES "")
endif()

if(DEFINED SEARCH_ROOTS)
    string(REPLACE "@@" ";" SEARCH_ROOTS "${SEARCH_ROOTS}")
else()
    set(SEARCH_ROOTS "")
endif()

set(runtime_file "")

# First try the exact path reported by the imported target.
if(EXISTS "${RUNTIME_CANDIDATE}")
    set(runtime_file "${RUNTIME_CANDIDATE}")
endif()

# If that path is invalid, search recursively.
if(NOT runtime_file)
    get_filename_component(candidate_name
        "${RUNTIME_CANDIDATE}"
        NAME
    )

    set(names_to_find
        "${candidate_name}"
        ${RUNTIME_NAMES}
    )

    list(REMOVE_DUPLICATES names_to_find)

    foreach(root IN LISTS SEARCH_ROOTS)
        if(NOT IS_DIRECTORY "${root}")
            continue()
        endif()

        foreach(name IN LISTS names_to_find)
            file(GLOB_RECURSE matches
                LIST_DIRECTORIES FALSE
                "${root}/${name}"
            )

            if(matches)
                list(GET matches 0 runtime_file)
                break()
            endif()
        endforeach()

        if(runtime_file)
            break()
        endif()
    endforeach()
endif()

if(NOT runtime_file)
    message(FATAL_ERROR
        "Could not find runtime library.\n"
        "Expected path: ${RUNTIME_CANDIDATE}\n"
        "Names searched: ${RUNTIME_NAMES}\n"
        "Roots searched: ${SEARCH_ROOTS}"
    )
endif()

file(MAKE_DIRECTORY "${RUNTIME_DESTINATION}")

execute_process(
    COMMAND "${CMAKE_COMMAND}" -E copy_if_different
        "${runtime_file}"
        "${RUNTIME_DESTINATION}"
    RESULT_VARIABLE copy_result
)

if(NOT copy_result EQUAL 0)
    message(FATAL_ERROR
        "Failed to copy '${runtime_file}' to "
        "'${RUNTIME_DESTINATION}'"
    )
endif()

message(STATUS
    "Copied runtime library '${runtime_file}' "
    "to '${RUNTIME_DESTINATION}'"
)