if (WIN32)
    set(NoSpherA2_EXE ${CMAKE_BINARY_DIR}/Src/NoSpherA2.exe)
else ()
    set(NoSpherA2_EXE ${CMAKE_BINARY_DIR}/Src/NoSpherA2)
endif ()

function(add_nos_test test_name test_args)
    # Parse keyworded args: ARGS (single value), COMPARE_NAME (single value)
    cmake_parse_arguments(PARSE_ARGV 2 TEST "" "COMPARE_NAME;COMPARE_GOOD;DIRECTORY" "")

    message(STATUS "Running test for ${test_name}!")

    if(NOT TEST_COMPARE_NAME)
        set(TEST_COMPARE_NAME "${test_name}.good")
    endif()
    if (NOT TEST_DIRECTORY)
        set(TEST_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}")
        set(TEST_BIN_DIR ${CMAKE_CURRENT_BINARY_DIR}/${test_name})
    else ()
        set(TEST_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIRECTORY}")
        set(TEST_BIN_DIR ${CMAKE_CURRENT_BINARY_DIR}/${TEST_DIRECTORY})
    endif ()
    if (NOT TEST_COMPARE_GOOD)
        set(TEST_COMPARE_GOOD "NoSpherA2.log")
    endif ()
    file(COPY ${TEST_SRC_DIR} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

    add_test(
        NAME ${test_name}
        WORKING_DIRECTORY ${TEST_BIN_DIR}
        COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/runTest.py
                "${TEST_COMPARE_NAME}" "${TEST_COMPARE_GOOD}" "${NoSpherA2_EXE} ${test_args}"
    )
endfunction()
