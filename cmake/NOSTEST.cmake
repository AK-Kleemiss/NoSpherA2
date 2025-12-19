function(add_nos_test TEST_NAME ARGS)
    cmake_parse_arguments(PARSE_ARGV 2 TEST "" "DIRECTORY;GOOD_LOG;ACTUAL_LOG" "")
    if (NOT TEST_ACTUAL_LOG)
        set(TEST_ACTUAL_LOG "NoSpherA2.log")
    endif ()
    if (NOT TEST_DIRECTORY)
        set(TEST_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME})
        set(TEST_BIN_DIR ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME})
    else ()
        set(TEST_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIRECTORY})
        set(TEST_BIN_DIR ${CMAKE_CURRENT_BINARY_DIR}/${TEST_DIRECTORY})
    endif ()

    file(COPY ${TEST_SRC_DIR} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    if (NOT TEST_GOOD_LOG)
        set(TEST_GOOD_LOG "${TEST_BIN_DIR}/${TEST_NAME}.good")
    endif ()
    add_test(
        NAME ${TEST_NAME}
        COMMAND ${Python3_EXECUTABLE}
                "${CMAKE_SOURCE_DIR}/tests/runTest.py"
                --exe $<TARGET_FILE:NoSpherA2>
                --actual "${TEST_ACTUAL_LOG}"
                --good "${TEST_GOOD_LOG}"
                --args="${ARGS}"
                --dir ${TEST_BIN_DIR}
        WORKING_DIRECTORY ${TEST_BIN_DIR}
    )
endfunction()