include(CTest)
find_package(Python COMPONENTS Interpreter REQUIRED)

if (WIN32)
    set(NoSpherA2_EXE ${CMAKE_BINARY_DIR}/Src/NoSpherA2.exe)
else ()
    set(NoSpherA2_EXE ${CMAKE_BINARY_DIR}/Src/NoSpherA2)
endif ()
FILE(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runTest.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

function(add_nos_test test_name)
    # Parse keyworded args: ARGS (single value), COMPARE_NAME (single value)
    cmake_parse_arguments(TEST "" "ARGS;COMPARE_NAME" "" ${ARGN})

    message(STATUS "Running test for ${test_name}!")

    if(NOT TEST_ARGS)
        set(TEST_ARGS "-cif ${test_name}.cif -dmin 0.5 -wfn ${test_name}.owf.fchk -acc 1 -no-date")
    endif()

    if(NOT TEST_COMPARE_NAME)
        set(TEST_COMPARE_NAME "${test_name}.good")
    endif()

    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${test_name} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

    add_test(
        NAME ${test_name}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_name}
        COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/runTest.py
                "${TEST_COMPARE_NAME}" "${NoSpherA2_EXE} ${TEST_ARGS}"
    )
endfunction()


add_nos_test(alanine_occ ARGS "-cif alanine.cif -dmin 0.5 -wfn alanine.owf.fchk -acc 1 -no-date")
add_nos_test(disorder ARGS "-cif thpp.cif -hkl thpp.hkl -mtc olex2/Wfn_job/Part_1/thpp.wfx 0.1 olex2/Wfn_job/Part_2/thpp.wfx 0.2 -mtc_mult 1 1 -mtc_charge 0 0 -acc 0 -no-date" COMPARE_NAME "disorder_THPP.good")

message(WARNING ${CMAKE_CURRENT_BINARY_DIR})
