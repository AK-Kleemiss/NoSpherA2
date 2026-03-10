include(CTest)

find_package(Python3 COMPONENTS Interpreter)

if (NOT OCC_DATA_PATH)
    set(OCC_DATA_PATH "${CMAKE_SOURCE_DIR}/occ/share")
endif ()

set(_NOS_TEST_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/run_test.py")
set(_NOS_TESTS
    alanine
    alanine_integrated_occ
    disorder_THPP
    fourier_transform
    fractal
    grown_water
    Hybrid_mode
    malbac_SF_ECP
    openBLAS
    properties
    reading_SALTED
    rubredoxin_cmtc
    SALTED
    sucrose_IAM
    sucrose_ptb
    sucrose_SF
    sucrose_twin
    wfn_reading
)

if (NOT (CMAKE_SYSTEM_NAME STREQUAL "Darwin"))
    list(APPEND _NOS_TESTS ri_fit)
endif ()

foreach(TEST_NAME IN LISTS _NOS_TESTS)
    add_test(
        NAME ${TEST_NAME}
        COMMAND ${Python3_EXECUTABLE} -m pytest ${_NOS_TEST_SCRIPT} -k ${TEST_NAME} -v
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
    set_tests_properties(${TEST_NAME} PROPERTIES
        ENVIRONMENT "NOS_EXE=$<TARGET_FILE:NoSpherA2>;OCC_DATA_PATH=${OCC_DATA_PATH}"
    )
endforeach()

if (RUN_FULL_TEST)
    foreach(TEST_NAME IN ITEMS fchk_conversion fourier_transform_full)
        add_test(
            NAME ${TEST_NAME}
            COMMAND ${Python3_EXECUTABLE} -m pytest ${_NOS_TEST_SCRIPT} -k ${TEST_NAME} -v
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        )
        set_tests_properties(${TEST_NAME} PROPERTIES
            ENVIRONMENT "NOS_EXE=$<TARGET_FILE:NoSpherA2>;OCC_DATA_PATH=${OCC_DATA_PATH};RUN_FULL_TEST=1"
        )
    endforeach()
endif ()
