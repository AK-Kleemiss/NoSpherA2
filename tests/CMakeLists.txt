include(CTest)

if (NOT OCC_DATA_PATH)
    set(OCC_DATA_PATH "${CMAKE_SOURCE_DIR}/occ/share")
endif ()

find_package(Pytest REQUIRED HINTS "${CMAKE_SOURCE_DIR}"
    PATHS "$ENV{HOME}/.local/share" "$ENV{APPDATA}"
)

set(_ENVIRONMENT
    "NOS_EXE=$<TARGET_FILE:NoSpherA2>"
    "OCC_DATA_PATH=${OCC_DATA_PATH}"
)

if (RUN_FULL_TEST)
    list(APPEND _ENVIRONMENT "RUN_FULL_TEST=1")
endif ()

set(_EXCLUDE_PATTERN "")
if (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(_EXCLUDE_PATTERN "ri_fit")
endif ()

pytest_discover_tests(
    NoSpherA2Tests
    TEST_PATHS "${CMAKE_CURRENT_SOURCE_DIR}/run_test.py"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    TRIM_FROM_NAME "^test"
    STRIP_PARAM_BRACKETS
    ENVIRONMENT ${_ENVIRONMENT}
    DEPENDS NoSpherA2
    EXTRA_ARGS "--color=yes"
)

if (_EXCLUDE_PATTERN)
    set_tests_properties(${_EXCLUDE_PATTERN} PROPERTIES DISABLED TRUE)
endif ()

# Show full pytest output (with colour) when a test fails
file(WRITE "${CMAKE_BINARY_DIR}/CTestCustom.cmake"
    "set(CTEST_OUTPUT_ON_FAILURE TRUE)\n"
)
