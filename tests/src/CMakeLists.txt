find_package(
    GTest
    CONFIG
    REQUIRED
)

file(
    GLOB_RECURSE NOSPHERA2_TEST_SOURCES
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx"
)

if(NOT NOSPHERA2_TEST_SOURCES)
    message(STATUS
        "No test source files found; NoSpherA2_Tests will not be built"
    )
    return()
endif()

add_executable(
    NoSpherA2_Tests
    ${NOSPHERA2_TEST_SOURCES}
)

target_compile_features(
    NoSpherA2_Tests
    PRIVATE
        cxx_std_20
)

target_link_libraries(
    NoSpherA2_Tests
    PRIVATE
        NoSpherA2::Core
        GTest::gtest
)

target_precompile_headers(
    NoSpherA2_Tests
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/pch.h"
)

set(
    NOSPHERA2_OCC_DATA_PATH
    "${PROJECT_SOURCE_DIR}/occ/share"
    CACHE PATH
    "Directory containing OCC runtime data"
)


include(GoogleTest)

gtest_discover_tests(
    NoSpherA2_Tests
    DISCOVERY_MODE PRE_TEST
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    PROPERTIES
        ENVIRONMENT
            "OCC_DATA_PATH=${NOSPHERA2_OCC_DATA_PATH}"
)

set_property(
    DIRECTORY
    APPEND
    PROPERTY TEST_INCLUDE_FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/SetIntegrationTestLocks.cmake"
)

set_target_properties(
    NoSpherA2_Tests
    PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY
            "${CMAKE_BINARY_DIR}/bin"
)

nosphera2_copy_runtime_libraries(NoSpherA2_Tests)

if(WIN32)
    if(NOT TARGET GTest::gtest)
        message(FATAL_ERROR
            "GTest::gtest target is unavailable; cannot copy the GTest runtime DLL"
        )
    endif()

    add_custom_command(
        TARGET NoSpherA2_Tests
        POST_BUILD
        COMMAND
            "${CMAKE_COMMAND}" -E copy_if_different
            "$<TARGET_FILE:GTest::gtest>"
            "$<TARGET_FILE_DIR:NoSpherA2_Tests>"
        COMMENT
            "Copying GTest runtime DLL for NoSpherA2_Tests"
        VERBATIM
    )
endif()
