file(
    GLOB NOSPHERA2_APP_SOURCES
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx"
)

add_executable(
    NoSpherA2
    ${NOSPHERA2_APP_SOURCES}
)

if(WIN32)
    target_sources(
        NoSpherA2
        PRIVATE
            "${PROJECT_SOURCE_DIR}/Windows/NoSpherA2/Resource.rc"
    )

    target_include_directories(
        NoSpherA2
        PRIVATE
            "${PROJECT_SOURCE_DIR}/Windows/NoSpherA2"
    )
endif()

nosphera2_enable_optimizations(NoSpherA2)

target_compile_features(
    NoSpherA2
    PRIVATE
        cxx_std_20
)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND NOT APPLE AND NOT WIN32)
    target_link_options(
        NoSpherA2
        PRIVATE
            -static-libstdc++
            -static-libgcc
    )
endif()

target_link_libraries(
    NoSpherA2
    PRIVATE
        NoSpherA2::Core
)

set_target_properties(
    NoSpherA2
    PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY
            "${CMAKE_BINARY_DIR}/bin"
)


nosphera2_copy_runtime_libraries(NoSpherA2)
