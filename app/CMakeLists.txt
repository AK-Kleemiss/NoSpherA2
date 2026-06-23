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

target_compile_features(
    NoSpherA2
    PRIVATE
        cxx_std_20
)

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