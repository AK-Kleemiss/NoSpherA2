file(
    GLOB_RECURSE NOSPHERA2_DLL_SOURCES
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx"
)

file(
    GLOB_RECURSE NOSPHERA2_DLL_HEADERS
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx"
)


add_library(
    NoSpherA2Dll
    SHARED
        ${NOSPHERA2_DLL_SOURCES}
        ${NOSPHERA2_DLL_HEADERS}
)

nosphera2_enable_optimizations(NoSpherA2Dll)

add_library(
    NoSpherA2::Dll
    ALIAS
    NoSpherA2Dll
)

target_compile_features(
    NoSpherA2Dll
    PUBLIC
        cxx_std_20
)

target_include_directories(
    NoSpherA2Dll
    PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
        "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/Src>"
)

target_link_libraries(
    NoSpherA2Dll
    PRIVATE
        NoSpherA2::Core
)

target_precompile_headers(
    NoSpherA2Dll
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/pch.h"
)

target_compile_definitions(
    NoSpherA2Dll
    PRIVATE
        NOSPHERA2_DLL_BUILD
)

set_target_properties(
    NoSpherA2Dll
    PROPERTIES
        OUTPUT_NAME "NoSpherA2"
        WINDOWS_EXPORT_ALL_SYMBOLS OFF
)

nosphera2_copy_runtime_libraries(NoSpherA2Dll)