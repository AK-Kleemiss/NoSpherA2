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
nosphera2_isa_check(NoSpherA2Dll DLL)

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
        BUILDING_DLL
)

set_target_properties(
    NoSpherA2Dll
    PROPERTIES
        OUTPUT_NAME "NoSpherA2"
        WINDOWS_EXPORT_ALL_SYMBOLS OFF
        #beside the NoSpherA2 launcher, which loads it from there, and in the packaged bin
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

#Linux: export the entry points only. Everything else default-visible is interposable in a .so, so GCC's
#LTO would not inline or internalise the core as it does in an executable; this also keeps the static
#BLAS/TBB/OpenMP copies out of the dynamic symbol table, as Olex2 holds its own
if(LINUX)
    target_link_options(NoSpherA2Dll PRIVATE "LINKER:--version-script=${CMAKE_CURRENT_SOURCE_DIR}/exports.map")
    set_property(TARGET NoSpherA2Dll APPEND PROPERTY LINK_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/exports.map")
endif()

nosphera2_copy_runtime_libraries(NoSpherA2Dll)

if(NOSPHERA2_USE_HIP)
    #see app/CMakeLists.txt: the kernels carry their own fat binaries, the C++ linker is all that is needed
    set_target_properties(NoSpherA2Dll PROPERTIES LINKER_LANGUAGE CXX)
endif()