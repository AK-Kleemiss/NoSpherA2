file(
    GLOB_RECURSE NOSPHERA2_CORE_SOURCES
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx"
)

list(
    FILTER NOSPHERA2_CORE_SOURCES
    EXCLUDE
    REGEX "[/\\\\](vec_nos_math|mat_nos_math)\\.cpp$"
)

file(
    GLOB_RECURSE NOSPHERA2_CORE_HEADERS
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp"
)

#Basis Set Generator:

set(
    NOSPHERA2_BASIS_SET_DIRECTORY
    "${PROJECT_SOURCE_DIR}/BasisSetGenerator/basis_sets"
)

set(
    NOSPHERA2_GENERATED_DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/generated"
)

set(
    NOSPHERA2_GENERATED_BASIS_SOURCE
    "${PROJECT_SOURCE_DIR}/Src/basis_data.cpp"
)

add_custom_command(
    OUTPUT
        "${NOSPHERA2_GENERATED_BASIS_SOURCE}"

    COMMAND
        "${CMAKE_COMMAND}" -E make_directory
        "${NOSPHERA2_GENERATED_DIRECTORY}"

    COMMAND
        "$<TARGET_FILE:BasisSetConverter>"
        "${NOSPHERA2_BASIS_SET_DIRECTORY}"
        "${PROJECT_SOURCE_DIR}/Src"

    DEPENDS
        BasisSetConverter

    WORKING_DIRECTORY
        "${PROJECT_SOURCE_DIR}/BasisSetGenerator"

    COMMENT
        "Generating basis_data.cpp"

    VERBATIM
)

#Generator end....

add_library(
    NoSpherA2Core
    STATIC
        ${NOSPHERA2_CORE_SOURCES}
        ${NOSPHERA2_CORE_HEADERS}
        "${NOSPHERA2_GENERATED_BASIS_SOURCE}"
)

add_library(
    NoSpherA2::Core
    ALIAS
    NoSpherA2Core
)

target_compile_features(
    NoSpherA2Core
    PUBLIC
        cxx_std_20
)

target_include_directories(
    NoSpherA2Core
    PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
        "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/Src>"
)

target_precompile_headers(
    NoSpherA2Core
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/pch.h"
)

if(APPLE)
    find_library(
        NOSPHERA2_ACCELERATE_FRAMEWORK
        Accelerate
        REQUIRED
    )
    target_link_libraries(
        NoSpherA2Core
        PRIVATE
            ${NOSPHERA2_ACCELERATE_FRAMEWORK}
        PUBLIC
            OpenMP::OpenMP_CXX
    )

    target_compile_definitions(
        NoSpherA2Core
        PRIVATE
            ACCELERATE_NEW_LAPACK
    )
else()
    target_link_libraries(
        NoSpherA2Core
        PUBLIC
            MKL::MKL
    )
endif()

target_link_libraries(
    NoSpherA2Core
    PUBLIC
        libcint::libcint
        occ_main
        featomic::static
)