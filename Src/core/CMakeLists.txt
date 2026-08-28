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

if(NOSPHERA2_USE_CUDA)
    list(APPEND NOSPHERA2_CORE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/sf_cuda.cu")
endif()

#Basis Set Generator:
set(
    NOSPHERA2_BASIS_SET_DIRECTORY
    "${PROJECT_SOURCE_DIR}/BasisSetGenerator/basis_sets"
)

set(
    NOSPHERA2_GENERATED_BASIS_SOURCE
    "${PROJECT_SOURCE_DIR}/Src/basis_data.cpp"
)

# The generated source depends on the CSVs it is generated FROM, not only on the
# converter that reads them.
#
# Without this the rule only ever fired when BasisSetConverter itself was rebuilt,
# so adding or editing a basis set left a stale basis_data.cpp in place - and it
# is .gitignore'd, so nothing in the repository showed a difference either. That
# is how this checkout ended up with 98 basis sets against 117 CSVs: STO-3G was
# missing, and every XCW test failed with "Basis set sto-3g not defined" until the
# file was deleted by hand.
#
# CONFIGURE_DEPENDS re-globs at build time, so a NEW csv triggers a reconfigure;
# naming the files in DEPENDS below covers an EDITED one. Both are needed - the
# glob alone would not notice a change to a file already in the list.
file(
    GLOB NOSPHERA2_BASIS_SET_FILES
    CONFIGURE_DEPENDS
    "${NOSPHERA2_BASIS_SET_DIRECTORY}/*.csv"
)

if(NOT NOSPHERA2_BASIS_SET_FILES)
    # Not fatal: the converter target would already be missing if the submodule
    # were absent. But silence here means a stale basis_data.cpp again, so say so.
    message(WARNING
        "No .csv basis sets found in ${NOSPHERA2_BASIS_SET_DIRECTORY} - "
        "basis_data.cpp will not regenerate when basis sets change. "
        "Is the BasisSetGenerator submodule checked out?"
    )
endif()

add_custom_command(
    OUTPUT
        "${NOSPHERA2_GENERATED_BASIS_SOURCE}"

    COMMAND
        "$<TARGET_FILE:BasisSetConverter>"
        "${NOSPHERA2_BASIS_SET_DIRECTORY}"
        "${PROJECT_SOURCE_DIR}/Src"

    DEPENDS
        BasisSetConverter
        ${NOSPHERA2_BASIS_SET_FILES}

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

nosphera2_enable_optimizations(NoSpherA2Core)

add_library(
    NoSpherA2::Core
    ALIAS
    NoSpherA2Core
)

target_compile_definitions(
    NoSpherA2Core
    PRIVATE
        NOSPHERA2_BASIS_DATA_SEPARATE
)

if(NOSPHERA2_USE_CUDA)
    target_compile_definitions(NoSpherA2Core PRIVATE NOSPHERA2_USE_CUDA)
    set_target_properties(NoSpherA2Core PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()

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
