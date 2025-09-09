include(ExternalProject)
#CPMAddPackage(
#    NAME Ccache.cmake
#    GITHUB_REPOSITORY TheLartians/Ccache.cmake
#    VERSION 1.2.5
#)
#
CPMAddPackage(
    NAME mdspan
    GITHUB_REPOSITORY "kokkos/mdspan"
    GIT_TAG "stable"
)

CPMAddPackage(
    NAME metatensor
    GITHUB_REPOSITORY "metatensor/metatensor"
    GIT_TAG "metatensor-core-v0.1.17"
    OPTIONS
        "CMAKE_BUILD_TYPE Release"
        "BUILD_SHARED_LIBS OFF"
        "METATENSOR_INSTALL_BOTH_STATIC_SHARED OFF"
)

set(metatensor_STATIC_LIB ${metatensor_BINARY_DIR}/metatensor-core/target/release/libmetatensor.a)
CPMAddPackage(
    NAME featomic
    URL "https://github.com/metatensor/featomic/releases/download/featomic-v0.6.1/featomic-cxx-0.6.1.tar.gz"
#    DOWNLOAD_ONLY YES
    OPTIONS
        "CMAKE_BUILD_TYPE Release"
        "FEATOMIC_FETCH_METATENSOR OFF"
        "BUILD_SHARED_LIBS OFF"
        "FEATOMIC_INSTALL_BOTH_STATIC_SHARED OFF"
)
if (featomic_ADDED)
    add_library(featomic_static INTERFACE IMPORTED GLOBAL)
    target_link_libraries(featomic_static INTERFACE
        ${metatensor_STATIC_LIB}
        "${featomic_BINARY_DIR}/target/release/libfeatomic.a"
    )
    target_include_directories(
        featomic_static INTERFACE
           "${featomic_SOURCE_DIR}/include"
    )
endif ()

CPMAddPackage(
    NAME cmakeLLVM
    URL "https://github.com/llvm/llvm-project/releases/download/llvmorg-21.1.0/cmake-21.1.0.src.tar.xz"
    DOWNLOAD_ONLY YES
)

CPMAddPackage(
    NAME openmpLLVM
    URL "https://github.com/llvm/llvm-project/releases/download/llvmorg-21.1.0/openmp-21.1.0.src.tar.xz"
    DOWNLOAD_ONLY YES
)

if (cmakeLLVM_ADDED AND openmpLLVM_ADDED)
    set(openmp_SOURCE_DIR "${CMAKE_BINARY_DIR}/_deps/openmpLLVM-src")
    file(MAKE_DIRECTORY ${openmp_SOURCE_DIR})
    file(COPY "${cmakeLLVM_SOURCE_DIR}/" DESTINATION "${openmp_SOURCE_DIR}/cmake")
    file(COPY "${openmpLLVM_SOURCE_DIR}/" DESTINATION "${openmp_SOURCE_DIR}/openmp")
    set(openmp_INSTALL_DIR ${CMAKE_BINARY_DIR}/_deps/openmpLLVM-install)
    file(MAKE_DIRECTORY ${openmp_INSTALL_DIR})
    file(MAKE_DIRECTORY ${openmp_INSTALL_DIR}/include)
    message(WARNING ${CMAKE_BUILD_PARALLEL_LEVEL})
    ExternalProject_Add(
        openmpLLVM
        SOURCE_DIR "${openmp_SOURCE_DIR}/openmp"
        BINARY_DIR "${CMAKE_BINARY_DIR}/_deps/openmpLLVM-build"
        CMAKE_ARGS
            -DLIBOMP_ENABLE_SHARED=OFF
            -DLIBOMP_FORTRAN_MODULES=ON
        INSTALL_COMMAND ${CMAKE_COMMAND} --install . --prefix ${openmp_INSTALL_DIR}
        BUILD_BYPRODUCTS ${openmp_INSTALL_DIR}/lib/libomp.a
    )
    add_library(openmp INTERFACE IMPORTED GLOBAL)
    add_dependencies(openmp openmpLLVM)
    target_link_libraries(openmp INTERFACE ${openmp_INSTALL_DIR}/lib/libomp.a)
    target_include_directories(
        openmp INTERFACE
           ${openmp_INSTALL_DIR}/include
    )
endif ()


CPMAddPackage(
    NAME openblas
    URL "https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.30/OpenBLAS-0.3.30.tar.gz"
    DOWNLOAD_ONLY YES
)

if (openblas_ADDED)
    ExternalProject_Add(
        OpenBLAS
        SOURCE_DIR ${openblas_SOURCE_DIR}
        DEPENDS openmp
        BUILD_IN_SOURCE ON
        BUILD_COMMAND make FC=gfortran USE_OPENMP=1 -j
        INSTALL_COMMAND make install PREFIX=${CMAKE_BINARY_DIR}/3rdparty/openblas

        BUILD_BYPRODUCTS ${CMAKE_BINARY_DIR}/3rdparty/openblas/lib/libopenblas.a
        LOG_BUILD ON
    )
    # Create imported target
    add_library(openblas STATIC IMPORTED GLOBAL)
    add_dependencies(openblas OpenBLAS)
    set_target_properties(openblas PROPERTIES
        IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/3rdparty/openblas/lib/libopenblas.a
    )
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/3rdparty/openblas/include)
    set(OPENBLAS_INCLUDE_DIR ${CMAKE_BINARY_DIR}/3rdparty/openblas/include PARENT_SCOPE)
    message(WARNING ${OPENBLAS_INCLUDE_DIR})
    target_include_directories(openblas INTERFACE ${CMAKE_BINARY_DIR}/3rdparty/openblas/include)

    set(DESTINATION_INCLUDE_DIR "${CMAKE_BINARY_DIR}/3rdparty/openblas/include")
    set(openblas_SOURCE_DIR ${CMAKE_BINARY_DIR}/_deps/openblas-src)

    message(WARNING "SOURCE DIR: ${openblas_SOURCE_DIR}")
    file(MAKE_DIRECTORY "${DESTINATION_INCLUDE_DIR}")
endif ()

