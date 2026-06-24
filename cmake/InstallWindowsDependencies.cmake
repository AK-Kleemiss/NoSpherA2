if(NOT WIN32)
    return()
endif()

add_custom_target(NoSpherA2WindowsDependencies ALL
    DEPENDS
        cint
        tbb
        occ
        BasisSetConverter
        cargo-build-featomic
)

install(
    TARGETS
        cint
        tbb
        spdlog
        gau2grid_static
        fmt
        dftd4
        xc
        scn
        _subprocess
        libecpint_static
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)

install(
    TARGETS
    occ
    occ_cg
    occ_core
    occ_crystal
    occ_descriptors
    occ_dft
    occ_disp
    occ_dma
    occ_driver
    occ_elastic_fit
    occ_geometry
    occ_gto
    occ_interaction
    occ_ints
    occ_io
    occ_isosurface
    occ_main
    occ_numint
    occ_opt
    occ_qm
    occ_sht
    occ_slater
    occ_solvent
    occ_xdm
    occ_xtb
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)


set(NOSPHERA2_EIGEN_SOURCE_DIR "")

install(
    DIRECTORY "${CMAKE_BINARY_DIR}/_deps/eigen3-src/Eigen/"
    DESTINATION include/Eigen
)
install(
    DIRECTORY "${CMAKE_BINARY_DIR}/_deps/eigen3-src/unsupported/Eigen/"
    DESTINATION include/unsupported/Eigen
)

install(
    DIRECTORY "${CMAKE_BINARY_DIR}/_deps/cli11-src/include/"
    DESTINATION include
    FILES_MATCHING
        PATTERN "*.h"
        PATTERN "*.hpp"
)

install(
    DIRECTORY "${CMAKE_BINARY_DIR}/_deps/fmt-src/include/"
    DESTINATION include
    FILES_MATCHING
        PATTERN "*.h"
        PATTERN "*.hpp"
)

install(
    DIRECTORY "${CMAKE_BINARY_DIR}/_deps/unordered_dense-src/include/"
    DESTINATION include
    FILES_MATCHING
        PATTERN "*.h"
        PATTERN "*.hpp"
)

install(
    DIRECTORY "${CMAKE_BINARY_DIR}/_deps/gemmi-src/include/"
    DESTINATION include
    FILES_MATCHING
        PATTERN "*.h"
        PATTERN "*.hpp"
)

install(
    DIRECTORY "${CMAKE_BINARY_DIR}/_deps/spdlog-src/include/"
    DESTINATION include
    FILES_MATCHING
        PATTERN "*.h"
        PATTERN "*.hpp"
)

install(
    DIRECTORY "${CMAKE_BINARY_DIR}/_deps/onetbb-src/include/"
    DESTINATION include
    FILES_MATCHING
        PATTERN "*.h"
        PATTERN "*.hpp"
)

install(
    DIRECTORY "${occ_SOURCE_DIR}/src/3rdparty/libecpint/"
    DESTINATION include
    FILES_MATCHING
        PATTERN "*.h"
        PATTERN "*.hpp"
)
install(
    DIRECTORY "${occ_SOURCE_DIR}/src/3rdparty/gau2grid/"
    DESTINATION include
    FILES_MATCHING
        PATTERN "*.h"
        PATTERN "*.hpp"
)

# libcint public headers
install(
    DIRECTORY "${cint_SOURCE_DIR}/include/"
    DESTINATION include
    FILES_MATCHING
        PATTERN "*.h"
        PATTERN "*.hpp"
)
install(
    DIRECTORY "${CMAKE_BINARY_DIR}/_deps/libcint/include/"
    DESTINATION include
    FILES_MATCHING
        PATTERN "*.h"
        PATTERN "*.hpp"
)

file(GLOB_RECURSE FEATOMIC_LIBS "${CMAKE_BINARY_DIR}/_deps/featomic/target/*/*.lib")
file(GLOB_RECURSE METATENSOR_LIBS "${CMAKE_BINARY_DIR}/_deps/metatensor-build/target/*/*.lib")
install(
    FILES
        ${FEATOMIC_LIBS}
        ${METATENSOR_LIBS}
    DESTINATION lib
)

# featomic public C API headers
install(
    DIRECTORY "${featomic_SOURCE_DIR}/include/"
    DESTINATION include
    FILES_MATCHING
        PATTERN "*.h"
        PATTERN "*.hpp"
)

install(
    DIRECTORY "${CMAKE_BINARY_DIR}/_deps/metatensor-src/include/"
    DESTINATION include
    FILES_MATCHING
        PATTERN "*.h"
        PATTERN "*.hpp"
)

install(
    DIRECTORY "${CMAKE_BINARY_DIR}/_deps/metatensor-build/include/"
    DESTINATION include
    FILES_MATCHING
        PATTERN "*.h"
        PATTERN "*.hpp"
)

install(
    DIRECTORY "${CMAKE_SOURCE_DIR}/occ/include/"
    DESTINATION include
)

install(
    FILES
        "${MICROMAMBA_ENV_PREFIX}/Library/bin/libiomp5md.dll"
    DESTINATION bin
)

install(
    DIRECTORY "${MICROMAMBA_ENV_PREFIX}/Library/include/"
    DESTINATION include
    FILES_MATCHING
        PATTERN "mkl*.h"
        PATTERN "mkl*.hpp"
)


install(
    FILES
        "${MICROMAMBA_ENV_PREFIX}/Library/lib/libiomp5md.lib"
        "${MICROMAMBA_ENV_PREFIX}/Library/lib/mkl_intel_lp64.lib"
        "${MICROMAMBA_ENV_PREFIX}/Library/lib/mkl_intel_thread.lib"
        "${MICROMAMBA_ENV_PREFIX}/Library/lib/mkl_core.lib"
        "${MICROMAMBA_ENV_PREFIX}/Library/lib/mkl_rt.lib"
        "${MICROMAMBA_ENV_PREFIX}/Library/lib/gtest.lib"
    DESTINATION lib
)