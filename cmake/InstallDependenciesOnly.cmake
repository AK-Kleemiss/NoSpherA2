enable_googletests()

add_custom_target(NoSpherA2WindowsDependencies ALL
    DEPENDS
        cint
        tbb
        occ_main
        gtest
        BasisSetConverter
        cargo-build-featomic
)

if(NOT WIN32)
    return()
endif()

install(
    TARGETS
        cint
        tbb
        spdlog
        gau2grid_static
        fmt
        xc
        scn
        _subprocess
        libecpint_static
        gtest
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)

install(
    TARGETS
    occ_cc
    occ_cg
    occ_core
    occ_correlation
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
    occ_mults
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

# occ's public headers include nlohmann/json.hpp (occ/dft/dft_method.h)
install(
    DIRECTORY "${nlohmann_json_SOURCE_DIR}/include/"
    DESTINATION include
)

# and xc.h (occ/dft/functional.h); xc_version.h is generated
install(
    FILES
        "${Libxc_SOURCE_DIR}/src/xc.h"
        "${Libxc_SOURCE_DIR}/src/xc_funcs.h"
        "${Libxc_SOURCE_DIR}/src/xc_funcs_removed.h"
        "${Libxc_BINARY_DIR}/xc_version.h"
    DESTINATION include
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

install(
    DIRECTORY
        "${googletest_SOURCE_DIR}/googletest/include/"
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

# The imported targets name the libraries of this build's RUST_BUILD_TARGET; a glob
# over target/ also picks up a stale host (x64) build left in an ARM64 tree
install(
    FILES
        $<TARGET_FILE:featomic::static>
        $<TARGET_FILE:metatensor::static>
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

if(NOSPHERA2_ARMPL)
    # armpl_lp64.lib (static, Release) or armpl_lp64.dll.lib (Debug, DLL to bin) in lib is what
    # switches the vcxprojs to ArmPL (NosArmPL in NoSpherA2_universal.props); the flang runtime
    # libraries ride along with the static one.
    install(FILES ${_armpl_libs} DESTINATION lib)
    if(NOSPHERA2_ARMPL_DLL)
        install(FILES "${NOSPHERA2_ARMPL_DLL}" DESTINATION bin)
    endif()
    install(DIRECTORY "${NOSPHERA2_ARMPL_INCLUDE_DIR}/" DESTINATION include)
elseif(NOSPHERA2_OPENBLAS)
    # back from ArmPL: drop the markers so the vcxprojs link OpenBLAS again
    install(CODE "file(REMOVE \"\${CMAKE_INSTALL_PREFIX}/lib/armpl_lp64.lib\" \"\${CMAKE_INSTALL_PREFIX}/lib/armpl_lp64.dll.lib\")")
    # The vcxprojs link openblas.lib whatever conda-forge named the import library
    install(FILES "${NOSPHERA2_OPENBLAS_DLL}" DESTINATION bin)
    install(FILES "${NOSPHERA2_OPENBLAS_IMPLIB}" DESTINATION lib RENAME openblas.lib)
    install(
        DIRECTORY "${NOSPHERA2_OPENBLAS_INCLUDE_DIR}/"
        DESTINATION include
        FILES_MATCHING
            PATTERN "cblas*.h"
            PATTERN "lapack*.h"
            PATTERN "openblas*.h"
            PATTERN "f77blas.h"
    )
else()
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
    DESTINATION lib
)
endif()