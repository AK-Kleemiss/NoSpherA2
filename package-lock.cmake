# CPM Package Lock
# This file should be committed to version control

# Ccache.cmake
CPMDeclarePackage(Ccache.cmake
  NAME Ccache.cmake
  VERSION 1.2.5
  GITHUB_REPOSITORY TheLartians/Ccache.cmake
)
# mdspan (unversioned)
 CPMDeclarePackage(mdspan
  NAME mdspan
  GIT_TAG "546d4dd63697c6a331554adb6fe650e09b690812"
  GITHUB_REPOSITORY kokkos/mdspan
)
# featomic
CPMDeclarePackage(featomic
  NAME featomic
  GITHUB_REPOSITORY "metatensor/featomic"
  GIT_TAG "featomic-v0.6.4"
  SOURCE_SUBDIR "featomic"
  OPTIONS
    "G Ninja"
    "CMAKE_BUILD_TYPE Release"
    "FEATOMIC_FETCH_METATENSOR ON"
    "BUILD_SHARED_LIBS OFF"
    "FEATOMIC_INSTALL_BOTH_STATIC_SHARED OFF"
    "FEATOMIC_USE_STATIC_METATENSOR ON"
)
# occ (unversioned)
# CPMDeclarePackage(occ
#  NAME occ
#  GIT_TAG libs_rebase2
#  GITHUB_REPOSITORY MilitaoLucas/occ
#  OPTIONS
#    "CMAKE_BUILD_TYPE Release"
#    "WITH_PYTHON_BINDINGS OFF"
#    "CMAKE_INSTALL_PREFIX /home/lucas/CLionProjects/NoSpherA2/Lib/occ_install"
#    "CMAKE_CXX_FLAGS -O3 -ffast-math"
#    "CMAKE_C_FLAGS -O3 -ffast-math"
#    "BLA_VENDOR Intel10_64lp"
#    "USE_SYSTEM_TBB OFF"
#    "TBB_PREFIX_PATH /home/lucas/CLionProjects/NoSpherA2/cmake-build-release-gcc/environments/mambaenv"
#    "BLAS_ROOT /home/lucas/CLionProjects/NoSpherA2/cmake-build-release-gcc/environments/mambaenv/lib/cmake/mkl"
#    "USE_QCINT OFF"
#    "ENABLE_HOST_OPT OFF"
#    "USE_FORTRAN OFF"
#    "BLA_STATIC ON"
#)
# fmt
CPMDeclarePackage(fmt
  NAME fmt
  VERSION 11.1.0
  GIT_TAG 11.1.0
  GITHUB_REPOSITORY fmtlib/fmt
)
# spdlog
CPMDeclarePackage(spdlog
  NAME spdlog
  VERSION 1.x
  GITHUB_REPOSITORY gabime/spdlog
  OPTIONS
    "SPDLOG_FMT_EXTERNAL ON"
)
# tomlplusplus
CPMDeclarePackage(tomlplusplus
  NAME tomlplusplus
  VERSION 3.4.0
  GITHUB_REPOSITORY marzer/tomlplusplus
)
# oneTBB
CPMDeclarePackage(oneTBB
  NAME oneTBB
  VERSION 2022.2.0
  GITHUB_REPOSITORY uxlfoundation/oneTBB
  OPTIONS
    "TBB_TEST OFF"
    "TBB_STRICT OFF"
    "TBB_DISABLE_HWLOC_AUTOMATIC_SEARCH OFF"
)
# unordered_dense
CPMDeclarePackage(unordered_dense
  NAME unordered_dense
  VERSION 4.5.0
  GITHUB_REPOSITORY martinus/unordered_dense
)
# scnlib
CPMDeclarePackage(scnlib
  NAME scnlib
  VERSION 4.0.1
  GITHUB_REPOSITORY eliaskosunen/scnlib
  OPTIONS
    "SIMDUTF_TESTS OFF"
    "ENABLE_FULL OFF"
)
# CLI11
CPMDeclarePackage(CLI11
  NAME CLI11
  VERSION 2.4.2
  GITHUB_REPOSITORY CLIUtils/CLI11
)
# nlohmann_json
CPMDeclarePackage(nlohmann_json
  NAME nlohmann_json
  VERSION 3.11.3
  GITHUB_REPOSITORY nlohmann/json
)
# eigen3
CPMDeclarePackage(eigen3
  NAME eigen3
  DOWNLOAD_ONLY YES
  URL
    "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip"
)
# gemmi
CPMDeclarePackage(gemmi
  NAME gemmi
  VERSION 0.6.5
  DOWNLOAD_ONLY YES
  GITHUB_REPOSITORY project-gemmi/gemmi
)
# dftd4_cpp
CPMDeclarePackage(dftd4_cpp
  NAME dftd4_cpp
  VERSION 2.2.0
  GIT_TAG main
  GITHUB_REPOSITORY peterspackman/cpp-d4
  OPTIONS
    "DFTD4_USE_EIGEN ON"
    "BUILD_SHARED_LIBS OFF"
)
# LBFGSpp (unversioned)
# CPMDeclarePackage(LBFGSpp
#  NAME LBFGSpp
#  GIT_TAG master
#  DOWNLOAD_ONLY YES
#  GITHUB_REPOSITORY yixuan/LBFGSpp
#)
# libcint
if(MSVC)
  set(LIBCINT_REPO "MilitaoLucas/libcint")
else()
  set(LIBCINT_REPO "peterspackman/libcint")
endif ()
CPMDeclarePackage(libcint
  NAME libcint
  VERSION 6.1.2
  GIT_TAG master
  GITHUB_REPOSITORY ${LIBCINT_REPO}
  OPTIONS
    "WITH_FORTRAN OFF"
    "WITH_CINT2_INTERFACE OFF"
    "ENABLE_STATIC ON"
    "BUILD_SHARED_LIBS OFF"
    "PYPZPX ON"
    "WITH_RANGE_COULOMB ON"
    "CMAKE_C_FLAGS -Wno-implicit-function-declaration -Wno-deprecated-non-prototype -D_GNU_SOURCE"
)
# Libxc
CPMDeclarePackage(Libxc
  NAME Libxc
  URL
    "https://gitlab.com/libxc/libxc/-/archive/6.2.2/libxc-6.2.2.tar.gz"
  OPTIONS
    "CMAKE_BUILD_TYPE Release"
    "BUILD_TESTING OFF"
    "ENABLE_XHOST OFF"
    "BUILD_FPIC OFF"
    "CMAKE_POLICY_VERSION_MINIMUM 3.5"
)
