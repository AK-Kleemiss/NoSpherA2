include(ExternalProject)
CPMAddPackage(
    NAME Ccache.cmake
    GITHUB_REPOSITORY TheLartians/Ccache.cmake
    VERSION 1.2.5
)

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
        "BUILD_SHARED_LIBS ON"
        "METATENSOR_INSTALL_BOTH_STATIC_SHARED ON"
)
message(WARNING "${CMAKE_INSTALL_LIBDIR}")
list(APPEND CMAKE_PREFIX_PATH "${metatensor_BINARY_DIR}")
CPMAddPackage(
    NAME featomic
    URL "https://github.com/metatensor/featomic/releases/download/featomic-v0.6.1/featomic-cxx-0.6.1.tar.gz"
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
    file(MAKE_DIRECTORY ${featomic_SOURCE_DIR}/.cargo)
    file(WRITE ${featomic_SOURCE_DIR}/.cargo/config.toml "paths = ['${metatensor_SOURCE_DIR}']")
endif ()
add_custom_target(featomic_target DEPENDS featomic_static featomic)
set(MKL_INSTALL_DIR ${CMAKE_SOURCE_DIR}/Lib/MKL)
message(WARNING "HERE: ${ONEAPI_DIRECTORY}")
message(WARNING "HERE: ${CMAKE_CURRENT_BINARY_DIR}")
if (${FAKE_HOME_TEST})
    set(HOME "/tmp/home")
else ()
    set(HOME $ENV{HOME})
endif ()

set(MKL_LINK "static")
if (NOT EXISTS ${ONEAPI_DIRECTORY})
    set(_MKL_ROOT "${MKL_INSTALL_DIR}/mkl/latest")
    set(_MKL_CONFIG "${_MKL_ROOT}/lib/cmake/mkl/MKLConfig.cmake")

    if(NOT EXISTS "${_MKL_CONFIG}")
        file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/Lib/MKL")
        set(_MKL_SH "${CMAKE_SOURCE_DIR}/Lib/intel-mkl-installer.sh")
        if (NOT EXISTS ${_MKL_SH})
            file(DOWNLOAD
              "https://registrationcenter-download.intel.com/akdlm/IRC_NAS/47c7d946-fca1-441a-b0df-b094e3f045ea/intel-onemkl-2025.2.0.629_offline.sh"
              "${_MKL_SH}"
              EXPECTED_HASH SHA256=7bfdde1379a9cd3c2784af4352931a0527a5c483c9da71e851a3f06e055af7c7
              SHOW_PROGRESS
              STATUS _DL_STATUS
            )
            list(GET _DL_STATUS 0 _DL_CODE)
            if(NOT _DL_CODE EQUAL 0)
              message(FATAL_ERROR "Failed to download oneMKL: ${_DL_STATUS}")
            endif()
        endif ()
        set(OLD_HOME $ENV{HOME})
        set(ENV{HOME} "/tmp/home")
        message(STATUS "INSTALLING MKL IN: ${MKL_INSTALL_DIR}")
        file(MAKE_DIRECTORY $ENV{HOME})
        execute_process(COMMAND chmod +x ${_MKL_SH})
        execute_process(
            COMMAND ${_MKL_SH} -a --silent --eula accept --install-dir "${MKL_INSTALL_DIR}"

          RESULT_VARIABLE _MKL_RC
        )
        FILE(REMOVE $ENV{HOME})
        set(ENV{HOME} ${OLD_HOME})
        if(NOT _MKL_RC EQUAL 0)
          message(FATAL_ERROR "oneMKL installer failed with code ${_MKL_RC}")
        endif()
    endif()
    list(APPEND CMAKE_PREFIX_PATH "${_MKL_ROOT}")
    find_package(MKL CONFIG REQUIRED
        HINTS "${_MKL_ROOT}"
        PATH_SUFFIXES "lib/cmake/mkl"
    )
    message(STATUS "Imported oneMKL targets: ${MKL_IMPORTED_TARGETS}")

else ()
    set(_ONEAPI_MKL_ROOT "${ONEAPI_DIRECTORY}/mkl/latest")
    find_package(MKL CONFIG REQUIRED
        HINTS "${_ONEAPI_MKL_ROOT}" "${ONEAPI_DIRECTORY}"
        PATH_SUFFIXES "lib/cmake/mkl"
    )
  message(STATUS "Imported oneMKL targets: ${MKL_IMPORTED_TARGETS}")
endif ()

