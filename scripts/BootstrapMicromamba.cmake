cmake_minimum_required(VERSION 3.25)

get_filename_component(
    NOSPHERA2_SOURCE_DIR
    "${CMAKE_CURRENT_LIST_DIR}/.."
    ABSOLUTE
)

include(
    "${CMAKE_CURRENT_LIST_DIR}/MicromambaEnvironment.cmake"
)

set(_mamba_root
    "${NOSPHERA2_SOURCE_DIR}/.mambaenv/root"
)

set(_mamba_bootstrap
    "${NOSPHERA2_SOURCE_DIR}/.mambaenv/bootstrap"
)

if(APPLE)
    set(_environment_file
        "${NOSPHERA2_SOURCE_DIR}/environment-macos.yaml"
    )

    setup_micromamba_environment(
        ENVIRONMENT_FILE "${_environment_file}"
        PLATFORM "osx-arm64"
        PREFIX "${NOSPHERA2_SOURCE_DIR}/.mambaenv/env-arm64"
        ROOT_PREFIX "${_mamba_root}"
        DOWNLOAD_DIRECTORY "${_mamba_bootstrap}"
        EXPORT_PREFIX_VARIABLE MICROMAMBA_ENV_ARM64_PREFIX
    )


    setup_micromamba_environment(
        ENVIRONMENT_FILE "${_environment_file}"
        PLATFORM "osx-64"
        PREFIX "${NOSPHERA2_SOURCE_DIR}/.mambaenv/env-x86_64"
        ROOT_PREFIX "${_mamba_root}"
        DOWNLOAD_DIRECTORY "${_mamba_bootstrap}"
        EXPORT_PREFIX_VARIABLE MICROMAMBA_ENV_X86_64_PREFIX
    )

    message(STATUS "Bootstrap complete")
    message(STATUS
        "macOS arm64 environment: "
        "${NOSPHERA2_SOURCE_DIR}/.mambaenv/env-arm64"
    )
    message(STATUS
        "macOS x86_64 environment: "
        "${NOSPHERA2_SOURCE_DIR}/.mambaenv/env-x86_64"
    )

else()
    set(_environment_file
        "${NOSPHERA2_SOURCE_DIR}/environment.yaml"
    )

    setup_micromamba_environment(
        ENVIRONMENT_FILE
            "${_environment_file}"
        PREFIX
            "${NOSPHERA2_SOURCE_DIR}/.mambaenv/env"
        ROOT_PREFIX
            "${_mamba_root}"
        DOWNLOAD_DIRECTORY
            "${_mamba_bootstrap}"
    )

    message(STATUS "Bootstrap complete")
    message(STATUS "Environment: ${MICROMAMBA_ENV_PREFIX}")
    message(STATUS "Configure with: cmake --preset release-xxxx")
    message(STATUS "Build with: cmake --build --preset release-xxxx")
endif()