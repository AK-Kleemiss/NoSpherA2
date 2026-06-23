cmake_minimum_required(VERSION 3.25)

get_filename_component(
    NOSPHERA2_SOURCE_DIR
    "${CMAKE_CURRENT_LIST_DIR}/.."
    ABSOLUTE
)

include(
    "${NOSPHERA2_SOURCE_DIR}/cmake/MicromambaEnvironment.cmake"
)

setup_micromamba_environment(
    ENVIRONMENT_FILE
        "${NOSPHERA2_SOURCE_DIR}/environment.yaml"
    PREFIX
        "${NOSPHERA2_SOURCE_DIR}/.mambaenv/env"
    ROOT_PREFIX
        "${NOSPHERA2_SOURCE_DIR}/.mambaenv/root"
    DOWNLOAD_DIRECTORY
        "${NOSPHERA2_SOURCE_DIR}/.mambaenv/bootstrap"
)

message(STATUS "Bootstrap complete")
message(STATUS "Environment: ${MICROMAMBA_ENV_PREFIX}")
message(STATUS "Configure with: cmake --preset release-xxxx")
message(STATUS "Build with: cmake --build --preset release-xxxx")