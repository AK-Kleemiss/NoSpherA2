cmake_minimum_required(VERSION 3.25)

get_filename_component(
    NOSPHERA2_SOURCE_DIR
    "${CMAKE_CURRENT_LIST_DIR}/.."
    ABSOLUTE
)

include(
    "${CMAKE_CURRENT_LIST_DIR}/MicromambaEnvironment.cmake"
)

include(
    "${CMAKE_CURRENT_LIST_DIR}/GpuToolkit.cmake"
)

# ON by default so a machine with a card gets the GPU path without being asked. Set it to OFF
# to keep the bootstrap to the packages every build needs - worth doing on a metered
# connection, since the CUDA packages are the largest thing here by a wide margin.
option(NOSPHERA2_BOOTSTRAP_GPU "Add a CUDA toolkit to the environment when an NVIDIA GPU is present" ON)
# AUTO looks at the driver. NVIDIA or AMD fetches that toolkit whether or not a card is
# present - for the CI runners and for a head node that builds for the compute nodes. The
# CUDA version is only used when one is fetched; the CI pins it so the artifact covers the
# same cards from one week to the next (see the arch list in the top-level CMakeLists.txt).
set(NOSPHERA2_BOOTSTRAP_GPU_VENDOR "AUTO" CACHE STRING "GPU toolkit to fetch: AUTO, NVIDIA or AMD")
set(NOSPHERA2_BOOTSTRAP_CUDA_VERSION "" CACHE STRING "cuda-version to pin when fetching the CUDA toolkit, e.g. 12.9")

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

    # macOS is deliberately excluded above: no CUDA, and no AMD compute stack either.
    if(NOSPHERA2_BOOTSTRAP_GPU)
        if(NOSPHERA2_BOOTSTRAP_GPU_VENDOR STREQUAL "AUTO")
            nosphera2_detect_gpu(_gpu_vendor)
        else()
            set(_gpu_vendor "${NOSPHERA2_BOOTSTRAP_GPU_VENDOR}")
            message(STATUS "GPU toolkit chosen by NOSPHERA2_BOOTSTRAP_GPU_VENDOR: ${_gpu_vendor}")
        endif()
        if(_gpu_vendor STREQUAL "NVIDIA")
            nosphera2_bootstrap_cuda_toolkit(
                PREFIX      "${MICROMAMBA_ENV_PREFIX}"
                ROOT_PREFIX "${MICROMAMBA_ROOT_PREFIX}"
                EXECUTABLE  "${MICROMAMBA_EXECUTABLE}"
                VERSION     "${NOSPHERA2_BOOTSTRAP_CUDA_VERSION}"
            )
        elseif(_gpu_vendor STREQUAL "AMD")
            nosphera2_find_rocm(_rocm_hipcc)
            if(_rocm_hipcc)
                message(STATUS "ROCm already installed: ${_rocm_hipcc}")
            else()
                nosphera2_bootstrap_rocm_toolkit(
                    PREFIX      "${MICROMAMBA_ENV_PREFIX}"
                    ROOT_PREFIX "${MICROMAMBA_ROOT_PREFIX}"
                    EXECUTABLE  "${MICROMAMBA_EXECUTABLE}"
                )
            endif()
        else()
            message(STATUS "No GPU driver detected; building for the CPU")
        endif()
    endif()

    message(STATUS "Bootstrap complete")
    message(STATUS "Environment: ${MICROMAMBA_ENV_PREFIX}")
    message(STATUS "Configure with: cmake --preset release-xxxx")
    message(STATUS "Build with: cmake --build --preset release-xxxx")
endif()