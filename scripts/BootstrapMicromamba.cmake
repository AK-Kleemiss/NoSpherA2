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
# ALL sets up the fat CUDA + HIP build the CI ships: the CUDA toolkit plus AMD's ROCm SDK
# wheels (hipBLAS included) in .mambaenv/rocm, unless a ROCm is installed already.
set(NOSPHERA2_BOOTSTRAP_GPU_VENDOR "AUTO" CACHE STRING "GPU toolkit to fetch: AUTO, NVIDIA, AMD or ALL")
set(NOSPHERA2_BOOTSTRAP_CUDA_VERSION "" CACHE STRING "cuda-version to pin when fetching the CUDA toolkit, e.g. 12.9")
set(NOSPHERA2_BOOTSTRAP_ROCM_VERSION "10.0.0" CACHE STRING "rocm[devel] wheel version for NOSPHERA2_BOOTSTRAP_GPU_VENDOR=ALL (the CI's)")

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
    # Windows on ARM. On an ARM64 host the one environment is win-arm64 (OpenBLAS, no MKL). On an
    # x64 host it cross-builds: host tools stay in env, the ARM64 libraries go to env-win-arm64 and
    # env's rustc gets the aarch64 standard library.
    if(WIN32 AND "$ENV{PROCESSOR_ARCHITECTURE}" STREQUAL "ARM64")
        set(_win_arm64_default ON)
    else()
        set(_win_arm64_default OFF)
    endif()
    option(NOSPHERA2_BOOTSTRAP_WIN_ARM64 "Set up for a Windows ARM64 build (native or cross from x64)" ${_win_arm64_default})

    cmake_host_system_information(RESULT _host_processor QUERY OS_PLATFORM)
    if(NOSPHERA2_BOOTSTRAP_WIN_ARM64 AND "$ENV{PROCESSOR_ARCHITECTURE}" STREQUAL "ARM64")
        set(_environment_file
            "${NOSPHERA2_SOURCE_DIR}/environment-win-arm64.yaml"
        )
    elseif(CMAKE_HOST_LINUX AND _host_processor MATCHES "aarch64|arm64")
        set(_environment_file
            "${NOSPHERA2_SOURCE_DIR}/environment-linux-aarch64.yaml"
        )
    else()
        set(_environment_file
            "${NOSPHERA2_SOURCE_DIR}/environment.yaml"
        )
    endif()

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

    if(NOSPHERA2_BOOTSTRAP_WIN_ARM64 AND NOT "$ENV{PROCESSOR_ARCHITECTURE}" STREQUAL "ARM64")
        setup_micromamba_environment(
            ENVIRONMENT_FILE "${NOSPHERA2_SOURCE_DIR}/environment-win-arm64.yaml"
            PLATFORM "win-arm64"
            PREFIX "${NOSPHERA2_SOURCE_DIR}/.mambaenv/env-win-arm64"
            ROOT_PREFIX "${_mamba_root}"
            DOWNLOAD_DIRECTORY "${_mamba_bootstrap}"
            EXPORT_PREFIX_VARIABLE MICROMAMBA_ENV_WIN_ARM64_PREFIX
        )
        # rust-std is noarch; it has to match env's rustc exactly
        execute_process(
            COMMAND "${MICROMAMBA_ENV_PREFIX}/Library/bin/rustc.exe" --version
            OUTPUT_VARIABLE _rustc_version
            COMMAND_ERROR_IS_FATAL ANY
        )
        string(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+" _rustc_version "${_rustc_version}")
        execute_process(
            COMMAND
                "${CMAKE_COMMAND}" -E env
                "MAMBA_ROOT_PREFIX=${MICROMAMBA_ROOT_PREFIX}"
                "${MICROMAMBA_EXECUTABLE}"
                install
                --yes
                --prefix "${MICROMAMBA_ENV_PREFIX}"
                -c conda-forge
                "rust-std-aarch64-pc-windows-msvc=${_rustc_version}"
            COMMAND_ERROR_IS_FATAL ANY
        )
        message(STATUS "Windows ARM64 libraries: ${MICROMAMBA_ENV_WIN_ARM64_PREFIX}")
    endif()

    # macOS is deliberately excluded above: no CUDA, and no AMD compute stack either.
    if(NOSPHERA2_BOOTSTRAP_GPU)
        if(NOSPHERA2_BOOTSTRAP_GPU_VENDOR STREQUAL "AUTO")
            nosphera2_detect_gpu(_gpu_vendor)
        else()
            set(_gpu_vendor "${NOSPHERA2_BOOTSTRAP_GPU_VENDOR}")
            message(STATUS "GPU toolkit chosen by NOSPHERA2_BOOTSTRAP_GPU_VENDOR: ${_gpu_vendor}")
        endif()
        if(_gpu_vendor STREQUAL "ALL")
            nosphera2_bootstrap_cuda_toolkit(
                PREFIX      "${MICROMAMBA_ENV_PREFIX}"
                ROOT_PREFIX "${MICROMAMBA_ROOT_PREFIX}"
                EXECUTABLE  "${MICROMAMBA_EXECUTABLE}"
                VERSION     "${NOSPHERA2_BOOTSTRAP_CUDA_VERSION}"
            )
            nosphera2_find_rocm(_rocm_hipcc)
            if(_rocm_hipcc)
                message(STATUS "ROCm already installed: ${_rocm_hipcc}")
            else()
                if(WIN32)
                    set(_env_python "${MICROMAMBA_ENV_PREFIX}/python.exe")
                else()
                    set(_env_python "${MICROMAMBA_ENV_PREFIX}/bin/python")
                endif()
                nosphera2_bootstrap_rocm_sdk(
                    PYTHON    "${_env_python}"
                    DIRECTORY "${NOSPHERA2_SOURCE_DIR}/.mambaenv/rocm"
                    VERSION   "${NOSPHERA2_BOOTSTRAP_ROCM_VERSION}"
                )
            endif()
        elseif(_gpu_vendor STREQUAL "NVIDIA")
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

    # Arm Performance Libraries (Linux aarch64, Windows ARM64). Arm's site refuses scripted
    # downloads, so the package is fetched by hand (https://learn.arm.com/install-guides/armpl/)
    # and named in NOSPHERA2_ARMPL_PACKAGE or put into .mambaenv/bootstrap. Installing it accepts
    # Arm's licence. It goes to .mambaenv/armpl, whose root.txt the top-level CMakeLists.txt reads.
    set(NOSPHERA2_ARMPL_PACKAGE "$ENV{NOSPHERA2_ARMPL_PACKAGE}" CACHE FILEPATH "Arm Performance Libraries package (.tar, .sh or .msi) to install")
    if(NOT NOSPHERA2_ARMPL_PACKAGE)
        file(GLOB _armpl_packages "${_mamba_bootstrap}/arm-performance-libraries_*")
        list(FILTER _armpl_packages INCLUDE REGEX "\\.(tar|sh|msi)$")
        # the Arm64EC .msi is for x64-compatible ARM code, not a native ARM64 link
        list(FILTER _armpl_packages EXCLUDE REGEX "Arm64EC")
        list(SORT _armpl_packages COMPARE NATURAL ORDER DESCENDING)
        list(POP_FRONT _armpl_packages NOSPHERA2_ARMPL_PACKAGE)
    endif()
    set(_armpl_dir "${NOSPHERA2_SOURCE_DIR}/.mambaenv/armpl")
    if(NOSPHERA2_ARMPL_PACKAGE AND NOT EXISTS "${_armpl_dir}/root.txt")
        message(STATUS "Installing ${NOSPHERA2_ARMPL_PACKAGE} into ${_armpl_dir}, which accepts Arm's licence")
        if(NOSPHERA2_ARMPL_PACKAGE MATCHES "\\.msi$")
            # administrative install: unpacks the files without touching Program Files or the registry
            file(TO_NATIVE_PATH "${NOSPHERA2_ARMPL_PACKAGE}" _armpl_msi)
            file(TO_NATIVE_PATH "${_armpl_dir}" _armpl_target)
            execute_process(
                COMMAND msiexec /a "${_armpl_msi}" /qn "TARGETDIR=${_armpl_target}" ACCEPT_EULA=1
                COMMAND_ERROR_IS_FATAL ANY
            )
        else()
            set(_armpl_installer "${NOSPHERA2_ARMPL_PACKAGE}")
            if(_armpl_installer MATCHES "\\.tar$")
                file(ARCHIVE_EXTRACT INPUT "${_armpl_installer}" DESTINATION "${_mamba_bootstrap}/armpl_package")
                file(GLOB_RECURSE _armpl_installer "${_mamba_bootstrap}/armpl_package/arm-performance-libraries_*.sh")
            endif()
            # it unpacks ~1 GB into TMPDIR, more than a small machine's tmpfs /tmp holds (Pi 4: 922 MB);
            # --force lets a rerun install over an attempt that died half way (no root.txt yet)
            execute_process(
                COMMAND "${CMAKE_COMMAND}" -E env "TMPDIR=${_mamba_bootstrap}"
                    bash "${_armpl_installer}" --accept --force --install-to "${_armpl_dir}"
                COMMAND_ERROR_IS_FATAL ANY
            )
        endif()
        file(GLOB_RECURSE _armpl_headers "${_armpl_dir}/armpl.h")
        list(FILTER _armpl_headers INCLUDE REGEX "/include/armpl\\.h$")
        list(POP_FRONT _armpl_headers _armpl_header)
        if(NOT _armpl_header)
            message(FATAL_ERROR "No include/armpl.h under ${_armpl_dir} after installing ${NOSPHERA2_ARMPL_PACKAGE}")
        endif()
        get_filename_component(_armpl_root "${_armpl_header}" DIRECTORY)
        get_filename_component(_armpl_root "${_armpl_root}" DIRECTORY)
        file(WRITE "${_armpl_dir}/root.txt" "${_armpl_root}")
        message(STATUS "Arm Performance Libraries: ${_armpl_root}")
    endif()

    message(STATUS "Bootstrap complete")
    message(STATUS "Environment: ${MICROMAMBA_ENV_PREFIX}")
    message(STATUS "Configure with: cmake --preset release-xxxx")
    message(STATUS "Build with: cmake --build --preset release-xxxx")
endif()