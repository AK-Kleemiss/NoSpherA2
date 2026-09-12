# Find out whether this machine has a GPU worth building for, and get a CUDA toolkit onto it
# if it does.
#
# The detection has to work with no toolkit installed - that is the whole point, since its
# answer decides whether to fetch one - so it looks for the *driver*. A driver is what a user
# with a graphics card actually has; the toolkit is the thing they have not got round to
# installing. nvcuda.dll and libcuda.so.1 ship with the NVIDIA driver and with nothing else.
#
# Both vendors are fetched, but they are not equally well served and it is worth knowing why
# before trusting the AMD side.
#
# CUDA is complete: conda-forge ships the compiler, the runtime and cuBLAS.
#
# ROCm is partial. conda-forge has hipcc and hip-devel for linux-64 only, and no hipBLAS at
# all - checked against the package index and then confirmed by actually creating an
# environment, which produced hip_runtime.h and no hipblas.h. AMD's own conda channel at
# repo.radeon.com/rocm/conda is empty (its repodata.json 404s), and everything else AMD
# distributes is an apt/rpm repository or a Windows installer wanting administrator rights.
# So on linux-64 we fetch what exists and build the three kernels that need no BLAS - which
# includes the scattering-factor transform, the one path that runs without being asked for -
# and the I tensor and the offloaded GEMM stay on the CPU. Everywhere else ROCm has to be
# installed by hand, and this detects it and uses it.

# NVIDIA, AMD or NONE, from the driver rather than from any toolkit.
function(nosphera2_detect_gpu out_var)
    set(_vendor "NONE")

    if(WIN32)
        if(EXISTS "$ENV{SystemRoot}/System32/nvcuda.dll")
            set(_vendor "NVIDIA")
        elseif(EXISTS "$ENV{SystemRoot}/System32/amdhip64.dll"
            OR EXISTS "$ENV{SystemRoot}/System32/amdhip64_6.dll")
            set(_vendor "AMD")
        endif()
    else()
        if(EXISTS "/proc/driver/nvidia/version")
            set(_vendor "NVIDIA")
        elseif(EXISTS "/sys/module/amdgpu")
            set(_vendor "AMD")
        else()
            foreach(_dir "/usr/lib/x86_64-linux-gnu" "/usr/lib64" "/usr/lib" "/lib/x86_64-linux-gnu")
                if(EXISTS "${_dir}/libcuda.so.1")
                    set(_vendor "NVIDIA")
                    break()
                endif()
            endforeach()
        endif()
    endif()

    # A container or an unusual driver layout can hide the library while nvidia-smi still
    # works, so fall back to asking it. Only reached when the cheap checks found nothing.
    if(_vendor STREQUAL "NONE")
        find_program(NOSPHERA2_NVIDIA_SMI NAMES nvidia-smi)
        mark_as_advanced(NOSPHERA2_NVIDIA_SMI)
        if(NOSPHERA2_NVIDIA_SMI)
            execute_process(
                COMMAND "${NOSPHERA2_NVIDIA_SMI}" -L
                RESULT_VARIABLE _smi_result
                OUTPUT_VARIABLE _smi_output
                ERROR_VARIABLE  _smi_error
                TIMEOUT 10
            )
            if(_smi_result EQUAL 0 AND _smi_output MATCHES "GPU [0-9]+:")
                set(_vendor "NVIDIA")
            endif()
        endif()
    endif()

    set(${out_var} "${_vendor}" PARENT_SCOPE)
endfunction()

# An existing ROCm, wherever the platform puts it. Checked before fetching anything, and on
# Windows it is the only way to get one.
function(nosphera2_find_rocm out_var)
    set(_roots "")
    if(DEFINED ENV{HIP_PATH})
        list(APPEND _roots "$ENV{HIP_PATH}")
    endif()
    if(DEFINED ENV{ROCM_PATH})
        list(APPEND _roots "$ENV{ROCM_PATH}")
    endif()
    list(APPEND _roots "/opt/rocm" "C:/Program Files/AMD/ROCm")
    foreach(_root ${_roots})
        foreach(_exe "bin/hipcc" "bin/hipcc.exe" "bin/hipcc.bat")
            if(EXISTS "${_root}/${_exe}")
                set(${out_var} "${_root}/${_exe}" PARENT_SCOPE)
                return()
            endif()
        endforeach()
        # The Windows SDK installs under a version directory, so look one level down
        file(GLOB _versioned "${_root}/*/bin/hipcc.exe" "${_root}/*/bin/hipcc.bat")
        if(_versioned)
            list(GET _versioned 0 _first)
            set(${out_var} "${_first}" PARENT_SCOPE)
            return()
        endif()
    endforeach()
    set(${out_var} "" PARENT_SCOPE)
endfunction()

# Where hipcc would live inside a micromamba environment.
function(nosphera2_env_hipcc env_prefix out_var)
    if(EXISTS "${env_prefix}/bin/hipcc")
        set(${out_var} "${env_prefix}/bin/hipcc" PARENT_SCOPE)
    else()
        set(${out_var} "" PARENT_SCOPE)
    endif()
endfunction()

# Add what conda-forge has of ROCm to the environment. linux-64 only, and no hipBLAS - see
# the note at the top of this file. Non-fatal for the same reason the CUDA one is.
function(nosphera2_bootstrap_rocm_toolkit)
    cmake_parse_arguments(GPU "" "PREFIX;ROOT_PREFIX;EXECUTABLE" "" ${ARGN})

    nosphera2_env_hipcc("${GPU_PREFIX}" _existing)
    if(_existing)
        message(STATUS "HIP toolchain already in the environment: ${_existing}")
        return()
    endif()

    if(NOT CMAKE_HOST_SYSTEM_NAME STREQUAL "Linux"
       OR NOT CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
        message(STATUS
            "AMD GPU detected. conda-forge packages hipcc for linux-64 only, so nothing can\n"
            "  be fetched here - install the HIP SDK and configure will pick it up.")
        return()
    endif()

    message(STATUS "AMD GPU found and no ROCm present - fetching the HIP toolchain (~1 GB)")
    execute_process(
        COMMAND
            "${CMAKE_COMMAND}" -E env
            "MAMBA_ROOT_PREFIX=${GPU_ROOT_PREFIX}"
            "${GPU_EXECUTABLE}"
            install
            --yes
            --prefix "${GPU_PREFIX}"
            -c conda-forge
            hipcc hip-devel
        RESULT_VARIABLE _rocm_result
        OUTPUT_VARIABLE _rocm_output
        ERROR_VARIABLE  _rocm_error
    )

    if(NOT _rocm_result EQUAL 0)
        message(STATUS
            "Could not add a HIP toolchain; building without the GPU path.\n"
            "  This is not an error - NoSpherA2 runs on the CPU.\n"
            "${_rocm_error}")
        return()
    endif()

    nosphera2_env_hipcc("${GPU_PREFIX}" _installed)
    if(_installed)
        message(STATUS
            "HIP toolchain ready: ${_installed}\n"
            "  conda-forge has no hipBLAS, so the I tensor and the offloaded GEMM stay on the\n"
            "  CPU. Install the full ROCm stack if you want those on the device too.")
    else()
        message(STATUS "HIP packages installed but no hipcc appeared; building without the GPU path")
    endif()
endfunction()

# Where nvcc would live inside a micromamba environment, which differs by platform.
function(nosphera2_env_nvcc env_prefix out_var)
    if(WIN32)
        set(_candidates "${env_prefix}/Library/bin/nvcc.exe" "${env_prefix}/bin/nvcc.exe")
    else()
        set(_candidates "${env_prefix}/bin/nvcc")
    endif()
    foreach(_c ${_candidates})
        if(EXISTS "${_c}")
            set(${out_var} "${_c}" PARENT_SCOPE)
            return()
        endif()
    endforeach()
    set(${out_var} "" PARENT_SCOPE)
endfunction()

# Add the CUDA compiler and the two libraries the kernels link against to an existing
# micromamba environment.
#
# Deliberately a separate step from creating that environment rather than three more lines in
# environment.yaml. The environment is what every build needs and its creation is fatal if it
# fails; this is optional, and on a machine whose driver is too old for any packaged toolkit
# the solve will fail. Folding it into the base environment would turn "your driver is old"
# into "you cannot build NoSpherA2 at all". Here a failure costs the GPU path and nothing else.
function(nosphera2_bootstrap_cuda_toolkit)
    cmake_parse_arguments(GPU "" "PREFIX;ROOT_PREFIX;EXECUTABLE" "" ${ARGN})

    nosphera2_env_nvcc("${GPU_PREFIX}" _existing)
    if(_existing)
        message(STATUS "CUDA toolkit already in the environment: ${_existing}")
        return()
    endif()

    # No version pin. Micromamba reports the installed driver as the __cuda virtual package
    # and conda-forge's cuda-version metapackage is constrained against it, so the solver
    # picks a toolkit this driver can actually run. Pinning here would override that and
    # hand the user a toolkit their driver refuses.
    message(STATUS "NVIDIA driver found and no CUDA toolkit present - fetching one (~600 MB)")
    execute_process(
        COMMAND
            "${CMAKE_COMMAND}" -E env
            "MAMBA_ROOT_PREFIX=${GPU_ROOT_PREFIX}"
            "${GPU_EXECUTABLE}"
            install
            --yes
            --prefix "${GPU_PREFIX}"
            -c conda-forge
            cuda-nvcc cuda-cudart-dev libcublas-dev
        RESULT_VARIABLE _cuda_result
        OUTPUT_VARIABLE _cuda_output
        ERROR_VARIABLE  _cuda_error
    )

    if(NOT _cuda_result EQUAL 0)
        message(STATUS
            "Could not add a CUDA toolkit; building without the GPU path.\n"
            "  This is not an error - NoSpherA2 runs on the CPU. A driver too old for any\n"
            "  packaged toolkit is the usual reason.\n"
            "${_cuda_error}")
        return()
    endif()

    nosphera2_env_nvcc("${GPU_PREFIX}" _installed)
    if(_installed)
        message(STATUS "CUDA toolkit ready: ${_installed}")
    else()
        message(STATUS "CUDA packages installed but no nvcc appeared; building without the GPU path")
    endif()
endfunction()
