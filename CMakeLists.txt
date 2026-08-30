cmake_minimum_required(VERSION 3.25)

project(NoSpherA2 LANGUAGES C CXX)

message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

option(NOSPHERA2_BUILD_DLL "Build the NoSpherA2 DLL" OFF)
option(NOSPHERA2_BUILD_TESTS "Build NoSpherA2 tests" OFF)
option(NOSPHERA2_DEPENDENCIES_ONLY "Build dependencies without building NoSpherA2" OFF)
option(NOSPHERA2_USE_CUDA "Build the CUDA scattering-factor path" OFF)
option(NOSPHERA2_USE_HIP "Build the HIP (AMD) scattering-factor path" OFF)
option(NOSPHERA2_GPU_AUTO "Pick the GPU backend from the hardware present" ON)
#For a binary that will be copied to other machines rather than run where it was built.
option(NOSPHERA2_CUDA_PORTABLE "Build CUDA code for every major architecture, not just this card" OFF)
#CUTLASS runs the I tensor GEMM at cuBLAS speed against the built-in kernel's ~57%. Off
#turns that call site back over to the built-in one, which is also what a HIP build uses.
option(NOSPHERA2_USE_CUTLASS "Use CUTLASS for the I tensor GEMM (CUDA only)" ON)

#Detection may switch a backend on, never off, and only runs while neither is already on.
#Both halves are load-bearing. It has to re-decide rather than skip once the variables exist,
#because option() caches them on the first configure and from the second onwards "the user
#asked for this" and "nobody has ever asked" are the same state - a build directory created
#before the toolkit existed would otherwise stay CPU-only however often it was reconfigured.
#But overwriting unconditionally is worse: it overrides an explicit -DNOSPHERA2_USE_CUDA=ON
#and so breaks building for a GPU on a machine that has none, which is the ordinary cluster
#arrangement of a GPU-less head node compiling for the compute nodes. Deferring to a backend
#that is already on covers both, since wanting to build for an absent GPU is exactly what
#detection cannot observe. Set NOSPHERA2_GPU_AUTO=OFF to pin the two flags by hand.
if(NOSPHERA2_GPU_AUTO AND NOT NOSPHERA2_USE_CUDA AND NOT NOSPHERA2_USE_HIP)
    include("${CMAKE_CURRENT_LIST_DIR}/scripts/GpuToolkit.cmake")
    nosphera2_detect_gpu(NOSPHERA2_GPU_VENDOR)
    if(NOSPHERA2_GPU_VENDOR STREQUAL "NVIDIA")
        set(NOSPHERA2_USE_CUDA ON CACHE BOOL "" FORCE)
    elseif(NOSPHERA2_GPU_VENDOR STREQUAL "AMD")
        set(NOSPHERA2_USE_HIP ON CACHE BOOL "" FORCE)
    endif()
    message(STATUS "GPU auto-detect: ${NOSPHERA2_GPU_VENDOR} (-DNOSPHERA2_GPU_AUTO=OFF to choose by hand)")
    #A missing compiler after auto-enabling is a status line below, not the warning someone
    #who named the backend deserves.
    set(_nosphera2_gpu_auto_enabled TRUE)
endif()

#The toolkit the bootstrap fetched, in preference to whatever is first on PATH. Without this
#the build silently follows PATH: on a machine carrying both, the binary can end up importing
#a cuBLAS that the machine it runs on does not have.
#
#Tested for usability rather than for being set. check_language caches CMAKE_CUDA_COMPILER as
#NOTFOUND when it fails, so a directory configured before the toolkit arrived holds a defined
#but useless value, and check_language will not look again. Gating on DEFINED here reproduced
#exactly the freeze this file already has one fix for: bootstrap reports a toolkit, configure
#reports no compiler, and no amount of reconfiguring changes its mind.
if(NOSPHERA2_USE_CUDA)
    #Drop a cached compiler that no longer resolves, so check_language searches again. This
    #covers a toolkit installed system-wide after the failed configure as well as a fetched one.
    if(DEFINED CMAKE_CUDA_COMPILER AND (NOT CMAKE_CUDA_COMPILER OR NOT EXISTS "${CMAKE_CUDA_COMPILER}"))
        unset(CMAKE_CUDA_COMPILER CACHE)
    endif()
    if(MICROMAMBA_ENV_PREFIX AND NOT CMAKE_CUDA_COMPILER)
        foreach(_c "${MICROMAMBA_ENV_PREFIX}/Library/bin/nvcc.exe" "${MICROMAMBA_ENV_PREFIX}/bin/nvcc")
            if(EXISTS "${_c}")
                set(CMAKE_CUDA_COMPILER "${_c}" CACHE FILEPATH "CUDA compiler from the micromamba environment" FORCE)
                message(STATUS "Using the bootstrapped CUDA compiler: ${_c}")
                break()
            endif()
        endforeach()
    endif()
endif()

if(NOSPHERA2_USE_CUDA AND NOSPHERA2_USE_HIP)
    message(FATAL_ERROR "NOSPHERA2_USE_CUDA and NOSPHERA2_USE_HIP are mutually exclusive")
endif()

#Nothing delay-loads a CUDA library any more, cuBLAS having been the only one. Deleting the
#code that set this was not enough to stop it being used: it was cached, so a build directory
#configured before that change kept handing the flag and delayimp to every link. Cleared
#rather than left to rot, the same way an unusable cached compiler is above.
unset(NOSPHERA2_CUDA_DELAYLOAD CACHE)

if(NOSPHERA2_USE_CUDA)
    include(CheckLanguage)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
        set(CMAKE_CUDA_STANDARD 17)
        set(CMAKE_CUDA_STANDARD_REQUIRED ON)
        #native builds for the card doing the build and nothing else, which is right for a
        #local build and wrong for one that gets copied.
        #
        #The portable list names minor architectures, not major ones. all-major would be the
        #obvious choice and is a trap: a cubin runs on any later minor of the same major, so
        #sm_80 does execute on an sm_86 or sm_89 card - but as sm_80 code. Measured on the I
        #tensor, an sm_80 cubin on an sm_89 card gives up 17% against one built for sm_89.
        #The consumer parts are exactly the ones that would inherit that, so they are named.
        #
        #Chosen from the compiler version rather than by filtering
        #CMAKE_CUDA_ARCHITECTURES_ALL. That list suffixes its entries with -real, so a plain
        #IN_LIST against bare numbers silently matches almost nothing, and in a directory
        #that has seen more than one toolkit it goes stale as well. The version is neither.
        #
        #CUDA 13 dropped Volta, so 70 appears only on a 12.x toolkit - which is what a V100
        #needs anyway - and Blackwell only on 13.
        if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
            if(NOSPHERA2_CUDA_PORTABLE)
                if(CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL 13)
                    set(CMAKE_CUDA_ARCHITECTURES 75 80 86 89 90 100 120)
                else()
                    set(CMAKE_CUDA_ARCHITECTURES 70 75 80 86 89 90)
                endif()
            else()
                set(CMAKE_CUDA_ARCHITECTURES native)
            endif()
        endif()
        #A pinned list compiles to those architectures and nothing else, so the binary has
        #no code at all for a card outside it and every GPU path quietly declines - which
        #looks exactly like having no GPU. Carrying PTX for the newest entry lets the driver
        #compile for anything newer on first run. It cannot help an older card: PTX is
        #compiled forward, never back, so this is a backstop and not a substitute for
        #NOSPHERA2_CUDA_PORTABLE. "native" is left alone, being the one case that cannot
        #arise.
        if(NOT CMAKE_CUDA_ARCHITECTURES MATCHES "native|all|all-major")
            set(_arch_sorted ${CMAKE_CUDA_ARCHITECTURES})
            list(SORT _arch_sorted COMPARE NATURAL)
            list(GET _arch_sorted -1 _arch_newest)
            string(REGEX MATCH "^[0-9]+" _arch_newest "${_arch_newest}")
            if(_arch_newest AND NOT "${_arch_newest}-virtual" IN_LIST CMAKE_CUDA_ARCHITECTURES)
                list(APPEND CMAKE_CUDA_ARCHITECTURES "${_arch_newest}-virtual")
            endif()
        endif()
        find_package(CUDAToolkit REQUIRED)
        #Nothing to delay-load any more. cuBLAS was the only dynamic CUDA library linked,
        #and the GEMMs that used it are ours now, so the binary imports no DLL that a
        #machine without CUDA could fail to provide.
        message(STATUS "CUDA enabled, arch ${CMAKE_CUDA_ARCHITECTURES}")

        #CUTLASS supplies the I tensor's GEMM, which is the only device GEMM on a hot path.
        #Headers only: SOURCE_SUBDIR names a directory that does not exist so that
        #MakeAvailable populates without configuring CUTLASS's own build, which would bring
        #its tests, examples and profiler with it. Nothing is linked and nothing ships - it
        #compiles into this binary, so it does not reintroduce a runtime dependency.
        if(NOSPHERA2_USE_CUTLASS)
            include(FetchContent)
            FetchContent_Declare(cutlass
                URL https://github.com/NVIDIA/cutlass/archive/refs/tags/v4.7.1.tar.gz
                SOURCE_SUBDIR headers-only-no-build)
            FetchContent_MakeAvailable(cutlass)
            if(EXISTS "${cutlass_SOURCE_DIR}/include/cutlass/cutlass.h")
                set(NOSPHERA2_CUTLASS_INCLUDE "${cutlass_SOURCE_DIR}/include"
                    CACHE INTERNAL "CUTLASS headers")
                message(STATUS "CUTLASS enabled from ${cutlass_SOURCE_DIR}")
            else()
                message(WARNING "CUTLASS was fetched but its headers are missing; "
                                "the I tensor falls back to the built-in GEMM")
            endif()
        endif()
    elseif(_nosphera2_gpu_auto_enabled)
        message(STATUS "NVIDIA GPU present but no CUDA compiler found; building for the CPU")
        set(NOSPHERA2_USE_CUDA OFF)
    else()
        message(WARNING "NOSPHERA2_USE_CUDA is ON but no CUDA compiler was found; building without it")
        set(NOSPHERA2_USE_CUDA OFF)
    endif()
endif()

if(NOSPHERA2_USE_HIP)
    include(CheckLanguage)
    check_language(HIP)
    if(CMAKE_HIP_COMPILER)
        enable_language(HIP)
        set(CMAKE_HIP_STANDARD 17)
        set(CMAKE_HIP_STANDARD_REQUIRED ON)
        if(NOT DEFINED CMAKE_HIP_ARCHITECTURES)
            set(CMAKE_HIP_ARCHITECTURES gfx90a;gfx942;gfx1100)
        endif()
        #ROCm's bin directory is not on PATH by default, so a HIP build fails to start with
        #STATUS_DLL_NOT_FOUND on a machine that has the GPU but not the SDK. Measured on the
        #AMD node: without this the process is rejected before main(). hipBLAS is no longer
        #among these, which is what lets every kernel build on a conda-forge ROCm - it
        #packages hipcc without hipBLAS, and that used to strand the I tensor on the CPU.
        if(MSVC)
            string(REGEX MATCH "^[0-9]+" NOSPHERA2_HIP_MAJOR "${CMAKE_HIP_COMPILER_VERSION}")
            if(NOT NOSPHERA2_HIP_MAJOR)
                set(NOSPHERA2_HIP_MAJOR 6)
            endif()
            set(NOSPHERA2_HIP_DELAYLOAD
                "/DELAYLOAD:amdhip64_${NOSPHERA2_HIP_MAJOR}.dll"
                CACHE INTERNAL "ROCm DLLs to delay-load")
        endif()
        message(STATUS "HIP enabled, arch ${CMAKE_HIP_ARCHITECTURES}")
    elseif(_nosphera2_gpu_auto_enabled)
        message(STATUS "AMD GPU present but no HIP compiler found; building for the CPU")
        set(NOSPHERA2_USE_HIP OFF)
    else()
        message(WARNING "NOSPHERA2_USE_HIP is ON but no HIP compiler was found; building without it")
        set(NOSPHERA2_USE_HIP OFF)
    endif()
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NO_TRIVIAL_INIT)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        add_compile_options(-ftrivial-auto-var-init=pattern)
    else()
        message(WARNING
            "NO_TRIVIAL_INIT is enabled, but "
            "${CMAKE_CXX_COMPILER_ID} does not support "
            "-ftrivial-auto-var-init=pattern"
        )
    endif()
endif()

if(WIN32)
    set(CMAKE_TRY_COMPILE_CONFIGURATION Release)
endif()

list(PREPEND CMAKE_MODULE_PATH
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake"
)

include(NoSpherA2Avx)
include(NoSpherA2Sanitizers)
include(ValidateMambaEnvironemt)
include(CopyRuntimeLibraries)
include(NoSpherA2Optimizations)
validate_micromamba_environment()


# ---------------------------------------------------------------------------
# libcint
# ---------------------------------------------------------------------------

set(WITH_FORTRAN OFF)
set(WITH_CINT2_INTERFACE OFF)
set(ENABLE_STATIC ON)
set(PYPZPX ON)
set(WITH_RANGE_COULOMB ON)
set(WITH_PYTHON_BINDINGS OFF)
set(BUILD_SHARED_LIBS OFF)

add_subdirectory(
    libcint
    "${CMAKE_BINARY_DIR}/_deps/libcint"
    EXCLUDE_FROM_ALL
)

# Patch libcint GCC ivdep pragmas for MSVC without modifying the submodule.
if(MSVC)
    set(_libcint_ivdep_pragma
        "#if defined(_MSC_VER)\n#pragma loop(ivdep)\n#elif defined(__GNUC__) || defined(__clang__)\n#pragma GCC ivdep\n#endif"
    )

    get_target_property(_cint_sources cint SOURCES)

    foreach(_libcint_ivdep_source IN ITEMS
        "src/cart2sph.c"
        "src/g1e_grids.c"
        "src/fblas.c"
    )
        set(_libcint_ivdep_in
            "${CMAKE_CURRENT_SOURCE_DIR}/libcint/${_libcint_ivdep_source}"
        )
        set(_libcint_ivdep_out
            "${CMAKE_CURRENT_BINARY_DIR}/patched/libcint/${_libcint_ivdep_source}"
        )

        file(READ "${_libcint_ivdep_in}" _libcint_ivdep_text)
        string(REPLACE
            "#pragma GCC ivdep"
            "${_libcint_ivdep_pragma}"
            _libcint_ivdep_text
            "${_libcint_ivdep_text}"
        )

        get_filename_component(_libcint_ivdep_out_dir
            "${_libcint_ivdep_out}"
            DIRECTORY
        )
        file(MAKE_DIRECTORY "${_libcint_ivdep_out_dir}")
        file(WRITE "${_libcint_ivdep_out}" "${_libcint_ivdep_text}")

        list(REMOVE_ITEM _cint_sources "${_libcint_ivdep_source}")
        list(APPEND _cint_sources "${_libcint_ivdep_out}")
    endforeach()

    set_target_properties(cint PROPERTIES SOURCES "${_cint_sources}")
endif()

if(NOT TARGET cint)
    message(FATAL_ERROR
        "The bundled libcint project did not define the target 'cint'"
    )
endif()

# libcint's current CMake target does not expose all source-tree public
# headers when consumed directly with add_subdirectory().
target_include_directories(
    cint
    PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/libcint/include>"
)

if(NOT TARGET libcint::libcint)
    add_library(libcint::libcint ALIAS cint)
endif()

# ---------------------------------------------------------------------------
# OCC
# ---------------------------------------------------------------------------

set(USE_SYSTEM_LIBCINT ON)
set(BUILD_FPIC ON)
set(WITH_PYTHON_BINDINGS OFF)
set(BUILD_SHARED_LIBS OFF)

add_subdirectory(
    occ
    "${CMAKE_BINARY_DIR}/_deps/occ"
    EXCLUDE_FROM_ALL
)

# ---------------------------------------------------------------------------
# mdspan
# ---------------------------------------------------------------------------

add_subdirectory(
    mdspan
    "${CMAKE_BINARY_DIR}/_deps/mdspan"
    EXCLUDE_FROM_ALL
)

# ---------------------------------------------------------------------------
# featomic
# ---------------------------------------------------------------------------

set(FEATOMIC_FETCH_METATENSOR ON CACHE BOOL
    "Fetch the matching metatensor-core source" FORCE
)
set(FEATOMIC_INSTALL_BOTH_STATIC_SHARED OFF CACHE BOOL
    "Build only the selected featomic library type" FORCE
)
set(FEATOMIC_USE_STATIC_METATENSOR ON CACHE BOOL
    "Link featomic against static metatensor-core" FORCE
)
set(BUILD_SHARED_LIBS OFF)

add_subdirectory(
    featomic/featomic
    "${CMAKE_BINARY_DIR}/_deps/featomic"
    EXCLUDE_FROM_ALL
)

# ---------------------------------------------------------------------------
# Intel MKL
# ---------------------------------------------------------------------------

set(MKL_LINK "static")
set(MKL_INTERFACE lp64)
set(MKL_THREADING intel_thread)
set(MKL_MPI intelmpi)
if(WIN32)
    set(
        MKL_DIR
        "${MICROMAMBA_ENV_PREFIX}/Library/lib/cmake/mkl"
        CACHE PATH
        "MKL CMake package directory"
    )
else()
    set(
        MKL_DIR
        "${MICROMAMBA_ENV_PREFIX}/lib/cmake/mkl"
        CACHE PATH
        "MKL CMake package directory"
    )
endif()

if(APPLE)
    find_library(
        NOSPHERA2_ACCELERATE_FRAMEWORK
        Accelerate
        REQUIRED
    )
    find_package(OpenMP REQUIRED COMPONENTS CXX)
else()
    find_package(MKL CONFIG REQUIRED)
endif()


# ---------------------------------------------------------------------------
# Basis Set Generation
# ---------------------------------------------------------------------------
add_subdirectory(
    BasisSetGenerator
    "${CMAKE_BINARY_DIR}/tools/BasisSetConverter"
)

# ---------------------------------------------------------------------------
# Function to include GoogleTest
# ---------------------------------------------------------------------------

function(enable_googletests)
    include(FetchContent)
    FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
    )
    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
endfunction()

# ---------------------------------------------------------------------------
# NoSpherA2 targets
# ---------------------------------------------------------------------------
if(NOSPHERA2_DEPENDENCIES_ONLY)
    enable_googletests()
    include(cmake/InstallDependenciesOnly.cmake)
    return()
endif()

add_subdirectory(Src/core)
add_subdirectory(app)

#Small --install target to copy the contents of bin. USAGE: cmake --install build/release-linux --prefix "$PWD/install-release"
include(GNUInstallDirs)
install(
    DIRECTORY
        "${CMAKE_BINARY_DIR}/bin/"
    DESTINATION
        "${CMAKE_INSTALL_BINDIR}"
    USE_SOURCE_PERMISSIONS
    COMPONENT Runtime
)


if(NOSPHERA2_BUILD_DLL)
    add_subdirectory(Src/dll)
endif()

include(CTest)

if(NOSPHERA2_BUILD_TESTS)
    add_subdirectory(tests/src)
endif()

if(APPLE)
    add_custom_target(
        NoSpherA2Universal

        COMMAND
            "${CMAKE_COMMAND}"
            "-DARM64_DIRECTORY=${PROJECT_SOURCE_DIR}/build/release-macos-arm64/bin"
            "-DX86_64_DIRECTORY=${PROJECT_SOURCE_DIR}/build/release-macos-x86_64/bin"
            "-DOUTPUT_DIRECTORY=${PROJECT_SOURCE_DIR}/build/release-macos-universal/bin"
            -P
            "${PROJECT_SOURCE_DIR}/cmake/CreateMacOSUniversal.cmake"

        COMMENT
            "Creating universal macOS binaries"

        VERBATIM
    )
endif()
