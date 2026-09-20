cmake_minimum_required(VERSION 3.25)

project(NoSpherA2 LANGUAGES C CXX)

message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

option(NOSPHERA2_BUILD_DLL "Build the NoSpherA2 DLL" OFF)
option(NOSPHERA2_BUILD_TESTS "Build NoSpherA2 tests" OFF)
option(NOSPHERA2_DEPENDENCIES_ONLY "Build dependencies without building NoSpherA2" OFF)
option(NOSPHERA2_PROFILE_SYMBOLS "Release code with a PDB, so VTune and the VS profiler resolve functions (scripts/profile_cpu.ps1)" OFF)
option(NOSPHERA2_COVERAGE "gcc/clang: instrument the NoSpherA2 targets for gcov (scripts/coverage.sh), dependencies stay plain" OFF)
option(NOSPHERA2_USE_CUDA "Build the CUDA scattering-factor path" OFF)
option(NOSPHERA2_USE_HIP "Build the HIP (AMD) scattering-factor path" OFF)
#Apple silicon only: the Metal I tensor path wants unified memory, which an Intel Mac has not.
#Detected rather than assumed: the arm64 slice of the universal binary carries the Metal path
#when the SDK the build runs against has the two frameworks, and is a CPU build otherwise, so
#a machine or SDK without them still builds. The frameworks are the only build-time need,
#the kernels being compiled at run time; whether a device is there is decided at run time
#too (itensor_metal.mm gives up when MTLCreateSystemDefaultDevice returns nil).
#The target architecture is CMAKE_OSX_ARCHITECTURES when a preset names one (the x86_64 preset
#on an M-series host is a cross build and gets no Metal), the host's otherwise.
set(_nosphera2_metal_default OFF)
set(_nosphera2_metal_why "not an Apple silicon build")
set(_nosphera2_target_arch "${CMAKE_OSX_ARCHITECTURES}")
if(NOT _nosphera2_target_arch)
    set(_nosphera2_target_arch "${CMAKE_SYSTEM_PROCESSOR}")
endif()
if(APPLE AND _nosphera2_target_arch MATCHES "arm64|aarch64")
    set(_nosphera2_metal_why "the SDK has no Metal or MetalPerformanceShaders framework")
    set(_sdk "${CMAKE_OSX_SYSROOT}")
    if(NOT _sdk OR NOT IS_DIRECTORY "${_sdk}")
        execute_process(COMMAND xcrun --sdk macosx --show-sdk-path
            OUTPUT_VARIABLE _sdk ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
    endif()
    if(IS_DIRECTORY "${_sdk}/System/Library/Frameworks/Metal.framework"
       AND IS_DIRECTORY "${_sdk}/System/Library/Frameworks/MetalPerformanceShaders.framework")
        set(_nosphera2_metal_default ON)
        set(_nosphera2_metal_why "Apple silicon, frameworks in ${_sdk}")
    endif()
endif()
option(NOSPHERA2_USE_METAL "Build the Metal XCW I tensor path (Apple silicon)" ${_nosphera2_metal_default})
if(APPLE)
    if(NOSPHERA2_USE_METAL)
        message(STATUS "Metal I tensor path: ON (${_nosphera2_metal_why})")
    elseif(_nosphera2_metal_default)
        message(STATUS "Metal I tensor path: OFF (switched off, NOSPHERA2_USE_METAL=OFF)")
    else()
        message(STATUS "Metal I tensor path: OFF (${_nosphera2_metal_why}; -DNOSPHERA2_USE_METAL=ON to force)")
    endif()
endif()
option(NOSPHERA2_GPU_AUTO "Pick the GPU backend from the hardware present" ON)
#For a binary that will be copied to other machines rather than run where it was built.
option(NOSPHERA2_CUDA_PORTABLE "Build CUDA code for every major architecture, not just this card" OFF)
#The AMD counterpart. HIP has no PTX, so a binary runs only on the ISAs it was compiled for
#and a portable one has to name every card it might meet; the list is in the HIP block below.
option(NOSPHERA2_HIP_PORTABLE "Build HIP code for every supported AMD architecture, not just this card" OFF)
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
    #A Visual Studio developer prompt can retain an older toolkit first on PATH after a
    #new CUDA Toolkit is installed. Prefer the newest system nvcc in that case: otherwise
    #check_language() rejects an old toolkit before it ever considers the new one.
    if(WIN32 AND NOT CMAKE_CUDA_COMPILER)
        file(GLOB _system_nvcc "$ENV{ProgramFiles}/NVIDIA GPU Computing Toolkit/CUDA/v*/bin/nvcc.exe")
        if(_system_nvcc)
            list(SORT _system_nvcc COMPARE NATURAL)
            list(GET _system_nvcc -1 _system_nvcc)
            set(CMAKE_CUDA_COMPILER "${_system_nvcc}" CACHE FILEPATH "Newest system CUDA compiler" FORCE)
            message(STATUS "Using the newest system CUDA compiler: ${_system_nvcc}")
        endif()
    endif()
endif()

#Both backends at once is a fat binary: each kernel source is compiled twice, once per
#backend, into its own namespace (Src/core/gpu_api.h), and Src/core/gpu_dispatch.cpp forwards
#the global entry points to whichever backend has a device. Neither runtime is linked - CUDA's
#is static, HIP's is opened by name (Src/core/hip_runtime_shim.cpp) - so the one binary runs
#on an NVIDIA machine, an AMD machine and a machine with no GPU at all. That is what the CI
#ships (.github/workflows/c-cpp_all.yml).

#Whether this file may choose CMAKE_<LANG>_ARCHITECTURES. A list given on the command line is
#the user's and is left alone. The trap is that enable_language() creates the cache entry
#itself when nobody has (sm_52 for nvcc, gfx906 or whatever rocm_agent_enumerator reports
#for HIP), so testing DEFINED after the language is enabled always says "the user set it"
#and the portable list never applies - the Linux HIP build came out as gfx906 only. Asked
#before the language is enabled, and a list this file wrote last time (remembered in
#NOSPHERA2_<LANG>_ARCH_MANAGED) counts as not given, so a reconfigure with the portable
#option toggled takes effect. A tree configured before this marker existed still carries
#CMake's default without the marker and looks user-set: -UCMAKE_<LANG>_ARCHITECTURES once.
macro(nosphera2_arch_managed lang)
    set(_nosphera2_${lang}_arch_managed OFF)
    if(NOT DEFINED CMAKE_${lang}_ARCHITECTURES
       OR (DEFINED NOSPHERA2_${lang}_ARCH_MANAGED
           AND "${CMAKE_${lang}_ARCHITECTURES}" STREQUAL "${NOSPHERA2_${lang}_ARCH_MANAGED}"))
        set(_nosphera2_${lang}_arch_managed ON)
    endif()
endmacro()
macro(nosphera2_set_arch lang)
    set(CMAKE_${lang}_ARCHITECTURES "${ARGN}" CACHE STRING "${lang} architectures" FORCE)
    set(NOSPHERA2_${lang}_ARCH_MANAGED "${ARGN}" CACHE INTERNAL "The list this project chose")
endmacro()

#Nothing delay-loads a GPU library any more: cuBLAS was the last CUDA one, and the HIP
#runtime is opened by name (Src/core/hip_runtime_shim.cpp) rather than delay-loaded. Deleting
#the code that set these was not enough to stop them being used: they were cached, so a build
#directory configured before that change kept handing the flag and delayimp to every link.
#Cleared rather than left to rot, the same way an unusable cached compiler is above.
unset(NOSPHERA2_CUDA_DELAYLOAD CACHE)
unset(NOSPHERA2_HIP_DELAYLOAD CACHE)

if(NOSPHERA2_USE_CUDA)
    #nvcc takes its host compiler from PATH unless told, and on a machine whose PATH carries
    #a different gcc than the one CMake picked for C++ the two halves of one binary would be
    #built by different compilers. Named before the language is enabled, since that is when
    #it is first used. Windows is left alone: nvcc finds cl.exe through the developer prompt.
    if(NOT WIN32 AND NOT DEFINED CMAKE_CUDA_HOST_COMPILER AND NOT DEFINED ENV{CUDAHOSTCXX})
        set(CMAKE_CUDA_HOST_COMPILER "${CMAKE_CXX_COMPILER}" CACHE FILEPATH "Host compiler for nvcc")
    endif()
    include(CheckLanguage)
    nosphera2_arch_managed(CUDA)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
        set(CMAKE_CUDA_STANDARD 17)
        set(CMAKE_CUDA_STANDARD_REQUIRED ON)
        #A distributed executable must not import cudart64_*.dll. The static CUDA runtime
        #opens the display-driver API only when a GPU path is used, so the ordinary CPU
        #fallback still starts on a machine with no NVIDIA software installed.
        set(CMAKE_CUDA_RUNTIME_LIBRARY Static)
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
        #needs anyway. Blackwell arrived with 12.8, so a 12.8 or 12.9 toolkit is the one
        #that covers everything from a V100 to an RTX 50 in a single binary, which is why
        #the CI release build pins 12.9 (scripts/GpuToolkit.cmake).
        if(_nosphera2_CUDA_arch_managed)
            if(NOSPHERA2_CUDA_PORTABLE)
                if(CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL 13)
                    set(_arch_list 75 80 86 89 90 100 120)
                elseif(CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL 12.8)
                    set(_arch_list 70 75 80 86 89 90 100 120)
                else()
                    set(_arch_list 70 75 80 86 89 90)
                endif()
            else()
                set(_arch_list native)
            endif()
            nosphera2_set_arch(CUDA ${_arch_list})
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
    #Two ways to compile the kernels. On Linux they are a HIP language for CMake and clang
    #handles the whole program. On Windows the program is built with MSVC and CMake refuses to
    #mix that with a clang HIP compiler ("mixes Clang and MSVC", Platform/Windows-Clang.cmake,
    #still there in 4.3); building all of NoSpherA2 with ROCm's clang is a different compiler
    #for the whole program, so instead the kernels are compiled the way the Visual Studio
    #projects do it (NoSpherA2_gpu.props): clang++ -x hip on each .cu into an MSVC-ABI object
    #that links like any other, in Src/core. This block finds that clang and settles the
    #architecture list for both routes. The compiler is taken from -DCMAKE_HIP_COMPILER, else
    #from HIP_PATH (the HIP SDK keeps it in bin, the ROCm pip wheels in lib/llvm/bin), else
    #from the CMAKE_PREFIX_PATH entries - never from PATH, where Visual Studio's own clang++
    #would be found first and cannot compile HIP.
    set(_nosphera2_hip_version "")
    if(MSVC)
        if(CMAKE_HIP_COMPILER)
            set(NOSPHERA2_HIP_CLANG "${CMAKE_HIP_COMPILER}" CACHE FILEPATH "ROCm clang++ that compiles the HIP kernels")
        else()
            set(_hip_clang_hints "")
            if(DEFINED ENV{HIP_PATH})
                file(TO_CMAKE_PATH "$ENV{HIP_PATH}" _hip_path)
                list(APPEND _hip_clang_hints "${_hip_path}/bin" "${_hip_path}/lib/llvm/bin")
            endif()
            foreach(_p ${CMAKE_PREFIX_PATH})
                list(APPEND _hip_clang_hints "${_p}/bin" "${_p}/lib/llvm/bin")
            endforeach()
            find_program(NOSPHERA2_HIP_CLANG clang++ HINTS ${_hip_clang_hints} NO_DEFAULT_PATH
                DOC "ROCm clang++ that compiles the HIP kernels")
        endif()
        if(NOSPHERA2_HIP_CLANG)
            execute_process(COMMAND "${NOSPHERA2_HIP_CLANG}" --version
                OUTPUT_VARIABLE _hip_clang_banner ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
            if(_hip_clang_banner MATCHES "clang version ([0-9]+\\.[0-9]+\\.[0-9]+)")
                set(_nosphera2_hip_version "${CMAKE_MATCH_1}")
            endif()
            message(STATUS "HIP kernels compiled by ${NOSPHERA2_HIP_CLANG} (clang ${_nosphera2_hip_version})")
        endif()
        #The only HIP platform on Windows; hip-config.cmake would otherwise ask hipconfig,
        #which the pip wheels do not ship.
        if(NOT DEFINED HIP_PLATFORM AND NOT DEFINED ENV{HIP_PLATFORM})
            set(HIP_PLATFORM "amd")
        endif()
        set(_nosphera2_hip_ready ${NOSPHERA2_HIP_CLANG})
    else()
        include(CheckLanguage)
        nosphera2_arch_managed(HIP)
        check_language(HIP)
        if(CMAKE_HIP_COMPILER)
            enable_language(HIP)
            set(CMAKE_HIP_STANDARD 17)
            set(CMAKE_HIP_STANDARD_REQUIRED ON)
            set(_nosphera2_hip_version "${CMAKE_HIP_COMPILER_VERSION}")
            #The HIP language would otherwise link libamdhip64 into every target that
            #carries a HIP object, directly or through a static library (hip-lang::device,
            #CMakeHIPInformation.cmake), and a binary that needs libamdhip64.so to start is
            #a binary that does not start without ROCm. NONE names no runtime library at
            #all; Src/core/hip_runtime_shim.cpp supplies the entry points and opens the real
            #library by name when there is one.
            set(CMAKE_HIP_RUNTIME_LIBRARY NONE)
        endif()
        set(_nosphera2_hip_ready ${CMAKE_HIP_COMPILER})
    endif()
    if(_nosphera2_hip_ready)
        #Unlike CUDA there is no intermediate code the driver could finish at run time: a
        #kernel built for gfx1030 does not run on a gfx1031, and the generic targets
        #(gfx10-3-generic and so on) demand a runtime newer than most machines carry. So the
        #portable list is every ISA ROCm still supports, one by one - GCN Vega, the three CDNA
        #generations and RDNA 1 to 4. The two newest generations are only known to the clang
        #of ROCm 6.3 and later, so the list is trimmed on an older toolkit rather than
        #failing the build. A card missing from the list is what sf_gpu_available() reports.
        if(MSVC)
            #No enable_language on this path, so nothing has created the entry yet.
            nosphera2_arch_managed(HIP)
        endif()
        if(_nosphera2_HIP_arch_managed)
            if(NOSPHERA2_HIP_PORTABLE)
                set(_arch_list
                    gfx900 gfx906 gfx908 gfx90a gfx942
                    gfx1010 gfx1012 gfx1030 gfx1031 gfx1032 gfx1034 gfx1035 gfx1036
                    gfx1100 gfx1101 gfx1102 gfx1103)
                if(_nosphera2_hip_version VERSION_GREATER_EQUAL 19)
                    list(APPEND _arch_list gfx1150 gfx1151 gfx1200 gfx1201)
                endif()
                if(_nosphera2_hip_version VERSION_GREATER_EQUAL 20)
                    list(APPEND _arch_list gfx950)
                endif()
            else()
                set(_arch_list gfx90a gfx942 gfx1100)
            endif()
            nosphera2_set_arch(HIP ${_arch_list})
        endif()
        #Nothing links the HIP runtime, on either route. ROCm's lib/bin directory is not on
        #the path of a machine that has the card but not the SDK, and a binary whose loader
        #wants libamdhip64.so / amdhip64_7.dll is rejected before main() (measured on the
        #AMD node as STATUS_DLL_NOT_FOUND). Src/core/hip_runtime_shim.cpp defines the entry
        #points the kernel objects import and opens the library by name on the first call,
        #answering "no device" when it is absent. The name carries the HIP runtime's major
        #version, which is neither the ROCm release nor the clang version (ROCm 6.x ships
        #amdhip64_6.dll, ROCm 7 to 10 amdhip64_7.dll); the shim takes it from HIP_VERSION_MAJOR
        #in the headers it was compiled against. hipBLAS is not needed at all, which is what
        #lets every kernel build on a ROCm that packages the compiler without it.
        if(MSVC)
            #hip-config runs this to name the compiler-rt builtins clang's host code may
            #need; left to itself it would ask cl.exe and get nothing.
            set(HIP_CXX_COMPILER "${NOSPHERA2_HIP_CLANG}")
            find_package(hip REQUIRED)
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
