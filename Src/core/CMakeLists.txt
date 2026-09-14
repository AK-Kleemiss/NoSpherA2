file(
    GLOB_RECURSE NOSPHERA2_CORE_SOURCES
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx"
)

list(
    FILTER NOSPHERA2_CORE_SOURCES
    EXCLUDE
    REGEX "[/\\\\](vec_nos_math|mat_nos_math)\\.cpp$"
)

#The two GPU host files are for particular builds only. The shim stands in for the HIP
#runtime wherever HIP kernels are compiled in; the dispatcher exists only when both backends
#are, since it defines the global entry points that a single backend's objects define
#themselves.
if(NOT NOSPHERA2_USE_HIP)
    list(FILTER NOSPHERA2_CORE_SOURCES EXCLUDE REGEX "[/\\\\]hip_runtime_shim\\.cpp$")
endif()
if(NOT (NOSPHERA2_USE_CUDA AND NOSPHERA2_USE_HIP))
    list(FILTER NOSPHERA2_CORE_SOURCES EXCLUDE REGEX "[/\\\\]gpu_dispatch\\.cpp$")
endif()

file(
    GLOB_RECURSE NOSPHERA2_CORE_HEADERS
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp"
)

#One source, two languages: nvcc takes it as CUDA, hipcc as HIP. On an MSVC host the HIP
#route compiles them by hand further down (NOSPHERA2_HIP_CLANG, see the top-level HIP block),
#so they are kept out of the target's own sources there.
#
#A fat build (both backends) has to compile each of them twice within one target, and a
#source file has one LANGUAGE. So the HIP compile of a fat build goes through a one-line
#.hip file in the build tree that includes the .cu: the extension makes it HIP to CMake, the
#include makes it the same text, and the depfile keeps a kernel edit rebuilding both halves.
set(NOSPHERA2_CORE_GPU_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/sf_gpu.cu"
    "${CMAKE_CURRENT_SOURCE_DIR}/itensor_gpu.cu"
    "${CMAKE_CURRENT_SOURCE_DIR}/salted_gpu.cu"
    "${CMAKE_CURRENT_SOURCE_DIR}/grid_gpu.cu"
    "${CMAKE_CURRENT_SOURCE_DIR}/blas_gpu.cu"
    "${CMAKE_CURRENT_SOURCE_DIR}/aux_density_gpu.cu")
if(NOSPHERA2_USE_CUDA)
    list(APPEND NOSPHERA2_CORE_SOURCES ${NOSPHERA2_CORE_GPU_SOURCES})
endif()
if(NOSPHERA2_USE_HIP AND NOT NOSPHERA2_HIP_CLANG)
    if(NOSPHERA2_USE_CUDA)
        foreach(_src ${NOSPHERA2_CORE_GPU_SOURCES})
            get_filename_component(_name "${_src}" NAME_WE)
            set(_hip "${CMAKE_CURRENT_BINARY_DIR}/hip/${_name}.hip")
            file(CONFIGURE OUTPUT "${_hip}" CONTENT "#include \"${_src}\"\n" @ONLY)
            list(APPEND NOSPHERA2_CORE_SOURCES "${_hip}")
        endforeach()
    else()
        list(APPEND NOSPHERA2_CORE_SOURCES ${NOSPHERA2_CORE_GPU_SOURCES})
        set_source_files_properties(${NOSPHERA2_CORE_GPU_SOURCES} PROPERTIES LANGUAGE HIP)
    endif()
endif()

#Basis Set Generator:
set(
    NOSPHERA2_BASIS_SET_DIRECTORY
    "${PROJECT_SOURCE_DIR}/BasisSetGenerator/basis_sets"
)

set(
    NOSPHERA2_GENERATED_BASIS_SOURCE
    "${PROJECT_SOURCE_DIR}/Src/basis_data.cpp"
)

# The generated source depends on the CSVs it is generated FROM, not only on the
# converter that reads them.
#
# Without this the rule only ever fired when BasisSetConverter itself was rebuilt,
# so adding or editing a basis set left a stale basis_data.cpp in place - and it
# is .gitignore'd, so nothing in the repository showed a difference either. That
# is how this checkout ended up with 98 basis sets against 117 CSVs: STO-3G was
# missing, and every XCW test failed with "Basis set sto-3g not defined" until the
# file was deleted by hand.
#
# CONFIGURE_DEPENDS re-globs at build time, so a NEW csv triggers a reconfigure;
# naming the files in DEPENDS below covers an EDITED one. Both are needed - the
# glob alone would not notice a change to a file already in the list.
file(
    GLOB NOSPHERA2_BASIS_SET_FILES
    CONFIGURE_DEPENDS
    "${NOSPHERA2_BASIS_SET_DIRECTORY}/*.csv"
)

if(NOT NOSPHERA2_BASIS_SET_FILES)
    # Not fatal: the converter target would already be missing if the submodule
    # were absent. But silence here means a stale basis_data.cpp again, so say so.
    message(WARNING
        "No .csv basis sets found in ${NOSPHERA2_BASIS_SET_DIRECTORY} - "
        "basis_data.cpp will not regenerate when basis sets change. "
        "Is the BasisSetGenerator submodule checked out?"
    )
endif()

add_custom_command(
    OUTPUT
        "${NOSPHERA2_GENERATED_BASIS_SOURCE}"

    COMMAND
        "$<TARGET_FILE:BasisSetConverter>"
        "${NOSPHERA2_BASIS_SET_DIRECTORY}"
        "${PROJECT_SOURCE_DIR}/Src"

    DEPENDS
        BasisSetConverter
        ${NOSPHERA2_BASIS_SET_FILES}

    WORKING_DIRECTORY
        "${PROJECT_SOURCE_DIR}/BasisSetGenerator"

    COMMENT
        "Generating basis_data.cpp"

    VERBATIM
)

#Generator end....

add_library(
    NoSpherA2Core
    STATIC
        ${NOSPHERA2_CORE_SOURCES}
        ${NOSPHERA2_CORE_HEADERS}
        "${NOSPHERA2_GENERATED_BASIS_SOURCE}"
)

nosphera2_enable_optimizations(NoSpherA2Core)

add_library(
    NoSpherA2::Core
    ALIAS
    NoSpherA2Core
)

target_compile_definitions(
    NoSpherA2Core
    PRIVATE
        NOSPHERA2_BASIS_DATA_SEPARATE
    # occ is built without OpenMP and parallelises with TBB; its header templates
    # instantiated in these OpenMP-compiled units would otherwise run Eigen's
    # OpenMP gemm inside every TBB worker. PUBLIC so the app and the tests match.
    PUBLIC
        EIGEN_DONT_PARALLELIZE
)

#Host code sees NOSPHERA2_USE_GPU and nothing else; which backend a kernel object is built
#for is told to that compile alone, by language, so that a fat build can hand the CUDA
#compiles one namespace and the HIP compiles the other (gpu_api.h).
if(NOSPHERA2_USE_CUDA OR NOSPHERA2_USE_HIP)
    target_compile_definitions(NoSpherA2Core PRIVATE NOSPHERA2_USE_GPU)
endif()
if(NOSPHERA2_USE_CUDA AND NOSPHERA2_USE_HIP)
    target_compile_definitions(NoSpherA2Core PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:NOSPHERA2_GPU_BACKEND_NS=nosphera2_cuda>
        $<$<COMPILE_LANGUAGE:HIP>:NOSPHERA2_GPU_BACKEND_NS=nosphera2_hip>)
    message(STATUS "Fat GPU build: CUDA and HIP kernels in one binary")
endif()

if(NOSPHERA2_USE_CUDA)
    target_compile_definitions(NoSpherA2Core PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:NOSPHERA2_USE_CUDA>)
    set_target_properties(NoSpherA2Core PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    #No cuBLAS to link, delay-load or ship: the GEMMs are ours. cudart is static. cuBLAS can
    #still be opened by name at run time under -gpu_cublas, which needs no more from the
    #build than the major version to look for.
    target_compile_definitions(NoSpherA2Core PRIVATE
        NOSPHERA2_CUDART_MAJOR=${CUDAToolkit_VERSION_MAJOR})
    if(NOSPHERA2_CUTLASS_INCLUDE)
        target_include_directories(NoSpherA2Core SYSTEM PRIVATE
            $<$<COMPILE_LANGUAGE:CUDA>:${NOSPHERA2_CUTLASS_INCLUDE}>)
        target_compile_definitions(NoSpherA2Core PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:NOSPHERA2_HAVE_CUTLASS>)
        #CUTLASS pulls in cuda/std, which refuses to compile under MSVC's traditional
        #preprocessor, and calls constexpr host functions from host-device ones.
        target_compile_options(NoSpherA2Core PRIVATE
            $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>
            $<$<COMPILE_LANGUAGE:CUDA>:-diag-suppress=20013,20015>)
        if(MSVC)
            target_compile_options(NoSpherA2Core PRIVATE
                $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/Zc:preprocessor>)
        endif()
    endif()
endif()

if(NOSPHERA2_USE_HIP)
    find_package(hip REQUIRED)
    target_compile_definitions(NoSpherA2Core PRIVATE $<$<COMPILE_LANGUAGE:HIP>:NOSPHERA2_USE_HIP>)
    #The runtime shim is ordinary C++ that needs the HIP headers for the prototypes it
    #defines - the headers only, not hip::host, whose link interface is the runtime library
    #this file exists to avoid.
    set_source_files_properties("${CMAKE_CURRENT_SOURCE_DIR}/hip_runtime_shim.cpp" PROPERTIES
        INCLUDE_DIRECTORIES "${hip_INCLUDE_DIR}"
        COMPILE_DEFINITIONS "__HIP_PLATFORM_AMD__=1")
    if(NOT WIN32)
        target_link_libraries(NoSpherA2Core PRIVATE ${CMAKE_DL_LIBS})
    endif()
    if(NOSPHERA2_HIP_CLANG)
        #MSVC host: each .cu becomes an object through ROCm's clang, with the same runtime
        #library as the rest of the build (-fms-runtime-lib is clang's /MD), and the object
        #joins the archive like one cl.exe produced. The depfile keeps a header edit
        #rebuilding the kernels under Ninja; the launcher is sccache when the CI sets it.
        set(_hip_arch_flags "")
        foreach(_arch ${CMAKE_HIP_ARCHITECTURES})
            list(APPEND _hip_arch_flags "--offload-arch=${_arch}")
        endforeach()
        set(_hip_ns_flags "")
        if(NOSPHERA2_USE_CUDA)
            set(_hip_ns_flags "-DNOSPHERA2_GPU_BACKEND_NS=nosphera2_hip")
        endif()
        foreach(_src ${NOSPHERA2_CORE_GPU_SOURCES})
            get_filename_component(_name "${_src}" NAME_WE)
            set(_obj "${CMAKE_CURRENT_BINARY_DIR}/hip/${_name}.obj")
            add_custom_command(
                OUTPUT "${_obj}"
                COMMAND ${CMAKE_HIP_COMPILER_LAUNCHER} "${NOSPHERA2_HIP_CLANG}"
                    -x hip ${_hip_arch_flags} -std=c++17
                    "$<IF:$<CONFIG:Debug>,-O0,-O3>" "$<$<CONFIG:Debug>:-g>"
                    "-fms-runtime-lib=$<IF:$<CONFIG:Debug>,dll_dbg,dll>"
                    -DNOSPHERA2_USE_GPU -DNOSPHERA2_USE_HIP -D__HIP_PLATFORM_AMD__=1 ${_hip_ns_flags}
                    -I "${CMAKE_CURRENT_SOURCE_DIR}" -I "${PROJECT_SOURCE_DIR}/Src"
                    -isystem "${hip_INCLUDE_DIR}"
                    -MD -MF "${_obj}.d" -c "${_src}" -o "${_obj}"
                DEPFILE "${_obj}.d"
                DEPENDS "${_src}"
                COMMENT "hip ${_name}.cu"
                VERBATIM COMMAND_EXPAND_LISTS)
            target_sources(NoSpherA2Core PRIVATE "${_obj}")
        endforeach()
        #No hip::host: its link interface is amdhip64.lib, and the point of the shim is that
        #nothing imports that DLL.
    endif()
    #On the CMake HIP language route the compiler finds the HIP headers next to itself and
    #defines the platform macro; hip::device would hand -x hip to every C++ source (which a
    #gcc host answers with "language hip not recognized") and hip::host would link the
    #runtime, so neither is used. CMAKE_HIP_RUNTIME_LIBRARY NONE in the top-level file keeps
    #the language itself from linking it.
endif()

target_compile_features(
    NoSpherA2Core
    PUBLIC
        cxx_std_20
)

target_include_directories(
    NoSpherA2Core
    PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
        "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/Src>"
)

target_precompile_headers(
    NoSpherA2Core
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/pch.h"
)

if(APPLE)
    find_library(
        NOSPHERA2_ACCELERATE_FRAMEWORK
        Accelerate
        REQUIRED
    )
    #Apple silicon runs the XCW I tensor on the GPU through Metal Performance Shaders. The
    #kernels are compiled at run time, so only the frameworks are needed here; an Intel Mac
    #has no unified memory and keeps the CPU path.
    if(NOSPHERA2_USE_METAL)
        enable_language(OBJCXX)
        target_sources(NoSpherA2Core PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/itensor_metal.mm")
        set_source_files_properties("${CMAKE_CURRENT_SOURCE_DIR}/itensor_metal.mm"
            PROPERTIES COMPILE_OPTIONS "-fobjc-arc" SKIP_PRECOMPILE_HEADERS ON)
        target_compile_definitions(NoSpherA2Core PRIVATE NOSPHERA2_USE_METAL)
        target_link_libraries(NoSpherA2Core PRIVATE "-framework Metal" "-framework MetalPerformanceShaders" "-framework Foundation")
    endif()
    target_link_libraries(
        NoSpherA2Core
        PRIVATE
            ${NOSPHERA2_ACCELERATE_FRAMEWORK}
        PUBLIC
            OpenMP::OpenMP_CXX
    )

    target_compile_definitions(
        NoSpherA2Core
        PRIVATE
            ACCELERATE_NEW_LAPACK
    )
else()
    target_link_libraries(
        NoSpherA2Core
        PUBLIC
            MKL::MKL
    )
endif()

target_link_libraries(
    NoSpherA2Core
    PUBLIC
        libcint::libcint
        occ_main
        featomic::static
)
