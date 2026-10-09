# AVX detection / compile options for NoSpherA2 + OCC.
include_guard(GLOBAL)

# Vector instruction policy. OCC and NoSpherA2 MUST be built with the same
# setting: Eigen object layouts differ between AVX and non-AVX builds.
#   NOS_AVX=ON   build with AVX
#   NOS_AVX=OFF  SSE4.2 only — the compatibility build, used by CI and for binaries handed out
#   unset        honor the NOS_AVX environment variable, else detect the host (AVX, then AVX2 + FMA)
#   NOS_AVX2=ON/OFF  force AVX2 + FMA on (implies AVX) or keep a detecting build at AVX
# The ISA is a build-time choice, never a run-time dispatch: a native build uses what its host has.
if(NOT DEFINED NOS_AVX)
    if(DEFINED ENV{NOS_AVX} AND NOT "$ENV{NOS_AVX}" STREQUAL "")
        set(NOS_AVX "$ENV{NOS_AVX}")
    else()
        set(NOS_AVX "AUTO")
    endif()
endif()

set(NOS_USE_AVX OFF)
set(NOS_USE_AVX2 OFF)
if(NOS_AVX STREQUAL "AUTO")
    if(NOT CMAKE_CROSSCOMPILING AND CMAKE_SYSTEM_PROCESSOR MATCHES "(x86_64|AMD64|amd64)")
        try_run(NOS_AVX_RUN_RESULT NOS_AVX_COMPILE_RESULT
            ${CMAKE_BINARY_DIR}/nos_avx_check
            ${CMAKE_CURRENT_SOURCE_DIR}/cmake/detect_avx.c)
        if(NOS_AVX_COMPILE_RESULT AND NOS_AVX_RUN_RESULT EQUAL 0)
            set(NOS_USE_AVX ON)
            if(NOT DEFINED NOS_AVX2)
                try_run(NOS_AVX2_RUN_RESULT NOS_AVX2_COMPILE_RESULT
                    ${CMAKE_BINARY_DIR}/nos_avx2_check
                    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/detect_avx.c
                    COMPILE_DEFINITIONS -DNOS_DETECT_AVX2)
                if(NOS_AVX2_COMPILE_RESULT AND NOS_AVX2_RUN_RESULT EQUAL 0)
                    set(NOS_USE_AVX2 ON)
                endif()
            endif()
        endif()
    endif()
elseif(NOS_AVX)
    set(NOS_USE_AVX ON)
endif()
if(NOS_AVX2)
    set(NOS_USE_AVX ON)
    set(NOS_USE_AVX2 ON)
endif()
message(STATUS "NOS_AVX=${NOS_AVX} -> building with AVX: ${NOS_USE_AVX}, AVX2+FMA: ${NOS_USE_AVX2}")

#C and C++ only: nvcc refuses -m flags, and the GPU sources get their host flags from it
if(LINUX)
    add_compile_options("$<$<COMPILE_LANGUAGE:C,CXX>:-msse2;-msse3;-msse4.1;-msse4.2>")
    if(NOS_USE_AVX)
        add_compile_options($<$<COMPILE_LANGUAGE:C,CXX>:-mavx>)
    endif()
endif()

#AVX2 is -mavx2 -mfma, i.e. 256-bit packed fp64 with a fused multiply-add. On Linux it has to be a
#LINK option as well as a compile option, and that is not a detail: the build is -flto=auto, so the
#machine code is generated at link time and per-file compile flags produced a byte-identical binary
#when this was tried per translation unit. Only where the host has it - an AVX2 instruction on a
#pre-Haswell host is SIGILL, not a slow path. i7-7700HQ, rubredoxin -cmtc: calc_SF's Fourier
#transform 17 % faster with the SF kernels built AVX2 (FLOWOFFICE 13 %).
if(NOS_USE_AVX2 AND LINUX)
    add_compile_options("$<$<COMPILE_LANGUAGE:C,CXX>:-mavx2;-mfma>")
    add_link_options(-mavx2 -mfma)
endif()

if(WIN32 AND CMAKE_CXX_COMPILER_ARCHITECTURE_ID STREQUAL "x64")
    if(NOS_USE_AVX2)
        add_compile_options($<$<COMPILE_LANGUAGE:CXX>:/arch:AVX2>)
    elseif(NOS_USE_AVX)
        add_compile_options($<$<COMPILE_LANGUAGE:CXX>:/arch:AVX>)
    else()
        #the Linux baseline above; MSVC's own default is SSE2
        add_compile_options($<$<COMPILE_LANGUAGE:CXX>:/arch:SSE4.2>)
    endif()
endif()
