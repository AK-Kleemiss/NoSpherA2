# AVX detection / compile options for NoSpherA2 + OCC.
include_guard(GLOBAL)

# Vector instruction policy. OCC and NoSpherA2 MUST be built with the same
# setting: Eigen object layouts differ between AVX and non-AVX builds.
#   NOS_AVX=ON   build with AVX
#   NOS_AVX=OFF  SSE only — maximum compatibility, used by CI
#   unset        honor the NOS_AVX environment variable, else detect the host
if(NOT DEFINED NOS_AVX)
    if(DEFINED ENV{NOS_AVX} AND NOT "$ENV{NOS_AVX}" STREQUAL "")
        set(NOS_AVX "$ENV{NOS_AVX}")
    else()
        set(NOS_AVX "AUTO")
    endif()
endif()

set(NOS_USE_AVX OFF)
if(NOS_AVX STREQUAL "AUTO")
    if(NOT CMAKE_CROSSCOMPILING AND CMAKE_SYSTEM_PROCESSOR MATCHES "(x86_64|AMD64|amd64)")
        try_run(NOS_AVX_RUN_RESULT NOS_AVX_COMPILE_RESULT
            ${CMAKE_BINARY_DIR}/nos_avx_check
            ${CMAKE_CURRENT_SOURCE_DIR}/cmake/detect_avx.c)
        if(NOS_AVX_COMPILE_RESULT AND NOS_AVX_RUN_RESULT EQUAL 0)
            set(NOS_USE_AVX ON)
        endif()
    endif()
elseif(NOS_AVX)
    set(NOS_USE_AVX ON)
endif()
message(STATUS "NOS_AVX=${NOS_AVX} -> building with AVX: ${NOS_USE_AVX}")

if(LINUX)
    add_compile_options(
        -msse2
        -msse3
        -msse4.1
        -msse4.2
    )
    if(NOS_USE_AVX)
        add_compile_options(-mavx)
    endif()
endif()

#NOS_AVX2=ON goes one step further: -mavx2 -mfma, i.e. 256-bit packed fp64 with a fused
#multiply-add. It has to be a LINK option as well as a compile option, and that is not a detail:
#the build is -flto=auto, so the machine code is generated at link time and per-file compile flags
#produced a byte-identical binary when this was tried per translation unit. Off by default - an AVX2
#instruction on a pre-Haswell host is SIGILL, not a slow path - and NOT wired into the AUTO
#detection above until the win is measured on a real case.
if(NOS_AVX2 AND LINUX)
    add_compile_options(-mavx2 -mfma)
    add_link_options(-mavx2 -mfma)
endif()

if(WIN32 AND NOS_USE_AVX)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:/arch:AVX>)
endif()
