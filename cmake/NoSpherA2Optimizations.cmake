include_guard(GLOBAL)

include(CheckIPOSupported)

check_ipo_supported(
    RESULT NOSPHERA2_IPO_SUPPORTED
    OUTPUT NOSPHERA2_IPO_ERROR
    LANGUAGES CXX
)

if(NOT NOSPHERA2_IPO_SUPPORTED)
    message(
        WARNING
        "IPO/LTO is unavailable: ${NOSPHERA2_IPO_ERROR}"
    )
endif()

#clang PGO: build with GENERATE, run NoSpherA2_Tests --gtest_filter=TomlIntegrationTests.*, llvm-profdata merge, rebuild with USE
option(NOSPHERA2_PGO_GENERATE "Instrument for clang profile-guided optimisation" OFF)
set(NOSPHERA2_PGO_USE "" CACHE FILEPATH "Merged .profdata for clang profile-guided optimisation")
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
    set(NOSPHERA2_LIBMVEC ON)
endif()
#lld-link is called directly; the objects name clang_rt.profile but not where it lives
if(NOSPHERA2_PGO_GENERATE AND MSVC AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    execute_process(COMMAND "${CMAKE_CXX_COMPILER}" /clang:-print-runtime-dir OUTPUT_VARIABLE NOSPHERA2_CLANG_RT_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

function(nosphera2_enable_optimizations target_name)
    if(NOT TARGET "${target_name}")
        message(FATAL_ERROR
            "Target '${target_name}' does not exist"
        )
    endif()

    # Link-time/interprocedural optimization.
    if(NOSPHERA2_IPO_SUPPORTED)
        set_property(
            TARGET "${target_name}"
            PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE
        )
    endif()

    if(MSVC)
        target_compile_options(
            "${target_name}"
            PRIVATE
                $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:
                    /fp:fast
                    /fp:except-
                    /Zc:inline
                >
                #clang-cl: /openmp:experimental emits __kmpc_* calls that MKL's libiomp5md serves, /Qpar it ignores.
                #its /fp:fast is -ffast-math, whose finite-math-only folds every isnan/isfinite away; MSVC's keeps them
                $<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:CXX,MSVC>>:/Qpar>
                $<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:CXX,Clang>>:/clang:-fno-finite-math-only /clang:-O3>
                $<$<AND:$<BOOL:${NOSPHERA2_PGO_GENERATE}>,$<COMPILE_LANG_AND_ID:CXX,Clang>>:/clang:-fprofile-generate>
                $<$<AND:$<BOOL:${NOSPHERA2_PGO_USE}>,$<COMPILE_LANG_AND_ID:CXX,Clang>>:/clang:-fprofile-use=${NOSPHERA2_PGO_USE} /clang:-Wno-profile-instr-unprofiled /clang:-Wno-profile-instr-out-of-date /clang:-Wno-backend-plugin>
                $<$<COMPILE_LANGUAGE:CXX>:/openmp:experimental>
                # unoptimised Debug objects outgrow the COFF section limit (C1128)
                $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CXX>>:/bigobj>
                $<$<AND:$<BOOL:${NOSPHERA2_PROFILE_SYMBOLS}>,$<COMPILE_LANGUAGE:CXX>>:/Zi>
        )
        if(NOSPHERA2_PROFILE_SYMBOLS)
            target_link_options("${target_name}" PRIVATE "$<HOST_LINK:/DEBUG>")
        endif()
        if(NOSPHERA2_CLANG_RT_DIR)
            target_link_options("${target_name}" PRIVATE "$<HOST_LINK:/LIBPATH:${NOSPHERA2_CLANG_RT_DIR}>")
        endif()

        #HOST_LINK, because a CUDA target links twice. Link options otherwise reach the
        #device link as well, where nvcc hands what it does not recognise to cl - and cl
        #reads /OPT:REF as its own /O flag followed by rubbish, warning once per character.
        #They were ignored there rather than misapplied, so this was noise rather than a
        #defect, but it buried real warnings and made every build look unclean.
        #ARM64 has no libiomp5; there vcomp is the OpenMP runtime
        if(NOT NOSPHERA2_OPENBLAS)
            target_link_options(
                "${target_name}"
                PRIVATE
                    "$<HOST_LINK:/NODEFAULTLIB:vcomp>"
                    "$<HOST_LINK:/NODEFAULTLIB:vcompd>"
            )
        endif()

        get_target_property(target_type "${target_name}" TYPE)

        if(NOT target_type STREQUAL "STATIC_LIBRARY")
            target_link_options(
                "${target_name}"
                PRIVATE
                    "$<HOST_LINK:$<$<CONFIG:Release>:/OPT:REF>>"
                    "$<HOST_LINK:$<$<CONFIG:Release>:/OPT:ICF>>"
                    "$<HOST_LINK:$<$<CONFIG:Release>:/INCREMENTAL:NO>>"
            )
        endif()

    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        # Put individual functions and data objects in separate sections.
        target_compile_options(
            "${target_name}"
            PRIVATE
                $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:
                    -ffunction-sections
                    -fdata-sections
                >
                #MSVC's /fp:fast without the NaN/Inf folding; compile only, so no crtfastmath.o sets FTZ/DAZ
                #process-wide in a host that loads the library. gcc, the test-only build, stays strict IEEE
                $<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:CXX,Clang,AppleClang>>:
                    -ffast-math
                    -fno-finite-math-only
                    #clang 18 lets -ffast-math flush denormals to zero in the code itself (get_bessel_ratio's DBL_MIN rescale gave NaN)
                    -fdenormal-fp-math=ieee
                >
                #vectorised exp/log/sin/pow: glibc's libmvec (inside libm), Apple's Accelerate (linked already)
                $<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:CXX,Clang>,$<BOOL:${NOSPHERA2_LIBMVEC}>>:-fveclib=libmvec>
                $<$<AND:$<CONFIG:Release>,$<COMPILE_LANG_AND_ID:CXX,Clang,AppleClang>,$<BOOL:${APPLE}>>:-fveclib=Accelerate>
                $<$<AND:$<BOOL:${NOSPHERA2_PGO_GENERATE}>,$<COMPILE_LANG_AND_ID:CXX,Clang,AppleClang>>:-fprofile-generate>
                $<$<AND:$<BOOL:${NOSPHERA2_PGO_USE}>,$<COMPILE_LANG_AND_ID:CXX,Clang,AppleClang>>:
                    -fprofile-use=${NOSPHERA2_PGO_USE}
                    -Wno-profile-instr-unprofiled
                    -Wno-profile-instr-out-of-date
                >
        )
        #only nosphera2_run/nosphera2_main (NOS_API) leave the shared library
        set_target_properties("${target_name}" PROPERTIES CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON)
        if(NOSPHERA2_PGO_GENERATE)
            target_link_options("${target_name}" PRIVATE -fprofile-generate)
        endif()

        # -fopenmp has to be a compile option on every target, static library included.
        # It used to sit inside the "not a static library" branch below, next to the link
        # options where it does not belong: NoSpherA2Core is a static library and holds
        # every omp pragma in the project, so on Linux they were all compiled away. The
        # executables got the flag and contain almost no parallel code, which is why this
        # was invisible. Measured on the cluster before the fix: the scattering-factor
        # transform ran at 0.5 GFLOP/s against 22.1 on a laptop, and OMP_NUM_THREADS=48 and
        # =1 gave the same runtime to within a second. MSVC was never affected because
        # /openmp:experimental above is applied unconditionally.
        #
        # The compile flag is deliberately not matched by -fopenmp at link time: that
        # makes gcc pull in libgomp, and the process then holds two OpenMP runtimes,
        # since MKL brings libiomp5. libiomp5 implements the GOMP_* entry points that
        # -fopenmp generates, so linking MKL is enough to resolve them and the build
        # keeps a single runtime, which is how it worked before the link flag was
        # added alongside the compile one.
        # AppleClang rejects a bare -fopenmp ("unsupported option '-fopenmp'"), and
        # even a clang that accepts it would not find omp.h or libomp, both of which
        # live in the Micromamba environment. The imported target carries what that
        # compiler needs - -Xpreprocessor -fopenmp, the include directory and the
        # libomp path, all set by the macOS presets - so use it there. Linking it to
        # a static library adds no library but does hand the target its compile
        # options and include directories, which is the half that matters here.
        if(APPLE)
            target_link_libraries(
                "${target_name}"
                PRIVATE
                    OpenMP::OpenMP_CXX
            )
        else()
            target_compile_options(
                "${target_name}"
                PRIVATE
                    $<$<COMPILE_LANGUAGE:CXX>:-fopenmp>
            )
        endif()

        if(NOSPHERA2_COVERAGE)
            target_compile_options("${target_name}" PRIVATE $<$<COMPILE_LANGUAGE:CXX>:--coverage>)
            target_link_options("${target_name}" PRIVATE --coverage)
        endif()

        get_target_property(target_type "${target_name}" TYPE)

        if(NOT target_type STREQUAL "STATIC_LIBRARY")
            if(APPLE)
                target_link_options(
                    "${target_name}"
                    PRIVATE
                        $<$<CONFIG:Release>:
                            LINKER:-dead_strip
                        >
                )
            else()
                target_link_options(
                    "${target_name}"
                    PRIVATE
                        $<$<CONFIG:Release>:
                            LINKER:--gc-sections
                        >
                )
                #lld-link already folds with /OPT:ICF; safe = only functions whose address is never taken
                #LINKER_ID, not CMAKE_LINKER_TYPE: CMake < 3.29 ignores the latter and links with ld.bfd
                if(CMAKE_CXX_COMPILER_LINKER_ID STREQUAL "LLD")
                    target_link_options("${target_name}" PRIVATE $<$<CONFIG:Release>:LINKER:--icf=safe>)
                endif()
            endif()
            # Android has no MKL to bring libiomp5: the NDK's libomp, linked in statically
            if(ANDROID)
                target_link_options("${target_name}" PRIVATE -fopenmp -static-openmp)
            endif()
        endif()
    endif()
endfunction()