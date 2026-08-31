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
                    /Qpar
                    /Zc:inline
                >
                $<$<COMPILE_LANGUAGE:CXX>:/openmp:experimental>
        )

        #HOST_LINK, because a CUDA target links twice. Link options otherwise reach the
        #device link as well, where nvcc hands what it does not recognise to cl - and cl
        #reads /OPT:REF as its own /O flag followed by rubbish, warning once per character.
        #They were ignored there rather than misapplied, so this was noise rather than a
        #defect, but it buried real warnings and made every build look unclean.
        target_link_options(
            "${target_name}"
            PRIVATE
                "$<HOST_LINK:/NODEFAULTLIB:vcomp>"
                "$<HOST_LINK:/NODEFAULTLIB:vcompd>"
        )

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
        )

        # -fopenmp has to be a compile option on every target, static library included.
        # It used to sit inside the "not a static library" branch below, next to the link
        # options where it does not belong: NoSpherA2Core is a static library and holds
        # every omp pragma in the project, so on Linux they were all compiled away. The
        # executables got the flag and contain almost no parallel code, which is why this
        # was invisible. Measured on the cluster before the fix: the scattering-factor
        # transform ran at 0.5 GFLOP/s against 22.1 on a laptop, and OMP_NUM_THREADS=48 and
        # =1 gave the same runtime to within a second. MSVC was never affected because
        # /openmp:experimental above is applied unconditionally.
        target_compile_options(
            "${target_name}"
            PRIVATE
                $<$<COMPILE_LANGUAGE:CXX>:-fopenmp>
        )

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
                target_link_libraries(
                    "${target_name}"
                    PRIVATE
                        OpenMP::OpenMP_CXX
                )
            else()
                target_link_options(
                    "${target_name}"
                    PRIVATE
                        $<$<CONFIG:Release>:
                            LINKER:--gc-sections
                        >
                        # Linking the runtime is the half that genuinely only applies to
                        # things that get linked; the compile flag moved out of here.
                        -fopenmp
                )
            endif()
        endif()
    endif()
endfunction()