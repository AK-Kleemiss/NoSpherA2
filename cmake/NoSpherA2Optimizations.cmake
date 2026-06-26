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
                /openmp:experimental
        )

        get_target_property(target_type "${target_name}" TYPE)

        if(NOT target_type STREQUAL "STATIC_LIBRARY")
            target_link_options(
                "${target_name}"
                PRIVATE
                    $<$<CONFIG:Release>:
                        /OPT:REF
                        /OPT:ICF
                        /INCREMENTAL:NO
                    >
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
            endif()
            target_link_options(
                "${target_name}"
                PRIVATE
                    OpenMP::OpenMP_CXX
                )
        endif()
    endif()
endfunction()