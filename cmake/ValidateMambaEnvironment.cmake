include_guard(GLOBAL)

function(validate_micromamba_environment)
    if(
        NOT DEFINED MICROMAMBA_ENV_PREFIX
        OR MICROMAMBA_ENV_PREFIX STREQUAL ""
    )
        message(FATAL_ERROR
            "MICROMAMBA_ENV_PREFIX is not set.\n"
            "Configure through a CMake preset."
        )
    endif()

    if(NOT EXISTS "${MICROMAMBA_ENV_PREFIX}/conda-meta/history")
        message(FATAL_ERROR
            "The selected micromamba environment is missing:\n"
            "  ${MICROMAMBA_ENV_PREFIX}\n\n"
            "Run first:\n"
            "  cmake -P scripts/BootstrapMicromamba.cmake"
        )
    endif()

    message(STATUS
        "Using micromamba environment: ${MICROMAMBA_ENV_PREFIX}"
    )
endfunction()