include_guard(GLOBAL)

include(MicromambaEnvironment)

function(nosphera2_setup_dependencies)
    setup_micromamba_environment(
        ENVIRONMENT_FILE
            "${PROJECT_SOURCE_DIR}/environment.yaml"
    )

    # CMakePresets.json normally provides this explicitly. This lookup also
    # validates that Python is usable from the environment.
    find_package(
        Python3 3.12
        REQUIRED
        COMPONENTS Interpreter
    )

    message(STATUS "Micromamba environment: ${MICROMAMBA_ENV_PREFIX}")
    message(STATUS "Python executable: ${Python3_EXECUTABLE}")
endfunction()
