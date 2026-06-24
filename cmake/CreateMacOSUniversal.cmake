cmake_minimum_required(VERSION 3.25)

foreach(_required IN ITEMS ARM64_DIRECTORY X86_64_DIRECTORY OUTPUT_DIRECTORY)
    if(NOT DEFINED "${_required}")
        message(FATAL_ERROR "${_required} was not provided")
    endif()
endforeach()

file(MAKE_DIRECTORY "${OUTPUT_DIRECTORY}")

function(nosphera2_lipo_file filename)
    set(_arm64_file "${ARM64_DIRECTORY}/${filename}")
    set(_x86_64_file "${X86_64_DIRECTORY}/${filename}")
    set(_output_file "${OUTPUT_DIRECTORY}/${filename}")

    if(NOT EXISTS "${_arm64_file}")
        message(FATAL_ERROR "Missing arm64 file: ${_arm64_file}")
    endif()

    if(NOT EXISTS "${_x86_64_file}")
        message(FATAL_ERROR "Missing x86_64 file: ${_x86_64_file}")
    endif()

    execute_process(
        COMMAND
            /usr/bin/lipo
            -create
            "${_arm64_file}"
            "${_x86_64_file}"
            -output
            "${_output_file}"
        COMMAND_ERROR_IS_FATAL ANY
    )

    execute_process(
        COMMAND
            /usr/bin/lipo
            -archs
            "${_output_file}"
        OUTPUT_VARIABLE _architectures
        OUTPUT_STRIP_TRAILING_WHITESPACE
        COMMAND_ERROR_IS_FATAL ANY
    )

    message(STATUS "${filename}: ${_architectures}")
endfunction()

nosphera2_lipo_file("NoSpherA2")
nosphera2_lipo_file("libtbb.12.dylib")
nosphera2_lipo_file("libomp.dylib")