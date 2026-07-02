# Sanitizer configuration for NoSpherA2 (Debug-only).
include_guard(GLOBAL)

if(
    NOT "${USE_SANITIZER}" STREQUAL ""
    AND NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Debug"
)
    message(
        FATAL_ERROR
        "Sanitizer should only be used in Debug mode. Without debug symbols it is really hard to know where the issue is."
    )
endif()

if(NOT "${USE_SANITIZER}" STREQUAL "" AND MSVC)
    message(
        FATAL_ERROR
        "Sanitizers configured via USE_SANITIZER are not supported with the MSVC toolchain in this project."
    )
endif()

if(USE_SANITIZER MATCHES "Address")
    add_compile_options(
        $<$<CONFIG:Debug>:-fsanitize=address,leak,undefined>
        $<$<CONFIG:Debug>:-fno-omit-frame-pointer>
        $<$<CONFIG:Debug>:-O0>
    )
    add_link_options($<$<CONFIG:Debug>:-fsanitize=address,leak,undefined>)
elseif(USE_SANITIZER MATCHES "Thread")
    add_compile_options(
        $<$<CONFIG:Debug>:-fsanitize=thread>
        $<$<CONFIG:Debug>:-fno-omit-frame-pointer>
        $<$<CONFIG:Debug>:-O0>
    )
    add_link_options($<$<CONFIG:Debug>:-fsanitize=thread>)
elseif(USE_SANITIZER MATCHES "Memory")
    add_compile_options(
        $<$<CONFIG:Debug>:-fsanitize=memory>
        $<$<CONFIG:Debug>:-fsanitize-memory-track-origins>
        $<$<CONFIG:Debug>:-fno-omit-frame-pointer>
        $<$<CONFIG:Debug>:-O0>
    )
    add_link_options($<$<CONFIG:Debug>:-fsanitize=memory>)
endif()
