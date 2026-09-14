file(
    GLOB NOSPHERA2_APP_SOURCES
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx"
)

add_executable(
    NoSpherA2
    ${NOSPHERA2_APP_SOURCES}
)

if(WIN32)
    target_sources(
        NoSpherA2
        PRIVATE
            "${PROJECT_SOURCE_DIR}/Windows/NoSpherA2/Resource.rc"
    )

    target_include_directories(
        NoSpherA2
        PRIVATE
            "${PROJECT_SOURCE_DIR}/Windows/NoSpherA2"
    )
endif()

nosphera2_enable_optimizations(NoSpherA2)

target_compile_features(
    NoSpherA2
    PRIVATE
        cxx_std_20
)

target_link_libraries(
    NoSpherA2
    PRIVATE
        NoSpherA2::Core
)

set_target_properties(
    NoSpherA2
    PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY
            "${CMAKE_BINARY_DIR}/bin"
)


nosphera2_copy_runtime_libraries(NoSpherA2)

if(NOSPHERA2_CUDA_DELAYLOAD)
    target_link_options(NoSpherA2 PRIVATE ${NOSPHERA2_CUDA_DELAYLOAD})
    target_link_libraries(NoSpherA2 PRIVATE delayimp)
endif()

if(NOSPHERA2_USE_HIP)
    #A static library with HIP objects makes CMake link this with the HIP compiler and
    #the --offload-arch list, and clang then unbundles every object on the line - a plain
    #host object comes out empty, so main() goes missing. The kernels carry their own
    #fat binaries (no separable compilation), so the C++ linker is all that is needed.
    set_target_properties(NoSpherA2 PROPERTIES LINKER_LANGUAGE CXX)
endif()

if(NOSPHERA2_HIP_DELAYLOAD)
    target_link_options(NoSpherA2 PRIVATE ${NOSPHERA2_HIP_DELAYLOAD})
    target_link_libraries(NoSpherA2 PRIVATE delayimp)
    #The failure hook has to be an object of this link, not of the core library - see the file.
    target_sources(NoSpherA2 PRIVATE "${PROJECT_SOURCE_DIR}/Src/hip_delayload_hook.cpp")
endif()
