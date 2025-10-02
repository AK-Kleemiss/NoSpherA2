add_executable(BasisSetConverter
    "${CMAKE_CURRENT_SOURCE_DIR}/basis_set_helper.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/BasisSetConverter.cpp"
)

target_include_directories(BasisSetConverter PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

add_library(basis_set_helper STATIC
    "${CMAKE_CURRENT_SOURCE_DIR}/basis_set_helper.cpp"
)
target_include_directories(basis_set_helper PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
if (WIN32)
    set(BSC_EXE "${CMAKE_CURRENT_BINARY_DIR}/BasisSetConverter.exe")
else ()
    set(BSC_EXE "${CMAKE_CURRENT_BINARY_DIR}/BasisSetConverter")
endif ()
add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/Src/auxiliary_basis.cpp
    COMMAND $<TARGET_FILE:BasisSetConverter> ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_BINARY_DIR}/Src
    DEPENDS BasisSetConverter
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating auxiliary_basis.cpp using BasisSetConverter"
    VERBATIM
)

add_custom_target(generate_auxiliary_basis
    DEPENDS ${CMAKE_BINARY_DIR}/Src/auxiliary_basis.cpp
)

