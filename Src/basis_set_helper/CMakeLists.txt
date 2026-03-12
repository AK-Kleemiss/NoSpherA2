add_executable(BasisSetConverter EXCLUDE_FROM_ALL
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
file(GLOB BASIS_SET_CSVS CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/basis_sets/*.csv"
)

add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/Src/auxiliary_basis.cpp
    COMMAND $<TARGET_FILE:BasisSetConverter> ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_BINARY_DIR}/Src
    DEPENDS BasisSetConverter ${BASIS_SET_CSVS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating auxiliary_basis.cpp using BasisSetConverter"
    VERBATIM
)

add_custom_target(generate_auxiliary_basis
    DEPENDS ${CMAKE_BINARY_DIR}/Src/auxiliary_basis.cpp
)

