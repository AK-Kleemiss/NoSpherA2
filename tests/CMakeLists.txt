include(CTest)
find_package(Python)

if (WIN32)
    set(NoSpherA2_EXE ${CMAKE_BINARY_DIR}/Src/NoSpherA2.exe)
else ()
    set(NoSpherA2_EXE ${CMAKE_BINARY_DIR}/Src/NoSpherA2)
endif ()

message(WARNING ${CMAKE_CURRENT_BINARY_DIR})
FILE(COPY ${CMAKE_CURRENT_SOURCE_DIR}/alanine_occ DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_test(
    NAME alanine_occ
    COMMAND ${Python_EXECUTABLE} runTest.py "NoSpherA2_EXE -cif alanine.cif -dmin 0.5 -wfn alanine.owf.fchk -acc 1 -no-date"
)
