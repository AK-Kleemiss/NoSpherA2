add_executable(convert convert.cpp)
target_link_libraries(convert libnosphera2 occ::libocc)
target_compile_options(convert PRIVATE "-march=native")
