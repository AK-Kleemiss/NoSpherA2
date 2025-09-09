cmake_minimum_required(VERSION 3.14)
project(NoSpherA2 CXX C)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(BUILD_DIR ${CMAKE_SOURCE_DIR}/build)
set(DEBUG_DIR ${CMAKE_SOURCE_DIR}/build_debug)
set(CMAKE_INSTALL_PREFIX $ENV{HOME})
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/CPM.cmake")
add_subdirectory(3rdparty EXCLUDE_FROM_ALL)
add_subdirectory(Src)
enable_testing()
add_subdirectory(tests)
