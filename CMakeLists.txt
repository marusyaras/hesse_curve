cmake_minimum_required(VERSION 3.17)
project(hesse_curve)

set(CMAKE_CXX_STANDARD 20)

add_executable(hesse_curve main.cpp)

target_include_directories(hesse_curve PRIVATE .)

find_package(PkgConfig REQUIRED)
pkg_check_modules(tommath REQUIRED IMPORTED_TARGET libtommath)
target_link_libraries(hesse_curve PRIVATE PkgConfig::tommath)

