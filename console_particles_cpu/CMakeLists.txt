cmake_minimum_required(VERSION 3.0)

project(particles VERSION 1.0 LANGUAGES CXX)

add_compile_options(
    -std=c++11
    -O3
    -Werror
    -Wall
    -Wextra
    -Wpedantic
    -Wcast-align
    -Wcast-qual
    -Wconversion
    -Wctor-dtor-privacy
    -Wenum-compare
    -Wfloat-equal
    -Wnon-virtual-dtor
    -Wold-style-cast
    -Woverloaded-virtual
    -Wredundant-decls
    -Wsign-conversion
    -Wsign-promo
)

set(SRC_LIST
    ./src/particle.cpp
    ./src/vec3.cpp
)

add_executable(particles src/main.cpp ${SRC_LIST})
add_executable(getcoulombforce src/getcoulombforce.cpp ${SRC_LIST})

target_include_directories(particles PRIVATE ./include/)
target_include_directories(getcoulombforce PRIVATE ./include/)
