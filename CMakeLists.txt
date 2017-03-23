project(ExtendedKF)

cmake_minimum_required (VERSION 3.5)

add_definitions(-std=c++0x)


FILE(GLOB_RECURSE LibFiles src/*)
add_custom_target(headers SOURCES ${LibFiles})

set(sources
    src/FusionEKF.cpp
    src/kalman_filter.cpp
    src/main.cpp
    src/tools.cpp)

add_executable(ExtendedKF ${sources})
