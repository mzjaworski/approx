cmake_minimum_required(VERSION 3.17)
project(approx VERSION 0.1)

add_library(approx INTERFACE)

target_compile_features(approx INTERFACE cxx_std_20)

target_include_directories(approx
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/src)
