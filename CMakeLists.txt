cmake_minimum_required(VERSION 3.17)
project(approx)

set(CMAKE_CXX_STANDARD 20)

add_executable(approx examples/example.cpp include/approx/approx.h src/approx/internals.h)