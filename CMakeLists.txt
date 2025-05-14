cmake_minimum_required(VERSION 3.29)
project(LBFGS)

set(CMAKE_CXX_STANDARD 11)

add_executable(
        LBFGS
        src/LBFGS.h
        main.cpp)
