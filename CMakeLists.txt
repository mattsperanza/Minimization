cmake_minimum_required(VERSION 3.29)
project(Minimization)

set(CMAKE_CXX_STANDARD 11)

add_executable(
        Minimization
        main.cpp
        src/LBFGS.cpp
        src/LBFGS.h
        src/Minimizer.cpp
        src/Minimizer.h
)
