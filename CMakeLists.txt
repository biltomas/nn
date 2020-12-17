cmake_minimum_required(VERSION 3.0)
project(pv021)

set(CMAKE_CXX_COMPILER "clang++")

if (CMAKE_CXX_COMPILER_ID MATCHES Clang OR ${CMAKE_CXX_COMPILER_ID} STREQUAL GNU)
    # using regular Clang or AppleClang or GCC
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -g -Wall -Wextra -pedantic -Wold-style-cast")
elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    # using Visual Studio C++
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++14 /W4 /permissive-")
    add_definitions("-D_CRT_SECURE_NO_WARNINGS")
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS on)

add_executable(pv021_nn main.cpp NeuralNetwork.cpp)
