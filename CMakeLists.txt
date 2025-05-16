cmake_minimum_required(VERSION 3.28)
project(software)

set(CMAKE_CXX_STANDARD 17)

add_executable(software main.cpp
        Graph.cpp
        Config.cpp
        Metrics.cpp
        AnynetGenerator.hpp
)