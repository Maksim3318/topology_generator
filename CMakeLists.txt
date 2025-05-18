cmake_minimum_required(VERSION 3.28)
project(topology_generator)

set(CMAKE_CXX_STANDARD 17)

add_executable(topology_generator
        main.cpp
        Graph.cpp
        Config.cpp
        Metrics.cpp
        AnynetGenerator.hpp
)