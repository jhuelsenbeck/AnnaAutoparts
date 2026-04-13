cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(PartitionDistanceUtility)

include_directories(${PROJECT_SOURCE_DIR})

file(GLOB SOURCES "${PROJECT_SOURCE_DIR}/*.cpp")

add_executable(partition_distance ${SOURCES})

target_include_directories(partition_distance PRIVATE ${PROJECT_SOURCE_DIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O3")
