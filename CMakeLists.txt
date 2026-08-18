cmake_minimum_required(VERSION 3.10)

project(Autoparts)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)

find_package(Eigen3 QUIET)

if(NOT Eigen3_FOUND)
    FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0
    )
    FetchContent_MakeAvailable(eigen)
endif()

file(GLOB_RECURSE SOURCES
    "${PROJECT_SOURCE_DIR}/src/*.cpp"
)

add_executable(autoparts ${SOURCES})

target_include_directories(autoparts PRIVATE
    ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(autoparts PRIVATE
    Eigen3::Eigen
)

target_compile_options(autoparts PRIVATE
    -march=native
    -O3
)