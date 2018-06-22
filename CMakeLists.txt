cmake_minimum_required(VERSION 3.0.0)
project(Civil-CFD)
set(CMAKE_BUILD_TYPE Release)
list(APPEND CMAKE_MODULE_PATH "cmake-modules")
include(CheckFunctionExists)


add_subdirectory(src/mesh3d)
add_subdirectory(src/solver3d)
add_subdirectory(src/app)
add_subdirectory(gui)