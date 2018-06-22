cmake_minimum_required(VERSION 3.0.0)
project(solver3d)
set(CMAKE_BUILD_TYPE Release)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../cmake-modules")
 
cmake_policy(SET CMP0053 OLD)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(PETSc)

if(${CMAKE_VERSION} VERSION_GREATER "3.10.0") 
    find_package(Iconv)
endif()
find_package(LibXml2)

find_package(ZLIB)

file(GLOB SOURCES "*.c")
#Generate the static library from the sources
add_library(solver3d STATIC ${SOURCES})

#Bring the headers into the project
target_include_directories(solver3d PRIVATE)
include_directories(SYSTEM ../mesh3d ${PETSC_INCLUDE_DIR} ${PETSC_INCLUDE_CONF} ${MPI_INCLUDE_PATH} ${LIBXML2_INCLUDE_DIR} ${Iconv_INCLUDE_DIR} ${ZLIB_INCLUDE_DIR})
 
install(TARGETS solver3d DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../../lib)