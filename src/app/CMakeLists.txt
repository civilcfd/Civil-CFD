cmake_minimum_required(VERSION 3.0.0)
project(app)
set(CMAKE_BUILD_TYPE Release)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../cmake-modules")
include(CheckFunctionExists)

if(MSVC)
    set(QHULL_USE_STATIC TRUE)
endif(MSVC)

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

find_package(Qhull)

file(GLOB SOURCES "*.c")
#Generate the static library from the sources
add_executable(solver3d-bin solver3d.c)
add_executable(mesh3d-bin mesh3d.c)
add_executable(inspect_cell inspect_cell.c)

#Bring the headers into the project
include_directories(SYSTEM ../mesh3d ../solver3d ${PETSC_INCLUDE_DIR} ${PETSC_INCLUDE_CONF} ${MPI_INCLUDE_PATH} ${LIBXML2_INCLUDE_DIR} ${Iconv_INCLUDE_DIR} ${ZLIB_INCLUDE_DIR}  ${QHULL_INCLUDE_DIRS})
 

CHECK_FUNCTION_EXISTS(pow RESULT)
if(NOT RESULT)
#  unset(RESULT)
#  list(APPEND CMAKE_REQUIRED_LIBRARIES m)
#  CHECK_FUNCTION_EXISTS(pow RESULT)
#  if(RESULT)
    target_link_libraries(solver3d-bin m)
    target_link_libraries(mesh3d-bin m)
    target_link_libraries(inspect_cell m)
#  else()
#    message(FATAL_ERROR "No pow() found")
#  endif()
endif()

if(WIN32)
  target_link_libraries(solver3d-bin ws2_32)
  target_link_libraries(mesh3d-bin ws2_32)
  target_link_libraries(inspect_cell ws2_32)
endif()

target_link_libraries(solver3d-bin solver3d)
target_link_libraries(solver3d-bin mesh3d)
target_link_libraries(solver3d-bin ${MPI_LIBRARIES})
target_link_libraries(solver3d-bin ${OPENMP_LIBRARIES})
target_link_libraries(solver3d-bin ${LIBXML2_LIBRARIES})
target_link_libraries(solver3d-bin ${Iconv_LIBRARIES})
target_link_libraries(solver3d-bin ${PETSC_LIBRARIES})
target_link_libraries(solver3d-bin ${ZLIB_LIBRARIES})

target_link_libraries(mesh3d-bin solver3d)
target_link_libraries(mesh3d-bin mesh3d)
target_link_libraries(mesh3d-bin ${MPI_LIBRARIES})
target_link_libraries(mesh3d-bin ${OPENMP_LIBRARIES})
target_link_libraries(mesh3d-bin ${LIBXML2_LIBRARIES})
target_link_libraries(mesh3d-bin ${Iconv_LIBRARIES})
target_link_libraries(mesh3d-bin ${PETSC_LIBRARIES})
target_link_libraries(mesh3d-bin ${ZLIB_LIBRARIES})
target_link_libraries(mesh3d-bin ${QHULL_LIBRARIES})

target_link_libraries(inspect_cell solver3d)
target_link_libraries(inspect_cell mesh3d)
target_link_libraries(inspect_cell ${MPI_LIBRARIES})
target_link_libraries(inspect_cell ${OPENMP_LIBRARIES})
target_link_libraries(inspect_cell ${LIBXML2_LIBRARIES})
target_link_libraries(inspect_cell ${Iconv_LIBRARIES})
target_link_libraries(inspect_cell ${PETSC_LIBRARIES})
target_link_libraries(inspect_cell ${ZLIB_LIBRARIES})


set_target_properties(solver3d-bin PROPERTIES OUTPUT_NAME solver3d)
set_target_properties(mesh3d-bin PROPERTIES OUTPUT_NAME mesh3d)

install(TARGETS solver3d-bin DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../../bin)
install(TARGETS mesh3d-bin DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../../bin)
install(TARGETS inspect_cell DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../../bin)