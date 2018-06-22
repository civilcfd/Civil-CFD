cmake_minimum_required (VERSION 3.0.0)
project (gui)
set(CMAKE_BUILD_TYPE Release)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake-modules")
include(CheckFunctionExists)

if(MSVC)
    set(QHULL_USE_STATIC TRUE)
endif(MSVC)

cmake_policy(SET CMP0053 OLD)

find_package(VTK 8.1.0 REQUIRED)
IF(VTK_FOUND)
  message("found VTK. Version:" ${VTK_VERSION}. VTK_DIR: ${VTK_DIR})
ELSE(VTK_FOUND)
  MESSAGE(FATAL_ERROR
    "Cannot build the executable without VTK. Please set the                                                                                                   
VTK_DIR")
ENDIF(VTK_FOUND)
include(${VTK_USE_FILE})

set(VTK_IGNORE_GLDRIVER_BUGS:BOOL=ON)

# Find the QtWidgets library
find_package(Qt5Widgets)
find_package(Qt5PrintSupport REQUIRED)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

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

if(WIN32)
  set (GUI "WIN32")
endif()

# Tell CMake to create the helloworld executable
add_executable(Civil-CFD ${GUI} main.cpp MainWindow.cpp MeshDisplay.cpp GeometryDisplay.cpp Render.cpp RenderDisplay.cpp MainWindow_Boundaries.cpp Boundary.cpp BoundaryDisplay.cpp Solver.cpp ResultList.cpp Visualize.cpp VisualizeDisplay.cpp Baffle.cpp MainWindow_Baffles.cpp qcustomplot.cpp Visualize3D.cpp Visualize3DDisplay.cpp About.cpp InspectDialog.cpp Simulation.cpp SolverDialog.h ../src/solver3d/solver.h ../src/mesh3d/mesh.h gui.qrc civilcfd.rc)

# Include directories
target_include_directories(Civil-CFD PRIVATE ../src/solver3d ../src/mesh3d)
include_directories(SYSTEM ${MPI_INCLUDE_PATH} ${LIBXML2_INCLUDE_DIR} ${ZLIB_INCLUDE_DIR} ${Iconv_INCLUDE_DIR})

if(WIN32)
  target_link_libraries(Civil-CFD ws2_32)
endif()

# Use the Widgets module from Qt 5.
target_link_libraries(Civil-CFD solver3d)
target_link_libraries(Civil-CFD mesh3d)
target_link_libraries(Civil-CFD Qt5::Widgets)
target_link_libraries(Civil-CFD Qt5::PrintSupport)
target_link_libraries(Civil-CFD ${VTK_LIBRARIES})
target_link_libraries(Civil-CFD ${MPI_LIBRARIES})
target_link_libraries(Civil-CFD ${OPENMP_LIBRARIES})
target_link_libraries(Civil-CFD ${LIBXML2_LIBRARIES})
target_link_libraries(Civil-CFD ${Iconv_LIBRARIES})
target_link_libraries(Civil-CFD ${PETSC_LIBRARIES})
target_link_libraries(Civil-CFD ${ZLIB_LIBRARIES})
target_link_libraries(Civil-CFD ${QHULL_LIBRARIES})

CHECK_FUNCTION_EXISTS(pow pow_exist)
if(NOT pow_exist)
#  unset(pow_exist)
#  list(APPEND CMAKE_REQUIRED_LIBRARIES m)
#  CHECK_FUNCTION_EXISTS(pow pow_exist)
#  if(pow_exist)
    target_link_libraries(Civil-CFD m)
#  else()
#    message(FATAL_ERROR "No pow() found")
#  endif()
endif()

install(TARGETS Civil-CFD DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../bin)