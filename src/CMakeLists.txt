cmake_minimum_required(VERSION 3.16)
project(Voxel2Tet LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(VOXEL2TET_EXTERNAL_INCLUDE_DIRS "")
set(VOXEL2TET_EXTERNAL_LIBS "")

# Look for required libraries

#### HDF5
find_package(HDF5 REQUIRED COMPONENTS CXX)
if (HDF5_INCLUDE_DIRS)
  list(APPEND VOXEL2TET_EXTERNAL_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
elseif (HDF5_INCLUDE_DIR)
  list(APPEND VOXEL2TET_EXTERNAL_INCLUDE_DIRS ${HDF5_INCLUDE_DIR})
endif ()
list(APPEND VOXEL2TET_EXTERNAL_LIBS ${HDF5_LIBRARIES})
message(STATUS "HDF5 root: ${HDF5_ROOT}")

#### ARMADILLO
find_package(Armadillo REQUIRED)
if (TARGET Armadillo::armadillo)
  list(APPEND VOXEL2TET_EXTERNAL_LIBS Armadillo::armadillo)
else ()
  if (ARMADILLO_INCLUDE_DIRS)
    list(APPEND VOXEL2TET_EXTERNAL_INCLUDE_DIRS ${ARMADILLO_INCLUDE_DIRS})
  endif ()
  list(APPEND VOXEL2TET_EXTERNAL_LIBS ${ARMADILLO_LIBRARIES})
endif ()

#### VTK
find_package(VTK REQUIRED NO_MODULE)
if (VTK_INCLUDE_DIRS)
  list(APPEND VOXEL2TET_EXTERNAL_INCLUDE_DIRS ${VTK_INCLUDE_DIRS})
else ()
  find_path(VTK_HEADERS_DIR vtkSmartPointer.h PATH_SUFFIXES vtk vtk-9.3 vtk-9.2 vtk-9.1 vtk-9.0)
  if (VTK_HEADERS_DIR)
    list(APPEND VOXEL2TET_EXTERNAL_INCLUDE_DIRS ${VTK_HEADERS_DIR})
  endif ()
endif ()
list(APPEND VOXEL2TET_EXTERNAL_LIBS ${VTK_LIBRARIES})

#### TetGen
set(TETGEN_NAMES libtetgen tetgen libtet libtet.so)
find_library(TETGEN_LIB NAMES ${TETGEN_NAMES} DOC "The TetGen library")
if (NOT TETGEN_LIB)
  message(FATAL_ERROR "TetGen library not found")
endif ()
message(STATUS "TetGen library found in ${TETGEN_LIB}")
list(APPEND VOXEL2TET_EXTERNAL_LIBS ${TETGEN_LIB})

find_path(TETGEN_H tetgen.h)
if (NOT TETGEN_H)
  message(FATAL_ERROR "Include file tetgen.h not found")
endif ()
message(STATUS "Include file tetgen.h found in ${TETGEN_H}")
list(APPEND VOXEL2TET_EXTERNAL_INCLUDE_DIRS ${TETGEN_H})

list(REMOVE_DUPLICATES VOXEL2TET_EXTERNAL_INCLUDE_DIRS)
list(REMOVE_DUPLICATES VOXEL2TET_EXTERNAL_LIBS)

add_subdirectory(lib)
add_subdirectory(bin)
