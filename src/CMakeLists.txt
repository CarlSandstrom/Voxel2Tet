cmake_minimum_required (VERSION 2.6)
project (Voxel2Tet CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -Wno-comment")
set(linkerlist Voxel2TetLib)

# Look for required libraries

#### hdf5 library

find_package(HDF5 COMPONENTS CXX)

if (HDF5_FOUND)
  include_directories(${HDF5_INCLUDE_DIR})
  set(hdf5Lib ${hdf5Lib} ${HDF5_LIBRARIES})
  set(linkerlist ${linkerlist} ${hdf5Lib})
  message(STATUS "HDF5 root: ${HDF5_ROOT}")
  message(STATUS "Includedir: ${HDF5_INCLUDE_DIR}")
else (HDF5_FOUND)
  message(ERROR "HDF5 library not found")
endif (HDF5_FOUND)

#### ARMADILLO

find_package(Armadillo REQUIRED)
if (ARMADILLO_FOUND)
  include_directories(${ARMADILLO_INCLUDE_DIRS})
  set(linkerlist ${linkerlist} ${ARMADILLO_LIBRARIES})
else (ARMADILLO_FOUND)  
  message(ERROR "Armadillo library not found")
endif (ARMADILLO_FOUND)

#### VTK

find_package(VTK REQUIRED NO_MODULE)
if (VTK_FOUND)
  include(${VTK_USE_FILE})
  set(linkerlist ${linkerlist} ${VTK_LIBRARIES})
else (VTK_FOUND)
  message(ERROR "VTK not found")
endif (VTK_FOUND)

#### TetGen

set(TETGEN_NAMES ${TETGEN_NAMES} libtetgen tetgen libtet libtet.so)
find_library(TETGEN_LIB NAMES ${TETGEN_NAMES}
            DOC "The TetGen library")
MESSAGE( STATUS "Tetgen library found in ${TETGEN_LIB}")
set(linkerlist ${linkerlist} ${TETGEN_LIB})

find_path(TETGEN_H tetgen.h)
MESSAGE( STATUS "Include file tetgen.h found in ${TETGEN_H}" )
include_directories (${TETGEN_H})

#find_library(TETGEN_LIBRARY libtet)
#set (TETGEN_LIB ${TETGEN_LIB}/libtet.a)
#set(TETGEN_LIB "/usr/lib/x86_64-linux-gnu/libtet.so")

#### Project files

include_directories(lib lib/Export lib/Import lib/Mesh lib/General)
add_subdirectory(lib)

get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
#message (STATUS "inc_dirs = ${inc_dirs}")

include_directories(bin)
add_subdirectory(bin)

#message(STATUS "LinkerList: ${linkerlist}")
