cmake_minimum_required (VERSION 2.6)
project (Voxel2Tet C CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")# -Wno-comment")
set(linkerlist Voxel2TetLib)

# Look for required libraries

#### hdf5 library

FIND_PACKAGE(HDF5 COMPONENTS CXX HL NO_MODULE REQUIRED static)
message(STATUS "C libraries: ${HDF5_LIBRARIES}")

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
endif (ARMADILLO_FOUND)

#### VTK

find_package(VTK REQUIRED NO_MODULE)
if (VTK_FOUND)
  include(${VTK_USE_FILE})
  set(linkerlist ${linkerlist} ${VTK_LIBRARIES})
endif (VTK_FOUND)

#### TetGen

find_path(TETGEN_H tetgen.h)
MESSAGE( STATUS "Include file tetgen.h found in ${TETGEN_H}" )
include_directories (${TETGEN_H})

find_library(TETGEN_LIB libtet.so)
MESSAGE( STATUS "Tetgen library found in ${TETGEN_LIB}")
set(linkerlist ${linkerlist} ${TETGEN_LIB})

#### Project files

include_directories(lib lib/Export lib/Import lib/Mesh lib/General)
add_subdirectory(lib)

message (STATUS "My files: ${lib}")

get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message (STATUS "inc_dirs = ${inc_dirs}")

include_directories(bin)
add_subdirectory(bin)

message(STATUS "LinkerList: ${linkerlist}")