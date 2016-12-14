cmake_minimum_required (VERSION 2.6)
project (Voxel2Tet)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")
set(linkerlist Voxel2TetLib)

# Look for required libraries

#### hdf5 library

#set (HDF5_USE_STATIC_LIBRARIES ON)
find_package(HDF5 COMPONENTS CXX)

if (HDF5_FOUND)
  option(WITH_HDF5 "Use HDF5 library" ON)
  include_directories(${HDF5_INCLUDE_DIRS})
  set(hdf5Lib ${hdf5Lib} ${HDF5_CXX_LIBRARIES})
  set(linkerlist ${linkerlist} ${hdf5Lib})
  #set(linkerlist ${linkerlist} "hdf5_cpp")
  message(STATUS "Libraries: ${hdf5Lib}")
  message(STATUS "LinkeList: ${linkerlist}")
else (HDF5_FOUND)
  message(FATAL_ERROR "HDF5 library is required to compile this project")
  option(WITH_HDF5 "Use HDF5 library" OFF)
endif (HDF5_FOUND)

if (WITH_HDF5)
  if (NOT HDF5_FOUND)
      message(SEND_ERROR "Cannot find HDF5 library")
  endif (NOT HDF5_FOUND)
  add_definitions(-DWITH_HDF5=1)
  add_definitions(-DH5_NO_DEPRECATED_SYMBOLS)
endif (WITH_HDF5)

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

get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message("inc_dirs = ${inc_dirs}")

include_directories(bin)
add_subdirectory(bin)

