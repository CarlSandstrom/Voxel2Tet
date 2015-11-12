cmake_minimum_required (VERSION 2.6)
project (Voxel2Tet)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -ffloat-store -std=c++11")

set(linkerlist Voxel2TetLib)

# Look for required libraries

#### hdf5 library

#find_library(HDF5_CXX_LIB hdf5_cpp)
find_package(HDF5)

if (HDF5_FOUND)
  option(WITH_HDF5 "Use HDF5 library" ON)
  include_directories(${HDF5_INCLUDE_DIRS})
  set(hdf5Lib ${hdf5Lib} ${HDF5_LIBRARIES})
  set(linkerlist ${linkerlist} ${hdf5Lib})
  set(linkerlist ${linkerlist} "hdf5_cpp")
  message(STATUS "${linkerlist}")
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

#### Armadillo library

find_package(Armadillo REQUIRED)
if (ARMADILLO_FOUND)
  include_directories(${ARMADILLO_INCLUDE_DIRS})
  set(linkerlist ${linkerlist} ${ARMADILLO_LIBRARIES})
endif (ARMADILLO_FOUND)

#### VTK

find_package(VTK REQUIRED)
if (VTK_FOUND)
  include(${VTK_USE_FILE})
  set(linkerlist ${linkerlist} ${VTK_LIBRARIES})
endif (VTK_FOUND)

#### Project files

include_directories(lib)
add_subdirectory(lib)

include_directories(Examples)
add_subdirectory(Examples)

