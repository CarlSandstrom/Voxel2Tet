set(Voxel2TetLibList
    Surface.cpp
    Voxel2Tet.cpp
    PhaseEdge.cpp
    Volume.cpp
)

set(subdirs Export General Import Mesh)

foreach(subdir ${subdirs})
	set(thisdir "${subdir}")
	#message(STATUS "Process directory ${thisdir}")
	add_subdirectory(${subdir})
	#message(STATUS "Add include path ${thisdir}")
	include_directories(${thisdir})
	#message(STATUS "LIB: Current list: ${Voxel2TetLibList}")
endforeach(subdir)

add_library(Voxel2TetLib ${Voxel2TetLibList})
