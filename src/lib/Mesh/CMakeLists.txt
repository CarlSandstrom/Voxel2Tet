set(MeshLibList
	${CMAKE_CURRENT_SOURCE_DIR}/MeshData.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/MeshManipulations.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/MeshComponents.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/MeshGenerator3D.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/TetGenCaller.cpp
)

set(Voxel2TetLibList ${Voxel2TetLibList} ${MeshLibList} PARENT_SCOPE)
