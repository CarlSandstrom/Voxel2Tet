set(MeshLibList
	${CMAKE_CURRENT_SOURCE_DIR}/MeshData.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/MeshManipulations.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/MeshComponents.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/MeshGenerator3D.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/TetGenCaller.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/Smoother.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/SpringSmoother.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/SpringSmootherPenalty.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/LaplaceSmoother.cpp
)

set(Voxel2TetLibList ${Voxel2TetLibList} ${MeshLibList} PARENT_SCOPE)
