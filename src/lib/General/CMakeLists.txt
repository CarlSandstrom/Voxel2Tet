set(GeneralLibList
	${CMAKE_CURRENT_SOURCE_DIR}/MiscFunctions.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/Options.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/VertexOctreeNode.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/TriTriIntersect.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/TimeStamp.cpp
)

set(Voxel2TetLibList ${Voxel2TetLibList} ${GeneralLibList} PARENT_SCOPE)

