set(ExportLibList
	${CMAKE_CURRENT_SOURCE_DIR}/Exporter.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/VTKExporter.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/OFFExporter.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/TetGenExporter.cpp
)

set(Voxel2TetLibList ${Voxel2TetLibList} ${ExportLibList} PARENT_SCOPE)


