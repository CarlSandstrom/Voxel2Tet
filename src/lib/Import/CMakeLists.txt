set(ImportLibList
	${CMAKE_CURRENT_SOURCE_DIR}/CallbackImporter.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/Importer.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Dream3DDataReader.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/VTKStructuredReader.cpp
)

set(Voxel2TetLibList ${Voxel2TetLibList} ${ImportLibList} PARENT_SCOPE)
