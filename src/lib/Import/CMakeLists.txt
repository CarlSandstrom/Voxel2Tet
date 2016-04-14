set(ImportLibList
	${CMAKE_CURRENT_SOURCE_DIR}/CallbackImporter.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/Importer.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/hdf5DataReader.cpp
)

set(Voxel2TetLibList ${Voxel2TetLibList} ${ImportLibList} PARENT_SCOPE)
