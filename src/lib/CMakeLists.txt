set(Voxel2TetLibList
    Options.cpp
    Importer.cpp
    hdf5DataReader.cpp
    VertexOctreeNode.cpp
    Exporter.cpp
    VTKExporter.cpp
)

add_library(Voxel2TetLib ${Voxel2TetLibList})
