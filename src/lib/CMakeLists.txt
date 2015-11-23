set(Voxel2TetLibList
    Options.cpp
    Importer.cpp
    hdf5DataReader.cpp
    VertexOctreeNode.cpp
    Exporter.cpp
    VTKExporter.cpp
    MeshData.cpp
    MeshComponents.cpp
    MiscFunctions.cpp
    Surface.cpp
    Voxel2Tet.cpp
    PhaseEdge.cpp
    CallbackImporter.cpp
    TetGenExporter.cpp
)

add_library(Voxel2TetLib ${Voxel2TetLibList})
