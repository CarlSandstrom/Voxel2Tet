#ifndef TETGENCALLER_H
#define TETGENCALLER_H

#include <algorithm>

// This is for telling tetgen.h to use library definitions
#define TETLIBRARY = 1

#include "tetgen.h"
#include "MeshGenerator3D.h"

namespace voxel2tet
{
/**
 * @brief The TetGenCaller class implement functionality to call the TetGen library for producing tetrahedral meshes.
 */
class TetGenCaller : public MeshGenerator3D
{
private:
    // Mapping for vertices. Key is Vertex.ID and value is the condensed ID (i.e. the ID where all unused vertices are removed)
    std::map<VertexType *, int> VertexMapFromID;
    std::map<int, VertexType *> VertexMapFromTetgen;

    void UpdateVertexMapping();

    void CopyMeshFromSelf(tetgenio *io);

    MeshData *CopyTetMesh(tetgenio *io);

public:

    /**
     * @brief Constructor
     */
    TetGenCaller();

    /**
     * @brief Create mesh and return pointer to MeshData object containing the new mesh
     * @return Pointer to MeshData object
     */
    virtual MeshData *Execute();

    /**
     * @brief Check if mesh is valid
     */
    virtual void TestMesh();
};
}
#endif // MESHGENERATOR3D_H
