#ifndef EXPORTER_H
#define EXPORTER_H

#include <vector>
#include <string>

#include "MeshComponents.h"
#include "MiscFunctions.h"

namespace voxel2tet
{
enum Exporter_FileTypes { FT_VTK, FT_Poly, FT_OFF, FT_OOFEM, FT_SIMPLE };

/**
 * @brief Abstract class for all classes exporting meshes.
 */
class Exporter
{
protected:
    /**
     * @brief Pointer to a list of pointers to the triangles in the mesh. This is typically a pointer to the list of triangles
     * in the MeshData object.
     */
    std :: vector< TriangleType * > *Triangles;

    /**
     * @brief Pointer to a list of pointers to the vertices in the mesh. This is typically a pointer to the list of vertices
     * in the MeshData object.
     */
    std :: vector< VertexType * > *Vertices;

    /**
     * @brief Pointer to a list of pointers to the edges in the mesh. This is typically a pointer to the list of edges
     * in the MeshData object.
     */
    std :: vector< EdgeType * > *Edges;

    /**
     * @brief Pointer to a list of pointers to the tetrahedrons in the mesh. This is typically a pointer to the list of tetrahedrons
     * in the MeshData object.
     */
    std :: vector< TetType * > *Tets;
public:
    /**
     * @brief Constructor of the class given lists of neccessary elements.
     * @param Triangles Pointer to a list of triangles
     * @param Vertices Pointer to a list of vertices
     * @param Edges Pointer to a list of edges
     * @param Tets Pointer to a list of tetrahedrons
     */
    Exporter(std :: vector< TriangleType * > *Triangles, std :: vector< VertexType * > *Vertices, std :: vector< EdgeType * > *Edges, std :: vector< TetType * > *Tets);

    /**
     * @brief Exports a surface mesh. I.e. interfaces between different material IDs. Tetrahedrons will be ignored.
     * @param Filename String. File name of target file.
     */
    virtual void WriteSurfaceData(std :: string Filename) = 0;

    /**
     * @brief Exports a volume mesh.
     * @param Filename String. File name of target file.
     */
    virtual void WriteVolumeData(std :: string Filename) = 0;
};
}

#endif // EXPORTER_H
