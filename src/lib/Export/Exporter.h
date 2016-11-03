#ifndef EXPORTER_H
#define EXPORTER_H

#include <vector>
#include <string>
#include <map>

#include "MeshComponents.h"
#include "MiscFunctions.h"

namespace voxel2tet
{
enum Exporter_FileTypes
{
    FT_VTK, FT_Poly, FT_OFF, FT_OOFEM, FT_SIMPLE, FT_ABAQUS, FT_ABAQUSPHON
};

/**
 * @brief Abstract class for all classes exporting meshes.
 */
class Exporter
{
protected:
    /**
     * List of all vertices used in mesh
     */
    std::vector<VertexType *> UsedVertices;

    /**
     * Updates member UsedVertices such that is contains all used vertices (Mesh.vertices contains all vertices, both used and unused).
     */
    void UpdateUsedVertices();

    /**
     * Maps my material ID to the destination material
     */
    std::map<int, int> MapSelfMaterials;

    /**
     * Update MapSelfMaterial
     */
    void UpdateMaterialsMapping();

    /**
     * Contains highest x, y and z coordinates
     */
    std::array<double, 3> MaxCoords;

    /**
     * Contains lowest x, y and z coordinates
     */
    std::array<double, 3> MinCoords;

    /**
     * Updates MaxCoords and MinCoords members
     */
    void UpdateMinMaxCoordinates();

    /**
     * Contains the nodes in maximum {x, y, z} direction. I.e. first element in array is nodes on the X+ surface etc.
     */
    std::array<std::vector<VertexType *>, 3> MaxNodes;

    /**
     * Contains the nodes in minimum {x, y, z} direction. I.e. first element in array is nodes on the X+ surface etc.
     */
    std::array<std::vector<VertexType *>, 3> MinNodes;

    /**
     * Updates MaxNodes and MinNodes
     */
    void UpdateMinMaxNodes();

    /**
     * Contains the elements with a side in maximum {x, y, z} direction. I.e. has a side in the X+ surface
     */
    std::array<std::vector<TetType *>, 3> MaxElements;

    /**
     * Tells which element side the elements in MaxElements are on the maximum surface
     */
    std::array<std::vector<int>, 3> MaxSide;

    /**
     * Contains the elements with a side in minimum {x, y, z} direction. I.e. has a side in the X+ surface
     */
    std::array<std::vector<TetType *>, 3> MinElements;

    /**
     * Tells which element side the elements in MinElements are on the maximum surface
     */
    std::array<std::vector<int>, 3> MinSide;

    /**
     * Updates MaxElements, MinElements, MaxSide and MinSide
     */
    void UpdateMinMaxElements();

    /**
     * @brief GrainSets is a vector of pairs. The pair, in turn, contains a vector of TetTypes and an associated integer which is the grain ID of those tets.
     * This is used to determin which set the grain belongs to
     */
    std::vector<std::pair<std::vector<TetType *>, int> *> GrainSets;

    /**
     * Updates GrainSets
     */
    void UpdateGrainSets();

    /**
     * @brief TriangleSets is a vector of pairs. The pair, in turn, contains a vector of TriangleTypes and an associated integer which is the interface ID of those triangles.
     * This is used to determin which set the triangles belongs to
     */
    std::vector<std::pair<std::vector<TriangleType *>, int> *> TriangleSets;

    /**
     * Updates TriangleSets
     */
    void UpdateTriangleSets();

    /**
     * @brief Pointer to a list of pointers to the triangles in the mesh. This is typically a pointer to the list of triangles
     * in the MeshData object.
     */
    std::vector<TriangleType *> *Triangles;

    /**
     * @brief Pointer to a list of pointers to the vertices in the mesh. This is typically a pointer to the list of vertices
     * in the MeshData object.
     */
    std::vector<VertexType *> *Vertices;

    /**
     * @brief Pointer to a list of pointers to the edges in the mesh. This is typically a pointer to the list of edges
     * in the MeshData object.
     */
    std::vector<EdgeType *> *Edges;

    /**
     * @brief Pointer to a list of pointers to the tetrahedrons in the mesh. This is typically a pointer to the list of tetrahedrons
     * in the MeshData object.
     */
    std::vector<TetType *> *Tets;
public:
    /**
     * @brief Constructor of the class given lists of neccessary elements.
     * @param Triangles Pointer to a list of triangles
     * @param Vertices Pointer to a list of vertices
     * @param Edges Pointer to a list of edges
     * @param Tets Pointer to a list of tetrahedrons
     */
    Exporter(std::vector<TriangleType *> *Triangles, std::vector<VertexType *> *Vertices,
             std::vector<EdgeType *> *Edges, std::vector<TetType *> *Tets);

    /**
     * @brief Exports a surface mesh. I.e. interfaces between different material IDs. Tetrahedrons will be ignored.
     * @param Filename String. File name of target file.
     */
    virtual void WriteSurfaceData(std::string Filename) = 0;

    /**
     * @brief Exports a volume mesh.
     * @param Filename String. File name of target file.
     */
    virtual void WriteVolumeData(std::string Filename) = 0;
};
}

#endif // EXPORTER_H
