#ifndef MESHDATA_H
#define MESHDATA_H

#include <vector>
#include <string>
#include <algorithm>

#include "MeshComponents.h"
#include "VertexOctreeNode.h"
#include "VTKExporter.h"
#include "SimpleExporter.h"
#include "TriTriIntersect.h"

namespace voxel2tet
{


enum FC_MESH {
    FC_OK, FC_FIXEDVERTEX, FC_NORMAL, FC_CHORD, FC_SMALLAREA, FC_AREACHANGETOOLARGE, FC_TOOMANYTRIANGLES,
    FC_WORSEMINANGLE, FC_VERTICESONDIFFERENTSHAPES, FC_TRIANGLESINTERSECT, FC_DUPLICATETRIANGLE, FC_INVALIDEDGE,
    FC_DIFFERENTSURFACES, FC_ANGLESNOTIMPROVED, FC_TOOLARGEERROR, FC_SMALLANGLE
};


/**
 * @brief The MeshData class supplies information and methods for handling mesh data.
 */
class MeshData
{
private:
    int TriangleCounter;
    int EdgeCounter;

public:

    /**
     * @brief Perform sanity check. Mostly for debugging purposes to make sure that the connectivity of the mesh is kept in order.
     */
    void DoSanityCheck();

    /**
     * @brief Bounding box of mesh
     */
    BoundingBoxType BoundingBox;

    /**
     * @brief List of all triangles in the mesh
     */
    std :: vector< TriangleType * >Triangles;

    /**
     * @brief List of all vertices in the mesh
     */
    std :: vector< VertexType * >Vertices;

    /**
     * @brief List of all edges in the mesh
     */
    std :: vector< EdgeType * >Edges;

    /**
     * @brief List of all tetrahedrons in the mesh
     */
    std :: vector< TetType * >Tets;

    /**
     * @brief VertexOctreeNode root object for all vertices. Uses octree algorithm for performance.
     */
    VertexOctreeNode *VertexOctreeRoot;

    /**
     * @brief Constructor
     * @param BoundingBox Bounding box of mesh
     */
    MeshData(BoundingBoxType BoundingBox);

    ~MeshData();

    /**
     * @brief Exports all Surface objects to file of preferred format
     * @param FileName Name of file
     * @param FileType Type of file
     */
    void ExportSurface(std :: string FileName, Exporter_FileTypes FileType);


    /**
     * @brief Exports all Volume objects to a file of preferred format
     * @param FileName Name of file
     * @param FileType Type of file
     */
    void ExportVolume(std :: string FileName, Exporter_FileTypes FileType);

    /**
     * @brief Adds an Edge object to Edges given vertex IDs
     *
     * AddEdge(std::vector<int>) creates a new Edge object and calls AddEdge(EdgeType *) which in turn checks
     * if the new edge already exists. If so, a reference to the existing EdgeType object is returned, otherwise
     * a reference to the new object is returned.
     *
     * @param VertexIDs Vertex indices
     * @return Reference to EdgeType object
     */
    EdgeType *AddEdge(std :: array< int, 2 >VertexIDs);

    // Adds and edge object to list and update vertices
    /**
     * @brief Adds an EdgeType object to Edges
     *
     * @param e Pointer EdgeType object to add
     * @return If the edge is already in the list, return a pointer to that edge, if not, return a pointer to the new edge.
     */
    EdgeType *AddEdge(EdgeType *e);

    /**
     * @brief Remove EdgeType object from Edges
     * @param e Pointer to EdgeType object to remove
     */
    void RemoveEdge(EdgeType *e);

    /**
     * @brief Removes Triangle from Triangles list
     * @param t Pointer to triangle object to remove
     */
    void RemoveTriangle(TriangleType *t);

    // Adds a triangle using coordinates
    /**
     * @brief Adds Triangle a TriangleType object to the Triangles list given coordinates
     * @param v0 Coordinates of vertex 0
     * @param v1 Coordinates of vertex 1
     * @param v2 Coordinates of vertex 2
     * @return If triangle already exists, a pointer to that object is returned. If not, a pointer to the newly created object is returned.
     */
    TriangleType *AddTriangle(std :: array< double, 3 >v0, std :: array< double, 3 >v1, std :: array< double, 3 >v2);

    /**
     * @brief Adds a TriangleType object to Triangles list given vertex IDs
     * @param VertexIDs Array of vertex IDs
     * @return If triangle already exists, a pointer to that object is returned. If not, a pointer to the newly created object is returned.
     */
    TriangleType *AddTriangle(std :: array< int, 3 >VertexIDs);

    /**
     * @brief Adds a TriangleType object to Triangles list given a pointer to an TriangleType object
     * @param NewTriangle Pointer to new triangle
     * @return If triangle already exists, a pointer to that object is returned. If not, a pointer to the newly created object is returned.
     */
    TriangleType *AddTriangle(TriangleType *NewTriangle);

    /**
     * @brief Adds a TetType object to the Tets list given vertex IDs
     * @param VertexIDs Array of vertex IDs
     * @return Returns pointer to the added TetType object
     */
    TetType *AddTetrahedron(std :: array< int, 4 >VertexIDs);

    /**
     * @brief Adds a TetType object to the Tets list given a pointer to a TetType object
     * @param NewTet Pointer to TetType object
     * @return Pointer to the added TetType object
     */
    TetType *AddTetrahedron(TetType *NewTet);

    /**
     * @brief Removes a TetType object from Tets list. NOT YET IMPLEMENTED.
     * @param t Pointer to TetType object
     */
    void RemoveTetragedron(TetType *t);

    /**
     * @brief GetTrianglesAround gives a vector containing all triangles with at least one vertex withing the sphere created by the coordinate c and distance r.
     * @param c [in] Center of sphere
     * @param r [in] Radius of sphere
     * @return List of triangles
     */
    std :: vector< TriangleType * >GetTrianglesAround(std :: array< double, 3 >c, double r);

    /**
     * @brief CheckSameOrientation tells if two neighbouring triangles are oriented in the same
     * way by comparing the order of the vertices on the shared edge
     * @param t1
     * @param t2
     * @return True if oriented in the same way, False otherwise
     */
    bool CheckSameOrientation(TriangleType *t1, TriangleType *t2);

    /**
     * @brief Checks if two triangle penetrate given two TriangleType objects.
     *
     * @param t1 Pointer to TriangleType object
     * @param t2 Pointer to TriangleType object
     * @return FC_MESH type determining the status of the intersection of the two triangles
     */
    FC_MESH CheckTrianglePenetration(TriangleType *t1, TriangleType *t2);

    /**
     * @brief Check if two triangle penetrate given the vertex IDs
     *
     * The procedure uses the intersection algorithms by Tomas Möller.
     *
     * @param t1 Array of three vertex IDs
     * @param t2 Array of three vertex IDs
     * @param sharedvertices Number of shared vertices
     * @return FC_MESH type determining the status of the intersection of the two triangles
     */
    FC_MESH CheckTrianglePenetration(std :: array< VertexType *, 3 >t1, std :: array< VertexType *, 3 >t2, int &sharedvertices);
};
}
#endif // MESHDATA_H
