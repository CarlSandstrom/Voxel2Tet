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
    FC_DIFFERENTSURFACES, FC_ANGLESNOTIMPROVED, FC_TOOLARGEERROR
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
     * @brief Adds an Edge object to EdgeList given vertex IDs
     *
     * AddEdge(std::vector<int>) creates a new Edge object and calls AddEdge(EdgeType *) which in turn checks
     * if the new edge already exists. If so, a reference to the existing EdgeType object is returned, otherwise
     * a reference to the new object is returned.
     *
     * @param VertexIDs Vertex indices
     * @return Reference to EdgeType object
     */
    EdgeType *AddEdge(std :: vector< int >VertexIDs);

    // Adds and edge object to list and update vertices
    EdgeType *AddEdge(EdgeType *e);

    void RemoveEdge(EdgeType *e);
    //TODO: Add function for removing edges

    void RemoveTriangle(TriangleType *t);

    // Adds a triangle using coordinates
    TriangleType *AddTriangle(std :: vector< double >n0, std :: vector< double >n1, std :: vector< double >n2);

    // Adds a triangle using vertex indices
    TriangleType *AddTriangle(std :: vector< int >VertexIDs);

    // Adds a triangle object to the list of triangles
    TriangleType *AddTriangle(TriangleType *NewTriangle);

    TetType *AddTetrahedron(std :: vector< int >VertexIDs);

    TetType *AddTetrahedron(TetType *NewTet);

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

    FC_MESH CheckTrianglePenetration(TriangleType *t1, TriangleType *t2);

    FC_MESH CheckTrianglePenetration(std :: array< VertexType *, 3 >t1, std :: array< VertexType *, 3 >t2, int &sharedvertices);
};
}
#endif // MESHDATA_H
