#ifndef MESHDATA_H
#define MESHDATA_H

#include<vector>
#include<string>
#include<algorithm>

#include "MeshComponents.h"
#include "VertexOctreeNode.h"
#include "VTKExporter.h"
#include "TriTriIntersect.h"

namespace voxel2tet
{

enum FC_MESH {FC_OK, FC_FIXEDVERTEX, FC_NORMAL, FC_CHORD, FC_SMALLAREA, FC_AREACHANGETOOLARGE, FC_TOOMANYTRIANGLES,
              FC_WORSEMINANGLE, FC_VERTICESONDIFFERENTSHAPES, FC_TRIANGLESINTERSECT, FC_DUPLICATETRIANGLE, FC_INVALIDEDGE,
             FC_DIFFERENTSURFACES, FC_ANGLESNOTIMPROVED};


/**
 * @brief The MeshData class supplies information and methods for adding vertices and triangles.
 */
class MeshData
{
private:
    int TriangleCounter;

public:
    void DoSanityCheck();
    BoundingBoxType BoundingBox;
    std :: vector <TriangleType*> Triangles;
    std :: vector <VertexType*> Vertices;
    std :: vector <EdgeType*> Edges;
    std :: vector <TetType*> Tets;

    // TODO: Add a AddVertex function
    VertexOctreeNode *VertexOctreeRoot;

    // Either with all vertices or with a bounding box
    MeshData(BoundingBoxType BoundingBox);

    ~MeshData();

    // Export surface mesh
    void ExportSurface(std::string FileName, Exporter_FileTypes FileType);

    // Export volume mesh
    void ExportVolume(std::string FileName, Exporter_FileTypes FileType);

    // Adds an edge using vertex indices
    EdgeType *AddEdge(std::vector<int> VertexIDs);

    // Adds and edge object to list and update vertices
    EdgeType *AddEdge(EdgeType *e);

    void RemoveEdge(EdgeType *e);
    //TODO: Add function for removing edges

    void RemoveTriangle(TriangleType *t);

    // Adds a triangle using coordinates
    TriangleType *AddTriangle(std::vector<double> n0, std::vector<double> n1, std::vector<double> n2);

    // Adds a triangle using vertex indices
    TriangleType *AddTriangle(std::vector<int> VertexIDs);

    // Adds a triangle object to the list of triangles
    TriangleType *AddTriangle(TriangleType *NewTriangle);

    TetType *AddTetrahedron(std::vector<int> VertexIDs);

    TetType *AddTetrahedron(TetType *NewTet);

    void RemoveTetragedron(TetType *t);

    /**
     * @brief CheckSameOrientation tells if two neighbouring triangles are oriented in the same
     * way by comparing the order of the vertices on the shared edge
     * @param t1
     * @param t2
     * @return True if oriented in the same way, False otherwise
     */
    bool CheckSameOrientation(TriangleType *t1, TriangleType *t2);

    FC_MESH CheckTrianglePenetration(TriangleType *t1, TriangleType *t2);

    FC_MESH CheckTrianglePenetration(std::array<VertexType *, 3> t1, std::array<VertexType *, 3> t2, int &sharedvertices);
};

}
#endif // MESHDATA_H
