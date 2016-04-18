#ifndef MESHDATA_H
#define MESHDATA_H

#include<vector>
#include<string>
#include<algorithm>

#include "MeshComponents.h"
#include "VertexOctreeNode.h"
#include "VTKExporter.h"

namespace voxel2tet
{

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
};

}
#endif // MESHDATA_H
