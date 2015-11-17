#ifndef MESHDATA_H
#define MESHDATA_H

#include<vector>
#include<string>

#include "MeshComponents.h"
#include "VertexOctreeNode.h"
#include "VTKExporter.h"

namespace voxel2tet
{

class MeshData
{
private:
    BoundingBoxType BoundingBox;
public:
    std :: vector <TriangleType*> Triangles;
    std :: vector <VertexType*> Vertices;
    std :: vector <EdgeType*> Edges;

    VertexOctreeNode *VertexOctreeRoot;

    // Either with all vertices or with a bounding box
    MeshData(BoundingBoxType BoundingBox);

    // Export mesh to VTK
    void ExportVTK(std::string FileName);

    // Adds an edge using vertex indices
    EdgeType *AddEdge(std::vector<int> VertexIDs);

    // Adds a triangle using coordinates
    TriangleType *AddTriangle(std::vector<double> n0, std::vector<double> n1, std::vector<double> n2);

    // Adds a triangle using vertex indices
    TriangleType *AddTriangle(std::vector<int> VertexIDs);

};

}
#endif // MESHDATA_H
