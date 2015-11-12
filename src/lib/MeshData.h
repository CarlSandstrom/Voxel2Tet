#ifndef MESHDATA_H
#define MESHDATA_H

#include<vector>
#include<string>

#include "MeshComponents.h"

namespace voxel2tet
{

class MeshData
{
private:
    std :: vector <TriangleType*> Triangles;
    std :: vector <VertexType*> Vertices;
    std :: vector <EdgeType*> Edges;
public:
    MeshData();

    // Export mesh to VTK
    void ExportVTK(std::string FileName);

    // Adds a triangle using coordinates
    int AddTriangle(double x, double y, double z);

};

}
#endif // MESHDATA_H
