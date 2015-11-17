#ifndef SURFACE_H
#define SURFACE_H

#include <string>

#include "MeshComponents.h"

namespace voxel2tet
{

class Surface
{
public:
    Surface(int Phase1, int Phase2);
    int Phases[2];
    std::vector <VertexType*> Vertices;
    std::vector <TriangleType*> Triangles;

    void AddVertex(VertexType* Vertex);
    void AddTriangle(TriangleType* Triangle);
};

}
#endif // SURFACE_H
