#ifndef SURFACE_H
#define SURFACE_H

#include <string>

#include "MeshComponents.h"
#include "PhaseEdge.h"

namespace voxel2tet
{

class Surface
{
private:
    Options *Opt;
public:
    std::vector <VertexType*> FixedVertices;

    Surface(int Phase1, int Phase2, Options *Opt);
    int Phases[2];
    std::vector <VertexType*> Vertices;
    std::vector <TriangleType*> Triangles;

    void AddVertex(VertexType* Vertex);
    void AddTriangle(TriangleType* Triangle);

    void Smooth();
};

}
#endif // SURFACE_H
