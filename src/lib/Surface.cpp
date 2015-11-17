#include "Surface.h"

namespace voxel2tet
{

Surface::Surface(int Phase1, int Phase2)
{
    this->Phases[0] = Phase1;
    this->Phases[1] = Phase2;
}

void Surface::AddVertex(VertexType *Vertex)
{
    for (auto v: this->Vertices) {
        if (v==Vertex) {
            return;
        }
    }
    this->Vertices.push_back(Vertex);
}

void Surface::AddTriangle(TriangleType *Triangle)
{
    for (auto t: this->Triangles) {
        if (t==Triangle) {
            return;
        }
    }
    this->Triangles.push_back(Triangle);

}



}
