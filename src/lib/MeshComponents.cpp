#include<algorithm>

#include "MeshComponents.h"

namespace voxel2tet
{

// VertexType

void VertexType :: AddTriangle(TriangleType *Triangle)
{
    for (TriangleType *T: this->Triangles) {
        if (T == Triangle) {
            return;
        }
    }
    this->Triangles.push_back(Triangle);
}

void VertexType :: RemoveTriangle(TriangleType *Triangle)
{
    this->Triangles.erase(std::remove(this->Triangles.begin(), this->Triangles.end(), Triangle), this->Triangles.end());
}

void VertexType :: AddEdge(EdgeType *Edge)
{
    for (EdgeType *E: this->Edges) {
        if (E == Edge) {
            return;
        }
    }
    this->Edges.push_back(Edge);
}

void VertexType :: RemoveEdge(EdgeType *Edge)
{
    this->Edges.erase(std::remove(this->Edges.begin(), this->Edges.end(), Edge), this->Edges.end());
}



}
