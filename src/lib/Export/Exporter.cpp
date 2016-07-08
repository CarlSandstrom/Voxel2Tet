#include "Exporter.h"

namespace voxel2tet
{

Exporter :: Exporter(std :: vector< TriangleType * > *Triangles, std :: vector< VertexType * > *Vertices, std :: vector< EdgeType * > *Edges, std :: vector< TetType * > *Tets)
{
    LOG("Create exporter for MeshData@%p\n", Triangles);
    this->Vertices = Vertices;
    this->Edges = Edges;
    this->Triangles = Triangles;
    this->Tets = Tets;
}
}
