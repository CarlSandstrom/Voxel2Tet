#include "Exporter.h"

namespace voxel2tet
{

Exporter::Exporter()
{

}

Exporter::Exporter(std::vector<TriangleType *> *Triangles, std::vector<VertexType *> *Vertices, std::vector<EdgeType *> *Edges)
{
    LOG("Create exporter for MeshData@%p\n", Triangles);
    this->Vertices = Vertices;
    this->Edges = Edges;
    this->Triangles = Triangles;
}

}
