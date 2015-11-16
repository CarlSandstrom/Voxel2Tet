#include "Exporter.h"
#include "MiscFunctions.h"

namespace voxel2tet
{

Exporter::Exporter()
{

}

Exporter::Exporter(MeshData *mesh)
{
    log("Create exporter for MeshData@%p\n", mesh);
    this->Vertices = std::addressof( mesh->Vertices );
    this->Edges = std::addressof( mesh->Edges );
    this->Triangles = std::addressof( mesh->Triangles );
}

}
