#ifndef TETGETEXPORTER_H
#define TETGETEXPORTER_H
#include <string>
#include <array>

#include "Exporter.h"
#include "MeshComponents.h"
#include "MiscFunctions.h"
#include "Exporter.h"

namespace voxel2tet
{
// TODO: Do we need this class?
class TetGenExporter : public Exporter
{
public:
    TetGenExporter();
    TetGenExporter(std :: vector< TriangleType * > *Triangles, std :: vector< VertexType * > *Vertices, std :: vector< EdgeType * > *Edges, std :: vector< TetType * > *Tets);
    virtual void WriteSurfaceData(std :: string Filename);
    virtual void WriteVolumeData(std :: string Filename) {}
};
}
#endif // TETGETEXPORTER_H
