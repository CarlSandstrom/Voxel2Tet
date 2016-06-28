#ifndef OOFEMEXPORTER_H
#define OOFEMEXPORTER_H

#include "Exporter.h"

namespace voxel2tet
{
class OOFEMExporter : public Exporter
{
public:
    OOFEMExporter();
    OOFEMExporter(std :: vector< TriangleType * > *Triangles, std :: vector< VertexType * > *Vertices, std :: vector< EdgeType * > *Edges, std :: vector< TetType * > *Tets);
    virtual void WriteSurfaceData(std :: string Filename) {}
    virtual void WriteVolumeData(std :: string Filename);
};
}

#endif // OOFEMEXPORTER_H
