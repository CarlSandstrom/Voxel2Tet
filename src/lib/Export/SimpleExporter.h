#ifndef SIMPLEEXPORTER_H
#define SIMPLEEXPORTER_H

#include "Exporter.h"

namespace voxel2tet
{

class SimpleExporter: public Exporter
{
public:
    SimpleExporter(std :: vector< TriangleType * > *Triangles, std :: vector< VertexType * > *Vertices, std :: vector< EdgeType * > *Edges, std :: vector< TetType * > *Tets);
    virtual void WriteSurfaceData(std :: string Filename);
    virtual void WriteVolumeData(std :: string Filename) {}

};

}
#endif // SIMPLEEXPORTER_H
