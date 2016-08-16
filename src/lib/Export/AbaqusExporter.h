#ifndef ABAQUSEXPORTER_H
#define ABAQUSEXPORTER_H

#include "Exporter.h"

namespace voxel2tet
{

class AbaqusExporter : public Exporter
{
public:
    AbaqusExporter(std :: vector< TriangleType * > *Triangles, std :: vector< VertexType * > *Vertices, std :: vector< EdgeType * > *Edges, std :: vector< TetType * > *Tets);
    virtual void WriteSurfaceData(std :: string Filename) {}
    virtual void WriteVolumeData(std :: string Filename);
};

}

#endif // ABAQUSEXPORTER_H
