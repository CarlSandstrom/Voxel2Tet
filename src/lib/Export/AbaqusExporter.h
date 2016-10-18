#ifndef ABAQUSEXPORTER_H
#define ABAQUSEXPORTER_H

#include "Exporter.h"

namespace voxel2tet
{

class AbaqusExporter : public Exporter
{
private:
    bool UsePhon;
public:
    AbaqusExporter(std :: vector< TriangleType * > *Triangles, std :: vector< VertexType * > *Vertices, std :: vector< EdgeType * > *Edges, std :: vector< TetType * > *Tets, bool UsePhon=false);
    virtual void WriteSurfaceData(std :: string Filename) {}
    virtual void WriteVolumeData(std :: string Filename);
};

}

#endif // ABAQUSEXPORTER_H
