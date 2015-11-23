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

class TetGenExporter: public Exporter
{
public:
    TetGenExporter();
    TetGenExporter(std :: vector <TriangleType*> *Triangles, std :: vector <VertexType*> *Vertices, std :: vector <EdgeType*> *Edges);
    virtual void WriteData(std::string Filename);
};

}
#endif // TETGETEXPORTER_H
