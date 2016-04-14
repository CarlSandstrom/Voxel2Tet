#ifndef OFFEXPORTER_H
#define OFFEXPORTER_H

#include <string>
#include <array>

#include "Exporter.h"
#include "MeshComponents.h"
#include "MiscFunctions.h"
#include "Exporter.h"

namespace voxel2tet
{

class OFFExporter: public Exporter
{
public:
    OFFExporter();
    OFFExporter(std :: vector <TriangleType*> *Triangles, std :: vector <VertexType*> *Vertices, std :: vector <EdgeType*> *Edges);
    virtual void WriteData(std::string Filename);
};

}

#endif // OFFEXPORTER_H
