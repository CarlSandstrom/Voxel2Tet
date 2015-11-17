#ifndef EXPORTER_H
#define EXPORTER_H

#include <vector>
#include <string>

#include "MeshComponents.h"
#include "MiscFunctions.h"

namespace voxel2tet
{

class Exporter
{
private:
public:
    Exporter();
    Exporter(std :: vector <TriangleType*> *Triangles, std :: vector <VertexType*> *Vertices, std :: vector <EdgeType*> *Edges);
    std :: vector <TriangleType*> *Triangles;
    std :: vector <VertexType*> *Vertices;
    std :: vector <EdgeType*> *Edges;
    virtual void WriteData(std::string Filename) {}
};

}

#endif // EXPORTER_H
