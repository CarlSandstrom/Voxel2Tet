#ifndef EXPORTER_H
#define EXPORTER_H

#include <vector>
#include <string>

//#include "MeshData.h"
#include "MeshComponents.h"

namespace voxel2tet
{

class MeshData;

class Exporter
{
private:
public:
    Exporter();
    Exporter(MeshData *mesh);
    std :: vector <TriangleType*> *Triangles;
    std :: vector <VertexType*> *Vertices;
    std :: vector <EdgeType*> *Edges;
    virtual void WriteData(std::string Filename) {}
};

}

#endif // EXPORTER_H
