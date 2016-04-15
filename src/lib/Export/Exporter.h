#ifndef EXPORTER_H
#define EXPORTER_H

#include <vector>
#include <string>

#include "MeshComponents.h"
#include "MiscFunctions.h"

namespace voxel2tet
{

enum Exporter_FileTypes {FT_VTK, FT_Poly, FT_OFF};

class Exporter
{
private:
public:
    Exporter();
    Exporter(std :: vector <TriangleType*> *Triangles, std :: vector <VertexType*> *Vertices, std :: vector <EdgeType*> *Edges);
    std :: vector <TriangleType*> *Triangles;
    std :: vector <VertexType*> *Vertices;
    std :: vector <EdgeType*> *Edges;
    virtual void WriteData(std::string Filename) = 0;
};

}

#endif // EXPORTER_H
