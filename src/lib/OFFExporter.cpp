#include "OFFExporter.h"
#include <iostream>
#include <fstream>
#include <iomanip>

namespace voxel2tet
{

OFFExporter::OFFExporter()
{

}

OFFExporter::OFFExporter(std :: vector <TriangleType*> *Triangles, std :: vector <VertexType*> *Vertices, std :: vector <EdgeType*> *Edges):
    Exporter (Triangles, Vertices, Edges)
{
    LOG("Create TetGenExporter object\n",0);
}

void OFFExporter::WriteData(std::string Filename)
{
    STATUS("Write .off file %s\n", Filename.c_str());
    std::ofstream OFFFile;
    OFFFile.open(Filename);

    // Write header
    OFFFile << "OFF\n" << this->Vertices->size() << "\t" << this->Triangles->size() << "\t0\n";

    // Write vertices
    for (auto v: *this->Vertices) {
        OFFFile << std::setiosflags(std::ios::fixed) << std::setprecision(15) << v->c[0]<< " " << v->c[1] << " " << v->c[2] << "\n";
    }

    // Write facets
    int j=0;

    for (auto t: *this->Triangles) {
        OFFFile << 3;
        std::array<int,3> VertexIDs;
        for (int i=0; i<3; i++) {
            VertexIDs[i] = t->Vertices[i]->ID;
            OFFFile << " " << VertexIDs[i];
        }
        OFFFile << "\t#" << j << "\n";
        j++;
    }


}

}
