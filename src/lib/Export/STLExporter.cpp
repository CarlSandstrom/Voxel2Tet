#include "STLExporter.h"

namespace voxel2tet
{

STLExporter::STLExporter(std::vector<TriangleType *> *Triangles, std::vector<VertexType *> *Vertices,
                         std::vector<EdgeType *> *Edges, std::vector<TetType *> *Tets) :
        Exporter(Triangles, Vertices, Edges, Tets)
{
    LOG("Create OFFExporter object\n", 0);
}

void STLExporter::WriteSurfaceData(std::string Filename, int VolumeID)
{
    STATUS("Write .stl file %s\n", Filename.c_str());
    std::ofstream STLFile;
    STLFile.open(Filename);
    STLFile << "solid smoothed\n";

    for (auto t : *this->Triangles) {

        if ((VolumeID==-1) | (t->NegNormalMatID==VolumeID) | (t->PosNormalMatID==VolumeID)){
            std::array<double, 3> Normal = t->GiveNormal();
            STLFile << "facet normal " << std::scientific << Normal[0] << " " << Normal[1] << " " << Normal[2] << "\n";
            STLFile << "\touter loop\n";
            for(int i=0; i<3; i++ ) {
                STLFile << "\t\tvertex " << t->Vertices[i]->get_c(0) << " " << t->Vertices[i]->get_c(1) << " " <<t->Vertices[i]->get_c(2) << "\n";
            }
            STLFile << "\tendloop\nendfacet\n";
        }
    }


    STLFile << "endsolid smoothed";

}

}
