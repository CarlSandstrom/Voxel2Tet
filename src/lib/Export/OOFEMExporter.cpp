#include "OOFEMExporter.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <map>
#include <vector>

namespace voxel2tet
{

OOFEMExporter::OOFEMExporter()
{
    LOG("Create TetGenExporter object\n",0);
}

OOFEMExporter::OOFEMExporter(std :: vector <TriangleType*> *Triangles, std :: vector <VertexType*> *Vertices, std :: vector <EdgeType*> *Edges, std :: vector <TetType*> *Tets):
    Exporter (Triangles, Vertices, Edges, Tets)
{
    LOG("Create TetGenExporter object\n",0);
}


void OOFEMExporter::WriteVolumeData(std::string Filename)
{
    STATUS("Write .in (oofem) file %s\n", Filename.c_str());
    std::ofstream OOFEMFile;
    OOFEMFile.open(Filename);

    // Prepare information
    std::vector<VertexType *> UsedVertices;
    for (TetType *t: *this->Tets) {
        for (VertexType *v: t->Vertices) {
            UsedVertices.push_back(v);
        }
    }

    std::sort(UsedVertices.begin(), UsedVertices.end());
    UsedVertices.erase( std::unique(UsedVertices.begin(), UsedVertices.end()), UsedVertices.end());

    int i=0;
    for (VertexType *v: UsedVertices) {
        v->tag = i;
        i++;
    }

    std::map<int, int> Self2OofemMaterials;

    for (TetType *t: *this->Tets) {
        if (Self2OofemMaterials.find(t->MaterialID) == Self2OofemMaterials.end()) {
            Self2OofemMaterials[t->MaterialID] = Self2OofemMaterials.size();
        }
    }

    // Write header

    OOFEMFile << Filename << ".out\n";
    OOFEMFile << "Exported from Voxel2Tet\n";
    OOFEMFile << "staticstructural nsteps 1 deltaT 1 lstype 3 smtype 7 rtolf 1.e-12 nmodules 1 maxiter 50\n";
    OOFEMFile << "vtkxml tstep_step 1 domain_all\n";
    OOFEMFile << "domain 3d\n";
    OOFEMFile << "OutputManager tstep_all dofman_all element_all\n";
    OOFEMFile << "ndofman " << UsedVertices.size() << " nelem " << this->Tets->size() << " ncrosssect 1 nmat " \
              << Self2OofemMaterials.size() << " nbc 3 nic 1 nltf 1 nset 3 nxfemman 0\n";

    // Write vertices

    for (size_t i = 0; i<this->Vertices->size(); i++) {
        OOFEMFile << "node " << i+1 << "\tcoords 3 " << UsedVertices.at(i)->get_c(0) << "\t" \
                  << UsedVertices.at(i)->get_c(1) << "\t" << UsedVertices.at(i)->get_c(2) << "\n";
    }

    // Write elements

    for (size_t i = 0; i<this->Tets->size(); i++) {
        TetType *t = this->Tets->at(i);
        OOFEMFile << "ltrspace " << i+1 << "\tnodes 4\t" << t->Vertices[0]->tag \
                  << "\t" << t->Vertices[1]->tag << "\t" << t->Vertices[2]->tag << "\n";
    }

}

}
