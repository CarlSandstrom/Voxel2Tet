#include "AbaqusExporter.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <map>
#include <vector>

namespace voxel2tet
{

AbaqusExporter::AbaqusExporter(std :: vector< TriangleType * > *Triangles, std :: vector< VertexType * > *Vertices, std :: vector< EdgeType * > *Edges, std :: vector< TetType * > *Tets) :
    Exporter(Triangles, Vertices, Edges, Tets)
{

}

void AbaqusExporter::WriteVolumeData(std :: string Filename)
{
    STATUS( "Write .inp (Abaqus) file %s\n", Filename.c_str() );
    std :: ofstream AbaqusFile;
    AbaqusFile.open(Filename);

    // Prepare information

    // Used vertices
    this->UpdateUsedVertices();

    // Materials
    this->UpdateMaterialsMapping();

    // Find node sets
    this->UpdateMinMaxCoordinates();

    // Max and min nodes are arrays of lists of vertices where the index of the array tells in which direction the vertex is located
    this->UpdateMinMaxNodes();

    // Find element boundaries
    this->UpdateMinMaxElements();

    AbaqusFile << "**\n";
    AbaqusFile << "** File created by Voxel2Tet\n";
    AbaqusFile << "**\n";
    AbaqusFile << "*HEADING\n\n";

    AbaqusFile << "**\n**\n**\n";

    AbaqusFile << "*NODE\n";

    AbaqusFile.precision(10);

    // Write vertices/nodes
    for ( size_t i = 0; i < this->Vertices->size(); i++ ) {
        AbaqusFile << i + 1 << ",\t" << UsedVertices.at(i)->get_c(0) << ",\t" \
                  << UsedVertices.at(i)->get_c(1) << ",\t" << UsedVertices.at(i)->get_c(2) << "\n";
    }

    AbaqusFile << "**\n** SOLID ELEMENTS\n**\n";
    for (auto a: Self2OofemMaterials) {
        AbaqusFile << "*ELEMENT, TYPE=C3D4, ELSET=SOLID_" << a.first << "\n";
        for (TetType *t: *this->Tets) {
            if (t->MaterialID == a.first) {
                AbaqusFile << "  " << t->ID+1 << ",\t" << t->Vertices[0]->tag+1 << ",\t" <<
                              t->Vertices[1]->tag+1 << ",\t" << t->Vertices[2]->tag+1 <<
                              ",\t" << t->Vertices[3]->tag+1 << "\n";
            }
        }
        LOG("%u\n", a.first);
    }

    AbaqusFile << "**\n** SECTION DATA\n**\n";
    for (auto a: Self2OofemMaterials) {
        AbaqusFile << "*SOLID SECTION, ELSET=SOLID_" << a.first << ", MATERIAL=MATERIAL_" << a.first << "\n";
    }

    AbaqusFile << "**\n** MATERIALS\n**\n";
    int i=1;
    for (auto a: Self2OofemMaterials) {
        AbaqusFile << "*MATERIAL, NAME=MATERIAL_" << a.first << "\n";
        AbaqusFile << "*DENSITY\n\t1,\n";
        AbaqusFile << "ELASTIC, TYPE=ISOTROPIC\n\t " << (200e9+i*10e9) << ",\t0.3\n";
        i++;
    }

    // Boundary set
    AbaqusFile << "** Complete node boundary set\n";
    AbaqusFile << "*NSET, NSET=All surface nodes\n";
    for ( int i = 0; i < 3; i++ ) {
        // Write max nodes
        int k = 0;
        for ( size_t j = 0; j < MaxNodes [ i ].size(); j++ ) {
            AbaqusFile << "\t" << MaxNodes [ i ].at(j)->tag + 1 << ",";
            if (k==10) {
                AbaqusFile << "\n";
                k=0;
            } else {
                k++;
            }
        }

        // Write min nodes
        k = 0;
        for ( size_t j = 0; j < MinNodes [ i ].size(); j++ ) {
            AbaqusFile << "\t" << MinNodes [ i ].at(j)->tag + 1 << ",";
            if (k==10) {
                AbaqusFile << "\n";
                k=0;
            } else {
                k++;
            }
        }
    }
    AbaqusFile << "\n";

    // Individual max and min sets in all directions
    std :: array< std :: string, 3 >MaxNames = { { "MaxX", "MaxY", "MaxZ" } };
    std :: array< std :: string, 3 >MinNames = { { "MinX", "MinY", "MinZ" } };

    for ( int i = 0; i < 3; i++ ) {
        int k=0;

        // Write max nodes
        AbaqusFile << "*NSET, NSET=" << MaxNames [ i ] << "\n";
        for ( size_t j = 0; j < MaxNodes [ i ].size(); j++ ) {
            AbaqusFile << "\t" << MaxNodes [ i ].at(j)->tag + 1;
            if (j!=MaxNodes[i].size()-1) { AbaqusFile << ","; }
            if (k==10) {
                AbaqusFile << "\n";
                k=0;
            } else {
                k++;
            }
        }
        AbaqusFile << "\n";

        // Write min nodes
        k=0;
        AbaqusFile << "*NSET, NSET=" << MinNames [ i ] << "\n";
        for ( size_t j = 0; j < MinNodes [ i ].size(); j++ ) {
            AbaqusFile << "\t" << MinNodes [ i ].at(j)->tag + 1;
            if (j!=MinNodes[i].size()-1) { AbaqusFile << ","; }
            if (k==10) {
                AbaqusFile << "\n";
                k=0;
            } else {
                k++;
            }
        }
        AbaqusFile << "\n";
    }

    // Element boundaries

    // Complete boundary by element sides
    int ecount = 0;
    for ( std :: vector< int >s : MaxSide ) {
        ecount = ecount + s.size();
    }
    for ( std :: vector< int >s : MinSide ) {
        ecount = ecount + s.size();
    }

    AbaqusFile << "** Complete boundary by element sides\n";
    AbaqusFile << "*SURFACE, NAME=AllSurfaces, TYPE=ELEMENT\n";
    for ( int k = 0; k < 2; k++ ) {
        std :: array< std :: vector< TetType * >, 3 > *Elements = ( k == 0 ) ? & MaxElements : & MinElements;
        std :: array< std :: vector< int >, 3 > *Sides = ( k == 0 ) ? & MaxSide : & MinSide;
        for ( int i = 0; i < 3; i++ ) {
            for ( size_t j = 0; j < Elements->at(i).size(); j++ ) {
                AbaqusFile << " " << Elements->at(i).at(j)->ID + 1 << ", S" << Sides->at(i).at(j) << "\n";
            }
        }
    }
    AbaqusFile << "\n";

    // Boundaries in all directions
    for ( int k = 0; k < 2; k++ ) {
        std :: array< std :: string, 3 > *Comments = ( k == 0 ) ? & MaxNames : & MinNames;
        std :: array< std :: vector< TetType * >, 3 > *Elements = ( k == 0 ) ? & MaxElements : & MinElements;
        std :: array< std :: vector< int >, 3 > *Sides = ( k == 0 ) ? & MaxSide : & MinSide;

        for ( int i = 0; i < 3; i++ ) {
            AbaqusFile << "*SURFACE, NAME=Surface" << Comments->at(i) << "\n";
            for ( size_t j = 0; j < Elements->at(i).size(); j++ ) {
                AbaqusFile << "\t" << Elements->at(i).at(j)->ID + 1 << ", S" << Sides->at(i).at(j) << "\n";
            }
            AbaqusFile << "\n";
        }
    }
}

}
