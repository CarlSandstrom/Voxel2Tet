#include "OFFExporter.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>

namespace voxel2tet
{

OFFExporter :: OFFExporter(std :: vector< TriangleType * > *Triangles, std :: vector< VertexType * > *Vertices, std :: vector< EdgeType * > *Edges, std :: vector< TetType * > *Tets) :
    Exporter(Triangles, Vertices, Edges, Tets)
{
    LOG("Create OFFExporter object\n", 0);
}

void OFFExporter :: WriteSurfaceData(std :: string Filename)
{
    STATUS( "Write .off file %s\n", Filename.c_str() );
    std :: ofstream OFFFile;
    OFFFile.open(Filename);

    // Prepare information
    std :: vector< VertexType * >UsedVertices;
    for ( TriangleType *t : *this->Triangles ) {
        for ( VertexType *v : t->Vertices ) {
            UsedVertices.push_back(v);
        }
    }

    std :: sort( UsedVertices.begin(), UsedVertices.end() );
    UsedVertices.erase( std :: unique( UsedVertices.begin(), UsedVertices.end() ), UsedVertices.end() );

    int i = 0;
    for ( VertexType *v : UsedVertices ) {
        v->tag = i;
        i++;
    }

    // Write header
    OFFFile << "OFF\n" << UsedVertices.size() << "\t" << this->Triangles->size() << "\t0\n";

    // Write vertices
    for ( auto v : UsedVertices ) {
        OFFFile << std :: setiosflags(std :: ios :: fixed) << std :: setprecision(15) << v->get_c(0) << " " << v->get_c(1) << " " << v->get_c(2) << "\n";
    }

    // Write facets
    int j = 0;

    for ( auto t : *this->Triangles ) {
        OFFFile << 3;
        std :: array< int, 3 >VertexIDs;
        for ( int i = 0; i < 3; i++ ) {
            VertexIDs [ i ] = t->Vertices [ i ]->tag;
            OFFFile << " " << VertexIDs [ i ];
        }
        OFFFile << "\t#" << j << "\n";
        j++;
    }
}
}
