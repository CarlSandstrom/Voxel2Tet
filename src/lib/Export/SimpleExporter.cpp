#include "SimpleExporter.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>

namespace voxel2tet
{

SimpleExporter :: SimpleExporter(std :: vector< TriangleType * > *Triangles, std :: vector< VertexType * > *Vertices, std :: vector< EdgeType * > *Edges, std :: vector< TetType * > *Tets)
    : Exporter(Triangles, Vertices, Edges, Tets)
{

}

void  SimpleExporter::WriteSurfaceData(std :: string Filename)
{
    STATUS( "Write .simple file %s\n", Filename.c_str() );
    std :: ofstream SimpleFile;
    SimpleFile.open(Filename);

    // Prepare information
    std :: vector< VertexType * >UsedVertices;
    for ( TriangleType *t : *this->Triangles ) {
        for ( VertexType *v : t->Vertices ) {
            UsedVertices.push_back(v);
        }
    }

    std :: sort( UsedVertices.begin(), UsedVertices.end(), SortByID<VertexType *> );
    UsedVertices.erase( std :: unique( UsedVertices.begin(), UsedVertices.end() ), UsedVertices.end() );

    int i = 0;
    for ( VertexType *v : UsedVertices ) {
        v->tag = i;
        i++;
    }

    // Write header
    SimpleFile << "# Simple mesh file\n# Vertices\n" << UsedVertices.size() << "\n";

    // Write vertices
    for ( auto v : UsedVertices ) {
        SimpleFile << std :: setiosflags(std :: ios :: fixed) << std :: setprecision(15) << v->ID << " " << v->get_c(0) << " " << v->get_c(1) << " " << v->get_c(2) << "\n";
    }

    // Write edges
    SimpleFile << "# Edges\n" << this->Edges->size() << "\n";

    for (EdgeType *e: *this->Edges) {
        SimpleFile << e->ID << "\t" << e->Vertices[0]->ID << "\t" << e->Vertices[1]->ID << "\n";
    }

    // Write triangles
    SimpleFile << "# Triangles\n" << this->Triangles->size() << "\n";

    for (TriangleType *t: *this->Triangles) {
        SimpleFile << t->ID << "\t" << t->Vertices[0]->ID << "\t" << t->Vertices[1]->ID << "\t"<< t->Vertices[2]->ID << "\n";
    }

}

}
