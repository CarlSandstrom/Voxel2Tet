#include "Exporter.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <map>
#include <vector>

namespace voxel2tet
{

Exporter :: Exporter(std :: vector< TriangleType * > *Triangles, std :: vector< VertexType * > *Vertices, std :: vector< EdgeType * > *Edges, std :: vector< TetType * > *Tets)
{
    LOG("Create exporter for MeshData@%p\n", Triangles);
    this->Vertices = Vertices;
    this->Edges = Edges;
    this->Triangles = Triangles;
    this->Tets = Tets;
}

void Exporter::UpdateUsedVertices()
{

    UsedVertices.clear();

    for ( TetType *t : *this->Tets ) {
        for ( VertexType *v : t->Vertices ) {
            UsedVertices.push_back(v);
        }
    }

    std :: sort( UsedVertices.begin(), UsedVertices.end() );
    UsedVertices.erase( std :: unique( UsedVertices.begin(), UsedVertices.end() ), UsedVertices.end() );


}

void Exporter::UpdateMaterialsMapping()
{
    for ( TetType *t : *this->Tets ) {
        if ( Self2OofemMaterials.find(t->MaterialID) == Self2OofemMaterials.end() ) {
            Self2OofemMaterials [ t->MaterialID ] = Self2OofemMaterials.size();
        }
    }
}

void Exporter::UpdateMinMaxCoordinates()
{
    MaxCoords = { { UsedVertices [ 0 ]->get_c(0), UsedVertices [ 0 ]->get_c(1), UsedVertices [ 0 ]->get_c(2) } };
    MinCoords = { { UsedVertices [ 0 ]->get_c(0), UsedVertices [ 0 ]->get_c(1), UsedVertices [ 0 ]->get_c(2) } };

    for ( VertexType *v : UsedVertices ) {
        for ( size_t i = 0; i < 3; i++ ) {
            if ( v->get_c(i) > MaxCoords [ i ] ) {
                MaxCoords [ i ] = v->get_c(i);
            }
            if ( v->get_c(i) < MinCoords [ i ] ) {
                MinCoords [ i ] = v->get_c(i);
            }
        }
    }

}

void Exporter::UpdateMinMaxNodes()
{
    double eps = 1e-8;

    for ( VertexType *v : UsedVertices ) {
        for ( int i = 0; i < 3; i++ ) {
            double cvalue = v->get_c(i);
            if ( fabs(cvalue - MaxCoords [ i ]) < eps ) {
                MaxNodes [ i ].push_back(v);
            }
            if ( fabs(cvalue - MinCoords [ i ]) < eps ) {
                MinNodes [ i ].push_back(v);
            }
        }
    }
}

void Exporter::UpdateMinMaxElements()
{

}

}
