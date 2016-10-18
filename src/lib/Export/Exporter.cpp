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

    std :: sort( UsedVertices.begin(), UsedVertices.end(), SortByID<VertexType *> );
    UsedVertices.erase( std :: unique( UsedVertices.begin(), UsedVertices.end() ), UsedVertices.end() );

    int i = 0;
    for ( VertexType *v : UsedVertices ) {
        v->tag = i;
        i++;
    }
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
    // Find element boundaries

    for ( TetType *t : *this->Tets ) {
        for ( int k = 0; k < 2; k++ ) {   // Test max/min
            std :: array< std :: vector< TetType * >, 3 > *XElements = ( k == 0 ) ? & MaxElements : & MinElements;
            std :: array< std :: vector< int >, 3 > *XSide = ( k == 0 ) ? & MaxSide : & MinSide;
            for ( int i = 0; i < 3; i++ ) {   // Test direction
                std :: vector< int >TheNodes;

                for ( int j = 0; j < 4; j++ ) {   // Test node
                    VertexType *v = t->Vertices [ j ];
                    double cvalue = v->get_c(i);
                    if (k==0) {
                        if ( fabs(cvalue - MaxCoords [ i ]) < 1e-8 ) {
                            TheNodes.push_back(j + 1);                                 // +1 to match the numbering in the elemenent manual
                        }
                    } else {
                        if ( fabs(cvalue - MinCoords [ i ]) < 1e-8 ) {
                            TheNodes.push_back(j + 1);                                 // +1 to match the numbering in the elemenent manual
                        }
                    }
                }

                if ( TheNodes.size() == 3 ) {
                    int Side = -1;
                    if ( ( TheNodes [ 0 ] == 1 ) & ( TheNodes [ 1 ] == 2 ) & ( TheNodes [ 2 ] == 3 ) ) {
                        Side = 1;
                    }
                    if ( ( TheNodes [ 0 ] == 1 ) & ( TheNodes [ 1 ] == 2 ) & ( TheNodes [ 2 ] == 4 ) ) {
                        Side = 2;
                    }
                    if ( ( TheNodes [ 0 ] == 2 ) & ( TheNodes [ 1 ] == 3 ) & ( TheNodes [ 2 ] == 4 ) ) {
                        Side = 3;
                    }
                    if ( ( TheNodes [ 0 ] == 1 ) & ( TheNodes [ 1 ] == 3 ) & ( TheNodes [ 2 ] == 4 ) ) {
                        Side = 4;
                    }
                    XSide->at(i).push_back(Side);
                    XElements->at(i).push_back(t);
                }
            }
        }
    }
}

void Exporter :: UpdateGrainSets()
{
    for (TetType *t: *this->Tets) {
        int GrainID = t->MaterialID;
        std::pair< std::vector<TetType *>, int > *GrainSet = NULL;
        for (auto p: this->GrainSets) {
            if (p->second == GrainID) {
                GrainSet = p;
                break;
            }
        }
        if (GrainSet == NULL) {
            GrainSet = new std::pair< std::vector<TetType *>, int >;
            GrainSet->second = GrainID;
            this->GrainSets.push_back(GrainSet);
        }
        GrainSet->first.push_back(t);
    }
}

void Exporter :: UpdateTriangleSets()
{
    for (TriangleType *t: *this->Triangles) {
        int InterfaceID = t->InterfaceID;
        std::pair< std::vector<TriangleType *>, int > *TriangleSet = NULL;
        for (auto p: this->TriangleSets) {
            if (p->second == InterfaceID) {
                TriangleSet = p;
                break;
            }
        }
        if (TriangleSet == NULL) {
            TriangleSet = new std::pair< std::vector<TriangleType *>, int >;
            TriangleSet->second = InterfaceID;
            this->TriangleSets.push_back(TriangleSet);
        }
        TriangleSet->first.push_back(t);
    }
}

}
