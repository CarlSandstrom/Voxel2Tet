#include <string>
#include <stdio.h>
#include <cmath>
#include <stdexcept>

#include "VertexOctreeNode.h"

namespace voxel2tet
{
VertexOctreeNode :: VertexOctreeNode(BoundingBoxType BoundingBox, std :: vector< VertexType * > *Vertices, int level)
{
    this->BoundingBox = BoundingBox;
    this->Vertices = Vertices;
    this->level = level;

    this->maxvertices = 20;
    this->eps = 1e-12;
}

VertexOctreeNode :: ~VertexOctreeNode()
{
    for ( auto c : this->children ) {
        delete c;
    }
}

VertexType *VertexOctreeNode :: FindVertexByCoords(double x, double y, double z)
{
    if ( this->children.size() == 0 ) {
        for ( auto VertexID : this->VertexIds ) {
            double d = std :: sqrt( std :: pow(this->Vertices->at(VertexID)->get_c(0) - x, 2) + std :: pow(this->Vertices->at(VertexID)->get_c(1) - y, 2) + std :: pow(this->Vertices->at(VertexID)->get_c(2) - z, 2) );
            if ( d < EPS ) {
                return this->Vertices->at(VertexID);
            }
        }
    } else {
        for ( auto c : this->children ) {
            if ( c->IsInBoundingBox(x, y, z) ) {
                return c->FindVertexByCoords(x, y, z);
            }
        }
    }
    return NULL;
}

int VertexOctreeNode :: AddVertex(double x, double y, double z)
{
    int newvertexid = -1;

    // If this is a leaf, locate the node and return the ID
    if ( this->children.size() == 0 ) {
        for ( auto VertexID : this->VertexIds ) {
            double d = std :: sqrt( std :: pow(this->Vertices->at(VertexID)->get_c(0) - x, 2) + std :: pow(this->Vertices->at(VertexID)->get_c(1) - y, 2) + std :: pow(this->Vertices->at(VertexID)->get_c(2) - z, 2) );
            if ( d < this->eps ) {
                return VertexID;
            }
        }
    }

    if ( ( this->children.size() == 0 ) & ( this->VertexIds.size() < std :: size_t(this->maxvertices) ) ) { // This is a leaf that is not full
        if ( this->IsInBoundingBox(x, y, z) == false ) {
            throw std :: out_of_range("Vertex is located outside the bounding box");
        }

        this->Vertices->push_back( new VertexType(x, y, z) );
        int VertexID = this->Vertices->size() - 1;
        this->Vertices->at(this->Vertices->size() - 1)->ID = VertexID;
        this->VertexIds.push_back(VertexID);

        return VertexID;
    } else if ( this->children.size() > 0 ) { // Has underlying nodes
        for ( VertexOctreeNode *child : this->children ) {
            if ( child->IsInBoundingBox(x, y, z) ) {
                newvertexid = child->AddVertex(x, y, z);
                return newvertexid;
            }
        }
    } else if ( this->VertexIds.size() == std :: size_t(this->maxvertices) ) { // This node is full. Split, then add to self
        this->split();
        newvertexid = this->AddVertex(x, y, z);
    }

    return newvertexid;
}

void VertexOctreeNode :: split()
{
    int newlevel = this->level + 1;
    double xmin = this->BoundingBox.minvalues [ 0 ];
    double xmax = this->BoundingBox.maxvalues [ 0 ];
    double ymin = this->BoundingBox.minvalues [ 1 ];
    double ymax = this->BoundingBox.maxvalues [ 1 ];
    double zmin = this->BoundingBox.minvalues [ 2 ];
    double zmax = this->BoundingBox.maxvalues [ 2 ];
    double xc = ( xmax + xmin ) / 2.0;
    double yc = ( ymax + ymin ) / 2.0;
    double zc = ( zmax + zmin ) / 2.0;

    BoundingBoxType b1;
    b1.maxvalues [ 0 ] = xc;
    b1.maxvalues [ 1 ] = yc;
    b1.maxvalues [ 2 ] = zc;
    b1.minvalues [ 0 ] = xmin;
    b1.minvalues [ 1 ] = ymin;
    b1.minvalues [ 2 ] = zmin;
    this->children.push_back( new VertexOctreeNode(b1, this->Vertices, newlevel) );

    b1.maxvalues [ 0 ] = xmax;
    b1.maxvalues [ 1 ] = yc;
    b1.maxvalues [ 2 ] = zc;
    b1.minvalues [ 0 ] = xc;
    b1.minvalues [ 1 ] = ymin;
    b1.minvalues [ 2 ] = zmin;
    this->children.push_back( new VertexOctreeNode(b1, this->Vertices, newlevel) );

    b1.maxvalues [ 0 ] = xmax;
    b1.maxvalues [ 1 ] = yc;
    b1.maxvalues [ 2 ] = zmax;
    b1.minvalues [ 0 ] = xc;
    b1.minvalues [ 1 ] = ymin;
    b1.minvalues [ 2 ] = zc;
    this->children.push_back( new VertexOctreeNode(b1, this->Vertices, newlevel) );

    b1.maxvalues [ 0 ] = xc;
    b1.maxvalues [ 1 ] = yc;
    b1.maxvalues [ 2 ] = zmax;
    b1.minvalues [ 0 ] = xmin;
    b1.minvalues [ 1 ] = ymin;
    b1.minvalues [ 2 ] = zc;
    this->children.push_back( new VertexOctreeNode(b1, this->Vertices, newlevel) );

    b1.maxvalues [ 0 ] = xc;
    b1.maxvalues [ 1 ] = ymax;
    b1.maxvalues [ 2 ] = zc;
    b1.minvalues [ 0 ] = xmin;
    b1.minvalues [ 1 ] = yc;
    b1.minvalues [ 2 ] = zmin;
    this->children.push_back( new VertexOctreeNode(b1, this->Vertices, newlevel) );

    b1.maxvalues [ 0 ] = xmax;
    b1.maxvalues [ 1 ] = ymax;
    b1.maxvalues [ 2 ] = zc;
    b1.minvalues [ 0 ] = xc;
    b1.minvalues [ 1 ] = yc;
    b1.minvalues [ 2 ] = zmin;
    this->children.push_back( new VertexOctreeNode(b1, this->Vertices, newlevel) );

    b1.maxvalues [ 0 ] = xmax;
    b1.maxvalues [ 1 ] = ymax;
    b1.maxvalues [ 2 ] = zmax;
    b1.minvalues [ 0 ] = xc;
    b1.minvalues [ 1 ] = yc;
    b1.minvalues [ 2 ] = zc;
    this->children.push_back( new VertexOctreeNode(b1, this->Vertices, newlevel) );

    b1.maxvalues [ 0 ] = xc;
    b1.maxvalues [ 1 ] = ymax;
    b1.maxvalues [ 2 ] = zmax;
    b1.minvalues [ 0 ] = xmin;
    b1.minvalues [ 1 ] = yc;
    b1.minvalues [ 2 ] = zc;
    this->children.push_back( new VertexOctreeNode(b1, this->Vertices, newlevel) );

    for ( auto VertexID : this->VertexIds ) {
        bool nodefound = false;
        for ( VertexOctreeNode *child : this->children ) {
            if ( child->IsInBoundingBox( this->Vertices->at(VertexID)->get_c(0), this->Vertices->at(VertexID)->get_c(1), this->Vertices->at(VertexID)->get_c(2) ) ) {
                child->VertexIds.push_back(VertexID);
                nodefound = true;
                break;
            }
        }
        if ( !nodefound ) {
            throw std :: out_of_range("Node not found");
        }
    }

    this->VertexIds.clear();
}

bool VertexOctreeNode :: IsInBoundingBox(double x, double y, double z)
{
    if ( ( x >= this->BoundingBox.minvalues [ 0 ] ) & ( y >= this->BoundingBox.minvalues [ 1 ] ) & ( z >= this->BoundingBox.minvalues [ 2 ] ) &
         ( x < this->BoundingBox.maxvalues [ 0 ] ) & ( y < this->BoundingBox.maxvalues [ 1 ] ) & ( z < this->BoundingBox.maxvalues [ 2 ] ) ) {
        return true;
    }
    return false;
}

std :: vector< VertexType * >VertexOctreeNode :: GiveVerticesWithinSphere(double x, double y, double z, double r)
{
    std :: vector< VertexType * >ResultList;

    if ( this->children.size() == 0 ) { // If this is a leaf, check all vertices
        for ( int VertexId : this->VertexIds ) {
            VertexType *v = this->Vertices->at(VertexId);
            double distance = std :: sqrt( ( v->get_c(0) - x ) * ( v->get_c(0) - x ) + ( v->get_c(1) - y ) * ( v->get_c(1) - y ) + ( v->get_c(2) - z ) * ( v->get_c(2) - z ) );
            if ( distance < r ) {
                ResultList.push_back(v);
            }
        }
    } else { // Check which children is, fully or partly, within the sphere
        double P [ 3 ] = {
            x, y, z
        };
        double d [ 3 ];

        for ( VertexOctreeNode *von : this->children ) {
            for ( int i : { 0, 1, 2 } ) {
                if ( P [ i ] > von->BoundingBox.maxvalues [ i ] ) {
                    d [ i ] = P [ i ] - von->BoundingBox.maxvalues [ i ];
                } else if ( P [ i ] < von->BoundingBox.minvalues [ i ] ) {
                    d [ i ] = von->BoundingBox.minvalues [ i ] - P [ i ];
                } else {
                    d [ i ] = 0.0;
                }
            }

            double d_scalar = std :: sqrt(d [ 0 ] * d [ 0 ] + d [ 1 ] * d [ 1 ] + d [ 2 ] * d [ 2 ]);

            if ( von->IsInBoundingBox(x, y, z) ) {
                d_scalar = 0.0;
            }

            if ( d_scalar <= r ) {
                std :: vector< VertexType * >ChildNodes = von->GiveVerticesWithinSphere(x, y, z, r);
                for ( VertexType *v : ChildNodes ) {
                    ResultList.push_back(v);
                }
            }
        }
    }
    return ResultList;
}

void VertexOctreeNode :: printself()
{
    std :: string tab;

    // Setup tabs for this node
    for ( int i = 0; i < this->level; i++ ) {
        tab = tab + "\t";
    }

    // Print bounding box
    printf("%s[[%f, %f, %f],[%f, %f, %f]]\n", tab.c_str(), this->BoundingBox.minvalues [ 0 ], this->BoundingBox.minvalues [ 1 ], this->BoundingBox.minvalues [ 2 ], this->BoundingBox.maxvalues [ 0 ], this->BoundingBox.maxvalues [ 1 ], this->BoundingBox.maxvalues [ 2 ]);

    // If node has children, print them, otherwise print vertices
    if ( this->children.size() > 0 ) {
        for ( VertexOctreeNode *child : this->children ) {
            child->printself();
        }
    } else {
        for ( auto VertexID : this->VertexIds ) {
            printf( "%s\t#%u: (%f, %f, %f)\n", tab.c_str(), VertexID, this->Vertices->at(VertexID)->get_c(0), this->Vertices->at(VertexID)->get_c(1), this->Vertices->at(VertexID)->get_c(2) );
        }
    }
}
}
