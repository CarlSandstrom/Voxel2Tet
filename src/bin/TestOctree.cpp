#include <vector>
#include <algorithm>
#include <cmath>
#include <stdlib.h>
#include <ctime>


#include <VertexOctreeNode.h>
#include <MeshComponents.h>

int main()
{
    voxel2tet :: BoundingBoxType bb;
    for ( int i = 0; i < 3; i++ ) {
        bb.maxvalues [ i ] = 1.0;
        bb.minvalues [ i ] = 0.0;
    }

    std :: vector< voxel2tet :: VertexType * >VertexList;
    voxel2tet :: VertexOctreeNode RootNode(bb, & VertexList, 0);

    bool UseFixed = false;

    printf("Populate octree...\n");

    if ( UseFixed ) {
        RootNode.AddVertex(0.1, 0.1, 0.1);
        RootNode.AddVertex(0.5, 0.5, 0.5);
        RootNode.AddVertex(0.4, 0.4, 0.4);
        RootNode.AddVertex(0.11, 0.15, 0.1);
        RootNode.AddVertex(0.13, 0.15, 0.1);
        RootNode.AddVertex(0.19, 0.12, 0.1);
        RootNode.AddVertex(0.11, 0.15, 0.21);
        RootNode.AddVertex(0.1, 0.115, 0.1);
        RootNode.AddVertex(0.1, 0.1, 0.7);
        RootNode.AddVertex(0.11, 0.125, 0.1);
        RootNode.AddVertex(0.9, 0.15, 0.1);
        RootNode.printself();
    } else {
        srand( time(NULL) );

        int n = 10000000;
        for ( int i = 0; i < n; i++ ) {
            double x = double( rand() ) / double( RAND_MAX );
            double y = double( rand() ) / double( RAND_MAX );
            double z = double( rand() ) / double( RAND_MAX );
            RootNode.AddVertex(x, y, z);
        }
    }

    double x = 0.5, y = 0.5, z = 0.5, r = .1;

    std :: clock_t start = std :: clock();
    std :: vector< voxel2tet :: VertexType * >SphereVerticesOctree = RootNode.GiveVerticesWithinSphere(x, y, z, r);
    double OctreeDuration = ( std :: clock() - start )  / ( double ) CLOCKS_PER_SEC;


    std :: vector< voxel2tet :: VertexType * >SphereVerticesBrute;

    start = std :: clock();
    for ( voxel2tet :: VertexType *v : VertexList ) {
        double d = sqrt( ( v->get_c(0) - x ) * ( v->get_c(0) - x ) + ( v->get_c(1) - y ) * ( v->get_c(1) - y ) + ( v->get_c(2) - z ) * ( v->get_c(2) - z ) );
        if ( d <= r ) {
            SphereVerticesBrute.push_back(v);
        }
    }
    double LinearDuration = ( std :: clock() - start ) / ( double ) CLOCKS_PER_SEC;

    printf("Octree: %f Linear: %f\n", OctreeDuration, LinearDuration);

    std :: sort( SphereVerticesBrute.begin(), SphereVerticesBrute.end() );
    std :: sort( SphereVerticesOctree.begin(), SphereVerticesOctree.end() );

    if ( SphereVerticesBrute.size() == SphereVerticesOctree.size() ) {
        for ( unsigned int i = 0; i < SphereVerticesBrute.size(); i++ ) {
            if ( SphereVerticesBrute.at(i) != SphereVerticesOctree.at(i) ) {
                printf("Item at %u differs\n", i);
            }
        }
    } else {
        printf("Different size\n");
    }
}
