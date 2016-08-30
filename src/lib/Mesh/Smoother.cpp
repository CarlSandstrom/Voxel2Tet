#include "Smoother.h"

namespace voxel2tet
{

std::vector<std::vector<VertexType *>> Smoother::GetConnectivityVector(std::vector<VertexType *> Vertices)
{
    std::vector<std::vector<VertexType *>> Connectivity;
    for (VertexType *v: Vertices) {
        std::vector<VertexType *> ConnectedVertices;

        // If v belongs to several phase edges, it is static and should not be effected by any other vertices
        if (!v->IsFixedVertex()) {
            // This will contain all connected vertices. Depending on type of v, some other vertices will be removed

            for (EdgeType *e: v->Edges) {
                if (!e->IsTransverse) {
                    for (VertexType *ve: e->Vertices) {
                        ConnectedVertices.push_back(ve);
                    }
                }
            }

            // Remove duplicate vertices
            std::sort(ConnectedVertices.begin(), ConnectedVertices.end());
            ConnectedVertices.erase( std :: unique( ConnectedVertices.begin(), ConnectedVertices.end() ), ConnectedVertices.end() );

            // Remove self
            ConnectedVertices.erase(std::remove(ConnectedVertices.begin(), ConnectedVertices.end(), v), ConnectedVertices.end());

            // If v is on a PhaseEdge, v should not be effected by vertices not belonging to the same PhaseEdge
            if (v->IsPhaseEdgeVertex()) {

                size_t i=0;
                while (i<ConnectedVertices.size()) {
                    VertexType *v2=ConnectedVertices[i];
                    if (!v2->IsPhaseEdgeVertex()) {
                        // v2 does not belong to any PhaseEdge
                        ConnectedVertices.erase(std::remove(ConnectedVertices.begin(), ConnectedVertices.end(), v2), ConnectedVertices.end());
                    } else {
                        // v2 belongs to a PhaseEdge. If it belongs to some other PhaseEdge, remove it.
                        bool BelongToSamePhaseEdge = false;
                        for (PhaseEdge *pe: v2->PhaseEdges) {
                            if (pe==v->PhaseEdges[0]) {
                                BelongToSamePhaseEdge = true;
                                break;
                            }
                        }
                        if (!BelongToSamePhaseEdge) {
                            ConnectedVertices.erase(std::remove(ConnectedVertices.begin(), ConnectedVertices.end(), v2), ConnectedVertices.end());
                        } else {
                            i++;
                        }
                    }
                }
            }

        }
        Connectivity.push_back(ConnectedVertices);
    }
    return Connectivity;
}

std :: vector< std :: pair< TriangleType *, TriangleType * > > Smoother::CheckPenetration(std :: vector< VertexType * > *Vertices, MeshManipulations *Mesh)
{
    std :: vector< std :: pair< TriangleType *, TriangleType * > >IntersectingTriangles;
    std :: vector< TriangleType * >Triangles;

    for ( VertexType *v : *Vertices ) {
        for ( TriangleType *t : v->Triangles ) {
            Triangles.push_back(t);
        }
    }
    std :: sort( Triangles.begin(), Triangles.end() );
    Triangles.erase( std :: unique( Triangles.begin(), Triangles.end() ), Triangles.end() );

    for ( TriangleType *t1 : Triangles ) {
        std :: array< double, 3 >c = t1->GiveCenterOfMass();
        double d = t1->GiveLongestEdgeLength();

        std :: vector< TriangleType * >NearTriangles = Mesh->GetTrianglesAround(c, Mesh->LongestEdgeLength);

        for ( TriangleType *t2 : NearTriangles ) {
            if ( t1 != t2 ) {
                if ( Mesh->CheckTrianglePenetration(t1, t2) ) {
                    //Mesh->CheckTrianglePenetration(t1, t2);
                    IntersectingTriangles.push_back({ t1, t2 });
                    LOG("Triangles %u and %u intersect!\n", t1->ID, t2->ID);
                    LOG( "t1(%u): (%f, %f, %f), (%f, %f, %f), (%f, %f, %f)\n", t1->ID, t1->Vertices.at(0)->get_c(0), t1->Vertices.at(0)->get_c(1), t1->Vertices.at(0)->get_c(2),
                         t1->Vertices.at(1)->get_c(0), t1->Vertices.at(1)->get_c(1), t1->Vertices.at(1)->get_c(2),
                         t1->Vertices.at(2)->get_c(0), t1->Vertices.at(2)->get_c(1), t1->Vertices.at(2)->get_c(2) );
                    LOG( "t2(%u): (%f, %f, %f), (%f, %f, %f), (%f, %f, %f)\n", t2->ID, t2->Vertices.at(0)->get_c(0), t2->Vertices.at(0)->get_c(1), t2->Vertices.at(0)->get_c(2),
                         t2->Vertices.at(1)->get_c(0), t2->Vertices.at(1)->get_c(1), t2->Vertices.at(1)->get_c(2),
                         t2->Vertices.at(2)->get_c(0), t2->Vertices.at(2)->get_c(1), t2->Vertices.at(2)->get_c(2) );
                }
            }
        }
    }
    return IntersectingTriangles;
}

void Smoother :: PullBackAtIntersections(std :: vector< VertexType * > Vertices, MeshManipulations *Mesh)
{

    // Check for intersecting triangles. If some triangles intersect, stiffen the structure in that area and re-smooth
    std :: vector< std :: pair< TriangleType *, TriangleType * > >IntersectingTriangles = CheckPenetration(&Vertices, Mesh);
    int intersecting_count = 0;
    // Pull back nodes until no intersecting triangles remain

    Mesh->ExportSurface(strfmt("/tmp/Intersection_step_%u.vtp", intersecting_count), FT_VTK);

    while ( IntersectingTriangles.size() > 0 ) {
        STATUS("Stiffen and re-smooth surface/edge due to %u intersecting triangles, iteration %u\n", IntersectingTriangles.size(), intersecting_count);

        // Stiffen spring
        std :: vector< VertexType * >TriangleVertices;
        for ( std :: pair< TriangleType *, TriangleType * >p : IntersectingTriangles ) {
            for ( VertexType *v : p.first->Vertices ) {
                TriangleVertices.push_back(v);
            }
            for ( VertexType *v : p.second->Vertices ) {
                TriangleVertices.push_back(v);
            }
        }

        std :: sort( TriangleVertices.begin(), TriangleVertices.end() );
        TriangleVertices.erase( std :: unique( TriangleVertices.begin(), TriangleVertices.end() ), TriangleVertices.end() );

        for ( VertexType *v : TriangleVertices ) {
            arma :: vec displacement = {
                v->get_c(0) - v->originalcoordinates [ 0 ], v->get_c(1) - v->originalcoordinates [ 1 ], v->get_c(2) - v->originalcoordinates [ 2 ]
            };
            arma :: vec newdisplacement = displacement * .9;

            for ( int i = 0; i < 3; i++ ) {
                v->set_c(v->originalcoordinates [ i ] + newdisplacement [ i ], i);
            }
        }

        Mesh->ExportSurface(strfmt("/tmp/Intersection_step_%u.vtp", intersecting_count+1), FT_VTK);

        intersecting_count++;
        IntersectingTriangles = CheckPenetration(&Vertices, Mesh);

    }

}

Smoother::Smoother()
{

}

}
