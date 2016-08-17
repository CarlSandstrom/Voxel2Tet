#include "MiscFunctions.h"
#include "MeshManipulations.h"
#include "PhaseEdge.h"

namespace voxel2tet
{
MeshManipulations :: MeshManipulations(BoundingBoxType BoundingBox) : MeshData(BoundingBox)
{
    TOL_FLIP_MAXAREACHANGE = 1e-2;

    TOL_COL_SMALLESTAREA = 1e-8;
    TOL_COL_MAXNORMALCHANGE = 10 * 2 * 3.1415 / 360;
    TOL_COL_CHORD_MAXNORMALCHANGE = 10 * 2 * 3.141593 / 360;
    TOL_COL_MINANGLE = 20*2*3.1415/360;

    TOL_FLIP_SMALLESTAREA = 1e-8;
    TOL_FLIP_MAXNORMALCHANGE = 10 * 2 * 3.141593 / 360;
    TOL_FLIP_MAXNORMALDIFFERENCE = 10 * 2 * 3.1415 / 360;

    TOL_COL_MAXVOLUMECHANGE = .5 * .5 * .5 * 2;
    TOL_COL_MAXERROR = .5 * .5 * .5;
}

void MeshManipulations :: SortEdgesByLength()
{
    std :: vector< std :: pair< double, EdgeType * > >EdgeLength;

    for ( EdgeType *e : this->Edges ) {
        EdgeLength.push_back( std :: make_pair(e->GiveLength(), e) );
    }

    std :: sort( EdgeLength.begin(), EdgeLength.end() );

    this->Edges.clear();
    for ( std :: pair< double, EdgeType * >epair : EdgeLength ) {
        this->Edges.push_back(epair.second);
    }
}

void MeshManipulations :: SortEdgesByMinArea()
{
    std :: vector< std :: pair< double, EdgeType * > >EdgeArea;

    for ( EdgeType *e : this->Edges ) {
        std :: vector< TriangleType * >ts = e->GiveTriangles();

        unsigned int i;
        double smallestarea = 0.0;

        for ( i = 0; i < ts.size(); i++ ) {
            if ( ts.at(i)->GiveArea() < smallestarea ) {
                smallestarea = ts.at(i)->GiveArea();
            }
        }

        EdgeArea.push_back( std :: make_pair(smallestarea, e) );
    }

    std :: sort( EdgeArea.begin(), EdgeArea.end() );

    this->Edges.clear();
    for ( std :: pair< double, EdgeType * >epair : EdgeArea ) {
        this->Edges.push_back(epair.second);
    }
}

FC_MESH MeshManipulations :: GetFlippedEdgeData(EdgeType *EdgeToFlip, EdgeType *NewEdge, std :: array< TriangleType *, 2 > *NewTriangles)
{
    LOG("Get flipped edge data for edge %u@%p\n", EdgeToFlip->ID, EdgeToFlip);

    std :: vector< TriangleType * >EdgeTriangles = EdgeToFlip->GiveTriangles();

    if ( EdgeTriangles.size() != 2 ) {
        LOG("Unable to flip edge. To many or only one triangle connected\n", 0);
        return FC_TOOMANYTRIANGLES;
    }

    // Ensure that both triangles are oriented in the same way
    if ( !this->CheckSameOrientation(EdgeTriangles [ 0 ], EdgeTriangles [ 1 ]) ) {
        EdgeTriangles [ 0 ]->FlipNormal();
    }

    // Find shared edge and not shared verticec
    std :: array< VertexType *, 2 >NewEdgeVertices;

    bool finished = false;

    for ( int i = 0; i < 3; i++ ) {
        for ( int j = 0; j < 3; j++ ) {
            if ( EdgeTriangles [ 0 ]->Vertices [ i ] == EdgeTriangles [ 1 ]->Vertices [ j ] ) {
                bool t0forward = EdgeTriangles [ 0 ]->Vertices [ ( i + 1 ) % 3 ] == EdgeTriangles [ 1 ]->Vertices [ ( j + 2 ) % 3 ];

                LOG("t0forward: %d\n", t0forward);

                int increment = ( t0forward ) ? 1 : 2;
                std :: array< VertexType *, 2 >t0edge = { { EdgeTriangles [ 0 ]->Vertices [ i ], EdgeTriangles [ 0 ]->Vertices [ ( i + increment ) % 3 ] } };

                // Produce new edge
                if ( t0forward ) {
                    NewEdgeVertices [ 0 ] = EdgeTriangles [ 0 ]->Vertices [ ( i + 2 ) % 3 ];
                    NewEdgeVertices [ 1 ] = EdgeTriangles [ 1 ]->Vertices [ ( j + 1 ) % 3 ];
                } else {
                    NewEdgeVertices [ 1 ] = EdgeTriangles [ 0 ]->Vertices [ ( i + 1 ) % 3 ];
                    NewEdgeVertices [ 0 ] = EdgeTriangles [ 1 ]->Vertices [ ( j + 2 ) % 3 ];
                }

                //                finished = true;
                //                break;

                // Produce new triangles
                TriangleType *new_t0, *new_t1;
                new_t0 = new TriangleType({ NewEdgeVertices [ 1 ], NewEdgeVertices [ 0 ], t0edge [ 0 ] });
                new_t1 = new TriangleType({ NewEdgeVertices [ 0 ], NewEdgeVertices [ 1 ], t0edge [ 1 ] });

                new_t0->PosNormalMatID = new_t1->PosNormalMatID = EdgeTriangles [ 0 ]->PosNormalMatID;
                new_t0->NegNormalMatID = new_t1->NegNormalMatID = EdgeTriangles [ 0 ]->NegNormalMatID;
                new_t0->InterfaceID = new_t1->InterfaceID = EdgeTriangles [ 0 ]->InterfaceID;

                NewTriangles->at(0) = new_t0;
                NewTriangles->at(1) = new_t1;

                NewEdge->Vertices [ 0 ] = NewEdgeVertices [ 0 ];
                NewEdge->Vertices [ 1 ] = NewEdgeVertices [ 1 ];

                /*                std::array<double, 3> OldNormal0=EdgeTriangles[0]->GiveUnitNormal();
                 *              std::array<double, 3> OldNormal1=EdgeTriangles[1]->GiveUnitNormal();
                 *              LOG("Old normal 0: [%f, %f, %f]\n", OldNormal0[0], OldNormal0[1], OldNormal0[2]);
                 *              LOG("Old normal 1: [%f, %f, %f]\n", OldNormal1[0], OldNormal1[1], OldNormal1[2]);
                 *
                 *              std::array<double, 3> NewNormal0=new_t0->GiveUnitNormal();
                 *              std::array<double, 3> NewNormal1=new_t1->GiveUnitNormal();
                 *              LOG("New normal 0: [%f, %f, %f]\n", NewNormal0[0], NewNormal0[1], NewNormal0[2]);
                 *              LOG("New normal 1: [%f, %f, %f]\n", NewNormal1[0], NewNormal1[1], NewNormal1[2]);*/

                return FC_OK;
            }
        }
        if ( finished ) {
            break;
        }
    }

    return FC_INVALIDEDGE;
}

FC_MESH MeshManipulations :: FlipEdge(EdgeType *Edge)
{

    LOG("Flip edge %u@%p (%u, %u)\n", Edge->ID, Edge, Edge->Vertices [ 0 ]->ID, Edge->Vertices [ 1 ]->ID);

    // this->DoSanityCheck();

    std :: vector< TriangleType * >EdgeTriangles = Edge->GiveTriangles();

    // Check if edge belongs to too many triangles
    if ( EdgeTriangles.size() != 2 ) {
        LOG("\tUnable to flip edge. To many or only one triangle connected\n", 0);
        for ( TriangleType *t : EdgeTriangles ) {
            LOG("\t\t%u\n", t->ID);
        }
        return FC_TOOMANYTRIANGLES;
    }

    // Check if triangels belongs to same surface
    if ( EdgeTriangles [ 0 ]->InterfaceID != EdgeTriangles [ 1 ]->InterfaceID ) {
        LOG("\tUnable to flip edge. Triangels belongs to different surfaces\n", 0);
        return FC_DIFFERENTSURFACES;
    }

    // Check angle between elements
    double Angle = ComputeAngleBetweenVectors( EdgeTriangles [ 0 ]->GiveUnitNormal(), EdgeTriangles [ 1 ]->GiveUnitNormal() );
    if ( std :: min( fabs(Angle), fabs(Angle - 3.141593) )  > TOL_FLIP_MAXNORMALDIFFERENCE ) {
        return FC_NORMAL;
    }

    std :: array< TriangleType *, 2 >NewTriangles;
    EdgeType NewEdge;

    FC_MESH FC;

    // Get new data for flip (i.e new triangles and new edge)
    FC = this->GetFlippedEdgeData(Edge, & NewEdge, & NewTriangles);
    if ( FC != FC_OK ) {
        return FC;
    }

    // Check if normal changes too much
    FC = this->CheckFlipNormal(& EdgeTriangles, NewTriangles);
    if ( FC != FC_OK ) {
        LOG("\tUnable to flip edge. Changes in normal direction prevents flipping\n", 0);
        for ( TriangleType *t : NewTriangles ) {
            delete t;
            t = NULL;
        }
        return FC;
    }

    // Ensure that minimal angle is improved if we continue
    double minAngleCurrent = std :: min( EdgeTriangles [ 0 ]->GiveSmallestAngle(), EdgeTriangles [ 1 ]->GiveSmallestAngle() );
    double minAngleNew = std :: min( NewTriangles [ 0 ]->GiveSmallestAngle(), NewTriangles [ 1 ]->GiveSmallestAngle() );

    if ( minAngleNew < minAngleCurrent ) {
        LOG("\tUnable to flip edge. New minimal angles worse than current (New: %f, Current: %f).\n", minAngleNew, minAngleCurrent);
        for ( TriangleType *t : NewTriangles ) {
            delete t;
            t = NULL;
        }
        return FC_WORSEMINANGLE;
    }

    // If the minimal angle does not change, prevent flipping as this has no value
    if ( fabs(minAngleNew - minAngleCurrent) < 1e-8 ) {
        LOG("\tUnable to flip edge. Flipping does not improve quality\n", minAngleNew, minAngleCurrent);
        for ( TriangleType *t : NewTriangles ) {
            delete t;
            t = NULL;
        }
        return FC_ANGLESNOTIMPROVED;
    }

    // Check change in area of region. Should be small
    double CurrentArea = EdgeTriangles.at(0)->GiveArea() + EdgeTriangles.at(1)->GiveArea();
    double NewArea = NewTriangles.at(0)->GiveArea() + NewTriangles.at(1)->GiveArea();

    if ( std :: fabs(CurrentArea - NewArea) > TOL_FLIP_MAXAREACHANGE ) {
        LOG("The combined area of the triangles changes too much. Prevent flipping.\n", 0);
        for ( TriangleType *t : NewTriangles ) {
            delete t;
            t = NULL;
        }
        return FC_AREACHANGETOOLARGE;
    }

    // Check if any of the new triangles already exist
    std :: vector< TriangleType * >ConnectedTriangles;
    for ( VertexType *v : NewEdge.Vertices ) {
        for ( TriangleType *t : v->Triangles ) {
            ConnectedTriangles.push_back(t);
        }
    }

    std :: sort( ConnectedTriangles.begin(), ConnectedTriangles.end() );

    for ( unsigned int i = 0; i < ConnectedTriangles.size() - 1; i++ ) {
        if ( ConnectedTriangles.at(i) == ConnectedTriangles.at(i + 1) ) {
            LOG("Edge already exists!\n", 0);
            return FC_INVALIDEDGE;
        }
    }

    // Check if any of the new triangles penetrates any of the old triangles nearby

    // Add all triangles nearby to NearTriangles
    std :: array< double, 3 >cm = NewEdge.GiveCenterPoint();
    std :: vector< TriangleType * >NearTriangles = this->GetTrianglesAround(cm, this->LongestEdgeLength);

    // Remove EdgeTriangles from NearTriangles
    for ( int i = 0; i < 2; i++ ) {
        NearTriangles.erase( std :: remove(NearTriangles.begin(), NearTriangles.end(), EdgeTriangles [ i ]), NearTriangles.end() );
    }

    // Check if new triangles penetrates existing triangles (except those that will be deleted of course)
    for ( TriangleType *t1 : NewTriangles ) {
        for ( TriangleType *t2 : NearTriangles ) {
            FC_MESH R = this->CheckTrianglePenetration(t1, t2);
            if ( R != FC_OK ) {
                LOG("Unable to flip edge. Will result in penetration\n", 0);
                return R;
            }
        }
    }


    LOG("Ok to flip edge\n", 0);

    // Update mesh data

    // Remove triangles from verices
    for ( int i : { 0, 1 } ) {
        for ( TriangleType *t : EdgeTriangles ) {
            Edge->Vertices [ i ]->RemoveTriangle(t);
        }
    }

    // Remove edge from vertices
    for ( int i : { 0, 1 } ) {
        Edge->Vertices [ i ]->RemoveEdge(Edge);
    }

    // Add new edge to vertices
    for ( int i : { 0, 1 } ) {
        Edge->Vertices [ i ] = NewEdge.Vertices [ i ];
        Edge->Vertices [ i ]->AddEdge(Edge);
    }

    // Add new triangles list (and thus also to vertices)
    for ( TriangleType *t : NewTriangles ) {
        this->AddTriangle(t);
    }

    // Remove edge from list
    //this->RemoveEdge(Edge);
    //this->Edges.erase(std::remove(this->Edges.begin(), this->Edges.end(), Edge), this->Edges.end());

    // 5. Free old triangles
    for ( TriangleType *t : EdgeTriangles ) {
        this->RemoveTriangle(t);
    }

    LOG("Edge %u@%p (%u, %u) is now flipped\n", Edge->ID, Edge, Edge->Vertices[0]->ID, Edge->Vertices[1]->ID);

    return FC_OK;
}

FC_MESH MeshManipulations :: CheckFlipNormal(std :: vector< TriangleType * > *OldTriangles, std :: array< TriangleType *, 2 >NewTriangles)
{
    double MaxAngle = 0.0;

    for ( TriangleType *OldTriangle : *OldTriangles ) {
        std :: array< double, 3 >OldNormal = OldTriangle->GiveUnitNormal();

        for ( unsigned int j = 0; j < NewTriangles.size(); j++ ) {
            TriangleType *NewTriangle = NewTriangles.at(j);
            double NewArea = NewTriangle->GiveArea();

            if ( NewArea < TOL_FLIP_SMALLESTAREA ) {
                return FC_SMALLAREA;
            }

            std :: array< double, 3 >NewNormal = NewTriangle->GiveUnitNormal();

            // Compute angle between new and old normal
            double angle1 = std :: acos(OldNormal [ 0 ] * NewNormal [ 0 ] + OldNormal [ 1 ] * NewNormal [ 1 ] + OldNormal [ 2 ] * NewNormal [ 2 ]);
            double angle2 = std :: acos( -( OldNormal [ 0 ] * NewNormal [ 0 ] + OldNormal [ 1 ] * NewNormal [ 1 ] + OldNormal [ 2 ] * NewNormal [ 2 ] ) );
            double MinAngle = std :: min(angle1, angle2);
            MaxAngle = std :: max(MaxAngle, MinAngle);
        }
    }

    if ( MaxAngle > TOL_FLIP_MAXNORMALCHANGE ) {
        return FC_NORMAL;
    } else {
        return FC_OK;
    }
}

FC_MESH MeshManipulations :: CollapseEdgeTest(std :: vector< TriangleType * > *TrianglesToSave, std :: vector< TriangleType * > *TrianglesToRemove, std :: vector< TriangleType * > *NewTriangles, EdgeType *EdgeToCollapse, int RemoveVertexIndex)
{
    int SaveVertexIndex = ( RemoveVertexIndex == 0 ) ? 1 : 0;

    // Check if NewTriangles are good replacements
    LOG("Check validity of new normals...\n", 0);
    FC_MESH FC;

    double error;
    FC = this->CheckCoarsenNormalImproved(TrianglesToSave, TrianglesToRemove, NewTriangles, error);
    if ( FC != FC_OK ) {
        return FC;
    }

    // Check if accumulated error will exceed limit
    std :: vector< VertexType * >eVertices;
    for ( TriangleType *t : *TrianglesToSave ) {
        for ( int i = 0; i < 3; i++ ) {
            eVertices.push_back(t->Vertices [ i ]);
        }
    }
    std :: sort( eVertices.begin(), eVertices.end() );
    eVertices.erase( std :: unique( eVertices.begin(), eVertices.end() ), eVertices.end() );

    double TotalError = 0;
    for ( VertexType *v : eVertices ) {
        TotalError = TotalError + v->error;
    }

    if ( TotalError > TOL_COL_MAXERROR ) {
        return FC_TOOLARGEERROR;
    }

    LOG("Check validity of new chord...\n", 0);
    FC = this->CheckCoarsenChord(EdgeToCollapse, EdgeToCollapse->Vertices [ RemoveVertexIndex ], EdgeToCollapse->Vertices [ SaveVertexIndex ]);
    if ( FC != FC_OK ) {
        return FC;
    }

    LOG("Check area of new triangles...\n", 0);
    // A too small area should be ok as long as it is smaller than the previous one
    for ( TriangleType *t : *NewTriangles ) {
        if ( t->GiveArea() < 1e-7 ) {
            LOG(" - Check failed\n", 0);
            return FC_SMALLAREA;
        }
        if (t->GiveSmallestAngle() < TOL_COL_MINANGLE ) {
            LOG(" - Check failed. To small angle on new triangle.", 0);
            return FC_SMALLANGLE;
        }
    }

    LOG("Check if new triangles intersect...\n", 0);

    // Collect all point close to the centerpoint of the edge to collapse. The, form a list of all triangles connected to those points and perform check on all triangles in that list (except with triangles to remove).
    std :: array< double, 3 >cp = EdgeToCollapse->GiveCenterPoint();
    std :: vector< TriangleType * >TrianglesNear;
    TrianglesNear = this->GetTrianglesAround(cp, this->LongestEdgeLength);

    for ( TriangleType *t : TrianglesNear ) {
        // Skip triangles that will be removed
        bool DoCheck = true;
        for ( TriangleType *remt : *TrianglesToRemove ) {
            if ( remt == t ) {
                DoCheck = false;
            }
        }

        if ( DoCheck ) {
            // For this test we need to change RemoveVertex --> SaveVertex
            std :: array< VertexType *, 3 >NearTriVertices = t->Vertices;
            for ( int i = 0; i < 3; i++ ) {
                if ( NearTriVertices [ i ] == EdgeToCollapse->Vertices [ RemoveVertexIndex ] ) {
                    NearTriVertices [ i ] = EdgeToCollapse->Vertices [ SaveVertexIndex ];
                    break;
                }
            }

            for ( TriangleType *newt : *NewTriangles ) {
                std :: array< VertexType *, 3 >t1 = {
                    NearTriVertices [ 0 ], NearTriVertices [ 1 ], NearTriVertices [ 2 ]
                };
                std :: array< VertexType *, 3 >t2 = {
                    newt->Vertices [ 0 ], newt->Vertices [ 1 ], newt->Vertices [ 2 ]
                };

                int sv;
                FC_MESH nf = CheckTrianglePenetration(t1, t2, sv);

                if ( nf == FC_DUPLICATETRIANGLE ) {
                    if ( newt->ID == -t->ID ) {
                        nf = FC_OK;
                    }
                }

                if ( nf != FC_OK ) {
                    return nf;
                }
            }
        }
    }

    // Add new error
    for ( VertexType *v : eVertices ) {
        v->error = v->error + error / ( ( double ) eVertices.size() );
    }

    return FC_OK;
}

FC_MESH MeshManipulations :: CollapseEdge(EdgeType *EdgeToCollapse, int RemoveVertexIndex, bool PerformTesting)
{

    LOG("Collapse edge %u@%p (%u, %u) by removing vertex %u\n", EdgeToCollapse->ID, EdgeToCollapse, EdgeToCollapse->Vertices [ 0 ]->ID,
            EdgeToCollapse->Vertices [ 1 ]->ID, EdgeToCollapse->Vertices [ RemoveVertexIndex ]->ID);

    // Cannot remove a fixed vertex
    if ( EdgeToCollapse->Vertices [ RemoveVertexIndex ]->IsFixedVertex() ) {
        return FC_FIXEDVERTEX;
    }

    int SaveVertexIndex = ( RemoveVertexIndex == 0 ) ? 1 : 0;

    // If the vertex to remove belongs to a PhaseEdge, ensure that the other vertex belongs to the same PhaseEdge
    if ( EdgeToCollapse->Vertices [ RemoveVertexIndex ]->PhaseEdges.size() > 0 ) {
        if ( EdgeToCollapse->Vertices [ SaveVertexIndex ]->PhaseEdges.size() == 0 ) {
            return FC_CHORD;
        }
    }

    VertexType *SaveVertex = EdgeToCollapse->Vertices [ SaveVertexIndex ];

    // Create new triangles. These are create by moving RemoveVertex to the other end of the edge and remove the 0-area triangles
    std :: vector< TriangleType * >TrianglesToRemove = EdgeToCollapse->GiveTriangles();

    LOG("Connected triangle IDs: %u, %u\n", TrianglesToRemove.at(0)->ID, TrianglesToRemove.at(1)->ID);
    std :: vector< TriangleType * >ConnectedTriangles = EdgeToCollapse->Vertices [ RemoveVertexIndex ]->Triangles;

    std :: sort( TrianglesToRemove.begin(), TrianglesToRemove.end() );
    std :: sort( ConnectedTriangles.begin(), ConnectedTriangles.end() );

    std :: vector< TriangleType * >TrianglesToSave;
    std :: set_difference( ConnectedTriangles.begin(), ConnectedTriangles.end(),
                           TrianglesToRemove.begin(), TrianglesToRemove.end(), std :: inserter( TrianglesToSave, TrianglesToSave.begin() ) );

    // NewTriangles is the updated subset of of TrianglesToSave with the removed vertex changed to the saved vertex
    std :: vector< TriangleType * >NewTriangles;
    for ( TriangleType *t : TrianglesToSave ) {
        t->UpdateNormal();
        TriangleType *NewTriangle = new TriangleType;

        // Copy data
        NewTriangle->InterfaceID = t->InterfaceID;
        NewTriangle->ID = -t->ID; // Minus to distinguish from existing triangles
        NewTriangle->Vertices = t->Vertices;
        NewTriangle->PosNormalMatID = t->PosNormalMatID;
        NewTriangle->NegNormalMatID = t->NegNormalMatID;

        NewTriangle->UpdateNormal();

        for ( int i = 0; i < 3; i++ ) {
            if ( t->Vertices [ i ] == EdgeToCollapse->Vertices [ RemoveVertexIndex ] ) {
                NewTriangle->Vertices [ i ] = EdgeToCollapse->Vertices [ SaveVertexIndex ];
            } else {
                NewTriangle->Vertices [ i ] = t->Vertices [ i ];
            }
        }

        NewTriangle->UpdateNormal();
        NewTriangles.push_back(NewTriangle);
    }


    FC_MESH FC;
    if ( PerformTesting ) {
        FC = this->CollapseEdgeTest(& TrianglesToSave, & TrianglesToRemove, & NewTriangles, EdgeToCollapse, RemoveVertexIndex);
    } else {
        FC = FC_OK;
    }

    if ( FC != FC_OK ) {
        for ( TriangleType *t : NewTriangles ) {
            delete t;
            t = NULL;
        }
        return FC;
    }

    LOG("Checks ok. Proceed\n", 0);

    // Update triangulation

    // Find edges to remove (all edges connected to RemoveVertex and in any triangle in TrianglesToRemove)
    std :: vector< EdgeType * >EdgesToRemove;
    std :: vector< EdgeType * >RemoveVertexEdges = EdgeToCollapse->Vertices.at(RemoveVertexIndex)->Edges;

    std :: vector< EdgeType * >TriangleToRemoveEdges;
    for ( TriangleType *t : TrianglesToRemove ) {
        for ( EdgeType *e : t->GiveEdges() ) {
            TriangleToRemoveEdges.push_back(e);
        }
    }

    std :: sort( RemoveVertexEdges.begin(), RemoveVertexEdges.end() );
    std :: sort( TriangleToRemoveEdges.begin(), TriangleToRemoveEdges.end() );

    std :: set_intersection( RemoveVertexEdges.begin(), RemoveVertexEdges.end(),
                             TriangleToRemoveEdges.begin(), TriangleToRemoveEdges.end(), std :: back_inserter(EdgesToRemove) );

    // Update edges
    std :: vector< EdgeType * > ConnectedEdges;
    std::set_difference(RemoveVertexEdges.begin(), RemoveVertexEdges.end(), EdgesToRemove.begin(), EdgesToRemove.end(), std::back_inserter(ConnectedEdges));

    // Ensure that we don't end up with copies edges, i.e. moves one edge onto another. This means that we "snap of" a volume
    std::vector<EdgeType *> SaveVertexEdges = EdgeToCollapse->Vertices[SaveVertexIndex]->Edges;

    for (EdgeType *se: SaveVertexEdges) {
        if (se!=EdgeToCollapse) {
            for (EdgeType *ce: ConnectedEdges) {

                // Construct updated edge
                EdgeType UpdatedConnectedEdge = *ce;
                for (size_t i=0; i<2; i++) {
                    if (UpdatedConnectedEdge.Vertices[i] == EdgeToCollapse->Vertices[RemoveVertexIndex]) {
                        UpdatedConnectedEdge.Vertices[i] = EdgeToCollapse->Vertices[SaveVertexIndex];
                    }
                }

                //Compare
                if ( ( (se->Vertices[0] == UpdatedConnectedEdge.Vertices[0]) & (se->Vertices[1] == UpdatedConnectedEdge.Vertices[1]) ) |
                     ( (se->Vertices[0] == UpdatedConnectedEdge.Vertices[1]) & (se->Vertices[1] == UpdatedConnectedEdge.Vertices[0]) ) ) {
                    return FC_INVALIDEDGE;
                }

            }
        }
    }

    for ( EdgeType *e : ConnectedEdges ) {
        if ( e != EdgeToCollapse ) {
            for ( int i : { 0, 1 } ) {
                // "Move" vertex
                if ( e->Vertices [ i ] == EdgeToCollapse->Vertices.at(RemoveVertexIndex) ) {
                    e->Vertices [ i ]->RemoveEdge(e);
                    e->Vertices [ i ] = EdgeToCollapse->Vertices.at(SaveVertexIndex);
                    e->Vertices [ i ]->AddEdge(e);
                }
            }
        }
    }

    // Remove triangles
    for ( TriangleType *t : ConnectedTriangles ) {
        this->RemoveTriangle(t);
    }

    // Remove edges
    for ( EdgeType *e : EdgesToRemove ) {
        this->RemoveEdge(e);
    }

    // Add new triangles
    for ( TriangleType *t : NewTriangles ) {
        this->AddTriangle(t);
    }

    // DoSanityCheck();

    ConnectedEdges = SaveVertex->Edges;
    bool edgeflipped = true;
    while ( edgeflipped ) {
        edgeflipped = false;
        size_t i = 0;
        while ( i < SaveVertex->Edges.size() ) {
            EdgeType *e = SaveVertex->Edges [ i ];
            if ( FlipEdge(e) != FC_OK ) {
                //FlipEdge(e);
                LOG("Failed to flip edge\n", 0);
                i++;
            } else {
                LOG("Edge flipped!\n", 0);
                edgeflipped = true;
            }
        }
    }

    return FC_OK;
}

FC_MESH MeshManipulations :: CheckCoarsenNormal(std :: vector< TriangleType * > *OldTriangles, std :: vector< TriangleType * > *NewTriangles)
{
    // Here, we use that the order of the old and new triangles are the same in both lists. I.e. NewTriangle[i] is the same triangle as OldTriangles[i] except that vD is replaced by vR.

    double MaxAngle = 0;

    for ( unsigned int i = 0; i < OldTriangles->size(); i++ ) {
        std :: array< double, 3 >OldNormal = OldTriangles->at(i)->GiveUnitNormal();
        std :: array< double, 3 >NewNormal = NewTriangles->at(i)->GiveUnitNormal();

        double angle1 = std :: acos(OldNormal [ 0 ] * NewNormal [ 0 ] + OldNormal [ 1 ] * NewNormal [ 1 ] + OldNormal [ 2 ] * NewNormal [ 2 ]);
        double angle2 = std :: acos( -( OldNormal [ 0 ] * NewNormal [ 0 ] + OldNormal [ 1 ] * NewNormal [ 1 ] + OldNormal [ 2 ] * NewNormal [ 2 ] ) );

        double oldarea  = OldTriangles->at(i)->GiveArea();
        double newarea  = NewTriangles->at(i)->GiveArea();

        if ( newarea < TOL_COL_SMALLESTAREA ) {
            return FC_SMALLAREA;
        }

        LOG("angle1=%f,\tangle2=%f\tOldArea=%f\tNewArea=%f\n", angle1, angle2, oldarea, newarea);

        if ( std :: min(angle1, angle2) > MaxAngle ) {
            MaxAngle = std :: min(angle1, angle2);
        }

        if ( angle1 > TOL_COL_MAXNORMALCHANGE ) {
            LOG("Should not be collapsed... (angle1=%f) \n", angle1);
            return FC_NORMAL;
        }
    }

    if ( MaxAngle > TOL_COL_MAXNORMALCHANGE ) {
        return FC_NORMAL;
    } else {
        return FC_OK;
    }
}

FC_MESH MeshManipulations :: CheckCoarsenNormalImproved(std :: vector< TriangleType * > *OldTriangles, std :: vector< TriangleType * > *TrianglesToRemove, std :: vector< TriangleType * > *NewTriangles, double &error)
{
    // This is very unfinished....

    double MaxAngle = 0;

    for ( unsigned int i = 0; i < OldTriangles->size(); i++ ) {
        std :: array< double, 3 >OldNormal = OldTriangles->at(i)->GiveUnitNormal();
        std :: array< double, 3 >NewNormal = NewTriangles->at(i)->GiveUnitNormal();

        double angle1 = std :: acos(OldNormal [ 0 ] * NewNormal [ 0 ] + OldNormal [ 1 ] * NewNormal [ 1 ] + OldNormal [ 2 ] * NewNormal [ 2 ]);
        double angle2 = std :: acos( -( OldNormal [ 0 ] * NewNormal [ 0 ] + OldNormal [ 1 ] * NewNormal [ 1 ] + OldNormal [ 2 ] * NewNormal [ 2 ] ) );

        double oldarea  = OldTriangles->at(i)->GiveArea();
        double newarea  = NewTriangles->at(i)->GiveArea();

        LOG("angle1=%f,\tangle2=%f\tOldArea=%f\tNewArea=%f\n", angle1, angle2, oldarea, newarea);

        if ( newarea < TOL_COL_SMALLESTAREA ) {
            return FC_SMALLAREA;
        }

        MaxAngle = std :: max( MaxAngle, std :: min( fabs(angle1), fabs(angle2) ) );

        if ( angle1 > 3.141592 / 2 ) {
            LOG("Should not be collapsed...\n", 0);
            return FC_NORMAL;
        }
    }

    if ( MaxAngle > 90 * 2 * 3.141593 / 360 ) {
        return FC_NORMAL;
    }

    std :: vector< TriangleType * >RemoveVolume = * OldTriangles;
    for ( TriangleType *t : *TrianglesToRemove ) {
        RemoveVolume.push_back(t);
    }

    double V = 0.0;
    double ContributionSign = 1.0;
    int Phase = RemoveVolume [ 0 ]->PosNormalMatID;

    for ( TriangleType *t : RemoveVolume ) {
        double v1x = t->Vertices [ 0 ]->get_c(0);
        double v1y = t->Vertices [ 0 ]->get_c(1);
        double v1z = t->Vertices [ 0 ]->get_c(2);

        double v2x = t->Vertices [ 1 ]->get_c(0);
        double v2y = t->Vertices [ 1 ]->get_c(1);
        double v2z = t->Vertices [ 1 ]->get_c(2);

        double v3x = t->Vertices [ 2 ]->get_c(0);
        double v3y = t->Vertices [ 2 ]->get_c(1);
        double v3z = t->Vertices [ 2 ]->get_c(2);

        ContributionSign = ( t->PosNormalMatID == Phase ) ? 1 : -1;

        double Vc = 1. / 6. * ( v3x * ( v1y * v2z - v2y * v1z ) + v3y * ( v2x * v1z - v1x * v2z ) + v3z * ( v1x * v2y - v1y * v2x ) ) * ContributionSign;
        V = V + Vc;
    }

    for ( TriangleType *t : *NewTriangles ) {
        double v1x = t->Vertices [ 0 ]->get_c(0);
        double v1y = t->Vertices [ 0 ]->get_c(1);
        double v1z = t->Vertices [ 0 ]->get_c(2);

        double v2x = t->Vertices [ 1 ]->get_c(0);
        double v2y = t->Vertices [ 1 ]->get_c(1);
        double v2z = t->Vertices [ 1 ]->get_c(2);

        double v3x = t->Vertices [ 2 ]->get_c(0);
        double v3y = t->Vertices [ 2 ]->get_c(1);
        double v3z = t->Vertices [ 2 ]->get_c(2);

        ContributionSign = ( t->PosNormalMatID == Phase ) ? 1 : -1;

        double Vc = 1. / 6. * ( v3x * ( v1y * v2z - v2y * v1z ) + v3y * ( v2x * v1z - v1x * v2z ) + v3z * ( v1x * v2y - v1y * v2x ) ) * ContributionSign;
        V = V - Vc;
    }

    error = fabs(V);

    if ( error > TOL_COL_MAXVOLUMECHANGE ) {
        return FC_AREACHANGETOOLARGE;
    } else {
        return FC_OK;
    }
}

FC_MESH MeshManipulations :: CheckCoarsenChord(EdgeType *EdgeToCollapse, VertexType *RemoveVertex, VertexType *SaveVertex)
{
    // If both vertices are located on an edge, proceed wih check. If not, collapsing is ok since this is not a chord.
    if ( !( SaveVertex->IsPhaseEdgeVertex() && RemoveVertex->IsPhaseEdgeVertex() ) ) {
        return FC_OK;
    }

    // If RemoveVertex is fixed, prevent collapsing
    if ( RemoveVertex->IsFixedVertex() ) {
        return FC_FIXEDVERTEX;
    }

    // If vertices are located on different PhaseEdges, collapsing is not ok

    // Here we know that RemoveVertex only contains one PhaseEdge and both vertices are located on the chord.
    PhaseEdge *RemoveVertexPhaseEdge = RemoveVertex->PhaseEdges.at(0);
    if ( std :: find(SaveVertex->PhaseEdges.begin(), SaveVertex->PhaseEdges.end(), RemoveVertexPhaseEdge) == SaveVertex->PhaseEdges.end() ) {
        return FC_VERTICESONDIFFERENTSHAPES;
    }

    // Check if chord changes too much...
    std :: array< VertexType *, 2 >NewEdge;
    std :: array< VertexType *, 2 >OtherEdge;


    // Find the other edge connected to RemoveVertex (the one connected to the same PhaseEdge)
    std :: vector< VertexType * >ConnectedVertices = RemoveVertexPhaseEdge->GiveVerticesConnectedToVertex(RemoveVertex);
    if ( ConnectedVertices.size() == 1 ) {
        // This should be ok if the edge is small. This simply removes the (very small) chord.
        // Note that this can imply duplicate triangles and it is quite cumbersome to solve this. Thus, we push this feature forward
        double EdgeLength = EdgeToCollapse->GiveLength();
        LOG("EdgeLength = %f\n", EdgeLength);
        if (EdgeLength>.25) {
            return FC_CHORD;
        }
    }
    VertexType *OtherVertex = ( ConnectedVertices [ 0 ] == SaveVertex ) ? ConnectedVertices [ 1 ] : ConnectedVertices [ 0 ];
    NewEdge = { { SaveVertex, OtherVertex } };
    OtherEdge = { { RemoveVertex, OtherVertex } };

    // Compute normalized vectors
    std :: array< double, 3 >NewNormal = ComputeNormalizedVector( NewEdge.at(0), NewEdge.at(1) );
    std :: array< std :: array< double, 3 >, 2 >Normals = { { ComputeNormalizedVector( OtherEdge.at(0), OtherEdge.at(1) ),
                                                              ComputeNormalizedVector( EdgeToCollapse->Vertices.at(0), EdgeToCollapse->Vertices.at(1) ) } };
    double MaxAngle = 0.0;

    for ( int i : { 0, 1 } ) {
        double Alpha = std :: acos(NewNormal.at(0) * Normals [ i ] [ 0 ] + NewNormal.at(1) * Normals [ i ] [ 1 ] + NewNormal.at(2) * Normals [ i ] [ 2 ]);
        Alpha = std :: min( std :: abs(Alpha), std :: abs(Alpha - 3.141592) );
        if ( Alpha > MaxAngle ) {
            MaxAngle = Alpha;
        }
    }

    if ( MaxAngle > TOL_COL_CHORD_MAXNORMALCHANGE ) {
        return FC_CHORD;
    }

    return FC_OK;
}

std :: vector< VertexType * >MeshManipulations :: FindIndependentSet()
{
    std :: vector< VertexType * >IndepSet;

    // Add all fixed vertices (Vertices shared between more than 2 phases)
    for ( VertexType *v : this->Vertices ) {
        if ( v->IsFixedVertex() ) {
            IndepSet.push_back(v);
        }
    }

    // Create vectors of vertices on edges and faces respectively
    std :: vector< VertexType * >PhaseEdgeVertices;
    std :: vector< VertexType * >FaceVertices;
    for ( VertexType *v : this->Vertices ) {
        if ( v->IsPhaseEdgeVertex() ) {
            PhaseEdgeVertices.push_back(v);
        } else {
            FaceVertices.push_back(v);
        }
    }


    for ( std :: vector< VertexType * >vl : { PhaseEdgeVertices, FaceVertices } ) {
        for ( VertexType *v : vl ) {
            if ( std :: find(IndepSet.begin(), IndepSet.end(), v) == IndepSet.end() ) {
                bool AddToSet = true;

                // Add v to IndepSet if no neighbour is in the set
                for ( EdgeType *e : v->Edges ) {
                    VertexType *w;
                    if ( e->Vertices [ 0 ] == v ) {
                        w = e->Vertices [ 1 ];
                    } else {
                        w = e->Vertices [ 0 ];
                    }
                    if ( std :: find(IndepSet.begin(), IndepSet.end(), w) != IndepSet.end() ) {
                        AddToSet = false;
                    }
                }

                if ( AddToSet ) {
                    IndepSet.push_back(v);
                }
            }
        }
    }

    return IndepSet;
}

int MeshManipulations :: FlipAll()
{
    this->UpdateLongestEdgeLength();

    int flipcount = 0;
    int i = 0;
    bool edgeflipped = true;
    while ( edgeflipped ) {
        edgeflipped = false;
        size_t j = 0;
        while ( j < this->Edges.size() ) {
            EdgeType *e = this->Edges [ j ];
            LOG("Flip edge iteration %u: edge @%p (%u, %u)\n", i, e, e->Vertices [ 0 ]->ID, e->Vertices [ 1 ]->ID);
            if ( this->FlipEdge(e) == FC_OK ) {
                flipcount++;
                edgeflipped = true;
            } else {
                j++;
            }
            i++;
        }
    }
    return flipcount;
}
void MeshManipulations :: UpdateLongestEdgeLength()
{
    this->LongestEdgeLength=0;
    for (EdgeType *e: this->Edges) {
        double ThisLength = e->GiveLength();
        this->LongestEdgeLength = std::max(ThisLength, this->LongestEdgeLength);
    }
}

void MeshManipulations :: CoarsenMesh()
{
    STATUS("Coarsen mesh\n", 0);

    bool CoarseningOccurs = true;
    int iter = 0;

#if TEST_MESH_FOR_EACH_COARSENING_ITERATION
    TetGenCaller Generator;
    Generator.Mesh = this;
#endif

#if EXPORT_MESH_COARSENING
    dooutputlogmesh(* this, "/tmp/Coarsening_%u.vtp", 0);
#endif

#if EXPORT_MESH_COARSENING
    int MeshIndex = 0;
#endif

    while ( CoarseningOccurs ) {
        CoarseningOccurs = false;

        // Find independent sets. Vertices belonging to IndepSet cannot be collapsed (according to de Cougnt)
        std :: vector< VertexType * >IndepSet = FindIndependentSet();
        for ( VertexType *v : this->Vertices ) {
            v->tag = 0;
        }
        for ( VertexType *v : IndepSet ) {
            v->tag = 1;
        }

        unsigned int i = 0;

        while ( i < this->Edges.size() ) {
            STATUS( "%c[2K\rCoarsening iteration %u, collapse edge %u (%u)", 27, iter, i, this->Edges.size() );
            fflush(stdout);

            EdgeType *e = this->Edges.at(i);
            // Try to collapse vertices on current edge
            std :: array< VertexType *, 2 >EdgeVertices = { { e->Vertices [ 0 ], e->Vertices [ 1 ] } };
            int vi = 0;
            for ( VertexType *v : EdgeVertices ) {
                // If vertex v is not in the set of independent vertices, try to collapse
                if ( std :: find(IndepSet.begin(), IndepSet.end(), v) == IndepSet.end() ) {
                    if ( this->CollapseEdge(e, vi) == FC_OK ) {
                        this->UpdateLongestEdgeLength();
                        CoarseningOccurs = true;
#if EXPORT_MESH_COARSENING
                        this->ExportSurface(strfmt("/tmp/Coarseningp_%u.simple", MeshIndex), FT_SIMPLE);
                        dooutputlogmesh(* this, "/tmp/Coarsening_%u.vtp", MeshIndex++);
#endif

#if TEST_MESH_FOR_EACH_COARSENING_ITERATION
                        Generator.TestMesh();
#endif
                        break;
                    }
                }
                vi++;
            }
            i++;
        }
        iter++;


#if SANITYCHECK == 1
        for ( TriangleType *t1 : this->Triangles ) {
            for ( TriangleType *t2 : this->Triangles ) {
                if ( t1 != t2 ) {
                    bool ispermutation = std :: is_permutation( t1->Vertices.begin(), t1->Vertices.end(), t2->Vertices.begin() );
                    if ( ispermutation ) {
                        STATUS("\nDuplicate triangles at iteration %u, failcount %u\n", iter, failcount);
                        throw 0;
                    }
                }
            }
        }
#endif
    }
    STATUS("\n", 0);

    bool DoesCollapse = true;
    int CleanupIteration = 0;
    int CleanupCount=0;

    std::vector<FC_MESH> Reasons;

/*    TOL_COL_CHORD_MAXNORMALCHANGE = TOL_COL_CHORD_MAXNORMALCHANGE*2;
    TOL_COL_MAXERROR = TOL_COL_MAXERROR*2;
    TOL_COL_MAXVOLUMECHANGE = TOL_COL_MAXVOLUMECHANGE*2; */

    STATUS("Clean up triangles with very small angles\n", 0);
    while (DoesCollapse) {
        DoesCollapse = false;
        Reasons.clear();

        size_t i=0;

        while (i<this->Triangles.size()) {
            STATUS( "%c[2K\rCleanup iteration %u, edge %u (%u)", 27, iter, i, this->Edges.size() );
            fflush(stdout);

            TriangleType *t=this->Triangles.at(i);

            int VertexIndex;
            double Angle = t->GiveSmallestAngle(&VertexIndex);
            if (Angle < 0.2) {
                // If angle is small, try to collapse the shortest edge
                std::array <EdgeType *, 3> tEdges = t->GiveEdges();
                std::sort(tEdges.begin(), tEdges.end(), CompareEdgeLength);

                EdgeType *e=tEdges[0];

                FC_MESH r = this->CollapseEdge(e, 0);
                if (!(r==FC_OK)) {
                    r = this->CollapseEdge(e, 1);
                }
                Reasons.push_back(r);
                if (r==FC_OK) {
                    DoesCollapse = true;
                    CleanupCount++;
                }
                LOG("Collapse result:%u\n", r);

            }

            i++;
        }
        CleanupIteration++;
        std::sort(Reasons.begin(), Reasons.end());
        this->FlipAll();
    }
    STATUS("\nCleanup removed %u triangles\n", CleanupCount);

}

}
