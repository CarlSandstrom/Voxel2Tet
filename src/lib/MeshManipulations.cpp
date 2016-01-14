#include "MeshManipulations.h"
#include "PhaseEdge.h"

namespace voxel2tet
{

MeshManipulations::MeshManipulations(BoundingBoxType BoundingBox) : MeshData(BoundingBox)
{

}

void MeshManipulations :: RemoveDegenerateTriangles()
{
    STATUS("Remove degenerate triangles...\n", 0);

    // Find triangles with zero-area
    for (unsigned int i=0; i<this->Triangles.size(); i++) {
        TriangleType *t = this->Triangles.at(i);

        int index;
        double LargestAngle = t->GiveLargestAngle(&index);

        LOG("Triangle #%i@%p  has largest angle %f\n", i, t, LargestAngle);
        if (LargestAngle > 3) {
            // Note that this can occur by either a very small angle and/or two vertices located very close to each other
            LOG("Triangle #%i@%p has a too large largest angle\n", i, t);

            // Find largest angle and the associated index

            // Find opposite edge and flip that edge
            std::array<EdgeType *, 3> Edges = t->GiveEdges();
            EdgeType *OppositeEdge = NULL;
            for (EdgeType *e: Edges) {
                if ( (e->Vertices[0]!=t->Vertices[index]) & (e->Vertices[1]!=t->Vertices[index])) {
                    OppositeEdge = e;
                    break;
                }
            }

            if (OppositeEdge==NULL) {
                LOG ("Something unexpected occured\n", 0);
            }

            if (!this->FlipEdge(OppositeEdge)) {
                STATUS("The longest edge of triangle %u failed to flip\n", i);
            }
        }
    }
}

bool MeshManipulations :: GetFlippedEdgeData(EdgeType *EdgeToFlip, EdgeType *NewEdge, std::array<TriangleType*, 2> *NewTriangles)
{
    LOG("Get flipped edge data for edge %p\n", EdgeToFlip);

    std::vector<TriangleType *> EdgeTriangles = EdgeToFlip->GiveTriangles();

    if (EdgeTriangles.size()!=2) {
        LOG("Unable to flip edge. To many or only one triangle connected\n", 0);
        return false;
    }

    for (int i=0; i<2; i++) std::sort(EdgeTriangles.at(i)->Vertices.begin(), EdgeTriangles.at(i)->Vertices.end());

    std::vector<VertexType *> NewEdgeVertices;

    std::set_symmetric_difference(EdgeTriangles.at(0)->Vertices.begin(), EdgeTriangles.at(0)->Vertices.end(),
                                  EdgeTriangles.at(1)->Vertices.begin(), EdgeTriangles.at(1)->Vertices.end(), std::back_inserter(NewEdgeVertices));

    // Create triangles
    for (int i=0; i<2; i++) {
        TriangleType* t=new TriangleType;
        if (i==0) {
            t->Vertices[0] = NewEdgeVertices[1];
            t->Vertices[1] = NewEdgeVertices[0];
            t->Vertices[2] = EdgeToFlip->Vertices[i];
        } else {
            t->Vertices[0] = NewEdgeVertices[0];
            t->Vertices[1] = NewEdgeVertices[1];
            t->Vertices[2] = EdgeToFlip->Vertices[i];
        }
        t->InterfaceID = EdgeTriangles.at(0)->InterfaceID;
        t->NegNormalMatID = -1; // TODO: Keep track of materials
        t->PosNormalMatID = -1;
        t->UpdateNormal();
        NewTriangles->at(i) = t;
    }

    // Create edge
    NewEdge->Vertices[0] = NewEdgeVertices[0];
    NewEdge->Vertices[1] = NewEdgeVertices[1];

    return true;

}

bool MeshManipulations :: FlipEdge(EdgeType *Edge)
{
    LOG ("\tFlip edge %p (%u, %u)\n", Edge, Edge->Vertices[0]->ID, Edge->Vertices[1]->ID);

    std::vector<TriangleType *> EdgeTriangles = Edge->GiveTriangles();
    if (EdgeTriangles.size()!=2) {
        LOG("Unable to flip edge. To many or only one triangle connected\n", 0);
        return false;
    }

    std::array<TriangleType*, 2> NewTriangles;
    EdgeType NewEdge;

    this->GetFlippedEdgeData(Edge, &NewEdge, &NewTriangles);

    if (!this->CheckFlipNormal(&EdgeTriangles, NewTriangles)) {
        LOG("Changes in normal direction prevents flipping\n", 0);
        for (TriangleType *t: NewTriangles) {
            delete t;
        }
        return false;
    }

    // Ensure that minimal angle is improved if we continue
    double minAngleCurrent = std::min(EdgeTriangles[0]->GiveSmallestAngle(), EdgeTriangles[1]->GiveSmallestAngle());
    double minAngleNew = std::min(NewTriangles[0]->GiveSmallestAngle(), NewTriangles[1]->GiveSmallestAngle());

    if (minAngleNew<minAngleCurrent) {
        LOG("New minimal angles worse than current (New: %f, Current: %f). Prevent flipping.\n", minAngleNew, minAngleCurrent);
        for (TriangleType *t: NewTriangles) {
            delete t;
        }
        return false;
    }

    // Check change in area of region. Should be small
    double CurrentArea = EdgeTriangles.at(0)->GiveArea()+EdgeTriangles.at(1)->GiveArea();
    double NewArea = NewTriangles.at(0)->GiveArea()+NewTriangles.at(1)->GiveArea();

    if (std::fabs(CurrentArea-NewArea)>1e-4) { // TODO: Use variable instead of fixed value
        LOG("The combined area of the triangles changes too much. Prevent flipping.\n", 0);
        for (TriangleType *t: NewTriangles) {
            delete t;
        }
        return false;
    }

    // Check if any of the new triangles already exist
    std::vector<TriangleType *> ConnectedTriangles;
    for (VertexType *v: NewEdge.Vertices) {
        for (TriangleType *t: v->Triangles) {
            ConnectedTriangles.push_back(t);
        }
    }

    std::sort(ConnectedTriangles.begin(), ConnectedTriangles.end());

    for (unsigned int i=0; i<ConnectedTriangles.size()-1; i++) {
        if (ConnectedTriangles.at(i) == ConnectedTriangles.at(i+1)) {
            LOG("Edge already exists!\n", 0);
            return false;
        }
    }


    LOG("Flip edge!\n", 0);
    // Update mesh data

    // Update triangles
/*    for (int i=0; i<2; i++) {
        for (int j=0; j<3; j++) {
            EdgeTriangles.at(i)->Vertices[j] = NewTriangles.at(i)->Vertices[j];
        }
    }*/
    // Remove triangles from verices
    for (int i: {0, 1}) {
        for (TriangleType *t: EdgeTriangles) {
            Edge->Vertices[i]->RemoveTriangle(t);
        }
    }

    // Remove edge from vertices
    for (int i: {0, 1}) Edge->Vertices[i]->RemoveEdge(Edge);



    // Add new edge to vertices
    for (int i: {0, 1}) {
        Edge->Vertices[i] = NewEdge.Vertices[i];
        Edge->Vertices[i]->AddEdge(Edge);
    }

    // Add new triangles list (and thus also to vertices)
    for (TriangleType *t: NewTriangles) {
        this->AddTriangle(t);
    }

    // 5. Free old triangles
    for (TriangleType *t: EdgeTriangles) {
        this->RemoveTriangle(t);
    }

    return true;

}

bool MeshManipulations :: CheckFlipNormal(std::vector<TriangleType*> *OldTriangles, std::array<TriangleType *, 2> NewTriangles)
{
    double MaxAngle = 0.0;

    for (TriangleType *OldTriangle: *OldTriangles) {
        for (unsigned int j=0; j<NewTriangles.size(); j++) {
            TriangleType *NewTriangle = NewTriangles.at(j);
            std::array<double, 3> OldNormal = OldTriangle->GiveUnitNormal(); // TODO: Move to first for loop
            std::array<double, 3> NewNormal = NewTriangle->GiveUnitNormal();
            double angle1 = std::acos( OldNormal[0]*NewNormal[0] + OldNormal[1]*NewNormal[1] + OldNormal[2]*NewNormal[2]);
            double angle2 = std::acos( -(OldNormal[0]*NewNormal[0] + OldNormal[1]*NewNormal[1] + OldNormal[2]*NewNormal[2]) );
            if (std::min(angle1, angle2)>MaxAngle) {
                MaxAngle = std::min(angle1, angle2);
            }
        }
    }

    if (MaxAngle > (10*2*3.1415/360))
        return false;
    else
        return true;

}

bool MeshManipulations :: CollapseEdge(EdgeType *EdgeToCollapse, int RemoveVertexIndex)
{
    // TODO: In its current setting, this procedure only checks if collapsing is ok from the current topology. It should compare to the original topology, otherwise degeneration can occur gradually

    LOG("Collapse edge %p (%u, %u) by removing vertex %u\n", EdgeToCollapse, EdgeToCollapse->Vertices[0]->ID,
            EdgeToCollapse->Vertices[1]->ID, RemoveVertexIndex);

    // Cannot remove a fixed vertex
    if (EdgeToCollapse->Vertices[RemoveVertexIndex]->IsFixedVertex()) return false;

    int SaveVertexIndex = (RemoveVertexIndex == 0) ? 1 : 0;

    VertexType *SaveVertex = EdgeToCollapse->Vertices[SaveVertexIndex];

    // Create new triangles. These are create by moving RemoveVertex to the other end of the edge and remove the 0-area triangles
    std::vector<TriangleType*> TrianglesToRemove = EdgeToCollapse->GiveTriangles();
    LOG("Connected triangle IDs: %u, %u\n", TrianglesToRemove.at(0)->ID, TrianglesToRemove.at(1)->ID);
    std::vector<TriangleType*> ConnectedTriangles = EdgeToCollapse->Vertices[RemoveVertexIndex]->Triangles;

    std::sort(TrianglesToRemove.begin(), TrianglesToRemove.end());
    std::sort(ConnectedTriangles.begin(), ConnectedTriangles.end());

    std::vector<TriangleType*> TrianglesToSave;
    std::set_difference(ConnectedTriangles.begin(), ConnectedTriangles.end(),
                        TrianglesToRemove.begin(), TrianglesToRemove.end(), std::inserter( TrianglesToSave, TrianglesToSave.begin() ) );

    std::vector<TriangleType*> NewTriangles;
    for (TriangleType* t: TrianglesToSave) {
        t->UpdateNormal();
        TriangleType* NewTriangle = new TriangleType;

        // TODO: Keep track of materials
        NewTriangle->InterfaceID = t->InterfaceID;
        NewTriangle->ID = -1;

        for (int i=0; i<3; i++) {
            NewTriangle->Vertices[i] = t->Vertices[i];
        }
        NewTriangle->UpdateNormal();

        for (int i=0; i<3; i++) {
            if (t->Vertices[i]==EdgeToCollapse->Vertices[RemoveVertexIndex]) {
                NewTriangle->Vertices[i] = EdgeToCollapse->Vertices[SaveVertexIndex];
            } else {
                NewTriangle->Vertices[i] = t->Vertices[i];
            }
        }

        NewTriangle->UpdateNormal();
        NewTriangles.push_back(NewTriangle);
    }

    bool DoCollapse = true;

    // Check if NewTriangles are good replacements
    LOG("Check validity of new normals...\n", 0);
    if (!this->CheckCoarsenNormal(&TrianglesToSave, &NewTriangles)) {
        LOG(" - Check failed\n", 0);
        DoCollapse = false;
    }

    LOG("Check validity of new chord...\n", 0);
    if (!this->CheckCoarsenChord(EdgeToCollapse, EdgeToCollapse->Vertices[RemoveVertexIndex], EdgeToCollapse->Vertices[SaveVertexIndex])) {
        LOG(" - Check failed\n", 0);
        DoCollapse = false;
    }

    LOG("Check area of new triangles...\n", 0);
    for (TriangleType *t: NewTriangles) {
        if (t->GiveArea()<1e-7) {
            LOG(" - Check failed\n", 0);
            DoCollapse = false;
        }
    }

    LOG("Check cmmbined area of new vs current triangles\n", 0);
    double NewArea=0, CurrentArea=0;
    for (TriangleType *t: NewTriangles) NewArea = NewArea + t->GiveArea();
    for (TriangleType *t: ConnectedTriangles) CurrentArea = CurrentArea + t->GiveArea();
    if (fabs(NewArea-CurrentArea) > 1e-2) {
        LOG("Area change significant. Prevent collapse\n", 0);
        DoCollapse = false;
    }

    if (!DoCollapse) {
        for (TriangleType *t: NewTriangles) {
            delete t;
        }
        return false;
    }

    LOG ("Checks ok. Proceed\n", 0);

    // Update triangulation

    // Find edges to remove (all edges connected to RemoveVertex and in any triangle in TrianglesToRemove)
    std::vector<EdgeType *> EdgesToRemove;
    std::vector<EdgeType *> RemoveVertexEdges = EdgeToCollapse->Vertices.at(RemoveVertexIndex)->Edges;

    std::vector<EdgeType *> TriangleToRemoveEdges;
    for (TriangleType *t: TrianglesToRemove) {
        for (EdgeType *e: t->GiveEdges()) {
            TriangleToRemoveEdges.push_back(e);
        }
    }

    std::sort(RemoveVertexEdges.begin(), RemoveVertexEdges.end());
    std::sort(TriangleToRemoveEdges.begin(), TriangleToRemoveEdges.end());

    std::set_intersection(RemoveVertexEdges.begin(), RemoveVertexEdges.end(),
                          TriangleToRemoveEdges.begin(), TriangleToRemoveEdges.end(), std::back_inserter(EdgesToRemove));

    // Update edges
    std::vector<EdgeType*> ConnectedEdges = EdgeToCollapse->Vertices.at(RemoveVertexIndex)->Edges;
    for (EdgeType *e: ConnectedEdges) {
        if (e!=EdgeToCollapse) {
            for (int i: {0, 1}) {
                // "Move" vertex
                if (e->Vertices[i]==EdgeToCollapse->Vertices.at(RemoveVertexIndex)) {
                    e->Vertices[i]->RemoveEdge(e);
                    e->Vertices[i] = EdgeToCollapse->Vertices.at(SaveVertexIndex);
                    e->Vertices[i]->AddEdge(e);
                }
            }
        }
    }

    // Remove triangles
    for (TriangleType* t: ConnectedTriangles) {
        this->RemoveTriangle(t);
    }

    // Remove edges
    for (EdgeType* e: EdgesToRemove) {
        this->RemoveEdge(e);
        //ConnectedEdges.erase( std::remove(ConnectedEdges.begin(), ConnectedEdges.end(), e), ConnectedEdges.end());
    }

    // Add new triangles
    for (TriangleType* t: NewTriangles) {
        this->AddTriangle(t);
    }

    ConnectedEdges = SaveVertex->Edges;
    for (EdgeType *e: ConnectedEdges) {
        FlipEdge(e);
    }

    return true;

}

bool MeshManipulations :: CheckCoarsenNormal(std::vector<TriangleType*> *OldTriangles, std::vector<TriangleType*> *NewTriangles)
{
    // Here, we use that the order of the old and new triangles are the same in both lists. I.e. NewTriangle[i] is the same triangle as OldTriangles[i] except that vD is replaced by vR.

    double MaxAngle = 0;

    for (unsigned int i=0; i<OldTriangles->size(); i++) {

        std::array<double, 3> OldNormal = OldTriangles->at(i)->GiveUnitNormal();
        std::array<double, 3> NewNormal = NewTriangles->at(i)->GiveUnitNormal();

        double angle1 = std::acos( OldNormal[0]*NewNormal[0] + OldNormal[1]*NewNormal[1] + OldNormal[2]*NewNormal[2]);
        double angle2 = std::acos( -(OldNormal[0]*NewNormal[0] + OldNormal[1]*NewNormal[1] + OldNormal[2]*NewNormal[2]) );

        double oldarea  = OldTriangles->at(i)->GiveArea();
        double newarea  = NewTriangles->at(i)->GiveArea();

        LOG("angle1=%f,\tangle2=%f\tOldArea=%f\tNewArea=%f\n", angle1, angle2, oldarea, newarea);

        if (std::min(angle1, angle2)>MaxAngle) {
            MaxAngle = std::min(angle1, angle2);
        }

        if (angle1 > (170*2*3.1415/360)) {
            LOG("Should not be collapsed...\n", 0);
            return false;
        }


    }

    if (MaxAngle > (10*2*3.1415/360))
        return false;
    else
        return true;
}

bool MeshManipulations :: CheckCoarsenChord(EdgeType *EdgeToCollapse, VertexType* RemoveVertex, VertexType* SaveVertex)
{
    // If both vertices are located on an edge, proceed wih check. If not, collapsing is ok.
    if (!(SaveVertex->IsEdgeVertex() && RemoveVertex->IsEdgeVertex())) {
        return true;
    }

    // If RemoveVertex is fixed, prevent collapsing
    if (RemoveVertex->IsFixedVertex()) {
        return false;
    }

    // If vertices are located on different PhaseEdges, collapsing is not ok

    // Here we know that RemoveVertex only contains one PhaseEdge and both vertices are located on edges.
    PhaseEdge* RemoveVertexPhaseEdge = RemoveVertex->PhaseEdges.at(0);
    bool SamePhaseEdge = false;
    for (PhaseEdge* pe: SaveVertex->PhaseEdges) {
        if (pe == RemoveVertexPhaseEdge) {
            SamePhaseEdge = true;
            break;
        }
    }

    if (!SamePhaseEdge) return false;

    // Check if chord changes too much...
    std::array<VertexType*, 2> NewEdge;
    std::array<VertexType*, 2> *OtherEdge = NULL;

    // Find the other edge connected to RemoveVertex (the one connected to the same PhaseEdge)
    for (std::array<VertexType*, 2> e: RemoveVertexPhaseEdge->EdgeSegments) {
        int OtherIndex = -1;
        if (e.at(0) == RemoveVertex) OtherIndex = 1;
        if (e.at(1) == RemoveVertex) OtherIndex = 0;

        if ( OtherIndex != -1 ) {
            // Is 'e' the edge that will be collapsed?
            if ( ( (e.at(0)==EdgeToCollapse->Vertices.at(0)) && (e.at(1)==EdgeToCollapse->Vertices.at(1)) ) ||
                 ( (e.at(1)==EdgeToCollapse->Vertices.at(0)) && (e.at(0)==EdgeToCollapse->Vertices.at(1)) ) ) {
                // This is the current edge. Pass
            } else {
                NewEdge.at(0) = SaveVertex;
                NewEdge.at(1) = e.at(OtherIndex);
                OtherEdge = &e;
                break;
            }
        }
    }

    if (OtherEdge==NULL) {
        LOG("OtherEdge not found. \n", 0);
        return false;
    }

    // Compute normalized vectors
    std::array<double, 3> NewNormal = ComputeNormalizedVector(NewEdge.at(0), NewEdge.at(1));
    std::array<std::array<double, 3>, 2> Normals = {ComputeNormalizedVector(OtherEdge->at(0), OtherEdge->at(1)),
                                                    ComputeNormalizedVector(EdgeToCollapse->Vertices.at(0), EdgeToCollapse->Vertices.at(1)) };
    double MaxAngle = 0.0;

    for (int i: {0, 1}) {
        double Alpha = std::acos(NewNormal.at(0)*Normals[i][0] + NewNormal.at(1)*Normals[i][1] + NewNormal.at(2)*Normals[i][2] );
        Alpha = std::min(std::abs(Alpha), std::abs(Alpha-3.141592));
        if (Alpha>MaxAngle) MaxAngle=Alpha;
    }

    if (MaxAngle > 45*3.141593/360) {
        return false;
    }

    return true;

}

std::vector<VertexType *> MeshManipulations :: FindIndependentSet()
{
    std::vector<VertexType *> IndepSet;

    // Add all fixed vertices
    for (VertexType *v: this->Vertices) {
        if (v->IsFixedVertex()) {
            IndepSet.push_back(v);
        }
    }

    // Create vectors of vertices on edges and faces
    std::vector<VertexType *> EdgeVertices;
    std::vector<VertexType *> FaceVertices;
    for (VertexType *v: this->Vertices) {
        if (v->IsEdgeVertex()) {
            EdgeVertices.push_back(v);
        } else {
            FaceVertices.push_back(v);
        }
    }

    for (std::vector<VertexType *> vl: {EdgeVertices, FaceVertices}) {
        for (VertexType *v: vl) {

            if (std::find(IndepSet.begin(), IndepSet.end(), v) == IndepSet.end()) {

                bool AddToSet = true;

                // Add v to IndepSet if no neighbour is in the set
                for (EdgeType *e: v->Edges) {
                    VertexType *w;
                    if (e->Vertices[0] == v) {
                        w = e->Vertices[1];
                    } else {
                        w = e->Vertices[0];
                    }
                    if (std::find(IndepSet.begin(), IndepSet.end(), w) != IndepSet.end()) AddToSet=false;
                }

                if (AddToSet) {
                    IndepSet.push_back(v);
                }

            }
        }
    }

    return IndepSet;
}

bool MeshManipulations :: FlipAll()
{
    int i=0;
    for (EdgeType *e: this->Edges) {
        LOG ("Flip edge iteration %u: edge @%p (%u, %u)\n", i, e, e->Vertices[0]->ID, e->Vertices[1]->ID);
        this->FlipEdge(e);
        std::ostringstream FileName;
        FileName << "/tmp/TestFlip_" << i << ".vtp";
        //this->ExportVTK( FileName.str() );
        i++;
    }
}

bool MeshManipulations :: CoarsenMeshImproved()
{

    STATUS("Coarsen mesh\n", 0);

    bool CoarseningOccurs=true;
    int iter=0;
    int failcount;

    while (CoarseningOccurs) {

        STATUS("%c[2K\rCoarsening iteration %u", 27, iter);
        fflush(stdout);

        CoarseningOccurs = false;
        std::vector<VertexType*> IndepSet = FindIndependentSet();
        failcount=0;
        int i=0;
        while (i<this->Edges.size()) {
            EdgeType *e = this->Edges.at(i);
            if (e->GiveLength()<1e10) { // TODO: Use argument here
                bool BothInSet = (std::find(IndepSet.begin(), IndepSet.end(), e->Vertices[0]) == IndepSet.end()) && (std::find(IndepSet.begin(), IndepSet.end(), e->Vertices[1]) == IndepSet.end());
                if (!BothInSet) {
                    int D = 0;
                    if (std::find(IndepSet.begin(), IndepSet.end(), e->Vertices[D]) != IndepSet.end()) {
                        D = 1;
                    }
                    if ( (iter==2616) & (failcount==731)) {
                        this->ExportVTK( "/tmp/beforecollapse.vtp" );
                    }
                    bool CoarseOk = this->CollapseEdge(e, D);
                    if (!CoarseOk) {
                        D = 1 ? 0 : 1;
                        CoarseOk = this->CollapseEdge(e, D);
                    }
                    if (CoarseOk) {
                        CoarseningOccurs = true;
                        iter++;
                    } else {
                        failcount++;
                    }

                }
            }
            i++;
        }
#if SANITYCHECK == 1
            for (TriangleType *t1: this->Triangles) {
                for (TriangleType *t2: this->Triangles) {
                    if (t1!=t2) {
                        bool ispermutation = std::is_permutation(t1->Vertices.begin(), t1->Vertices.end(), t2->Vertices.begin());
                        if (ispermutation) {
                            STATUS ("\nDuplicate triangles at iteration %u, failcount %u\n", iter, failcount);
                            throw 0;
                        }
                    }
                }
            }
#endif
    }
    STATUS("\n",0);

    //return true;

    for (TriangleType *t: this->Triangles) {
        double A = t->GiveArea();
        if (A<1e-3) {
            STATUS("Triangle %u has a too small area (A=%f)\n", t->ID, A);
            for (EdgeType *e: t->GiveEdges()) {
                bool CollapseOk = this->CollapseEdge(e, 0);
                if (!CollapseOk) {
                    CollapseOk = this->CollapseEdge(e, 1);
                }
                if (CollapseOk) {
                    break;
                } else {
                    STATUS("Collapse failed\n", 0);
                }
            }
        }
    }

    return true;
}

bool MeshManipulations :: CoarsenMesh()
{

    this->ExportVTK("/tmp/TestCoarsen_0.vtp");
    int iter=1;
    bool EdgeCollapsed = true;
    while (EdgeCollapsed) {
        EdgeCollapsed=false;
        int failcount=0;

        LOG ("Start iteration %u ===================== \n", iter);

//        int i=0;
//        for (EdgeType* e: this->Edges) {
//            if ( (e->Vertices.at(0)->ID==8) & (e->Vertices.at(1)->ID==1)) {
//                LOG("Edge found at index %u! ***********\n ", i);
//                break;
//            }
//            i++;
//        }

        for (EdgeType* e: this->Edges) {

            LOG("Collapse iteration %u (%u, %u), failcount=%u\n", iter, e->Vertices.at(0)->ID, e->Vertices.at(1)->ID, failcount);

            if (!this->CollapseEdge(e,0)) {
                EdgeCollapsed=this->CollapseEdge(e,1);
            } else {
                EdgeCollapsed=true;
            }

            if (EdgeCollapsed) {
                LOG("Collapse success!\n", iter, failcount);
                std::ostringstream FileName;
                FileName << "/tmp/TestCoarsen_" << iter << ".vtp";
                //this->ExportVTK( FileName.str() );
                iter++;
                break;
            } else {
                LOG ("Collapse failed...\n", 0);
                failcount++;
            }
        }
    }
    return true;
}

}
