#include "MiscFunctions.h"
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

FC_MESH MeshManipulations :: GetFlippedEdgeData(EdgeType *EdgeToFlip, EdgeType *NewEdge, std::array<TriangleType*, 2> *NewTriangles)
{
    LOG("Get flipped edge data for edge %p\n", EdgeToFlip);

    std::vector<TriangleType *> EdgeTriangles = EdgeToFlip->GiveTriangles();

    if (EdgeTriangles.size()!=2) {
        LOG("Unable to flip edge. To many or only one triangle connected\n", 0);
        return FC_TOOMANYTRIANGLES;
    }

    for (int i=0; i<2; i++) std::sort(EdgeTriangles.at(i)->Vertices.begin(), EdgeTriangles.at(i)->Vertices.end());

    std::vector<VertexType *> NewEdgeVertices;

    std::set_symmetric_difference(EdgeTriangles.at(0)->Vertices.begin(), EdgeTriangles.at(0)->Vertices.end(),
                                  EdgeTriangles.at(1)->Vertices.begin(), EdgeTriangles.at(1)->Vertices.end(), std::back_inserter(NewEdgeVertices));

    if (NewEdgeVertices.size()==0) {
        STATUS("Edge triangles are same!\n", 0);
        throw 0;
    }

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

    return FC_OK;

}

FC_MESH MeshManipulations :: FlipEdge(EdgeType *Edge)
{
    LOG ("\tFlip edge %p (%u, %u)\n", Edge, Edge->Vertices[0]->ID, Edge->Vertices[1]->ID);

    std::vector<TriangleType *> EdgeTriangles = Edge->GiveTriangles();
    if (EdgeTriangles.size()!=2) {
        LOG("Unable to flip edge. To many or only one triangle connected\n", 0);
        return FC_TOOMANYTRIANGLES;
    }

    std::array<TriangleType*, 2> NewTriangles;
    EdgeType NewEdge;

    FC_MESH FC;

    // Check if flipping is an ok action
    FC = this->GetFlippedEdgeData(Edge, &NewEdge, &NewTriangles);
    if (FC != FC_OK) {
        return FC;
    }

    FC = this->CheckFlipNormal(&EdgeTriangles, NewTriangles);
    if (FC != FC_OK) {
        LOG("Changes in normal direction prevents flipping\n", 0);
        for (TriangleType *t: NewTriangles) {
            delete t;
        }
        return FC;
    }

    // Ensure that minimal angle is improved if we continue
    double minAngleCurrent = std::min(EdgeTriangles[0]->GiveSmallestAngle(), EdgeTriangles[1]->GiveSmallestAngle());
    double minAngleNew = std::min(NewTriangles[0]->GiveSmallestAngle(), NewTriangles[1]->GiveSmallestAngle());

    if (minAngleNew<minAngleCurrent) {
        LOG("New minimal angles worse than current (New: %f, Current: %f). Prevent flipping.\n", minAngleNew, minAngleCurrent);
        for (TriangleType *t: NewTriangles) {
            delete t;
        }
        return FC_WORSEMINANGLE;
    }

    // Check change in area of region. Should be small
    double CurrentArea = EdgeTriangles.at(0)->GiveArea()+EdgeTriangles.at(1)->GiveArea();
    double NewArea = NewTriangles.at(0)->GiveArea()+NewTriangles.at(1)->GiveArea();

    if (std::fabs(CurrentArea-NewArea)>1e-4) { // TODO: Use variable instead of fixed value
        LOG("The combined area of the triangles changes too much. Prevent flipping.\n", 0);
        for (TriangleType *t: NewTriangles) {
            delete t;
        }
        return FC_SMALLAREA;
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
            throw 0;
        }
    }

    LOG("Flip edge!\n", 0);

    // Update mesh data

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

    return FC_OK;

}

FC_MESH MeshManipulations :: CheckFlipNormal(std::vector<TriangleType*> *OldTriangles, std::array<TriangleType *, 2> NewTriangles)
{
    double MaxAngle = 0.0;

    for (TriangleType *OldTriangle: *OldTriangles) {
        for (unsigned int j=0; j<NewTriangles.size(); j++) {
            TriangleType *NewTriangle = NewTriangles.at(j);
            double NewArea = NewTriangle->GiveArea();
            // TODO: Use variable instead
            if (NewArea < 1e-8) {
                return FC_SMALLAREA;
            }
            std::array<double, 3> OldNormal = OldTriangle->GiveUnitNormal(); // TODO: Move to first for loop
            std::array<double, 3> NewNormal = NewTriangle->GiveUnitNormal();

            // Compute angle between new and old normal
            double angle1 = std::acos( OldNormal[0]*NewNormal[0] + OldNormal[1]*NewNormal[1] + OldNormal[2]*NewNormal[2]);
            double angle2 = std::acos( -(OldNormal[0]*NewNormal[0] + OldNormal[1]*NewNormal[1] + OldNormal[2]*NewNormal[2]) );
/*            if (std::min(angle1, angle2)>MaxAngle) {
                MaxAngle = std::min(angle1, angle2);
            }*/
            MaxAngle = std::max(MaxAngle, angle1);
        }
    }

    if (MaxAngle > (10*2*3.1415/360))
        return FC_NORMAL;
    else
        return FC_OK;

}

FC_MESH MeshManipulations :: CollapseEdgeTest(std::vector<TriangleType *> *TrianglesToSave, std::vector<TriangleType *> *TrianglesToRemove, std::vector<TriangleType *> *NewTriangles, EdgeType *EdgeToCollapse, int RemoveVertexIndex)
{
    int SaveVertexIndex = (RemoveVertexIndex == 0) ? 1 : 0;

    // Check if NewTriangles are good replacements
    LOG("Check validity of new normals...\n", 0);
    FC_MESH FC;

    FC = this->CheckCoarsenNormal(TrianglesToSave, NewTriangles);
    if (FC != FC_OK) {
        return FC;
    }

    LOG("Check validity of new chord...\n", 0);
    FC = this->CheckCoarsenChord(EdgeToCollapse, EdgeToCollapse->Vertices[RemoveVertexIndex], EdgeToCollapse->Vertices[SaveVertexIndex]);
    if (FC != FC_OK) {
        return FC;
    }

    LOG("Check area of new triangles...\n", 0);
    for (TriangleType *t: *NewTriangles) {
        if (t->GiveArea()<1e-7) { // TODO: Use variable
            LOG(" - Check failed\n", 0);
            return FC_SMALLAREA;
        }
    }

    LOG("Check if new triangles intersect...\n", 0);
    LOG("Remove vertex: %u\n", EdgeToCollapse->Vertices[RemoveVertexIndex]->ID);

    // Collect all point close to the centerpoint of the edge to collapse. The, form a list of all triangles connected to those points and perform check on all triangles in that list (except with triangles to remove).
    std::array<double, 3> cp = EdgeToCollapse->GiveCenterPoint();
    std::vector<VertexType *> VerticesNear = this->VertexOctreeRoot->GiveVerticesWithinSphere(cp[0], cp[1], cp[2], EdgeToCollapse->GiveLength());
    std::vector<TriangleType *> TrianglesNear;

    for (VertexType *v: VerticesNear) {
        for (TriangleType *t: v->Triangles) {
            TrianglesNear.push_back(t);
        }
    }

    std::sort(TrianglesNear.begin(), TrianglesNear.end());
    TrianglesNear.erase( std::unique(TrianglesNear.begin(), TrianglesNear.end()), TrianglesNear.end());

    for (TriangleType *t: TrianglesNear) {
        bool DoCheck = true;

        // Possible skip neighbors of new triangles
        for (TriangleType *remt: *TrianglesToRemove) { // Skip triangles that will be removed
            if (remt==t) DoCheck = false;
        }

        if (DoCheck) {
            // For this test we need to change RemoveVertex --> SaveVertex
            std::array<VertexType *, 3> NearTriVertices = t->Vertices;
            for (int i=0; i<3; i++) {
                if (NearTriVertices[i] == EdgeToCollapse->Vertices[RemoveVertexIndex]) {
                    NearTriVertices[i] = EdgeToCollapse->Vertices[SaveVertexIndex];
                    break;
                }
            }

            double V0[3] = {NearTriVertices[0]->get_c(0),NearTriVertices[0]->get_c(1),NearTriVertices[0]->get_c(2)};
            double V1[3] = {NearTriVertices[1]->get_c(0),NearTriVertices[1]->get_c(1),NearTriVertices[1]->get_c(2)};
            double V2[3] = {NearTriVertices[2]->get_c(0),NearTriVertices[2]->get_c(1),NearTriVertices[2]->get_c(2)};
            for (TriangleType *newt: *NewTriangles) {

                // Skip if triangles sharen one ore more vertices
                int sharedvertices=0;
                for (VertexType *v1: t->Vertices) {
                    if (v1 == EdgeToCollapse->Vertices[RemoveVertexIndex]) v1 = EdgeToCollapse->Vertices[SaveVertexIndex];
                    for (VertexType *v2: newt->Vertices) {
                        if (v1==v2) {
                            sharedvertices ++;
                            break;
                        }
                    }
                }

                if (sharedvertices==0) {
                    // I am a bit unsure if this is correct. Naturally two triangles with one shared vertex can intersect, but will we ever have that situation?
                    // Anyway, best would be to check if edges not members of the other triangle penetrates the surface since there is a "bug" in the used algorithm that gives a "false" true if vertices are shared
                    double U0[3] = {newt->Vertices[0]->get_c(0),newt->Vertices[0]->get_c(1),newt->Vertices[0]->get_c(2)};
                    double U1[3] = {newt->Vertices[1]->get_c(0),newt->Vertices[1]->get_c(1),newt->Vertices[1]->get_c(2)};
                    double U2[3] = {newt->Vertices[2]->get_c(0),newt->Vertices[2]->get_c(1),newt->Vertices[2]->get_c(2)};
                    int intersects = tri_tri_intersect(V0, V1, V2, U0, U1, U2);
                    if (intersects==1) {
                        return FC_TRIANGLESINTERSECT;
                    } else {
                        LOG("Triangles does not intersect\n", 0);
                    }
                } else if (sharedvertices==3) {

                    LOG ("Duplicate triangle, newt.id=%i, t.id=%i\n", newt->ID, t->ID);
                    if (t->ID!=-newt->ID) {
                        return FC_DUPLICATETRIANGLE;
                    }
                }
            }
        }
    }

    return FC_OK;
}

FC_MESH MeshManipulations :: CollapseEdge(EdgeType *EdgeToCollapse, int RemoveVertexIndex, bool PerformTesting)
{
    // TODO: In its current setting, this procedure only checks if collapsing is ok from the current topology. It should compare to the original topology, otherwise degeneration can occur gradually

    LOG("Collapse edge %p (%u, %u) by removing vertex %u\n", EdgeToCollapse, EdgeToCollapse->Vertices[0]->ID,
            EdgeToCollapse->Vertices[1]->ID, EdgeToCollapse->Vertices[RemoveVertexIndex]->ID);

    // Cannot remove a fixed vertex
    if (EdgeToCollapse->Vertices[RemoveVertexIndex]->IsFixedVertex()) return FC_FIXEDVERTEX;

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

    // NewTriangles is the updated subset of of TrianglesToSave with the removed vertex changed to the saved vertex
    std::vector<TriangleType*> NewTriangles;
    for (TriangleType* t: TrianglesToSave) {
        t->UpdateNormal();
        TriangleType* NewTriangle = new TriangleType;

        // TODO: Keep track of materials
        NewTriangle->InterfaceID = t->InterfaceID;
        NewTriangle->ID = -t->ID; // Minus to distinguish from existing triangles

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


    FC_MESH FC;
    if (PerformTesting) {
        FC = this->CollapseEdgeTest(&TrianglesToSave, &TrianglesToRemove, &NewTriangles, EdgeToCollapse, RemoveVertexIndex );
    } else {
        FC = FC_OK;
    }

    if (FC != FC_OK) {
        for (TriangleType *t: NewTriangles) {
            delete t;
        }
        return FC;
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

#if SANITYCHECK == 1
    //this->DoSanityCheck();
#endif

    return FC_OK;

}

FC_MESH MeshManipulations :: CheckCoarsenNormal(std::vector<TriangleType*> *OldTriangles, std::vector<TriangleType*> *NewTriangles)
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

        if (newarea < 1e-8) {
            return FC_SMALLAREA; // TODO: Change to variable tol
        }

        LOG("angle1=%f,\tangle2=%f\tOldArea=%f\tNewArea=%f\n", angle1, angle2, oldarea, newarea);

        if (std::min(angle1, angle2)>MaxAngle) {
            MaxAngle = std::min(angle1, angle2);
        }

        if (angle1 > (170*2*3.1415/360)) {
            LOG("Should not be collapsed...\n", 0);
            return FC_NORMAL;
        }


    }

    if (MaxAngle > (15*2*3.1415/360))
        return FC_NORMAL;
    else
        return FC_OK;
}

FC_MESH MeshManipulations :: CheckCoarsenChord(EdgeType *EdgeToCollapse, VertexType* RemoveVertex, VertexType* SaveVertex)
{
    // If both vertices are located on an edge, proceed wih check. If not, collapsing is ok since this is not a chord.
    if (!(SaveVertex->IsPhaseEdgeVertex() && RemoveVertex->IsPhaseEdgeVertex())) {
        return FC_OK;
    }

    // If RemoveVertex is fixed, prevent collapsing
    if (RemoveVertex->IsFixedVertex()) {
        return FC_FIXEDVERTEX;
    }

    // If vertices are located on different PhaseEdges, collapsing is not ok

    // Here we know that RemoveVertex only contains one PhaseEdge and both vertices are located on the chord.
    PhaseEdge* RemoveVertexPhaseEdge = RemoveVertex->PhaseEdges.at(0);
    if ( std::find(SaveVertex->PhaseEdges.begin(), SaveVertex->PhaseEdges.end(), RemoveVertexPhaseEdge) == SaveVertex->PhaseEdges.end() ) {
        return FC_VERTICESONDIFFERENTSHAPES;
    }

    // Check if chord changes too much...
    std::array<VertexType*, 2> NewEdge;
    std::array<VertexType*, 2> OtherEdge;


    // Find the other edge connected to RemoveVertex (the one connected to the same PhaseEdge)
    std::vector<VertexType*> ConnectedVertices = RemoveVertexPhaseEdge->GiveVerticesConnectedToVertex(RemoveVertex);
    if (ConnectedVertices.size()==1) {
        // This should be ok if the edge is small. This simply removes the (very small) chord.
        // Note that this can imply duplicate triangles and it is quite cumbersome to solve this. Thus, we push this feature forward
        return FC_CHORD;
    }
    VertexType *OtherVertex = (ConnectedVertices[0] == SaveVertex) ? ConnectedVertices[1] : ConnectedVertices[0];
    NewEdge={SaveVertex, OtherVertex};
    OtherEdge = {RemoveVertex, OtherVertex};

    // Compute normalized vectors
    std::array<double, 3> NewNormal = ComputeNormalizedVector(NewEdge.at(0), NewEdge.at(1));
    std::array<std::array<double, 3>, 2> Normals = {ComputeNormalizedVector(OtherEdge.at(0), OtherEdge.at(1)),
                                                    ComputeNormalizedVector(EdgeToCollapse->Vertices.at(0), EdgeToCollapse->Vertices.at(1)) };
    double MaxAngle = 0.0;

    for (int i: {0, 1}) {
        double Alpha = std::acos(NewNormal.at(0)*Normals[i][0] + NewNormal.at(1)*Normals[i][1] + NewNormal.at(2)*Normals[i][2] );
        Alpha = std::min(std::abs(Alpha), std::abs(Alpha-3.141592));
        if (Alpha>MaxAngle) MaxAngle=Alpha;
    }

    if (MaxAngle > 45*3.141593/360) { // TODO: Use variable
        return FC_CHORD;
    }

    return FC_OK;

}

std::vector<VertexType *> MeshManipulations :: FindIndependentSet()
{
    std::vector<VertexType *> IndepSet;

    // Add all fixed vertices (Vertices shared between more than 2 phases)
    for (VertexType *v: this->Vertices) {
        if (v->IsFixedVertex()) {
            IndepSet.push_back(v);
        }
    }

    // Create vectors of vertices on edges and faces respectively
    std::vector<VertexType *> PhaseEdgeVertices;
    std::vector<VertexType *> FaceVertices;
    for (VertexType *v: this->Vertices) {
        if (v->IsPhaseEdgeVertex()) {
            PhaseEdgeVertices.push_back(v);
        } else {
            FaceVertices.push_back(v);
        }
    }


    for (std::vector<VertexType *> vl: {PhaseEdgeVertices, FaceVertices}) {
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

int MeshManipulations::FlipAll()
{
    int i=0;
    int flipcount = 0;
    for (EdgeType *e: this->Edges) {
        LOG ("Flip edge iteration %u: edge @%p (%u, %u)\n", i, e, e->Vertices[0]->ID, e->Vertices[1]->ID);
        if (this->FlipEdge(e)) flipcount++;
        std::ostringstream FileName;
        FileName << "/tmp/TestFlip_" << i << ".vtp";
        i++;
    }
    return flipcount;
}

bool MeshManipulations :: CoarsenMeshImproved()
{

    STATUS("Coarsen mesh\n", 0);

    bool CoarseningOccurs=true;
    int iter=0;
    int failcount;

    dooutputlogmesh(*this, "/tmp/Coarsening_%u.vtp", 0);

    while (CoarseningOccurs) {


        CoarseningOccurs = false;
        std::vector<VertexType*> IndepSet = FindIndependentSet();
        failcount=0;
        unsigned int i=0;

        std::vector<std::pair<double, EdgeType *>> EdgeLength;

        for (EdgeType *e: this->Edges) {
            EdgeLength.push_back(std::make_pair (e->GiveLength(), e));
        }

        std::sort(EdgeLength.begin(), EdgeLength.end());

        this->Edges.clear();
        for (std::pair<double, EdgeType *> epair: EdgeLength) {
            this->Edges.push_back(epair.second);
        }

        while (i<this->Edges.size()) {

            STATUS("%c[2K\rCoarsening iteration %u, failcount %u", 27, iter, failcount);
            fflush(stdout);

            EdgeType *e = this->Edges.at(i);
            if (e->GiveLength()<1e10) { // TODO: Use argument here
                bool BothInSet = (std::find(IndepSet.begin(), IndepSet.end(), e->Vertices[0]) == IndepSet.end()) && (std::find(IndepSet.begin(), IndepSet.end(), e->Vertices[1]) == IndepSet.end());
                if (!BothInSet) {
                    int D = 0;
                    if (std::find(IndepSet.begin(), IndepSet.end(), e->Vertices[D]) != IndepSet.end()) {
                        D = 1;
                    }

                    bool CoarseOk = (this->CollapseEdge(e, D) == FC_OK);

                    if (!CoarseOk) {
                        D = 1 ? 0 : 1;
                        CoarseOk = (this->CollapseEdge(e, D) == FC_OK);
                    }

                    if (CoarseOk) {
                        CoarseningOccurs = true;
                        dooutputlogmesh(*this, "/tmp/Coarsening_%u.vtp", iter+1);
#if SANITYCHECK == 1
//                        this->DoSanityCheck();
#endif

#if TEST_MESH_FOR_EACH_COARSENING_ITERATION
                        TetGenCaller Generator;
                        Generator.Mesh = this;
                        Generator.TestMesh();
#endif
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

    this->ExportSurface("/tmp/TestCoarsen_0.vtp", FT_VTK);
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
