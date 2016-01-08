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
        t->Vertices[0] = NewEdgeVertices[0];
        t->Vertices[1] = NewEdgeVertices[1];
        t->Vertices[2] = EdgeToFlip->Vertices[i];
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
    LOG ("\tFlip edge %p\n", Edge);

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

    // Update mesh data

    // Update triangles
    for (int i=0; i<2; i++) {
        for (int j=0; j<3; j++) {
            EdgeTriangles.at(i)->Vertices[j] = NewTriangles.at(i)->Vertices[j];
        }
    }

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
/*    for (TriangleType *t: NewTriangles) {
        for (VertexType *v: t->Vertices) {
            v->AddTriangle(t);
        }
    }*/

    for (TriangleType *t: NewTriangles) {
        this->AddTriangle(t);
    }

    // 5. Free old triangles
    for (TriangleType *t: EdgeTriangles) {
        this->RemoveTriangle(t);
        delete t;
    }

    return true;

}

bool MeshManipulations :: CheckFlipNormal(std::vector<TriangleType*> *OldTriangles, std::array<TriangleType *, 2> NewTriangles)
{
    double MaxAngle = 0.0;

    for (TriangleType *OldTriangle: *OldTriangles) {
        for (unsigned int j=0; j<NewTriangles.size(); j++) {
            TriangleType *NewTriangle = NewTriangles.at(j);
            std::array<double, 3> OldNormal = OldTriangle->GiveUnitNormal();
            std::array<double, 3> NewNormal = NewTriangle->GiveUnitNormal();
            double angle1 = std::acos( OldNormal[0]*NewNormal[0] + OldNormal[1]*NewNormal[1] + OldNormal[2]*NewNormal[2]);
            double angle2 = std::acos( -(OldNormal[0]*NewNormal[0] + OldNormal[1]*NewNormal[1] + OldNormal[2]*NewNormal[2]) );
            if (std::min(angle1, angle2)>MaxAngle) {
                MaxAngle = std::min(angle1, angle2);
            }
        }
    }

    if (MaxAngle > (90*3.1415/360))
        return false;
    else
        return true;

}

bool MeshManipulations :: CollapseEdge(EdgeType *EdgeToCollapse, int RemoveVertexIndex)
{
    // TODO: In its current setting, this procedure only checks if collapsing is ok from the current topology. It should compare to the original topology, otherwise degeneration can occur gradually

    // Cannot remove a fixed vertex
    LOG("Collapse edge %p (%u, %u) by removing vertex %u\n", EdgeToCollapse, EdgeToCollapse->Vertices[0]->ID,
            EdgeToCollapse->Vertices[1]->ID, RemoveVertexIndex);
    if (EdgeToCollapse->Vertices[RemoveVertexIndex]->IsFixedVertex()) return false;

    int SaveVertexIndex = (RemoveVertexIndex == 0) ? 1 : 0;

    // Create new triangles. These are create by moving RemoveVertex to the other end of the edge and remove the 0-area triangles
    std::vector<TriangleType*> TrianglesToRemove = EdgeToCollapse->GiveTriangles();
    std::vector<TriangleType*> ConnectedTriangles = EdgeToCollapse->Vertices[RemoveVertexIndex]->Triangles;

    std::sort(TrianglesToRemove.begin(), TrianglesToRemove.end());
    std::sort(ConnectedTriangles.begin(), ConnectedTriangles.end());

    std::vector<TriangleType*> TrianglesToSave;
    std::set_difference(ConnectedTriangles.begin(), ConnectedTriangles.end(),
                        TrianglesToRemove.begin(), TrianglesToRemove.end(), std::inserter( TrianglesToSave, TrianglesToSave.begin() ) );

    std::vector<TriangleType*> NewTriangles;
    for (TriangleType* t: TrianglesToSave) {
        TriangleType* NewTriangle = new TriangleType;

        // TODO: Keep track of materials
        NewTriangle->InterfaceID = t->InterfaceID;
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

    // Check if NewTriangles are good replacements
    if (!this->CheckCoarsenNormal(&TrianglesToSave, &NewTriangles)) {
        return false;
    }

    if (!this->CheckCoarsenChord(EdgeToCollapse, EdgeToCollapse->Vertices[RemoveVertexIndex], EdgeToCollapse->Vertices[SaveVertexIndex])) {
        return false;
    }

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
                if (e->Vertices[i]==EdgeToCollapse->Vertices.at(RemoveVertexIndex)) {
                    e->Vertices[i] = EdgeToCollapse->Vertices.at(SaveVertexIndex);
                }
            }
        }
    }

    // Remove triangles
    for (TriangleType* t: ConnectedTriangles) {
        for (VertexType *v: t->Vertices) {
            v->RemoveTriangle(t);
        }
        this->Triangles.erase(std::remove(this->Triangles.begin(), this->Triangles.end(), t), this->Triangles.end());
        delete t;
    }

    // Remove edges
    for (EdgeType* e: EdgesToRemove) {
        this->Edges.erase((std::remove(this->Edges.begin(), this->Edges.end(), e)));
        delete e;
    }


    // Add new triangles
    for (TriangleType* t: NewTriangles) {
        this->Triangles.push_back(t);
        for (VertexType *v: t->Vertices) {
            v->AddTriangle(t);
        }
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
        if (std::min(angle1, angle2)>MaxAngle) {
            MaxAngle = std::min(angle1, angle2);
        }

    }

    if (MaxAngle > (20*3.1415/360))
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

bool MeshManipulations :: CoarsenMesh()
{

    this->ExportVTK("/tmp/TestCoarsen_0.vtp");
    int iter=1;
    bool EdgeCollapsed = true;
    while (EdgeCollapsed) {
        EdgeCollapsed=false;
        int failcount=0;
        for (EdgeType* e: this->Edges) {

            LOG("Collapse iteration %u, failcount=%u\n", iter, failcount);

            if (!this->CollapseEdge(e,0)) {
                EdgeCollapsed=this->CollapseEdge(e,1);
            } else {
                EdgeCollapsed=true;
            }

            if (EdgeCollapsed) {
                LOG("Collapse iteration %u succeeded after %u failed tries\n", iter, failcount);
                std::ostringstream FileName;
                FileName << "/tmp/TestCoarsen_" << iter << ".vtp";
                this->ExportVTK( FileName.str() );
                iter++;
                break;
            } else {
                failcount++;
            }
        }
    }
    return true;
}

}
