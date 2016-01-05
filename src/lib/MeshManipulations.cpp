#include "MeshManipulations.h"

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
                STATUS("The longest edge of triangle %u failed to flipped\n", i);
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

    // Remove edge from vertices
    for (int i: {0, 1}) Edge->Vertices[i]->RemoveEdge(Edge);

    // Add new edge to vertices
    for (int i: {0, 1}) {
        Edge->Vertices[i] = NewEdge.Vertices[i];
        Edge->Vertices[i]->AddEdge(Edge);
    }

    // 5. Free old triangles
    for (TriangleType *t: NewTriangles) {
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
    // Cannot remove a fixed vertex
    if (EdgeToCollapse->Vertices[RemoveVertexIndex]->FixedVertex) return false;

    int SaveVertexIndex = (RemoveVertexIndex == 0) ? 1 : 0;

    // TODO: In the Python code, is the topology check necessary?

    // Create new triangles. These are create by moving RemoveVertex to the other end of the edge and remove the 0-area triangles
    std::vector<TriangleType*> TrianglesToRemove = EdgeToCollapse->GiveTriangles();
    std::vector<TriangleType*> ConnectedTriangle = EdgeToCollapse->Vertices[RemoveVertexIndex]->Triangles;

    std::sort(TrianglesToRemove.begin(), TrianglesToRemove.end());
    std::sort(ConnectedTriangle.begin(), ConnectedTriangle.end());

    std::vector<TriangleType*> TrianglesToSave;
    std::set_difference(ConnectedTriangle.begin(), ConnectedTriangle.end(),
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

    // Check if NewTriangles is a good approximation




    for (TriangleType* t: TrianglesToRemove) {
        this->Triangles.erase(std::remove(this->Triangles.begin(), this->Triangles.end(), t), this->Triangles.end());
    }
    for (TriangleType* t: NewTriangles) {
        this->Triangles.push_back(t);
    }

    TrianglesToSave.clear();
}

}
