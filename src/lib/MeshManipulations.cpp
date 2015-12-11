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

            if (!this->FlipEdge(OppositeEdge)) {
                STATUS("The longest edge of triangle %u failed to flipped\n", i);
            }
        }
    }
}

bool MeshManipulations :: FlipEdge(EdgeType *Edge)
{
    LOG ("\tFlip edge %p\n", Edge);

    std::vector<TriangleType *> EdgeTriangles = Edge->GiveTriangles();

    if (EdgeTriangles.size()!=2) {
        LOG("Unable to flip edge. To many or only one triangle connected\n", 0);
        return false;
    }

    std::vector<VertexType *> NewEdge;

    for (int i=0; i<2; i++) std::sort(EdgeTriangles.at(i)->Vertices.begin(), EdgeTriangles.at(i)->Vertices.end());

    std::set_symmetric_difference(EdgeTriangles.at(0)->Vertices.begin(), EdgeTriangles.at(0)->Vertices.end(),
                                  EdgeTriangles.at(1)->Vertices.begin(), EdgeTriangles.at(1)->Vertices.end(), std::back_inserter(NewEdge));

    // Update triangles
    for (int i=0; i<2; i++) {
        TriangleType *t = EdgeTriangles.at(i);
        for (int j=0; j<3; j++) {
            t->Vertices[j]->RemoveTriangle(t);
        }
    }

    for (int i=0; i<2; i++) {
        TriangleType *t = EdgeTriangles.at(i);
        t->Vertices[0] = NewEdge[0];
        t->Vertices[1] = NewEdge[1];
        t->Vertices[2] = Edge->Vertices[i];
        for (int j=0; j<3; j++) {
            t->Vertices[j]->AddTriangle(t);
        }
    }

    // Update edge
    for (int i: {0, 1}) Edge->Vertices[i]->RemoveEdge(Edge);

    for (int i: {0, 1}) {
        Edge->Vertices[i] = NewEdge[i];
        NewEdge[i]->AddEdge(Edge);
    }

    return true;

}

}
