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
        double Area=t->GiveArea();
        // LOG("Triangle %p (#%i) has area %f\n", t, i, Area);
        if (Area<1e-5) {
            // Note that this can occur by either a very small angle and/or two vertices located very close to each other
            LOG("Triangle %p has a too small area\n", t);

            // Find largest angle and the associated index
            int index;
            t->GiveLargestAngle(&index);

            // Find opposite edge and flip that edge
            std::array<EdgeType *, 3> Edges = t->GiveEdges();
            EdgeType *OppositeEdge = NULL;
            for (EdgeType *e: Edges) {
                if ( (e->Vertices[0]!=t->Vertices[index]) & (e->Vertices[1]!=t->Vertices[index])) {
                    OppositeEdge = e;
                    break;
                }
            }

            this->FlipEdge(OppositeEdge);

        }
    }
}

bool MeshManipulations :: FlipEdge(EdgeType *Edge)
{
    LOG ("\tFlip edge %p\n", Edge);

    std::vector<TriangleType *> EdgeTriangles = Edge->GiveTriangles();

    if (EdgeTriangles.size()!=2) {
        LOG("Unable to flip edge\n", 0);
        return false;
    }



    return true;


}

}
