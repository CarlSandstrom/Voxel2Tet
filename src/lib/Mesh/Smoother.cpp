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


            // Finally, ensure that the referenced vertices exists in the list
/*            size_t i=0;
            while (i<ConnectedVertices.size()) {
                if (std::find(Vertices.begin(), Vertices.end(), ConnectedVertices[i]) == Vertices.end()) {
                    // Connected vertex does not exists in list. Remove.
                    ConnectedVertices.erase(std::remove(ConnectedVertices.begin(), ConnectedVertices.end(), ConnectedVertices[i]), ConnectedVertices.end());
                } else {
                    i++;
                }
            }*/
        }
        Connectivity.push_back(ConnectedVertices);
    }
    return Connectivity;
}

Smoother::Smoother()
{

}

}
