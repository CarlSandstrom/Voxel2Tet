#include "SmoothingPipeline.h"

#include <algorithm>
#include <array>
#include <vector>

#include "Importer.h"
#include "MeshManipulations.h"
#include "MiscFunctions.h"
#include "PhaseEdge.h"
#include "Smoother.h"
#include "Surface.h"

namespace voxel2tet
{
SmoothingPipeline::SmoothingPipeline(Importer *Imp, MeshManipulations *Mesh, Smoother *EdgeSmoother,
                                     std::vector<Surface *> *Surfaces, std::vector<PhaseEdge *> *PhaseEdges,
                                     double eps)
{
    this->Imp = Imp;
    this->Mesh = Mesh;
    this->EdgeSmoother = EdgeSmoother;
    this->Surfaces = Surfaces;
    this->PhaseEdges = PhaseEdges;
    this->eps = eps;
}

void SmoothingPipeline::SmoothEdges() const
{
    //TODO: This does far too much stuff. A lot of thing has moved to the Smoother class. Clean this up.
    STATUS("Smooth edges (simultaneously)\n", 0);

    // We want to sort several lists according to one list. (http://stackoverflow.com/questions/1723066/c-stl-custom-sorting-one-vector-based-on-contents-of-another)
    struct VertexConnectivity
    {
        VertexType *v;
        std::vector<VertexType *> Connections;
        std::array<bool, 3> FixedDirections;
    };

    struct by_vertexptr
    {
        bool operator()(VertexConnectivity *a, VertexConnectivity *b)
        { return (a->v < b->v); }
    };

    std::vector<VertexConnectivity *> VertexConnections;

    // Collect all vertices on phase edges
    for (PhaseEdge *p : *this->PhaseEdges) {
        std::vector<std::vector<VertexType *> > Connections;
        std::vector<std::array<bool, 3> > FixedDirectionsList;

        std::vector<VertexType *> EdgeVertices = p->GetFlatListOfVertices();
        p->GiveTopologyLists(&Connections, &FixedDirectionsList);

        for (unsigned int i = 0; i < EdgeVertices.size(); i++) {
            VertexType *v = EdgeVertices.at(i);
            VertexConnectivity *vc = new VertexConnectivity;
            vc->v = v;
            vc->Connections.insert(vc->Connections.begin(), Connections.at(i).begin(), Connections.at(i).end());
            vc->FixedDirections = FixedDirectionsList.at(i);
            VertexConnections.push_back(vc);
        }
    }

    // Sort list by VertexType address and merge information of duplicate vertices
    std::sort(VertexConnections.begin(), VertexConnections.end(), by_vertexptr());

    unsigned int i = 0;
    while (i < VertexConnections.size()) {
        // Compare element i to element i+1. If the elements points to the same vertex, merge connections and remove element i+1
        if (i < (VertexConnections.size() - 1)) {
            if (VertexConnections.at(i)->v == VertexConnections.at(i + 1)->v) {
                VertexConnectivity *thisvc = VertexConnections.at(i);
                VertexConnectivity *nextvc = VertexConnections.at(i + 1);

                thisvc->Connections.insert(thisvc->Connections.end(), nextvc->Connections.begin(),
                                           nextvc->Connections.end());

                // Uniqueify connections
                std::sort(thisvc->Connections.begin(), thisvc->Connections.end());
                std::vector<VertexType *>::iterator it;
                it = std::unique(thisvc->Connections.begin(), thisvc->Connections.end());
                if (it != thisvc->Connections.end()) {
                    thisvc->Connections.erase(it);
                }

                // Erase next connection
                VertexConnections.erase(VertexConnections.begin() + i + 1);
            } else {
                i++;
            }
        } else {
            i++;
        }
    }

    // Build connectivity and FixedDirectionsLists
    std::vector<VertexType *> VertexList;
    std::vector<std::vector<VertexType *> > Connections;
    std::vector<bool> FixedDirectionsList;

    for (unsigned int i = 0; i < VertexConnections.size(); i++) {
        VertexType *v = VertexConnections.at(i)->v;

        VertexList.push_back(v);
        Connections.push_back(VertexConnections.at(i)->Connections);

        // Determine which directions are locked TODO: This should be done elsewhere
        for (int j = 0; j < 3; j++) {
            if ((v->get_c(j) > (this->Imp->GiveBoundingBox().maxvalues[j] - this->eps)) |
                    (v->get_c(j) < (this->Imp->GiveBoundingBox().minvalues[j] + this->eps))) {
                v->Fixed[j] = true;
            } else {
                v->Fixed[j] = false;
            }
        }
        FixedDirectionsList.push_back(false);
    }

    for (VertexConnectivity *v : VertexConnections) {
        delete v;
    }

    this->EdgeSmoother->Smooth(VertexList, this->Mesh);
}

void SmoothingPipeline::SmoothSurfaces() const
{
    STATUS("Smooth surfaces\n", 0);

    int i = 1;
    for (Surface *s : *this->Surfaces) {
        STATUS("Smoothing surface %i (%i)\n", i, this->Surfaces->size());
        s->Smooth(this->Mesh);
        i++;
    }
}
}
