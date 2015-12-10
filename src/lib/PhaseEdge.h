#ifndef PHASEEDGE_H
#define PHASEEDGE_H
#include <vector>
#include <algorithm>

#include <MeshComponents.h>
#include <MiscFunctions.h>
#include "Options.h"
#include "MeshData.h"

namespace voxel2tet
{

class PhaseEdge
{
private:
    Options *Opt;
public:
    void SortAndFixBrokenEdge(std::vector<PhaseEdge*> *FixedEdges);
    void SplitAtVertex(VertexType *Vertex, std::vector<PhaseEdge*> *SplitEdges);
    void GiveTopologyLists(std::vector<std::vector<VertexType *>> *Connections, std::vector<std::array<bool,3>> *FixedDirections);
    void Smooth(MeshData *Mesh);

    bool IsClosed();

    std::vector<VertexType*> GetFlatListOfVertices();

    std::vector<VertexType*> FixedVertices;
    std::vector<std::vector<VertexType*>> EdgeSegments; // TODO:Should be vector of arrays
    std::vector<int> Phases;
    PhaseEdge(Options *Opt);
};

}
#endif // PHASEEDGE_H
