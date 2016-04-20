#ifndef PHASEEDGE_H
#define PHASEEDGE_H
#include <vector>
#include <array>
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
public:
    Options *Opt; // TODO: Should be private

    void SortAndFixBrokenEdge(std::vector<PhaseEdge*> *FixedEdges);
    bool SplitAtVertex(VertexType *Vertex, std::vector<PhaseEdge*> *SplitEdges);
    void GiveTopologyLists(std::vector<std::vector<VertexType *>> *Connections, std::vector<std::array<bool,3>> *FixedDirections);
    std::vector<VertexType *> GiveVerticesConnectedToVertex(VertexType *v);
    void Smooth(MeshData *Mesh);

    bool IsClosed();

    std::vector<VertexType*> GetFlatListOfVertices();

    std::vector<VertexType*> FixedVertices;
    std::vector<std::array<VertexType*, 2>> EdgeSegments; // TODO:Should be vector of arrays
    std::vector<int> Phases;
    PhaseEdge(Options *Opt);
};

}
#endif // PHASEEDGE_H
