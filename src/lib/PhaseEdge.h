#ifndef PHASEEDGE_H
#define PHASEEDGE_H
#include <vector>
#include <algorithm>

#include <MeshComponents.h>
#include <MiscFunctions.h>

namespace voxel2tet
{

class PhaseEdge
{
public:
    void SortAndFixBrokenEdge(std::vector<PhaseEdge*> *FixedEdges);
    std::vector<VertexType*> GetFlatListOfVertices();

    std::vector<std::vector<VertexType*>> EdgeSegments;
    std::vector<int> Phases;
    PhaseEdge();
};

}
#endif // PHASEEDGE_H
