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
private:
public:
    void SortAndFixBrokenEdge(std::vector<PhaseEdge*> *FixedEdges);
    void SplitAtVertex(VertexType *Vertex, std::vector<PhaseEdge*> *SplitEdges);
    void Smooth();

    std::vector<VertexType*> GetFlatListOfVertices();

    std::vector<VertexType*> FixedVertices;
    std::vector<std::vector<VertexType*>> EdgeSegments;
    std::vector<int> Phases;
    PhaseEdge();
};

}
#endif // PHASEEDGE_H
