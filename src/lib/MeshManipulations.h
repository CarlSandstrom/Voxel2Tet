#ifndef MESHMANIPULATIONS_H
#define MESHMANIPULATIONS_H

#include <vector>
#include <array>
#include <algorithm>
#include <armadillo>
#include <math.h>

#include "MeshData.h"
#include "MiscFunctions.h"
#include "MeshComponents.h"

namespace voxel2tet
{

class MeshManipulations: public MeshData
{
private:
    bool GetFlippedEdgeData(EdgeType *EdgeToFlip, EdgeType *NewEdge, std::array<TriangleType*, 2> *NewTriangles);



public:
    MeshManipulations(BoundingBoxType BoundingBox);

    void RemoveDegenerateTriangles();

    bool FlipEdge(EdgeType *Edge);
    bool CheckFlipNormal(std::vector<TriangleType*> *OldTriangles, std::array<TriangleType*, 2> NewTriangles);

    bool CollapseEdge(EdgeType *EdgeToCollapse, int RemoveVertexIndex);
    bool CheckCoarsenNormal(std::vector<TriangleType*> *OldTriangles, std::vector<TriangleType*> *NewTriangles);
    bool CheckCoarsenChord(EdgeType *EdgeToCollapse, VertexType* RemoveVertex, VertexType* SaveVertex);

    bool CoarsenMesh();

};

}
#endif // MESHMANIPULATIONS_H
