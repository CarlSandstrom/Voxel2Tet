#ifndef MESHMANIPULATIONS_H
#define MESHMANIPULATIONS_H

#include <vector>
#include <array>
#include <algorithm>

#include "MeshData.h"
#include "MiscFunctions.h"
#include "MeshComponents.h"

namespace voxel2tet
{

class MeshManipulations: public MeshData
{
public:
    MeshManipulations(BoundingBoxType BoundingBox);

    void RemoveDegenerateTriangles();

    bool FlipEdge(EdgeType *Edge);
};

}
#endif // MESHMANIPULATIONS_H
