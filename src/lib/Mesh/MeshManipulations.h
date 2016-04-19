#ifndef MESHMANIPULATIONS_H
#define MESHMANIPULATIONS_H

#include <vector>
#include <array>
#include <algorithm>
#include <armadillo>
#include <math.h>
#include <utility>
#include <map>

#include "MeshData.h"
#include "MiscFunctions.h"
#include "MeshComponents.h"

namespace voxel2tet
{

enum FC_MESH {FC_OK, FC_FIXEDVERTEX, FC_NORMAL, FC_CHORD, FC_SMALLAREA, FC_TOOMANYTRIANGLES, FC_WORSEMINANGLE, FC_VERTICESONDIFFERENTSHAPES};

class MeshManipulations: public MeshData
{
private:
    /**
     * @brief GetFlippedEdgeData Produce new edge and new triangles for a flipped edge
     * @param EdgeToFlip Vector to flip
     * @param NewEdge New edge (out)
     * @param NewTriangles NewTriangles (out)
     * @return  Status
     */
    FC_MESH GetFlippedEdgeData(EdgeType *EdgeToFlip, EdgeType *NewEdge, std::array<TriangleType*, 2> *NewTriangles);

public:
    MeshManipulations(BoundingBoxType BoundingBox);

    void RemoveDegenerateTriangles();

    FC_MESH FlipEdge(EdgeType *Edge);

    /* Check if flipping is permitted. For return codes, see CollapseEdge
     *
     */
    FC_MESH CheckFlipNormal(std::vector<TriangleType*> *OldTriangles, std::array<TriangleType*, 2> NewTriangles);

    /* Collapses an edge. Return values
     *
     *  0   Edge Collapsed
     *  -1  Failes due to fixed vertex
     *  -2  Failed due to a too large change in normal
     *  -3  Failed due to a too large change in chord
     *  -4  Failed due to a too small area of new triangle(s)
     *
     */
    FC_MESH CollapseEdge(EdgeType *EdgeToCollapse, int RemoveVertexIndex, bool PerformTesting = true);
    FC_MESH CollapseEdgeTest(std::vector<TriangleType *> *TrianglesToSave, std::vector<TriangleType *> *NewTriangles, EdgeType *EdgeToCollapse, int RemoveVertexIndex);
    FC_MESH CheckCoarsenNormal(std::vector<TriangleType*> *OldTriangles, std::vector<TriangleType*> *NewTriangles);
    FC_MESH CheckCoarsenChord(EdgeType *EdgeToCollapse, VertexType* RemoveVertex, VertexType* SaveVertex);

    bool CoarsenMesh();
    bool CoarsenMeshImproved();
    std::vector<VertexType *> FindIndependentSet();

    int FlipAll();

};

}
#endif // MESHMANIPULATIONS_H
