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
#include "TetGenCaller.h"
#include "TriTriIntersect.h"

namespace voxel2tet
{

enum FC_MESH {FC_OK, FC_FIXEDVERTEX, FC_NORMAL, FC_CHORD, FC_SMALLAREA, FC_AREACHANGETOOLARGE, FC_TOOMANYTRIANGLES,
              FC_WORSEMINANGLE, FC_VERTICESONDIFFERENTSHAPES, FC_TRIANGLESINTERSECT, FC_DUPLICATETRIANGLE, FC_INVALIDEDGE};

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

    void SortEdgesByLength();
    void SortEdgesByMinArea();

public:
    /**
     * @brief TOL_MAXAREACHANGE determines the maximum change in area of triangles associated with the collapse of n edge
     */
    double TOL_MAXAREACHANGE;

    /**
     * @brief TOL_COL_SMALLESTAREA Smallest allowed area for triangles in edge collapsing
     */
    double TOL_COL_SMALLESTAREA;

    /**
     * @brief TOL_COL_MAXNORMALCHANGE Maximum allowed change in normal angle during collapse (rad)
     */
    double TOL_COL_MAXNORMALCHANGE;

    /**
     * @brief TOL_COL_CHORD_MAXNORMALCHANGE Maximum change in normal for elements connected to a chord during collapsing
     */
    double TOL_COL_CHORD_MAXNORMALCHANGE;

    /**
     * @brief TOL_SMALLESTAREA determines the smallest allowed area of a triangle in edge flipping
     */
    double TOL_FLIP_SMALLESTAREA;

    /**
     * @brief TOL_FLIP_MAXNORMALCHANGE Maximum change in angle for formal in flipping (rad)
     */
    double TOL_FLIP_MAXNORMALCHANGE;

    /**
     * @brief TOL_FLIP_MAXNORMALDIFFERENCE Maximum difference in normal angle of original triangles for allowing flipping (rad)
     */
    double TOL_FLIP_MAXNORMALDIFFERENCE;

    MeshManipulations(BoundingBoxType BoundingBox);

    void RemoveDegenerateTriangles();

    FC_MESH FlipEdge(EdgeType *Edge);

    /* Check if flipping is permitted. For return codes, see CollapseEdge
     *
     */
    FC_MESH CheckFlipNormal(std::vector<TriangleType*> *OldTriangles, std::array<TriangleType*, 2> NewTriangles);

    /* Collapses an edge. Returns FC_MESH value
     *
     */
    FC_MESH CollapseEdge(EdgeType *EdgeToCollapse, int RemoveVertexIndex, bool PerformTesting = true);
    FC_MESH CollapseEdgeTest(std::vector<TriangleType *> *TrianglesToSave, std::vector<TriangleType *> *TrianglesToRemove, std::vector<TriangleType *> *NewTriangles, EdgeType *EdgeToCollapse, int RemoveVertexIndex);
    FC_MESH CheckCoarsenNormal(std::vector<TriangleType*> *OldTriangles, std::vector<TriangleType*> *NewTriangles);
    FC_MESH CheckCoarsenChord(EdgeType *EdgeToCollapse, VertexType* RemoveVertex, VertexType* SaveVertex);

    bool CoarsenMesh();
    bool CoarsenMeshImproved();
    std::vector<VertexType *> FindIndependentSet();

    int FlipAll();

};

}
#endif // MESHMANIPULATIONS_H
