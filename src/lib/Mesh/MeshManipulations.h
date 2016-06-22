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

    /**
     * @brief TOL_COL_MAXVOLUMECHANGE Maximum change in volume for one edge collapse
     */
    double TOL_COL_MAXVOLUMECHANGE;

    /**
     * @brief TOL_COL_MAXERROR Largest accumulated error in nodes for collapsing
     */
    double TOL_COL_MAXERROR;

    MeshManipulations(BoundingBoxType BoundingBox);

    void RemoveDegenerateTriangles();

    /**
     * @brief Flips edge if permitted
     * @param Edge [in] Edge to be flipped
     * @return Value depends on success of flip
     */
    FC_MESH FlipEdge(EdgeType *Edge);

    /**
     * @brief Check if flipping is permitted. For return codes, see CollapseEdge
     * @return Returns wether flipping is permitted
     */
    FC_MESH CheckFlipNormal(std::vector<TriangleType*> *OldTriangles, std::array<TriangleType*, 2> NewTriangles);

    /**
     * @brief Collapses an edge if possible. Returns FC_MESH value
     * @return Returns wether collapsing was successfull or not
     */
    FC_MESH CollapseEdge(EdgeType *EdgeToCollapse, int RemoveVertexIndex, bool PerformTesting = true);

    /**
     * @brief Determines wether collapsing an edge is permitted
     * @param TrianglesToSave [in] Triangles to keep
     * @param TrianglesToRemove [in] Triangles to remove from existing triangulation
     * @param NewTriangles [in] New triangles
     * @param EdgeToCollapse [in] Edge to collapse
     * @param RemoveVertexIndex [in] Index in EdgeToCollapse to the vertex about to be collapsed
     * @return
     */
    FC_MESH CollapseEdgeTest(std::vector<TriangleType *> *TrianglesToSave, std::vector<TriangleType *> *TrianglesToRemove, std::vector<TriangleType *> *NewTriangles, EdgeType *EdgeToCollapse, int RemoveVertexIndex);

    /**
     * @brief Check if the change in normals of existing triangles are small enough to allow edge collapse
     */
    FC_MESH CheckCoarsenNormal(std::vector<TriangleType*> *OldTriangles, std::vector<TriangleType*> *NewTriangles);

    /**
     * @brief Check if edge collapse results in a too large loss of volume.
     * @param OldTriangles
     * @param TrianglesToRemove
     * @param NewTriangles
     * @param error [out] Accumulated error
     * @return
     */
    FC_MESH CheckCoarsenNormalImproved(std::vector<TriangleType*> *OldTriangles, std::vector<TriangleType*> *TrianglesToRemove, std::vector<TriangleType*> *NewTriangles, double &error);

    FC_MESH CheckCoarsenChord(EdgeType *EdgeToCollapse, VertexType* RemoveVertex, VertexType* SaveVertex);

    bool CoarsenMesh();
    bool CoarsenMeshImproved();
    std::vector<VertexType *> FindIndependentSet();

    int FlipAll();

};

}
#endif // MESHMANIPULATIONS_H
