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

/**
 * @brief The MeshManipulations class supplies methods for manipulating the mesh.
 */
class MeshManipulations : public MeshData
{
private:

    double LongestEdgeLength;
    void UpdateLongestEdgeLength();

    /**
     * @brief GetFlippedEdgeData Produce new edge and new triangles for a flipped edge
     * @param EdgeToFlip Vector to flip
     * @param NewEdge New edge (out)
     * @param NewTriangles NewTriangles (out)
     * @return  Status
     */
    FC_MESH GetFlippedEdgeData(EdgeType *EdgeToFlip, EdgeType *NewEdge, std :: array< TriangleType *, 2 > *NewTriangles);

    void SortEdgesByLength();
    void SortEdgesByMinArea();

public:

    /**
     * @brief Determines the maximum change in area of triangles associated with a flip of an edge
     */
    double TOL_FLIP_MAXAREACHANGE;

    /**
     * @brief Smallest allowed area for triangles in edge collapsing
     */
    double TOL_COL_SMALLESTAREA;

    /**
     * @brief Maximum allowed change in normal angle during collapse (rad)
     */
    double TOL_COL_MAXNORMALCHANGE;

    /**
     * @brief Maximum change in normal for elements connected to a chord during collapsing
     */
    double TOL_COL_CHORD_MAXNORMALCHANGE;

    /**
     * @brief Determines the smallest allowed area of a triangle in edge flipping
     */
    double TOL_FLIP_SMALLESTAREA;

    /**
     * @brief Maximum change in angle for formal in flipping (rad)
     */
    double TOL_FLIP_MAXNORMALCHANGE;

    /**
     * @brief Maximum difference in normal angle of original triangles for allowing flipping (rad)
     */
    double TOL_FLIP_MAXNORMALDIFFERENCE;

    /**
     * @brief Maximum change in volume for one edge collapse
     */
    double TOL_COL_MAXVOLUMECHANGE;

    /**
     * @brief Largest accumulated error in nodes for collapsing
     */
    double TOL_COL_MAXERROR_ACCUMULATED;

    double TOL_COL_MINANGLE;

    /**
     * @brief Constructor
     * @param BoundingBox Bounding box of mesh
     */
    MeshManipulations(BoundingBoxType BoundingBox);

    /**
     * @brief Flips edge if permitted
     * @param Edge [in] Edge to be flipped
     * @return Value depends on success of flip
     */
    FC_MESH FlipEdge(EdgeType *Edge, bool SkipIntersectionCheck = false);

    /**
     * @brief Check if flipping is permitted. For return codes, see CollapseEdge
     * @return Returns wether flipping is permitted
     */
    FC_MESH CheckFlipNormal(std :: vector< TriangleType * > *OldTriangles, std :: array< TriangleType *, 2 >NewTriangles);

    /**
     * @brief Collapses an edge if possible. Returns FC_MESH value
     * @param EdgeToCollapse Pointer to the edge to collapse
     * @param RemoveVertexIndex Index of the vertex to remove. This is the index wrt the edge, i.e. the index is either 0 or 1.
     * @param PerformTesting
     * @return Returns wether collapsing was successfull or not
     */
    FC_MESH CollapseEdge(EdgeType *EdgeToCollapse, int RemoveVertexIndex, bool PerformTesting = true);

    /**
     * @brief Determines wether collapsing an edge is permitted
     * @param TrianglesToSave [in] Triangles containing the vertex to be removed, but will still have an area after collapsing
     * @param TrianglesToRemove [in] Triangles to remove from existing triangulation
     * @param NewTriangles [in] New triangles. Basically same triangles as in TrianglesToSave but where the vertex to be removed has been exchanged to the vertex to keep.
     * @param EdgeToCollapse [in] Edge to collapse
     * @param RemoveVertexIndex [in] Index in EdgeToCollapse to the vertex about to be collapsed
     * @return FC_MESH type stating if the collapse was a success
     */
    FC_MESH CollapseEdgeTest(std :: vector< TriangleType * > *TrianglesToSave, std :: vector< TriangleType * > *TrianglesToRemove, std :: vector< TriangleType * > *NewTriangles, EdgeType *EdgeToCollapse, int RemoveVertexIndex);

    /**
     * @brief Check if the change in normals of existing triangles are small enough to allow edge collapse
     *
     * The method compares the normals of triangles in OldTriangles and NewTriangles and computes the angle between then.
     *
     * @param OldTriangles Vector of pointers to TriangleTypes
     * @param NewTriangles Vector of pointers to TriangleTypes
     * @return FC_MESH type stating if collapsing is ok or not
     */
    FC_MESH CheckCoarsenNormal(std :: vector< TriangleType * > *OldTriangles, std :: vector< TriangleType * > *NewTriangles);

    /**
     * @brief Check if edge collapse results in a too large loss of volume or a too large change in normal.
     *
     * The method compares the normals of the triangle in OldTriangles and NewTriangles and computes the angle between them. If
     * the angle is too large, collapsing is prohibited.
     *
     * It also computes the change in volume the collapse would imply and stops the collapse if the error or the accumulated error is to large.
     *
     * @param OldTriangles
     * @param TrianglesToRemove
     * @param NewTriangles
     * @param error [out] Accumulated error
     * @return FC_MESH type stating if collapsing is ok or not
     */
    FC_MESH CheckCoarsenNormalImproved(std :: vector< TriangleType * > *OldTriangles, std :: vector< TriangleType * > *TrianglesToRemove, std :: vector< TriangleType * > *NewTriangles, double &error);

    /**
     * @brief Check if change in chord is small enough to allow collapsing
     * @param EdgeToCollapse
     * @param RemoveVertex
     * @param SaveVertex
     * @return FC_MESH type stating if collapsing is ok or not
     */
    FC_MESH CheckCoarsenChord(EdgeType *EdgeToCollapse, VertexType *RemoveVertex, VertexType *SaveVertex);

    /**
     * @brief Coarsen surface mesh
     */
    void CoarsenMesh();

    /**
     * @brief Find the set of independet vertices.
     * @return Vector of pointers to VertexType objects
     */
    std :: vector< VertexType * >FindIndependentSet();

    /**
     * @brief Perform flipping of edges until no more edges can be flipped. This implies better quality of the mesh.
     * @return Number of flips performed
     */
    int FlipAll(bool SkipIntersectionCheck = false);

    int CleanupTetrahedrals();

    bool RemoveTetrahedron(int index);

    bool RemoveTetrahedron(TetType *t);
};
}
#endif // MESHMANIPULATIONS_H
