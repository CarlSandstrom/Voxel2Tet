#ifndef PHASEEDGE_H
#define PHASEEDGE_H
#include <vector>
#include <array>
#include <algorithm>
#include <MeshComponents.h>
#include <MiscFunctions.h>

#include "SpringSmoother.h"
#include "Options.h"
#include "MeshData.h"

namespace voxel2tet
{
/**
 * @brief The PhaseEdge class contains information and means of manipulation of the curves formd
 * where three or more materials meet.
 *
 * A phase edge is the curve formed where three or more materials meet. When the voxel data is
 * analyzed, PhaseEdge objects are created for each set of combinations of materials. For example,
 * one PhaseEdge object for the curves where materials 1, 2 and 3 meet and one for the curve where
 * materials 2, 3 and 4 meet. Note that at this stage, the curve in the PhaseEdge objects can very
 * well be segments spread over the volume.
 */
class PhaseEdge
{
private:
    /**
     * @brief Opt hold the Options object created from the command line arguments
     */
    Options *Opt;

    /**
     * @brief FixedVertices holds all vertices that are fixed.
     *
     * Fixed vertices are vertices where a PhaseEdge ends, i.e. where it meets one or more new
     * PhaseEdge objects.
     */
    std :: vector< VertexType * >FixedVertices;

public:

    /**
     * @brief EdgeSegments hold all edge segments for this edge
     */
    std :: vector< std :: array< VertexType *, 2 > >EdgeSegments;

    /**
     * @brief Pointer to smoother object for this phase edge.
     */
    SpringSmoother *EdgeSmoother;

    /**
     * @brief Constructor
     * @param Opt Input. An Options object for communicating the command line parameters.
     * @param EdgeSmoother Input. Smoother object to the phase edge.
     */
    PhaseEdge(Options *Opt, SpringSmoother *EdgeSmoother);

    /**
     * @brief SortAndFixBrokenEdge Identifies all separate phase edges within this PhaseEdge object and outputs all (new) internally connected phase edges.
     * @param FixedEdges Output. Holds all new, internally connected, phase edges.
     */
    void SortAndFixBrokenEdge(std :: vector< PhaseEdge * > *FixedEdges);

    /**
     * @brief Creates a list of connected vertices for each vertex in the object as well as an array describing
     * in which directions a vertex is allowed to move.
     *
     * @param Connections Output. List of connected vertices for each vertex in PhaseEdge object. Each item corresponds
     *                    to the vertex of the same index.
     * @param FixedDirections Output. List of fixed directions for each vertex. Each item corresponds to the vertex with
     *                    the same index.
     */
    void GiveTopologyLists(std :: vector< std :: vector< VertexType * > > *Connections, std :: vector< std :: array< bool, 3 > > *FixedDirections);

    /**
     * @brief Returns the vertices connected to vertex \p v.
     * @param v Input. Vertex to which connected vertices are sought.
     * @return Returns vertices connected to \p v.
     */
    std :: vector< VertexType * >GiveVerticesConnectedToVertex(VertexType *v);

    /**
     * @brief Perform smoothing on this PhaseEdge.
     * @param Mesh
     */
    void Smooth(MeshData *Mesh);

    /**
     * @brief Add an edge segment to the PhaseEdge
     * @param v1 Input. Vertex
     * @param v2 Input. Vertex
     */
    void AddPhaseEdgeSegment(VertexType *v1, VertexType *v2);

    /**
     * @brief Determines if the PhaseEdge is closed or not, i.e. if it has ends
     * @return Wither the curve is closed or not.
     */
    bool IsClosed();

    /**
     * @brief Returns a list of all vertices in the order they appear on the PhaseEdge.
     * @return List of vertices
     */
    std :: vector< VertexType * >GetFlatListOfVertices();

    /**
     * @brief List of all phase IDs ni contact with this phase edge
     */
    std :: vector< int >Phases;

    void LogPhaseEdge();
};
}
#endif // PHASEEDGE_H
