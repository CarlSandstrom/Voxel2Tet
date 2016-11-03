#ifndef SMOOTHER_H
#define SMOOTHER_H

#include <vector>

#include "MeshComponents.h"
#include "MeshData.h"
#include "MeshManipulations.h"

namespace voxel2tet
{

/**
 * @brief Abstract class for all smoothers
 */
class Smoother
{
    friend std::ostream &operator<<(std::ostream &stream, const Smoother &ThisSmoother)
    {
        stream << ThisSmoother.DoOutput();
        return stream;
    }

public:
    /**
     * @brief constructor
     */
    Smoother();

    /**
     * @brief Creates a vector where each element is a vector of vertices connected to the corresponding node in the Vertices vector.
     * Static vertices has no connections and vertices on a PhaseEdge does not have any connection with vertices on surfaces. However,
     * vertices on the surfaces has connections to vertices on PhaseEdges.
     * @param Vertices
     * @return
     */
    std::vector<std::vector<VertexType *>> GetConnectivityVector(std::vector<VertexType *> Vertices);

    /**
     * @brief Smooth all given vertices.
     *
     * The smoothing depends on the connectivity and in wich directions the vertices are fixed
     *
     * @param Vertices Vector of pointer to objects of VertexType
     * @param Mesh For exporting (debugging purposes)
     */
    virtual void Smooth(std::vector<VertexType *> Vertices, MeshData *Mesh = NULL) = 0;

    /**
     * Checks for penetrating triangles
     * @param Vertices All triangles connected to any vertex in Vertices will be checked for penetration
     * @param Mesh Mesh object containing the triangles
     * @return List of intersecting triangle pairs
     */
    std::vector<std::pair<TriangleType *, TriangleType *> >
    CheckPenetration(std::vector<VertexType *> *Vertices, MeshManipulations *Mesh);

    /**
     * Performs pullback of vertices connected to intersecting triangles
     * @param Vertices List of vertices. All triangles connected to any vertex in Vertices will be checked for penetration
     * @param Mesh Mesh object
     */
    void PullBackAtIntersections(std::vector<VertexType *> Vertices, MeshManipulations *Mesh);

protected:
    /**
     * Returns a string for output in .stat file
     * @return
     */
    virtual std::string DoOutput() const = 0;
};

}

#endif // SMOOTHER_H
