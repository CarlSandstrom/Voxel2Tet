#ifndef SMOOTHER_H
#define SMOOTHER_H

#include <vector>

#include "MeshComponents.h"
#include "MeshData.h"

namespace voxel2tet
{

/**
 * @brief Abstract class for all smoothers
 */
class Smoother
{
    friend std::ostream &operator<<(std::ostream &stream, const Smoother &Smoother) {return stream; }
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
     * The smoothing depends on the connectivity and in wich directions the verteces are fixed
     *
     * @param Vertices Vector of pointer to objects of VertexType
     * @param Fixed Vextor of array of 3 indicating in which directions the vertices are fixed. Should be the same length as Vertices
     * @param Connections Gives the connections of each vertex. Should be the same length as Vertices
     * @param Mesh For exporting (debugging purposes)
     */
    virtual void Smooth(std :: vector< VertexType * >Vertices, MeshData *Mesh = NULL) = 0;

    std :: vector< std :: pair< TriangleType *, TriangleType * > >CheckPenetration(std :: vector< VertexType * > *Vertices, MeshData *Mesh);

    void PullBackAtIntersections(std :: vector< VertexType * > Vertices, MeshData *Mesh);

};

}

#endif // SMOOTHER_H
