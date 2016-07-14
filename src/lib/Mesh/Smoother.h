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
public:
    /**
     * @brief constructor
     */
    Smoother();

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
    virtual void Smooth(std :: vector< VertexType * >Vertices, std :: vector< bool >Fixed,
                      std :: vector< std :: vector< VertexType * > >Connections,
                      MeshData *Mesh = NULL) = 0;

};

}

#endif // SMOOTHER_H
