#ifndef SMOOTHER_H
#define SMOOTHER_H

#include <vector>

#include "MeshComponents.h"
#include "MeshData.h"

namespace voxel2tet
{

class Smoother
{
public:
    Smoother();

    virtual void Smooth(std :: vector< VertexType * >Vertices, std :: vector< bool >Fixed,
                      std :: vector< std :: vector< VertexType * > >Connections,
                      MeshData *Mesh = NULL) = 0;

};

}

#endif // SMOOTHER_H
