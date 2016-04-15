#ifndef MESHGENERATOR3D_H
#define MESHGENERATOR3D_H

#include "MeshData.h"

namespace voxel2tet
{

enum MeshGenerator_Type {MG_TETGEN};

class MeshGenerator3D
{
public:
    MeshGenerator3D();
    MeshData *Mesh;

    virtual void Execute()=0;
};

}
#endif // MESHGENERATOR3D_H
