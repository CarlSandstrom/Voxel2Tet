#ifndef TETGENCALLER_H
#define TETGENCALLER_H

#include "tetgen.h"
#include "MeshGenerator3D.h"

namespace voxel2tet
{

// Class for calling TetGen 3D mesh generator
class TetGenCaller: public MeshGenerator3D
{
public:
    TetGenCaller();
};

}
#endif // MESHGENERATOR3D_H
