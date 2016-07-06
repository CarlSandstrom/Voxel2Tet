#ifndef MESHGENERATOR3D_H
#define MESHGENERATOR3D_H

#include "MeshData.h"

namespace voxel2tet
{
enum MeshGenerator_Type { MG_TETGEN };

/**
 * @brief Abstract class for 3D volume mesh generators
 */
class MeshGenerator3D
{
public:
    /**
     * @brief Constructor
     */
    MeshGenerator3D();

    /**
     * @brief Pointer to a MeshData object containing the surface mesh
     */
    MeshData *Mesh;

    /**
     * @brief Calling this member function will execute the 3D mesh genetrator and return a tetrahedral mesh
     * @return Pointer MeshData object
     */
    virtual MeshData *Execute() = 0;
};
}
#endif // MESHGENERATOR3D_H
