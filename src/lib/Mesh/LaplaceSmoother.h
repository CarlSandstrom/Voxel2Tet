#ifndef LaplaceSmoother_H
#define LaplaceSmoother_H

#include "SpringSmoother.h"

namespace voxel2tet
{

/**
 * Class for Laplace smoothing
 */
class LaplaceSmoother : public SpringSmoother
{
public:
    /**
     * Constructor for the Laplace smoothing class
     * @param VoxelCharLength
     * @param c
     * @param alpha
     * @param c_factor
     * @param compute_c
     * @return
     */
    LaplaceSmoother(double VoxelCharLength, double c, double alpha, double c_factor, bool compute_c = false);

    ~LaplaceSmoother()
    {}

    virtual void Smooth(std::vector<VertexType *> Vertices, MeshData *Mesh = NULL);
};

}
#endif // LaplaceSmoother_H
