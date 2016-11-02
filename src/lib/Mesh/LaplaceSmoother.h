#ifndef LaplaceSmoother_H
#define LaplaceSmoother_H

#include "SpringSmoother.h"

namespace voxel2tet
{

class LaplaceSmoother : public SpringSmoother
{
public:
    LaplaceSmoother(double VoxelCharLength, double c, double alpha, double c_factor, bool compute_c = false);

    ~LaplaceSmoother()
    {}

    virtual void Smooth(std::vector<VertexType *> Vertices, MeshData *Mesh = NULL);
};

}
#endif // LaplaceSmoother_H
