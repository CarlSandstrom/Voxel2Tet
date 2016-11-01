#ifndef SPRINGSMOOTHERPENALTY_H
#define SPRINGSMOOTHERPENALTY_H

#include "SpringSmoother.h"

namespace voxel2tet
{

class SpringSmootherPenalty : public SpringSmoother
{
public:
    SpringSmootherPenalty(double VoxelCharLength, double c, double alpha, double c_factor, bool compute_c=false);
    ~SpringSmootherPenalty() {}
    virtual void Smooth(std :: vector< VertexType * >Vertices, MeshData *Mesh = NULL);
protected:
    std::string DoOutput() const;
};

}
#endif // SPRINGSMOOTHERPENALTY_H
