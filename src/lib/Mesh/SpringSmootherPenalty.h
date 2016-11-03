#ifndef SPRINGSMOOTHERPENALTY_H
#define SPRINGSMOOTHERPENALTY_H

#include "SpringSmoother.h"

namespace voxel2tet
{

/**
 * Performs smoothing of voxel data. The smoothing is the same as in SpringSmoother but with the difference that the
 * volume is preserved using the penalty method (Experimental)
 */
class SpringSmootherPenalty : public SpringSmoother
{
public:

    /**
     * @brief Constructor given data for the spring smoothing method
     * @param VoxelCharLength Charachteristic length of a voxel
     * @param c Constant
     * @param alpha Constant
     * @param c_factor If c is computed, use a factor c_factor of the VoxelCharLength when determining the equilibrium
     * @param compute_c Determines if c should be used from input or computed
     */
    SpringSmootherPenalty(double VoxelCharLength, double c, double alpha, double c_factor, bool compute_c = false);

    ~SpringSmootherPenalty()
    {}

    virtual void Smooth(std::vector<VertexType *> Vertices, MeshData *Mesh = NULL);

protected:
    std::string DoOutput() const;
};

}
#endif // SPRINGSMOOTHERPENALTY_H
