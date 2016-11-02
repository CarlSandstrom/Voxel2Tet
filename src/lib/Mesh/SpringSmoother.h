#ifndef SPRINGSMOOTHER_H
#define SPRINGSMOOTHER_H

#include <vector>
#include <math.h>
#include <sstream>
#include <armadillo>

#include "MeshComponents.h"
#include "MeshData.h"
#include "TetGenCaller.h"
#include "Options.h"
#include "Smoother.h"

namespace voxel2tet
{

class SmootherConfiguration
{

};

/**
 * @brief Provides functionality for smoothing a set of connected vertices given their relation and connectivity.
 *
 * This particular smoothing class uses spring attached between all connected vertices. For each vertex, an additional
 * (non-linear) spring is also connected to its original location. This limits the decrease in volume of the objects.
 *
 */
class SpringSmoother : public Smoother
{
protected:
    double c;
    double c_factor;
    double alpha;
    double charlength;

    arma::vec
    ComputeOutOfBalance(std::vector<arma::vec3> ConnectionCoords, arma::vec3 xc, arma::vec3 x0, double alpha, double c);

    arma::mat
    ComputeAnalyticalTangent(std::vector<arma::vec3> ConnectionCoords, arma::vec xc, arma::vec x0, double alpha,
                             double c);

private:

    double Compute_c(double l, double alpha);

    arma::mat
    ComputeNumericalTangent(std::vector<arma::vec3> ConnectionCoords, arma::vec xc, arma::vec x0, double alpha,
                            double c);

protected:
    std::string DoOutput() const;

public:

    /**
     * @brief Constructor given data for the spring smoothing method
     * @param VoxelCharLength Charachteristic length of a voxel
     * @param c Constant
     * @param alpha Constant
     * @param c_factor If c is computed, use a factor c_factor of the VoxelCharLength when determining the equilibrium
     * @param compute_c Determines if c should be used from input or computed
     */
    SpringSmoother(double VoxelCharLength, double c, double alpha, double c_factor, bool compute_c = false);

    ~SpringSmoother()
    {}

    void Smooth(std::vector<VertexType *> Vertices, MeshData *Mesh = NULL);
};

}
#endif
