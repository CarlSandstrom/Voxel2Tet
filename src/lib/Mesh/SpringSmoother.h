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
    /**
     * Constant that determines how far away from the original position a vertex can move
     */
    double c;

    /**
     * Used for determining c such that the largest distance is a factor of the side length of a voxel
     */
    double c_factor;

    /**
     * Determines how hard movement outside the range determined by c should be
     */
    double alpha;

    /**
     * Characteristic length of a voxel
     */
    double charlength;

    /**
     * Compute out-of-balance vector for a vertex
     * @param ConnectionCoords List of coordinates for connected vertices
     * @param xc Vertex current coordinate
     * @param x0 Original position for vertex
     * @param alpha See SpringSmoother::alpha
     * @param c See SpringSmoother::c
     * @return
     */
    virtual arma::vec
    ComputeOutOfBalance(std::vector<arma::vec3> ConnectionCoords, arma::vec3 xc, arma::vec3 x0, double alpha, double c);

    /**
     * Computes the tangent for the out-of-balance force.
     * @param ConnectionCoords List of coordinates for connected vertices
     * @param xc Vertex current coordinate
     * @param x0 Original position for vertex
     * @param alpha See SpringSmoother::alpha
     * @param c See SpringSmoother::c
     * @return Tangent
     */
    virtual arma::mat
    ComputeAnalyticalTangent(std::vector<arma::vec3> ConnectionCoords, arma::vec xc, arma::vec x0, double alpha,
                             double c);

    /**
     * Computes the tangent for the out-of-balance force using numerical differentiation
     * @param ConnectionCoords List of coordinates for connected vertices
     * @param xc Vertex current coordinate
     * @param x0 Original position for vertex
     * @param alpha See SpringSmoother::alpha
     * @param c See SpringSmoother::c
     * @return Tangent
     */
    virtual arma::mat
    ComputeNumericalTangent(std::vector<arma::vec3> ConnectionCoords, arma::vec xc, arma::vec x0, double alpha,
                            double c);
private:

    /**
     * Computes c such that the voxel can move a distance l from its original position
     * @param l See SpringSmoother::charlength
     * @param alpha See SpringSmoother::alpha
     * @return Value of c
     */
    double Compute_c(double l, double alpha);


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

    /**
     * Perform smoothing on Vertices
     * @param Vertices List of vertices that should be smoothed
     * @param Mesh Pointer to Mesh object.
     */
    void Smooth(std::vector<VertexType *> Vertices, MeshData *Mesh = NULL);
};

}
#endif
