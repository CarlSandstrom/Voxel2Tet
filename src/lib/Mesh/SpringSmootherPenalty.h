#ifndef SPRINGSMOOTHERPENALTY_H
#define SPRINGSMOOTHERPENALTY_H

#include "SpringSmoother.h"
#include "Volume.h"
#include <set>

namespace voxel2tet
{

/**
 * Performs smoothing of voxel data. The smoothing is the same as in SpringSmoother but with the difference that the
 * volume is preserved using the penalty method (Experimental)
 */
class SpringSmootherPenalty : public SpringSmoother
{
private:
    std::vector<Volume *> *Volumes;
    std::map<VertexType*, std::set<Volume*>> VertexVolumes;
    std::map<Volume *, double> VolumePenalty;
    std::map<Volume *, double> VolumeZero;
    std::map<Volume *, double> VolumeCurrent;
public:

    /**
     * @brief Constructor given data for the spring smoothing method
     * @param VoxelCharLength Charachteristic length of a voxel
     * @param c Constant
     * @param alpha Constant
     * @param c_factor If c is computed, use a factor c_factor of the VoxelCharLength when determining the equilibrium
     * @param compute_c Determines if c should be used from input or computed
     */
    SpringSmootherPenalty(double VoxelCharLength, double c, double alpha, double c_factor, std::vector<Volume *> *, bool compute_c = false);

    ~SpringSmootherPenalty()
    {}

    virtual void Smooth(std::vector<VertexType *> Vertices, MeshData *Mesh = NULL);

protected:
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
    ComputeOutOfBalance(std::vector<arma::vec3> ConnectionCoords, arma::vec3 xc, arma::vec3 x0, double alpha, double c, VertexType *v);

    virtual arma::mat ComputeNumericalTangent(std::vector<arma::vec3> ConnectionCoords, arma::vec xc, arma::vec x0,
                                                      double alpha, double c, VertexType *v);

    std::string DoOutput() const;
};

}
#endif // SPRINGSMOOTHERPENALTY_H
