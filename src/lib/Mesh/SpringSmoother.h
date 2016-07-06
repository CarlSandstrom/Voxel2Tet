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
    friend std::ostream &operator<<(std::ostream &stream, const SpringSmoother &Smoother);
private:

    double c;
    double c_factor;
    double alpha;
    double charlength;

    double Compute_c(double l, double alpha);
    arma :: vec ComputeOutOfBalance(std :: vector< std :: array< double, 3 > >ConnectionCoords, arma :: vec xc, arma :: vec x0, double alpha, double c);
    arma :: mat ComputeNumericalTangent(std :: vector< std :: array< double, 3 > >ConnectionCoords, arma :: vec xc, arma :: vec x0, double alpha, double c);
    arma :: mat ComputeAnalyticalTangent(std :: vector< std :: array< double, 3 > >ConnectionCoords, arma :: vec xc, arma :: vec x0, double alpha, double c);
    arma :: mat ComputeAnalyticalTangentGlobal(std :: vector< std :: array< double, 3 > >ConnectionCoords, arma :: vec xc, arma :: vec x0, double alpha, double c);

    void SpringSmoothGlobal(std :: vector< VertexType * >Vertices, std :: vector< bool >Fixed,
                            std :: vector< std :: vector< VertexType * > >Connections,
                            double c, double alpha, double charlength, bool Automatic_c = false,
                            MeshData *Mesh = NULL);

    std :: vector< std :: pair< TriangleType *, TriangleType * > >CheckPenetration(std :: vector< VertexType * > *Vertices, MeshData *Mesh);

public:

    SpringSmoother(double VoxelCharLength, double c, double alpha, double c_factor, bool compute_c=false);
    ~SpringSmoother() {}


    void Smooth(std :: vector< VertexType * >Vertices, std :: vector< bool >Fixed,
                      std :: vector< std :: vector< VertexType * > >Connections,
                      MeshData *Mesh = NULL);
};

}
#endif
