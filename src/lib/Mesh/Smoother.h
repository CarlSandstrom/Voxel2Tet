#ifndef SMOOTHER_H
#define SMOOTHER_H

#include <vector>
#include <math.h>
#include <sstream>
#include <armadillo>

#include "MeshComponents.h"
#include "MeshData.h"
#include "TetGenCaller.h"
#include "Options.h"

namespace voxel2tet
{

class SmootherConfiguration
{

};

class SmootherClass
{
private:

    double c;
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
    SmootherClass(double VoxelCharLength, double c, double alpha, double c_factor, bool compute_c=false);
    ~SmootherClass() {}

    void SpringSmooth(std :: vector< VertexType * >Vertices, std :: vector< bool >Fixed,
                      std :: vector< std :: vector< VertexType * > >Connections,
                      MeshData *Mesh = NULL);
};

}
#endif
