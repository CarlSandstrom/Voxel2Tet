#ifndef SMOOTHER_H
#define SMOOTHER_H

#include <vector>
#include <math.h>
#include <sstream>
#include <armadillo>

#include "MeshComponents.h"
#include "MeshData.h"
#include "TetGenCaller.h"

namespace voxel2tet
{

double Compute_c(double l, double alpha);

arma::vec ComputeOutOfBalance(std::vector<std::array<double, 3> > ConnectionCoords, arma::vec xc, arma::vec x0, double alpha, double c);

arma::mat ComputeNumericalTangent(std::vector<std::array<double, 3> > ConnectionCoords, arma::vec xc, arma::vec x0, double alpha, double c);

arma::mat ComputeAnalyticalTangent(std::vector<std::array<double, 3> > ConnectionCoords, arma::vec xc, arma::vec x0, double alpha, double c);

void SpringSmooth(std::vector<VertexType*> Vertices, std::vector<bool> Fixed,
                  std::vector<std::vector<VertexType*>> Connections,
                  double c, double alpha, double charlength, bool Automatic_c=false,
                  MeshData *Mesh=NULL);

}
#endif
