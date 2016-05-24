#include <vector>
#include <math.h>
#include <sstream>
#include <armadillo>

#include "MeshComponents.h"
#include "MeshData.h"
#include "TetGenCaller.h"

namespace voxel2tet
{

//void ComputeOutOfBalance(std::vector<std::array<double, 3>> ConnectionCoords, std::array<double,3> x0, std::array<double,3> xc, double alpha, double c);
arma::vec ComputeOutOfBalance(std::vector<std::array<double, 3> > ConnectionCoords, arma::vec xc, arma::vec x0, double alpha, double c);

arma::mat ComputeNumericalTangent(std::vector<std::array<double, 3> > ConnectionCoords, arma::vec xc, arma::vec x0, double alpha, double c);

arma::mat ComputeAnalyticalTangent(std::vector<std::array<double, 3> > ConnectionCoords, arma::vec xc, arma::vec x0, double alpha, double c);

void SpringSmooth(std::vector<VertexType*> Vertices, std::vector<bool> Fixed,
                  std::vector<std::vector<VertexType*>> Connections, double K, MeshData *Mesh=NULL);


}
