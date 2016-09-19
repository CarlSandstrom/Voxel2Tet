#include "LaplaceSmoother.h"

namespace voxel2tet
{

LaplaceSmoother::LaplaceSmoother(double VoxelCharLength, double c, double alpha, double c_factor, bool compute_c) : SpringSmoother( VoxelCharLength, c, alpha, c_factor, compute_c)
{

}

void LaplaceSmoother :: Smooth(std :: vector< VertexType * >Vertices, MeshData *Mesh)
{
    double MAXCHANGE = 1e-3 * charlength;

    std::vector<std::vector<VertexType *>> Connections = this->GetConnectivityVector(Vertices);

    // Create vectors for current and previous positions for all involved vertices (even those vertices connected to a vertex in Vertices vector)
    std :: map<VertexType *, arma::vec3 >OriginalPositions;
    std :: map<VertexType *, arma::vec3 >CurrentPositions;
    std :: map<VertexType *, arma::vec3 >PreviousPositions;

    for (std::vector<VertexType *> VertexList: Connections) {
        for (VertexType *v: VertexList) {
            OriginalPositions[v] = v->get_c_vec();
            CurrentPositions[v] = v->get_c_vec();
            PreviousPositions[v] = v->get_c_vec();
        }
    }

    for (VertexType *v: Vertices) {
        OriginalPositions[v] = v->get_c_vec();
        CurrentPositions[v] = v->get_c_vec();
        PreviousPositions[v] = v->get_c_vec();
    }

    bool intersecting = true;
    this->CheckPenetration(&Vertices, (MeshManipulations*) Mesh);

    while (intersecting) {

        double deltamax = 1e8;
        intersecting = false;

        while (deltamax > MAXCHANGE) {

            deltamax = 0.0;
            size_t i=0;

            for (VertexType *v: Vertices) {

                std::vector<arma::vec3> ConnectionCoords;
                for (VertexType *cv: Connections[i]) {
                    ConnectionCoords.push_back(cv->get_c_vec());
                }

                arma::vec3 xc = CurrentPositions[v];
                arma::vec3 x0 = OriginalPositions[v];

                arma::vec3 NewPosition = {0, 0, 0};

                for (VertexType *cv: Connections[i]) {
                    NewPosition = NewPosition + 1.0 / double ( Connections[i].size() ) * cv->get_c_vec();
                }
                
                for (int i=0; i<3; i++) {
                    if (!v->Fixed[i]) {
                        CurrentPositions[v][i] = NewPosition[i];
                        v->set_c(CurrentPositions[v][i], i);
                    }
                }

                double delta = arma::norm(CurrentPositions[v]-PreviousPositions[v]);

                for (int i=0; i<3; i++) {
                    if (!v->Fixed[i]) {
                        PreviousPositions[v][i] = CurrentPositions[v][i];
                    }
                }

                arma::vec3 R = ComputeOutOfBalance(ConnectionCoords, xc, x0, alpha, c);

                deltamax = std::max(delta, deltamax);
                i++;
            }
        }
    }
}

}
