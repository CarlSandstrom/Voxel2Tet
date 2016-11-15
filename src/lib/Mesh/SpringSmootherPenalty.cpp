#include "SpringSmootherPenalty.h"

namespace voxel2tet
{

SpringSmootherPenalty::SpringSmootherPenalty(double VoxelCharLength, double c, double alpha, double c_factor,
                                             bool compute_c) : SpringSmoother(VoxelCharLength, c, alpha, c_factor,
                                                                              compute_c)
{

}

void SpringSmootherPenalty::Smooth(std::vector<VertexType *> Vertices, MeshData *Mesh)
{

    STATUS("Smooth surface\n", 0);

    double MAXCHANGE = 1e-4 * charlength;

    std::vector<std::vector<VertexType *>> Connections = this->GetConnectivityVector(Vertices);

    // Create vectors for current and previous positions for all involved vertices (even those vertices connected to a vertex in Vertices vector)
    std::map<VertexType *, arma::vec3> OriginalPositions;
    std::map<VertexType *, arma::vec3> CurrentPositions;
    std::map<VertexType *, arma::vec3> PreviousPositions;

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

    this->CheckPenetration(&Vertices, (MeshManipulations *) Mesh);

    double deltamax = 1e8;
    int iter=0;

    while (deltamax > MAXCHANGE) {

        deltamax = 0.0;
        size_t j = 0;

        for (VertexType *v: Vertices) {

            std::vector<arma::vec3> ConnectionCoords;
            for (VertexType *cv: Connections[j]) {
                ConnectionCoords.push_back(cv->get_c_vec());
            }

            arma::vec3 xc = CurrentPositions[v];
            arma::vec3 x0 = OriginalPositions[v];
            arma::vec3 R = ComputeOutOfBalance(ConnectionCoords, xc, x0, alpha, c);

            double err = arma::norm(R);

            while (err > 1e-5) {
                arma::mat K = ComputeAnalyticalTangent(ConnectionCoords, xc, x0, alpha, c);
                arma::vec d = -arma::solve(K, R);
                xc = xc + d;
                R = ComputeOutOfBalance(ConnectionCoords, xc, x0, alpha, c);
                err = arma::norm(R);
            }


            for (int i = 0; i < 3; i++) {
                if (!v->Fixed[i]) {
                    CurrentPositions[v][i] = xc[i];
                    v->set_c(CurrentPositions[v][i], i);
                }
            }

            double delta = arma::norm(CurrentPositions[v] - PreviousPositions[v]);

            for (int i = 0; i < 3; i++) {
                if (!v->Fixed[i]) {
                    PreviousPositions[v][i] = CurrentPositions[v][i];
                }
            }

            deltamax = std::max(delta, deltamax);
            j++;
        }
        STATUS("%c[2K\r\tIteration %u end with deltamax=%f\r", 27, iter, deltamax);
        iter++;
    }
    STATUS("\n", 0);

}

std::string SpringSmootherPenalty::DoOutput() const
{
    std::string stream;
    stream = "\talpha = " + std::to_string(alpha) + ", " +
             "c = " + std::to_string(c) + ", c_factor = " + std::to_string(c_factor) + "\n";
    return stream;
}

}
