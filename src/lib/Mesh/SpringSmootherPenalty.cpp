#include "SpringSmootherPenalty.h"

namespace voxel2tet
{

SpringSmootherPenalty::SpringSmootherPenalty(double VoxelCharLength, double c, double alpha, double c_factor, std::vector<Volume *> *Volumes,
                                             bool compute_c) : SpringSmoother(VoxelCharLength, c, alpha, c_factor,
                                                                              compute_c)
{
    this->Volumes = Volumes;
}

arma::vec SpringSmootherPenalty::ComputeOutOfBalance(std::vector<arma::vec3> ConnectionCoords, arma::vec3 xc, arma::vec3 x0, double alpha, double c, VertexType *v)
{
    arma::vec F = {
        0., 0., 0.
    };

    // Compute nonlinear part of force
    arma::vec n0;

    double d0 = arma::norm(x0 - xc);
    if (d0 < 1e-8) {
        n0 = {
            0., 0., 0.
        };
    } else {
        n0 = (x0 - xc) / d0;
    }

    F = (exp(pow(d0 / c, alpha)) - 1) * n0;

    // Compute linear part of force
    for (unsigned int i = 0; i < ConnectionCoords.size(); i++) {
        arma::vec xi = {
            ConnectionCoords.at(i)[0], ConnectionCoords.at(i)[1], ConnectionCoords.at(i)[2]
        };
        arma::vec nj;
        double dj = arma::norm(xi - xc);
        if (dj < 1e-8) {
            nj = {
                0., 0., 0.,
            };
        } else {
            nj = (xi - xc) / dj;
        }
        F = F + dj * nj / ConnectionCoords.size();
    }

    // Compute volume part

    for (Volume *V: this->VertexVolumes[v]) {
        double CurrentVolume = VolumeCurrent[V];
        double Coeff = this->VolumePenalty[V]*(this->VolumeZero[V]-CurrentVolume);

        // Compute derivative
        double eps = 1e-8;
        arma::vec dVdx = {0.,0.,0.};

/*        std::vector<TriangleType *> VertexTriangles = v->Triangles;

        for (TriangleType *t: VertexTriangles) {
            std::vector<VertexType *> Vertices;

            for (int i=0; i<3; i++) {
                if (t->Vertices[i] != v) {
                    Vertices.push_back(t->Vertices[i]);
                }
            }

            arma::vec tContribution = {Vertices[0]->get_c(1)*Vertices[1]->get_c(2)-Vertices[1]->get_c(1)*Vertices[0]->get_c(2),
                                       -Vertices[0]->get_c(0)*Vertices[1]->get_c(2)+Vertices[1]->get_c(0)*Vertices[0]->get_c(2),
                                       Vertices[0]->get_c(0)*Vertices[1]->get_c(1)-Vertices[1]->get_c(0)*Vertices[0]->get_c(1)};

            if (V->Phase==t->NegNormalMatID) {
                tContribution = -tContribution;
            }

            //LOG("tContribution = [%f, %f, %f]\n", tContribution[0], tContribution[1], tContribution[2]);
            // Compare MatID of volume with PosNormMatID
            F = F + Coeff*tContribution;

        }*/


        for (int i=0; i<3; i++) {
            double oldvalue = v->get_c(i);
            v->set_c(oldvalue+eps, i);
            double Veps = V->ComputeVolume();
            dVdx[i] = (Veps - CurrentVolume) / eps;
            if (dVdx[i]>1000000) {
                LOG ("Fishy...\n", 0);
            }
            v->set_c(oldvalue, i);
        }
        F = F + dVdx*Coeff;
        // LOG("Coeff = %f, dVdx=[%f, %f, %f]\n", Coeff, dVdx[0], dVdx[1], dVdx[2]);


    }

    LOG ("F = [%f, %f, %f]\n", F[0], F[1], F[2]);

/*    for (int i=0; i<3; i++) {
        if (v->Fixed[i]) {
            F[i] = 0.0;
        }
    }*/

    return F;
}

arma::mat SpringSmootherPenalty::ComputeNumericalTangent(std::vector<arma::vec3> ConnectionCoords, arma::vec xc, arma::vec x0,
                                                         double alpha, double c, VertexType *v)
{
    double eps = 1e-10;
    arma::mat Tangent = arma::zeros<arma::mat>(3, 3);

    arma::vec Fval = ComputeOutOfBalance(ConnectionCoords, xc, x0, alpha, c, v);

    for (int i = 0; i < 3; i++) {
        arma::vec xi = xc;
        xi[i] = xi[i] + eps;
        double oldc = v->get_c(i);
        v->set_c(xi[i]+eps, i);
        arma::vec Fvali = ComputeOutOfBalance(ConnectionCoords, xi, x0, alpha, c, v);
        arma::vec dF = (Fvali - Fval) / eps;
        Tangent.col(i) = dF;
        v->set_c(oldc, i);
    }

    return Tangent;
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

    // Create mapping Vertices <-> Volumes
    VertexVolumes.clear();
    for (Volume* V: *this->Volumes) {
        for (Surface* S: V->Surfaces) {
            for (VertexType* v: S->Vertices) {
                VertexVolumes[v].insert(V);
            }
        }
    }

    // Compute all original volumes
    VolumeZero.clear();
    for (Volume *V: *this->Volumes) {
        VolumeZero[V] = V->ComputeVolume();
    }

    this->CheckPenetration(&Vertices, (MeshManipulations *) Mesh);


    for (double p=0; p<100; p++) {

        // Initialize penalties
        VolumePenalty.clear();
        for (Volume *V: *this->Volumes) {
            VolumePenalty[V] = p*100;
        }

        double deltamax = 1e8;
        int iter=0;

        while (deltamax > MAXCHANGE) {

            deltamax = 0.0;
            size_t j = 0;

            // Compute current volumes
            for (Volume *V: *this->Volumes) {
                VolumeCurrent[V] = V->ComputeVolume();
            }

            for (VertexType *v: Vertices) {

                std::vector<arma::vec3> ConnectionCoords;
                for (VertexType *cv: Connections[j]) {
                    ConnectionCoords.push_back(cv->get_c_vec());
                }

                arma::vec3 xc = CurrentPositions[v];
                arma::vec3 x0 = OriginalPositions[v];
                arma::vec3 R = ComputeOutOfBalance(ConnectionCoords, xc, x0, alpha, c, v);

                double err = arma::norm(R);

                while (err > 1e-5) {
                    arma::mat K = ComputeNumericalTangent(ConnectionCoords, xc, x0, alpha, c, v);
                    arma::vec d = -arma::solve(K, R);
                    xc = xc + d;
                    R = ComputeOutOfBalance(ConnectionCoords, xc, x0, alpha, c, v);
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

}

std::string SpringSmootherPenalty::DoOutput() const
{
    std::string stream;
    stream = "\talpha = " + std::to_string(alpha) + ", " +
            "c = " + std::to_string(c) + ", c_factor = " + std::to_string(c_factor) + "\n";
    return stream;
}

}
