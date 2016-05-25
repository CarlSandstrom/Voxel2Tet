#include <omp.h>
#include"Smoother.h"

namespace voxel2tet
{

arma::vec ComputeOutOfBalance(std::vector<std::array<double, 3> > ConnectionCoords, arma::vec xc, arma::vec x0, double alpha, double c)
{

    arma::vec F = {0., 0., 0.};

    // Compute nonlinear part of force
    arma::vec n0;

    double d0 = arma::norm(x0-xc);
    if (d0<1e-8) {
        n0 = {0.,0.,0.};
    } else {
        n0 = (x0-xc)/d0;
    }

    F = (exp(pow(d0/c, alpha))-1) * n0;

    // Compute linear part of force
    for (unsigned int i=0; i<ConnectionCoords.size(); i++) {
        arma::vec xi = {ConnectionCoords.at(i)[0], ConnectionCoords.at(i)[1], ConnectionCoords.at(i)[2]};
        arma::vec nj;
        double dj = arma::norm(xi-xc);
        if (dj<1e-8) {
            nj = {0., 0., 0.,};
        } else {
            nj = (xi-xc)/dj;
        }
        F = F + dj*nj/ConnectionCoords.size();
    }

    return F;

}

arma::mat  ComputeNumericalTangent(std::vector<std::array<double, 3> > ConnectionCoords, arma::vec xc, arma::vec x0, double alpha, double c)
{

    double eps=1e-10;
    arma::mat Tangent=arma::zeros<arma::mat>(3,3);

    arma::vec Fval=ComputeOutOfBalance(ConnectionCoords, xc, x0, alpha, c);

    for (int i=0; i<3; i++) {
        arma::vec xi = xc;
        xi[i] = xi[i] + eps;
        arma::vec Fvali=ComputeOutOfBalance(ConnectionCoords, xi, x0, alpha, c);
        arma::vec dF = (Fvali-Fval)/eps;
        Tangent.col(i) = dF;
    }

    return Tangent;

}

arma::mat ComputeAnalyticalTangent(std::vector<std::array<double, 3> > ConnectionCoords, arma::vec xc, arma::vec x0, double alpha, double c)
{

    arma::mat Tangent=arma::zeros<arma::mat>(3,ConnectionCoords.size()*3+3);
    arma::mat TangentSelf=arma::zeros<arma::mat>(3,3);

    arma::vec a0 = x0-xc;
    double d0 = arma::norm(a0);
    arma::vec n0;

    // delta x_i part

    // Non-linear part
    // We run into numerical trouble if d0=0. However, in the case of d0=0, everyting nonlinear is zero...
    arma::vec Dexp;
    arma::mat Dn0;
    if (d0 < 1e-8) {
        n0 = {0.,0.,0.};
        Dexp = {0.,0.,0.};
        Dn0 = -arma::zeros<arma::mat>(3,3);

    } else {
        n0 = a0/d0;
        Dexp = -std::exp(std::pow(d0/c, alpha))*alpha/c*pow(d0/c, alpha-1)*n0;
        Dn0 = -arma::eye<arma::mat>(3,3)/d0 + a0*a0.t()/pow(d0,3);
    }
    TangentSelf = n0*Dexp.t() + Dn0 * ( std::exp(std::pow(d0/c, alpha)) - 1);


    // Linear part
    TangentSelf = TangentSelf - arma::eye<arma::mat>(3,3);

    // Assemble self part to tangent
    Tangent(arma::span(0,2),arma::span(0,2)) = TangentSelf;

    // delta x_j part
    for (unsigned int i=0; i<ConnectionCoords.size(); i++) {
        for (int j=0; j<3; j++) {
            Tangent(j, 3*(i+1)+j) = 1.0/ConnectionCoords.size();
        }
    }

    return Tangent;
}

void SpringSmooth(std::vector<VertexType*> Vertices, std::vector<bool> Fixed, std::vector<std::vector<VertexType*>> Connections, double K, voxel2tet::MeshData *Mesh)
{
    int MAX_ITER_COUNT=10000;
    double NEWTON_TOL=1e-10;
    double alpha = 4;

    // Create vectors for current and previous positions
    std::vector<std::array<double, 3>> CurrentPositions;
    std::vector<std::array<int, 3>> dofids;

    int ndof = 0;

    for (unsigned int i=0; i<Vertices.size(); i++) {
        std::array<double, 3> cp;
        std::array<int, 3> dofid = {-1, -1, -1};
        for (int j=0; j<3; j++) {
            cp.at(j) = Vertices.at(i)->get_c(j);
            if (!Fixed[i]) {
                if (!Vertices.at(i)->Fixed[j]) {
                    dofid[j]=ndof;
                    ndof++;
                }
            }
        }
        CurrentPositions.push_back(cp);
        dofids.push_back(dofid);
    }

    // Create vector-vector for accessing neightbours
    std::vector <std::vector <int>> ConnectionVertexIndex;
    for (std::vector<VertexType*> n: Connections) {
        ConnectionVertexIndex.push_back({});
        for (VertexType *v: n) {
            int VertexIndex = std::distance(Vertices.begin(), std::find(Vertices.begin(), Vertices.end(), v));
            ConnectionVertexIndex.at(ConnectionVertexIndex.size()-1).push_back(VertexIndex);
        }
    }

#if EXPORT_SMOOTHING_ANIMATION == 1
    std::ostringstream FileName;
    if (Mesh!=NULL) {
        FileName << "/tmp/Smoothing" << 0 << ".vtp";
        Mesh->ExportSurface( FileName.str(), FT_VTK );
    }
#endif

    int itercount = 0;
    double err = 1e8;

    arma::sp_mat Kff(ndof, ndof);
    arma::vec Rf = arma::ones<arma::vec>(ndof);

    while ((itercount < MAX_ITER_COUNT) & (err>1e-6)) {

        err = 0.0;
        Kff.zeros();

        for (size_t i=0; i<Vertices.size(); i++) {
            if (!Fixed[i]) {

                std::vector<std::array<double, 3>> ConnectionCoords;
                std::vector<VertexType*> MyConnections = Connections.at(i);

                for (unsigned k=0; k<MyConnections.size(); k++) {
                    ConnectionCoords.push_back({CurrentPositions.at(ConnectionVertexIndex.at(i).at(k))[0],
                                                CurrentPositions.at(ConnectionVertexIndex.at(i).at(k))[1],
                                                CurrentPositions.at(ConnectionVertexIndex.at(i).at(k))[2]});
                }

                arma::vec xc = {CurrentPositions.at(i)[0], CurrentPositions.at(i)[1], CurrentPositions.at(i)[2]};
                arma::vec x0 = {Vertices.at(i)->get_c(0), Vertices.at(i)->get_c(1), Vertices.at(i)->get_c(2)};

                // Assemble out-of-balance vector
                arma::vec R = ComputeOutOfBalance(ConnectionCoords, xc, x0, alpha, K);

                for (int j=0; j<3; j++) {
                    if (dofids[i][j]!=-1) { // Not stationary
                        Rf[dofids[i][j]] = R[j];
                    }
                }

                // Assemble to tangent
                std::vector <int> ConnectionIndices = ConnectionVertexIndex.at(i);
                arma::mat T = ComputeAnalyticalTangent(ConnectionCoords, xc, x0, alpha, K);

                std::vector <int> Rows = {dofids[i][0], dofids[i][1], dofids[i][2]};
                std::vector <int> Cols = {dofids[i][0], dofids[i][1], dofids[i][2]};

                for (int index: ConnectionIndices) {
                    for (int j=0; j<3; j++) {
                        Cols.push_back(dofids[index][j]);
                    }
                }

                for (int j=0; j<3; j++) {
                    if (Rows[j]!=-1) {
                        for (int k=0; k<Cols.size(); k++) {
                            if (Cols[k]!=-1) {
                                Kff(Rows[j], Cols[k]) = Kff(Rows[j], Cols[k]) + T(j, k);
                            }
                        }
                    }
                }
            }
        }

        err = arma::norm(Rf);
        //Kff.print("Kff = ");

        arma::vec delta = arma::spsolve(Kff, -Rf);
        if (err>1e-10) {
            double maxdelta = arma::max(delta);
            if (maxdelta > 0.066667) {
                delta = delta * 0.066667/maxdelta;
            }
        }

        //delta.print("update = ");

        // Update current positions
        for (size_t i=0; i<Vertices.size(); i++) {
            for (int j=0; j<3; j++) {
                if (dofids[i][j]!=-1) {
                    CurrentPositions.at(i)[j] = CurrentPositions.at(i)[j] + delta[dofids[i][j]];
                }
            }
        }

        STATUS("Iteration %i end with err=%f\n", itercount, err);

        itercount ++;

#if EXPORT_SMOOTHING_ANIMATION == 1
        // ************************** DEBUG STUFF
        // Update vertices
        for (unsigned int i=0; i<Vertices.size(); i++) {
            VertexType *v=Vertices.at(i);
            for (int j=0; j<3; j++) {
                if (!v->Fixed[j]) {
                    v->set_c(CurrentPositions.at(i)[j], j);
                }
            }
        }

        if (Mesh!=NULL) {
            FileName.str(""); FileName.clear();
            FileName << "/tmp/Smoothing" << itercount++ << ".vtp";
            Mesh->ExportSurface(FileName.str(), FT_VTK);
        }
        // ************************** /DEBUG STUFF
#endif

#if TEST_MESH_FOR_EACH_SMOOTHING
        TetGenCaller Tetgen;
        Tetgen.Mesh = Mesh;
        Tetgen.TestMesh();
#endif
    }

    if (itercount > 999) {
        STATUS("WARNING: Smoothing did not converge\n", 0);
    }

    // Update vertices
    for (unsigned int i=0; i<Vertices.size(); i++) {
        for (int j=0; j<3; j++) {
            Vertices.at(i)->set_c(CurrentPositions.at(i)[j], j); // TODO: Use array to improve performance
        }
    }

}

}
