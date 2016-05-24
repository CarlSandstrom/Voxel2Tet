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
        F = F + dj*nj;
    }

    return F;

}

arma::mat  ComputeNumericalTangent(std::vector<std::array<double, 3> > ConnectionCoords, arma::vec xc, arma::vec x0, double alpha, double c)
{

    double eps=1e-8;
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

void SpringSmooth(std::vector<VertexType*> Vertices, std::vector<bool> Fixed, std::vector<std::vector<VertexType*>> Connections, double K, voxel2tet::MeshData *Mesh)
{
    int MAX_ITER_COUNT=10000;
    double NEWTON_TOL=1e-10;
    double alpha = 4;

    // Create vectors for current and previous positions
    std::vector<std::array<double, 3>> CurrentPositions;
    std::vector<std::array<double, 3>> PreviousPositions;

    for (unsigned int i=0; i<Vertices.size(); i++) {
        std::array<double, 3> cp;
        std::array<double, 3> pp;
        for (int j=0; j<3; j++) {
            cp.at(j) = pp.at(j) = Vertices.at(i)->get_c(j);
        }
        CurrentPositions.push_back(cp);
        PreviousPositions.push_back(pp);
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
    int threadcount = 1;
    double deltamax = 1e8;
    int deltamaxnode;

#ifdef OPENMP
    threadcount = omp_get_max_threads();
#endif

    while ((itercount < MAX_ITER_COUNT) & (deltamax>1e-6)) {

        deltamax=0.0;
        int deltamaxnodes[threadcount];
        double deltamaxvalues[threadcount];

        for (int i=0; i<threadcount; i++) {
            deltamaxnodes[i]=0;
            deltamaxvalues[i]=0.0;
        }

        //#pragma omp parallel default(shared)
        {
            int threadid=0;
#ifdef OPENMP
            threadid=omp_get_thread_num();
#endif

            //#pragma omp for schedule(static, 100)
            for (size_t i=0; i<Vertices.size(); i++) {
                if (!Fixed[i]) {

                    std::vector<std::array<double, 3>> ConnectionCoords;
                    std::vector<VertexType*> MyConnections = Connections.at(i);

                    for (unsigned k=0; k<MyConnections.size(); k++) {
                        ConnectionCoords.push_back({PreviousPositions.at(ConnectionVertexIndex.at(i).at(k))[0],
                                                    PreviousPositions.at(ConnectionVertexIndex.at(i).at(k))[1],
                                                    PreviousPositions.at(ConnectionVertexIndex.at(i).at(k))[2]});
                    }

                    arma::vec xc = {PreviousPositions.at(i)[0], PreviousPositions.at(i)[1], PreviousPositions.at(i)[2]};
                    arma::vec x0 = {Vertices.at(i)->get_c(0), Vertices.at(i)->get_c(1), Vertices.at(i)->get_c(2)};

                    // Compute out-of-balance vector
                    arma::vec R = ComputeOutOfBalance(ConnectionCoords, xc, x0, alpha, K);
                    double err = arma::norm(R);
                    int iter = 0;

                    // Find equilibrium
                    while ((err>1e-5) & (iter<1000)) {
                        arma::mat T = ComputeNumericalTangent(ConnectionCoords, xc, x0, alpha, K);
                        arma::vec delta = -arma::solve(T, R);
                        xc = xc + 0.1*delta;
                        R = ComputeOutOfBalance(ConnectionCoords, xc, x0, alpha, K);
                        err = arma::norm(R);
                        iter ++;
                    }

                    // If too many iterations, throw an exception and investigate...
                    if (iter>999) {
                        throw(0);
                    }

                    // Update current position
                    for (int j=0; j<3; j++) {
                        if (!Vertices.at(i)->Fixed[j]) {
                            CurrentPositions.at(i)[j] = xc[j];
                        }
                    }

                    // Update maximum delta
                    arma::vec d={PreviousPositions.at(i)[0]-CurrentPositions.at(i)[0], PreviousPositions.at(i)[1]-CurrentPositions.at(i)[1], PreviousPositions.at(i)[2]-CurrentPositions.at(i)[2]};
                    if (arma::norm(d)>deltamaxvalues[threadid]) {
                        deltamaxvalues[threadid] = arma::norm(d);
                        deltamaxnodes[threadid] = i;
                    }

                }
            }
            // Update previous positions
            for (size_t i=0; i<Vertices.size(); i++) {
                for (int j=0; j<3; j++) {
                    if (!Vertices.at(i)->Fixed[j]) {
                        PreviousPositions.at(i)[j] = CurrentPositions.at(i)[j];
                    }
                }
            }

        } // End of OpenMP section

        deltamaxnode=0;
        for (int i=0; i<threadcount; i++) {
            if (deltamaxvalues[i]>deltamax) {
                deltamax = deltamaxvalues[i];
                deltamaxnode = deltamaxnodes[i];
            }
        }

        STATUS("Iteration %i end with deltamax=%f at node %i\n", itercount, deltamax, deltamaxnode);

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
