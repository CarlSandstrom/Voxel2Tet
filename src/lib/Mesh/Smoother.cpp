#include"Smoother.h"

namespace voxel2tet
{

void SpringSmooth(std::vector<VertexType*> Vertices, std::vector<std::array<bool,3>> FixedDirections, std::vector<std::vector<VertexType*>> Connections, double K, voxel2tet::MeshData *Mesh)
{
    int MAX_ITER_COUNT=10000;
    double NEWTON_TOL=1e-10;

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


    int itercount = 0;
    double deltamax = 1e8;

#if EXPORT_SMOOTHING_ANIMATION == 1
    std::ostringstream FileName;
    if (Mesh!=NULL) {
        FileName << "/tmp/Smoothing" << itercount << ".vtp";
        Mesh->ExportSurface( FileName.str(), FT_VTK );
    }
#endif

    while ((itercount < MAX_ITER_COUNT) & (deltamax>1e-4)){

        deltamax=0.0;
        int deltamaxnode = -1;
        unsigned int i;

#pragma omp parallel default(shared), private(i)
        {
#pragma omp for
        for (i=0; i<Vertices.size(); i++) {
            std::array<double, 3> NewCoords = {0,0,0};

            std::vector<VertexType*> MyConnections = Connections.at(i);
            for (unsigned int k=0; k<MyConnections.size(); k++) {
                // We need the index of MyConnection such that we can use CurrentPosition (as CurrentPosition contains the updated coordinates)

                //int VertexIndex = std::distance(Vertices.begin(), std::find(Vertices.begin(), Vertices.end(), MyConnections.at(k)));

                for (int m=0; m<3; m++) {
                    NewCoords.at(m) = NewCoords.at(m) + CurrentPositions.at(ConnectionVertexIndex.at(i).at(k))[m] / double(MyConnections.size());
                }
            }

            for (int j=0; j<3; j++) {
                if (!FixedDirections.at(i)[j]) {
                    CurrentPositions.at(i)[j] = NewCoords[j];
                }
            }

            // Pull back
            double d;
            std::array<double, 3> delta, unitdelta;
            for (int j=0; j<3; j++) delta[j]=CurrentPositions.at(i)[j]-Vertices.at(i)->originalcoordinates[j];

            double d0 = sqrt( pow(delta[0], 2) + pow(delta[1], 2) + pow(delta[2], 2) );
            if ((d0>1e-8) & (K!=0.0)) {
                for (int j=0; j<3; j++) unitdelta[j] = delta[j]/d0;

                double F = d0*1.0;

                int SolveIterations = 0;

                d=0;
                double fval=F-d*exp(d*d/K);

                while (std::fabs(fval)>NEWTON_TOL) {
                    double update = fval/(std::exp(d*d/K)*(1+2*d*d/K));
                    d=d+update;
                    fval=F-d*exp(d*d/K);
                    SolveIterations++;
                }

                if (SolveIterations == 100) {
                    STATUS("WARNING: Newton iterations did not converge\n", 0);
                    d=0.0;
                }

            }

            // Update current position with new delta
            for (int j=0; j<3; j++) {
                if (!FixedDirections.at(i)[j]) {
                    CurrentPositions.at(i)[j] = Vertices.at(i)->originalcoordinates[j] + unitdelta[j]*d;
                }
            }

            // Compute difference from previous step
            double df = sqrt( pow(CurrentPositions.at(i)[0]-PreviousPositions.at(i)[0],2) +
                    pow(CurrentPositions.at(i)[1]-PreviousPositions.at(i)[1],2) +
                    pow(CurrentPositions.at(i)[2]-PreviousPositions.at(i)[2],2) );

            if (df>deltamax) {
                deltamaxnode = i;
                deltamax=df;
            }

            // Update previous positions
            for (unsigned int j=0; j<3; j++) {
                if (!FixedDirections.at(i)[j]) {
                    PreviousPositions.at(i)[j]=CurrentPositions.at(i)[j];
                }
            }


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
                if (!FixedDirections.at(i)[j]) {
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

            if (!FixedDirections.at(i)[j]) {
                Vertices.at(i)->set_c(CurrentPositions.at(i)[j], j); // TODO: Use array to improve performance
            }
        }
    }

}

}
