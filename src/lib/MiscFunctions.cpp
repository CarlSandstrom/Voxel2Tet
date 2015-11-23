#include <vector>
#include <array>
#include <stdio.h>
#include <algorithm>
#include <stdarg.h>
#include <MeshComponents.h>

namespace voxel2tet {

void dolog(const char *functionname, const char *fmt, ...)
{
#if LOGOUTPUT > 0
        va_list argp;
        va_start(argp, fmt);
        printf("%s:\t", functionname);
        vfprintf(stdout, fmt, argp);
        va_end(argp);
#endif
}

void dooutputstat(const char *fmt, ...)
{
        va_list argp;
        va_start(argp, fmt);
        vfprintf(stdout, fmt, argp);
        va_end(argp);
}

template <typename T>
std::vector<int> FindSubsetIndices(std::vector<T> Container, std::vector<T> Subset)
{
    std::vector <int> Indices;
    for (auto s: Subset) {
        int pos = std::distance(Container.begin(), std::find(Container.begin(), Container.end(), s));
        Indices.push_back(pos);
    }
    return Indices;
}

void SpringSmooth(std::vector<VertexType*> Vertices, std::vector<std::array<bool,3>> FixedDirections, std::vector<std::vector<VertexType*>> Connections, double K)
{
    std::vector<std::array<double, 3>> CurrentPositions;
    std::vector<std::array<double, 3>> PreviousPositions;

    for (unsigned int i=0; i<Vertices.size(); i++) {
        std::array<double, 3> cp;
        std::array<double, 3> pp;
        for (int j=0; j<3; j++) {
            cp.at(j) = pp.at(j) = Vertices.at(i)->c[j];
        }
        CurrentPositions.push_back(cp);
        PreviousPositions.push_back(pp);
    }

    int itercount = 0;
    double deltamax = 1e12;
    while ((itercount < 1000) & (deltamax>1e-8)){

        deltamax=0.0;

        for (unsigned int i=0; i<Vertices.size(); i++) {
            std::array<double, 3> NewCoords = {0,0,0};

            std::vector<VertexType*> MyConnections = Connections.at(i);
            for (unsigned int k=0; k<MyConnections.size(); k++) {
                // We need the index of MyConnection such that we can use CurrentPosition (as CurrentPosition contains the updated coordinates)

                int VertexIndex = std::distance(Vertices.begin(), std::find(Vertices.begin(), Vertices.end(), MyConnections.at(k)));

                for (int m=0; m<3; m++) {
                    VertexType* cv = MyConnections.at(k);
                    NewCoords.at(m) = NewCoords.at(m) + CurrentPositions.at(VertexIndex)[m] / double(MyConnections.size());
                }
            }

            for (int j=0; j<3; j++) {
                if (!FixedDirections.at(i)[j]) {
                    CurrentPositions.at(i)[j] = NewCoords[j];
                }
            }

            // Pull back
            std::array<double, 3> delta, unitdelta;
            for (int j=0; j<3; j++) delta[j]=CurrentPositions.at(i)[j]-Vertices.at(i)->c[j];

            double d0 = sqrt( pow(delta[0], 2) + pow(delta[1], 2) + pow(delta[2], 2) );
            if ((d0>1e-8) & (K!=0.0)) {
                for (int j=0; j<3; j++) unitdelta[j] = delta[j]/d0;

                double F = d0*1.0;
                double d = d0;

                double change=1e8;
                while (change>1e-8) {
                    double NewDelta = F*1.0/exp(pow(d , 2) / K);

                    change = fabs(d-NewDelta);
                    d = NewDelta;
                }

                for (int j=0; j<3; j++) CurrentPositions.at(i)[j] = Vertices.at(i)->c[j] + unitdelta[j]*d;
            }

            // Compute difference from previous step
            double d = sqrt( pow(CurrentPositions.at(i)[0]-PreviousPositions.at(i)[0],2) +
                             pow(CurrentPositions.at(i)[1]-PreviousPositions.at(i)[1],2) +
                             pow(CurrentPositions.at(i)[2]-PreviousPositions.at(i)[2],2) );
            if (d>deltamax) deltamax=d;

            // Update previous positions
            for (unsigned int j=0; j<3; j++) {
                if (!FixedDirections.at(i)[j]) {
                    PreviousPositions.at(i)[j]=CurrentPositions.at(i)[j];
                }
            }


        }
        itercount ++;
    }

    // Update vertices
    for (unsigned int i=0; i<Vertices.size(); i++) {
        for (int j=0; j<3; j++) {
            if (!FixedDirections.at(i)[j]) {
                Vertices.at(i)->c[j] = CurrentPositions.at(i)[j];
            }
        }
    }

}




template std::vector<int> FindSubsetIndices(std::vector<VertexType*>, std::vector<VertexType*>);

}