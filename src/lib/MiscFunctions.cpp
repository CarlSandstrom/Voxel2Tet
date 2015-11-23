#include <vector>
#include <array>
#include <stdio.h>
#include <algorithm>
#include <stdarg.h>
#include <MeshComponents.h>

namespace voxel2tet {

void dolog(const char *functionname, const char *fmt, ...)
{
#if LOGOUTPUT > -10
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

void SpringSmooth(std::vector<VertexType*> Vertices, std::vector<VertexType*> FixedVertices,
                  std::vector<std::array<bool,3>> FixedDirections, std::vector<std::vector<VertexType*>> Connections, double K)
{
    std::vector<std::array<double, 3>> CurrentPositions;
    std::vector<std::array<double, 3>> PreviousPositions;

    std::vector<int> FixedVerticesIndices = FindSubsetIndices(Vertices, FixedVertices);

    for (unsigned int i=0; i<Vertices.size(); i++) {
        std::array<double, 3> cp;
        std::array<double, 3> pp;
        for (int j=0; j<3; j++) {
            cp.at(j) = pp.at(j) = Vertices.at(i)->c[j];
        }
        CurrentPositions.push_back(cp);
        PreviousPositions.push_back(pp);
    }

    int itercount =0;
    while (itercount < 10) {
        for (unsigned int i=0; i<Vertices.size(); i++) {
            std::array<double, 3> NewCoords = {0,0,0};

            std::vector<VertexType*> MyConnections = Connections.at(i);
            for (unsigned int k=0; k<MyConnections.size(); k++) {
                for (int m=0; m<3; m++) {
                    VertexType* cv = MyConnections.at(k);
                    NewCoords.at(m) = NewCoords.at(m) + cv->c[m] / double(MyConnections.size());
                }
            }

            for (int j=0; j<3; j++) {
                CurrentPositions.at(i)[j] = NewCoords[j];
            }

        }
        itercount ++;
    }

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
