#include <vector>
#include <array>
#include <stdio.h>
#include <algorithm>
#include <stdarg.h>
#include <MeshComponents.h>
#include <sstream>

#include "MeshData.h"

namespace voxel2tet {

void dolog(const char *functionname, const char *fmt, ...)
{
#if LOGOUTPUT == 1
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

void dooutputlogmesh(MeshData &Mesh, char *filename, ...)
{
    va_list argp;
    va_start(argp, filename);
    char buffer[255];
    vsprintf(buffer, filename, argp);
    Mesh.ExportSurface(buffer, FT_VTK);
    va_end(argp);
}

template <typename T>
std::vector<int> FindSubsetIndices(std::vector<T> Container, std::vector<T> Subset)
{
    std::vector <int> Indices;
    for (auto s: Subset) {
        unsigned int pos = std::distance(Container.begin(), std::find(Container.begin(), Container.end(), s));
        if (pos<Container.size()) {
            Indices.push_back(pos);
        }
    }
    return Indices;
}

std::array<double,3> ComputeNormalizedVector(VertexType* v1, VertexType* v2)
{
    std::array<double, 3> Normal;
    for (int i=0; i<3; i++) Normal.at(i) = v1->get_c(i)-v2->get_c(i);
    double NewNormalLength = sqrt(Normal.at(0)*Normal.at(0) + Normal.at(1)*Normal.at(1) +Normal.at(2)*Normal.at(2));
    for (int i=0; i<3; i++) Normal.at(i) = Normal.at(i) / NewNormalLength;
    return Normal;
}



template std::vector<int> FindSubsetIndices(std::vector<VertexType*>, std::vector<VertexType*>);

}
