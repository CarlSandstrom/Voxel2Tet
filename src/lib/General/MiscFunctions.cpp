#include <vector>
#include <array>
#include <stdio.h>
#include <algorithm>
#include <stdarg.h>
#include <MeshComponents.h>
#include <sstream>
#ifdef OPENMP
#include <omp.h>
#endif

#include "MeshData.h"

namespace voxel2tet {

void dolog(const char *functionname, const char *fmt, ...)
{
#if LOGOUTPUT == 1
    va_list argp;
    va_start(argp, fmt);
#ifdef OPENMP
    printf("%s[%u]:\t", functionname, omp_get_thread_num());
#else
    printf("%s:\t", functionname);
#endif
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

std::vector<std::string> SplitString (std::string Text, char Delimiter)
{
    int i=0;
    std::vector<std::string> StringList;
    std::string CurrentItem;

    while (i<Text.length()) {
        if (Text.at(i)!=Delimiter) {
            CurrentItem += Text.at(i);
        } else {
            StringList.push_back(CurrentItem);
            CurrentItem.clear();
        }
        i++;
    }

    if (CurrentItem.length()>0) {
        StringList.push_back(CurrentItem);
        CurrentItem.clear();
    }

    return StringList;
}

std::array<double,3> ComputeNormalizedVector(VertexType* v1, VertexType* v2)
{
    std::array<double, 3> Normal;
    for (int i=0; i<3; i++) Normal.at(i) = v1->get_c(i)-v2->get_c(i);
    double NewNormalLength = sqrt(Normal.at(0)*Normal.at(0) + Normal.at(1)*Normal.at(1) +Normal.at(2)*Normal.at(2));
    for (int i=0; i<3; i++) Normal.at(i) = Normal.at(i) / NewNormalLength;
    return Normal;
}

double ComputeAngleBetweenVectors(std::array<double, 3> v1, std::array<double, 3> v2)
{
    return std::acos( v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
}


template std::vector<int> FindSubsetIndices(std::vector<VertexType*>, std::vector<VertexType*>);

}
