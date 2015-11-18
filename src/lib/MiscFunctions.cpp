#include <vector>
#include <stdio.h>
#include <algorithm>
#include <stdarg.h>

namespace voxel2tet {

void dolog(const char *functionname, const char *fmt, ...)
{
#if LOGOUTPUT > -1
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

/*template <typename T>
std::vector<int> FindSubsetIndices(std::vector<T> Container, std::vector<T> Subset)
{
    std::vector <int> Indices;
    for (auto s: Subset) {
        int pos = std::distance(Container.begin(), std::find(Container.begin(), Container.end(), s));
        Indices.push_back(pos);
    }
    return Indices;
}


template typename T<int>; */


}
