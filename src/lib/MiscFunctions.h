#ifndef MISCFUNCTIONS_H_
#define MISCFUNCTIONS_H_

#define LOGOUTPUT 1
#define STATOUTPUT 1

#define LOG(format, args...) dolog (__FUNCTION__, format, args)
#define STATUS(format, args...) dooutputstat (format, args)

#define EPS 1.0e-8

#include<vector>

namespace voxel2tet
{

void dolog(const char *functionname, const char *fmt, ...);
void dooutputstat(const char *fmt, ...);

template <typename T> std::vector<int> FindSubsetIndices(std::vector<T> Container, std::vector<T> Subset);

}

#endif /* MISCFUNCTIONS_H_ */
