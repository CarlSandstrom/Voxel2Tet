#ifndef MISCFUNCTIONS_H_
#define MISCFUNCTIONS_H_

#define LOGOUTPUT 1
#define STATOUTPUT 1

#define LOG(format, args...) dolog (__FUNCTION__, format, args)
#define STATUS(format, args...) dooutputstat (format, args)

#define EPS 1.0e-8

#include <vector>
#include <array>
#include "MeshComponents.h"

namespace voxel2tet
{

void dolog(const char *functionname, const char *fmt, ...);
void dooutputstat(const char *fmt, ...);

void SpringSmooth(std::vector<VertexType*> Vertices, std::vector<VertexType*> FixedVertices,
                  std::vector<std::array<bool,3>> FixedDirections, std::vector<std::vector<VertexType*>> Connections, double K);

template <typename T> std::vector<int> FindSubsetIndices(std::vector<T> Container, std::vector<T> Subset);

}

#endif /* MISCFUNCTIONS_H_ */
