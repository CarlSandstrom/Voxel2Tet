#ifndef MISCFUNCTIONS_H_
#define MISCFUNCTIONS_H_

#define LOGOUTPUT 0
#define STATOUTPUT 1
#define SANITYCHECK 0

#define TEST_MESH_FOR_EACH_SMOOTHING 0
#define EXPORT_SMOOTHING_ANIMATION 0

#define TEST_MESH_FOR_EACH_COARSENING_ITERATION 0
#define EXPORT_MESH_COARSENING 0

#define TEST_MESH_BETWEEN_STEPS_TETGEN 0

#define SMOOTH_EDGES_INDIVIDUALLY 0 // Is this used?

#ifdef _MSC_VER 
#define strncasecmp _strnicmp
#define strcasecmp _stricmp
#endif

#ifdef _MSC_VER
#define LOG(format, args, ...) dolog(__FUNCTION__, format, args)
#define STATUS(format, args, ...) dooutputstat(format, args)
#else
#define LOG(format, args ...) dolog(__FUNCTION__, format, args)
#define STATUS(format, args ...) dooutputstat(format, args)
#endif

#define EPS 1.0e-8

#include <vector>
#include <array>
#include <string>
#include <memory>

#include "MeshComponents.h"


namespace voxel2tet
{
#define RADIANS(DEG) DEG * 2 * 3.141593 / 360
#define DEGREES(RAD) RAD * 360 / 2 / 3.141593

class MeshData;

void dolog(const char *functionname, const char *fmt, ...);

void dooutputstat(const char *fmt, ...);

void dooutputlogmesh(MeshData &Mesh, char *filename, ...);

std::vector<std::string> SplitString(std::string Text, char Delimiter);

template<typename ... Args>
std::string strfmt(const std::string &format, Args ... args)
{
    size_t size = snprintf(nullptr, 0, format.c_str(), args ...) + 1;   // Extra space for '\0'
    std::unique_ptr<char[]> buf(new char[size]);
    snprintf(buf.get(), size, format.c_str(), args ...);
    return std::string(buf.get(), buf.get() + size - 1); // We don't want the '\0' inside
}

std::array<double, 3> ComputeNormalizedVector(VertexType *v1, VertexType *v2);

double ComputeAngleBetweenVectors(std::array<double, 3> v1, std::array<double, 3> v2);

bool CompareTriangleAngles(TriangleType *t1, TriangleType *t2);

bool CompareEdgeLength(EdgeType *e1, EdgeType *e2);

// Operator overloading in order to perform sorting of Vertices and triangles according to ID
template<typename T>
bool SortByID(T obj1, T obj2);

template<typename T>
bool EqualByID(T obj1, T obj2);

template<typename T>
std::vector<int> FindSubsetIndices(std::vector<T> Container, std::vector<T> Subset);

}

#endif /* MISCFUNCTIONS_H_ */
