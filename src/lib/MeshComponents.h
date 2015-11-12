#ifndef MESHCOMPONENTS_H
#define MESHCOMPONENTS_H

#include<vector>

namespace voxel2tet
{

class TriangleType;
class EdgeType;

class VertexType {
public:
    VertexType (double x, double y, double z) {this->c[0] = x; this->c[1] = y; this->c[2] = z; }
    double c[3];
    std::vector <TriangleType*> Triangles;
    std::vector <EdgeType*> Edges;
};

class EdgeType {
public:
    VertexType Vertices[2];
};

class TriangleType {
public:
    EdgeType Edges[3];
};

}

#endif // MESHCOMPONENTS_H
