#ifndef MESHCOMPONENTS_H
#define MESHCOMPONENTS_H

#include<vector>

namespace voxel2tet
{

class TriangleType;
class EdgeType;

class VertexType {
private:

public:
    VertexType (double x, double y, double z) {this->c[0] = x; this->c[1] = y; this->c[2] = z; }
    double c[3];
    std::vector <TriangleType*> Triangles;
    std::vector <EdgeType*> Edges;

    // Adds a triangle to vertex triangle list and ensures uniqueness on list
    void AddTriangle(TriangleType *Triangle);

    // Removes a triangle from triangle list
    void RemoveTriangle(TriangleType *Triangle);

    // Adds an edge to vertex edge list and ensures uniqueness on list
    void AddEdge(EdgeType *Edge);

    // Removes an edge from edge list
    void RemoveEdge(EdgeType *Edge);
};

class EdgeType {
public:
    VertexType *Vertices[2];
};

class TriangleType {
public:
    int InterfaceID;
    int PosNormalMatID;
    int NegNormalMatID;
    VertexType *Vertices[3];
};

}

#endif // MESHCOMPONENTS_H
