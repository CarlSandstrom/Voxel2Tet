#ifndef MESHCOMPONENTS_H
#define MESHCOMPONENTS_H

#include<vector>
#include<array>

#include "math.h"

namespace voxel2tet
{

class TriangleType;
class EdgeType;

class VertexType {
private:

public:
    int ID;

    VertexType (double x, double y, double z);
    double c[3];
    double originalcoordinates[3];
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

    // Finds all neighbouring vertices
    std::vector <VertexType*> FetchNeighbouringVertices();

};


class EdgeType {
public:
    std::array<VertexType *, 2> Vertices;

    std::vector<TriangleType *> GiveTriangles();
};


class TriangleType {
private:
    std::array<double, 3> GiveVector(int node);
public:
    int InterfaceID;
    int PosNormalMatID;
    int NegNormalMatID;
    std::array<VertexType *, 3> Vertices;

    std::array<EdgeType *, 3> GiveEdges();

    // Compute area of triangle
    double GiveArea();

    double GiveLargestAngle(int *index=NULL);
};

}

#endif // MESHCOMPONENTS_H
