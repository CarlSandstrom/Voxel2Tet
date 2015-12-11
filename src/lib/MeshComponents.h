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
    std::array<double,3> c;

public:
    int ID;

    VertexType (double x, double y, double z);
    double originalcoordinates[3];

    void set_c(std::array<double,3> newc);
    void set_c(double c, int index);
    std::array<double, 3> get_c();
    double get_c(int index);


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
    // Give vector for edge between node and node+1
    std::array<double, 3> GiveEdgeVector(int node);

    // Normal of element
    std::array<double, 3> Normal;
public:
    int InterfaceID;
    int PosNormalMatID;
    int NegNormalMatID;

    std::array<VertexType *, 3> Vertices;

    std::array<EdgeType *, 3> GiveEdges();

    // Compute area of triangle
    double GiveArea();

    double GiveLargestAngle(int *index=NULL);

    std::array<double, 3> GiveNormal() {return Normal;}
    void UpdateNormal();

};

}

#endif // MESHCOMPONENTS_H
