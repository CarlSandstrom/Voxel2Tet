#ifndef MESHCOMPONENTS_H
#define MESHCOMPONENTS_H

#include<vector>
#include<array>

#include "math.h"

namespace voxel2tet
{

class TriangleType;
class EdgeType;
class PhaseEdge;

class VertexType {
private:
    std::array<double,3> c;

public:
    /**
     * @brief Specifies the value of the 'c' constant in the non-linear spring
     */
    double c_constant;

    int ID;

    int tag;

    std::vector<PhaseEdge*> PhaseEdges;

    void AddPhaseEdge(PhaseEdge*);

    /**
     * @brief IsFixedVertex tells if the vertex is (always) fixed. A vertex that is the end of a phase edge is always fixed.
     * @return If the vertex is fixed or not
     */
    bool IsFixedVertex() {return PhaseEdges.size()>1;}

    std::array<bool, 3> Fixed;
    bool IsPhaseEdgeVertex() {return PhaseEdges.size()>0;}

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

    double GiveLength();
    std::array<double, 3> GiveCenterPoint();
};


class TriangleType {
private:
    // Give vector for edge between node and node+1
    std::array<double, 3> GiveEdgeVector(int node);

    // Normal of element
    std::array<double, 3> Normal;
public:
    TriangleType () {}
    TriangleType (std::array<VertexType *, 3> Vertices);
    int ID;
    int InterfaceID;
    int PosNormalMatID;
    int NegNormalMatID;

    std::array<VertexType *, 3> Vertices;

    EdgeType * GiveEdge(int Index);
    std::array<EdgeType *, 3> GiveEdges();

    // Compute area of triangle
    double GiveArea();

    // Compute signed area of triangle
    double GiveSignedArea();

    std::array<double, 3> GiveCenterOfMass();

    double GiveLargestAngle(int *index=NULL);
    double GiveSmallestAngle(int *index=NULL);

    std::array<double, 3> GiveNormal() {return Normal;}
    std::array<double, 3> GiveUnitNormal();

    void UpdateNormal();

    // Reverse the order of the vertices to flip the normal
    void FlipNormal();

};


class TetType {
private:
public:
    int ID;
    int MaterialID;
    std::array<VertexType *, 4> Vertices;
    std::array<double, 3> GiveCenterOfMass();
};

}

#endif // MESHCOMPONENTS_H
