#ifndef MESHCOMPONENTS_H
#define MESHCOMPONENTS_H

#include <vector>
#include <array>
#define ARMA_USE_CXX11
#include <armadillo>
#include "math.h"

namespace voxel2tet
{
class TriangleType;

class EdgeType;

class PhaseEdge;

/**
 * @brief The VertexType class provides information relevant to one vertex in the mesh.
 *
 * The VertexType object knows which triangles and edges it is connected to. The TriangleType
 * and EdgeType objects however, does only know which vertices they are connected to. In order
 * for a triangle to access the connected edges, it har to ask its vertices. The same goes for
 * the edges looking for connected triangles.
 */
class VertexType
{
private:
    std::array<double, 3> c;

public:
    /**
     * @brief Specifies the value of the 'c' constant in the non-linear spring
     */
    double c_constant;

    /**
     * @brief The error accumulated due to mesh coarsening on this vertex
     */
    double error;

    /**
     * @brief ID of current vertex. Mainly for debugging purposes.
     */
    int ID;

    /**
     * @brief Arbitrary tag on vertex.
     */
    int tag;

    /**
     * @brief List of Phase edges the vertex is connected to.
     */
    std::vector<PhaseEdge *> PhaseEdges;

    /**
     * @brief Add a phase edge to PhaseEdges
     * @param pe Input. Pointer to phase edge object to add.
     */
    void AddPhaseEdge(PhaseEdge *pe);

    /**
     * @brief IsFixedVertex tells if the vertex is (always) fixed. A vertex that is the end of a phase edge is always fixed.
     * @return If the vertex is fixed or not
     */
    bool IsFixedVertex()
    {
        this->PhaseEdges.erase(std::unique(this->PhaseEdges.begin(), this->PhaseEdges.end()), this->PhaseEdges.end());
        return PhaseEdges.size() > 1;
    }

    /**
     * @brief Determines in which directions this vertex is fixed during smoothing.
     */
    std::array<bool, 3> Fixed;

    /**
     * @brief Determines if this vertex belongs to an phase edge.
     * @return True or false depending on if the vertex belongs to a phase edge.
     */
    bool IsPhaseEdgeVertex()
    { return PhaseEdges.size() > 0; }

    /**
     * @brief Constructor
     * @param x X coordinate
     * @param y Y coodrinate
     * @param z Z coordinate
     */
    VertexType(double x, double y, double z);

    /**
     * @brief Contains the original coordinates before smoothing.
     */
    double originalcoordinates[3];

    /**
     * @brief Update coordinates of this vertex
     * @param newc Array containing new coordinate information
     */
    void set_c(std::array<double, 3> newc);

    /**
     * @brief Update coordinates of this vertex
     * @param c Coordinate value
     * @param index Index of coordinate, 0 is X, 1 is Y and 2 is Z.
     */
    void set_c(double c, int index);

    /**
     * @brief Retrieve coordinates of current vertex
     * @return Coordinates
     */
    std::array<double, 3> get_c();

    /**
     * @brief Retrieve coordinate by index
     * @param index Index of coordinate. 0 is X, 1 is Y and 2 is Z.
     * @return Coordinate value
     */
    double get_c(int index);

    /**
     * Gives coordinate as arma::vec
     * @return Coordinate
     */
    arma::vec get_c_vec();

    /**
     * @brief List of triangles connected to this vertex
     */
    std::vector<TriangleType *> Triangles;

    /**
     * @brief List of edges connected to this vertex
     */
    std::vector<EdgeType *> Edges;

    /**
     * @brief Adds a triangle to vertex triangle list and ensures uniqueness on list
     * @param Triangle Pointer to triangle to add
     */
    void AddTriangle(TriangleType *Triangle);

    /**
     * @brief Removes triangle from list. Note, no memory is freed.
     * @param Triangle Pointer to triangle to remove
     */
    void RemoveTriangle(TriangleType *Triangle);

    /**
     * @brief Adds an edge to vertex edge list and ensures uniqueness on list
     * @param Edge Pointer to edge to add
     */
    void AddEdge(EdgeType *Edge);

    /**
     * @brief Removes an edge from edge list
     * @param Edge Edge to be removed. Note, no memory is freed
     */
    void RemoveEdge(EdgeType *Edge);

    /**
     * @brief Produce a list of all neighbouring vertices
     * @return List of neighbouring vertices
     */
    std::vector<VertexType *> FetchNeighbouringVertices();

};

/**
 * @brief Provides information of one edge.
 *
 * Here, an edge referres to an edge of a triangle which can be shared among several triangles. An edge is
 * defined by it's two end points.
 *
 * The Edge object knows which two vertices it is defined by. It does not, however, know which triangles are
 * connected. For that information it relies on the vertices, which in turn, knows.
 *
 */
class EdgeType
{
public:
    /**
     * True if the edge is a transverse, i.e. a diagonal of a square. Edges that are transverse should no cound as a spring during smoothing
     */
    bool IsTransverse;

    /**
     * @brief ID of edge. Mainly for debugging purposes.
     */
    int ID;

    /**
     * @brief Array of two vertices describing the edge.
     */
    std::array<VertexType *, 2> Vertices;

    /**
     * @brief Produces a list of triangles connected to this edge.
     *
     * The triangles connected to this edge is computed as
     *
     * \f$ T(V_1) \cap T(V_2) \f$
     *
     * where \f$ T(V_\alpha) \f$ is the set of all triangles connected to vertex \f$V_\alpha\f$.
     *
     * @return List of triangles
     */
    std::vector<TriangleType *> GiveTriangles();

    /**
     * @brief Computes the length of this edge.
     * @return Length of edge
     */
    double GiveLength();

    /**
     * @brief Compute center point of edge.
     * @return Center point of edge
     */
    std::array<double, 3> GiveCenterPoint();
};

/**
 * @brief The TriangleType class provides information for a triangle.
 *
 * The Triangle objects is defined by three coordinates. Although it's orientation is determined by the order in which
 * the coordinates are given, we can determine if the normal is pointing inwards or outwards of the surface volume it is
 * part of, by comparing the PosNormalMatID and NegNormalMatID members to the material ID of the volume.
 *
 * The Triangle object does not know which Edge objects it is connected to but relies on its vertices to find out.
 */
class TriangleType
{
private:
    // Give vector for edge between node and node+1
    std::array<double, 3> GiveEdgeVector(int node);

    // Normal of element
    std::array<double, 3> Normal;
public:

    /**
     * @brief Constructor
     */
    TriangleType()
    {}

    /**
     * @brief Constructor
     * @param Vertices VertexType objects determining the corners of the triangle
     */
    TriangleType(std::array<VertexType *, 3> Vertices);

    /**
     * @brief ID of Triangle object. Mostly for debugging purposes.
     */
    int ID;

    /**
     * @brief ID of the interface the triangle is part of.
     */
    int InterfaceID;

    /**
     * @brief Material ID of the material in the positive normal direction
     */
    int PosNormalMatID;

    /**
     * @brief Material ID of the material in the negative normal direction
     */
    int NegNormalMatID;

    /**
     * @brief Array of vertices defining the triangle
     */
    std::array<VertexType *, 3> Vertices;

    /**
     * @brief Returns the Edge object located at edge index
     * @param Index Index of edge to be retrieved
     * @return Pointer to Edge object
     */
    EdgeType *GiveEdge(int Index);

    /**
     * @brief Returns array of pointer to Edge objects defining the triangle
     *
     * The edges are found as
     *
     * \f$ \cap_{i=1}^3 E(V_i) \f$
     *
     * where \f$ E(V_i) \f$ are all edges connected to vertex \f$ i \f$ of the triangle.
     *
     * @return Array of pointers to Edge objects
     */
    std::array<EdgeType *, 3> GiveEdges();

    /**
     * @brief Computes the longest edge length of the triangle.
     * @return Length of longest edge
     */
    double GiveLongestEdgeLength();

    /**
     * @brief Computes the (positive) area of the triangle
     * @return Area
     */
    double GiveArea();

    /**
     * @brief Computes the (signed) area of the triangle
     * @return Area
     */
    double GiveSignedArea();

    /**
     * @brief Computes center of mass for triangle
     * @return Coordinate
     */
    std::array<double, 3> GiveCenterOfMass();

    /**
     * @brief Computes the largest inner angle of the triangle.
     * @param index Output. Pointer to integer. If set, contains the index of the corner of the largest angle.
     * @return Largest inner angle
     */
    double GiveLargestAngle(int *index = NULL);

    /**
     * @brief Computes smallest inner angle of the triangle.
     * @param index Output. Pointer to integer. If set, contains the index of the corner of the smallest angle.
     * @return Smallest inner angle.
     */
    double GiveSmallestAngle(int *index = NULL);

    /**
     * @brief Gives normal of triangle.
     *
     * For performance, the normal is not computed every time and requires updating as coordinates or orientation change.
     *
     * @return Array of doubles describing the normal
     */
    std::array<double, 3> GiveNormal()
    { return Normal; }

    /**
     * @brief Gives normalized normal. See GiveNormal().
     * @return Normalized normal
     */
    std::array<double, 3> GiveUnitNormal();

    /**
     * @brief Updates normal
     */
    void UpdateNormal();

    /**
     * @brief Change orientation of triangle by reordering the vertizes. Also updates PosNormalMatID and NegNormalMatID.
     */
    void FlipNormal();
};

/**
 * @brief The TetType class contains information on a tetrahedral element
 */
class TetType
{
private:
    int FaceVertices[4][3] = {{0, 2, 1},
                              {0, 1, 3},
                              {1, 2, 3},
                              {0, 3, 2}};
    int FacePairs[6][2] = {{0, 1},
                           {0, 2},
                           {0, 3},
                           {1, 3},
                           {1, 2},
                           {2, 3}};
    int EdgeVertices[6][2] = {{0, 1},
                              {1, 2},
                              {2, 0},
                              {0, 3},
                              {1, 3},
                              {2, 3}};

public:
    /**
     * @brief ID of element. Mostly for debugging purposes.
     */
    int ID;

    /**
     * @brief Material ID
     */
    int MaterialID;

    /**
     * @brief Vertices determininmg the tetrahedral element.
     */
    std::array<VertexType *, 4> Vertices;

    /**
     * @brief Returns center of mass
     * @return Coordinate
     */
    std::array<double, 3> GiveCenterOfMass();

    /**
     * @brief GiveSmallestAngle
     * @param angle Pointer to integer containing index of the edge neighboring the smallest dihedral angle
     * @return Smallest dihedral angle
     */
    double GiveSmallestDihedralAngle(int &angle);

    /**
     * Computes the angle between two planes.
     * @param index Index of the pair of faces as given in MeshComponent::FacePairs
     * @return Angle
     */
    double GiveDihedralAngle(int index);

    /**
     * Computes normal of face
     * @param index Index of face
     * @return Normal vector
     */
    arma::vec GiveNormalOfFace(int index);

    /**
     * Compute volume of element
     * @return Volume
     */
    double GiveVolume();

    /**
     * Fins the index of the shortest edge
     * @return Index of shortest edge
     */
    int GiveShortestEdgeIndex();

};
}

#endif // MESHCOMPONENTS_H
