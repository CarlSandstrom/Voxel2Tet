#ifndef SURFACE_H
#define SURFACE_H

#include <string>

#include "MeshComponents.h"
#include "PhaseEdge.h"
#include "MiscFunctions.h"
#include "Smoother.h"
#include "armadillo"

namespace voxel2tet
{

/**
 * @brief The Surface class holds information and means for manipulation of the surface separating
 * two materials.
 *
 * One surface object is created for each pair of connected materials. The surface is described using
 * triangles.
 */
class Surface
{
private:
    /**
     * @brief Pointer to an Options object that contains all command line arguments.
     */
    Options *Opt;
public:

    /**
     * @brief Constructor for the Surface class
     * @param Phase1 Phase index on one side of the surface
     * @param Phase2 Phase index on the other side of the surface
     * @param Opt Options object for command line arguments and other settings
     */
    Surface(int Phase1, int Phase2, Options *Opt);

    /**
     * @brief Holds phase indices for phases on either side of the surface
     */
    std::array<int, 2> Phases;

    /**
     * @brief Vertices List of vertices on the surface
     */
    std::vector <VertexType*> Vertices;

    /**
     * @brief Triangles List of triangles on the surface
     */
    std::vector <TriangleType*> Triangles;

    /**
     * @brief PhaseEdges List of phase edges on the surface
     */
    std::vector <PhaseEdge*> PhaseEdges;

    /**
     * @brief AddVertex Adds the vertex \p Vertex to the surface
     * @param Vertex Vertex to be added to the Surface object
     */
    void AddVertex(VertexType* Vertex);

    /**
     * @brief AddTriangle Adds a triangle to the surface
     * @param Triangle Triangle to be added
     */
    void AddTriangle(TriangleType* Triangle);

    /**
     * @brief Perform smoothing on this Surface object
     * @param Mesh Mesh to be smoothed
     * @param c Non-linear spring constant
     */
    void Smooth(MeshData *Mesh, double c, double alpha, double charlength, bool Automatic_c=false);

    /**
     * @brief Compute area of surface
     * @return Area of Surface object
     */
    double ComputeArea();

    /**
     * @brief Reorient triangles such that the normal is consistent on the whole surface.
     * This is done by comparing the phases on the negative and positive normal sides on
     * the triangle rather than comparing the order of the nodes.
     */
    void ReorientTriangles();

    /**
     * @brief Computes the surface integral of x.n to be used when computing the volume
     * @return Integral of x.n
     */
    double ComputeIntegral_nx();

};

}
#endif // SURFACE_H
