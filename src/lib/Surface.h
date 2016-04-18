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

class Surface
{
private:
    Options *Opt;
public:
    std::vector <VertexType*> FixedVertices;

    Surface(int Phase1, int Phase2, Options *Opt);
    int Phases[2];
    std::vector <VertexType*> Vertices;
    std::vector <TriangleType*> Triangles;
    std::vector <PhaseEdge*> PhaseEdges;

    /* Moves all internal nodes as if it were a truss structure. This is done in order
     * to prevent triangles from intersecting during the smoothing phase.
     */
    void MoveAsTrussStructure();

    void AddVertex(VertexType* Vertex);
    void AddTriangle(TriangleType* Triangle);

    void Smooth();
};

}
#endif // SURFACE_H
