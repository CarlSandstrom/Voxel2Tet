#include "Surface.h"

namespace voxel2tet
{
Surface :: Surface(int Phase1, int Phase2, Options *Opt, Smoother *Smoother)
{
    this->Phases [ 0 ] = Phase1;
    this->Phases [ 1 ] = Phase2;
    this->Opt = Opt;
    this->SurfaceSmooth = Smoother;
}

void Surface :: AddVertex(VertexType *Vertex)
{
    for ( auto v : this->Vertices ) {
        if ( v == Vertex ) {
            return;
        }
    }
    this->Vertices.push_back(Vertex);
}

void Surface :: AddTriangle(TriangleType *Triangle)
{
    for ( auto t : this->Triangles ) {
        if ( t == Triangle ) {
            return;
        }
    }
    this->Triangles.push_back(Triangle);
}

void Surface :: Smooth(MeshData *Mesh)
{
    // Collect all of the surface vertices and then remove all vertices that belongs tho an EdgePhase
    std::vector<VertexType *> VerticesToSmooth = this->Vertices;

    size_t i=0;
    while(i<VerticesToSmooth.size()) {
        if (VerticesToSmooth[i]->PhaseEdges.size()>0) {
            VerticesToSmooth.erase(std::remove(VerticesToSmooth.begin(), VerticesToSmooth.end(), VerticesToSmooth[i]), VerticesToSmooth.end());
        } else {
            i++;
        }
    }
    this->SurfaceSmooth->Smooth(VerticesToSmooth, Mesh);
}

double Surface :: ComputeArea()
{
    double Area = 0;
    for ( TriangleType *t : this->Triangles ) {
        Area = Area + t->GiveArea();
    }
    return Area;
}

void Surface :: ReorientTriangles()
{
    // Reorient all triangles using first triangle in list as reference
    int PosRef = this->Triangles [ 0 ]->PosNormalMatID;
    for ( unsigned int i = 1; i < this->Triangles.size(); i++ ) {
        TriangleType *t = this->Triangles [ i ];
        if ( t->PosNormalMatID != PosRef ) {
            t->FlipNormal();
        }
    }
}

double Surface :: ComputeIntegral_nx() // I think this does not currently work.
{
    double result = 0;
    for ( TriangleType *t : this->Triangles ) {
        double da = t->GiveArea();
        std :: array< double, 3 >cm = t->GiveCenterOfMass();
        std :: array< double, 3 >n = t->GiveUnitNormal();
        result = result + da * ( cm [ 0 ] * n [ 0 ] + cm [ 1 ] * n [ 1 ] + cm [ 2 ] * n [ 2 ] );
    }

    return result;
}
}
