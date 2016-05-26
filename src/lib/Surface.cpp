#include "Surface.h"

namespace voxel2tet
{

Surface::Surface(int Phase1, int Phase2, Options *Opt)
{
    this->Phases[0] = Phase1;
    this->Phases[1] = Phase2;
    this->Opt = Opt;
}

void Surface::AddVertex(VertexType *Vertex)
{
    for (auto v: this->Vertices) {
        if (v==Vertex) {
            return;
        }
    }
    this->Vertices.push_back(Vertex);
}

void Surface::AddTriangle(TriangleType *Triangle)
{
    for (auto t: this->Triangles) {
        if (t==Triangle) {
            return;
        }
    }
    this->Triangles.push_back(Triangle);
}

void Surface::Smooth(MeshData *Mesh, double c, double alpha, double charlength, bool Automatic_c)
{

    std::vector<VertexType *> SmoothVertices;
    std::vector<std::vector<VertexType *>> Connections;
    std::vector<bool> FixedList;

    // Create connections matrix
    for (VertexType *v: this->Vertices) {
        std::vector <VertexType*> NeighbouringVertices = v->FetchNeighbouringVertices();
        std::sort (NeighbouringVertices.begin(), NeighbouringVertices.end());
        std::vector <VertexType*> ConnectedVertices;

        // Create list of indices of connected vertices
        std::set_intersection(NeighbouringVertices.begin(), NeighbouringVertices.end(),
                              this->Vertices.begin(), this->Vertices.end(), back_inserter(ConnectedVertices));

        Connections.push_back(ConnectedVertices);
        SmoothVertices.push_back(v);

        // Lock phase edges since they are already smoothed
        if (v->IsPhaseEdgeVertex()) {
            FixedList.push_back(true);
        } else {
            FixedList.push_back(false);
        }

    }

    SpringSmooth(this->Vertices, FixedList, Connections, c, alpha, charlength, Automatic_c, Mesh);

}

double Surface::ComputeArea()
{
    double Area = 0;
    for (TriangleType *t: this->Triangles) {
        Area=Area+t->GiveArea();
    }
    return Area;
}

void Surface::ReorientTriangles()
{

    // Use first triangle in list as reference
    int PosRef = this->Triangles.at(0)->PosNormalMatID;

    for (unsigned int i=1; i<this->Triangles.size(); i++) {
        TriangleType *t=this->Triangles[i];
        if (t->PosNormalMatID!=PosRef) {
            std::array<double, 3> cm = t->GiveCenterOfMass();
            VertexType *v = t->Vertices[0];
            t->Vertices[0] = t->Vertices[1];
            t->Vertices[1] = v;

            t->PosNormalMatID = t->NegNormalMatID;
            t->NegNormalMatID = PosRef;
        }
    }
}

double Surface::ComputeIntegral_nx()
{
    double result = 0;
    for (TriangleType *t: this->Triangles) {
        double da = t->GiveArea();
        std::array<double, 3> cm = t->GiveCenterOfMass();
        std::array<double, 3> n = t->GiveNormal();
        result = result + da*(cm[0]*n[0] + cm[1]*n[1] + cm[2]*n[2]);
    }

    return result;
}

}
