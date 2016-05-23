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

void Surface::Smooth(MeshData *Mesh)
{

    double K = this->Opt->GiveDoubleValue("spring_const");

    std::vector<VertexType *> SmoothVertices;
    std::vector<std::vector<VertexType *>> Connections;
    std::vector<bool> FixedList;

    // Create connections matrix
    //for (unsigned int i=0; i<this->Vertices.size(); i++) {
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

    SpringSmooth(this->Vertices, FixedList, Connections, K, Mesh);

}

double Surface::ComputeArea()
{
    double Area = 0;
    for (TriangleType *t: this->Triangles) {
        Area=Area+t->GiveArea();
    }
    return Area;
}

}
