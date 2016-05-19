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
    std::vector<std::array<bool,3>> FixedDirectionsList;

    // Create connections matrix
    for (unsigned int i=0; i<this->Vertices.size(); i++) {
        if (this->Vertices.at(i)->PhaseEdges.size()==0) {
            std::vector <VertexType*> NeighbouringVertices = this->Vertices.at(i)->FetchNeighbouringVertices();
            std::sort (NeighbouringVertices.begin(), NeighbouringVertices.end());
            std::vector <VertexType*> ConnectedVertices;

            // Create list of indices of connected vertices
            std::set_intersection(NeighbouringVertices.begin(), NeighbouringVertices.end(),
                                  this->Vertices.begin(), this->Vertices.end(), back_inserter(ConnectedVertices));

            Connections.push_back(ConnectedVertices);
            SmoothVertices.push_back(this->Vertices.at(i));
            FixedDirectionsList.push_back(this->Vertices.at(i)->Fixed);
        }
    }

    SpringSmooth(SmoothVertices, FixedDirectionsList, Connections, K, Mesh);

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
