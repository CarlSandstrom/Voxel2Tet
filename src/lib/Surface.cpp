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

void Surface :: MoveAsTrussStructure()
{
    // Find all nodes that has moved. These will constitute the set of fixed nodes
    std::vector<VertexType*> FixedNodes;

    for (auto v: this->Vertices) {
        for (int i=0; i<3; i++) {
            if (v->c[i]!=v->originalcoordinates[i]) {
                FixedNodes.push_back(v);
                break;
            }
        }
    }

}


void Surface::Smooth()
{

    double K = this->Opt->GiveDoubleValue("spring_const");

    std::vector<std::vector<VertexType *>> Connections;
    std::vector<std::array<bool,3>> FixedDirectionsList;

    // Create connections matrix
    for (unsigned int i=0; i<this->Vertices.size(); i++) {
        // Find connected vertices
        std::vector <VertexType*> NeighbouringVertices = this->Vertices.at(i)->FetchNeighbouringVertices();
        std::sort (NeighbouringVertices.begin(), NeighbouringVertices.end());
        std::vector <VertexType*> ConnectedVertices;

        // Create list of indices of connected vertices
        std::set_intersection(NeighbouringVertices.begin(), NeighbouringVertices.end(),
                              this->Vertices.begin(), this->Vertices.end(), back_inserter(ConnectedVertices));

        Connections.push_back(ConnectedVertices);

        std::array<bool,3> FixedDirections;
        if (std::find(this->FixedVertices.begin(), this->FixedVertices.end(), this->Vertices.at(i))==this->FixedVertices.end()  ) {
            FixedDirections = {false, false, false};
        } else {
            FixedDirections = {true, true, true};
        }

        FixedDirectionsList.push_back(FixedDirections);

    }

    SpringSmooth(this->Vertices, FixedDirectionsList, Connections, K);

}

}
