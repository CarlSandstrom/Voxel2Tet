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

void Surface::Smooth()
{

    double K = this->Opt->GiveDoubleValue("spring_const");

    std::vector<std::array<double, 3>> CurrentPositions;
    std::vector<std::array<double, 3>> PreviousPositions;

    std::vector<int> FixedVerticesIndices = FindSubsetIndices(this->Vertices, this->FixedVertices);

    for (unsigned int i=0; i<this->Vertices.size(); i++) {
        LOG("C = (%f, %f, %f)\n", this->Vertices.at(i)->c[0], this->Vertices.at(i)->c[1], this->Vertices.at(i)->c[2]);
        std::array<double, 3> cp;
        std::array<double, 3> pp;
        for (int j=0; j<3; j++) {
            cp.at(j) = pp.at(j) = this->Vertices.at(i)->c[j];
        }
        CurrentPositions.push_back(cp);
        PreviousPositions.push_back(pp);
    }

    int itercount=0;
    while (itercount < 100) {

        for (unsigned int i=0; i<this->Vertices.size(); i++) {
            // Move only if the vertex is not in the FixedVerticesIndecis list
            if (std::find(FixedVerticesIndices.begin(), FixedVerticesIndices.end(), i)==FixedVerticesIndices.end()) {
                // Find connected vertices
                std::vector <VertexType*> NeighbouringVertices = this->Vertices.at(i)->FetchNeighbouringVertices();
                std::sort (NeighbouringVertices.begin(), NeighbouringVertices.end());
                std::vector <VertexType*> ConnectedVertices;

                // Create list of indices of connected vertices
                std::set_intersection(NeighbouringVertices.begin(), NeighbouringVertices.end(),
                                      this->Vertices.begin(), this->Vertices.end(), back_inserter(ConnectedVertices));

                std::vector<int> ConnectedVerticesIndices = FindSubsetIndices(this->Vertices, ConnectedVertices);

                std::array<double, 3> NewCoords = {0,0,0};

                for (unsigned int j=0; j<ConnectedVertices.size(); j++) {
                    for (int k=0; k<3; k++) {
                        int index = ConnectedVerticesIndices.at(j);
                        NewCoords.at(k) = NewCoords.at(k) + PreviousPositions.at(index)[k] / double(ConnectedVertices.size());
                    }
                }

                for (int j=0; j<3; j++) {
                    CurrentPositions.at(i)[j] = NewCoords[j];
                }

                // Pull back
                std::array<double, 3> delta, unitdelta;
                for (int j=0; j<3; j++) delta[j]=CurrentPositions.at(i)[j]-this->Vertices.at(i)->c[j];

                double d0 = sqrt( pow(delta[0], 2) + pow(delta[1], 2) + pow(delta[2], 2) );
                for (int j=0; j<3; j++) unitdelta[j] = delta[j]/d0;

                double F = d0*1.0;
                double d = d0;

                double change=1e8;
                while (change>1e-8) {
                    double NewDelta = F*1.0/exp(pow(d , 2) / K);
                    change = fabs(d-NewDelta);
                    d = NewDelta;
                }

                for (int j=0; j<3; j++) CurrentPositions.at(i)[j] = this->Vertices.at(i)->c[j] + unitdelta[j]*d;

            }
        }

        // Update previous positions
        for (unsigned int j=0; j<PreviousPositions.size(); j++) {
            PreviousPositions.at(j)=CurrentPositions.at(j);
        }
        itercount ++;

    }

    //Update vertices
    for (unsigned int i=0; i<this->Vertices.size(); i++) {
        for (int j=0; j<3; j++) {
            this->Vertices.at(i)->c[j] =CurrentPositions.at(i)[j];
        }
    }
}

}
