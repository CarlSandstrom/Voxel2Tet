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
    // Find all nodes on PhaseEdges. These will constitute the set of fixed nodes
    std::vector<VertexType*> FixedNodes;

    for (PhaseEdge* p: this->PhaseEdges) {
        std::vector<VertexType*> PhaseEdgeVertices = p->GetFlatListOfVertices();
        FixedNodes.insert(FixedNodes.end(), PhaseEdgeVertices.begin(), PhaseEdgeVertices.end());
    }

    std::sort(FixedNodes.begin(), FixedNodes.end());

    // The complement to fixed nodes is all internal (free) nodes
    std::sort(this->Vertices.begin(), this->Vertices.end()); // TODO: We can improve performance by sorting somewhere else and assume that this.Vertices is sorted
    std::vector<VertexType*> FreeNodes;
    std::set_difference(this->Vertices.begin(), this->Vertices.end(), FixedNodes.begin(), FixedNodes.end(), std::inserter(FreeNodes, FreeNodes.begin()));

    if (FreeNodes.size()==0) return;

    // Setup solution vector
    int sizeF=FreeNodes.size()*3;
    arma::mat uF(sizeF, 1);

    // Setup prescribed vector
    int sizeC=FixedNodes.size()*3;
    arma::mat uC(sizeC, 1);

    // Set u to the difference between current and original position on all fixed nodes
    for (unsigned int i=0; i<FixedNodes.size(); i++) {
        for (int j=0; j<3; j++) {
            double delta = FixedNodes.at(i)->c[j]-FixedNodes.at(i)->originalcoordinates[j];
            uC.at(i*3+j) = -delta;
        }
    }
    //std::cout <<"uC=["<< uC<<"]\n";

    // Setup stiffness matrices. Store all values contineously.
    double k=1;
    arma::mat Kff(sizeF, sizeF, arma::fill::zeros);
    arma::mat Kfc(sizeF, sizeC, arma::fill::zeros);


    for (unsigned int i=0; i<FreeNodes.size(); i++) {
        VertexType *v=FreeNodes.at(i);

        // Fetch connected vertices
        std::vector<VertexType*> Connections = v->FetchNeighbouringVertices();

        // Setup Kff
        std::vector<int> FreeConnectionIndices = FindSubsetIndices(FreeNodes, Connections);

        for (int freeid: FreeConnectionIndices) {
            for (int j=0; j<3; j++) {
                Kff.at(i*3+j, freeid*3+j) = Kff.at(i*3+j, freeid*3+j) - k;
                Kff.at(i*3+j, i*3+j) = Kff.at(i*3+j, i*3+j) + k;
            }
        }

        // Setup Kfc
        std::vector<int> FixedConnectionIndices = FindSubsetIndices(FixedNodes, Connections);
        for (int FixedNode: FixedConnectionIndices) {
            for (int j=0; j<3; j++) {
                Kfc.at(i*3+j, FixedNode*3+j) = Kfc.at(i*3+j, FixedNode*3+j) - k;
                Kff.at(i*3+j, i*3+j) = Kff.at(i*3+j, i*3+j) + k;
            }
        }

        FreeConnectionIndices.clear();
        FixedConnectionIndices.clear();
    }


    // Compute RHS
    arma::mat KFCuC(sizeF, 1);
    KFCuC = Kfc*uC;
    uF = arma::solve(Kff, KFCuC);

    for (unsigned int i=0; i< FreeNodes.size(); i++) {
        VertexType *v = FreeNodes.at(i);
        for (int j=0; j<3; j++) {
            v->c[j] = v->c[j] + uF.at(3*i+j);
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
