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
    double *uF = new double[sizeF];

    // Setup prescribed vector
    int sizeC=FixedNodes.size()*3;
    double *uC = new double[sizeC];
/*
    // Set u to the difference between current and original position
    for (unsigned int i=0; i<FixedNodes.size(); i++) {
        printf("%u, ", i);
        for (int j=0; j<3; j++) {
            double delta = FixedNodes.at(i)->c[j]-FixedNodes.at(i)->originalcoordinates[j];
            uC[i*3+j] = delta;
            printf("VertexID %u, pos %u: %f\n", i, j, delta);
        }
    }

    printf("uC=[");
    for (int i=0; i<sizeC; i++) {
        printf("%f, ", uC[i]);
    }
    printf("]\n");
*/
    // Setup stiffness matrices. Store all values contineously.
    double k=1;
    double *Kff;
    Kff = new double [sizeF*sizeF];

    for (int i=0; i<sizeF*sizeF; i++) {
        Kff[i]=0.0;
    }

    double *Kfc;
    Kfc = new double [sizeF*sizeC];

    for (int i=0; i<sizeF*sizeC; i++) {
        Kfc[i]=0.0;
    }

    for (unsigned int i=0; i<FreeNodes.size(); i++) {
        VertexType *v=FreeNodes.at(i);

        for (int j=0; j<3; j++) {
            //Kff[(i*3+j)*sizeF+i*3+j]=Kff[(i*3+j)*sizeF+i*3+j] + k;
        }

        // Fetch connected vertices
        std::vector<VertexType*> Connections = v->FetchNeighbouringVertices();

        // Setup Kff
        std::vector<int> FreeConnectionIndices = FindSubsetIndices(FreeNodes, Connections);

        for (int freeid: FreeConnectionIndices) {
            for (int j=0; j<3; j++) {
                Kff[(i*3+j) + (freeid*3 + j)*sizeF] = Kff[(i*3+j) + (freeid*3 + j)*sizeF] - k;
                Kff[(i*3+j)*sizeF+i*3+j]=Kff[(i*3+j)*sizeF+i*3+j] + k;
            }
        }

        // Setup Kfc
        std::vector<int> FixedConnectionIndices = FindSubsetIndices(FixedNodes, Connections);
        for (int FixedNode: FixedConnectionIndices) {
            for (int j=0; j<3; j++) {
                Kfc[(i*3+j)+(FixedNode*3+j)*sizeF]=Kfc[(i*3+j)+(FixedNode*3+j)*sizeF]-k;
                Kff[(i*3+j)*sizeF+i*3+j]=Kff[(i*3+j)*sizeF+i*3+j] + k;
            }
        }

        FreeConnectionIndices.clear();
        FixedConnectionIndices.clear();
    }
    /*

    // Log Kff
    printf("Kff=[");
    for (int i = 0; i< sizeF;i++) {
        for (int j = 0; j< sizeF;j++) {
            printf("%f, ", Kff[i*sizeF+j]);
        }
        printf("\n");
    }
    printf("]\n");

    // Log Kfc
    printf("Kfc=[");
    for (int i = 0; i< sizeF;i++) {
        for (int j = 0; j< sizeC;j++) {
            printf("%f, ", Kfc[i*sizeC+j]);
        }
        printf("\n");
    }
    printf("]\n");
*/

    // Compute RHS
    double *KFCuC = new double[sizeF];
    char *transA = "N";
    char *transB = "N";
    double onedbl = -1.0;
    int oneint = 1;
    double zerodbl = 0.0;
    dgemm_(transA, transB, &sizeF, &oneint, &sizeC, &onedbl, Kfc, &sizeF, uC, &sizeC, &zerodbl, KFCuC, &sizeF);

    // Solve linear system of equations

    int *piv = new int[sizeF];
    int info;

    dgesv_(&sizeF, &oneint, Kff, &sizeF, piv, KFCuC, &sizeF, &info);

    for (unsigned int i=0; i< FreeNodes.size(); i++) {
        VertexType *v = FreeNodes.at(i);
        for (int j=0; j<3; j++) {
            v->c[j] = v->c[j] + KFCuC[3*i+j];
        }
    }

    /*
    printf("KFCuC=[");
    for (int i=0; i<sizeF; i++) {
        printf("%f, ", KFCuC[i]);
    }
    printf("]\n");

    printf("nodes = np.array([[");
    for (auto v: this->Vertices) {
        for (int j=0; j<3; j++) {
            printf("%f, ", v->originalcoordinates[j]);
        }
        printf("], [");
    }
    printf("]]\n");

    printf("topology=np.array([[");
    for (auto v: this->Vertices) {
        std::vector<VertexType*> Connections = v->FetchNeighbouringVertices();
        std::vector<int> indices = FindSubsetIndices(this->Vertices, Connections);
        for (auto i: indices) {
            printf("%i, ", i);
        }
        printf("],[");
    }
    printf("]]\n");

    printf("dirichlet=[[");
    std::vector<int> dnodes = FindSubsetIndices(this->Vertices, FixedNodes);
    for (unsigned int i=0; i<FixedNodes.size(); i++) {
        printf("%u, ", dnodes[i]);
        for (int j=0; j<3; j++) {
            printf("%f, ", FixedNodes.at(i)->c[j]-FixedNodes.at(i)->originalcoordinates[j] );
        }
        printf("],[");
    }
    printf("]]\n");

    */
    delete Kff;

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
