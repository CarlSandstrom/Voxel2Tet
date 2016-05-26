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

    for (unsigned int i=0; i<this->Triangles.size(); i++) {

        // See if any other triangles has the same edge. If so, it is oriented in another way
        TriangleType *t1 = this->Triangles[i];
        std::array<std::array<VertexType *, 2>, 3> edges1 = {{{t1->Vertices[0], t1->Vertices[1]}, {t1->Vertices[1], t1->Vertices[2]}, {t1->Vertices[2], t1->Vertices[0]}}};

        t1->UpdateNormal();
        std::array<double,3> n=t1->GiveNormal();
        LOG("n[%u] = [%f, %f, %f]\n", t1->ID, n[0], n[1], n[2]);

        for (unsigned int j=i; j<this->Triangles.size(); j++) {
            bool EdgeFound = false;
            if (i!=j) {
                TriangleType *t2 = this->Triangles[j];
                std::array<std::array<VertexType *, 2>, 3> edges2 = {{{t2->Vertices[0], t2->Vertices[1]}, {t2->Vertices[1], t2->Vertices[2]}, {t2->Vertices[2], t2->Vertices[0]}}};
                for (int k=0; k<3; k++) {
                    for (int l=0; l<3; l++) {
                        if ( (edges1[k][0]==edges2[l][0]) & (edges1[k][1]==edges2[l][1])) {
                            // Orientation is not consistent with other triangle defined earlier.

                            // Reorient
                            VertexType *v = t2->Vertices[1];
                            t2->Vertices[1] = t2->Vertices[0];
                            t2->Vertices[0] = v;

                            // Update phases on either positive and negative side
                            int temp = t2->PosNormalMatID;
                            t2->PosNormalMatID = t2->NegNormalMatID;
                            t2->NegNormalMatID = temp;
                            EdgeFound = true;
                            break;
                        }
                    }
                    if (EdgeFound) break;
                }
            }
            if (EdgeFound) break;
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
