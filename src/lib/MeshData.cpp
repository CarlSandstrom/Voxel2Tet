#include <vector>
#include "MeshData.h"

namespace voxel2tet
{


MeshData::MeshData(BoundingBoxType BoundingBox)
{
    this->BoundingBox = BoundingBox;
    this->VertexOctreeRoot = new VertexOctreeNode(this->BoundingBox, &this->Vertices, 0);
}

void MeshData :: ExportVTK(std::string FileName)
{

}

EdgeType *MeshData :: AddEdge(std::vector<int> VertexIDs)
{
    printf("Create edge from Vertex IDs %u and %u\n", VertexIDs.at(0), VertexIDs.at(1));

    // Check if edge already exists
    int ThisVertexID = VertexIDs.at(0);
    int OtherVertexID = VertexIDs.at(1);
    VertexType *ThisVertex = this->Vertices.at(ThisVertexID);
    VertexType *OtherVertex = this->Vertices.at(OtherVertexID);

    std::vector <EdgeType*> Edges = this->Vertices.at(ThisVertexID)->Edges;
    for (auto Edge: Edges) {
        if ( ( (Edge->Vertices[0] == ThisVertex) & (Edge->Vertices[1] == OtherVertex) ) | ( (Edge->Vertices[1] == ThisVertex) & (Edge->Vertices[0] == OtherVertex) )) {
            return Edge;
        }
    }

    EdgeType *NewEdge = new EdgeType;
    NewEdge->Vertices[0] = ThisVertex;
    NewEdge->Vertices[1] = OtherVertex;
    this->Edges.push_back(NewEdge);
    return NewEdge;
}

TriangleType *MeshData :: AddTriangle(std::vector<double> n0, std::vector<double> n1, std::vector<double> n2)
{
    // Insert vertices and create a triangle using the indices returned
    int VertexIDs[3];

    VertexIDs[0] = this->VertexOctreeRoot->AddVertex(n0[0], n0[1], n0[2]);
    VertexIDs[1] = this->VertexOctreeRoot->AddVertex(n1[0], n1[1], n1[2]);
    VertexIDs[2] = this->VertexOctreeRoot->AddVertex(n2[0], n2[1], n2[2]);

    return this->AddTriangle({VertexIDs[0], VertexIDs[1], VertexIDs[2]});
}

TriangleType *MeshData :: AddTriangle(std::vector<int> VertexIDs)
{
    TriangleType *NewTriangle = new TriangleType;

    for (int i=0; i<3; i++) {
        if (i<2) {
            this->AddEdge({VertexIDs.at(i), VertexIDs.at(i+1)});
        } else {
            this->AddEdge({VertexIDs.at(i), VertexIDs.at(0)});
        }
    }

    for (int i=0; i<3; i++) {
        this->Vertices.at(VertexIDs.at(i))->AddTriangle(NewTriangle);
    }
    this->Triangles.push_back(NewTriangle);
    return NewTriangle;

}

}
