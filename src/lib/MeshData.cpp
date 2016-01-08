#include <vector>
#include "MeshData.h"
#include "MiscFunctions.h"
#include "TetGenExporter.h"
#include "OFFExporter.h"

namespace voxel2tet
{


MeshData::MeshData(BoundingBoxType BoundingBox)
{
    this->BoundingBox = BoundingBox;
    this->VertexOctreeRoot = new VertexOctreeNode(this->BoundingBox, &this->Vertices, 0);
    this->TriangleCounter = 0;
}

MeshData::~MeshData()
{
    for (unsigned int i=0; i<this->Edges.size(); i++) {
        delete this->Edges.at(i);
    }

    for (auto t: this->Triangles) {
        delete t;
    }

    for (auto v: this->Vertices) delete v;

    delete this->VertexOctreeRoot;
}

void MeshData :: ExportVTK(std::string FileName)
{
    STATUS ("Export to %s\n", FileName.c_str());
    VTKExporter exporter(&this->Triangles, &this->Vertices, &this->Edges);
    exporter.WriteData(FileName);
}

void MeshData::ExportTetgen(std::string FileName)
{
    TetGenExporter exporter(&this->Triangles, &this->Vertices, &this->Edges);
    exporter.WriteData(FileName);
}

void MeshData::ExportOFF(std::string FileName)
{
    OFFExporter exporter(&this->Triangles, &this->Vertices, &this->Edges);
    exporter.WriteData(FileName);
}

EdgeType *MeshData :: AddEdge(std::vector<int> VertexIDs)
{

    // Check if edge already exists
    int ThisVertexID = VertexIDs.at(0);
    int OtherVertexID = VertexIDs.at(1);
    VertexType *ThisVertex = this->Vertices.at(ThisVertexID);
    VertexType *OtherVertex = this->Vertices.at(OtherVertexID);

    std::vector <EdgeType*> Edges = this->Vertices.at(ThisVertexID)->Edges;
    for (auto Edge: Edges) {
        if ( ( (Edge->Vertices[0] == ThisVertex) & (Edge->Vertices[1] == OtherVertex) ) | ( (Edge->Vertices[1] == ThisVertex) & (Edge->Vertices[0] == OtherVertex) )) {
            LOG("Edge %p already exists\n", Edge);
            return Edge;
        }
    }

    EdgeType *NewEdge = new EdgeType;
    LOG("Create edge from Vertex IDs %u and %u: %p\n", VertexIDs.at(0), VertexIDs.at(1), NewEdge);

    // Update edge
    NewEdge->Vertices[0] = ThisVertex;
    NewEdge->Vertices[1] = OtherVertex;
    this->Edges.push_back(NewEdge);

    // Update vertices
    ThisVertex->AddEdge(NewEdge);
    OtherVertex->AddEdge(NewEdge);
    return NewEdge;
}

void MeshData :: RemoveTriangle(TriangleType *t)
{
    for (VertexType *v: t->Vertices) {
        v->RemoveTriangle(t);
    }
    this->Triangles.erase(std::remove(this->Triangles.begin(), this->Triangles.end(), t), this->Triangles.end());
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
    LOG ("Create triangle %p from vertices (%u, %u, %u)@(%p, %p, %p)\n", NewTriangle, VertexIDs.at(0), VertexIDs.at(1), VertexIDs.at(2),
         this->Vertices.at(VertexIDs.at(0)), this->Vertices.at(VertexIDs.at(1)), this->Vertices.at(VertexIDs.at(2)));

    for (int i=0; i<3; i++) {
        if (i<2) {
            this->AddEdge({VertexIDs.at(i), VertexIDs.at(i+1)});
        } else {
            this->AddEdge({VertexIDs.at(i), VertexIDs.at(0)});
        }
    }

    // Update vertices and triangles with references to each other
    for (int i=0; i<3; i++) {
        this->Vertices.at(VertexIDs.at(i))->AddTriangle(NewTriangle);
        NewTriangle->Vertices[i]=this->Vertices.at(VertexIDs.at(i));
    }

    return AddTriangle(NewTriangle);

}

TriangleType *MeshData :: AddTriangle(TriangleType *NewTriangle)
{
    NewTriangle->UpdateNormal();
    NewTriangle->ID = TriangleCounter;
    for (VertexType *v: NewTriangle->Vertices) {
        v->AddTriangle(NewTriangle);
    }
    TriangleCounter++;
    this->Triangles.push_back(NewTriangle);
    return NewTriangle;

}

}
