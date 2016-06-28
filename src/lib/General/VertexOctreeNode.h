#ifndef VERTEXOCTREENODE_H
#define VERTEXOCTREENODE_H

#include <vector>

#include <Importer.h>
#include "MeshComponents.h"
#include "MiscFunctions.h"

namespace voxel2tet
{
class VertexOctreeNode
{
private:
    BoundingBoxType BoundingBox;
    int level;
    int maxvertices;
    double eps;

    // List of indices pointing to this->Vertices contained within this node
    // TODO: Is it neccessary to use indices? Why not pionters?
    std :: vector< int >VertexIds;

    // Splits node into eight nodes
    void split();

    // Determines wether a coordinate is located within this node
    bool IsInBoundingBox(double x, double y, double z);
public:
    VertexOctreeNode(BoundingBoxType BoundingBox, std :: vector< VertexType * > *Vertices, int level);
    ~VertexOctreeNode();

    // Adds vertex to this node or somwhere down the hierarchy
    int AddVertex(double x, double y, double z);

    // List of vertices. Note: This is the complete list, not only the vertices in the current node
    std :: vector< VertexType * > *Vertices;

    // List of child nodes
    std :: vector< VertexOctreeNode * >children;

    // Find vertex object by coordinates
    VertexType *FindVertexByCoords(double x, double y, double z);

    std :: vector< VertexType * >GiveVerticesWithinSphere(double x, double y, double z, double r);

    void printself();
};
}
#endif // VERTEXOCTREENODE_H
