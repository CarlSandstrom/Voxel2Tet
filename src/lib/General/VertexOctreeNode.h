#ifndef VERTEXOCTREENODE_H
#define VERTEXOCTREENODE_H

#include <vector>

#include <Importer.h>
#include "MeshComponents.h"
#include "MiscFunctions.h"

namespace voxel2tet
{

/**
 * @brief Provides octree functionality for accessing and storing verices to improve performance.
 *
 * For information on octrees, cf https://en.wikipedia.org/wiki/Octree
 *
 */
class VertexOctreeNode
{
private:
    BoundingBoxType BoundingBox;
    int level;
    int maxvertices;
    double eps;

    // List of indices pointing to this->Vertices contained within this node
    // TODO: Is it neccessary to use indices? Why not pionters?
    std::vector<int> VertexIds;

    // Splits node into eight nodes
    void split();

    // Determines wether a coordinate is located within this node
    bool IsInBoundingBox(double x, double y, double z);

public:

    /**
     * @brief Contructor.
     * @param BoundingBox Bounding box for this node. If root node, this is equal to the bounding box of the complete structure.
     * @param Vertices Pointer to a list for vertices. This is the list of vertices used henceforth.
     * @param level Level of this node. 0 is for root. Other levels are taken care of by its parent.
     */
    VertexOctreeNode(BoundingBoxType BoundingBox, std::vector<VertexType *> *Vertices, int level);

    ~VertexOctreeNode();

    /**
     * @brief Adds a vertex at a specified coordinate to the structure.
     *
     * If a vertex at the spcified coordinate already exists, the index of that vertex is returned, otherwise, the index of the new vertex is returned.
     *
     * If the current node is full, it is split into eight sub-nodes and all vertices are passed on down to the underlying nodes.
     *
     * @param x X coordinate
     * @param y Y coordinate
     * @param z Z coordinate
     * @return Index of vertex
     */
    int AddVertex(double x, double y, double z);

    /**
     * @brief List of vertices.
     *
     * This is the complete list of vertices, not only for the vertices in this node.
     */
    std::vector<VertexType *> *Vertices;

    /**
     * @brief List of nodes owned by this node
     */
    std::vector<VertexOctreeNode *> children;

    /**
     * @brief Find a vertex in the structure by coordinate. If a vertex at the specified coordinate exists, return a pointer to that vertex object.
     * @param x X coordinate
     * @param y Y coordinate
     * @param z Z coordinate
     * @return Pointer to vertex object
     */
    VertexType *FindVertexByCoords(double x, double y, double z);

    /**
     * @brief Produce a list of vertices which are located within a sphere
     * @param x Center X coordinate of sphere
     * @param y Center Y coordinate of sphere
     * @param z Center Z coordinate of sphere
     * @param r Radius of sphere
     * @return List of vertices
     */
    std::vector<VertexType *> GiveVerticesWithinSphere(double x, double y, double z, double r);

    /**
     * @brief Print information from this node and its children. For each level, a tab is added in front of the information to illustrate the hierachy.
     */
    void printself();
};
}
#endif // VERTEXOCTREENODE_H
