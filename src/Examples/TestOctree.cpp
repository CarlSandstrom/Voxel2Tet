#include<vector>

#include<VertexOctreeNode.h>
#include<Importer.h>

int main()
{
    voxel2tet::BoundingBoxType bb;
    for (int i = 0; i<3; i++) {
        bb.maxvalues[i]=1.0;
        bb.minvalues[i]=0.0;
    }

    std::vector <voxel2tet::DoubleTriplet> VertexList;
    voxel2tet::VertexOctreeNode RootNode(bb, &VertexList, 0);
    RootNode.AddVertex(0.1, 0.1, 0.1);
    RootNode.AddVertex(0.11, 0.15, 0.1);
    RootNode.AddVertex(0.13, 0.15, 0.1);
    RootNode.AddVertex(0.19, 0.12, 0.1);
    RootNode.AddVertex(0.11, 0.15, 0.21);
    RootNode.AddVertex(0.1, 0.115, 0.1);
    RootNode.AddVertex(0.1, 0.1, 0.7);
    RootNode.AddVertex(0.11, 0.125, 0.1);
    RootNode.AddVertex(0.9, 0.15, 0.1);
    RootNode.printself();

}
