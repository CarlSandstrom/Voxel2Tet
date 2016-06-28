#include <iostream>
#include <string>
#include <stdio.h>

#include "Voxel2Tet.h"
#include "Exporter.h"

int main(int argc, char *argv[])
{
    std :: map< std :: string, std :: string >DefaultOptions;
    voxel2tet :: Options *Options = new voxel2tet :: Options(argc, argv, DefaultOptions, {});

    voxel2tet::Voxel2TetClass v2t(Options);
    v2t.LoadFile("/home/carl/dev/testprojects/FindFaultyCC/testfile_2.0_1.dream3d");

    int v467id = v2t.Mesh->VertexOctreeRoot->AddVertex(0.236127331852912903,	0.86643218994140625, 1.34924781322479248);
    int v190id = v2t.Mesh->VertexOctreeRoot->AddVertex(0, 0.789500176906585693, 1.88603127002716064);
    int v184id = v2t.Mesh->VertexOctreeRoot->AddVertex(0,	0.862375915050506592,	1.4135744571685791);
    int v183id = v2t.Mesh->VertexOctreeRoot->AddVertex(0.107607662677764893,	0.886709988117218018,	1.37026321887969971);
    int v135id = v2t.Mesh->VertexOctreeRoot->AddVertex(0.0842119902372360229,	0.532879352569580078,	1.6611332893371582);
    int v132id = v2t.Mesh->VertexOctreeRoot->AddVertex(0,	0.638925135135650635,	1.49731588363647461);

    voxel2tet::TriangleType *t1 = v2t.Mesh->AddTriangle({{v467id, v190id, v184id}});
    voxel2tet::TriangleType *t2 = v2t.Mesh->AddTriangle({{v183id, v135id, v132id}});

    v2t.Mesh->CheckTrianglePenetration(t1, t2);

    v2t.ExportSurface("/tmp/test.vtp", voxel2tet::FT_VTK);

}
