#include "AbaqusExporter.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <map>
#include <vector>

namespace voxel2tet
{

AbaqusExporter::AbaqusExporter(std :: vector< TriangleType * > *Triangles, std :: vector< VertexType * > *Vertices, std :: vector< EdgeType * > *Edges, std :: vector< TetType * > *Tets) :
    Exporter(Triangles, Vertices, Edges, Tets)
{

}

void AbaqusExporter::WriteVolumeData(std :: string Filename)
{
    STATUS( "Write .in (oofem) file %s\n", Filename.c_str() );
    std :: ofstream AbaqusFile;
    AbaqusFile.open(Filename);


}

}
