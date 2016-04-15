#include <iostream>
#include <string>
#include <stdio.h>

#include "Options.h"
#include "Voxel2Tet.h"

int main( int argc, char *argv[] )
{
    std::map <std::string, std::string> DefaultOptions;
    voxel2tet::Options *Options = new voxel2tet::Options(argc, argv, DefaultOptions);

    voxel2tet::Voxel2Tet v2t(Options);

    v2t.LoadFile(Options->GiveStringValue("i"));
    v2t.Process();
    v2t.Tetrahedralize();

    v2t.ExportSurface("/tmp/Voxel2Tet/TetGenTest.poly", voxel2tet::FT_Poly);
    v2t.ExportSurface("/tmp/Voxel2Tet/TetGenTest.off", voxel2tet::FT_OFF);
}
