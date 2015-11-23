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
    v2t.LoadData();
    v2t.Process();
    v2t.Mesh->ExportVTK("/tmp/test0.vtp");
    //v2t.Mesh->ExportVTK("/tmp/test1.vtp");
}
