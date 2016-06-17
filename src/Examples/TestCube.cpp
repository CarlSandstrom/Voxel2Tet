#include <iostream>
#include <string>
#include <stdio.h>

#include "Options.h"
#include "Voxel2Tet.h"
#include "CallbackImporter.h"
#include "TetGenExporter.h"
#include "OFFExporter.h"

int GiveMaterialIDByCoordinateCube(double x, double y, double z)
{

    if ((x<=1.) & (x>=0.) & (y<=1.) & (y>=0.) & (z<=1.) & (z>=0)) {
        return 1;
    } else {
        return 0;
    }

}

int main( int argc, char *argv[] )
{
    std::map <std::string, std::string> DefaultOptions;
    voxel2tet::Options *Options = new voxel2tet::Options(argc, argv, DefaultOptions, {});

    voxel2tet::Voxel2TetClass v2t(Options);

    double spacing=0.1;
    double length=1.0;
    int dimensions= std::ceil(length/spacing);

    v2t.LoadCallback(&GiveMaterialIDByCoordinateCube, {{0,0,0}}, {{spacing, spacing, spacing}}, {{dimensions, dimensions, dimensions}});
    v2t.Process();
    v2t.ExportSurface("/tmp/Voxel2Tet/TetGenTest.off", voxel2tet::FT_OFF);


}
