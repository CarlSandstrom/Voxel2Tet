#include <iostream>
#include <string>
#include <stdio.h>

#include "Options.h"
#include "Voxel2Tet.h"
#include "CallbackImporter.h"
#include "TetGenExporter.h"
#include "OFFExporter.h"

int GiveMaterialIDByCoordinate(double x, double y, double z)
{

    if (sqrt(pow(x-.5, 2) + pow(y-.5, 2) + pow(z-.5, 2)) < .5) {
        return 1;
    } else {
        return 2;
    }

}

int main( int argc, char *argv[] )
{
    std::map <std::string, std::string> DefaultOptions;
    voxel2tet::Options *Options = new voxel2tet::Options(argc, argv, DefaultOptions);

    voxel2tet::Voxel2Tet v2t(Options);
    v2t.LoadCallback(&GiveMaterialIDByCoordinate, {0,0,0}, {.1, .1, .1}, {10, 10, 10});
    v2t.Process();
    v2t.ExportTetGenFile("/tmp/Voxel2Tet/TetGenTest.poly");
    v2t.ExportOFF("/tmp/Voxel2Tet/TetGenTest.off");


}
