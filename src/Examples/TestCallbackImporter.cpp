#include <iostream>
#include <string>
#include <stdio.h>

#include "Options.h"
#include "Voxel2Tet.h"
#include "CallbackImporter.h"
#include "TetGenExporter.h"
#include "OFFExporter.h"

int GiveMaterialIDByCoordinateSphere(double x, double y, double z)
{

    double r=.25;

    if (sqrt(pow(x-.5, 2) + pow(y-.5, 2) + pow(z-.5, 2)) < r) {
        return 1;
    } else {
        return 2;
    }

}

int GiveMaterialIDByCoordinateTwoSpheres(double x, double y, double z)
{

    double r=.25;
    double xc1=0.3, yc1=0.3, zc1=0.3;
    double xc2=0.3, yc2=0.3, zc2=0.7;

    bool InSphere1 = (sqrt(pow(x-xc1, 2) + pow(y-yc1, 2) + pow(z-zc1, 2)) < r);
    bool InSphere2 = (sqrt(pow(x-xc2, 2) + pow(y-yc2, 2) + pow(z-zc2, 2)) < r);

    if (InSphere1 & InSphere2) {
        return 0;
    } else if (InSphere1) {
        return 2;
    } else if (InSphere2) {
        return 3;
    } else {
        return 0;
    }

}

int main( int argc, char *argv[] )
{

    std::map <std::string, std::string> DefaultOptions;
    voxel2tet::Options *Options = new voxel2tet::Options(argc, argv, DefaultOptions);

    voxel2tet::Voxel2Tet v2t(Options);

    double spacing=0.05;
    double length=1.0;
    int dimensions= std::ceil(length/spacing);

    v2t.LoadCallback(&GiveMaterialIDByCoordinateSphere, {0,0,0}, {spacing, spacing, spacing}, {dimensions, dimensions, dimensions});
    v2t.Process();
    v2t.ExportSurface("/tmp/Voxel2Tet/TetGenTest.off", voxel2tet::FT_OFF);

}
