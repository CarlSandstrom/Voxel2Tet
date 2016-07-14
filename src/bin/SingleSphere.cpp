#include <iostream>
#include <string>
#include <stdio.h>

#include "Options.h"
#include "Voxel2Tet.h"
#include "CallbackImporter.h"
#include "TetGenExporter.h"
#include "OFFExporter.h"

/**
 * @brief Function being called from the Voxel2Tet class. The geometry of a sphere with radius .25
 * contained inside a cube, defined by coordinates (0,0,0) and (1,1,1) is implicitly given. The material
 * inside the sphere has ID 1 and the material inside the cube has ID 2.
 * @param x X coordinate
 * @param y Y coordinate
 * @param z Z coordinate
 * @return Material ID at coordinate
 */
int GiveMaterialIDByCoordinateSphere(double x, double y, double z)
{
    double r = .25;

    if ( sqrt( pow(x - .5, 2) + pow(y - .5, 2) + pow(z - .5, 2) ) < r ) {
        return 1;
    } else {
        return 2;
    }
}

int main(int argc, char *argv[])
{
    std :: map< std :: string, std :: string >DefaultOptions;
    voxel2tet :: Options *Options = new voxel2tet :: Options(argc, argv, DefaultOptions, {});

    voxel2tet :: Voxel2TetClass v2t(Options);

    double spacing = 0.025;
    double length = 1.0;
    int dimensions = std :: ceil(length / spacing);

    v2t.LoadCallback(& GiveMaterialIDByCoordinateSphere, { { 0, 0, 0 } }, { { spacing, spacing, spacing } }, { { dimensions, dimensions, dimensions } });
    v2t.Process();
    v2t.Tetrahedralize();
    v2t.ExportSurface("/tmp/SingleSphere0.vtp", voxel2tet :: FT_VTK);
    v2t.ExportSurface("/tmp/SingleSphere1.vtp", voxel2tet :: FT_VTK);
    v2t.ExportVolume("/tmp/SingleSphere.vtu", voxel2tet :: FT_VTK);
    int i=0;
    i++;

}
