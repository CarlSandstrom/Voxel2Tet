#include <iostream>

#include "Options.h"
#include "Voxel2Tet.h"
#include "TetGenExporter.h"
#include "OFFExporter.h"

/**
 * @brief Function being called by Voxel2Tet. The geometry of a cube, defined by the corners (0,0,0)
 * and (1,1,1) is implicitly given. If the coordinate is located inside the cube, 1 is returned
 * otherwise a 0 is returned.
 *
 * @param x X coordinate
 * @param y Y coordinate
 * @param z Z coordinate
 * @return Material ID at coordinate.
 *
 */
int GiveMaterialIDByCoordinateCube(double x, double y, double z)
{
    if ((x <= 1.) & (x >= 0.) & (y <= 1.) & (y >= 0.) & (z <= 1.) & (z >= 0)) {
        return 1;
    } else {
        return 0;
    }
}

int main(int argc, char *argv[])
{
    std::map<std::string, std::string> DefaultOptions;

    // The Options object is created with the command line arguments, thus, we can use the same
    // arguments as for the Voxel2Tet binary.
    voxel2tet::Options *Options = new voxel2tet::Options(argc, argv, DefaultOptions, {});

    voxel2tet::Voxel2TetClass v2t(Options);

    // Set dimensions
    double spacing = 0.4;
    double length = 1.0;
    int dimensions = std::ceil(length / spacing);

    v2t.LoadCallback(&GiveMaterialIDByCoordinateCube, {{0, 0, 0}}, {{spacing, spacing, spacing}},
                     {{dimensions, dimensions, dimensions}});
    v2t.Process();
    v2t.Tetrahedralize();
    v2t.ExportSurface("/tmp/Cube.vtp", voxel2tet::FT_VTK);
    v2t.ExportVolume("/tmp/Cube.vtu", voxel2tet::FT_VTK);

}
