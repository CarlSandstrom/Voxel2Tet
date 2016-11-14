#include <iostream>

#include "Options.h"
#include "Voxel2Tet.h"

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
    double r = .4;

    if (sqrt(pow(x - .50, 2) + pow(y - .50, 2) + pow(z - .50, 2)) < r) {
        return 1;
    } else {
        return 2;
    }
}

int main(int argc, char *argv[])
{
    std::map<std::string, std::string> DefaultOptions;
    voxel2tet::Options *Options = new voxel2tet::Options(argc, argv, DefaultOptions, {});

    Options->SetKey("spring_c_factor", 2.0);
    Options->SetKey("spring_alpha", 4.0);

    voxel2tet::Voxel2TetClass v2t(Options);

    int numberofvoxels = 20;
    double length = 1.0;
    double spacing = length/double(numberofvoxels);
    int dimensions = std::ceil(length / spacing);

    v2t.LoadCallback(&GiveMaterialIDByCoordinateSphere, {{0, 0, 0}}, {{spacing, spacing, spacing}},
                     {{dimensions, dimensions, dimensions}});
    v2t.Process();
    v2t.Tetrahedralize();
    v2t.ExportSurface("/tmp/SingleSphere.vtp", voxel2tet::FT_VTK);
    v2t.ExportVolume("/tmp/SingleSphere.in", voxel2tet::FT_ABAQUS);
    v2t.ExportAllVolumes();

}
