#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "Options.h"
#include "Voxel2Tet.h"
#include "CallbackImporter.h"
#include "TetGenExporter.h"
#include "OFFExporter.h"

std::vector<std::array<double, 3>> Coordinates;

/**
 * @brief Callback function being called from Voxel2Tet. Describes the geometry implicitly by returning
 * the ID of the material in the given coordinate.
 * @param x X coordinate
 * @param y Y coordinate
 * @param z Z coordinate
 * @return Material ID
 */
int GiveMaterialIDByCoordinateMultiSphere(double x, double y, double z)
{
    double r = .3;

    for (std::array<double, 3> C: Coordinates) {
        if ( sqrt( pow(x - C[0], 2) + pow(y - C[1], 2) + pow(z - C[2], 2) ) < r ) {
            return 1;
        }
    }

    return 2;

}

int main(int argc, char *argv[])
{
    std :: map< std :: string, std :: string >DefaultOptions;
    voxel2tet :: Options *Options = new voxel2tet :: Options(argc, argv, DefaultOptions, {});

    voxel2tet :: Voxel2TetClass v2t(Options);

    double spacing = 0.05;
    std::array<double, 3> length = {3.0, 2.0, 1.0};
    std::array<int, 3> dimensions;
    for (int i=0; i<3; i++) {
        dimensions[i] = std::ceil(length[i] / spacing);
    }

    // Generate random spheres
    int NumberOfSpheres = 20;

    srand(time(NULL));

    for (int i=0; i<NumberOfSpheres; i++) {
        std::array<double, 3> cp;
        for (int j=0; j<3; j++) {
            cp[j] = (double (std::rand()) / RAND_MAX) * length[j];
        }
        Coordinates.push_back(cp);
    }

    v2t.LoadCallback(& GiveMaterialIDByCoordinateMultiSphere, { { 0, 0, 0 } }, { { spacing, spacing, spacing } }, dimensions);
    v2t.Process();
    v2t.Tetrahedralize();

    v2t.ExportSurface("/tmp/MultiSphere.vtp", voxel2tet :: FT_VTK);
    v2t.ExportVolume("/tmp/MultiSphere.vtu", voxel2tet :: FT_VTK);

}
