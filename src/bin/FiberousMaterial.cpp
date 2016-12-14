#include <iostream>
#include <armadillo>

#include "Options.h"
#include "Voxel2Tet.h"

std::vector<std::array<double, 3>> Coordinates0;
std::vector<std::array<double, 3>> Coordinates1;

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
    double r = .1;

    arma::vec x0 = {x, y, z};

    for (size_t i = 0; i < Coordinates0.size(); i++) {
        arma::vec x1 = {Coordinates0[i][0], Coordinates0[i][1], Coordinates0[i][2]};
        arma::vec x2 = {Coordinates1[i][0], Coordinates1[i][1], Coordinates1[i][2]};
        double nominator = arma::norm(arma::cross(x2 - x1, x1 - x0));
        double denominator = arma::norm(x2 - x1);
        double d = nominator / denominator;
        if (d < r) return 1;
    }

    return 2;

}

int main(int argc, char *argv[])
{
    std::map<std::string, std::string> DefaultOptions;
    voxel2tet::Options *Options = new voxel2tet::Options(argc, argv, DefaultOptions, {});

    Options->SetKey("spring_alpha", 10.0);

    voxel2tet::Voxel2TetClass v2t(Options);

    double spacing = 0.025;
    std::array<double, 3> length = {{1.0, 1.0, 3.0}};
    std::array<int, 3> dimensions;
    for (int i = 0; i < 3; i++) {
        dimensions[i] = std::ceil(length[i] / spacing);
    }

    // Generate random coordinates on each opposite side of the box
    int NumberOfFibers = 5;

    srand(time(NULL));

    for (int i = 0; i < NumberOfFibers; i++) {
        std::array<double, 3> c0;
        std::array<double, 3> c1;
        c0[2] = 0;
        c1[2] = length[2];

        for (int j = 0; j < 2; j++) {
            c0[j] = (double(std::rand()) / RAND_MAX) * length[j];
            c1[j] = (double(std::rand()) / RAND_MAX) * length[j];
        }
        Coordinates0.push_back(c0);
        Coordinates1.push_back(c1);
    }

    v2t.LoadCallback(&GiveMaterialIDByCoordinateMultiSphere, {{0, 0, 0}}, {{spacing, spacing, spacing}}, dimensions);
    v2t.Process();
    v2t.Tetrahedralize();

    v2t.ExportSurface("FiberNetwork.vtp", voxel2tet::FT_VTK);
    v2t.ExportVolume("FiberNetwork.vtu", voxel2tet::FT_VTK);

}
