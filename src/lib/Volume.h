#ifndef VOLUME_H
#define VOLUME_H

#include <vector>

#include "Surface.h"
#include "MeshComponents.h"
#include "MiscFunctions.h"

namespace voxel2tet
{

class Volume
{
public:
    Volume();

    Volume(int Phase);

    std::vector<Surface*> Surfaces;

    std::vector<TriangleType*> GiveTriangles();

    double ComputeVolume();

    int Phase;
};

}
#endif // VOLUME_H
