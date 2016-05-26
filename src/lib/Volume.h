#ifndef VOLUME_H
#define VOLUME_H

#include <vector>

#include "Surface.h"
#include "MiscFunctions.h"

namespace voxel2tet
{

class Volume
{
public:
    Volume();

    Volume(int Phase);

    std::vector<Surface*> Surfaces;

    int Phase;
};

}
#endif // VOLUME_H
