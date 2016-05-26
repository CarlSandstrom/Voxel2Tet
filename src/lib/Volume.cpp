#include "Volume.h"

namespace voxel2tet
{

Volume::Volume()
{
    this->Phase = -1;
}

Volume::Volume(int Phase)
{
    this->Phase = Phase;
}

}
