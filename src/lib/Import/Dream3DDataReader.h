#ifndef CCBUILDERDATAREADER_H
#define CCBUILDERDATAREADER_H

#include <string>

#include "H5Cpp.h"
#include "Importer.h"

namespace voxel2tet {
class Dream3DDataReader : public Importer
{
private:
    H5 :: Group *VoxelDataContainer;
    std :: string DataContainer;
    std :: string MaterialGroup;
public:
    Dream3DDataReader();
    Dream3DDataReader(std :: string DataContainer, std :: string MaterialGroup);
    void LoadFile(std :: string FileName);
};
}

#endif // CCBUILDERDATAREADER_H
