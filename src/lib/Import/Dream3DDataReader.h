#ifndef CCBUILDERDATAREADER_H
#define CCBUILDERDATAREADER_H

#include <string>

#include "H5Cpp.h"
#include "Importer.h"

namespace voxel2tet {

/**
 * @brief The Dream3DDataReader class reads data from a Dream3D file.
 */
class Dream3DDataReader : public Importer
{
private:
    H5 :: Group *VoxelDataContainer;
    std :: string DataContainer;
    std :: string MaterialGroup;
public:
    Dream3DDataReader();

    /**
     * @brief Constructor with special parameters.
     *
     * Since the Dream 3D format differs some between versions, options for giving the names of the DataContainer
     * and the DataGroup locations are provided.
     *
     * @param DataContainer Name of DataContainer
     * @param MaterialGroup Name of Group
     */
    Dream3DDataReader(std :: string DataContainer, std :: string MaterialGroup);
    void LoadFile(std :: string FileName);
};
}

#endif // CCBUILDERDATAREADER_H
