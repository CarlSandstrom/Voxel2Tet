#ifndef HDF5DATAREADER_H
#define HDF5DATAREADER_H

#include <string>

#include "H5Cpp.h"
#include "Importer.h"

namespace voxel2tet {

class hdf5DataReader: public Importer
{
private:
    H5::Group *VoxelDataContainer;
    double spacing_data[3];
    double origin_data[3];
    int dimensions_data[3];
    BoundingBoxType BoundingBox;
    int *GrainIdsData;
public:
    hdf5DataReader();
    void LoadFile(std::string FileName);
    virtual int GiveMaterialIDByCoordinate(double x, double y, double z);
    virtual int GiveMaterialIDByIndex(int xi, int yi, int zi);
    virtual void GiveSpacing(double spacing[3]);
    virtual BoundingBoxType GiveBoundingBox();
    virtual void GiveDimensions(int dimensions[3]);
    virtual void GiveCoordinateByIndices(int xi, int yi, int zi, DoubleTriplet Coordinate);
    virtual void GiveOrigin(double origin[3]);
};

}

#endif // HDF5DATAREADER_H
