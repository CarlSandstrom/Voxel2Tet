#ifndef IMPORTER_H
#define IMPORTER_H

#include <vector>
#include <array>
#include <cmath>

#include "Options.h"

namespace voxel2tet {

typedef int IntTriplet[3];

typedef struct {
    double c[3];
} DoubleTriplet;

typedef struct {
    std::array<double,3> maxvalues;
    std::array<double,3> minvalues;
} BoundingBoxType;

class Importer
{
protected:
    double spacing_data[3];
    double origin_data[3];
    int dimensions_data[3];
    BoundingBoxType BoundingBox;
    int *GrainIdsData;

public:
    virtual void LoadFile(std::string FileName) = 0;
    virtual int GiveMaterialIDByCoordinate(double x, double y, double z);
    virtual int GiveMaterialIDByIndex(int xi, int yi, int zi);
    virtual void GiveSpacing(double spacing[3]);
    virtual BoundingBoxType GiveBoundingBox();
    virtual void GiveDimensions(int dimensions[3]);
    virtual void GiveCoordinateByIndices(int xi, int yi, int zi, DoubleTriplet Coordinate);
    virtual void GiveOrigin(double origin[3]);
};

}

#endif // IMPORTER_H
