#ifndef IMPORTER_H
#define IMPORTER_H

#include <vector>

namespace voxel2tet {

typedef int IntTriplet[3];

typedef struct {
    double c[3];
} DoubleTriplet;

typedef struct {
    double maxvalues[3];
    double minvalues[3];
} BoundingBoxType;

class Importer
{
private:
    std :: vector <double> Spacing;
public:
    Importer();
    virtual int GiveMaterialIDByCoordinate(double x, double y, double z) = 0;
    virtual int GiveMaterialIDByIndex(int xi, int yi, int zi) = 0;
    virtual void GiveSpacing(double spacing[3]) = 0;
    virtual BoundingBoxType GiveBoundingBox() = 0;
    virtual void GiveDimensions(int dimensions[3]) = 0;
    virtual void GiveCoordinateByIndices(int xi, int yi, int zi, DoubleTriplet Coordinate) = 0;
    virtual void GiveOrigin(double origin[3]) = 0;
};

}

#endif // IMPORTER_H
