#ifndef CALLBACKIMPORTER_H
#define CALLBACKIMPORTER_H

#include <array>
#include "Importer.h"
#include "MiscFunctions.h"

namespace voxel2tet
{

typedef int (*cbMaterialIDByCoordinate) (double, double, double);

class CallbackImporter: public Importer
{
private:
    double eps=1e-6;
    cbMaterialIDByCoordinate MaterialByCoordinate;
    std::array<double,3> spacing_data;
    std::array<double,3> origin_data;
    std::array<int,3> dimensions_data;
    BoundingBoxType BoundingBox;

public:
    CallbackImporter(cbMaterialIDByCoordinate MaterialIdByCoordinateCallback, std::array<double,3> Origin, std::array<double,3> Spacing, std::array<int,3> Dimensions);
    void LoadFile(std::string FileName) {STATUS("This class is made for loading voxeldata by evaluating functions. Loading of files are not possible\n", 0); }
    virtual int GiveMaterialIDByCoordinate(double x, double y, double z);
    virtual int GiveMaterialIDByIndex(int xi, int yi, int zi);
    virtual void GiveSpacing(double spacing[3]);
    virtual BoundingBoxType GiveBoundingBox();
    virtual void GiveDimensions(int dimensions[3]);
    virtual void GiveCoordinateByIndices(int xi, int yi, int zi, DoubleTriplet Coordinate);
    virtual void GiveOrigin(double origin[3]);
};

}
#endif // CALLBACKIMPORTER_H
