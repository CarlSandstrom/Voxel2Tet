#ifndef VTKSTRUCTUREDREADER_H
#define VTKSTRUCTUREDREADER_H

#include <iostream>
#include <fstream>
#include <string.h>

#include "Importer.h"
#include "MiscFunctions.h"

namespace voxel2tet
{

class VTKStructuredReader : public Importer
{
private:
    std::string VersionInfo;
    std::string Title;
    std::string DataName;
    std::string TableName;

    double spacing_data[3];
    double origin_data[3];
    int dimensions_data[3];
    int celldata;

    int *GrainIdsData;

    BoundingBoxType BoundingBox;

public:
    VTKStructuredReader();
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

#endif // VTKSTRUCTUREDREADER_H
