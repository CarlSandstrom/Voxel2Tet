#ifndef IMPORTER_H
#define IMPORTER_H

#include <vector>
#include <array>
#include <cmath>

#include "Options.h"

namespace voxel2tet {

typedef int IntTriplet [ 3 ];

typedef struct {
    double c [ 3 ];
} DoubleTriplet;

typedef struct {
    std :: array< double, 3 >maxvalues;
    std :: array< double, 3 >minvalues;
} BoundingBoxType;

/**
 * @brief Abstract class for all Importer objects.
 */
class Importer
{
protected:
    double spacing_data [ 3 ];
    double origin_data [ 3 ];
    int dimensions_data [ 3 ];
    BoundingBoxType BoundingBox;
    int *GrainIdsData;

public:

    /**
     * @brief Loads data from a file containing a voxel representation.
     *
     * As the file is loaded, internal variables for bounding boxes and such are updated. The data is stored in a private variable
     * and is accessed by GiveMaterialIDByCoordinate and GiveMaterialIDByIndex.
     *
     * @param FileName
     */
    virtual void LoadFile(std :: string FileName) = 0;

    /**
     * @brief Returns the identifier of the material located at coordinate (x, y, z).
     * @param x X coordinate
     * @param y Y coordinate
     * @param z Z coordinate
     * @return
     */
    virtual int GiveMaterialIDByCoordinate(double x, double y, double z);
    virtual int GiveMaterialIDByIndex(int xi, int yi, int zi);
    virtual void GiveSpacing(double spacing [ 3 ]);
    virtual BoundingBoxType GiveBoundingBox();
    virtual void GiveDimensions(int dimensions [ 3 ]);
    virtual void GiveCoordinateByIndices(int xi, int yi, int zi, DoubleTriplet Coordinate);
    virtual void GiveOrigin(double origin [ 3 ]);
};
}

#endif // IMPORTER_H
