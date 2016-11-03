#ifndef IMPORTER_H
#define IMPORTER_H

#include <vector>
#include <array>
#include <cmath>

#include "Options.h"

namespace voxel2tet
{

typedef int IntTriplet[3];

/**
  * Type for array of three doubles. Typically a coordinate.
  */
typedef struct
{
    /**
     * @brief Coordinate values
     */
    double c[3];
} DoubleTriplet;

/**
  * Holds the coordinates for the maximum (x, y, z) and the minimum (x, y, z)
  */
typedef struct
{
    /**
     * @brief Max values for (x, y, z)
     */
    std::array<double, 3> maxvalues;
    /**
     * @brief Min values for (x, y, z)
     */
    std::array<double, 3> minvalues;
} BoundingBoxType;

/**
  * Holds the number of voxels for the maximum (x, y, z) and the minimum (x, y, z)
  */
typedef struct
{
    /**
     * @brief Max number of voxels for (x, y, z)
     */
    std::array<int, 3> maxvalues;
    /**
     * @brief Min number of voxels for (x, y, z)
     */
    std::array<int, 3> minvalues;
} VoxelBoundingBoxType;

/**
 * @brief Abstract class for all Importer objects.
 */
class Importer
{
protected:
    /**
     * @brief Dimensions of a voxel
     */
    double spacing_data[3];

    /**
     * @brief Origin of voxel data
     */
    double origin_data[3];

    /**
     * @brief Number of voxels in each dimension
     */
    int dimensions_data[3];

    /**
     * @brief Bounding box
     */
    BoundingBoxType BoundingBox;

    /**
     * @brief Voxel data
     *
     * The data is stored in an array of integers where each integer is the material ID at that index. The order of the data is
     * X, Y, Z. I.e. the material id at index (xi, yi, zi) is given by
     *
     * MatID(xi, yi, zi) = zi*dimensions_data[0]*dimensions_data[1]+yi*dimensions_data[0]+xi
     *
     */
    int *GrainIdsData;

public:

    Importer()
    { UseCutOut = false; }

    /**
     * Tells if only a cutout of the data should be converted
     */
    bool UseCutOut;

    /**
     * Contains the voxel coordinates of the cutout
     */
    VoxelBoundingBoxType CutOut;

    /**
     * @brief Loads data from a file containing a voxel representation.
     *
     * As the file is loaded, internal variables for bounding boxes and such are updated. The data is stored in a private variable
     * and is accessed by GiveMaterialIDByCoordinate and GiveMaterialIDByIndex.
     *
     * @param FileName
     */
    virtual void LoadFile(std::string FileName) = 0;

    /**
     * @brief Returns the identifier of the material located at coordinate (x, y, z).
     * @param x X coordinate
     * @param y Y coordinate
     * @param z Z coordinate
     * @return Material ID
     */
    virtual int GiveMaterialIDByCoordinate(double x, double y, double z);

    /**
     * @brief Returns the identifier of the material located at index (xi, xy, zy).
     * @param xi Index in X direction
     * @param yi Index in Y direction
     * @param zi Index in Z direction
     * @return Material ID
     */
    virtual int GiveMaterialIDByIndex(int xi, int yi, int zi);

    /**
     * @brief Returns the dimensions of one voxel
     * @param spacing Array of 3
     */
    virtual void GiveSpacing(double spacing[3]);

    /**
     * @brief Returns bounding box for voxel data
     * @return Bounding box
     */
    virtual BoundingBoxType GiveBoundingBox();

    /**
     * @brief Returns the number of voxels in each dimension.
     * @param dimensions Array of 3 integers
     */
    virtual void GiveDimensions(int dimensions[3]);

    /**
     * @brief Convert indices to coordinate
     * @param xi Input. Index in X direction
     * @param yi Input. Index in Y direction
     * @param zi Input. Index in Z direction
     * @param Coordinate Output. Coordinate
     */
    virtual void GiveCoordinateByIndices(int xi, int yi, int zi, DoubleTriplet Coordinate);

    /**
     * Converts a coordinate to the index in the list of voxels
     * @param x
     * @param y
     * @param z
     * @param Indices
     */
    virtual void GiveIndicesByCoordinate(double x, double y, double z, IntTriplet Indices);

    /**
     * @brief Returns origin of voxel data
     * @param origin Coordinate. Array of 3 doubles
     */
    virtual void GiveOrigin(double origin[3]);
};
}

#endif // IMPORTER_H
