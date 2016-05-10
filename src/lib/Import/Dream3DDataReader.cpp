#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<cmath>

#include "Dream3DDataReader.h"
#include "MiscFunctions.h"

namespace voxel2tet {

Dream3DDataReader::Dream3DDataReader()
{
    this->DataContainer = "DataContainer";
    this->MaterialGroup = "GrainIds";
}

Dream3DDataReader::Dream3DDataReader(std::string DataContainer, std::string MaterialGroup) : Importer()
{
    this->DataContainer = DataContainer;
    this->MaterialGroup = MaterialGroup;
}

void Dream3DDataReader :: LoadFile(std::string FileName)
{

    LOG ("Open file %s\n", FileName.c_str());
    H5::H5File file( FileName, H5F_ACC_RDONLY );
    H5::Group DataContainers;

    DataContainers = H5 :: Group ( file.openGroup("DataContainers") );
    VoxelDataContainer = new H5 :: Group ( DataContainers.openGroup(this->DataContainer) );


    // Load dimensional/spatial data
    VoxelDataContainer->openGroup("_SIMPL_GEOMETRY").openDataSet("SPACING").read( spacing_data, H5::PredType::NATIVE_DOUBLE);
    VoxelDataContainer->openGroup("_SIMPL_GEOMETRY").openDataSet("ORIGIN").read( origin_data, H5::PredType::NATIVE_DOUBLE);
    VoxelDataContainer->openGroup("_SIMPL_GEOMETRY").openDataSet("DIMENSIONS").read( dimensions_data, H5::PredType::NATIVE_INT);

    for (int i=0; i<3; i++) this->BoundingBox.minvalues[i] = this->origin_data[i];
    for (int i=0; i<3; i++) this->BoundingBox.maxvalues[i] = this->origin_data[i]+this->dimensions_data[i]*this->spacing_data[i];

    // Load voxel data.
    H5 :: DataSet GrainIds = VoxelDataContainer->openGroup("CellData").openDataSet(this->MaterialGroup);
    H5 :: DataSpace space = GrainIds.getSpace();

    int Ndims = space.getSimpleExtentNdims();
    hsize_t dims[Ndims];

    space.getSimpleExtentDims(dims);

    GrainIdsData = (int*) malloc(sizeof(int)*dims[0]*dims[1]*dims[2]*dims[3]);
    GrainIds.read( GrainIdsData, H5::PredType::NATIVE_INT);

}

int Dream3DDataReader :: GiveMaterialIDByCoordinate(double x, double y, double z)
{
    double coords[3]={x, y, z};
    int indices[3];
    for (int i=0; i<3; i++) {
        double intcoord = (coords[i]-this->origin_data[i])/this->spacing_data[i];
        indices[i] =  floor(intcoord);
    }
    return this->GiveMaterialIDByIndex(indices[0], indices[1], indices[2]);
}

int Dream3DDataReader :: GiveMaterialIDByIndex(int xi, int yi, int zi)
{
    if (xi == -1) {
        return -1;
    } else if (xi == this->dimensions_data[0]) {
        return -2;
    } else if (yi == -1) {
        return -3;
    } else if (yi == this->dimensions_data[1]) {
        return -4;
    } else if (zi == -1) {
        return -5;
    } else if (zi == this->dimensions_data[2]) {
        return -6;
    }

    int index = zi*this->dimensions_data[1]*this->dimensions_data[0] + yi*this->dimensions_data[0] + xi;
    return this->GrainIdsData[index];

}

void Dream3DDataReader :: GiveSpacing(double spacing[3])
{
    for (int i=0; i<3; i++) {spacing[i] = this->spacing_data[i];}
}

BoundingBoxType Dream3DDataReader::GiveBoundingBox()
{
    BoundingBoxType BoundingBox;
    for (int i=0; i<3; i++) {
        BoundingBox.maxvalues[i] = this->BoundingBox.maxvalues[i];
        BoundingBox.minvalues[i] = this->BoundingBox.minvalues[i];
    };
    return BoundingBox;
}

void Dream3DDataReader :: GiveDimensions(int dimensions[3])
{
    for (int i=0; i<3; i++) {dimensions[i] = this->dimensions_data[i];}
}

void Dream3DDataReader :: GiveCoordinateByIndices(int xi, int yi, int zi, DoubleTriplet Coordinate)
{
    Coordinate.c[0] = xi*this->spacing_data[0] + this->origin_data[0];
    Coordinate.c[1] = yi*this->spacing_data[1] + this->origin_data[1];
    Coordinate.c[2] = zi*this->spacing_data[2] + this->origin_data[2];
}

void Dream3DDataReader :: GiveOrigin(double origin[3])
{
    for (int i=0; i<3; i++) {origin[i] = this->origin_data[i];}
}

}

