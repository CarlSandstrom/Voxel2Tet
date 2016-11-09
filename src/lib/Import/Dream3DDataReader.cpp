#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#include "Dream3DDataReader.h"
#include "MiscFunctions.h"

namespace voxel2tet
{
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

void Dream3DDataReader::LoadFile(std::string FileName)
{
    LOG("Open file %s\n", FileName.c_str());
    H5::H5File *file = NULL;

    try {
        file = new H5::H5File (FileName, H5F_ACC_RDONLY);
    } catch (...) {
        STATUS("Cound not open input file %s\n", FileName.c_str());
        exit(-1);
    }
    H5::Group DataContainers;

    DataContainers = H5::Group(file->openGroup("DataContainers"));
    VoxelDataContainer = new H5::Group(DataContainers.openGroup(this->DataContainer));


    // Load dimensional/spatial data
    VoxelDataContainer->openGroup("_SIMPL_GEOMETRY").openDataSet("SPACING").read(spacing_data,
                                                                                 H5::PredType::NATIVE_DOUBLE);
    VoxelDataContainer->openGroup("_SIMPL_GEOMETRY").openDataSet("ORIGIN").read(origin_data,
                                                                                H5::PredType::NATIVE_DOUBLE);
    VoxelDataContainer->openGroup("_SIMPL_GEOMETRY").openDataSet("DIMENSIONS").read(dimensions_data,
                                                                                    H5::PredType::NATIVE_INT);

    for (int i = 0; i < 3; i++) {
        this->BoundingBox.minvalues[i] = this->origin_data[i];
    }
    for (int i = 0; i < 3; i++) {
        this->BoundingBox.maxvalues[i] = this->origin_data[i] + this->dimensions_data[i] * this->spacing_data[i];
    }

    // Load voxel data.
    H5::DataSet GrainIds = VoxelDataContainer->openGroup("CellData").openDataSet(this->MaterialGroup);
    H5::DataSpace space = GrainIds.getSpace();

    int Ndims = space.getSimpleExtentNdims();
    hsize_t dims[Ndims];

    space.getSimpleExtentDims(dims);

    int DataLength = 1;
    for (int i = 0; i < Ndims; i++) {
        DataLength = DataLength * dims[i];
    }


    GrainIdsData = (int *) malloc(sizeof(int) * DataLength);
    GrainIds.read(GrainIdsData, H5::PredType::NATIVE_INT);
}
}
