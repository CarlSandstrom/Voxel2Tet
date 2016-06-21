#include "Importer.h"

namespace voxel2tet {

void Importer::GiveDimensions(int dimensions[3])
{
    for (int i=0; i<3; i++) {dimensions[i] = this->dimensions_data[i];}
}

void Importer::GiveSpacing(double spacing[3])
{
    for (int i=0; i<3; i++) {spacing[i] = this->spacing_data[i];}
}

void Importer::GiveOrigin(double origin[3])
{
    for (int i=0; i<3; i++) {origin[i] = this->origin_data[i];}
}

BoundingBoxType Importer::GiveBoundingBox()
{
    return this->BoundingBox;
}

int Importer::GiveMaterialIDByIndex(int xi, int yi, int zi)
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

int Importer::GiveMaterialIDByCoordinate(double x, double y, double z)
{
    double coords[3]={x, y, z};
    int indices[3];
    for (int i=0; i<3; i++) {
        double intcoord = (coords[i]-this->origin_data[i])/this->spacing_data[i];
        indices[i] = floor(intcoord);
    }
    return this->GiveMaterialIDByIndex(indices[0], indices[1], indices[2]);

}

void Importer::GiveCoordinateByIndices(int xi, int yi, int zi, DoubleTriplet Coordinate)
{
    Coordinate.c[0] = xi*this->spacing_data[0] + this->origin_data[0];
    Coordinate.c[1] = yi*this->spacing_data[1] + this->origin_data[1];
    Coordinate.c[2] = zi*this->spacing_data[2] + this->origin_data[2];
}


}

