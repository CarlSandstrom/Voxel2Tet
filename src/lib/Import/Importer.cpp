#include "Importer.h"

namespace voxel2tet {
void Importer :: GiveDimensions(int dimensions [ 3 ])
{
    for ( int i = 0; i < 3; i++ ) {
        if (this->UseCutOut) {
            dimensions [ i ] = this->CutOut.maxvalues[i]-this->CutOut.minvalues[i]+1;
        } else {
            dimensions [ i ] = this->dimensions_data[i];
        }
    }
}

void Importer :: GiveSpacing(double spacing [ 3 ])
{
    for ( int i = 0; i < 3; i++ ) {
        spacing [ i ] = this->spacing_data [ i ];
    }
}

void Importer :: GiveOrigin(double origin [ 3 ])
{
    for ( int i = 0; i < 3; i++ ) {
        origin [ i ] = this->origin_data [ i ];
    }
    if (UseCutOut) {
        for (int i=0; i<3; i++) origin[i]=this->origin_data[i]+this->CutOut.minvalues[i]*spacing_data[i];
    }
}

BoundingBoxType Importer :: GiveBoundingBox()
{
    if (this->UseCutOut) {
        IntTriplet PsuedoDimensions;
        this->GiveDimensions(PsuedoDimensions);
        BoundingBoxType bb;
        for (int i=0; i<3; i++) bb.minvalues[i] = (this->CutOut.minvalues[i]-1)*spacing_data[i];
        for (int i=0; i<3; i++) bb.maxvalues[i] = (this->CutOut.maxvalues[i]+1)*spacing_data[i];
        return bb;
    } else {
        return this->BoundingBox;
    }
}

int Importer :: GiveMaterialIDByIndex(int xi, int yi, int zi)
{

    IntTriplet PseudoDimensions;
    this->GiveDimensions(PseudoDimensions);


    if ( xi == -1 ) {
        return -1;
    } else if ( xi >= PseudoDimensions [ 0 ] ) {
        return -2;
    } else if ( yi == -1 ) {
        return -3;
    } else if ( yi == PseudoDimensions [ 1 ] ) {
        return -4;
    } else if ( zi == -1 ) {
        return -5;
    } else if ( zi == PseudoDimensions [ 2 ] ) {
        return -6;
    }

    // If using Cut out, xi, yi and zi has to be shifted
    if (this->UseCutOut) {
        xi=xi+this->CutOut.minvalues[0];
        yi=yi+this->CutOut.minvalues[1];
        zi=zi+this->CutOut.minvalues[2];
    }

    int index = zi * this->dimensions_data [ 1 ] * this->dimensions_data [ 0 ] + yi * this->dimensions_data [ 0 ] + xi;
    return this->GrainIdsData [ index ];
}

int Importer :: GiveMaterialIDByCoordinate(double x, double y, double z)
{
    double coords [ 3 ] = {
        x, y, z
    };
    int indices [ 3 ];
    for ( int i = 0; i < 3; i++ ) {
        double intcoord = ( coords [ i ] - this->origin_data [ i ] ) / this->spacing_data [ i ];
        indices [ i ] = floor(intcoord);
    }
    return this->GiveMaterialIDByIndex(indices [ 0 ], indices [ 1 ], indices [ 2 ]);
}

void Importer :: GiveIndicesByCoordinate(double x, double y, double z, IntTriplet Indices)
{
    Indices[0] = floor( (x-this->origin_data[0]) / this->spacing_data[0] );
    Indices[1] = floor( (y-this->origin_data[1]) / this->spacing_data[1] );
    Indices[2] = floor( (z-this->origin_data[2]) / this->spacing_data[2] );
}

void Importer :: GiveCoordinateByIndices(int xi, int yi, int zi, DoubleTriplet Coordinate)
{
    Coordinate.c [ 0 ] = xi * this->spacing_data [ 0 ] + this->origin_data [ 0 ];
    Coordinate.c [ 1 ] = yi * this->spacing_data [ 1 ] + this->origin_data [ 1 ];
    Coordinate.c [ 2 ] = zi * this->spacing_data [ 2 ] + this->origin_data [ 2 ];
}
}
