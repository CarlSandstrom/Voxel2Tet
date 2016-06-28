#include "CallbackImporter.h"

namespace  voxel2tet {
CallbackImporter :: CallbackImporter(cbMaterialIDByCoordinate MaterialIdByCoordinateCallback, std :: array< double, 3 >Origin, std :: array< double, 3 >Spacing, std :: array< int, 3 >Dimensions)
{
    this->MaterialByCoordinate = MaterialIdByCoordinateCallback;
    this->spacing_data = Spacing;
    this->origin_data = Origin;
    this->dimensions_data = Dimensions;

    this->BoundingBox.minvalues = Origin;
    for ( int i = 0; i < 3; i++ ) {
        this->BoundingBox.maxvalues [ i ] = Origin [ i ] + Spacing [ i ] * Dimensions [ i ];
    }
}

int CallbackImporter :: GiveMaterialIDByCoordinate(double x, double y, double z)
{
    if ( x < ( this->BoundingBox.minvalues [ 0 ] - eps ) ) {
        return -1;
    } else if ( x > ( this->BoundingBox.maxvalues [ 0 ] + eps ) ) {
        return -2;
    } else if ( y < ( this->BoundingBox.minvalues [ 1 ] - eps ) ) {
        return -3;
    } else if ( y > ( this->BoundingBox.maxvalues [ 1 ] + eps ) ) {
        return -4;
    } else if ( z < ( this->BoundingBox.minvalues [ 2 ] - eps ) ) {
        return -5;
    } else if ( z > ( this->BoundingBox.maxvalues [ 2 ] + eps ) ) {
        return -6;
    }

    return this->MaterialByCoordinate(x, y, z);
}

int CallbackImporter :: GiveMaterialIDByIndex(int xi, int yi, int zi)
{
    double x = this->origin_data [ 0 ] + this->spacing_data [ 0 ] * ( xi + .5 );
    double y = this->origin_data [ 1 ] + this->spacing_data [ 1 ] * ( yi + .5 );
    double z = this->origin_data [ 2 ] + this->spacing_data [ 2 ] * ( zi + .5 );
    return this->GiveMaterialIDByCoordinate(x, y, z);
}

void CallbackImporter :: GiveSpacing(double spacing [ 3 ])
{
    for ( int i = 0; i < 3; i++ ) {
        spacing [ i ] = this->spacing_data [ i ];
    }
}

BoundingBoxType CallbackImporter :: GiveBoundingBox()
{
    BoundingBoxType BoundingBox;
    for ( int i = 0; i < 3; i++ ) {
        BoundingBox.maxvalues [ i ] = this->BoundingBox.maxvalues [ i ];
        BoundingBox.minvalues [ i ] = this->BoundingBox.minvalues [ i ];
    }
    ;
    return BoundingBox;
}

void CallbackImporter :: GiveDimensions(int dimensions [ 3 ])
{
    for ( int i = 0; i < 3; i++ ) {
        dimensions [ i ] = this->dimensions_data [ i ];
    }
}

void CallbackImporter :: GiveCoordinateByIndices(int xi, int yi, int zi, DoubleTriplet Coordinate)
{
    Coordinate.c [ 0 ] = xi * this->spacing_data [ 0 ] + this->origin_data [ 0 ];
    Coordinate.c [ 1 ] = yi * this->spacing_data [ 1 ] + this->origin_data [ 1 ];
    Coordinate.c [ 2 ] = zi * this->spacing_data [ 2 ] + this->origin_data [ 2 ];
}

void CallbackImporter :: GiveOrigin(double origin [ 3 ])
{
    for ( int i = 0; i < 3; i++ ) {
        origin [ i ] = this->origin_data [ i ];
    }
}
}
