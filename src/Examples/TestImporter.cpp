#include <iostream>
#include <string>
#include <stdio.h>

#include "hdf5DataReader.h"

int main( int argc, char *argv[] )
{
    voxel2tet :: hdf5DataReader DataReader = voxel2tet :: hdf5DataReader();
    DataReader.LoadFile("/home/carl/dev/Voxel2Tet/ExampleInput/substructure.hdf5");
    double dt[3];
    DataReader.GiveSpacing(dt);
    for (int i=0; i<3; i++ ) {
        printf("%f\n", dt[i]);
        dt[i] = -9.1;
    }

}
