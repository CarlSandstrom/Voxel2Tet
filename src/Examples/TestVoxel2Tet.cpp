#include <iostream>
#include <string>
#include <stdio.h>
#include <ctime>

#include "Options.h"
#include "Voxel2Tet.h"

int main( int argc, char *argv[] )
{

    std::map <std::string, std::string> DefaultOptions;
    voxel2tet::Options *Options = new voxel2tet::Options(argc, argv, DefaultOptions, {"input"});

    voxel2tet::Voxel2Tet v2t(Options);

    printf("Voxel2Tet - ");
    printf("Copyright Carl Sandström, 2015-2016\n\n");

    v2t.LoadFile(Options->GiveStringValue("input"));

    clock_t tstart = clock();
    v2t.Process();
    clock_t tstopprocess = clock();

    v2t.Tetrahedralize();
    clock_t tstoptet = clock();

    v2t.ExportAll();

    double ProcessTime = double(tstopprocess-tstart) / CLOCKS_PER_SEC;
    double TetTime = double(tstoptet-tstopprocess) / CLOCKS_PER_SEC;

    printf("Time for mesh process: %f\n", ProcessTime);
    printf("Time for tetrahedralization: %f\n", TetTime);

}
