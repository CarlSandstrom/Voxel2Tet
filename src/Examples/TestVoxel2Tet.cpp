#include <iostream>
#include <string>
#include <stdio.h>
#include <ctime>

#include "Options.h"
#include "Voxel2Tet.h"

int main( int argc, char *argv[] )
{

    clock_t tstart=clock();

    std::map <std::string, std::string> DefaultOptions;
    voxel2tet::Options *Options = new voxel2tet::Options(argc, argv, DefaultOptions, {});

    printf("Voxel2Tet - ");
    printf("Copyright Carl Sandström, 2015-2016\n");

    voxel2tet::Voxel2TetClass v2t(Options);

    if (Options->has_key("help")) {
        v2t.PrintHelp();
        exit(0);
    }

    Options->AddRequiredKey("input");

    v2t.LoadFile(Options->GiveStringValue("input"));

    v2t.Process();

    v2t.Tetrahedralize();

    v2t.ExportAll();

    clock_t tstop=clock();

    double TotalTime = double(tstop-tstart) / CLOCKS_PER_SEC;

    printf("Total time: %fs\n", TotalTime);

}
