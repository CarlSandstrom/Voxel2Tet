#include <iostream>
#include <ctime>
#include <exception>

#include "Options.h"
#include "Voxel2Tet.h"

int main(int argc, char *argv[])
{
    try {
        std::map<std::string, std::string> DefaultOptions;
        voxel2tet::Options Options(argc, argv, DefaultOptions, {});

        printf("Voxel2Tet - ");
        printf("Copyright Carl Sandström, 2015-2016\n");

        voxel2tet::Voxel2TetClass v2t(&Options);

        if (Options.has_key("help") | (argc == 1)) {
            v2t.PrintHelp();
            return 0;
        }


        v2t.LoadFile(Options.GiveStringValue("input"));


        clock_t tstart = clock();
        v2t.Process();
        clock_t tstopprocess = clock();

        v2t.ExportAllSurfaces();
        v2t.ExportAllPhases();

        v2t.Tetrahedralize();

        v2t.ExportAllVolumes();
        v2t.ExportStatistics();

        clock_t tstoptet = clock();


        double ProcessTime = double(tstopprocess - tstart) / CLOCKS_PER_SEC;
        double TetTime = double(tstoptet - tstopprocess) / CLOCKS_PER_SEC;

        printf("Time for mesh process: %f\n", ProcessTime);
        printf("Time for tetrahedralization: %f\n", TetTime);
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
