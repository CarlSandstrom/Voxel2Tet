#include <iostream>
#include <string>
#include <stdio.h>

#include "Options.h"

int main( int argc, char *argv[] )
{
    std::map <std::string, std::string> DefaultOptions;
    DefaultOptions["IntVal"] = "1999";
    voxel2tet::Options O(argc, argv, DefaultOptions);

    int v = O.GiveIntegerValue("IntVal");
    printf ("IntVal = %u\n", v);

    double d = O.GiveDoubleValue("IntVal");
    printf ("IntVal = %5.5f\n", d);


}


