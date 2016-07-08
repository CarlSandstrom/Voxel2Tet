#include <stdio.h>

#include "TriTriIntersect.h"

int main(int, char *[])
{
    double V0 [ 3 ] = {
        0, 0, 0
    };
    double V1 [ 3 ] = {
        1, 0, 0
    };
    double V2 [ 3 ] = {
        0, 1, 0
    };

    double P [ 3 ] = {
        0.5, .5, 0.01
    };

    point_in_tri(V0, V1, V2, P);
}
