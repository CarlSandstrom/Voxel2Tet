#include <vector>
#include <stdio.h>
#include <algorithm>
#include <stdarg.h>

#include "MeshComponents.h"

namespace voxel2tet {

void dolog(const char *functionname, const char *fmt, ...)
{
#if LOGOUTPUT > -1
        va_list argp;
        va_start(argp, fmt);
        printf("%s:\t", functionname);
        vfprintf(stdout, fmt, argp);
        va_end(argp);
#endif
}

void dooutputstat(const char *fmt, ...)
{
        va_list argp;
        va_start(argp, fmt);
        vfprintf(stdout, fmt, argp);
        va_end(argp);
}

}
