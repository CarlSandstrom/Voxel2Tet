#ifndef MISCFUNCTIONS_H_
#define MISCFUNCTIONS_H_

#define LOGOUTPUT 1
#define STATOUTPUT 1

#define log(format, args...) dolog (__FUNCTION__, format, args)
#define outputstat(format, args...) dooutputstat (format, args)

namespace voxel2tet
{

void dolog(const char *functionname, const char *fmt, ...);
void dooutputstat(const char *fmt, ...);

}

#endif /* MISCFUNCTIONS_H_ */
