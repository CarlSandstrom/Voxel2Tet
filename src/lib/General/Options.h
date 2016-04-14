#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <vector>
#include <iostream>
#include <string>
#include <map>

namespace voxel2tet{

typedef std::map <std::string, std::string> ValueMap;

class Options
{
private:
    ValueMap OptionMap;
    ValueMap DefaultMap;
public:
    Options( int argc, char *argv[], ValueMap DefaultValues );

    void AddDefaultMap(std::string keyname, std::string value);

    bool has_key(std::string keyname);

    std::string GiveStringValue(std::string keyname);
    int GiveIntegerValue(std::string keyname);
    double GiveDoubleValue(std::string keyname);
    bool GiveBooleanValue(std::string keyname);

    int GiveIntegerList(std::string keyname);
};

}

#endif /* OPTIONS_H_ */
