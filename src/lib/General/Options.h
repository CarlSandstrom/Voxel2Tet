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
    std::vector<std::string> RequiredFields;
public:
    Options( int argc, char *argv[], ValueMap DefaultValues, std::vector<std::string> RequiredFields );

    void AddDefaultMap(std::string keyname, std::string value);
    void AddRequiredKey(std::string keyname);

    void CheckRequiredFields();

    bool has_key(std::string keyname);

    std::string GiveStringValue(std::string keyname);
    int GiveIntegerValue(std::string keyname);
    double GiveDoubleValue(std::string keyname);
    bool GiveBooleanValue(std::string keyname);

    int GiveIntegerList(std::string keyname);
};

}

#endif /* OPTIONS_H_ */
