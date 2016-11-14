#include <cstdlib>
#include <sstream>

#include "Options.h"
#include "MiscFunctions.h"

namespace voxel2tet
{
Options::Options(int argc, char *argv[], ValueMap DefaultValues, std::vector<std::string> RequiredFields)
{
    for (std::string s : RequiredFields) {
        this->RequiredFields.push_back(s);
    }

    // Parse command line arguments and simply store them in a map
    for (int i = 1; i < argc; i++) {
        if (char(argv[i][0]) == '-') {
            bool NextIsFlag = false;

            if (i == (argc - 1)) {
                NextIsFlag = true;
            } else if (char(argv[i + 1][0]) == '-') {
                NextIsFlag = true;
            }

            std::string FlagName;
            int j = 1;
            while (argv[i][j] != '\0') {
                FlagName += argv[i][j];
                j++;
            }

            std::string Value;
            if (!NextIsFlag) {
                Value = argv[i + 1];
            } else {
                Value = "1";
            }
            this->OptionMap[FlagName] = Value;

            //std :: cout << "FlagName: " << FlagName << ", Value:" << Value <<"\n";
        }
    }
    this->DefaultMap = DefaultValues;

    this->CheckRequiredFields();
}

void Options::CheckRequiredFields()
{
    // Check for required fields
    bool FieldsOk = true;
    for (std::string s : this->RequiredFields) {
        if (!this->has_key(s)) {
            STATUS("Input parameter \"-%s\" is required. Use \"-help\" for a list of avalible option\n", s.c_str());
            FieldsOk = false;
        }
    }
    if (!FieldsOk) {
        exit(-1);
    }
}

void Options::AddDefaultMap(std::string keyname, std::string value)
{
    this->DefaultMap[keyname] = value;
}

void Options::AddDefaultMap(std::string keyname, double value)
{
    this->AddDefaultMap(keyname, strfmt("%f", value));
}

void Options::AddRequiredKey(std::string keyname)
{
    this->RequiredFields.push_back(keyname);
    this->CheckRequiredFields();
}

void Options::SetKey(std::string keyname, std::string value)
{
    this->OptionMap[keyname] = value;
}

void Options::SetKey(std::string keyname, double value)
{
    this->OptionMap[keyname] = strfmt("%f", value);
}

bool Options::has_key(std::string keyname)
{
    if (this->OptionMap.find(keyname) == this->OptionMap.end()) {
        return false;
    }
    return true;
}

std::string Options::GiveStringValue(std::string keyname)
{
    std::string Value;
    if (this->has_key(keyname)) {
        Value = this->OptionMap[keyname];
    } else if (this->DefaultMap.find(keyname) != this->DefaultMap.end()) {
        Value = this->DefaultMap[keyname];
    }
    return Value;
}

int Options::GiveIntegerValue(std::string keyname)
{
    std::string Value = this->GiveStringValue(keyname);
    int iValue = std::atoi(Value.c_str());
    return iValue;
}

double Options::GiveDoubleValue(std::string keyname)
{
    std::string Value = this->GiveStringValue(keyname);
    double dValue = std::atof(Value.c_str());
    return dValue;
}

bool Options::GiveBooleanValue(std::string keyname)
{
    std::string Value = this->GiveStringValue(keyname);
    double iValue = std::atoi(Value.c_str());
    return bool(iValue);
}

std::vector<int> Options::GiveIntegerList(std::string keyname)
{
    // Not implemented
    std::string List = this->GiveStringValue(keyname);
    std::vector<int> IntList;
    if ((List[0] == '[') & (List[List.size() - 1] == ']')) {

        std::stringstream SubStringStream(List.substr(1, List.size() - 1));
        std::string StrItem;

        while (std::getline(SubStringStream, StrItem, ' ')) {
            IntList.push_back(std::atoi(StrItem.c_str()));
        }
    } else {
        STATUS("Argument for switch '%s' must be of type list (list enclosed in brackets [])\n", keyname.c_str());
    }
    return IntList;
}
}
