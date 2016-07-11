#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <vector>
#include <iostream>
#include <string>
#include <map>

namespace voxel2tet {
typedef std :: map< std :: string, std :: string >ValueMap;

/**
 * @brief The Options class supplies functionality for parsing, storing and accessing command line option.
 *
 * Functions for using defult values are provided. These are stored in the DefaultMap object. If the default values (DefaultMap)
 * are also given in the OptionMap object, the OptionMap has precedence.
 *
 * The user can also specify which field are required.
 */
class Options
{
private:
    ValueMap OptionMap;
    ValueMap DefaultMap;
    std :: vector< std :: string >RequiredFields;
public:

    /**
     * @brief Constructor.
     * @param argc Number of charachters in command line
     * @param argv Array of chars in command line.
     * @param DefaultValues Mapping object for storing default values.
     * @param RequiredFields
     */
    Options(int argc, char *argv[], ValueMap DefaultValues, std :: vector< std :: string >RequiredFields);

    /**
     * @brief Add a default value for a key.
     * @param keyname Name of key
     * @param value Default value of key
     */
    void AddDefaultMap(std :: string keyname, std :: string value);

    /**
     * \copydoc Options::AddDefaultMap
     */
    void AddDefaultMap(std :: string keyname, double value);

    /**
     * @brief Add a required key.
     * @param keyname Name of required key
     */
    void AddRequiredKey(std :: string keyname);

    /**
     * @brief Checks if all required keys are specified. If not, an exception is raised.
     */
    void CheckRequiredFields();

    /**
     * @brief Determines if a give key is specified in OptionsMap. Note that default values are not checked.
     * @param keyname Name of key.
     * @return True if key exists, false otherwise
     */
    bool has_key(std :: string keyname);

    /**
     * @brief Give string value of key
     * @param keyname Name of key
     * @return String value of key.
     */
    std :: string GiveStringValue(std :: string keyname);

    /**
     * @brief Give integer value of key
     * @param keyname Name of key
     * @return Integer value of key
     */
    int GiveIntegerValue(std :: string keyname);

    /**
     * @brief Give double value of key
     * @param keyname Name of key
     * @return Double value of key
     */
    double GiveDoubleValue(std :: string keyname);

    /**
     * @brief Give boolean value of key
     * @param keyname Name of key
     * @return Boolean value of key
     */
    bool GiveBooleanValue(std :: string keyname);

    /**
     * @brief Give list of integer corresponding to key
     * @param keyname Name of key
     * @return List of integers
     */
    std::vector<int> GiveIntegerList(std :: string keyname);
};
}

#endif /* OPTIONS_H_ */
