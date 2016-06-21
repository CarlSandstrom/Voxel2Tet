#ifndef VTKSTRUCTUREDREADER_H
#define VTKSTRUCTUREDREADER_H

#include <iostream>
#include <fstream>
#include <string.h>

#include "Importer.h"
#include "MiscFunctions.h"

namespace voxel2tet
{

class VTKStructuredReader : public Importer
{
private:
    std::string VersionInfo;
    std::string Title;
    std::string DataName;
    std::string TableName;
    int celldata;
public:
    VTKStructuredReader();
    void LoadFile(std::string FileName);
};

}

#endif // VTKSTRUCTUREDREADER_H
