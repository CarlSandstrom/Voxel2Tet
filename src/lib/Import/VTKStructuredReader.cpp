#include "VTKStructuredReader.h"

namespace voxel2tet
{


VTKStructuredReader::VTKStructuredReader()
{

}

void VTKStructuredReader::LoadFile(std::string FileName)
{

    std::ifstream Input;
    Input.open(FileName, std::ios::in);

    int linecount=0;
    if (Input.is_open()) {
        std::string line;
        while(std::getline(Input, line)) {
            if(linecount==0) {
                this->VersionInfo = line;
            } else if(linecount==1) {
                this->Title = line;
            } else if(linecount==2) {
                if (strcasecmp(line.c_str(), "ASCII")!=0) {
                    STATUS("Can only open VTK ASCII file right now\n", 0);
                }
            } else {
                std::vector<std::string> Strings = SplitString(line, ' ');
                if (Strings.size()>0) {
                    if (strcasecmp(Strings[0].c_str(), "DATASET") == 0) {
                        if (!strcasecmp(Strings[1].c_str(), "STRUCTURED_POINTS") == 0) {
                            STATUS ("Can only handle VTK files with dataset STRUCTURED_POINTS\n", 0);
                        }
                    } else if (strcasecmp(Strings[0].c_str(), "DIMENSIONS") == 0) {
                        this->dimensions_data[0] = std::stoi (Strings[1]);
                        this->dimensions_data[1] = std::stoi (Strings[2]);
                        this->dimensions_data[2] = std::stoi (Strings[3]);
                    } else if (strcasecmp(Strings[0].c_str(), "ORIGIN") == 0) {
                        this->origin_data[0] = std::stoi (Strings[1]);
                        this->origin_data[1] = std::stoi (Strings[2]);
                        this->origin_data[2] = std::stoi (Strings[3]);
                    } else if (strcasecmp(Strings[0].c_str(), "SPACING") == 0) {
                        this->spacing_data[0] = std::stoi (Strings[1]);
                        this->spacing_data[1] = std::stoi (Strings[2]);
                        this->spacing_data[2] = std::stoi (Strings[3]);
                    } else if (strcasecmp(Strings[0].c_str(), "CELL_DATA") == 0) {
                        this->celldata = std::stoi (Strings[1]);
                    } else if (strcasecmp(Strings[0].c_str(), "SCALARS") == 0) {

                    } else if (strcasecmp(Strings[0].c_str(), "LOOKUP_TABLE") == 0) {

                    } else {
                        STATUS("Token %s not recognized\n", Strings[0]);
                    }
                }
            }
            linecount++;
        }
    }

}

int VTKStructuredReader::GiveMaterialIDByCoordinate(double x, double y, double z)
{

}

int VTKStructuredReader::GiveMaterialIDByIndex(int xi, int yi, int zi)
{

}

void VTKStructuredReader::GiveSpacing(double spacing[3])
{

}

BoundingBoxType VTKStructuredReader::GiveBoundingBox()
{

}

void VTKStructuredReader::GiveDimensions(int dimensions[3])
{

}

void VTKStructuredReader::GiveCoordinateByIndices(int xi, int yi, int zi, DoubleTriplet Coordinate)
{

}

void VTKStructuredReader::GiveOrigin(double origin[3])
{

}

}
