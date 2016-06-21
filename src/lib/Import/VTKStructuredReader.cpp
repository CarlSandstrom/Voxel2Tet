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

    if (!Input.is_open()) {
        LOG("Cound not open input file %s\n", FileName);
        exit(-1);
    }

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
                        this->dimensions_data[0] = std::stoi (Strings[1])-1;
                        this->dimensions_data[1] = std::stoi (Strings[2])-1;
                        this->dimensions_data[2] = std::stoi (Strings[3])-1;
                    } else if (strcasecmp(Strings[0].c_str(), "ORIGIN") == 0) {
                        this->origin_data[0] = std::stof (Strings[1]);
                        this->origin_data[1] = std::stof (Strings[2]);
                        this->origin_data[2] = std::stof (Strings[3]);
                    } else if (strcasecmp(Strings[0].c_str(), "SPACING") == 0) {
                        this->spacing_data[0] = std::stof (Strings[1]);
                        this->spacing_data[1] = std::stof (Strings[2]);
                        this->spacing_data[2] = std::stof (Strings[3]);
                        for (int i=0; i<3; i++) this->BoundingBox.minvalues[i] = this->origin_data[i];
                        for (int i=0; i<3; i++) this->BoundingBox.maxvalues[i] = this->origin_data[i]+this->dimensions_data[i]*this->spacing_data[i];
                    } else if (strcasecmp(Strings[0].c_str(), "CELL_DATA") == 0) {
                        this->celldata = std::stoi (Strings[1]);
                        this->GrainIdsData = (int*) malloc(sizeof(int)*this->celldata);
                    } else if (strcasecmp(Strings[0].c_str(), "SCALARS") == 0) {
                        this->DataName = Strings[1];
                    } else if (strcasecmp(Strings[0].c_str(), "LOOKUP_TABLE") == 0) {
                        this->TableName = Strings[1];
                        int scount=0;
                        while (true) {
                            std::getline(Input, line);
                            std::vector<std::string> TableValues = SplitString(line, ' ');
                            for (std::string s: TableValues) {
                                if (s.length()>0) {
                                    int Value = std::stoi (s);
                                    this->GrainIdsData[scount] = Value;
                                    //printf("%u: %u\n", scount, Value);
                                    scount++;
                                }
                            }
                            if (scount==this->celldata) break;
                        }
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
    double coords[3]={x, y, z};
    int indices[3];
    for (int i=0; i<3; i++) {
        double intcoord = (coords[i]-this->origin_data[i])/this->spacing_data[i];
        indices[i] =  floor(intcoord);
    }
    return this->GiveMaterialIDByIndex(indices[0], indices[1], indices[2]);

}

int VTKStructuredReader::GiveMaterialIDByIndex(int xi, int yi, int zi)
{
    if (xi == -1) {
        return -1;
    } else if (xi == this->dimensions_data[0]) {
        return -2;
    } else if (yi == -1) {
        return -3;
    } else if (yi == this->dimensions_data[1]) {
        return -4;
    } else if (zi == -1) {
        return -5;
    } else if (zi == this->dimensions_data[2]) {
        return -6;
    }

    int index = zi*this->dimensions_data[1]*this->dimensions_data[0] + yi*this->dimensions_data[0] + xi;
    return this->GrainIdsData[index];

}

void VTKStructuredReader::GiveSpacing(double spacing[3])
{
    for (int i=0; i<3; i++) {spacing[i] = this->spacing_data[i];}
}

BoundingBoxType VTKStructuredReader::GiveBoundingBox()
{
    return this->BoundingBox;
}

void VTKStructuredReader::GiveDimensions(int dimensions[3])
{
    for (int i=0; i<3; i++) {dimensions[i] = this->dimensions_data[i];}
}

void VTKStructuredReader::GiveCoordinateByIndices(int xi, int yi, int zi, DoubleTriplet Coordinate)
{
    Coordinate.c[0] = xi*this->spacing_data[0] + this->origin_data[0];
    Coordinate.c[1] = yi*this->spacing_data[1] + this->origin_data[1];
    Coordinate.c[2] = zi*this->spacing_data[2] + this->origin_data[2];

}

void VTKStructuredReader::GiveOrigin(double origin[3])
{
    for (int i=0; i<3; i++) {origin[i] = this->origin_data[i];}
}

}
