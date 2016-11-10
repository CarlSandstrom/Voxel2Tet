#include "VTKStructuredReader.h"

namespace voxel2tet
{
VTKStructuredReader::VTKStructuredReader()
{}

void VTKStructuredReader::LoadFile(std::string FileName)
{
    std::ifstream Input;
    Input.open(FileName, std::ios::in);

    bool IsCellData = false;

    if (!Input.is_open()) {
        STATUS("Cound not open input file %s\n", FileName.c_str());
        exit(-1);
    }

    int linecount = 0;
    if (Input.is_open()) {
        std::string line;
        while (std::getline(Input, line)) {
            if (linecount == 0) {
                this->VersionInfo = line;
            } else if (linecount == 1) {
                this->Title = line;
            } else if (linecount == 2) {
                if (strcasecmp(line.c_str(), "ASCII") != 0) {
                    STATUS("Can only open VTK ASCII file right now\n", 0);
                }
            } else {
                std::vector<std::string> Strings = SplitString(line, ' ');
                if (Strings.size() > 0) {
                    if (strcasecmp(Strings[0].c_str(), "DATASET") == 0) {
                        if ((!strcasecmp(Strings[1].c_str(), "STRUCTURED_POINTS")) == 0) {
                            STATUS("Can only handle VTK files with dataset STRUCTURED_POINTS\n", 0);
                        }
                    } else if (strcasecmp(Strings[0].c_str(), "DIMENSIONS") == 0) {
                        this->dimensions_data[0] = std::stoi(Strings[1]) - 1;
                        this->dimensions_data[1] = std::stoi(Strings[2]) - 1;
                        this->dimensions_data[2] = std::stoi(Strings[3]) - 1;
                    } else if (strcasecmp(Strings[0].c_str(), "ORIGIN") == 0) {
                        this->origin_data[0] = std::stof(Strings[1]);
                        this->origin_data[1] = std::stof(Strings[2]);
                        this->origin_data[2] = std::stof(Strings[3]);
                    } else if (strcasecmp(Strings[0].c_str(), "SPACING") == 0) {
                        this->spacing_data[0] = std::stof(Strings[1]);
                        this->spacing_data[1] = std::stof(Strings[2]);
                        this->spacing_data[2] = std::stof(Strings[3]);
                        for (int i = 0; i < 3; i++) {
                            this->BoundingBox.minvalues[i] = this->origin_data[i];
                        }
                        for (int i = 0; i < 3; i++) {
                            this->BoundingBox.maxvalues[i] =
                                    this->origin_data[i] + this->dimensions_data[i] * this->spacing_data[i];
                        }
                    } else if (strcasecmp(Strings[0].c_str(), "CELL_DATA") == 0) {
                        IsCellData = true;
                        this->celldata = std::stoi(Strings[1]);
                        this->GrainIdsData = (int *) malloc(sizeof(int) * this->celldata);
                    } else if (strcasecmp(Strings[0].c_str(), "POINT_DATA") == 0) {
                        IsCellData = false;
                        this->celldata = std::stoi(Strings[1]);
                        this->GrainIdsData = (int *) malloc(sizeof(int) * this->celldata);
                    } else if (strcasecmp(Strings[0].c_str(), "SCALARS") == 0) {
                        this->DataName = Strings[1];
                    } else if (strcasecmp(Strings[0].c_str(), "LOOKUP_TABLE") == 0) {
                        this->TableName = Strings[1];
                        int scount = 0;
                        while (true) {
                            std::getline(Input, line);
                            std::vector<std::string> TableValues = SplitString(line, ' ');
                            for (std::string s : TableValues) {
                                if (s.length() > 0) {
                                    int Value = std::stoi(s);
                                    this->GrainIdsData[scount] = Value;
                                    //printf("%u: %u\n", scount, Value);
                                    scount++;
                                }
                            }
                            if (scount == this->celldata) {
                                break;
                            }
                        }
                    } else {
                        STATUS("Token %s not recognized\n", Strings[0].c_str());
                    }
                }
            }
            linecount++;
        }
    }

    // If we are given points, move origin and change dimensions
    if (!IsCellData) {
        for (int i=0; i<3; i++) {
            this->dimensions_data[i] = this->dimensions_data[i] + 1;
            this->origin_data[i] = this->origin_data[i] - this->spacing_data[i]/2.0;
        }
    }

}
}
