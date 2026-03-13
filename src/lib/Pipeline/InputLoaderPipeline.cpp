#include "InputLoaderPipeline.h"

#include <cstring>
#include <strings.h>
#include <vector>

#include "Dream3DDataReader.h"
#include "VTKStructuredReader.h"

namespace voxel2tet
{
InputLoaderPipeline::InputLoaderPipeline(Options *Opt)
{
    this->Opt = Opt;
}

Importer *InputLoaderPipeline::LoadFile(const std::string &FileName) const
{
    STATUS("Load file %s\n", FileName.c_str());

    // Find extension for input file
    char *ext;
    ext = strrchr((char *) FileName.c_str(), '.') + 1;

    Importer *Import;

    if ((strcasecmp(ext, "dream3d") == 0) | (strcasecmp(ext, "hdf5") == 0)) {
        Import = new Dream3DDataReader(this->Opt->GiveStringValue("DataContainer"),
                                       this->Opt->GiveStringValue("MaterialId"));
    } else if (strcasecmp(ext, "vtk") == 0) {
        Import = new VTKStructuredReader();
    } else {
        STATUS("File extension %s not recognized\n", ext);
        exit(-1);
    }

    Import->LoadFile(FileName);
    this->ApplyVoxelCutOut(Import);

    return Import;
}

Importer *InputLoaderPipeline::LoadCallback(cbMaterialIDByCoordinate MaterialIDByCoordinate,
                                            std::array<double, 3> origin, std::array<double, 3> spacing,
                                            std::array<int, 3> dimensions) const
{
    STATUS("Setup callback functions\n", 0);
    return new CallbackImporter(MaterialIDByCoordinate, origin, spacing, dimensions);
}

void InputLoaderPipeline::ApplyVoxelCutOut(Importer *Import) const
{
    // If specified, use cut-out argument
    if (this->Opt->has_key("voxelcutout")) {
        std::vector<int> VoxelCutOut = this->Opt->GiveIntegerList("voxelcutout");
        if (VoxelCutOut.size() != 6) {
            STATUS("VoxelCutOut should be specified as a bounding box consisting of integers for the indices", 0);
            exit(-1);
        }
        Import->UseCutOut = true;
        Import->CutOut.minvalues[0] = VoxelCutOut[0];
        Import->CutOut.minvalues[1] = VoxelCutOut[1];
        Import->CutOut.minvalues[2] = VoxelCutOut[2];
        Import->CutOut.maxvalues[0] = VoxelCutOut[3];
        Import->CutOut.maxvalues[1] = VoxelCutOut[4];
        Import->CutOut.maxvalues[2] = VoxelCutOut[5];
    }
}
}
