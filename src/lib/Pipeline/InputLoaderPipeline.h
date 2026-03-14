#ifndef INPUTLOADERPIPELINE_H
#define INPUTLOADERPIPELINE_H

#include <array>
#include <memory>
#include <string>

#include "CallbackImporter.h"
#include "Importer.h"
#include "Options.h"

namespace voxel2tet
{
class InputLoaderPipeline
{
public:
    explicit InputLoaderPipeline(Options *Opt);

    std::unique_ptr<Importer> LoadFile(const std::string &FileName) const;

    std::unique_ptr<Importer> LoadCallback(cbMaterialIDByCoordinate MaterialIDByCoordinate,
                                           std::array<double, 3> origin, std::array<double, 3> spacing,
                                           std::array<int, 3> dimensions) const;

private:
    void ApplyVoxelCutOut(Importer *Import) const;

    Options *Opt;
};
}

#endif // INPUTLOADERPIPELINE_H
