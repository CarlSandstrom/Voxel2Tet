#include "InputLoaderPipeline.h"

namespace voxel2tet
{
InputLoaderPipeline::InputLoaderPipeline(LoadFileDelegate LoadFileStep, LoadCallbackDelegate LoadCallbackStep)
{
    this->LoadFileStep = LoadFileStep;
    this->LoadCallbackStep = LoadCallbackStep;
}

void InputLoaderPipeline::LoadFile(const std::string &FileName) const
{
    this->LoadFileStep(FileName);
}

void InputLoaderPipeline::LoadCallback(cbMaterialIDByCoordinate MaterialIDByCoordinate, std::array<double, 3> origin,
                                       std::array<double, 3> spacing, std::array<int, 3> dimensions) const
{
    this->LoadCallbackStep(MaterialIDByCoordinate, origin, spacing, dimensions);
}
}
