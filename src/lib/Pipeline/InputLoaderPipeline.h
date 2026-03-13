#ifndef INPUTLOADERPIPELINE_H
#define INPUTLOADERPIPELINE_H

#include <array>
#include <functional>
#include <string>

#include "CallbackImporter.h"

namespace voxel2tet
{
class InputLoaderPipeline
{
public:
    typedef std::function<void(const std::string &)> LoadFileDelegate;
    typedef std::function<void(cbMaterialIDByCoordinate, std::array<double, 3>,
                               std::array<double, 3>, std::array<int, 3>)> LoadCallbackDelegate;

    InputLoaderPipeline(LoadFileDelegate LoadFileStep, LoadCallbackDelegate LoadCallbackStep);

    void LoadFile(const std::string &FileName) const;

    void LoadCallback(cbMaterialIDByCoordinate MaterialIDByCoordinate, std::array<double, 3> origin,
                      std::array<double, 3> spacing, std::array<int, 3> dimensions) const;

private:
    LoadFileDelegate LoadFileStep;
    LoadCallbackDelegate LoadCallbackStep;
};
}

#endif // INPUTLOADERPIPELINE_H
