#ifndef EXPORTPIPELINE_H
#define EXPORTPIPELINE_H

#include <functional>

namespace voxel2tet
{
class ExportPipeline
{
public:
    typedef std::function<void(void)> StepDelegate;

    ExportPipeline(StepDelegate ExportSurfacesStep, StepDelegate ExportVolumesStep, StepDelegate ExportPhasesStep);

    void ExportSurfaces() const;

    void ExportVolumes() const;

    void ExportPhases() const;

private:
    StepDelegate ExportSurfacesStep;
    StepDelegate ExportVolumesStep;
    StepDelegate ExportPhasesStep;
};
}

#endif // EXPORTPIPELINE_H
