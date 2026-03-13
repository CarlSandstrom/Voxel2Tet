#include "ExportPipeline.h"

namespace voxel2tet
{
ExportPipeline::ExportPipeline(StepDelegate ExportSurfacesStep, StepDelegate ExportVolumesStep,
                               StepDelegate ExportPhasesStep)
{
    this->ExportSurfacesStep = ExportSurfacesStep;
    this->ExportVolumesStep = ExportVolumesStep;
    this->ExportPhasesStep = ExportPhasesStep;
}

void ExportPipeline::ExportSurfaces() const
{
    this->ExportSurfacesStep();
}

void ExportPipeline::ExportVolumes() const
{
    this->ExportVolumesStep();
}

void ExportPipeline::ExportPhases() const
{
    this->ExportPhasesStep();
}
}
