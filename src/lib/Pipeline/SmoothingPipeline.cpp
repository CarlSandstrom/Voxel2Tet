#include "SmoothingPipeline.h"

namespace voxel2tet
{
SmoothingPipeline::SmoothingPipeline(StepDelegate SmoothEdgesStep, StepDelegate SmoothSurfacesStep)
{
    this->SmoothEdgesStep = SmoothEdgesStep;
    this->SmoothSurfacesStep = SmoothSurfacesStep;
}

void SmoothingPipeline::SmoothEdges() const
{
    this->SmoothEdgesStep();
}

void SmoothingPipeline::SmoothSurfaces() const
{
    this->SmoothSurfacesStep();
}
}
