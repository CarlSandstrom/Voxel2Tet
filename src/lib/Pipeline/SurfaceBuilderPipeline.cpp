#include "SurfaceBuilderPipeline.h"

namespace voxel2tet
{
SurfaceBuilderPipeline::SurfaceBuilderPipeline(StepDelegate FindSurfacesStep, StepDelegate FindEdgesStep)
{
    this->FindSurfacesStep = FindSurfacesStep;
    this->FindEdgesStep = FindEdgesStep;
}

void SurfaceBuilderPipeline::FindSurfaces() const
{
    this->FindSurfacesStep();
}

void SurfaceBuilderPipeline::FindEdges() const
{
    this->FindEdgesStep();
}
}
