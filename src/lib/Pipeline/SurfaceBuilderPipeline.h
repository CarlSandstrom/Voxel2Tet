#ifndef SURFACEBUILDERPIPELINE_H
#define SURFACEBUILDERPIPELINE_H

#include <functional>

namespace voxel2tet
{
class SurfaceBuilderPipeline
{
public:
    typedef std::function<void(void)> StepDelegate;

    SurfaceBuilderPipeline(StepDelegate FindSurfacesStep, StepDelegate FindEdgesStep);

    void FindSurfaces() const;

    void FindEdges() const;

private:
    StepDelegate FindSurfacesStep;
    StepDelegate FindEdgesStep;
};
}

#endif // SURFACEBUILDERPIPELINE_H
