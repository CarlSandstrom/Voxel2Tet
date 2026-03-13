#ifndef SMOOTHINGPIPELINE_H
#define SMOOTHINGPIPELINE_H

#include <functional>

namespace voxel2tet
{
class SmoothingPipeline
{
public:
    typedef std::function<void(void)> StepDelegate;

    SmoothingPipeline(StepDelegate SmoothEdgesStep, StepDelegate SmoothSurfacesStep);

    void SmoothEdges() const;

    void SmoothSurfaces() const;

private:
    StepDelegate SmoothEdgesStep;
    StepDelegate SmoothSurfacesStep;
};
}

#endif // SMOOTHINGPIPELINE_H
