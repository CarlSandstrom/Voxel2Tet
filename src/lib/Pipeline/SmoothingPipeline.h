#ifndef SMOOTHINGPIPELINE_H
#define SMOOTHINGPIPELINE_H

#include <vector>

namespace voxel2tet
{
class Importer;
class MeshManipulations;
class PhaseEdge;
class Smoother;
class Surface;

class SmoothingPipeline
{
public:
    SmoothingPipeline(Importer *Imp, MeshManipulations *Mesh, Smoother *EdgeSmoother,
                      std::vector<Surface *> *Surfaces, std::vector<PhaseEdge *> *PhaseEdges, double eps);

    void SmoothEdges() const;

    void SmoothSurfaces() const;

private:
    Importer *Imp;
    MeshManipulations *Mesh;
    Smoother *EdgeSmoother;
    std::vector<Surface *> *Surfaces;
    std::vector<PhaseEdge *> *PhaseEdges;
    double eps;
};
}

#endif // SMOOTHINGPIPELINE_H
