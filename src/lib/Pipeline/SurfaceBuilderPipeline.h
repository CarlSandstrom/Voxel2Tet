#ifndef SURFACEBUILDERPIPELINE_H
#define SURFACEBUILDERPIPELINE_H

#include <vector>

namespace voxel2tet
{
class Importer;
class MeshManipulations;
class Options;
class PhaseEdge;
class Smoother;
class Surface;
class VertexType;
class Volume;

class SurfaceBuilderPipeline
{
public:
    SurfaceBuilderPipeline(Importer *Imp, MeshManipulations *Mesh, Options *Opt, Smoother *SurfaceSmoother,
                           Smoother *EdgeSmoother, std::vector<Surface *> *Surfaces,
                           std::vector<Volume *> *Volumes, std::vector<PhaseEdge *> *PhaseEdges,
                           bool TreatZeroAsVoid);

    void FindSurfaces() const;

    void FindEdges() const;

private:
    void UpdateSurfaces() const;

    void AddSurfaceSquare(std::vector<int> VertexIDs, std::vector<int> phases, int normalphase) const;

    PhaseEdge *AddPhaseEdge(std::vector<VertexType *> EdgeSegment, std::vector<int> Phases) const;

    Importer *Imp;
    MeshManipulations *Mesh;
    Options *Opt;
    Smoother *SurfaceSmoother;
    Smoother *EdgeSmoother;
    std::vector<Surface *> *Surfaces;
    std::vector<Volume *> *Volumes;
    std::vector<PhaseEdge *> *PhaseEdges;
    bool TreatZeroAsVoid;
};
}

#endif // SURFACEBUILDERPIPELINE_H
