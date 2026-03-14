#ifndef TETPIPELINE_H
#define TETPIPELINE_H

#include <array>
#include <memory>
#include <vector>

namespace voxel2tet
{
class Importer;
class MeshManipulations;
class Volume;

class TetPipeline
{
public:
    TetPipeline(Importer *Imp, std::unique_ptr<MeshManipulations> *Mesh, std::vector<Volume *> *Volumes, double eps);

    void Tetrahedralize() const;

private:
    Volume *FindVolumeContainingPoint(std::array<double, 3> P) const;

    void UpdateFixed() const;

    Importer *Imp;
    std::unique_ptr<MeshManipulations> *Mesh;
    std::vector<Volume *> *Volumes;
    double eps;
};
}

#endif // TETPIPELINE_H
