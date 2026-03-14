#ifndef EXPORTPIPELINE_H
#define EXPORTPIPELINE_H

#include <memory>
#include <string>
#include <vector>

#include "Exporter.h"

namespace voxel2tet
{
class MeshData;
class MeshManipulations;
class Options;
class Volume;

class ExportPipeline
{
public:
    ExportPipeline(Options *Opt, MeshManipulations *Mesh, std::vector<Volume *> *Volumes);

    void SetMesh(MeshManipulations *Mesh);

    void SetVolumes(std::vector<Volume *> *Volumes);

    void ExportSurface(std::string FileName, Exporter_FileTypes FileType, int VolumeID = -1) const;

    void ExportVolume(std::string FileName, Exporter_FileTypes FileType) const;

    void ExportPhases(std::string FileName, Exporter_FileTypes FileType) const;

    void ExportSurfaces() const;

    void ExportVolumes() const;

    void ExportPhases() const;

    static void ExportSurface(MeshData &Mesh, std::string FileName, Exporter_FileTypes FileType, int VolumeID = -1);

    static void ExportVolume(MeshData &Mesh, std::string FileName, Exporter_FileTypes FileType);

private:
    static std::unique_ptr<Exporter> CreateSurfaceExporter(MeshData &Mesh, Exporter_FileTypes FileType);

    static std::unique_ptr<Exporter> CreateVolumeExporter(MeshData &Mesh, Exporter_FileTypes FileType);

    Options *Opt;
    MeshManipulations *Mesh;
    std::vector<Volume *> *Volumes;
};
}

#endif // EXPORTPIPELINE_H
