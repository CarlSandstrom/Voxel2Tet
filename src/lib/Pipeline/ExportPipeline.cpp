#include "ExportPipeline.h"

#include <sstream>

#include "AbaqusExporter.h"
#include "MeshData.h"
#include "MeshManipulations.h"
#include "MiscFunctions.h"
#include "OFFExporter.h"
#include "OOFEMExporter.h"
#include "Options.h"
#include "STLExporter.h"
#include "SimpleExporter.h"
#include "TetGenExporter.h"
#include "VTKExport.h"
#include "Volume.h"

namespace voxel2tet
{
ExportPipeline::ExportPipeline(Options *Opt, MeshManipulations *Mesh, std::vector<Volume *> *Volumes)
{
    this->Opt = Opt;
    this->Mesh = Mesh;
    this->Volumes = Volumes;
}

void ExportPipeline::SetMesh(MeshManipulations *Mesh)
{
    this->Mesh = Mesh;
}

void ExportPipeline::SetVolumes(std::vector<Volume *> *Volumes)
{
    this->Volumes = Volumes;
}

void ExportPipeline::ExportSurface(std::string FileName, Exporter_FileTypes FileType, int VolumeID) const
{
    if (this->Mesh == nullptr) {
        LOG("ExportSurface called before mesh was initialized\n", 0);
        return;
    }

    ExportPipeline::ExportSurface(*this->Mesh, FileName, FileType, VolumeID);
}

void ExportPipeline::ExportVolume(std::string FileName, Exporter_FileTypes FileType) const
{
    if (this->Mesh == nullptr) {
        LOG("ExportVolume called before mesh was initialized\n", 0);
        return;
    }

    ExportPipeline::ExportVolume(*this->Mesh, FileName, FileType);
}

void ExportPipeline::ExportPhases(std::string FileName, Exporter_FileTypes FileType) const
{
    if (this->Volumes == nullptr) {
        LOG("ExportPhases called without volumes list\n", 0);
        return;
    }

    for (Volume *v : *this->Volumes) {
        std::stringstream ThisFileName;
        ThisFileName << FileName << "." << v->Phase << ".stl";
        this->ExportSurface(ThisFileName.str(), FileType, v->Phase);
    }
}

void ExportPipeline::ExportSurfaces() const
{
    if (this->Opt->GiveBooleanValue("exportvtksurface")) {
        this->ExportSurface(strfmt("%s.surface.vtp", this->Opt->GiveStringValue("output").c_str()), FT_VTK);
    }

    if (this->Opt->GiveBooleanValue("exportstlsurface")) {
        this->ExportSurface(strfmt("%s.surface.stl", this->Opt->GiveStringValue("output").c_str()), FT_STL);
    }

    if (this->Opt->GiveBooleanValue("exportoff")) {
        this->ExportSurface(strfmt("%s.surface.off", this->Opt->GiveStringValue("output").c_str()), FT_OFF);
    }
}

void ExportPipeline::ExportVolumes() const
{
    if (this->Opt->GiveBooleanValue("exportvtkvolume")) {
        this->ExportVolume(strfmt("%s.volume_0.vtu", this->Opt->GiveStringValue("output").c_str()), FT_VTK);
    }

    if (this->Opt->GiveBooleanValue("exportoofem")) {
        this->ExportVolume(strfmt("%s.in", this->Opt->GiveStringValue("output").c_str()), FT_OOFEM);
    }

    if (this->Opt->GiveBooleanValue("exportabaqus")) {
        this->ExportVolume(strfmt("%s.inp", this->Opt->GiveStringValue("output").c_str()), FT_ABAQUS);
    }

    if (this->Opt->GiveBooleanValue("exportabaqusphon")) {
        this->ExportVolume(strfmt("%s.inp", this->Opt->GiveStringValue("output").c_str()), FT_ABAQUSPHON);
    }
}

void ExportPipeline::ExportPhases() const
{
    if (this->Opt->GiveBooleanValue("exportstlphases")) {
        this->ExportPhases(this->Opt->GiveStringValue("output").c_str(), FT_STL);
    }
}

void ExportPipeline::ExportSurface(MeshData &Mesh, std::string FileName, Exporter_FileTypes FileType, int VolumeID)
{
    STATUS("Export surface to %s\n", FileName.c_str());
    std::unique_ptr<Exporter> exporter = ExportPipeline::CreateSurfaceExporter(Mesh, FileType);
    exporter->WriteSurfaceData(FileName, VolumeID);
}

void ExportPipeline::ExportVolume(MeshData &Mesh, std::string FileName, Exporter_FileTypes FileType)
{
    STATUS("Export volume to %s\n", FileName.c_str());
    std::unique_ptr<Exporter> exporter = ExportPipeline::CreateVolumeExporter(Mesh, FileType);
    if (exporter == nullptr) {
        STATUS("No export filter found for %s\n", FileName.c_str());
        return;
    }
    exporter->WriteVolumeData(FileName);
}

std::unique_ptr<Exporter> ExportPipeline::CreateSurfaceExporter(MeshData &Mesh, Exporter_FileTypes FileType)
{
    switch (FileType) {
        case FT_OFF:
            return std::unique_ptr<Exporter>(new OFFExporter(&Mesh.Triangles, &Mesh.Vertices, &Mesh.Edges, &Mesh.Tets));
        case FT_Poly:
            return std::unique_ptr<Exporter>(new TetGenExporter(&Mesh.Triangles, &Mesh.Vertices, &Mesh.Edges, &Mesh.Tets));
        case FT_VTK:
            return std::unique_ptr<Exporter>(new VTKExporter(&Mesh.Triangles, &Mesh.Vertices, &Mesh.Edges, &Mesh.Tets));
        case FT_SIMPLE:
            return std::unique_ptr<Exporter>(new SimpleExporter(&Mesh.Triangles, &Mesh.Vertices, &Mesh.Edges, &Mesh.Tets));
        case FT_STL:
            return std::unique_ptr<Exporter>(new STLExporter(&Mesh.Triangles, &Mesh.Vertices, &Mesh.Edges, &Mesh.Tets));
        default:
            LOG("File type not supported for surface export\n", 0);
            throw(0);
    }
}

std::unique_ptr<Exporter> ExportPipeline::CreateVolumeExporter(MeshData &Mesh, Exporter_FileTypes FileType)
{
    switch (FileType) {
        case FT_OFF:
            return std::unique_ptr<Exporter>(new OFFExporter(&Mesh.Triangles, &Mesh.Vertices, &Mesh.Edges, &Mesh.Tets));
        case FT_Poly:
            return std::unique_ptr<Exporter>(new TetGenExporter(&Mesh.Triangles, &Mesh.Vertices, &Mesh.Edges, &Mesh.Tets));
        case FT_VTK:
            return std::unique_ptr<Exporter>(new VTKExporter(&Mesh.Triangles, &Mesh.Vertices, &Mesh.Edges, &Mesh.Tets));
        case FT_OOFEM:
            return std::unique_ptr<Exporter>(new OOFEMExporter(&Mesh.Triangles, &Mesh.Vertices, &Mesh.Edges, &Mesh.Tets));
        case FT_ABAQUS:
            return std::unique_ptr<Exporter>(
                    new AbaqusExporter(&Mesh.Triangles, &Mesh.Vertices, &Mesh.Edges, &Mesh.Tets, false));
        case FT_ABAQUSPHON:
            return std::unique_ptr<Exporter>(
                    new AbaqusExporter(&Mesh.Triangles, &Mesh.Vertices, &Mesh.Edges, &Mesh.Tets, true));
        default:
            return nullptr;
    }
}
}
