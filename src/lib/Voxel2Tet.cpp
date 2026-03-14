#include <algorithm>
#include <vector>
#include <iterator>
#include <time.h>
#include <iostream>
#include <fstream>

#include "Voxel2Tet.h"
#include "TetGenCaller.h"
#include "SpringSmoother.h"
#include "SpringSmootherPenalty.h"

namespace voxel2tet
{
Voxel2TetClass::Voxel2TetClass(Options *Opt)
{
    this->Opt = Opt;

    this->InputLoader = std::unique_ptr<InputLoaderPipeline>(new InputLoaderPipeline(this->Opt));
    this->SurfaceBuilder = nullptr;
    this->SmootherPipeline = nullptr;
    this->TetrahedralizationPipeline = nullptr;
    this->OutputPipeline = std::unique_ptr<ExportPipeline>(new ExportPipeline(
            [this]() {
                this->ExportAllSurfacesImpl();
            },
            [this]() {
                this->ExportAllVolumesImpl();
            },
            [this]() {
                this->ExportAllPhasesImpl();
            }));

    // Set defult options

    // Smoothing options
    this->Opt->AddDefaultMap("spring_c_factor", "1");
    this->Opt->AddDefaultMap("spring_alpha", "4");

    this->Opt->AddDefaultMap("edge_spring_c_factor", "1");
    this->Opt->AddDefaultMap("edge_spring_alpha", "4");

    this->Opt->AddDefaultMap("nocoarsening", "0");

    // Dream3D options
    this->Opt->AddDefaultMap("DataContainer", "VoxelDataContainer");
    this->Opt->AddDefaultMap("MaterialId", "GrainIds");

    // Export filters
    this->Opt->AddDefaultMap("exportvtksurface", "1");
    this->Opt->AddDefaultMap("exportvtkvolume", "1");
    this->Opt->AddDefaultMap("exportoff", "0");
    this->Opt->AddDefaultMap("exportoofem", "0");
    this->Opt->AddDefaultMap("exportabaqus", "0");
    this->Opt->AddDefaultMap("exportsteps", "0");

    //Treat zero as void
    if (this->Opt->has_key("treatzeroasvoid")) {
        TreatZeroAsVoid = true;
    } else {
        TreatZeroAsVoid = false;
    }


    // Input/output
    if (!this->Opt->has_key("output")) {
        std::string inputname = this->Opt->GiveStringValue("input");
        size_t lastindex = inputname.find_last_of(".");
        std::string name = inputname.substr(0, lastindex);
        this->Opt->AddDefaultMap("output", name.c_str());
    }

    SmoothSimultaneously = true;

    LOG("Starting Voxel2Tet\n", 0);
}

Voxel2TetClass::~Voxel2TetClass()
{
    for (unsigned int i = 0; i < PhaseEdges.size(); i++) {
        delete this->PhaseEdges.at(i);
    }

    for (auto s : this->Surfaces) {
        delete s;
    }

    delete this->Mesh;

    //delete this->Imp;

    //for (auto p: this->PhaseEdges) delete p;
}

void Voxel2TetClass::PrintHelp()
{
    printf("\nCommand:");
    printf("\n\tVoxel2Tet -input filename [options]");
    printf("\nArguments:");
    printf("\n\tfilename\t name of file containing voxel data");
    printf("\n\toptions");
    printf("\n\t\t-output filename");
    printf("\n\t\t\tFollowed by the name of the output file (without extension)");
    printf("\n\t\t-exportvtksurface");
    printf("\n\t\t\tExport the final surface in VTK format");
    printf("\n\t\t-exportvtkvolume");
    printf("\n\t\t\tExport the final volume in VTK format");
    printf("\n\t\t-exportoff");
    printf("\n\t\t\tExport the funal surface in .OFF format");
    printf("\n\t\t-exportoofem");
    printf("\n\t\t\tExport final volume as input file for OOFEM");
    printf("\n\t\t-exportabaqus");
    printf("\n\t\t\tExport final volume as input file for Abaqus");
    printf("\n\t\t-exportsteps");
    printf("\n\t\t\tExport the result of each step in VTK format (mainly for debugging purposes)");
    printf("\n\t\t-datacontainer name");
    printf("\n\t\t\t(Dream3D input) Name of data group, default 'VoxelDataContainer'");
    printf("\n\t\t-materialid name");
    printf("\n\t\t\t(Dream3D input) Field containing an identifier for the phase, default 'GrainIds'");
    printf("\n\t\t-voxelcutout arg");
    printf("\n\t\t\tOnly consider the subset of the input contained within the boundingbox defined by arg. Here, arg=\"[xmin ymin zmin xmax ymax zmax]\" (include citations and brackets)");
    printf("\n\t\t-treatzeroasvoid");
    printf("\n\t\t\tTreats a material with ID 0 as void. By default, this is considered a solid.");
    printf("\n\t\t-nocoarsening");
    printf("\n\t\t\tSkip coarsening step.");
    printf("\n\t\t-spring_c value");
    printf("\n\t\t\tDetermines how far a vertex can be displaced before penalized. Default depends on voxel size");
    printf("\n\t\t-spring_c_factor value");
    printf("\n\t\t\tChanges the value of spring_c to fit a multiple value of the characteristic length of a voxel. Default is 1.");
    printf("\n\t\t-spring_alpha value");
    printf("\n\t\t\tDetermines penalization for vertices moving too far from the original position. Default is 2.");
    printf("\n\t\t-edge_spring_c value");
    printf("\n\t\t\tSee spring_c flag.");
    printf("\n\t\t-edge_spring_factor value");
    printf("\n\t\t\tSee spring_factor flag. Default is 1.");
    printf("\n\t\t-edge_spring_alpha value");
    printf("\n\t\t\tSee spring_alpha flag. Default is 3.");
    printf("\n\t\t-TOL_FLIP_MAXAREACHANGE value");
    printf("\n\t\t\tLargest change in area due to a flip of a shared edge. Default is 1e-2.");
    printf("\n\t\t-TOL_FLIP_SMALLESTAREA value");
    printf("\n\t\t\tSmallest allowed area of a triangle after a flip. Default is 1e-8.");
    printf("\n\t\t-TOL_FLIP_MAXNORMALCHANGE value");
    printf("\n\t\t\tLargest change in direction of normal of triangles due to a flip (radians). Default is 0.35.");
    printf("\n\t\t-TOL_FLIP_MAXNORMALDIFFERENCE value");
    printf("\n\t\t\tMaximum angle of normals of two triangle to allow flipping (radians). Default is 0.175");
    printf("\n\t\t-TOL_COL_SMALLESTAREA value");
    printf("\n\t\t\tSmallest allowed area of a triangle after collapsing. Default is 1e-8.");
    printf("\n\t\t-TOL_COL_MINANGLE value");
    printf("\n\t\t\tSmallest inner angle of a triangle after collapsing (radians). Default is 0.35.");
    printf("\n\t\t-TOL_COL_MAXNORMALCHANGE value");
    printf("\n\t\t\tMaximum change in normal direction due to a collapse (radians). Default is 0.35.");
    printf("\n\t\t-TOL_COL_CHORD_MAXNORMALCHANGE value");
    printf("\n\t\t\tMaximum change in normal of a chord due to a collapse (radians). Default is 0.175.");
    printf("\n\t\t-TOL_COL_MAXVOLUMECHANGE value");
    printf("\n\t\t\tMaximum change in volume due to one collapse. Default is 2 voxels.");
    printf("\n\t\t-TOL_COL_MAXERROR value");
    printf("\n\t\t\tLargest accumulated error in a vertex. Default is 10 voxels.");
    printf("\n\t\t-TOL_COL_MAXVOLUMECHANGE_FACTOR value");
    printf("\n\t\t\tChanges the default value of TOL_COL_MAXVOLUMECHANGE such that the value is value times the size of a voxel. Default is 1.");
    printf("\n\t\t-TOL_COL_MAXERROR_FACTOR value");
    printf("\n\t\t\tChanges the default value of TOL_COL_MAXERROR such that the value is value times the size of a voxel. Default is 1.");
    printf("\n\n");
}

void Voxel2TetClass::LoadCallback(cbMaterialIDByCoordinate MaterialIDByCoordinate, std::array<double, 3> origin,
                                  std::array<double, 3> spacing, std::array<int, 3> dimensions)
{
    this->Imp = this->InputLoader->LoadCallback(MaterialIDByCoordinate, origin, spacing, dimensions);
    FinalizeLoad();
}

void Voxel2TetClass::LoadFile(std::string Filename)
{
    this->Imp = this->InputLoader->LoadFile(Filename);
    FinalizeLoad();
}

void Voxel2TetClass::FinalizeLoad()
{

    int dim[3];
    double cellspace[3];
    this->Imp->GiveSpacing(cellspace);
    this->Imp->GiveDimensions(dim);

    STATUS("\tVoxel dimensions are %f * %f * %f\n", cellspace[0], cellspace[1], cellspace[2]);
    STATUS("\tNumber of voxels are %u * %u * %u = %u\n", dim[0], dim[1], dim[2], dim[0] * dim[1] * dim[2]);

    this->EdgeSmoother = 0;

    if (this->Opt->has_key("smoothpenalty")) {
        SmoothSimultaneously = true;

        if (this->Opt->has_key("spring_c")) {
            this->SurfaceSmoother = new SpringSmootherPenalty(cellspace[0], Opt->GiveDoubleValue("spring_c"),
                    Opt->GiveDoubleValue("spring_alpha"),
                    Opt->GiveDoubleValue("spring_c_factor"), &this->Volumes, false);
        } else {
            this->SurfaceSmoother = new SpringSmootherPenalty(cellspace[0], Opt->GiveDoubleValue("spring_c"),
                    Opt->GiveDoubleValue("spring_alpha"),
                    Opt->GiveDoubleValue("spring_c_factor"), &this->Volumes, true);
        }

    } else {
        // Setup smoothing classes
        if (this->Opt->has_key("spring_c")) {
            this->SurfaceSmoother = new SpringSmoother(cellspace[0], Opt->GiveDoubleValue("spring_c"),
                    Opt->GiveDoubleValue("spring_alpha"),
                    Opt->GiveDoubleValue("spring_c_factor"), false);
        } else {
            this->SurfaceSmoother = new SpringSmoother(cellspace[0], Opt->GiveDoubleValue("spring_c"),
                    Opt->GiveDoubleValue("spring_alpha"),
                    Opt->GiveDoubleValue("spring_c_factor"), true);
        }

        if (this->Opt->has_key("edge_spring_c")) {
            this->EdgeSmoother = new SpringSmoother(cellspace[0], Opt->GiveDoubleValue("edge_spring_c"),
                    Opt->GiveDoubleValue("edge_spring_alpha"),
                    Opt->GiveDoubleValue("edge_spring_c_factor"), false);
        } else {
            this->EdgeSmoother = new SpringSmoother(cellspace[0], Opt->GiveDoubleValue("edge_spring_c"),
                    Opt->GiveDoubleValue("edge_spring_alpha"),
                    Opt->GiveDoubleValue("edge_spring_c_factor"), true);
        }
        SmoothSimultaneously = false;
    }

    // Setup tolearances in options
    if (!this->Opt->has_key("TOL_FLIP_MAXAREACHANGE")) {
        this->Opt->AddDefaultMap("TOL_FLIP_MAXAREACHANGE", 99);
    }
    if (!this->Opt->has_key("TOL_COL_SMALLESTAREA")) {
        this->Opt->AddDefaultMap("TOL_COL_SMALLESTAREA", 1e-8);
    }
    if (!this->Opt->has_key("TOL_COL_MINANGLE")) {
        this->Opt->AddDefaultMap("TOL_COL_MINANGLE", 10 * 2 * 3.1415 / 360);
    }

    if (!this->Opt->has_key("TOL_COL_MAXNORMALCHANGE")) {
        this->Opt->AddDefaultMap("TOL_COL_MAXNORMALCHANGE", 15 * 2 * 3.1415 / 360);
    }
    if (!this->Opt->has_key("TOL_COL_CHORD_MAXNORMALCHANGE")) {
        this->Opt->AddDefaultMap("TOL_COL_CHORD_MAXNORMALCHANGE", 15 * 2 * 3.1415 / 360);
    }
    if (!this->Opt->has_key("TOL_FLIP_SMALLESTAREA")) {
        this->Opt->AddDefaultMap("TOL_FLIP_SMALLESTAREA", 1e-8);
    }
    if (!this->Opt->has_key("TOL_FLIP_MAXNORMALCHANGE")) {
        this->Opt->AddDefaultMap("TOL_FLIP_MAXNORMALCHANGE", 15 * 2 * 3.1415 / 360);
    }
    if (!this->Opt->has_key("TOL_FLIP_MAXNORMALDIFFERENCE")) {
        this->Opt->AddDefaultMap("TOL_FLIP_MAXNORMALDIFFERENCE", 15 * 2 * 3.1415 / 360);
    }

    if (!this->Opt->has_key("TOL_COL_MAXVOLUMECHANGE_FACTOR")) {
        this->Opt->AddDefaultMap("TOL_COL_MAXVOLUMECHANGE_FACTOR", 5);
    }

    if (!this->Opt->has_key("TOL_COL_MAXERROR_ACCUMULATED_FACTOR")) {
        this->Opt->AddDefaultMap("TOL_COL_MAXERROR_ACCUMULATED_FACTOR", 10);
    }

    if (!this->Opt->has_key("TOL_COL_MAXVOLUMECHANGE")) {
        this->Opt->AddDefaultMap("TOL_COL_MAXVOLUMECHANGE", cellspace[0] * cellspace[1] * cellspace[2] *
                this->Opt->GiveDoubleValue(
                    "TOL_COL_MAXVOLUMECHANGE_FACTOR"));
    }
    if (!this->Opt->has_key("TOL_COL_MAXERROR_ACCUMULATED")) {
        this->Opt->AddDefaultMap("TOL_COL_MAXERROR_ACCUMULATED", cellspace[0] * cellspace[1] * cellspace[2] *
                this->Opt->GiveDoubleValue(
                    "TOL_COL_MAXERROR_ACCUMULATED_FACTOR"));
    }

    // Setup bounding box
    BoundingBoxType bb;

    bb = this->Imp->GiveBoundingBox();

    for (int i = 0; i < 3; i++) {
        bb.maxvalues[i] = bb.maxvalues[i] + cellspace[i];
        bb.minvalues[i] = bb.minvalues[i] - cellspace[i];
    }

    // Create mesh managing object
    Mesh = new MeshManipulations(bb);

    this->Mesh->TOL_COL_MAXVOLUMECHANGE = this->Opt->GiveDoubleValue("TOL_COL_MAXVOLUMECHANGE");
    this->Mesh->TOL_COL_MAXERROR_ACCUMULATED = this->Opt->GiveDoubleValue("TOL_COL_MAXERROR_ACCUMULATED");

    this->SurfaceBuilder = std::unique_ptr<SurfaceBuilderPipeline>(new SurfaceBuilderPipeline(
            this->Imp, this->Mesh, this->Opt, this->SurfaceSmoother, this->EdgeSmoother,
            &this->Surfaces, &this->Volumes, &this->PhaseEdges, this->TreatZeroAsVoid));
    this->SmootherPipeline = std::unique_ptr<SmoothingPipeline>(new SmoothingPipeline(
            this->Imp, this->Mesh, this->EdgeSmoother, &this->Surfaces, &this->PhaseEdges, this->eps));
    this->TetrahedralizationPipeline = std::unique_ptr<TetPipeline>(new TetPipeline(
            this->Imp, &this->Mesh, &this->Volumes, this->eps));

}

void Voxel2TetClass::UpdateSurfaces()
{
    for (Surface *s : this->Surfaces) {
        s->Triangles.clear();
    }

    for (TriangleType *t : this->Mesh->Triangles) {
        this->Surfaces.at(t->InterfaceID)->Triangles.push_back(t);
    }
}

void Voxel2TetClass::ExportSurface(std::string FileName, Exporter_FileTypes FileType)
{
    this->Mesh->ExportSurface(FileName, FileType);
}

void Voxel2TetClass::ExportPhases(std::string FileName, Exporter_FileTypes FileType)
{
    for (Volume *v: this->Volumes) {
        std::stringstream ThisFileName;
        ThisFileName << FileName << "." << v->Phase << ".stl";
        this->Mesh->ExportSurface(ThisFileName.str(), FileType, v->Phase);
    }
}

void Voxel2TetClass::Tetrahedralize()
{
    Timer.StartTimer("Tetrahedralize");
    this->TetrahedralizationPipeline->Tetrahedralize();
    Timer.StopTimer();
}

void Voxel2TetClass::ExportVolume(std::string FileName, Exporter_FileTypes FileType)
{
    this->Mesh->ExportVolume(FileName, FileType);
}

void Voxel2TetClass::UpdateFixed()
{
    for (VertexType *v: this->Mesh->Vertices) {
        for (int j = 0; j < 3; j++) {

            // If the node is located on the surface in the j:th direction, fix it in that direction
            if ((v->get_c(j) > (this->Imp->GiveBoundingBox().maxvalues[j] - eps)) |
                    (v->get_c(j) < (this->Imp->GiveBoundingBox().minvalues[j] + eps))) {
                v->Fixed[j] = true;
            } else {
                v->Fixed[j] = false;
            }

        }

        if (v->IsFixedVertex()) {
            for (int i = 0; i < 3; i++) v->Fixed[i] = true;
        }

    }
}

double Voxel2TetClass::GetListOfVolumes(std::vector<double> &VolumeList, std::vector<int> &PhaseList)
{
    VolumeList.clear();
    PhaseList.clear();

    double TotalVolume = 0.0;
    for (size_t i = 0; i < Volumes.size(); i++) {
        double v = this->Volumes.at(i)->ComputeVolume();
        VolumeList.push_back(v);
        PhaseList.push_back(this->Volumes.at(i)->Phase);
        TotalVolume = +v;
    }
    return TotalVolume;
}

Volume *Voxel2TetClass::FindVolumeContainingPoint(std::array<double, 3> P)
{
    for (Volume *v : this->Volumes) {
        if (v->IsPointInside(P)) {
            return v;
        }
    }
    return NULL;
}

void Voxel2TetClass::Process()
{
    STATUS("Start smoothing process...\n", 0);

    int outputindex = 0;

    // Find surface
    Timer.StartTimer("Find surfaces");
    this->SurfaceBuilder->FindSurfaces();
    Timer.StopTimer();

    // Compute volumes enclosed by surfaces
    Timer.StartTimer("Compute volumes");

    double TotalVolume = GetListOfVolumes(CurrentVolumes, PhaseList);
    LOG("Total volume: %f\n", TotalVolume);
    PhaseVolumes.push_back(CurrentVolumes);
    Timer.StopTimer();

    if (this->Opt->GiveBooleanValue("exportsteps")) {
        //        this->Mesh->ExportSurface(
        //                strfmt("%s_step_%u.simple", this->Opt->GiveStringValue("output").c_str(), outputindex), FT_SIMPLE);
        this->Mesh->ExportSurface(
                    strfmt("%s_step_%u.vtp", this->Opt->GiveStringValue("output").c_str(), outputindex++), FT_VTK);
    }

    Timer.StartTimer("Find edges");
    this->SurfaceBuilder->FindEdges();
    Timer.StopTimer();

    Timer.StartTimer("Smooth edges");
    double Spacing[3];
    this->Imp->GiveSpacing(Spacing);

    if (this->SmoothSimultaneously) {
        UpdateFixed();
        this->SurfaceSmoother->Smooth(this->Mesh->Vertices, this->Mesh);
        if (this->Opt->GiveBooleanValue("exportsteps")) {
            //            this->Mesh->ExportSurface(
            //                    strfmt("%s_step_%u.simple", this->Opt->GiveStringValue("output").c_str(), outputindex), FT_SIMPLE);
            this->Mesh->ExportSurface(
                        strfmt("%s_step_%u.vtp", this->Opt->GiveStringValue("output").c_str(), outputindex++), FT_VTK);
        }
        GetListOfVolumes(CurrentVolumes, PhaseList);
        PhaseVolumes.push_back(CurrentVolumes);

    } else {
        this->SmootherPipeline->SmoothEdges();
        Timer.StopTimer();

#if TEST_MESH_BETWEEN_STEPS_TETGEN == 1
        TetGenCaller Generator;
        Generator.Mesh = this->Mesh;
        Generator.TestMesh();
#endif

        GetListOfVolumes(CurrentVolumes, PhaseList);
        PhaseVolumes.push_back(CurrentVolumes);

        if (this->Opt->GiveBooleanValue("exportsteps")) {
            //            this->Mesh->ExportSurface(
            //                    strfmt("%s_step_%u.simple", this->Opt->GiveStringValue("output").c_str(), outputindex),
            //                    FT_SIMPLE);
            this->Mesh->ExportSurface(
                        strfmt("%s_step_%u.vtp", this->Opt->GiveStringValue("output").c_str(), outputindex++), FT_VTK);
        }

        Timer.StartTimer("Smooth surfaces");
        this->SmootherPipeline->SmoothSurfaces();
        Timer.StopTimer();

        GetListOfVolumes(CurrentVolumes, PhaseList);
        PhaseVolumes.push_back(CurrentVolumes);

        if (this->Opt->GiveBooleanValue("exportsteps")) {
            //            this->Mesh->ExportSurface(
            //                    strfmt("%s_step_%u.simple", this->Opt->GiveStringValue("output").c_str(), outputindex),
            //                    FT_SIMPLE);
            this->Mesh->ExportSurface(
                        strfmt("%s_step_%u.vtp", this->Opt->GiveStringValue("output").c_str(), outputindex++), FT_VTK);
        }

    }

#if TEST_MESH_BETWEEN_STEPS_TETGEN == 1
    Generator.TestMesh();
#endif

    // Perform flipping of edges where it will increase mesh quality
    this->Mesh->FlipAll(false);

    if (this->Opt->GiveBooleanValue("exportsteps")) {
        //        this->Mesh->ExportSurface(
        //                strfmt("%s_step_%u.simple", this->Opt->GiveStringValue("output").c_str(), outputindex), FT_SIMPLE);
        this->Mesh->ExportSurface(
                    strfmt("%s_step_%u.vtp", this->Opt->GiveStringValue("output").c_str(), outputindex++), FT_VTK);
    }

    this->UpdateSurfaces();

    // Check for intersections and pull back vertices in areas where intersections occur
    this->SurfaceSmoother->PullBackAtIntersections(this->Mesh->Vertices, this->Mesh);

    this->Mesh->FlipAll(false);
    this->UpdateSurfaces();

    GetListOfVolumes(CurrentVolumes, PhaseList);
    PhaseVolumes.push_back(CurrentVolumes);

    if (this->Opt->GiveBooleanValue("exportsteps")) {
        //        this->Mesh->ExportSurface(
        //                strfmt("%s_step_%u.simple", this->Opt->GiveStringValue("output").c_str(), outputindex), FT_SIMPLE);
        this->Mesh->ExportSurface(
                    strfmt("%s_step_%u.vtp", this->Opt->GiveStringValue("output").c_str(), outputindex++), FT_VTK);
    }

    // Perform coarsening

    if (!this->Opt->GiveBooleanValue("nocoarsening")) {
        Timer.StartTimer("Coarsen mesh");
        this->Mesh->CoarsenMesh();
        Timer.StopTimer();

        this->UpdateSurfaces();

        GetListOfVolumes(CurrentVolumes, PhaseList);
        PhaseVolumes.push_back(CurrentVolumes);

        if (this->Opt->GiveBooleanValue("exportsteps")) {
            //        this->Mesh->ExportSurface(
            //                strfmt("%s_step_%u.simple", this->Opt->GiveStringValue("output").c_str(), outputindex), FT_SIMPLE);
            this->Mesh->ExportSurface(
                        strfmt("%s_step_%u.vtp", this->Opt->GiveStringValue("output").c_str(), outputindex++), FT_VTK);
        }
        this->Mesh->FlipAll();
        this->UpdateSurfaces();
    }

    GetListOfVolumes(CurrentVolumes, PhaseList);
    PhaseVolumes.push_back(CurrentVolumes);

    if (this->Opt->GiveBooleanValue("exportsteps")) {
        //        this->Mesh->ExportSurface(
        //                strfmt("%s_step_%u.simple", this->Opt->GiveStringValue("output").c_str(), outputindex), FT_SIMPLE);
        this->Mesh->ExportSurface(
                    strfmt("%s_step_%u.vtp", this->Opt->GiveStringValue("output").c_str(), outputindex++), FT_VTK);
    }

#if EXPORT_MESH_COARSENING == 1
    dooutputlogmesh(* this->Mesh, ( char * ) "/tmp/finalcoarsening%u.vtp", 0);
#endif
    this->Timer.PrintTable();
}

void Voxel2TetClass::ExportAllSurfaces()
{
    this->OutputPipeline->ExportSurfaces();
}

void Voxel2TetClass::ExportAllSurfacesImpl()
{
    // Surfaces
    if (this->Opt->GiveBooleanValue("exportvtksurface")) {
        this->Mesh->ExportSurface(strfmt("%s.surface.vtp", this->Opt->GiveStringValue("output").c_str()), FT_VTK);
    }

    if (this->Opt->GiveBooleanValue("exportstlsurface")) {
        this->Mesh->ExportSurface(strfmt("%s.surface.stl", this->Opt->GiveStringValue("output").c_str()), FT_STL);
    }

    if (this->Opt->GiveBooleanValue("exportoff")) {
        this->Mesh->ExportSurface(strfmt("%s.surface.off", this->Opt->GiveStringValue("output").c_str()), FT_OFF);
    }
}

void Voxel2TetClass::ExportAllVolumes()
{
    this->OutputPipeline->ExportVolumes();
}

void Voxel2TetClass::ExportAllVolumesImpl()
{
    // Volumes
    if (this->Opt->GiveBooleanValue("exportvtkvolume")) {
        this->Mesh->ExportVolume(strfmt("%s.volume_0.vtu", this->Opt->GiveStringValue("output").c_str()), FT_VTK);
        //this->Mesh->ExportVolume(strfmt( "%s.volume_1.vtu", this->Opt->GiveStringValue("output").c_str() ), FT_VTK);
    }

    if (this->Opt->GiveBooleanValue("exportoofem")) {
        this->Mesh->ExportVolume(strfmt("%s.in", this->Opt->GiveStringValue("output").c_str()), FT_OOFEM);
    }

    if (this->Opt->GiveBooleanValue("exportabaqus")) {
        this->Mesh->ExportVolume(strfmt("%s.inp", this->Opt->GiveStringValue("output").c_str()), FT_ABAQUS);
    }
    if (this->Opt->GiveBooleanValue("exportabaqusphon")) {
        this->Mesh->ExportVolume(strfmt("%s.inp", this->Opt->GiveStringValue("output").c_str()), FT_ABAQUSPHON);
    }
}

void Voxel2TetClass::ExportAllPhases()
{
    this->OutputPipeline->ExportPhases();
}

void Voxel2TetClass::ExportAllPhasesImpl()
{
    if (this->Opt->GiveBooleanValue("exportstlphases")) {
        this->ExportPhases(this->Opt->GiveStringValue("output").c_str(), FT_STL);
    }
}

void Voxel2TetClass::ExportStatistics()
{
    int dimensions[3];
    this->Imp->GiveDimensions(dimensions);
    double origin[3];
    this->Imp->GiveOrigin(origin);
    double spacing[3];
    this->Imp->GiveSpacing(spacing);

    std::string Filename = strfmt("%s.stat", this->Opt->GiveStringValue("output").c_str());
    STATUS("Write statistics output file %s\n", Filename.c_str());
    std::ofstream StatFile;
    StatFile.open(Filename);
    StatFile << "Voxel2Tet statistics\n";
    StatFile << "====================\n\n";
    StatFile << "Input file:\t" << this->Opt->GiveStringValue("input") << "\n";
    StatFile << "\tDimensions: " << dimensions[0] << "*" << dimensions[1] << "*" << dimensions[2] << " voxels\n";
    StatFile << "\tOrigin:     " << origin[0] << "*" << origin[1] << "*" << origin[2] << "\n";
    StatFile << "\tSpacing:    " << spacing[0] << "*" << spacing[1] << "*" << spacing[2] << "\n";
    StatFile << "Output file:\t" << this->Opt->GiveStringValue("output") << "\n";

    StatFile << "\nConstants used during smoothing\n";
    StatFile << "-------------------------------\n";
    StatFile << "Smoothing algorithm:\n";
    StatFile << "Surface smoother:\n";
    StatFile << *this->SurfaceSmoother;
    if (this->EdgeSmoother != 0) {
        StatFile << "Edge smoother:\n";
        StatFile << *this->EdgeSmoother;
    }

    StatFile << "\nMesh coarsening:\n";
    StatFile << "TOL_MAXAREACHANGE = " << this->Opt->GiveStringValue("TOL_MAXAREACHANGE") << "\n";
    StatFile << "TOL_COL_SMALLESTAREA = " << this->Opt->GiveStringValue("TOL_COL_SMALLESTAREA") << "\n";
    StatFile << "TOL_COL_MAXNORMALCHANGE = " << this->Opt->GiveStringValue("TOL_COL_MAXNORMALCHANGE") << "\n";
    StatFile << "TOL_COL_CHORD_MAXNORMALCHANGE = " << this->Opt->GiveStringValue("TOL_COL_CHORD_MAXNORMALCHANGE")
             << "\n";
    StatFile << "TOL_FLIP_SMALLESTAREA = " << this->Opt->GiveStringValue("TOL_FLIP_SMALLESTAREA") << "\n";
    StatFile << "TOL_FLIP_MAXNORMALCHANGE = " << this->Opt->GiveStringValue("TOL_FLIP_MAXNORMALCHANGE") << "\n";
    StatFile << "TOL_FLIP_MAXNORMALDIFFERENCE = " << this->Opt->GiveStringValue("TOL_FLIP_MAXNORMALDIFFERENCE")
             << "\n";
    StatFile << "TOL_COL_MAXVOLUMECHANGE = " << this->Opt->GiveStringValue("TOL_COL_MAXVOLUMECHANGE") << "\n";
    StatFile << "TOL_COL_MAXERROR_ACCUMULATED = " << this->Opt->GiveStringValue("TOL_COL_MAXERROR_ACCUMULATED")
             << "\n";
    StatFile << "TOL_COL_MAXVOLUMECHANGE_FACTOR = " << this->Opt->GiveStringValue("TOL_COL_MAXVOLUMECHANGE_FACTOR")
             << "\n";
    StatFile << "TOL_COL_MAXERROR_ACCUMULATED_FACTOR = "
             << this->Opt->GiveStringValue("TOL_COL_MAXERROR_ACCUMULATED_FACTOR") << "\n";

    StatFile << "\nMesh\n----\n";
    StatFile << "Number of input voxels: " << dimensions[0] * dimensions[1] * dimensions[2] << "\n";
    StatFile << "Number of input nodes: " << (dimensions[0] + 1) * (dimensions[1] + 1) * (dimensions[2] + 1)
            << "\n";
    StatFile << "Number of output elements: " << this->Mesh->Tets.size() << "\n";
    StatFile << "Number of output nodes: " << this->Mesh->Vertices.size() << "\n";

    StatFile << "\nVolume change (One column per phase):\nPhase\t";

    for (int p : PhaseList) {
        StatFile << p << "\t";
    }

    double Vtot = 0;
    double E = 0;

    StatFile << "\nVoxel\t";
    for (size_t j = 0; j < PhaseVolumes.at(0).size(); j++) {
        StatFile << PhaseVolumes.at(0).at(j) << "\t";
        Vtot = Vtot + PhaseVolumes.at(0).at(j);
    }

    StatFile << "\nSmooth\t";
    for (size_t j = 0; j < PhaseVolumes.at(1).size(); j++) {
        StatFile << PhaseVolumes.at(1).at(j) << "\t";
    }

    StatFile << "\nTet  \t";
    for (size_t j = 0; j < PhaseVolumes.at(PhaseVolumes.size() - 1).size(); j++) {
        StatFile << PhaseVolumes.at(PhaseVolumes.size() - 1).at(j) << "\t";
    }

    StatFile << "\n%    \t";

    for (size_t j = 0; j < PhaseVolumes.at(PhaseVolumes.size() - 1).size(); j++) {
        double change = PhaseVolumes.at(PhaseVolumes.size() - 1).at(j) / PhaseVolumes.at(0).at(j);
        StatFile << (change - 1) * 100 << "\t";
    }

    for (size_t i = 0; i < PhaseVolumes.at(0).size(); i++) {
        double ei = fabs(PhaseVolumes.at(0).at(i) - PhaseVolumes.at(PhaseVolumes.size() - 1).at(i)) *
                PhaseVolumes.at(0).at(i);
        E = E + ei / Vtot;
    }

    StatFile << "\n\nWeighted error: " << E;

    StatFile << "\n\nTime\n----\n";

    std::vector<std::pair<double, std::string> > TimeTable = this->Timer.GetTable();
    for (std::pair<double, std::string> stamp : TimeTable) {
        StatFile << stamp.second << "\t" << stamp.first << "s\n";
    }
}
}
