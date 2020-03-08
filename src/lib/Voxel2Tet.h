#ifndef VOXEL2TET_H
#define VOXEL2TET_H

#include <vector>
#include <string>

#include "Options.h"

#include "Importer.h"
#include "CallbackImporter.h"
#include "Dream3DDataReader.h"
#include "VTKStructuredReader.h"

#include "MiscFunctions.h"
#include "MeshComponents.h"
#include "MeshManipulations.h"
#include "Surface.h"
#include "Volume.h"
#include "MiscFunctions.h"
#include "PhaseEdge.h"
#include "MeshGenerator3D.h"
#include "TimeStamp.h"
#include "Smoother.h"

namespace voxel2tet
{
/**
 * @brief The main class of the library. It supplies functions for loading and exporting data through one function, starting the smoothing process and more overall functions.
 */
class Voxel2TetClass
{
private:
    Options *Opt;
    Importer *Imp;
    std::vector<Surface *> Surfaces;
    std::vector<Volume *> Volumes;
    std::vector<PhaseEdge *> PhaseEdges;

    void FindSurfaces();

    void FindEdges();

    void SmoothEdgesSimultaneously();

    void SmoothSurfaces();

    void AddSurfaceSquare(std::vector<int> VoxelIDs, std::vector<int> phases, int normalphase);

    PhaseEdge *AddPhaseEdge(std::vector<VertexType *> EdgeSegment, std::vector<int> Phases);

    void FinalizeLoad();

    void UpdateSurfaces();

    double eps = 1e-6;

    TimeStamp Timer;

    std::vector<std::vector<double> > PhaseVolumes;
    std::vector<double> CurrentVolumes;
    std::vector<int> PhaseList;

    Smoother *SurfaceSmoother;
    Smoother *EdgeSmoother;
    bool SmoothSimultaneously;

    void UpdateFixed();

public:

    /**
     * @brief Constructor
     * @param Opt Input. Pointer to Options object.
     */
    Voxel2TetClass(Options *Opt);

    ~Voxel2TetClass();

    /**
     * Determines if a 0 is a void or solid
     */
    bool TreatZeroAsVoid = false;

    /**
     * @brief MeshManipulations object for accessing and modifying the mesh.
     */
    MeshManipulations *Mesh;

    /**
     * @brief Prints help message for switch usages
     */
    void PrintHelp();

    /**
     * @brief GetListOfVolumes return a std::vector of volumes for each phase and a vector for the phase
     * @param VolumeList [Out] List of volumes
     * @param PhaseList [Out] List of phases
     * @return Total volume
     */
    double GetListOfVolumes(std::vector<double> &VolumeList, std::vector<int> &PhaseList);

    /**
     * @brief Exports all surfaces to file
     * @param FileName File name of output file
     * @param FileType File type
     */
    void ExportSurface(std::string FileName, Exporter_FileTypes FileType);

    /**
     * @brief Exports all phase boundaries to file
     * @param FileName File name of output file
     * @param FileType File type
     */
    void ExportPhases(std::string FileName, Exporter_FileTypes FileType);

    /**
     * @brief Perform tetrahedralization.
     */
    void Tetrahedralize();

    /**
     * @brief Exports volumes to file
     * @param FileName File name of output file
     * @param FileType File type
     */
    void ExportVolume(std::string FileName, Exporter_FileTypes FileType);

    /**
     * @brief Loads a file
     * @param FileName File name of file to load
     */
    void LoadFile(std::string FileName);

    /**
     * @brief Loads data using a callback function.
     * @param MaterialIDByCoordinate Pointer to function wich returns the material ID at a given point
     * @param origin Origin of data
     * @param spacing Side length of voxels
     * @param dimensions Number of voxels in each dimension
     */
    void LoadCallback(cbMaterialIDByCoordinate MaterialIDByCoordinate, std::array<double, 3> origin,
                      std::array<double, 3> spacing, std::array<int, 3> dimensions);

    /**
     * @brief FindVolumeContainingPoint finds the volume that contains the point P
     * @param P [in] Point
     * @return Pointer to volume object that contains P or NULL if none found.
     */
    Volume *FindVolumeContainingPoint(std::array<double, 3> P);


    /**
     * @brief Performs smoothing process for loaded data. The process includes both smoothing and coarsening of the mesh.
     */
    void Process();

    /**
     * @brief Exports surfaces according to flags in Opt member
     */
    void ExportAllSurfaces();

    /**
     * @brief Exports volumes according to flags in Opt member
     */
    void ExportAllVolumes();

    /**
     * @brief Exports phse boundaries to separate files
     */
    void ExportAllPhases();

    /**
     * @brief Exports statistics on current (finished) job
     */
    void ExportStatistics();
};
}

#endif // VOXEL2TET_H
