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
    std::vector <Surface*> Surfaces;
    std::vector <Volume*> Volumes;
    std::vector <PhaseEdge*> PhaseEdges;
    void FindSurfaces();
    void FindEdges();

    void SmoothEdgesIndividually(double c, double alpha, double charlength, bool Automatic_c=false);
    void SmoothEdgesSimultaneously(double c, double alpha, double charlength, bool Automatic_c=false);
    void SmoothSurfaces(double c, double alpha, double charlength, bool Automatic_c=false);
    void SmoothAllAtOnce(double c, double alpha, double charlength, bool Automatic_c=false);

    void AddSurfaceSquare(std::vector<int> VoxelIDs, std::vector<int> phases, int normalphase);
    PhaseEdge* AddPhaseEdge(std::vector<VertexType*> EdgeSegment, std::vector<int> Phases);

    void FinalizeLoad();

    template <typename T>
    std::vector<int> FindSubsetIndices(std::vector<T> Container, std::vector<T> Subset);

    void UpdateSurfaces();

    double eps=1e-6;

    double spring_alpha;
    double spring_c;
    double edgespring_alpha;
    double edgespring_c;

    bool auto_c;

    TimeStamp Timer;

public:
    Voxel2TetClass(Options *Opt);
    ~Voxel2TetClass();

    MeshManipulations *Mesh;

    void PrintHelp();

    /**
     * @brief GetListOfVolumes return a std::vector of volumes for each phase and a vector for the phase
     * @param VolumeList [Out] List of volumes
     * @param PhaseList [Out] List of phases
     * @return Total volume
     */
    double GetListOfVolumes(std::vector<double> &VolumeList, std::vector<int> &PhaseList);

    void ExportSurface(std::string FileName, Exporter_FileTypes FileType);

    void Tetrahedralize();
    void ExportVolume(std::string FileName, Exporter_FileTypes FileType);

    void LoadFile(std::string FileName);
    void LoadCallback(cbMaterialIDByCoordinate MaterialIDByCoordinate, std::array<double, 3> origin, std::array<double, 3> spacing, std::array<int, 3> dimensions);

    /**
     * @brief FindVolumeContainingPoint finds the volume that contains the point P
     * @param P [in] Point
     * @return Pointer to volume object that contains P or NULL if none found.
     */
    Volume *FindVolumeContainingPoint (std::array<double, 3> P);

    void LoadData();
    void Process();

    /**
     * @brief Exports surfaces according to flags in Opt member
     */
    void ExportAllSurfaces();

    /**
     * @brief Exports volumes according to flags in Opt member
     */
    void ExportAllVolumes();
};

}

#endif // VOXEL2TET_H
