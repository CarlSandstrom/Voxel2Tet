#ifndef VOXEL2TET_H
#define VOXEL2TET_H

#include <vector>
#include <string>

#include "Options.h"
#include "Importer.h"
#include "CallbackImporter.h"
#include "hdf5DataReader.h"
#include "MiscFunctions.h"
#include "MeshComponents.h"
#include "MeshData.h"
#include "Surface.h"
#include "MiscFunctions.h"
#include "PhaseEdge.h"

namespace voxel2tet
{

class Voxel2Tet
{
private:
    Options *Opt;
    Importer *Imp;
    std::vector <Surface*> Surfaces;
    std::vector <PhaseEdge*> PhaseEdges;
    void FindSurfaces();
    void FindEdges();

    void SmoothEdges();
    void SmoothSurfaces();

    void AddSurfaceSquare(std::vector<int> VoxelIDs, std::vector<int> phases, int normalphase);
    void AddPhaseEdge(std::vector<VertexType*> EdgeSegment, std::vector<int> Phases);

    void FinalizeLoad();

    template <typename T>
    std::vector<int> FindSubsetIndices(std::vector<T> Container, std::vector<T> Subset);

public:
    Voxel2Tet(Options *Opt);
    ~Voxel2Tet();

    MeshData *Mesh;

    void LoadFile(std::string FileName);
    void LoadCallback(cbMaterialIDByCoordinate MaterialIDByCoordinate, std::array<double, 3> origin, std::array<double, 3> spacing, std::array<int, 3> dimensions);

    void LoadData();
    void Process();
};

}

#endif // VOXEL2TET_H
