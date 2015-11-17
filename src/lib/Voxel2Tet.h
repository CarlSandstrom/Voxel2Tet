#ifndef VOXEL2TET_H
#define VOXEL2TET_H

#include <vector>
#include <string>

#include "Options.h"
#include "Importer.h"
#include "hdf5DataReader.h"
#include "MiscFunctions.h"
#include "MeshComponents.h"
#include "MeshData.h"
#include "Surface.h"

namespace voxel2tet
{

class Voxel2Tet
{
private:
    Options *Opt;
    Importer *Imp;
    std::vector <Surface*> Surfaces;
    void FindSurfaces();
    void FindEdges();
    void AddSurfaceSquare(std::vector<int> VoxelIDs, std::vector<int> phases, int normalphase);
public:
    Voxel2Tet(Options *Opt);
    MeshData *Mesh;

    void LoadFile(std::string FileName);
    void LoadData();
    void Process();
};

}

#endif // VOXEL2TET_H
