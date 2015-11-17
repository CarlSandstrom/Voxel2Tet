#include "Voxel2Tet.h"
#include "Importer.h"
#include "hdf5DataReader.h"

namespace voxel2tet
{

Voxel2Tet::Voxel2Tet(Options *Opt)
{
    this->Opt = Opt;
}

void Voxel2Tet::LoadFile(std::string Filename)
{
    STATUS ("Load file %s\n", Filename.c_str());
    hdf5DataReader *DataReader = new hdf5DataReader();
    DataReader->LoadFile(Filename);
    this->Imp = DataReader;
}

void Voxel2Tet::FindSurfaces()
{

    STATUS("Find surfaces\n", 0);

    int dim[3];
    double spacing[3], origin[3];
    this->Imp->GiveDimensions(dim);
    this->Imp->GiveSpacing(spacing);
    this->Imp->GiveOrigin(origin);

    std::vector <double> signs = {1, -1};

    for (int i=0; i<dim[0]; i++) {
        for (int j=0; j<dim[1]; j++) {
            for (int k=0; k<dim[2]; k++) {

                // Specifies directions in which we will look for different materials
                std::vector <std::vector <double>> testdirections = {{1,0,0}, {0,1,0}, {0,0,1}};
                // If we the adjacent material is of other type, we will create a square by varying the coordinates marked '1' in vdirections
                std::vector <std::vector <double>> vdirections = {{0,1,1}, {1,0,1}, {1,1,0}};
                // vindex is the indices of the coordinates
                std::vector <std::vector <int>> vindex = {{1,2},{0,2},{0,1}};

                // If we are on a boundary, we need to check what is outside of that boundary
                if (i==0) {
                    testdirections.push_back({-1,0,0});
                    vdirections.push_back({0,1,1});
                    vindex.push_back({1,2});
                }
                if (j==0) {
                    testdirections.push_back({0,-1,0});
                    vdirections.push_back({1,0,1});
                    vindex.push_back({0,2});
                }
                if (k==0) {
                    testdirections.push_back({0,0,-1});
                    vdirections.push_back({1,1,0});
                    vindex.push_back({0,1});
                }

                int ThisPhase = this->Imp->GiveMaterialIDByIndex(i, j, k);
                int NeighboringPhase;
                bool SamePhase;

                // Check material in each direction
                for (unsigned int m=0; m<testdirections.size(); m++) {
                    int testi = testdirections.at(m).at(0) + i;
                    int testj = testdirections.at(m).at(1) + j;
                    int testk = testdirections.at(m).at(2) + k;

                    // If comparing inside the domain, simply compare
                    if ( (testi>=0) & (testj>=0) & (testk>=0) & (testi<dim[0]) & (testj<dim[1]) & (testk<dim[2])) {
                        NeighboringPhase = this->Imp->GiveMaterialIDByIndex(testi, testj, testk);
                        SamePhase = (ThisPhase == NeighboringPhase);
                    } else {
                        // If we are comparing with the outside, take into account that a we might have void (i.e. 0) in both voxels
                        if (ThisPhase!=0) {
                            SamePhase = false;
                            NeighboringPhase = this->Imp->GiveMaterialIDByIndex(testi, testj, testk);
                        } else {
                            // Void-to-void connection
                            SamePhase = true;
                            NeighboringPhase = ThisPhase;
                        }
                    }

                    if (!SamePhase) {   // Add surface. I.e. a square separating the two voxels

                        // Compute centre off square
                        double c[3];
                        c[0] = (double(i) + double(testdirections.at(m).at(0))/2.0) * spacing[0] + origin[0] + spacing[0]/2.0;
                        c[1] = (double(j) + double(testdirections.at(m).at(1))/2.0) * spacing[1] + origin[1] + spacing[1]/2.0;
                        c[2] = (double(k) + double(testdirections.at(m).at(2))/2.0) * spacing[2] + origin[2] + spacing[2]/2.0;

                        // Compute coordinate of corner point
                        double delta[3];
                        delta[0] = spacing[0]*vdirections.at(m)[0] / 2.0;
                        delta[1] = spacing[1]*vdirections.at(m)[1] / 2.0;
                        delta[2] = spacing[2]*vdirections.at(m)[2] / 2.0;

                        std::vector <int> VoxelIDs;

                        for (auto s1: signs) {
                            for (auto s2: signs) {
                                double newvertex[3];
                                newvertex[0] = c[0];
                                newvertex[1] = c[1];
                                newvertex[2] = c[2];

                                newvertex[vindex.at(m).at(0)] = newvertex[vindex.at(m).at(0)] + s1*delta[vindex.at(m)[0]];
                                newvertex[vindex.at(m).at(1)] = newvertex[vindex.at(m).at(1)] + s2*delta[vindex.at(m)[1]];

                                int id = Mesh->VertexOctreeRoot->AddVertex(newvertex[0], newvertex[1], newvertex[2]);
                                LOG ("Corner (id=%u) at (%f, %f, %f)\n", id, newvertex[0], newvertex[1], newvertex[2]);
                                VoxelIDs.push_back(id);
                            }
                        }
                        AddSurfaceSquare(VoxelIDs, {ThisPhase, NeighboringPhase}, NeighboringPhase);
                    }
                }
            }
        }
    }
}

void Voxel2Tet :: FindEdges()
{

    STATUS ("Find edges\n", 0);

    std::vector <VertexType*> SharedVertices;

    // Find all vertices that are shared among the surfaces
    for (auto surface1: this->Surfaces) {
        for (auto surface2: this->Surfaces) {

            if (surface1==surface2) break;

            for (int i=0; i<2; i++) {
                int mat1 = surface1->Phases[i];
                for (int j=0; j<2; j++) {
                    int mat2 = surface2->Phases[i];
                    LOG("Surfaces %p and %p shares phases\n", surface1, surface2);
                    LOG("   surface1->Phases = [%i, %i]\n", surface1->Phases[0], surface1->Phases[1]);
                    LOG("   surface2->Phases = [%i, %i]\n", surface2->Phases[0], surface2->Phases[1]);
                }
            }
        }
    }
}

void Voxel2Tet :: AddSurfaceSquare(std::vector<int> VoxelIDs, std::vector<int> phases, int normalphase)
{
    // Check is surface exists
    Surface *ThisSurface = NULL;
    int SurfaceID;
    for (unsigned int i=0; i<this->Surfaces.size(); i++) {
        if ( ( (this->Surfaces.at(i)->Phases[0] == phases.at(0) ) & ( this->Surfaces.at(i)->Phases[1] == phases.at(1) ) ) |
             ( (this->Surfaces.at(i)->Phases[0] == phases.at(1) ) & ( this->Surfaces.at(i)->Phases[1] == phases.at(0) ) ) ) {
            ThisSurface = this->Surfaces.at(i);
            SurfaceID = i;
            break;
        }
    }

    // If not, create it and add it to the list
    if (ThisSurface==NULL) {
        ThisSurface = new Surface(phases.at(0), phases.at(1));
        this->Surfaces.push_back(ThisSurface);
        SurfaceID = this->Surfaces.size()-1;
    }

    // Create square (i.e. two triangles)
    TriangleType *triangle0, *triangle1;
    triangle0 = Mesh->AddTriangle({VoxelIDs.at(0), VoxelIDs.at(1), VoxelIDs.at(2)});
    triangle1 = Mesh->AddTriangle({VoxelIDs.at(1), VoxelIDs.at(3), VoxelIDs.at(2)});
    triangle0->InterfaceID = SurfaceID;
    triangle1->InterfaceID = SurfaceID;
    triangle0->PosNormalMatID = normalphase;
    triangle1->PosNormalMatID = normalphase;

    // Update surface
    ThisSurface->AddTriangle(triangle0);
    ThisSurface->AddTriangle(triangle0);
    for (int i=0; i<3; i++) {
        ThisSurface->AddVertex(this->Mesh->Vertices.at(VoxelIDs.at(i)));
    }

}

void Voxel2Tet::LoadData()
{
    STATUS ("Load data\n", 0);
    if (this->Opt->has_key("i")) {
        LoadFile(this->Opt->GiveStringValue("i"));
    }
    BoundingBoxType bb;
    double spacing[3];

    bb = this->Imp->GiveBoundingBox();
    this->Imp->GiveSpacing(spacing);

    for (int i=0; i<3; i++) {
        bb.maxvalues[i] = bb.maxvalues[i] + spacing[i];
        bb.minvalues[i] = bb.minvalues[i] - spacing[i];
    }
    Mesh = new MeshData(bb);
}

void Voxel2Tet::Process()
{
    STATUS ("Proccess content\n", 0);
    this->FindSurfaces();

    this->FindEdges();

}

}
