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
    hdf5DataReader *DataReader = new hdf5DataReader();
    DataReader->LoadFile(Filename);
    this->Imp = DataReader;
}

void Voxel2Tet::FindSurfaces()
{

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
                std::vector <std::vector <double>> vdirections = {{1,0,0}, {0,1,0}, {0,0,1}};
                // vindex is the indices of the coordinates
                std::vector <std::vector <int>> vindex = {{1,2},{0,1},{0,1}};

                // If we are on a boundary, we need to check what is outside of that boundary
                if (i==0) {
                    testdirections.push_back({-1,0,0});
                    vdirections.push_back({0,1,1});
                    vindex.push_back({1,2});
                } else if (i==1) {
                    testdirections.push_back({-1,0,0});
                    vdirections.push_back({0,1,1});
                    vindex.push_back({1,2});
                } else if (i==2) {
                    testdirections.push_back({-1,0,0});
                    vdirections.push_back({0,1,1});
                    vindex.push_back({1,2});
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
                        log("c=%f, %f, %f", c[0], c[1], c[2]);

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
                                for (int n=0; n<2; n++) {
                                    newvertex[vindex.at(m).at(n)] = newvertex[vindex.at(m).at(n)] + ;
                                }
                                VoxelIDs.push_back(newvertex[0]);
                            }
                        }
                    }
                }
            }
        }
    }
}

void Voxel2Tet::LoadData()
{
    if (this->Opt->has_key("i")) {
        LoadFile(this->Opt->GiveStringValue("i"));
    }
}

void Voxel2Tet::Process()
{
    log ("Proccess content\n", 0);
    this->FindSurfaces();
}

}
