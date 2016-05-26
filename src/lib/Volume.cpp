#include "Volume.h"

namespace voxel2tet
{

Volume::Volume()
{
    this->Phase = -1;
}

Volume::Volume(int Phase)
{
    this->Phase = Phase;
}


std::vector<TriangleType*> Volume::GiveTriangles()
{
    std::vector<TriangleType*> VolTriangles;
    for (Surface *s: this->Surfaces) {
        for (TriangleType *t: s->Triangles) {
            VolTriangles.push_back(t);
        }
    }
    return VolTriangles;
}

double Volume::ComputeVolume()
{
    // Reorient triangles such that the normal is consisten over the surface
    std::vector<TriangleType*> VolTriangles = this->GiveTriangles();

    for (Surface *s: this->Surfaces) {

        s->ReorientTriangles();
        double nx = s->ComputeIntegral_nx();
        LOG("nx = %f\n", nx);

    }

}

}
