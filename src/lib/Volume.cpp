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

    double V = 0.0;

    // Use method described in http://n-e-r-v-o-u-s.com/blog/?p=4415
    for (Surface *s: this->Surfaces) {
        s->ReorientTriangles();
        int PosPhase = s->Triangles[0]->PosNormalMatID;
        double ContributionSign = (PosPhase == this->Phase ? -1.0 : 1.0);
        for (TriangleType *t: s->Triangles) {
            double v1x=t->Vertices[0]->get_c(0);
            double v1y=t->Vertices[0]->get_c(1);
            double v1z=t->Vertices[0]->get_c(2);

            double v2x=t->Vertices[1]->get_c(0);
            double v2y=t->Vertices[1]->get_c(1);
            double v2z=t->Vertices[1]->get_c(2);

            double v3x=t->Vertices[2]->get_c(0);
            double v3y=t->Vertices[2]->get_c(1);
            double v3z=t->Vertices[2]->get_c(2);

            double Vc = 1./6.*(v3x*(v1y*v2z-v2y*v1z) + v3y*(v2x*v1z-v1x*v2z) + v3z*(v1x*v2y-v1y*v2x))*ContributionSign;
            V = V + Vc;
        }
    }

    return V;

}

}
