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


bool Volume::IsPointInside(std::array<double, 3> P)
{
    // Use Ray-Casting method

    // Find a point that definitely is located outside of the volume
    std::array<double, 3> MinPoint = {this->Surfaces[0]->Vertices[0]->get_c(0), this->Surfaces[0]->Vertices[0]->get_c(1), this->Surfaces[0]->Vertices[0]->get_c(2)};
    std::vector<VertexType*> Vertices = this->GiveVertices();

    for (VertexType *v: Vertices) {
        for (int i=0; i<3; i++) {
            if (v->get_c(i)<MinPoint[i]) MinPoint[i]=v->get_c(i);
        }
    }

    for (int i=0; i<3; i++) {
        MinPoint[i] = MinPoint[i] - 10.0-double(i); // The last thing here is to avoid trouble where the line Minpoint-P crosses the edge of triangles
    }

    double O[3] = {MinPoint[0], MinPoint[1], MinPoint[2]};


    double Pa[3] = {P[0], P[1], P[2]};
    double D[3] = {Pa[0]-O[0], Pa[1]-O[1], Pa[2]-O[2]};

    int ti = 0;
    int count = 0;

    for (Surface *s: this->Surfaces) {
        s->ReorientTriangles();
        for (TriangleType *t: s->Triangles) {
            double V0[3] = {t->Vertices[0]->get_c(0), t->Vertices[0]->get_c(1), t->Vertices[0]->get_c(2)};
            double V1[3] = {t->Vertices[1]->get_c(0), t->Vertices[1]->get_c(1), t->Vertices[1]->get_c(2)};
            double V2[3] = {t->Vertices[2]->get_c(0), t->Vertices[2]->get_c(1), t->Vertices[2]->get_c(2)};

            float param;
            if (triangle_ray_intersection(V0, V1, V2, O, D, &param)) {
                //printf ("Intersects with triangle. t=%f\n", param );
                if ( (param > 0.0) & (param <= 1.0)) {
                    ti++;
                }
            }
            count++;

        }
    }

    return (ti % 2)!=0;


}

std::vector<VertexType*> Volume::GiveVertices()
{
    std::vector<VertexType*> VolVertices;
    for(Surface *s: this->Surfaces) {
        for (VertexType *v: s->Vertices) {
            VolVertices.push_back(v);
        }
    }
    return VolVertices;
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
