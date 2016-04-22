#ifndef TRITRIINTERSECT_H
#define TRITRIINTERSECT_H

int coplanar_tri_tri(double N[3],double V0[3],double V1[3],double V2[3], double U0[3],double U1[3],double U2[3]);
int tri_tri_intersect(double V0[3],double V1[3],double V2[3], double U0[3],double U1[3],double U2[3]);

bool point_in_tri(double V0[3],double V1[3],double V2[3], double P[3]);
bool tri_tri_intersect_shared_edge(double s0[3], double s1[3], double u0[3], double u1[3]);

#endif // TRITRIINTERSECT_H
